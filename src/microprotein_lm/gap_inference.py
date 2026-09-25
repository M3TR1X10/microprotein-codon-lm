"""Frozen causal-model inference for prospective codon-gap experiments.

The external cache preserves the original model and absolute nucleotide positions.
No hidden target or right flank is accepted by the proposal methods.
"""
from dataclasses import dataclass
import math

import torch
from torch.nn import functional as F

from .tokenization import Tokenizer


@dataclass
class State:
    logits: torch.Tensor
    layers: tuple
    length: int
    offset_tokens: int

    def select(self, indices):
        indices = torch.as_tensor(indices, dtype=torch.long, device=self.logits.device)
        return State(self.logits.index_select(0, indices),
                     tuple((k.index_select(0, indices), v.index_select(0, indices))
                           for k, v in self.layers), self.length, self.offset_tokens)


class GapDecoder:
    """Codon-synchronous, unconstrained decoding of a fixed pretrained model."""

    def __init__(self, model, mode):
        self.model = model.eval()
        self.tokenizer = Tokenizer(mode)
        self.codons = Tokenizer('codon').vocabulary
        if model.config.vocab_size != len(self.tokenizer.vocabulary):
            raise ValueError('Model and tokenizer disagree')
        if model.config.token_width != self.tokenizer.width:
            raise ValueError('Model and token width disagree')
        self.device = next(model.parameters()).device

    @torch.inference_mode()
    def advance(self, tokens, state=None, offset_nt=0):
        """Return cached state and all new logits; retain absolute CDS coordinates."""
        xids = torch.as_tensor(tokens, dtype=torch.long, device=self.device)
        if xids.ndim != 2 or not xids.shape[0] or not xids.shape[1]:
            raise ValueError('Expected a nonempty rectangular batch')
        if int(xids.min()) < 0 or int(xids.max()) >= self.model.config.vocab_size:
            raise ValueError('Token outside vocabulary')
        if type(offset_nt) is not int or offset_nt < 0 or offset_nt % 3:
            raise ValueError('Context must begin at a nonnegative complete-codon boundary')
        old_length = state.length if state else 0
        offset = state.offset_tokens if state else offset_nt // self.tokenizer.width
        if state and offset_nt and offset != offset_nt // self.tokenizer.width:
            raise ValueError('Cached coordinate offset cannot change')
        end = offset + old_length + xids.shape[1]
        if end > self.model.config.max_tokens:
            raise ValueError('Requested input exceeds trained positional capacity')
        if state and xids.shape[0] != state.logits.shape[0]:
            raise ValueError('Cache and token batch sizes differ')
        x = self.model.embedding(xids) + self.model.position[offset + old_length:end]
        new_layers = []
        for i, block in enumerate(self.model.blocks):
            h = block.ln1(x)
            b, t, d = h.shape
            heads = block.attn.heads
            q, k, v = [z.view(b, t, heads, d // heads).transpose(1, 2)
                       for z in block.attn.qkv(h).chunk(3, dim=-1)]
            if state:
                k = torch.cat((state.layers[i][0], k), dim=2)
                v = torch.cat((state.layers[i][1], v), dim=2)
            new_layers.append((k, v))
            if old_length:
                allowed = torch.arange(k.shape[2], device=self.device)[None, :] <= (
                    old_length + torch.arange(t, device=self.device)[:, None])
                y = F.scaled_dot_product_attention(q, k, v, attn_mask=allowed,
                                                  dropout_p=0.0, is_causal=False)
            else:
                y = F.scaled_dot_product_attention(q, k, v, dropout_p=0.0,
                                                  is_causal=True)
            y = y.transpose(1, 2).contiguous().view(b, t, d)
            x = x + block.attn.proj(y)
            x = x + block.mlp(block.ln2(x))
        logits = self.model.head(self.model.norm(x))
        return State(logits[:, -1], tuple(new_layers), old_length + xids.shape[1], offset), logits

    def prefix(self, rna, offset_nt=0):
        if len(rna) % 3:
            raise ValueError('Prefix must contain whole codons')
        return self.advance([self.tokenizer.encode(rna)], offset_nt=offset_nt)[0]

    @torch.inference_mode()
    def next_codons(self, state):
        """Exact 64-way conditional distribution, including stops under either code.

        In base mode sum three conditional log probabilities over the full 4^3
        tree. Return the two-base branch cache for efficient selected continuation.
        """
        first = F.log_softmax(state.logits.double(), dim=-1)
        if self.tokenizer.mode == 'codon':
            return first, state
        b = state.logits.shape[0]
        ids = torch.arange(4, device=self.device).repeat(b)
        one, _ = self.advance(ids[:, None], state.select(torch.arange(b, device=self.device).repeat_interleave(4)))
        second = F.log_softmax(one.logits.double(), dim=-1).reshape(b, 4, 4)
        two, _ = self.advance(torch.arange(4, device=self.device).repeat(b * 4)[:, None],
                              one.select(torch.arange(b * 4, device=self.device).repeat_interleave(4)))
        third = F.log_softmax(two.logits.double(), dim=-1).reshape(b, 4, 4, 4)
        joint = first[:, :, None, None] + second[:, :, :, None] + third
        return joint.reshape(b, 64), two

    def append_codons(self, branch_state, parents, codon_ids):
        parents = torch.as_tensor(parents, dtype=torch.long, device=self.device)
        codon_ids = torch.as_tensor(codon_ids, dtype=torch.long, device=self.device)
        if self.tokenizer.mode == 'codon':
            selected, ids = parents, codon_ids
        else:
            selected, ids = parents * 16 + codon_ids // 4, codon_ids % 4
        return self.advance(ids[:, None], branch_state.select(selected))[0]

    @torch.inference_mode()
    def proposals(self, prefix, gap_codons, width=1, offset_nt=0):
        """Fixed-width codon beam; width=1 is joint-codon MAP greedy.

        Ties follow AUCG vocabulary order. Every sequence has the requested length;
        premature stops remain measurable outcomes and never terminate decoding.
        """
        if type(gap_codons) is not int or gap_codons < 1:
            raise ValueError('Gap must contain a positive number of codons')
        if type(width) is not int or not 1 <= width <= 64:
            raise ValueError('Beam width must lie in 1..64')
        total_input_nt = offset_nt + len(prefix) + 3 * gap_codons - self.tokenizer.width
        if total_input_nt // self.tokenizer.width > self.model.config.max_tokens:
            raise ValueError('Gap exceeds the trained positional capacity')
        state = self.prefix(prefix, offset_nt)
        scores = torch.zeros(1, dtype=torch.float64, device=self.device)
        sequences = ['']
        paths = [()]
        for step in range(gap_codons):
            probabilities, branches = self.next_codons(state)
            flat = (scores[:, None] + probabilities).flatten()
            values = flat.tolist()
            # Parent rank is a score order, not lexical order. Tie-break on the
            # entire codon-index path, including when parent scores differ.
            ranked = sorted(range(len(values)), key=lambda j: (
                -values[j], paths[j // 64] + (j % 64,)))[:min(width, len(flat))]
            chosen = torch.tensor(ranked, dtype=torch.long, device=self.device)
            parents, codon_ids = chosen // 64, chosen % 64
            sequences = [sequences[p] + self.codons[c] for p, c in zip(parents.tolist(), codon_ids.tolist())]
            paths = [paths[p] + (c,) for p, c in zip(parents.tolist(), codon_ids.tolist())]
            scores = flat[chosen]
            if step + 1 < gap_codons:
                state = self.append_codons(branches, parents, codon_ids)
        return [{'rna': s, 'log_probability': float(v)} for s, v in zip(sequences, scores.tolist())]

    @torch.inference_mode()
    def continuation_log_probability(self, prefix, continuation, offset_nt=0):
        """Teacher-forced probability of a supplied continuation, without lookahead."""
        prefix_ids = self.tokenizer.encode(prefix)
        targets = self.tokenizer.encode(continuation)
        if len(prefix) % 3 or len(continuation) % 3:
            raise ValueError('Both boundaries must be codon aligned')
        inputs = prefix_ids + targets[:-1]
        _, logits = self.advance([inputs], offset_nt=offset_nt)
        logp = F.log_softmax(logits[0, len(prefix_ids) - 1:].double(), dim=-1)
        return float(logp[torch.arange(len(targets), device=self.device),
                          torch.tensor(targets, device=self.device)].sum())

    def rerank(self, prefix, candidates, right_flank, offset_nt=0):
        """Score a fixed prefix-only bank by log P(gap,right|prefix).

        Right-flank likelihoods are never normalized into a purported posterior.
        No true gap is accepted. Candidate generation has already finished.
        """
        if not candidates or len({len(c['rna']) for c in candidates}) != 1:
            raise ValueError('Reranking requires a nonempty equal-length bank')
        if not candidates[0]['rna'] or len(candidates[0]['rna']) % 3:
            raise ValueError('Candidates must contain complete codons')
        if len({c['rna'] for c in candidates}) != len(candidates):
            raise ValueError('Candidate bank must be distinct')
        scored = []
        for candidate in candidates:
            if not math.isfinite(candidate['log_probability']):
                raise ValueError('Nonfinite proposal score')
            self.tokenizer.encode(candidate['rna'])
            suffix_logp = self.continuation_log_probability(prefix + candidate['rna'], right_flank, offset_nt) if right_flank else 0.0
            scored.append({**candidate, 'right_log_probability': suffix_logp,
                           'joint_log_score': candidate['log_probability'] + suffix_logp})
        return sorted(scored, key=lambda c: -c['joint_log_score'])
