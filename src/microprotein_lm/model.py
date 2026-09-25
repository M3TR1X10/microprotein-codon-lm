"""Small decoder model with Hybrid RoPE + Frame Positional Encodings
and Multi-Task Pretraining Objectives (Autoregressive, FIM, Span Corruption, Frame).

See references/SLM-LICENSE and docs/scaffold.md for attribution and changes.
"""
import math
from dataclasses import dataclass, field
from typing import Optional, Dict, Tuple, Any

import torch
from torch import nn
from torch.nn import functional as F


@dataclass
class ModelConfig:
    vocab_size: int
    max_tokens: int
    token_width: int
    d_model: int = 48
    n_heads: int = 4
    n_layers: int = 2
    dropout: float = 0.1
    use_rope: bool = False
    use_frame_encoding: bool = True
    num_frames: int = 4
    fim_rate: float = 0.5
    span_mask_rate: float = 0.15

    def __post_init__(self):
        if self.vocab_size not in (4, 64) or self.token_width not in (1, 3):
            raise ValueError('Only fixed base and codon vocabularies are supported')
        if (self.vocab_size, self.token_width) not in ((4, 1), (64, 3)):
            raise ValueError('Vocabulary and token width disagree')
        if self.d_model % self.n_heads or self.d_model % 2 or self.d_model < 2:
            raise ValueError('d_model must be even and divisible by n_heads')
        if self.n_layers < 1 or self.max_tokens < 1 or not 0 <= self.dropout < 1:
            raise ValueError('Invalid model dimensions or dropout')
        if self.num_frames < 1:
            raise ValueError('num_frames must be >= 1')


class RotaryPositionalEncoding(nn.Module):
    """Rotary Position Embedding (RoPE) for sequence positions."""

    def __init__(self, dim: int, max_tokens: int = 2048, base: float = 10000.0, token_width: int = 1):
        super().__init__()
        self.dim = dim
        self.base = base
        self.token_width = token_width
        inv_freq = 1.0 / (base ** (torch.arange(0, dim, 2).float() / dim))
        self.register_buffer("inv_freq", inv_freq, persistent=False)
        self._build_cache(max_tokens)

    def _build_cache(self, max_tokens: int):
        pos = torch.arange(max_tokens, dtype=torch.float32) * self.token_width
        freqs = torch.outer(pos, self.inv_freq)
        emb = torch.cat((freqs, freqs), dim=-1)
        self.register_buffer("cos_cached", emb.cos(), persistent=False)
        self.register_buffer("sin_cached", emb.sin(), persistent=False)

    def forward(self, seq_len: int, device: torch.device) -> Tuple[torch.Tensor, torch.Tensor]:
        if seq_len > self.cos_cached.shape[0]:
            self._build_cache(seq_len)
        return self.cos_cached[:seq_len].to(device), self.sin_cached[:seq_len].to(device)


def apply_rope(x: torch.Tensor, cos: torch.Tensor, sin: torch.Tensor) -> torch.Tensor:
    """Applies Rotary Position Embedding to query/key tensors.
    x shape: (batch, n_heads, seq_len, head_dim)
    cos, sin shape: (seq_len, head_dim)
    """
    cos = cos.unsqueeze(0).unsqueeze(0)  # (1, 1, seq_len, head_dim)
    sin = sin.unsqueeze(0).unsqueeze(0)  # (1, 1, seq_len, head_dim)

    d = x.shape[-1]
    x1 = x[..., : d // 2]
    x2 = x[..., d // 2 :]
    rotated = torch.cat((-x2, x1), dim=-1)
    return (x * cos) + (rotated * sin)


class Attention(nn.Module):
    def __init__(self, c):
        super().__init__()
        self.qkv = nn.Linear(c.d_model, 3 * c.d_model)
        self.proj = nn.Linear(c.d_model, c.d_model)
        self.heads = c.n_heads
        self.head_dim = c.d_model // c.n_heads
        self.dropout = c.dropout
        self.residual_dropout = nn.Dropout(c.dropout)
        self.use_rope = getattr(c, 'use_rope', True)

        if self.use_rope:
            rope_dim = self.head_dim if self.head_dim % 2 == 0 else self.head_dim - 1
            self.rope = RotaryPositionalEncoding(
                dim=rope_dim,
                max_tokens=getattr(c, 'max_tokens', 2048),
                token_width=getattr(c, 'token_width', 1)
            )

    def forward(self, x):
        b, t, d = x.shape
        q, k, v = [z.view(b, t, self.heads, self.head_dim).transpose(1, 2) for z in self.qkv(x).chunk(3, dim=-1)]

        if self.use_rope:
            cos, sin = self.rope(t, x.device)
            rope_dim = cos.shape[-1]
            if rope_dim == self.head_dim:
                q = apply_rope(q, cos, sin)
                k = apply_rope(k, cos, sin)
            else:
                q_rope, q_pass = q[..., :rope_dim], q[..., rope_dim:]
                k_rope, k_pass = k[..., :rope_dim], k[..., rope_dim:]
                q_rope = apply_rope(q_rope, cos, sin)
                k_rope = apply_rope(k_rope, cos, sin)
                q = torch.cat((q_rope, q_pass), dim=-1)
                k = torch.cat((k_rope, k_pass), dim=-1)

        y = F.scaled_dot_product_attention(
            q, k, v, is_causal=True,
            dropout_p=self.dropout if self.training else 0.0
        )
        return self.residual_dropout(self.proj(y.transpose(1, 2).contiguous().view(b, t, d)))


class Block(nn.Module):
    def __init__(self, c):
        super().__init__()
        self.ln1 = nn.LayerNorm(c.d_model)
        self.ln2 = nn.LayerNorm(c.d_model)
        self.attn = Attention(c)
        self.mlp = nn.Sequential(
            nn.Linear(c.d_model, 4 * c.d_model),
            nn.GELU(),
            nn.Linear(4 * c.d_model, c.d_model),
            nn.Dropout(c.dropout)
        )

    def forward(self, x):
        x = x + self.attn(self.ln1(x))
        return x + self.mlp(self.ln2(x))


class Decoder(nn.Module):
    def __init__(self, c):
        super().__init__()
        self.config = c
        self.embedding = nn.Embedding(c.vocab_size, c.d_model)
        
        # Fixed positions measured in nucleotides; retained for backward compatibility.
        position = torch.arange(c.max_tokens).unsqueeze(1) * c.token_width
        frequency = torch.exp(torch.arange(0, c.d_model, 2) * (-math.log(10000) / c.d_model))
        pe = torch.zeros(c.max_tokens, c.d_model)
        pe[:, 0::2] = torch.sin(position * frequency)
        pe[:, 1::2] = torch.cos(position * frequency)
        self.register_buffer('position', pe)

        self.use_frame_encoding = getattr(c, 'use_frame_encoding', True)
        self.num_frames = getattr(c, 'num_frames', 4)
        self.frame_embedding = nn.Embedding(self.num_frames, c.d_model)
        self.frame_head = nn.Linear(c.d_model, self.num_frames)

        self.dropout = nn.Dropout(c.dropout)
        self.blocks = nn.ModuleList([Block(c) for _ in range(c.n_layers)])
        self.norm = nn.LayerNorm(c.d_model)
        self.head = nn.Linear(c.d_model, c.vocab_size, bias=False)
        self.head.weight = self.embedding.weight
        self.apply(self._init)

    def _init(self, module):
        if isinstance(module, (nn.Linear, nn.Embedding)):
            nn.init.normal_(module.weight, std=0.02)
        if isinstance(module, nn.Linear) and module.bias is not None:
            nn.init.zeros_(module.bias)

    def load_state_dict(self, state_dict: Dict[str, Any], strict: bool = True, assign: bool = False):
        if strict:
            has_fe = any('frame_embedding' in k or 'frame_head' in k for k in state_dict.keys())
            if not has_fe:
                res = super().load_state_dict(state_dict, strict=False, assign=assign)
                expected_missing = {'frame_embedding.weight', 'frame_head.weight', 'frame_head.bias'}
                actual_missing = set(res.missing_keys) - expected_missing
                if actual_missing or res.unexpected_keys:
                    raise RuntimeError(f"Error(s) in loading state_dict for Decoder: missing {actual_missing}, unexpected {res.unexpected_keys}")
                return res
        return super().load_state_dict(state_dict, strict=strict, assign=assign)

    def forward(self, x: torch.Tensor, frame_offsets: Optional[torch.Tensor] = None, return_hidden: bool = False):
        if x.shape[1] > self.config.max_tokens:
            raise ValueError('Sequence exceeds configured biological context')
        
        h = self.embedding(x)
        if hasattr(self, 'position') and self.position is not None and self.position.shape[0] >= x.shape[1]:
            h = h + self.position[:x.shape[1]]

        if self.use_frame_encoding and frame_offsets is not None and hasattr(self, 'frame_embedding'):
            f_ids = torch.where((frame_offsets >= 0) & (frame_offsets < 3), frame_offsets, self.num_frames - 1)
            h = h + self.frame_embedding(f_ids)

        x_out = self.dropout(h)
        for block in self.blocks:
            x_out = block(x_out)
        hidden = self.norm(x_out)
        logits = self.head(hidden)

        if return_hidden:
            return logits, hidden
        return logits

    def predict_frame(self, hidden: torch.Tensor) -> torch.Tensor:
        """Predicts reading frame logits from hidden representations."""
        return self.frame_head(hidden)

    def forward_multitask(
        self,
        x: torch.Tensor,
        frame_offsets: Optional[torch.Tensor] = None,
        fim_x: Optional[torch.Tensor] = None,
        span_x: Optional[torch.Tensor] = None,
    ) -> Dict[str, torch.Tensor]:
        """Runs forward passes for standard sequence and multi-task variants.

        Returns dict of logits:
            - 'logits': Autoregressive sequence logits (B, T, vocab_size)
            - 'hidden': Hidden representations (B, T, d_model)
            - 'frame_logits': Frame prediction logits (B, T, num_frames)
            - 'fim_logits': Fill-in-the-middle logits (if fim_x provided)
            - 'span_logits': Span corruption logits (if span_x provided)
        """
        logits, hidden = self.forward(x, frame_offsets=frame_offsets, return_hidden=True)
        frame_logits = self.predict_frame(hidden)

        res = {
            'logits': logits,
            'hidden': hidden,
            'frame_logits': frame_logits,
        }

        if fim_x is not None:
            res['fim_logits'] = self.forward(fim_x, frame_offsets=frame_offsets)

        if span_x is not None:
            res['span_logits'] = self.forward(span_x, frame_offsets=frame_offsets)

        return res


def summed_loss(logits, targets):
    return F.cross_entropy(logits.reshape(-1, logits.shape[-1]), targets.reshape(-1),
                           ignore_index=-100, reduction='sum')


# --- Multi-Task Objective Helpers & Loss Calculations ---

def compute_ar_loss(logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
    """Autoregressive next-token cross-entropy loss."""
    if logits.shape[1] > 1 and targets.shape[1] == logits.shape[1]:
        # Causal next-token alignment
        shift_logits = logits[:, :-1, :].contiguous()
        shift_targets = targets[:, 1:].contiguous()
        return summed_loss(shift_logits, shift_targets)
    return summed_loss(logits, targets)


def create_fim_batch(
    tokens: torch.Tensor,
    fim_rate: float = 0.5,
    psm_ratio: float = 0.5,
    pad_token_id: int = 0
) -> Tuple[torch.Tensor, torch.Tensor]:
    """Generates Fill-In-the-Middle (FIM) augmented inputs and target sequences.

    Rearranges sequence into SPM ([SUF] Suffix [PRE] Prefix [MID] Middle)
    or PSM ([PRE] Prefix [SUF] Suffix [MID] Middle).
    """
    b, t = tokens.shape
    device = tokens.device
    fim_tokens = tokens.clone()
    fim_targets = torch.full_like(tokens, -100)

    for i in range(b):
        if torch.rand(1).item() > fim_rate or t <= 4:
            continue

        # Split sequence into Prefix, Middle, Suffix
        split1 = torch.randint(1, t - 2, (1,)).item()
        split2 = torch.randint(split1 + 1, t - 1, (1,)).item()

        prefix = tokens[i, :split1]
        middle = tokens[i, split1:split2]
        suffix = tokens[i, split2:]

        is_psm = torch.rand(1).item() < psm_ratio
        if is_psm:
            # PSM: Prefix, Suffix, Middle
            new_seq = torch.cat([prefix, suffix, middle], dim=0)
            target_start = prefix.shape[0] + suffix.shape[0]
        else:
            # SPM: Suffix, Prefix, Middle
            new_seq = torch.cat([suffix, prefix, middle], dim=0)
            target_start = suffix.shape[0] + prefix.shape[0]

        fim_tokens[i, :new_seq.shape[0]] = new_seq
        fim_targets[i, target_start:new_seq.shape[0]] = middle

    return fim_tokens, fim_targets


def compute_fim_loss(logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
    """Computes Fill-In-the-Middle (FIM) loss."""
    return summed_loss(logits, targets)


def create_span_corrupted_batch(
    tokens: torch.Tensor,
    span_mask_rate: float = 0.15,
    mean_span_len: int = 3,
    mask_token_id: int = 0
) -> Tuple[torch.Tensor, torch.Tensor]:
    """Generates Span Corruption inputs and targets.

    Randomly masks contiguous spans of tokens with mask_token_id.
    Targets set masked positions to original token IDs, and unmasked positions to -100.
    """
    b, t = tokens.shape
    corrupted = tokens.clone()
    targets = torch.full_like(tokens, -100)

    for i in range(b):
        num_to_mask = int(t * span_mask_rate)
        if num_to_mask < 1:
            continue

        masked_count = 0
        while masked_count < num_to_mask:
            span_len = min(int(torch.poisson(torch.tensor([float(mean_span_len)])).item()) + 1, num_to_mask - masked_count)
            start_pos = torch.randint(0, max(1, t - span_len), (1,)).item()
            end_pos = start_pos + span_len

            corrupted[i, start_pos:end_pos] = mask_token_id
            targets[i, start_pos:end_pos] = tokens[i, start_pos:end_pos]
            masked_count += span_len

    return corrupted, targets


def compute_span_loss(logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
    """Computes Span Corruption loss over masked positions."""
    return summed_loss(logits, targets)


def compute_frame_loss(frame_logits: torch.Tensor, frame_targets: torch.Tensor, num_frames: int = 4) -> torch.Tensor:
    """Computes Frame Prediction cross-entropy loss.

    Maps unassigned/non-coding frame offset -1 to ignore index -100 or num_frames - 1.
    """
    mapped_targets = torch.where(
        (frame_targets >= 0) & (frame_targets < 3),
        frame_targets,
        torch.tensor(-100, device=frame_targets.device, dtype=frame_targets.dtype)
    )
    return F.cross_entropy(
        frame_logits.reshape(-1, num_frames),
        mapped_targets.reshape(-1),
        ignore_index=-100,
        reduction='sum'
    )


@dataclass
class MultiTaskPretrainingLoss:
    """Unified Multi-Task Pretraining Objective container."""
    w_ar: float = 1.0
    w_fim: float = 1.0
    w_span: float = 1.0
    w_frame: float = 1.0

    def __call__(
        self,
        model: Decoder,
        x: torch.Tensor,
        targets: Optional[torch.Tensor] = None,
        frame_offsets: Optional[torch.Tensor] = None,
        fim_x: Optional[torch.Tensor] = None,
        fim_targets: Optional[torch.Tensor] = None,
        span_x: Optional[torch.Tensor] = None,
        span_targets: Optional[torch.Tensor] = None,
    ) -> Dict[str, torch.Tensor]:
        """Calculates multi-task pretraining loss components and total weighted loss."""
        if targets is None:
            targets = x

        outputs = model.forward_multitask(x=x, frame_offsets=frame_offsets, fim_x=fim_x, span_x=span_x)

        ar_loss = compute_ar_loss(outputs['logits'], targets)

        fim_loss = torch.tensor(0.0, device=x.device)
        if 'fim_logits' in outputs and fim_targets is not None:
            fim_loss = compute_fim_loss(outputs['fim_logits'], fim_targets)

        span_loss = torch.tensor(0.0, device=x.device)
        if 'span_logits' in outputs and span_targets is not None:
            span_loss = compute_span_loss(outputs['span_logits'], span_targets)

        frame_loss = torch.tensor(0.0, device=x.device)
        if frame_offsets is not None:
            frame_loss = compute_frame_loss(outputs['frame_logits'], frame_offsets, num_frames=model.num_frames)

        total_loss = (
            self.w_ar * ar_loss +
            self.w_fim * fim_loss +
            self.w_span * span_loss +
            self.w_frame * frame_loss
        )

        return {
            'loss': total_loss,
            'ar_loss': ar_loss,
            'fim_loss': fim_loss,
            'span_loss': span_loss,
            'frame_loss': frame_loss,
        }

