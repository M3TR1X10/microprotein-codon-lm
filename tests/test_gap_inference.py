"""Inference correctness fixtures; these invented strings are never biological data."""
import math

import pytest
import torch

from microprotein_lm.gap_inference import GapDecoder
from microprotein_lm.model import Decoder, ModelConfig


@pytest.fixture(autouse=True)
def single_thread():
    old = torch.get_num_threads()
    torch.set_num_threads(1)
    yield
    torch.set_num_threads(old)


def decoder(mode='base', max_nt=60):
    torch.manual_seed(117)
    width = 1 if mode == 'base' else 3
    model = Decoder(ModelConfig(4 if width == 1 else 64, max_nt // width - 1,
                                width, d_model=16, n_heads=2, n_layers=2, dropout=.2))
    return GapDecoder(model, mode)


@pytest.mark.parametrize('mode', ['base', 'codon'])
def test_cache_matches_original_model_and_chunked_causality(mode):
    dec = decoder(mode)
    ids = dec.tokenizer.encode('AUGAAUCCGGAUUAC')
    tokens = torch.tensor([ids])
    with torch.inference_mode():
        expected = dec.model(tokens)
        _, actual = dec.advance(tokens)
        state, first = dec.advance(tokens[:, :2])
        state, second = dec.advance(tokens[:, 2:4], state)
        _, third = dec.advance(tokens[:, 4:], state)
    torch.testing.assert_close(actual, expected, rtol=2e-5, atol=2e-6)
    torch.testing.assert_close(torch.cat([first, second, third], dim=1), expected,
                               rtol=2e-5, atol=2e-6)


@pytest.mark.parametrize('mode', ['base', 'codon'])
def test_offset_uses_original_absolute_coordinates(mode):
    dec = decoder(mode)
    tokens = torch.tensor([dec.tokenizer.encode('AAUCCG')])
    offset_nt = 9
    offset = offset_nt // dec.tokenizer.width
    with torch.inference_mode():
        x = dec.model.embedding(tokens) + dec.model.position[offset:offset + tokens.shape[1]]
        for block in dec.model.blocks:
            x = block(x)
        expected = dec.model.head(dec.model.norm(x))
        _, actual = dec.advance(tokens, offset_nt=offset_nt)
    torch.testing.assert_close(actual, expected, rtol=2e-5, atol=2e-6)
    dec.model.position.zero_()
    a, _ = dec.advance(tokens, offset_nt=0)
    b, _ = dec.advance(tokens, offset_nt=offset_nt)
    torch.testing.assert_close(a.logits, b.logits, rtol=0, atol=0)


@pytest.mark.parametrize('mode', ['base', 'codon'])
def test_exact_codon_distribution_matches_exhaustive_teacher_forcing(mode):
    dec = decoder(mode)
    logp, branches = dec.next_codons(dec.prefix('AUG'))
    expected = torch.tensor([dec.continuation_log_probability('AUG', c)
                             for c in dec.codons], dtype=torch.float64)
    torch.testing.assert_close(logp[0], expected, rtol=1e-5, atol=1e-6)
    assert float(logp.exp().sum()) == pytest.approx(1, abs=1e-12)
    selected = [0, 7, 16, 63]
    state = dec.append_codons(branches, [0] * len(selected), selected)
    for i, codon in enumerate(selected):
        direct = dec.prefix('AUG' + dec.codons[codon])
        torch.testing.assert_close(state.logits[i], direct.logits[0], rtol=2e-5, atol=2e-6)


@pytest.mark.parametrize('mode', ['base', 'codon'])
def test_greedy_and_beam_use_true_joint_codon_scores(mode):
    dec = decoder(mode)
    greedy = dec.proposals('AUG', 2, width=1)
    assert len(greedy) == 1 and len(greedy[0]['rna']) == 6
    expected = dec.continuation_log_probability('AUG', greedy[0]['rna'])
    assert greedy[0]['log_probability'] == pytest.approx(expected, abs=2e-6)
    all_codons = dec.proposals('AUG', 1, width=64)
    assert {c['rna'] for c in all_codons} == set(dec.codons)
    assert all_codons[0]['rna'] == dec.proposals('AUG', 1)[0]['rna']
    assert sum(math.exp(c['log_probability']) for c in all_codons) == pytest.approx(1, abs=1e-12)
    assert {'UAA', 'UAG', 'UGA', 'AGA', 'AGG'} <= {c['rna'] for c in all_codons}
    beam = dec.proposals('AUG', 3, width=8)
    assert len({c['rna'] for c in beam}) == 8
    assert all(len(c['rna']) == 9 for c in beam)
    for candidate in beam:
        assert candidate['log_probability'] == pytest.approx(
            dec.continuation_log_probability('AUG', candidate['rna']), abs=3e-6)


@pytest.mark.parametrize('mode', ['base', 'codon'])
def test_reranking_is_fixed_bank_joint_likelihood_without_truth(mode):
    dec = decoder(mode)
    candidates = dec.proposals('AUG', 2, width=8)
    original = [dict(c) for c in candidates]
    no_suffix = dec.rerank('AUG', candidates, '')
    assert [c['rna'] for c in no_suffix] == [c['rna'] for c in candidates]
    ranked = dec.rerank('AUG', candidates, 'CCGAUU')
    assert candidates == original
    assert {c['rna'] for c in ranked} == {c['rna'] for c in candidates}
    for c in ranked:
        assert c['joint_log_score'] == pytest.approx(
            dec.continuation_log_probability('AUG', c['rna'] + 'CCGAUU'), abs=4e-6)
    assert [c['joint_log_score'] for c in ranked] == sorted(
        [c['joint_log_score'] for c in ranked], reverse=True)


@pytest.mark.parametrize('mode', ['base', 'codon'])
def test_last_output_can_reach_capacity_without_consuming_it(mode):
    dec = decoder(mode, max_nt=12)
    assert len(dec.proposals('AUG', 3)[0]['rna']) == 9
    with pytest.raises(ValueError, match='capacity'):
        dec.proposals('AUG', 4)
    with pytest.raises(ValueError, match='boundary'):
        dec.prefix('AUG', offset_nt=1)
    with pytest.raises(ValueError, match='whole codons'):
        dec.prefix('AU')


def test_zero_model_ties_are_deterministic_and_stops_not_suppressed():
    dec = decoder('codon')
    with torch.no_grad():
        for parameter in dec.model.parameters():
            parameter.zero_()
    result = dec.proposals('AUG', 2, width=8)
    assert [c['rna'] for c in result] == ['AAA' + c for c in dec.codons[:8]]
    assert all(c['log_probability'] == pytest.approx(-2 * math.log(64)) for c in result)


@pytest.mark.parametrize('mode', ['base', 'codon'])
def test_generation_does_not_mutate_any_checkpoint_tensor(mode):
    dec = decoder(mode)
    before = {k: v.clone() for k, v in dec.model.state_dict().items()}
    candidates = dec.proposals('AUG', 2, width=8)
    dec.rerank('AUG', candidates, 'GAA')
    assert not dec.model.training
    for key, value in dec.model.state_dict().items():
        torch.testing.assert_close(before[key], value, rtol=0, atol=0)
    assert all(p.grad is None for p in dec.model.parameters())


@pytest.mark.parametrize('mode', ['base', 'codon'])
def test_altering_future_tokens_cannot_change_earlier_logits(mode):
    dec = decoder(mode)
    a = dec.tokenizer.encode('AUGAAACCC')
    b = dec.tokenizer.encode('AUGAAAGGG')
    _, la = dec.advance([a])
    _, lb = dec.advance([b])
    cutoff = 6 // dec.tokenizer.width
    torch.testing.assert_close(la[:, :cutoff], lb[:, :cutoff], rtol=0, atol=0)


def test_full_path_tie_break_does_not_depend_on_parent_score_order(monkeypatch):
    dec = decoder('codon')
    calls = 0
    def distributions(state):
        nonlocal calls
        calls += 1
        values = torch.full((state.logits.shape[0], 64), -100., dtype=torch.float64)
        if calls == 1:
            values[0, 0], values[0, 1] = -2., -1.  # AAU parent ranks first
        else:
            values[0, 0], values[1, 0] = -2., -1.  # totals tie at -3
        return values, state
    monkeypatch.setattr(dec, 'next_codons', distributions)
    assert [c['rna'] for c in dec.proposals('AUG', 2, width=2)] == ['AAAAAA', 'AAUAAA']
