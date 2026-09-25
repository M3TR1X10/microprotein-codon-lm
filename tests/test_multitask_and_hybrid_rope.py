"""Tests for Multi-Task Pretraining Objectives (Autoregressive, FIM, Span Corruption, Frame)
and Hybrid RoPE + Frame Positional Encodings in microprotein_lm.model.
"""
import torch
import pytest
from microprotein_lm.model import (
    ModelConfig,
    Decoder,
    RotaryPositionalEncoding,
    apply_rope,
    compute_ar_loss,
    create_fim_batch,
    compute_fim_loss,
    create_span_corrupted_batch,
    compute_span_loss,
    compute_frame_loss,
    MultiTaskPretrainingLoss,
)


def test_model_config_multitask_fields():
    """Validates multi-task and positional encoding configuration defaults and constraints."""
    cfg = ModelConfig(vocab_size=4, max_tokens=64, token_width=1)
    assert cfg.use_rope is False
    assert cfg.use_frame_encoding is True
    assert cfg.num_frames == 4
    assert cfg.fim_rate == 0.5
    assert cfg.span_mask_rate == 0.15

    cfg_rope = ModelConfig(vocab_size=4, max_tokens=64, token_width=1, use_rope=True)
    assert cfg_rope.use_rope is True

    with pytest.raises(ValueError, match="num_frames must be >= 1"):
        ModelConfig(vocab_size=4, max_tokens=64, token_width=1, num_frames=0)


def test_rope_encoding_and_application():
    """Tests Rotary Position Embedding generation and application."""
    rope = RotaryPositionalEncoding(dim=8, max_tokens=32, token_width=1)
    cos, sin = rope(16, torch.device('cpu'))

    assert cos.shape == (16, 8)
    assert sin.shape == (16, 8)

    q = torch.randn(2, 4, 16, 8)
    q_rot = apply_rope(q, cos, sin)
    assert q_rot.shape == q.shape
    # Norm per vector should be preserved under rotation
    assert torch.allclose(torch.norm(q, dim=-1), torch.norm(q_rot, dim=-1), atol=1e-5)


def test_decoder_hybrid_rope_frame_forward():
    """Tests Decoder forward pass with both RoPE and Frame Positional Encodings."""
    cfg = ModelConfig(vocab_size=64, max_tokens=30, token_width=3, d_model=16, n_heads=2, n_layers=1, use_rope=True)
    model = Decoder(cfg)

    x = torch.randint(0, 64, (2, 10))
    frame_offsets = torch.tensor([
        [0, 1, 2, 0, 1, 2, 0, 1, 2, -1],
        [-1, -1, 0, 1, 2, 0, 1, 2, 0, 1]
    ], dtype=torch.long)

    logits = model(x, frame_offsets=frame_offsets)
    assert logits.shape == (2, 10, 64)

    # Test return_hidden
    logits_h, hidden = model(x, frame_offsets=frame_offsets, return_hidden=True)
    assert logits_h.shape == (2, 10, 64)
    assert hidden.shape == (2, 10, 16)


def test_decoder_frame_prediction():
    """Tests reading frame classifier head output."""
    cfg = ModelConfig(vocab_size=4, max_tokens=40, token_width=1, d_model=16, n_heads=2, n_layers=1)
    model = Decoder(cfg)

    x = torch.randint(0, 4, (2, 12))
    _, hidden = model(x, return_hidden=True)
    frame_logits = model.predict_frame(hidden)

    assert frame_logits.shape == (2, 12, 4)


def test_autoregressive_loss():
    """Tests autoregressive next-token cross-entropy loss."""
    logits = torch.randn(2, 10, 64)
    targets = torch.randint(0, 64, (2, 10))

    loss = compute_ar_loss(logits, targets)
    assert isinstance(loss, torch.Tensor)
    assert loss.dim() == 0
    assert loss.item() > 0


def test_fill_in_the_middle_objective():
    """Tests Fill-In-the-Middle (FIM) batch creation and loss calculation."""
    torch.manual_seed(42)
    tokens = torch.randint(1, 64, (4, 20))
    fim_tokens, fim_targets = create_fim_batch(tokens, fim_rate=1.0, psm_ratio=0.5)

    assert fim_tokens.shape == tokens.shape
    assert fim_targets.shape == tokens.shape

    cfg = ModelConfig(vocab_size=64, max_tokens=30, token_width=3, d_model=16, n_heads=2, n_layers=1)
    model = Decoder(cfg)
    fim_logits = model(fim_tokens)

    loss = compute_fim_loss(fim_logits, fim_targets)
    assert isinstance(loss, torch.Tensor)
    assert loss.item() >= 0


def test_span_corruption_objective():
    """Tests Span Corruption batch creation and loss calculation."""
    torch.manual_seed(42)
    tokens = torch.randint(1, 4, (4, 20))
    corrupted, targets = create_span_corrupted_batch(tokens, span_mask_rate=0.2, mean_span_len=3, mask_token_id=0)

    assert corrupted.shape == tokens.shape
    assert targets.shape == tokens.shape

    cfg = ModelConfig(vocab_size=4, max_tokens=30, token_width=1, d_model=16, n_heads=2, n_layers=1)
    model = Decoder(cfg)
    span_logits = model(corrupted)

    loss = compute_span_loss(span_logits, targets)
    assert isinstance(loss, torch.Tensor)
    assert loss.item() >= 0


def test_frame_loss_computation():
    """Tests reading frame cross-entropy loss computation with unassigned offsets."""
    frame_logits = torch.randn(2, 6, 4)
    frame_targets = torch.tensor([
        [0, 1, 2, 0, 1, -1],
        [-1, 0, 1, 2, 0, 1]
    ], dtype=torch.long)

    loss = compute_frame_loss(frame_logits, frame_targets, num_frames=4)
    assert isinstance(loss, torch.Tensor)
    assert loss.item() > 0


def test_multitask_pretraining_loss_full_pass():
    """Tests end-to-end MultiTaskPretrainingLoss calculation and backward gradient flow."""
    cfg = ModelConfig(vocab_size=64, max_tokens=40, token_width=3, d_model=16, n_heads=2, n_layers=1)
    model = Decoder(cfg)

    x = torch.randint(1, 64, (2, 15))
    frame_offsets = torch.tensor([
        [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, -1],
        [-1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1]
    ], dtype=torch.long)

    fim_x, fim_targets = create_fim_batch(x, fim_rate=1.0)
    span_x, span_targets = create_span_corrupted_batch(x, span_mask_rate=0.2)

    loss_fn = MultiTaskPretrainingLoss(w_ar=1.0, w_fim=1.0, w_span=1.0, w_frame=0.5)

    loss_dict = loss_fn(
        model=model,
        x=x,
        targets=x,
        frame_offsets=frame_offsets,
        fim_x=fim_x,
        fim_targets=fim_targets,
        span_x=span_x,
        span_targets=span_targets,
    )

    assert 'loss' in loss_dict
    assert 'ar_loss' in loss_dict
    assert 'fim_loss' in loss_dict
    assert 'span_loss' in loss_dict
    assert 'frame_loss' in loss_dict

    total_loss = loss_dict['loss']
    assert total_loss.requires_grad
    total_loss.backward()

    # Check gradients exist on embeddings and heads
    assert model.embedding.weight.grad is not None
    assert model.frame_head.weight.grad is not None


def test_legacy_checkpoint_strict_loading():
    """Ensures legacy state_dict without frame_embedding/frame_head loads cleanly with strict=True."""
    cfg = ModelConfig(vocab_size=4, max_tokens=20, token_width=1, d_model=16, n_heads=2, n_layers=1)
    model = Decoder(cfg)

    state = model.state_dict()
    # Simulate legacy checkpoint by removing new keys
    legacy_state = {k: v for k, v in state.items() if not k.startswith('frame_embedding') and not k.startswith('frame_head')}

    new_model = Decoder(cfg)
    new_model.load_state_dict(legacy_state, strict=True)
