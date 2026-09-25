import pytest
from microprotein_lm.tokenization import (
    BioTokenizer,
    Tokenizer,
    PAD_TOKEN,
    UNK_TOKEN,
    BOS_TOKEN,
    EOS_TOKEN,
    MASK_TOKEN,
    FRAME_0_TOKEN,
    FRAME_1_TOKEN,
    FRAME_2_TOKEN,
    CONTROL_TOKENS,
    IUPAC_DEGENERATE_BASES,
)


def test_biotokenizer_control_tokens_and_vocab_sizes():
    bt_base = BioTokenizer(mode="base", include_special_tokens=True, include_degenerate=True)
    # 8 control + 4 standard + 3 degenerate = 15
    assert len(bt_base.vocabulary) == 15
    for tok in CONTROL_TOKENS:
        assert tok in bt_base.vocabulary
        assert bt_base.is_control_token(tok)
        assert bt_base.is_control_token(bt_base.token_to_id[tok])

    assert bt_base.pad_id == 0
    assert bt_base.unk_id == 1
    assert bt_base.bos_id == 2
    assert bt_base.eos_id == 3
    assert bt_base.mask_id == 4
    assert bt_base.frame_0_id == 5
    assert bt_base.frame_1_id == 6
    assert bt_base.frame_2_id == 7


def test_biotokenizer_iupac_degenerate_bases_base_mode():
    bt_base = BioTokenizer(mode="base")
    for base in ["N", "R", "Y"]:
        assert base in bt_base.vocabulary
        assert bt_base.is_degenerate(base)
        assert bt_base.is_degenerate(bt_base.token_to_id[base])

    # Encode sequence with degenerate bases and DNA T
    encoded = bt_base.encode("ATGNRY")  # T -> U
    decoded = bt_base.decode(encoded)
    assert decoded == "AUGNRY"


def test_biotokenizer_iupac_degenerate_bases_codon_mode():
    bt_codon = BioTokenizer(mode="codon")
    # 8 control + 64 standard + 279 degenerate triplets = 351
    assert len(bt_codon.vocabulary) == 351
    assert "AUG" in bt_codon.vocabulary
    assert "NNN" in bt_codon.vocabulary
    assert "UAR" in bt_codon.vocabulary

    assert bt_codon.is_degenerate("NNN")
    assert bt_codon.is_degenerate("UAR")
    assert not bt_codon.is_degenerate("AUG")

    encoded = bt_codon.encode("AUGNNNUAR")
    decoded = bt_codon.decode(encoded)
    assert decoded == "AUGNNNUAR"


def test_biotokenizer_special_token_encoding_options():
    bt_codon = BioTokenizer(mode="codon")
    encoded = bt_codon.encode("AUGUUUUAA", add_bos=True, add_eos=True, frame=1)
    decoded_full = bt_codon.decode(encoded)
    decoded_skipped = bt_codon.decode(encoded, skip_special_tokens=True)

    assert decoded_full == "[FRAME_1][BOS]AUGUUUUAA[EOS]"
    assert decoded_skipped == "AUGUUUUAA"


def test_biotokenizer_expand_degenerate():
    bt_base = BioTokenizer(mode="base")
    assert sorted(bt_base.expand_degenerate("N")) == ["A", "C", "G", "U"]
    assert sorted(bt_base.expand_degenerate("R")) == ["A", "G"]
    assert sorted(bt_base.expand_degenerate("Y")) == ["C", "U"]

    bt_codon = BioTokenizer(mode="codon")
    assert sorted(bt_codon.expand_degenerate("UAR")) == ["UAA", "UAG"]
    assert sorted(bt_codon.expand_degenerate("CAY")) == ["CAC", "CAU"]
    assert len(bt_codon.expand_degenerate("NNN")) == 64


def test_legacy_tokenizer_backward_compatibility():
    t_base = Tokenizer("base")
    assert len(t_base.vocabulary) == 4
    assert t_base.vocabulary == ["A", "U", "C", "G"]

    encoded = t_base.encode("AUGUUUUAA")
    assert t_base.decode(encoded) == "AUGUUUUAA"

    for invalid in ["", "ATG", "ANN", "aug"]:
        with pytest.raises(ValueError):
            t_base.encode(invalid)

    t_codon = Tokenizer("codon")
    assert len(t_codon.vocabulary) == 64
    with pytest.raises(ValueError):
        t_codon.encode("AUGU")
