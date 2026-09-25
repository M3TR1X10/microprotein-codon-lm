import pytest
from microprotein_lm.extended_window import WindowConfig, ExtendedWindowExtractor


def test_window_config_defaults():
    cfg = WindowConfig()
    assert cfg.utr5_len == 50
    assert cfg.kozak_upstream == 6
    assert cfg.kozak_downstream == 4
    assert cfg.utr3_len == 30


def test_extract_window_basic():
    extractor = ExtendedWindowExtractor(WindowConfig(utr5_len=10, kozak_upstream=6, kozak_downstream=4, utr3_len=10))
    # Full transcript: 10 nt 5' UTR + AUG GCA UAA (9 nt CDS) + 10 nt 3' UTR = 29 nt
    utr5 = "C" * 10
    cds = "AUGGCAUAA"
    utr3 = "G" * 10
    transcript = utr5 + cds + utr3

    res = extractor.extract_window(transcript, tis_start=10, cds_len=9)
    assert res['utr5'] == utr5
    assert res['cds'] == cds
    assert res['utr3'] == utr3
    assert res['extended_seq'] == transcript
    assert res['tis_offset_in_extended'] == 10
    assert len(res['segment_ids']) == len(transcript)
    # Check Kozak window tagging (0 = 5' UTR, 1 = Kozak, 2 = CDS, 3 = 3' UTR)
    # tis_start is 10. Kozak range in extended_seq is [10-6 : 10+4] = [4:14]
    assert res['segment_ids'][0] == 0  # 5' UTR outside Kozak
    assert res['segment_ids'][4] == 1  # Start of Kozak in 5' UTR
    assert res['segment_ids'][9] == 1  # Last Kozak position in 5' UTR
    assert res['segment_ids'][10] == 1 # TIS start (Kozak)
    assert res['segment_ids'][13] == 1 # Kozak in CDS
    assert res['segment_ids'][14] == 2 # CDS outside Kozak
    assert res['segment_ids'][-1] == 3 # 3' UTR


def test_extract_window_short_transcript_padding():
    extractor = ExtendedWindowExtractor(WindowConfig(utr5_len=50, utr3_len=30, pad_char='A'))
    # Short transcript with only 5 nt 5' UTR, 6 nt CDS, 5 nt 3' UTR
    transcript = "GGGGG" + "AUGCUA" + "CCCCC"
    res = extractor.extract_window(transcript, tis_start=5, cds_len=6)
    assert len(res['utr5']) == 50
    assert res['utr5'].startswith("A" * 45)
    assert len(res['utr3']) == 30
    assert res['utr3'].endswith("A" * 25)
    assert len(res['extended_seq']) == 50 + 6 + 30
