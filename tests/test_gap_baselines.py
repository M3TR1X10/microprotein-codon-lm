import math

import pytest

from microprotein_lm.gap_baselines import GapBaselines, METHODS


def baseline():
    return GapBaselines([{'family': 'FIXTURE', 'tax_id': 1, 'rna': seq}
                         for seq in ('AUGAAACCCUAA', 'AUGAAACCGUAA')])


def test_uniform_is_exactly_two_bits_per_base_and_frozen_tie_order():
    b = baseline()
    result = b.score_and_generate('uniform', 'AUG', 'AAACCC', 1, {'family': 'X', 'tax_id': 2})
    assert result['gap_bits_per_base'] == pytest.approx(2)
    assert result['prediction'] == 'AAAAAA'
    assert result['gap_log_probability'] == pytest.approx(-2 * math.log(64))


def test_baselines_normalize_and_use_total_prior_two():
    b = baseline()
    for method in METHODS:
        p, _ = b.distribution(method, 'AUG', 1, {'family': 'FIXTURE', 'tax_id': 1})
        assert p.sum() == pytest.approx(1)
        assert p.min() > 0
    p, _ = b.distribution('position', 'AUG', 1, {})
    assert p[b.index['AAA']] == pytest.approx((2 + 2 / 64) / 4)


def test_missing_metadata_falls_back_only_to_training_counts():
    b = baseline()
    _, fallback = b.distribution('privileged_family_taxon_position', 'AUG', 1,
                                  {'family': 'UNSEEN', 'tax_id': 99})
    assert fallback == 'pooled_position'
    _, fallback = b.distribution('privileged_family_taxon_position', 'AUG', 50,
                                  {'family': 'UNSEEN', 'tax_id': 99})
    assert fallback == 'pooled_unigram'
    _, fallback = b.distribution('bigram', 'GGG', 1, {})
    assert fallback == 'pooled_unigram'


def test_greedy_path_cannot_depend_on_hidden_truth():
    b = baseline()
    record = {'family': 'FIXTURE', 'tax_id': 1}
    for method in METHODS:
        a = b.score_and_generate(method, 'AUG', 'AAACCC', 1, record)
        z = b.score_and_generate(method, 'AUG', 'GGGAAA', 1, record)
        assert a['prediction'] == z['prediction']
    a = b.score_and_generate('bigram', 'AUG', 'AAACCC', 1, record)
    assert a['scoring_fallbacks'] == {}
