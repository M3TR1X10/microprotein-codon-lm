import pytest
from microprotein_lm.gap_metrics import reconstruction, TrainingPositionSupport, translate_gap


def test_translation_uses_correct_internal_genetic_code():
    assert translate_gap('UGAAGAAGGAUA', 1) == '*RRI'
    assert translate_gap('UGAAGAAGGAUA', 2) == 'W**M'
    assert translate_gap('GUG', 2) == 'V'  # internal alternative start is not methionine


def test_reconstruction_keeps_exact_span_and_synonymous_error_distinct():
    result = reconstruction('AAUGAAUGA', 'AACGAAAGA', 2)
    assert result['exact_span'] == 0
    assert result['correct_codons'] == 1
    assert result['correct_amino_acids'] == 2
    assert result['synonymous_errors'] == 1
    assert result['has_premature_stop'] == 1
    assert result['first_error_codon'] == 1
    assert result['codon_matches'] == [0, 1, 0]
    assert result['amino_acid_matches'] == [1, 1, 0]
    with pytest.raises(ValueError):
        reconstruction('AAA', 'AAAAAA', 1)


def test_minority_support_is_not_saturated_variable_position_mask():
    rows = [{'family': 'X', 'tax_id': 1, 'rna': r} for r in
            ['AUGAAAUAA', 'AUGAAAUAA', 'AUGAACUAA']]
    support = TrainingPositionSupport(rows)
    minority = support.describe(rows[0], 1, 'AACUAG', 'AACUAA')
    assert minority['codon_support_labels'] == ['observed_minority', 'unobserved_codon']
    assert minority['training_support']['observed_minority']['correct_codons'] == 1
    other = support.describe({'family': 'X', 'tax_id': 2}, 1, 'AAA', 'AAA')
    assert other['family_seen'] and not other['taxon_seen']
    assert other['codon_support_labels'] == ['unsupported_position']


def test_tied_training_modes_are_all_majority():
    support = TrainingPositionSupport([{'family': 'X', 'tax_id': 1, 'rna': r}
                                        for r in ['AUGAAA', 'AUGAAC']])
    for truth in ['AAA', 'AAC']:
        result = support.describe({'family': 'X', 'tax_id': 1}, 1, truth, truth)
        assert result['codon_support_labels'] == ['majority']
