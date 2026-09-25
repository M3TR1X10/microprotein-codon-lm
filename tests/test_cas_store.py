import hashlib
import json
import pytest
from microprotein_lm.cas_store import CASStore, ObservationLedger, CanonicalSequences


def test_cas_store(tmp_path):
    store_dir = tmp_path / "raw_cas"
    cas = CASStore(store_dir)

    payload = "<?xml version='1.0'?><INSDSeq>TEST_PAYLOAD</INSDSeq>"
    sha256_hash = cas.store(payload)

    expected_hash = hashlib.sha256(payload.encode('utf-8')).hexdigest()
    assert sha256_hash == expected_hash
    assert cas.exists(sha256_hash) is True
    assert cas.retrieve_text(sha256_hash) == payload

    # Re-storing should be idempotent
    assert cas.store(payload) == sha256_hash


def test_observation_ledger(tmp_path):
    ledger_path = tmp_path / "raw_cas" / "observation_ledger.jsonl"
    ledger = ObservationLedger(ledger_path)

    obs1 = ledger.record_observation("NP_001.1", "https://ebi.ac.uk/api/xml/NP_001.1", "abc123sha256", {"gene": "GENE1"})
    obs2 = ledger.record_observation("NP_002.1", "https://ebi.ac.uk/api/xml/NP_002.1", "def456sha256", {"gene": "GENE2"})

    all_obs = ledger.get_observations()
    assert len(all_obs) == 2

    filtered = ledger.get_observations(accession="NP_001.1")
    assert len(filtered) == 1
    assert filtered[0]['accession'] == "NP_001.1"
    assert filtered[0]['raw_sha256'] == "abc123sha256"


def test_canonical_sequences(tmp_path):
    canon_path = tmp_path / "processed" / "canonical_sequences.jsonl"
    canon = CanonicalSequences(canon_path)

    seq_rna = "AUGGCUUAA"
    seq_sha256 = hashlib.sha256(seq_rna.encode('utf-8')).hexdigest()

    rec1 = canon.register_sequence(
        sequence_sha256=seq_sha256,
        rna=seq_rna,
        protein="MA",
        initiation_codon="AUG",
        stop_codon="UAA",
        metadata={"gene": "TEST1"},
        source_accessions=["ACC_1"],
        observation_ids=["obs_1"]
    )

    assert rec1['sequence_sha256'] == seq_sha256
    assert rec1['source_accessions'] == ["ACC_1"]

    # Upsert with second accession sharing the exact same sequence
    rec2 = canon.register_sequence(
        sequence_sha256=seq_sha256,
        rna=seq_rna,
        protein="MA",
        initiation_codon="AUG",
        stop_codon="UAA",
        metadata={"gene": "TEST1"},
        source_accessions=["ACC_2"],
        observation_ids=["obs_2"]
    )

    assert rec2['sequence_sha256'] == seq_sha256
    assert rec2['source_accessions'] == ["ACC_1", "ACC_2"]
    assert rec2['observation_ids'] == ["obs_1", "obs_2"]

    retrieved = canon.get_sequence(seq_sha256)
    assert retrieved is not None
    assert retrieved['source_accessions'] == ["ACC_1", "ACC_2"]
