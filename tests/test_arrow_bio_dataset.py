import hashlib
import json
import numpy as np
import pytest
import torch
import pyarrow as pa
import pyarrow.ipc as ipc
import pyarrow.parquet as pq

from microprotein_lm.data import ArrowBioDataset, Corpus, prepare, compute_directory_hash
from microprotein_lm.io import write_json, sha256
from microprotein_lm.tokenization import Tokenizer


@pytest.fixture
def synthetic_records():
    rnas = ['AUGUGAUAA', 'AUGCCCUAG', 'AUGAAACCCUAA']
    records = []
    for rna in rnas:
        seq_hash = hashlib.sha256(rna.encode()).hexdigest()
        records.append({
            'sequence_sha256': seq_hash,
            'rna': rna,
            'protein': 'TEST',
            'tax_id': 9606,
            'family': 'ATP8',
            'translation_table': 2,
        })
    return records


@pytest.fixture
def prepared_root(tmp_path, synthetic_records):
    root = tmp_path / 'data'
    root.mkdir()
    cohort_path = root / 'cohort.jsonl'
    cohort_path.write_text(''.join(json.dumps(r) + '\n' for r in synthetic_records), encoding='utf-8')
    write_json(root / 'manifest.json', {'cohort_sha256': sha256(cohort_path)})
    prepare(root, shard_size=2)
    return root


def test_arrow_bio_dataset_from_records_and_shards(tmp_path, synthetic_records):
    dataset_dir = tmp_path / 'arrow_ds'
    ds = ArrowBioDataset.from_records(synthetic_records, dataset_dir, mode='codon', shard_size=2)
    
    assert len(ds) == len(synthetic_records)
    assert (dataset_dir / 'tokens.arrow').exists()
    assert (dataset_dir / 'shards').exists()
    
    shards = list((dataset_dir / 'shards').glob('*.parquet'))
    assert len(shards) == 2  # 3 records with shard_size=2 -> 2 shards
    
    # Verify zero-copy sequence extraction
    tokenizer = Tokenizer('codon')
    for i, rec in enumerate(synthetic_records):
        expected_tokens = np.asarray(tokenizer.encode(rec['rna']), dtype=np.int64)
        np.testing.assert_array_equal(ds.sequence(i), expected_tokens)
        
        record_meta = ds.get_record(i)
        assert record_meta['sequence_sha256'] == rec['sequence_sha256']
        assert record_meta['rna'] == rec['rna']
        assert record_meta['family'] == rec['family']


def test_arrow_bio_dataset_prepare_and_corpus_compatibility(prepared_root):
    base_ds = Corpus(prepared_root, 'base')
    codon_ds = Corpus(prepared_root, 'codon')
    
    assert len(base_ds) == 3
    assert len(codon_ds) == 3
    
    # Check Arrow table and Parquet dataset availability
    assert base_ds.arrow_table is not None
    assert codon_ds.arrow_table is not None
    assert base_ds.parquet_dataset is not None
    assert codon_ds.parquet_dataset is not None

    # Test batch generation
    x_b, y_b = base_ds.batch([0, 1, 2])
    x_c, y_c = codon_ds.batch([0, 1, 2])
    
    assert x_b.shape[0] == 3
    assert x_c.shape[0] == 3
    assert (y_b[:, :2] == -100).all()  # Base mode masks first 2 target tokens


def test_arrow_bio_dataset_checksum_tampering(prepared_root):
    arrow_file = prepared_root / 'base' / 'tokens.arrow'
    arrow_file.write_bytes(arrow_file.read_bytes() + b'\x00')
    
    with pytest.raises(ValueError, match='checksum'):
        ArrowBioDataset(prepared_root / 'base', mode='base')

    shards_dir = prepared_root / 'codon' / 'shards'
    shard_files = list(shards_dir.glob('*.parquet'))
    shard_files[0].write_bytes(shard_files[0].read_bytes() + b'\x00')
    
    with pytest.raises(ValueError, match='checksum'):
        ArrowBioDataset(prepared_root / 'codon', mode='codon')


def test_standalone_arrow_ipc_and_parquet_utilities(tmp_path):
    schema = pa.schema([
        ('sequence_sha256', pa.string()),
        ('rna', pa.string()),
        ('length', pa.int32())
    ])
    table = pa.Table.from_pydict({
        'sequence_sha256': ['a', 'b', 'c'],
        'rna': ['AUGUGA', 'AUGUAA', 'AUGUAG'],
        'length': [6, 6, 6]
    }, schema=schema)
    
    arrow_path = tmp_path / 'test.arrow'
    hash_arrow = ArrowBioDataset.to_arrow_ipc(table, arrow_path)
    assert arrow_path.exists()
    assert hash_arrow == sha256(arrow_path)
    
    shards_dir = tmp_path / 'shards'
    hash_shards = ArrowBioDataset.to_parquet_shards(table, shards_dir, shard_size=2)
    assert shards_dir.exists()
    assert hash_shards == compute_directory_hash(shards_dir)
