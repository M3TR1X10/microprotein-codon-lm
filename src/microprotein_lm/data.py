import json
from pathlib import Path
from typing import List, Dict, Any, Union, Optional

import numpy as np
import pyarrow as pa
import pyarrow.ipc as ipc
import pyarrow.parquet as pq
import pyarrow.dataset as ds
import torch

from .io import read_json, write_json, sha256
from .tokenization import Tokenizer


def compute_directory_hash(directory_path: Path) -> str:
    """Compute deterministic SHA-256 over files in a directory."""
    import hashlib
    hasher = hashlib.sha256()
    for file_path in sorted(Path(directory_path).glob('**/*')):
        if file_path.is_file():
            hasher.update(file_path.name.encode('utf-8'))
            hasher.update(file_path.read_bytes())
    return hasher.hexdigest()


class ArrowBioDataset:
    """Apache Arrow & Parquet Sharded Storage Dataset with zero-copy memory mapping."""

    def __init__(self, root: Union[str, Path], mode: str = 'codon'):
        self.mode = mode
        self.root = Path(root)
        self.path = self.root if (self.root / 'metadata.json').exists() else self.root / mode
        if not (self.path / 'metadata.json').exists():
            raise FileNotFoundError(f"Metadata not found in {self.path}")

        self.metadata = read_json(self.path / 'metadata.json')
        self.tokenizer = Tokenizer(mode)

        if self.metadata.get('split') != 'train' or self.metadata.get('vocabulary') != self.tokenizer.vocabulary:
            raise ValueError('Unexpected split or vocabulary')

        # Validate checksums for integrity
        for name, key in [('tokens.bin', 'tokens_sha256'), ('offsets.npy', 'offsets_sha256')]:
            if (self.path / name).exists() and key in self.metadata:
                if sha256(self.path / name) != self.metadata[key]:
                    raise ValueError('Encoded data checksum mismatch')

        if (self.path / 'tokens.arrow').exists() and 'arrow_sha256' in self.metadata:
            if sha256(self.path / 'tokens.arrow') != self.metadata['arrow_sha256']:
                raise ValueError('Encoded data checksum mismatch')

        if (self.path / 'shards').exists() and 'parquet_shards_sha256' in self.metadata:
            if compute_directory_hash(self.path / 'shards') != self.metadata['parquet_shards_sha256']:
                raise ValueError('Encoded data checksum mismatch')

        # Load binary tokens & offsets if present
        if (self.path / 'tokens.bin').exists() and (self.path / 'offsets.npy').exists():
            self.tokens = np.memmap(self.path / 'tokens.bin', dtype=np.uint8, mode='r')
            self.offsets = np.load(self.path / 'offsets.npy', allow_pickle=False)
            if self.offsets[0] != 0 or self.offsets[-1] != len(self.tokens) or np.any(np.diff(self.offsets) <= 1):
                raise ValueError('Invalid sequence boundaries')
        else:
            self.tokens = None
            self.offsets = None

        # Memory-map Arrow IPC Table if available
        self.arrow_mmap = None
        self.arrow_table = None
        if (self.path / 'tokens.arrow').exists():
            self.arrow_mmap = pa.memory_map(str(self.path / 'tokens.arrow'), 'r')
            reader = ipc.RecordBatchFileReader(self.arrow_mmap)
            self.arrow_table = reader.read_all()

        # Load Parquet Dataset if available
        self.parquet_dataset = None
        if (self.path / 'shards').exists():
            self.parquet_dataset = ds.dataset(str(self.path / 'shards'), format='parquet')

    def __len__(self) -> int:
        if self.offsets is not None:
            return len(self.offsets) - 1
        if self.arrow_table is not None:
            return self.arrow_table.num_rows
        if self.parquet_dataset is not None:
            return self.parquet_dataset.to_table().num_rows
        return 0

    def sequence(self, index: int) -> np.ndarray:
        idx = int(index)
        if self.arrow_table is not None and 'tokens' in self.arrow_table.column_names:
            tokens_val = self.arrow_table.column('tokens')[idx]
            return np.asarray(tokens_val.values.to_numpy(), dtype=np.int64)
        if self.tokens is not None and self.offsets is not None:
            return self.tokens[self.offsets[idx]:self.offsets[idx + 1]].astype(np.int64)
        raise RuntimeError("No token data source available in dataset")

    def get_record(self, index: int) -> Dict[str, Any]:
        idx = int(index)
        if self.arrow_table is not None:
            return {col: self.arrow_table.column(col)[idx].as_py() for col in self.arrow_table.column_names}
        if self.parquet_dataset is not None:
            table = self.parquet_dataset.to_table()
            return {col: table.column(col)[idx].as_py() for col in table.column_names}
        raise RuntimeError("No metadata table available in dataset")

    def batch(self, indices: List[int], device: str = 'cpu'):
        sequences = [self.sequence(int(i)) for i in indices]
        length = max(len(s) for s in sequences) - 1
        x = torch.zeros((len(indices), length), dtype=torch.long)
        y = torch.full_like(x, -100)
        for i, seq in enumerate(sequences):
            x[i, :len(seq) - 1] = torch.from_numpy(seq[:-1])
            y[i, :len(seq) - 1] = torch.from_numpy(seq[1:])
            if self.mode == 'base':
                y[i, :2] = -100  # both models condition on the same first codon
        return x.to(device), y.to(device)

    @classmethod
    def to_arrow_ipc(cls, table: pa.Table, output_path: Union[str, Path]) -> str:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with pa.OSFile(str(output_path), 'wb') as sink:
            with ipc.new_file(sink, table.schema) as writer:
                writer.write_table(table)
        return sha256(output_path)

    @classmethod
    def to_parquet_shards(cls, table: pa.Table, output_dir: Union[str, Path], shard_size: int = 1000) -> str:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        n_rows = len(table)
        for start in range(0, max(1, n_rows), shard_size):
            end = min(start + shard_size, n_rows)
            shard_table = table.slice(start, end - start)
            shard_path = output_dir / f"shard_{start // shard_size:04d}.parquet"
            pq.write_table(shard_table, shard_path, compression='SNAPPY')
        return compute_directory_hash(output_dir)

    @classmethod
    def from_records(cls, records: List[Dict[str, Any]], output_dir: Union[str, Path], mode: str = 'codon', shard_size: int = 1000):
        tokenizer = Tokenizer(mode)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        seq_ids = [r['sequence_sha256'] for r in records]
        rnas = [r['rna'] for r in records]
        encoded = [tokenizer.encode(rna) for rna in rnas]
        lengths = [len(x) for x in encoded]

        schema_fields = [
            ('sequence_sha256', pa.string()),
            ('rna', pa.string()),
            ('tokens', pa.list_(pa.uint8())),
            ('length', pa.int32()),
        ]
        extra_keys = [k for k in records[0].keys() if k not in ('sequence_sha256', 'rna')] if records else []
        for k in extra_keys:
            val = records[0][k]
            if isinstance(val, int):
                schema_fields.append((k, pa.int64()))
            elif isinstance(val, float):
                schema_fields.append((k, pa.float64()))
            else:
                schema_fields.append((k, pa.string()))

        schema = pa.schema(schema_fields)

        columns = {
            'sequence_sha256': seq_ids,
            'rna': rnas,
            'tokens': encoded,
            'length': lengths,
        }
        for k in extra_keys:
            columns[k] = [str(r[k]) if not isinstance(r[k], (int, float)) else r[k] for r in records]

        table = pa.Table.from_pydict(columns, schema=schema)

        shards_dir = output_dir / 'shards'
        shards_hash = cls.to_parquet_shards(table, shards_dir, shard_size=shard_size)

        arrow_path = output_dir / 'tokens.arrow'
        arrow_hash = cls.to_arrow_ipc(table, arrow_path)

        arrays = [np.asarray(x, dtype=np.uint8) for x in encoded]
        tokens_bin = output_dir / 'tokens.bin'
        offsets_npy = output_dir / 'offsets.npy'
        np.concatenate(arrays).tofile(tokens_bin)
        np.save(offsets_npy, np.asarray([0, *np.cumsum(lengths)], dtype=np.int64))

        metadata = {
            'mode': mode,
            'vocabulary': tokenizer.vocabulary,
            'vocab_size': len(tokenizer.vocabulary),
            'sequences': len(records),
            'split': 'train',
            'tokens_sha256': sha256(tokens_bin),
            'offsets_sha256': sha256(offsets_npy),
            'arrow_sha256': arrow_hash,
            'parquet_shards_sha256': shards_hash,
            'tokens': sum(lengths),
            'max_sequence_tokens': max(lengths),
            'predicted_region': 'All nucleotides after first complete codon, including terminal stop',
            'sequence_ids': seq_ids
        }
        write_json(output_dir / 'metadata.json', metadata)
        return cls(output_dir, mode=mode)


class Corpus(ArrowBioDataset):
    """Backward-compatible Corpus wrapper around ArrowBioDataset."""
    pass


def prepare(processed='data/processed', shard_size=1000):
    root = Path(processed)
    manifest = read_json(root / 'manifest.json')
    cohort_jsonl = root / 'cohort.jsonl'
    cohort_parquet = root / 'cohort.parquet'

    if cohort_jsonl.exists():
        if sha256(cohort_jsonl) != manifest['cohort_sha256']:
            raise ValueError('Cohort changed after selection')
        records = [json.loads(line) for line in cohort_jsonl.read_text(encoding='utf-8').splitlines()]
    elif cohort_parquet.exists():
        table = pq.read_table(cohort_parquet)
        records = table.to_pylist()
    else:
        raise FileNotFoundError(f"Cohort file not found in {root}")

    if not records or len({r['sequence_sha256'] for r in records}) != len(records):
        raise ValueError('Cohort must be nonempty and contain unique sequences')

    for mode in ('base', 'codon'):
        tokenizer = Tokenizer(mode)
        arrays = [np.asarray(tokenizer.encode(r['rna']), dtype=np.uint8) for r in records]
        lengths = [len(x) for x in arrays]
        if min(lengths) < (4 if mode == 'base' else 2):
            raise ValueError('Need a first codon and at least one target codon')
        path = root / mode
        path.mkdir(exist_ok=True)

        seq_ids = [r['sequence_sha256'] for r in records]
        rnas = [r['rna'] for r in records]
        encoded = [tokenizer.encode(r['rna']) for r in records]

        schema = pa.schema([
            ('sequence_sha256', pa.string()),
            ('rna', pa.string()),
            ('tokens', pa.list_(pa.uint8())),
            ('length', pa.int32()),
        ])

        table = pa.Table.from_pydict({
            'sequence_sha256': seq_ids,
            'rna': rnas,
            'tokens': encoded,
            'length': lengths,
        }, schema=schema)

        shards_dir = path / 'shards'
        shards_hash = ArrowBioDataset.to_parquet_shards(table, shards_dir, shard_size=shard_size)

        arrow_path = path / 'tokens.arrow'
        arrow_hash = ArrowBioDataset.to_arrow_ipc(table, arrow_path)

        np.concatenate(arrays).tofile(path / 'tokens.bin')
        np.save(path / 'offsets.npy', np.asarray([0, *np.cumsum(lengths)], dtype=np.int64))

        write_json(path / 'metadata.json', {
            'mode': mode,
            'vocabulary': tokenizer.vocabulary,
            'vocab_size': len(tokenizer.vocabulary),
            'sequences': len(records),
            'cohort_sha256': manifest['cohort_sha256'],
            'split': 'train',
            'tokens_sha256': sha256(path / 'tokens.bin'),
            'offsets_sha256': sha256(path / 'offsets.npy'),
            'arrow_sha256': arrow_hash,
            'parquet_shards_sha256': shards_hash,
            'tokens': sum(lengths),
            'max_sequence_tokens': max(lengths),
            'predicted_region': 'All nucleotides after first complete codon, including terminal stop',
            'sequence_ids': seq_ids
        })


