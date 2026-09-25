"""Decoupled Relational Provenance & Content-Addressable Storage (CAS).

Manages CAS storage of raw XML payloads (data/raw_cas/<sha256>.xml),
observation ledger recording raw acquisition receipts, and canonical biological
sequence records decoupling observations from canonical entities.
"""

import hashlib
import json
from pathlib import Path
from typing import Dict, List, Optional, Union

from .io import now


class CASStore:
    """Content-Addressable Storage (CAS) for raw payloads."""

    def __init__(self, storage_dir: Union[str, Path] = "data/raw_cas"):
        self.storage_dir = Path(storage_dir)
        self.storage_dir.mkdir(parents=True, exist_ok=True)

    def store(self, content: Union[str, bytes]) -> str:
        """Saves payload content to data/raw_cas/<sha256>.xml and returns sha256 hash."""
        if isinstance(content, str):
            payload_bytes = content.encode('utf-8')
        else:
            payload_bytes = content

        sha256_hash = hashlib.sha256(payload_bytes).hexdigest()
        file_path = self.get_path(sha256_hash)
        if not file_path.exists():
            temp_path = file_path.with_suffix('.tmp')
            temp_path.write_bytes(payload_bytes)
            temp_path.replace(file_path)

        return sha256_hash

    def retrieve(self, sha256_hash: str) -> bytes:
        """Retrieves raw bytes for a given SHA256 hash."""
        file_path = self.get_path(sha256_hash)
        if not file_path.exists():
            raise FileNotFoundError(f"CAS payload with SHA256 hash {sha256_hash} not found at {file_path}")
        return file_path.read_bytes()

    def retrieve_text(self, sha256_hash: str, encoding: str = 'utf-8') -> str:
        """Retrieves raw string text for a given SHA256 hash."""
        return self.retrieve(sha256_hash).decode(encoding)

    def exists(self, sha256_hash: str) -> bool:
        """Returns True if the CAS payload exists."""
        return self.get_path(sha256_hash).exists()

    def get_path(self, sha256_hash: str) -> Path:
        """Returns file path for SHA256 hash."""
        return self.storage_dir / f"{sha256_hash}.xml"


class ObservationLedger:
    """Ledger recording raw observation instances and acquisition receipts."""

    def __init__(self, ledger_path: Union[str, Path] = "data/raw_cas/observation_ledger.jsonl"):
        self.ledger_path = Path(ledger_path)
        self.ledger_path.parent.mkdir(parents=True, exist_ok=True)

    def record_observation(
        self,
        accession: str,
        source_url: str,
        raw_sha256: str,
        metadata: Optional[Dict] = None
    ) -> Dict:
        """Appends an observation entry to the ledger.

        Args:
            accession: Target sequence or protein accession ID.
            source_url: URL endpoint used for fetching.
            raw_sha256: CAS SHA256 hash of the retrieved payload.
            metadata: Optional additional metadata dict.

        Returns:
            Recorded observation dict.
        """
        obs_id = f"obs_{raw_sha256[:12]}_{accession}"
        entry = {
            "observation_id": obs_id,
            "accession": accession,
            "source_url": source_url,
            "raw_sha256": raw_sha256,
            "timestamp_utc": now(),
            "metadata": metadata or {},
        }

        with self.ledger_path.open('a', encoding='utf-8') as f:
            f.write(json.dumps(entry, sort_keys=True) + '\n')

        return entry

    def get_observations(self, accession: Optional[str] = None) -> List[Dict]:
        """Reads all recorded observations, optionally filtering by accession."""
        if not self.ledger_path.exists():
            return []

        results = []
        with self.ledger_path.open('r', encoding='utf-8') as f:
            for line in f:
                if not line.strip():
                    continue
                record = json.loads(line)
                if accession is None or record.get('accession') == accession:
                    results.append(record)
        return results


class CanonicalSequences:
    """Store for canonical biological sequence records."""

    def __init__(self, canonical_path: Union[str, Path] = "data/processed/canonical_sequences.jsonl"):
        self.canonical_path = Path(canonical_path)
        self.canonical_path.parent.mkdir(parents=True, exist_ok=True)

    def register_sequence(
        self,
        sequence_sha256: str,
        rna: str,
        protein: str,
        initiation_codon: str,
        stop_codon: str,
        metadata: Dict,
        source_accessions: List[str],
        observation_ids: Optional[List[str]] = None
    ) -> Dict:
        """Upserts a canonical sequence record.

        Args:
            sequence_sha256: SHA256 hash of canonical RNA sequence.
            rna: Canonical RNA sequence.
            protein: Translated protein sequence.
            initiation_codon: Explicit initiation codon (e.g. 'AUG').
            stop_codon: Explicit stop codon (e.g. 'UAA').
            metadata: Associated metadata dict.
            source_accessions: List of source accessions sharing this canonical sequence.
            observation_ids: Optional list of linked observation IDs.

        Returns:
            The canonical sequence record dict.
        """
        existing_records = self.list_all()
        by_hash = {r['sequence_sha256']: r for r in existing_records}

        if sequence_sha256 in by_hash:
            rec = by_hash[sequence_sha256]
            acc_set = set(rec.get('source_accessions', [])) | set(source_accessions)
            obs_set = set(rec.get('observation_ids', [])) | set(observation_ids or [])
            rec['source_accessions'] = sorted(list(acc_set))
            rec['observation_ids'] = sorted(list(obs_set))
            rec['metadata'].update(metadata)
            rec['updated_utc'] = now()
            result = rec
        else:
            result = {
                "sequence_sha256": sequence_sha256,
                "rna": rna,
                "protein": protein,
                "initiation_codon": initiation_codon,
                "stop_codon": stop_codon,
                "source_accessions": sorted(list(set(source_accessions))),
                "observation_ids": sorted(list(set(observation_ids or []))),
                "metadata": metadata or {},
                "created_utc": now(),
                "updated_utc": now(),
            }
            by_hash[sequence_sha256] = result

        # Write back sorted by sequence_sha256
        sorted_records = sorted(by_hash.values(), key=lambda r: r['sequence_sha256'])
        with self.canonical_path.open('w', encoding='utf-8') as f:
            for r in sorted_records:
                f.write(json.dumps(r, sort_keys=True) + '\n')

        return result

    def get_sequence(self, sequence_sha256: str) -> Optional[Dict]:
        """Retrieves a canonical sequence by its SHA256 hash."""
        records = self.list_all()
        for r in records:
            if r['sequence_sha256'] == sequence_sha256:
                return r
        return None

    def list_all(self) -> List[Dict]:
        """Lists all canonical sequence records."""
        if not self.canonical_path.exists():
            return []

        records = []
        with self.canonical_path.open('r', encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    records.append(json.loads(line))
        return records
