"""Canonical JSON bytes for portable secondary manifests and identity hashes."""
import json
from pathlib import Path


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes((json.dumps(value, indent=2, sort_keys=True) + '\n').encode('utf-8'))
