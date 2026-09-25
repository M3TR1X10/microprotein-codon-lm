"""Restore the recorded public-source snapshot, failing on content changes."""
from pathlib import Path
from microprotein_lm.acquire import Download
from microprotein_lm.io import read_json, sha256

download = Download('data/raw')
registry = read_json('reports/source-registry.json')
for entry in registry['sources']:
    filename = entry['file']
    relative = Path(filename)
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('Unsafe snapshot path')
    path = download.get(filename,entry['url'])
    if sha256(path) != entry['sha256']:
        raise ValueError(f'Upstream content changed for {filename}; restore from your archived raw snapshot')
print(f"Verified {len(registry['sources'])} original source files")
