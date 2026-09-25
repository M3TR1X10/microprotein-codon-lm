"""Restore the recorded public raw snapshot; refuse changed upstream bytes."""
from pathlib import Path

from microprotein_lm.acquire import Download
from microprotein_lm.io import read_json, sha256


if __name__ == '__main__':
    registry = read_json('reports/secondary/source-registry.json')
    download = Download('data/raw')
    for index, entry in enumerate(registry['sources'], 1):
        relative = Path(entry['file'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Unsafe raw snapshot path')
        path = download.get(relative.as_posix(), entry['url'])
        if sha256(path) != entry['sha256']:
            raise ValueError(f'Upstream content changed: {relative}; restore your archived raw snapshot')
        if index % 25 == 0:
            print(f"Verified {index}/{len(registry['sources'])} source files", flush=True)
    print('Raw snapshot verified. Run scripts/build_secondary_data.py to reconstruct frozen cohorts.')
