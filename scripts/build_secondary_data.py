"""Acquire the declared panel, then freeze all eligible pools and matched views."""
import argparse
from pathlib import Path
from microprotein_lm.io import read_json
from microprotein_lm.secondary_acquire import acquire_secondary
from microprotein_lm.secondary_data import prepare_secondary
from microprotein_lm.secondary_genomes import enrich_genomes

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', default='.')
    args = parser.parse_args()
    root = Path(args.root)
    report = prepare_secondary(root, enrich_genomes(root, acquire_secondary(root)))
    print({name: value['n_sequences'] for name, value in report['views'].items()}, flush=True)
