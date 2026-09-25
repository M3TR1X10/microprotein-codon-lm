"""Acquire disjoint observed CDS without training or scoring a model."""
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from microprotein_lm.gap_acquire import acquire_gap

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--include-outer-parents', action='store_true')
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    acquire_gap(ROOT, args.include_outer_parents, args.freeze)
