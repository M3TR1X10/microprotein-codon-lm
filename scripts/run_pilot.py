"""Run all six prespecified pilots, serially, without selecting a winner."""
import argparse
from pathlib import Path
from microprotein_lm.io import read_json, write_json
from microprotein_lm.train import train

parser = argparse.ArgumentParser()
parser.add_argument('--config',default='configs/pilot.json')
parser.add_argument('--processed',default='data/processed')
parser.add_argument('--runs',default='runs/pilot')
args = parser.parse_args()
config = read_json(args.config)
summaries = []
for seed in config['seeds']:
    for mode in ('base','codon'):
        output = Path(args.runs)/f'{mode}-seed{seed}'
        summaries.append(train(config,mode,seed,args.processed,output))
        write_json('reports/pilot-results.json',summaries)
