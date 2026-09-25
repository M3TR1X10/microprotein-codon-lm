import argparse
import json

from .io import read_json


def main():
    parser = argparse.ArgumentParser(description='Training-only microprotein experiments')
    sub = parser.add_subparsers(dest='command',required=True)
    p = sub.add_parser('discover',help='Download the candidate inventory with source manifests')
    p.add_argument('--raw',default='data/raw')
    p = sub.add_parser('select',help='Fetch, quality-check and rank all candidates')
    p.add_argument('--raw',default='data/raw')
    p.add_argument('--processed',default='data/processed')
    p.add_argument('--min-identity',type=float,default=0.95)
    for name in ('prepare','baselines'):
        p = sub.add_parser(name)
        p.add_argument('--processed',default='data/processed')
    p = sub.add_parser('train')
    p.add_argument('--config',default='configs/pilot.json')
    p.add_argument('--mode',choices=['base','codon'],required=True)
    p.add_argument('--seed',type=int,default=17)
    p.add_argument('--processed',default='data/processed')
    p.add_argument('--output',required=True)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--stop-after',type=int,help='Stop after this many total updates; retain resumable state')
    args = parser.parse_args()
    if args.command == 'discover':
        from .acquire import discover
        discover(args.raw)
    elif args.command == 'select':
        from .quality import select
        print(json.dumps(select(args.raw,args.processed,args.min_identity),indent=2))
    elif args.command == 'prepare':
        from .data import prepare
        prepare(args.processed)
    elif args.command == 'baselines':
        from .baselines import training_baselines
        print(json.dumps(training_baselines(args.processed),indent=2))
    elif args.command == 'train':
        from .train import train
        train(read_json(args.config),args.mode,args.seed,args.processed,args.output,args.resume,args.stop_after)


if __name__ == '__main__':
    main()
