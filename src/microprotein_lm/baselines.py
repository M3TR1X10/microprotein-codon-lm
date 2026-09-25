"""Resubstitution diagnostics fitted and scored only on the training corpus."""
import math
import numpy as np

from .data import Corpus
from .io import write_json


def training_baselines(root='data/processed'):
    report = {'scope':'TRAINING resubstitution only; these are not generalization estimates','arms':{}}
    for mode in ('base','codon'):
        corpus = Corpus(root,mode)
        vocab = len(corpus.tokenizer.vocabulary)
        triples = []
        for i in range(len(corpus)):
            seq = corpus.sequence(i)
            first = 3 if mode == 'base' else 1
            triples.extend((pos,int(seq[pos-1]),int(seq[pos])) for pos in range(first,len(seq)))
        counts = np.zeros(vocab)
        bigram = np.zeros((vocab,vocab))
        positional = np.zeros((max(p for p,_,_ in triples)+1,vocab))
        for pos,prev,token in triples:
            counts[token] += 1
            bigram[prev,token] += 1
            positional[pos,token] += 1
        smooth = 0.5
        distributions = {
            'unigram':(counts+smooth)/(counts.sum()+smooth*vocab),
            'bigram':(bigram+smooth)/(bigram.sum(axis=1,keepdims=True)+smooth*vocab),
            'position':(positional+smooth)/(positional.sum(axis=1,keepdims=True)+smooth*vocab)}
        metrics = {'uniform':{'bits_per_base':math.log2(vocab)/corpus.tokenizer.width}}
        for name,p in distributions.items():
            nll, correct = 0.,0
            for pos,prev,token in triples:
                dist = p if name == 'unigram' else p[prev if name == 'bigram' else pos]
                nll -= math.log(float(dist[token]))
                correct += int(dist.argmax() == token)
            metrics[name] = {'bits_per_base':nll/(len(triples)*corpus.tokenizer.width*math.log(2)),
                'token_accuracy':correct/len(triples)}
        report['arms'][mode] = {'target_bases':len(triples)*corpus.tokenizer.width, 'metrics':metrics}
    write_json('reports/training-baselines.json',report)
    return report
