import hashlib
import json
from copy import deepcopy

import numpy as np
import pytest
import torch
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, BeforePosition
from Bio.SeqRecord import SeqRecord

from microprotein_lm.data import Corpus, prepare
from microprotein_lm.io import write_json, sha256
from microprotein_lm.model import Decoder, ModelConfig
from microprotein_lm.quality import validate_record
from microprotein_lm.tokenization import Tokenizer
from microprotein_lm.train import lr_at, train


def fixture_record(dna='ATGTGATAA', protein='MW'):
    record = SeqRecord(Seq(dna),id='TEST.1')
    record.features = [SeqFeature(FeatureLocation(0,len(dna)),type='source',qualifiers={'db_xref':['taxon:9606']}),
        SeqFeature(FeatureLocation(0,len(dna)),type='CDS',qualifiers={
            'protein_id':['TEST.1'],'codon_start':['1'],'transl_table':['2'],'translation':[protein]})]
    metadata = {'protein_id':'TEST.1','parent_accession':'TEST_PARENT','base_count':str(len(dna)),
        'sequence_md5':hashlib.md5(dna.encode()).hexdigest()}
    reference = {'gene':'MT-EXAMPLE','uniprot':'TEST','reference_protein':protein}
    return record,metadata,reference


def test_vocabularies_are_exact_and_round_trip():
    for mode,size in [('base',4),('codon',64)]:
        t = Tokenizer(mode)
        assert len(set(t.vocabulary)) == size
        text = 'AUGUUUUAA'
        assert t.decode(t.encode(text)) == text
        for invalid in ['', 'ATG', 'ANN', 'aug']:
            with pytest.raises(ValueError):
                t.encode(invalid)
    with pytest.raises(ValueError):
        Tokenizer('codon').encode('AUGU')
    assert all(c in Tokenizer('codon').vocabulary for c in ['UAA','UAG','UGA','AGA','AGG'])


def test_mitochondrial_code_and_complete_cds():
    r,m,ref = fixture_record()
    out,why = validate_record(r,m,ref)
    assert why is None and out['protein'] == 'MW'  # UGA is tryptophan in table 2
    r.features[1].qualifiers['transl_table'] = ['1']
    assert validate_record(r,m,ref)[1] == 'unexpected_genetic_code'


@pytest.mark.parametrize('dna,reason',[
    ('ATGNNNTAA','ambiguous_bases'), ('ATGTGATA','incomplete_codon'),
    ('ATGTAATAA','invalid_complete_cds'), ('ATGTGATGG','invalid_complete_cds')])
def test_invalid_cds_rejected(dna,reason):
    assert validate_record(*fixture_record(dna))[1] == reason


def test_partial_taxon_and_translation_checks():
    r,m,ref = fixture_record()
    r.features[1].location = FeatureLocation(BeforePosition(0),len(r))
    assert validate_record(r,m,ref)[1] == 'partial_location'
    r,m,ref = fixture_record()
    r.features[0].qualifiers['db_xref'] = ['taxon:10090']
    assert validate_record(r,m,ref)[1] == 'wrong_taxon'
    r,m,ref = fixture_record()
    r.features[1].qualifiers['translation'] = ['MM']
    assert validate_record(r,m,ref)[1] == 'translation_mismatch'


@pytest.fixture
def corpus_root(tmp_path):
    # Synthetic sequences are SOFTWARE TEST FIXTURES, never research data.
    root = tmp_path/'data'
    root.mkdir()
    records = [{'rna':r,'sequence_sha256':hashlib.sha256(r.encode()).hexdigest()}
        for r in ['AUGUGAUAA','AUGCCCUAG','AUGAAACC CUAA'.replace(' ','')]]
    path = root/'cohort.jsonl'
    path.write_text(''.join(json.dumps(r)+'\n' for r in records),encoding='utf-8')
    write_json(root/'manifest.json',{'cohort_sha256':sha256(path)})
    prepare(root)
    return root


def test_paired_targets_and_boundaries(corpus_root):
    base,codon = Corpus(corpus_root,'base'),Corpus(corpus_root,'codon')
    xb,yb = base.batch([0,1,2])
    xc,yc = codon.batch([0,1,2])
    assert (yb != -100).sum() == 3*(yc != -100).sum()
    assert (yb[:,:2] == -100).all()
    for i in range(3):
        b = base.tokenizer.decode(yb[i][yb[i] != -100].tolist())
        c = codon.tokenizer.decode(yc[i][yc[i] != -100].tolist())
        assert b == c
    assert (yb[0,8:] == -100).all()
    assert base.metadata['sequence_ids'] == codon.metadata['sequence_ids']


def test_causal_attention_cannot_read_future():
    torch.manual_seed(5)
    model = Decoder(ModelConfig(4,12,1,d_model=16,n_heads=2,n_layers=1,dropout=0)).eval()
    x = torch.tensor([[0,1,2,3,0,1]])
    changed = x.clone()
    changed[:,4:] = 3
    assert torch.allclose(model(x)[:,:4],model(changed)[:,:4],atol=1e-7)


def tiny_config():
    return {'stage':'training_only','model':{'d_model':16,'n_heads':2,'n_layers':1,'dropout':0.1},
        'max_nt_context':15,'batch_size':2,'gradient_accumulation':2,'updates':4,
        'warmup_updates':1,'learning_rate':0.001,'min_learning_rate':0.0001,'weight_decay':0.01,
        'gradient_clip':1.,'diagnostic_interval':2,'checkpoint_interval':2,
        'device':'cpu','precision':'float32','deterministic':True,'cpu_threads':1}


def test_schedule_uses_optimizer_updates():
    config = tiny_config()
    values = [lr_at(i,config) for i in range(config['updates'])]
    assert values[-1] == pytest.approx(config['min_learning_rate'])
    assert min(values) > 0 and max(values) <= config['learning_rate']
    config['min_learning_rate'] = 1.
    with pytest.raises(ValueError):
        lr_at(0,config)


def test_checkpoint_resume_is_exact_and_rejects_changed_config(corpus_root,tmp_path):
    config = tiny_config()
    full,part = tmp_path/'full',tmp_path/'part'
    train(config,'codon',7,corpus_root,full)
    train(config,'codon',7,corpus_root,part,stop_after=2)
    train(config,'codon',7,corpus_root,part,resume=True)
    a = torch.load(full/'checkpoint.pt',weights_only=True)
    b = torch.load(part/'checkpoint.pt',weights_only=True)
    assert a['bases_seen'] == b['bases_seen']
    for key in a['model']:
        assert torch.equal(a['model'][key],b['model'][key]),key
    changed = deepcopy(config)
    changed['learning_rate'] = 0.002
    with pytest.raises(ValueError,match='identical'):
        train(changed,'codon',7,corpus_root,part,resume=True)


def test_encoded_data_tampering_is_rejected(corpus_root):
    path = corpus_root/'base'/'tokens.bin'
    path.write_bytes(path.read_bytes()+b'\x00')
    with pytest.raises(ValueError,match='checksum'):
        Corpus(corpus_root,'base')
