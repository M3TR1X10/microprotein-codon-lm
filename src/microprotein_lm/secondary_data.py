"""Frozen, equal-count views of the observed ATP synthase sequence pool."""
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from .io import read_json, sha256
from .secondary_io import write_json


def ranked(rows, seed):
    return sorted(rows, key=lambda r: hashlib.sha256(
        f"{seed}:{r['family']}:{r['sequence_sha256']}".encode()).hexdigest())


def round_robin(rows, count, seed):
    """Balance available species without copying a sequence to fill a quota."""
    groups = defaultdict(list)
    for row in ranked(rows, seed):
        groups[int(row['tax_id'])].append(row)
    result = []
    while len(result) < count:
        added = False
        for tax in sorted(groups, key=lambda t: hashlib.sha256(f'{seed}:taxon:{t}'.encode()).hexdigest()):
            if groups[tax] and len(result) < count:
                result.append(groups[tax].pop(0))
                added = True
        if not added:
            break
    return result


def family_quotas(human, count):
    available = Counter(r['family'] for r in human if r['family'] != 'ATP8')
    quotas = {family: 0 for family in sorted(available)}
    budget = min(sum(available.values()), count // 2)
    while sum(quotas.values()) < budget:
        for family in quotas:
            if quotas[family] < available[family] and sum(quotas.values()) < budget:
                quotas[family] += 1
    quotas['ATP8'] = count - sum(quotas.values())
    return {key: value for key, value in quotas.items() if value}


def select_views(rows, biology):
    seed = biology['selection_seed']
    mammals = set(biology['mammals'])
    nonmammals = set(biology['nonmammalian_vertebrates'])
    if len({r['sequence_sha256'] for r in rows}) != len(rows):
        raise ValueError('Input pool must already be globally RNA-deduplicated')
    if any(int(r['tax_id']) not in mammals | nonmammals for r in rows):
        raise ValueError('Unexpected organism in pool')
    human = [r for r in rows if int(r['tax_id']) == 9606]
    human_atp8 = ranked([r for r in human if r['family'] == 'ATP8'], seed)
    count = len(human_atp8)
    if count < 3:
        raise ValueError('Insufficient human ATP8 sequences for the declared experiment')
    quotas = family_quotas(human, count)
    if len(quotas) < 2:
        raise ValueError('No usable related human family')
    hc, mc, vc = [], [], []
    for family, quota in sorted(quotas.items()):
        hr = ranked([r for r in human if r['family'] == family], seed)
        hchosen = hr[:quota]
        common_count = min(quota, max(1, round(quota / 3)))
        common = hchosen[:common_count]
        other_mammals = [r for r in rows if r['family'] == family
                         and int(r['tax_id']) in mammals - {9606}]
        other_vertebrates = [r for r in rows if r['family'] == family
                            and int(r['tax_id']) in nonmammals]
        remaining = quota - common_count
        mchosen = round_robin(other_mammals, remaining, seed)
        missing = remaining - len(mchosen)
        mammal = common + mchosen + hchosen[common_count:common_count + missing]
        # Replace a fixed fraction of actual nonhuman mammal selections, rather
        # than silently changing family quotas when nonmammal references are absent.
        vadded = round_robin(other_vertebrates, (len(mchosen) + 1) // 2, seed)
        # Keep a balanced subset, including at least two mammal taxa if available.
        retained = round_robin(mchosen, len(mchosen) - len(vadded), seed)
        vertebrate = common + retained + vadded + hchosen[common_count:common_count + missing]
        if not len(hchosen) == len(mammal) == len(vertebrate) == quota:
            raise ValueError(f'Family quota cannot be filled: {family}')
        hc.extend(hchosen)
        mc.extend(mammal)
        vc.extend(vertebrate)
    if len({int(r['tax_id']) for r in mc} - {9606}) < 2:
        raise ValueError('Mammal view requires at least two other mammal species')
    if not ({int(r['tax_id']) for r in vc} & nonmammals):
        raise ValueError('Vertebrate view has no nonmammalian sequence')
    if len({int(r['tax_id']) for r in vc} & (mammals - {9606})) < 2:
        raise ValueError('Vertebrate view requires at least two other mammal species')
    views = {'human_atp8': human_atp8, 'human_complex': hc,
             'mammal_complex': mc, 'vertebrate_complex': vc,
             'full_human': human,
             'full_mammal': [r for r in rows if int(r['tax_id']) in mammals],
             'full_vertebrate': rows}
    # A common hash order makes the selected dataset independent of accession order.
    return {name: ranked(view, seed) for name, view in views.items()}, quotas


def describe(rows):
    return {'n_sequences': len(rows), 'n_unique_peptides': len({r['protein'] for r in rows}),
            'family_counts': dict(sorted(Counter(r['family'] for r in rows).items())),
            'taxon_counts': dict(sorted(Counter(str(r['tax_id']) for r in rows).items())),
            'family_taxon_counts': dict(sorted(Counter(f"{r['family']}:{r['tax_id']}" for r in rows).items())),
            'code_counts': dict(Counter(str(r['translation_table']) for r in rows)),
            'total_nt': sum(len(r['rna']) for r in rows),
            'target_nt': sum(len(r['rna']) - 3 for r in rows),
            'length_nt_counts': dict(sorted(Counter(str(len(r['rna'])) for r in rows).items()))}


def prepare_secondary(root, rows):
    root = Path(root)
    biology_path = root / 'configs/secondary/biology.json'
    biology = read_json(biology_path)
    views, quotas = select_views(rows, biology)
    manifests = {}
    for name, view in views.items():
        directory = root / 'data/processed/secondary' / name
        directory.mkdir(parents=True, exist_ok=True)
        cohort = directory / 'cohort.jsonl'
        # Keep large accession histories once in the hashed enriched pool. The
        # training reader needs biological labels/RNA, not tens of thousands of
        # duplicate-source dictionaries occupying RAM on every model run.
        compact = [{**{key: r[key] for key in ('rna', 'sequence_sha256', 'family',
                    'tax_id', 'translation_table', 'protein')},
                    'taxa': r.get('taxa', [r['tax_id']])} for r in view]
        content = ''.join(json.dumps(r, sort_keys=True) + '\n' for r in compact)
        content_hash = hashlib.sha256(content.encode()).hexdigest()
        if cohort.exists() and sha256(cohort) != content_hash:
            raise ValueError(f'Frozen cohort differs; create a new experiment version: {cohort}')
        if not cohort.exists():
            cohort.write_bytes(content.encode())
        manifest = {'stage': 'training_only', 'dataset': name,
                    'cohort_sha256': sha256(cohort), 'biology_config_sha256': sha256(biology_path),
                    'selection_source_sha256': sha256(Path(__file__)),
                    'selection_io_sha256': sha256(Path(__file__).with_name('secondary_io.py')),
                    'protocol_sha256': sha256(root / 'docs/secondary-protocol.md'),
                    'eligible_pool_sha256': sha256(root / 'data/raw/secondary/genomes/enriched-records.jsonl'),
                    'provenance_file': 'data/raw/secondary/genomes/enriched-records.jsonl',
                    'provenance_join_key': 'sequence_sha256',
                    'acquisition_report_sha256': sha256(root / 'reports/secondary/acquisition.json'),
                    'genome_acquisition_report_sha256': sha256(root / 'reports/secondary/genome-acquisition.json'),
                    'selection_seed': biology['selection_seed'], **describe(view)}
        if (directory / 'manifest.json').exists() and read_json(directory / 'manifest.json') != manifest:
            raise ValueError(f'Frozen selection metadata differs: {directory}')
        write_json(directory / 'manifest.json', manifest)
        manifests[name] = manifest
    primary = ['human_atp8', 'human_complex', 'mammal_complex', 'vertebrate_complex']
    sets = {name: {r['sequence_sha256'] for r in views[name]} for name in primary}
    report = {'stage': 'training_only', 'family_quotas': quotas,
              'views': manifests,
              'overlap_distinct_rna': {f'{a}__{b}': len(sets[a] & sets[b])
                                      for a in primary for b in primary},
              'interpretation': 'Equal counts, not equal lengths or independent biological replicates. Taxonomic comparisons replace sequences; family coverage is incomplete.'}
    write_json(root / 'reports/secondary/cohorts.json', report)
    return report
