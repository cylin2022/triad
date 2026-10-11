"""Curate verified tumor SNV vaccination assays and lock study-separated splits."""
import collections
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from curate import AA, POSITIVE, rejection

STUDIES = {
    '31515258': ('train', 'peptide_vaccination', None),
    '33879241': ('validation', 'DNA_or_SLP_vaccination', 'E0771'),
    '36950115': ('test', 'BMDC_vaccination_restimulation', 'LLC1'),
    '40586324': ('test', 'pooled_peptide_vaccination', 'CMT93'),
}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def canonical(seq):
    return isinstance(seq, str) and 8 <= len(seq) <= 11 and set(seq) <= AA


def distance(a, b):
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(cur[-1] + 1, prev[j] + 1, prev[j-1] + (x != y)))
        prev = cur
    return prev[-1]


def close(a, b):
    return abs(len(a)-len(b)) <= .2*max(len(a), len(b)) and distance(a, b) <= .2*max(len(a), len(b))


def source_tokens(detail):
    tokens = []
    for k in ['related_object__source_molecule_iri', 'related_object__molecule_parent_iri']:
        s = detail.get(k)
        if s:
            tokens.append(re.sub(r'\.\d+$', '', s.rstrip('/').split('/')[-1]).upper())
    name = detail.get('related_object__source_molecule')
    if name:
        tokens.append('name:' + re.sub(r'\s+', ' ', name).strip().lower())
    return ';'.join(sorted(set(tokens)))


def overlap(a, b):
    if set(a['Source_Tokens'].split(';')) - {''} & (set(b['Source_Tokens'].split(';')) - {''}):
        return 'same_source_protein'
    if any(close(x, y) for x in [a['Peptide'], a['WT']] for y in [b['Peptide'], b['WT']]):
        return 'mutant_or_WT_sequence_similarity'
    return None


def select_record(raw, detail, smith):
    reason = rejection(raw)
    if reason != 'eligible':
        return reason, None
    pmid = str(raw['pubmed_id'])
    if pmid not in STUDIES:
        return 'outside_reviewed_studies', None
    mut = raw['linear_sequence']
    if detail.get('assay_antigen__name') != mut:
        return 'assay_tests_WT_or_other_antigen', None
    if detail.get('adoptive_transfer__flag') != 'N' or detail.get('effector_cell__tcr_name'):
        return 'adoptive_or_defined_TCR', None
    if detail.get('second_in_vivo_process__process_type'):
        return 'second_in_vivo_process', None
    if detail.get('first_in_vivo_process__process_type') not in ['Administration in vivo', 'Vaccination']:
        return 'no_vaccination', None
    immunogen = detail.get('first_immunogen__name') or ''
    if mut not in immunogen or not set(immunogen) <= AA:
        return 'not_mutant_peptide_or_encoded_sequence_priming', None
    wt = detail.get('related_object__name')
    if not canonical(wt) or len(wt) != len(mut) or sum(x != y for x,y in zip(wt,mut)) != 1:
        return 'not_verified_single_substitution_WT_pair', None
    # Require a mutation relation except Smith's supplementary table independently
    # confirms the sequence pair (including entries missing IEDB relation labels).
    if pmid != '31515258' and detail.get('related_object__epitope_relation') not in ['in-frame neo-epitope', 'unspecified neo-epitope']:
        return 'missing_mutation_relation', None
    split, context, tumor = STUDIES[pmid]
    allele = raw['mhc_allele_name']
    if pmid == '31515258':
        matches = smith[(smith.Antigen == mut) & (smith.Reference == wt) & (smith.MHC_haplotype == allele)]
        if len(matches) != 1 or matches.iloc[0].Tumor_model not in ['B16F10','BBN963','MB49','UPPL1541']:
            return 'not_confirmed_somatic_mouse_tumor_Smith_table', None
        tumor = matches.iloc[0].Tumor_model
    if detail['assay__qualitative_measurement'] != raw['qualitative_measure']:
        raise ValueError('Live detail label changed relative to snapshot')
    return 'eligible', {
        'Peptide': mut, 'WT': wt, 'MHC': allele.replace('H2-', 'H-2-'),
        'Label': int(raw['qualitative_measure'] in POSITIVE), 'PMID': pmid,
        'Reference_ID': str(raw['reference_id']), 'Tumor_Model': tumor,
        'Context': context, 'Split': split, 'Assay_ID': raw['tcell_id'],
        'Source_Tokens': source_tokens(detail),
        'Mutation_Name': detail.get('epitope__reference_name'),
        'Restriction_Evidence': raw['mhc_allele_evidence'],
    }


def prepare():
    out = ROOT / 'locked'
    out.mkdir(exist_ok=True)
    if (out / 'split_lock.json').exists():
        raise FileExistsError('Split already locked; use a new version directory for a new experiment')
    rawpath = ROOT.parent / 'data/iedb_mouse_kb_db.jsonl'
    manifest = json.loads((ROOT.parent/'data/manifest.json').read_text())
    if sha(rawpath) != manifest['sha256']:
        raise ValueError('Raw snapshot checksum mismatch')
    raw = [json.loads(l) for l in rawpath.open()]
    detailpath = ROOT/'sources/iedb_tumor_details.json'
    details = {r['assay_id']: r for r in json.loads(detailpath.read_text())}
    detailmanifest = json.loads((ROOT/'sources/details_manifest.json').read_text())
    if sha(detailpath) != detailmanifest['sha256']:
        raise ValueError('Details checksum mismatch')
    smith = pd.read_excel(ROOT/'sources/smith_supp.xlsx')
    selected, audit = [], []
    for r in raw:
        if str(r['pubmed_id']) in STUDIES:
            reason, item = select_record(r, details[r['tcell_id']], smith)
            audit.append({'Assay_ID':r['tcell_id'], 'PMID':r['pubmed_id'], 'Peptide':r['linear_sequence'], 'Reason':reason})
            if item:
                selected.append(item)
    pd.DataFrame(audit).to_csv(out/'assay_curation.csv', index=False)
    groups = collections.defaultdict(list)
    for item in selected:
        groups[(item['PMID'], item['Peptide'], item['MHC'])].append(item)
    collapsed, conflicts = [], []
    for _, rr in sorted(groups.items()):
        if len({r['Label'] for r in rr}) != 1 or len({r['WT'] for r in rr}) != 1:
            conflicts.extend(rr)
            continue
        item = dict(rr[0])
        item['Assay_IDs'] = ';'.join(str(r['Assay_ID']) for r in rr)
        item.pop('Assay_ID')
        item['Source_Tokens'] = ';'.join(sorted({t for r in rr for t in r['Source_Tokens'].split(';') if t}))
        collapsed.append(item)
    df = pd.DataFrame(collapsed)
    df.to_csv(out/'before_purging.csv', index=False)
    pd.DataFrame(conflicts).to_csv(out/'conflicts.csv', index=False)
    splits = {s: df[df.Split == s].to_dict('records') for s in ['train','validation','test']}
    purged = []
    for lower, higher in [('validation','test'), ('train','test'), ('train','validation')]:
        keep = []
        for r in splits[lower]:
            hit = next(((h, overlap(r,h)) for h in splits[higher] if overlap(r,h)), None)
            if hit:
                h, why = hit
                purged.append({**r, 'Reason':why, 'Against_Split':higher, 'Against_Peptide':h['Peptide'], 'Against_PMID':h['PMID']})
            else:
                keep.append(r)
        splits[lower] = keep
    pd.DataFrame(purged).to_csv(out/'purged.csv', index=False)
    overlap_audit = {}
    for a,b in [('train','validation'),('train','test'),('validation','test')]:
        hits = [(x['Peptide'], y['Peptide']) for x in splits[a] for y in splits[b] if overlap(x,y)]
        ref_hits = {r['PMID'] for r in splits[a]} & {r['PMID'] for r in splits[b]}
        tumor_hits = {r['Tumor_Model'] for r in splits[a]} & {r['Tumor_Model'] for r in splits[b]}
        if hits or ref_hits or tumor_hits:
            raise ValueError(f'Cross-split overlap: {a}/{b}')
        overlap_audit[a+'__'+b] = {'sequence_or_source':len(hits), 'studies':len(ref_hits), 'tumor_models':len(tumor_hits)}
    stats = {}
    for name, rr in splits.items():
        part = pd.DataFrame(rr, columns=df.columns)
        part.to_csv(out/(name+'.csv'), index=False)
        stats[name] = {'n':len(part), 'positive':int(part.Label.sum()), 'negative':int((part.Label==0).sum()),
                       'studies': sorted(part.PMID.unique()), 'tumor_models':sorted(part.Tumor_Model.unique()),
                       'by_allele':part.groupby('MHC').Label.agg(['count','sum']).to_dict('index')}
    lock = {'locked_utc':datetime.now(timezone.utc).isoformat(), 'user_scope':'tumor_neoantigen',
            'new_model_fit_has_started':False, 'prior_baseline_exposure':True,
            'statistics':stats, 'overlap_audit':overlap_audit,
            'curation_counts':dict(collections.Counter(x['Reason'] for x in audit)),
            'conflicting_assays':len(conflicts), 'purged_pairs':len(purged),
            'inputs':{str(p.relative_to(ROOT.parent)):sha(p) for p in [rawpath,detailpath,ROOT/'sources/smith_supp.xlsx']},
            'files':{p.name:sha(p) for p in out.glob('*.csv')},
            'protocol_sha256':sha(ROOT/'PROTOCOL.md'), 'prepare_sha256':sha(Path(__file__))}
    (out/'split_lock.json').write_text(json.dumps(lock,indent=2)+'\n')
    print(json.dumps(lock,indent=2))


if __name__ == '__main__':
    prepare()
