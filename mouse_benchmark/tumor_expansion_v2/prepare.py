"""Curate expanded mouse tumor neoantigen cohorts and lock study-separated splits."""
import collections
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent
sys.path.insert(0, str(ROOT.parent))
from curate import AA, POSITIVE, rejection

STUDIES_V1 = {
    '31515258': ('train', 'peptide_vaccination', None),
    '33879241': ('validation', 'DNA_or_SLP_vaccination', 'E0771'),
    '36950115': ('test_frozen_prior', 'BMDC_vaccination_restimulation', 'LLC1'),
    '40586324': ('test_frozen_prior', 'pooled_peptide_vaccination', 'CMT93'),
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
            cur.append(min(cur[-1] + 1, prev[j] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return prev[-1]


def close(a, b):
    return abs(len(a) - len(b)) <= 0.2 * max(len(a), len(b)) and distance(a, b) <= 0.2 * max(len(a), len(b))


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
    a_tokens = set(a.get('Source_Tokens', '').split(';')) - {''}
    b_tokens = set(b.get('Source_Tokens', '').split(';')) - {''}
    if a_tokens & b_tokens:
        return 'same_source_protein'
    if any(close(x, y) for x in [a['Peptide'], a['WT']] for y in [b['Peptide'], b['WT']]):
        return 'mutant_or_WT_sequence_similarity'
    return None


def select_record_iedb(raw, detail, smith):
    reason = rejection(raw)
    if reason != 'eligible':
        return reason, None
    pmid = str(raw['pubmed_id'])
    if pmid not in STUDIES_V1 and pmid != '31292145':
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
    if not canonical(wt) or len(wt) != len(mut) or sum(x != y for x, y in zip(wt, mut)) != 1:
        return 'not_verified_single_substitution_WT_pair', None
    if pmid not in ['31515258', '31292145'] and detail.get('related_object__epitope_relation') not in [
        'in-frame neo-epitope', 'unspecified neo-epitope'
    ]:
        return 'missing_mutation_relation', None

    allele = raw['mhc_allele_name']
    if pmid == '31515258':
        split, context, tumor = STUDIES_V1[pmid]
        matches = smith[(smith.Antigen == mut) & (smith.Reference == wt) & (smith.MHC_haplotype == allele)]
        if len(matches) != 1 or matches.iloc[0].Tumor_model not in ['B16F10', 'BBN963', 'MB49', 'UPPL1541']:
            return 'not_confirmed_somatic_mouse_tumor_Smith_table', None
        tumor = matches.iloc[0].Tumor_model
    elif pmid == '31292145':
        split, context, tumor = 'train_expansion', 'pooled_peptide_vaccination', 'MC-38'
    else:
        split, context, tumor = STUDIES_V1[pmid]

    return 'eligible', {
        'Peptide': mut,
        'WT': wt,
        'MHC': allele.replace('H2-', 'H-2-'),
        'Label': int(raw['qualitative_measure'] in POSITIVE),
        'PMID': pmid,
        'Reference_ID': str(raw['reference_id']),
        'Tumor_Model': tumor,
        'Context': context,
        'Split': split,
        'Assay_ID': raw['tcell_id'],
        'Source_Tokens': source_tokens(detail),
        'Mutation_Name': detail.get('epitope__reference_name', mut),
        'Restriction_Evidence': raw['mhc_allele_evidence'],
    }


def prepare():
    out = ROOT / 'locked'
    out.mkdir(exist_ok=True)
    if (out / 'split_lock.json').exists():
        raise FileExistsError('Split already locked; use a new version directory for a new experiment')

    rawpath = ROOT.parent / 'data/iedb_mouse_kb_db.jsonl'
    manifest = json.loads((ROOT.parent / 'data/manifest.json').read_text())
    if sha(rawpath) != manifest['sha256']:
        raise ValueError('Raw snapshot checksum mismatch')
    raw = [json.loads(l) for l in rawpath.open()]

    detailpath = ROOT / 'sources/iedb_tumor_details.json'
    details = {r['assay_id']: r for r in json.loads(detailpath.read_text())}
    smith = pd.read_excel(ROOT / 'sources/smith_supp.xlsx')

    # Curate IEDB records (Smith 2019, E0771, LLC1, CMT93, Fernando 2019)
    selected, audit = [], []
    for r in raw:
        pmid_str = str(r['pubmed_id'])
        if pmid_str in STUDIES_V1 or pmid_str == '31292145':
            reason, item = select_record_iedb(r, details[r['tcell_id']], smith)
            audit.append({
                'Assay_ID': r['tcell_id'],
                'PMID': r['pubmed_id'],
                'Peptide': r['linear_sequence'],
                'Reason': reason
            })
            if item:
                selected.append(item)

    # Collapse IEDB replicates
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
    df_iedb = pd.DataFrame(collapsed)

    # Curate Capietto et al. 2020 (PMC7144530 / PMID 31940002)
    excel_path = ROOT / 'sources/JEM_20190179_DataS1.xlsx'
    df_jem = pd.read_excel(excel_path, sheet_name='immunogenicity_modeling-dat_mod')
    mouse_jem = df_jem[df_jem['assay_org'] == 'mouse'].copy()

    capietto_c57bl6 = []
    capietto_balbc = []
    for _, row in mouse_jem.iterrows():
        mut = str(row['epitope']).strip()
        wt = str(row['wtCounterpart']).strip()
        allele = str(row['allele']).strip()
        label = int(row['response'] == 'CD8+')
        tumor = str(row['samples_CNAME']).strip()
        gene = str(row['antigen_gene_symbol']).strip()
        mut_id = str(row['mutation']).strip()
        aa_change = str(row['mutation_aa_change']).strip()
        record = {
            'Peptide': mut,
            'WT': wt,
            'MHC': allele,
            'Label': label,
            'PMID': '31940002',
            'Reference_ID': 'PMC7144530',
            'Tumor_Model': tumor,
            'Context': 'SLP_vaccination_polyIC_antiCD40',
            'Split': 'test_capietto_holdout',
            'Source_Tokens': f'name:{gene.lower()}',
            'Mutation_Name': f'{gene}_{aa_change}',
            'Restriction_Evidence': 'MHC binding prediction',
            'Assay_IDs': f'JEM_{mut_id}',
        }
        if allele in ['H-2-Kb', 'H-2-Db']:
            capietto_c57bl6.append(record)
        else:
            capietto_balbc.append(record)

    df_cap_b6 = pd.DataFrame(capietto_c57bl6)
    df_cap_balbc = pd.DataFrame(capietto_balbc)

    # Combine training candidates: Smith (100) + Fernando (13)
    train_v1 = df_iedb[df_iedb.Split == 'train'].copy()
    train_fernando = df_iedb[df_iedb.Split == 'train_expansion'].copy()
    train_fernando['Split'] = 'train'
    df_train = pd.concat([train_v1, train_fernando], ignore_index=True)

    df_val = df_iedb[df_iedb.Split == 'validation'].copy()
    df_test_prior = df_iedb[df_iedb.Split == 'test_frozen_prior'].copy()

    # Purge sequence overlap between Capietto holdout and train/val/prior_test
    # Note: 6 MC-38 peptides in Capietto are identical to Fernando MC-38 peptides.
    # To maintain 100% strict data leakage isolation, purge those from Capietto holdout.
    all_dev_records = df_train.to_dict('records') + df_val.to_dict('records')
    purged_capietto = []
    kept_capietto = []
    for r in df_cap_b6.to_dict('records'):
        hit = next(((d, overlap(r, d)) for d in all_dev_records if overlap(r, d)), None)
        if hit:
            d, why = hit
            purged_capietto.append({**r, 'Reason': why, 'Against_Peptide': d['Peptide'], 'Against_PMID': d['PMID']})
        else:
            kept_capietto.append(r)

    df_cap_holdout = pd.DataFrame(kept_capietto)
    pd.DataFrame(purged_capietto).to_csv(out / 'purged_capietto_from_holdout.csv', index=False)

    # Check zero overlap between holdout test splits and training
    for name, holdout_df in [('test_frozen_prior', df_test_prior), ('test_capietto_holdout', df_cap_holdout)]:
        hits = [(x['Peptide'], y['Peptide']) for x in df_train.to_dict('records') for y in holdout_df.to_dict('records') if overlap(x, y)]
        if hits:
            raise ValueError(f'Cross-split leakage detected between train and {name}: {hits}')

    # Save locked CSV tables
    df_train.to_csv(out / 'train.csv', index=False)
    df_val.to_csv(out / 'validation.csv', index=False)
    df_test_prior.to_csv(out / 'test_frozen_prior.csv', index=False)
    df_cap_holdout.to_csv(out / 'test_capietto_holdout.csv', index=False)
    df_cap_b6.to_csv(out / 'test_capietto_all.csv', index=False)
    df_cap_balbc.to_csv(out / 'exploratory_balbc.csv', index=False)

    # Also save combined training cohort for large-negative training experiment (Smith + Fernando + Capietto MC-38)
    combined_train = pd.concat([df_train, df_cap_holdout], ignore_index=True)
    combined_train.to_csv(out / 'train_combined_capietto.csv', index=False)

    pd.DataFrame(audit).to_csv(out / 'assay_curation.csv', index=False)
    pd.DataFrame(conflicts).to_csv(out / 'conflicts.csv', index=False)

    # Compute statistics
    stats = {
        'train': {
            'n': len(df_train),
            'positive': int(df_train.Label.sum()),
            'negative': int((df_train.Label == 0).sum()),
            'prevalence': float(df_train.Label.mean()),
            'studies': sorted(df_train.PMID.unique()),
            'tumor_models': sorted(df_train.Tumor_Model.unique()),
        },
        'validation': {
            'n': len(df_val),
            'positive': int(df_val.Label.sum()),
            'negative': int((df_val.Label == 0).sum()),
            'prevalence': float(df_val.Label.mean()),
            'studies': sorted(df_val.PMID.unique()),
            'tumor_models': sorted(df_val.Tumor_Model.unique()),
        },
        'test_frozen_prior': {
            'n': len(df_test_prior),
            'positive': int(df_test_prior.Label.sum()),
            'negative': int((df_test_prior.Label == 0).sum()),
            'prevalence': float(df_test_prior.Label.mean()),
            'studies': sorted(df_test_prior.PMID.unique()),
            'tumor_models': sorted(df_test_prior.Tumor_Model.unique()),
        },
        'test_capietto_holdout': {
            'n': len(df_cap_holdout),
            'positive': int(df_cap_holdout.Label.sum()),
            'negative': int((df_cap_holdout.Label == 0).sum()),
            'prevalence': float(df_cap_holdout.Label.mean()),
            'studies': sorted(df_cap_holdout.PMID.unique()),
            'tumor_models': sorted(df_cap_holdout.Tumor_Model.unique()),
        },
        'exploratory_balbc': {
            'n': len(df_cap_balbc),
            'positive': int(df_cap_balbc.Label.sum()),
            'negative': int((df_cap_balbc.Label == 0).sum()),
            'prevalence': float(df_cap_balbc.Label.mean()),
            'studies': sorted(df_cap_balbc.PMID.unique()),
            'tumor_models': sorted(df_cap_balbc.Tumor_Model.unique()),
        },
    }

    lock = {
        'locked_utc': datetime.now(timezone.utc).isoformat(),
        'experiment': 'mouse_tumor_neoantigen_expansion_v2',
        'statistics': stats,
        'purged_from_capietto_holdout': len(purged_capietto),
        'inputs': {
            'iedb_mouse_kb_db.jsonl': sha(rawpath),
            'iedb_tumor_details.json': sha(detailpath),
            'smith_supp.xlsx': sha(ROOT / 'sources/smith_supp.xlsx'),
            'JEM_20190179_DataS1.xlsx': sha(excel_path),
        },
        'files': {p.name: sha(p) for p in out.glob('*.csv')},
        'protocol_sha256': sha(ROOT / 'PROTOCOL.md'),
        'prepare_sha256': sha(Path(__file__)),
    }
    (out / 'split_lock.json').write_text(json.dumps(lock, indent=2) + '\n')
    print('Locked splits generated successfully:')
    print(json.dumps(stats, indent=2))


if __name__ == '__main__':
    prepare()
