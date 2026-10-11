"""Fit fixed research models only after verifying immutable split hashes."""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from prepare import ROOT, close, overlap, sha
sys.path.insert(0, str(ROOT.parent))
from predict import add_scores
import netmhcpan_api as api

SEED = 20261010
BINDING = ['log_EL_rank', 'log_BA_rank', 'log_affinity', 'length', 'is_Kb']
FEATURES = BINDING + ['WT_minus_mutant_log_affinity', 'WT_minus_mutant_log_EL_rank',
                      'hydrophobic_fraction', 'aromatic_fraction', 'charge_fraction',
                      'proline_fraction', 'mutation_position_fraction', 'hydrophobic_change']


def verify_lock():
    lock = json.loads((ROOT/'locked/split_lock.json').read_text())
    for name, expected in lock['files'].items():
        if sha(ROOT/'locked'/name) != expected:
            raise ValueError('Locked split changed: ' + name)
    for name, expected in lock['inputs'].items():
        if sha(ROOT.parent/name) != expected:
            raise ValueError('Input changed: ' + name)
    if sha(ROOT/'PROTOCOL.md') != lock['protocol_sha256'] or sha(ROOT/'prepare.py') != lock['prepare_sha256']:
        raise ValueError('Frozen protocol/curation code changed')
    return lock


def feature_frame(df):
    """Only sequence, mutation, allele and predicted binding enter the model."""
    hydrophobic = set('AVILMFWY')
    positions = [next(i for i,(a,b) in enumerate(zip(mut,wt)) if a!=b) for mut,wt in zip(df.Peptide,df.WT)]
    f = pd.DataFrame(index=df.index)
    f['log_EL_rank'] = np.log10(1 + df.Rank_EL)
    f['log_BA_rank'] = np.log10(1 + df.Rank_BA)
    f['log_affinity'] = np.log10(df.Affinity_nM.clip(1,50000))
    f['length'] = df.Peptide.str.len()
    f['is_Kb'] = (df.MHC == 'H-2-Kb').astype(int)
    f['WT_minus_mutant_log_affinity'] = np.log10(df.WT_Affinity_nM.clip(1,50000)) - f.log_affinity
    f['WT_minus_mutant_log_EL_rank'] = np.log10(1+df.WT_Rank_EL) - f.log_EL_rank
    f['hydrophobic_fraction'] = df.Peptide.map(lambda s: sum(a in hydrophobic for a in s)/len(s))
    f['aromatic_fraction'] = df.Peptide.map(lambda s: sum(a in 'FWY' for a in s)/len(s))
    f['charge_fraction'] = df.Peptide.map(lambda s: (sum(a in 'KR' for a in s)-sum(a in 'DE' for a in s))/len(s))
    f['proline_fraction'] = df.Peptide.map(lambda s: s.count('P')/len(s))
    f['mutation_position_fraction'] = [(p+1)/len(s) for p,s in zip(positions,df.Peptide)]
    f['hydrophobic_change'] = [int(m[p] in hydrophobic)-int(w[p] in hydrophobic) for m,w,p in zip(df.Peptide,df.WT,positions)]
    if list(f.columns) != FEATURES or not np.isfinite(f.to_numpy()).all():
        raise ValueError('Invalid features')
    return f


def models():
    return {
        'LR_binding': (make_pipeline(StandardScaler(), LogisticRegression(C=.1,max_iter=5000,random_state=SEED)), BINDING),
        'LR_tumor_primary': (make_pipeline(StandardScaler(), LogisticRegression(C=.1,max_iter=5000,random_state=SEED)), FEATURES),
        'GB_tumor_exploratory': (GradientBoostingClassifier(n_estimators=100,learning_rate=.05,max_depth=2,
                                                          min_samples_leaf=10,subsample=1,random_state=SEED), FEATURES),
    }


def predictions(cohort, out):
    rawdir = out/'netmhcpan_raw'
    rawdir.mkdir(exist_ok=True)
    seqs = sorted(set(cohort.Peptide) | set(cohort.WT))
    api.TMPDIR = rawdir/'tmp'
    api.TMPDIR.mkdir(exist_ok=True)
    runner = api.NetMHCpanRunner(timeout=1800)
    frames, commands = [], []
    for allele in ['H-2-Kb','H-2-Db']:
        p = rawdir/(allele+'.pep')
        p.write_text('\n'.join(seqs)+'\n')
        cmd = runner._build_cmd(str(p), '-p', [allele], True, False, False, False)
        print(f'Predicting {len(seqs)} mutant/WT sequences: {allele}',flush=True)
        raw = runner._run(cmd)
        (rawdir/(allele+'.txt')).write_text(raw)
        frame = api.parse_output(raw,include_ba=True)
        if len(frame)!=len(seqs) or set(frame.Peptide)!=set(seqs) or set(frame.MHC)!={allele}:
            raise ValueError('Incomplete netMHCpan predictions')
        frames.append(frame)
        commands.append(cmd)
    pred = pd.concat(frames,ignore_index=True)
    pred.to_csv(out/'binding_predictions.csv',index=False)
    result = cohort.merge(pred,on=['Peptide','MHC'],how='left',validate='many_to_one')
    wt = pred[['Peptide','MHC','Rank_EL','Rank_BA','Affinity_nM']].rename(columns={
        'Peptide':'WT','Rank_EL':'WT_Rank_EL','Rank_BA':'WT_Rank_BA','Affinity_nM':'WT_Affinity_nM'})
    result = result.merge(wt,on=['WT','MHC'],how='left',validate='many_to_one')
    if len(result)!=len(cohort) or result[['Rank_EL','WT_Rank_EL','WT_Affinity_nM']].isna().any().any():
        raise ValueError('Incomplete mutant/WT joins')
    return add_scores(result),commands


def train():
    lock = verify_lock()
    out = ROOT/'run'
    out.mkdir(exist_ok=True)
    if (out/'model_manifest.json').exists() or (out/'evaluation_started.json').exists():
        raise FileExistsError('This experiment already has fitted models; do not iterate against its test set')
    cohort = pd.concat([pd.read_csv(ROOT/'locked'/(s+'.csv'),dtype={'PMID':str,'Reference_ID':str})
                        for s in ['train','validation','test']],ignore_index=True)
    scored,commands = predictions(cohort,out)
    feats = feature_frame(scored)
    feats.to_csv(out/'features.csv',index=False)
    train_mask = scored.Split == 'train'
    train_df = scored[train_mask]
    if train_df.Label.nunique()!=2:
        raise ValueError('Training requires both classes')
    fitted = {}
    for name,(model,cols) in models().items():
        model.fit(feats.loc[train_mask,cols],train_df.Label)
        joblib.dump(model,out/(name+'.joblib'))
        fitted[name]=(model,cols)
    coef = fitted['LR_tumor_primary'][0].named_steps['logisticregression'].coef_[0]
    pd.DataFrame({'feature':FEATURES,'standardized_coefficient':coef}).to_csv(out/'primary_coefficients.csv',index=False)
    # Write the fitted checkpoint BEFORE calculating external scores.
    manifest = {
        'fitted_utc':datetime.now(timezone.utc).isoformat(), 'seed':SEED,
        'lock_sha256':sha(ROOT/'locked/split_lock.json'), 'train_code_sha256':sha(Path(__file__)),
        'protocol_sha256':lock['protocol_sha256'], 'training_rows':int(train_mask.sum()),
        'primary_model':'LR_tumor_primary','hyperparameter_search':False,'validation_used_for_fitting':False,
        'sklearn':sklearn.__version__,'features':FEATURES,'binding_features':BINDING,
        'model_files':{name+'.joblib':sha(out/(name+'.joblib')) for name in fitted},
        'commands':commands, 'netmhcpan_version':api.VERSION_FILE.read_text().strip(),
        'predictor_hashes':{str(p):sha(p) for p in [api.BINARY,api.SYN_NOCONTEXT,api.HLA_PSEUDO]},
    }
    (out/'model_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    # Fixed leave-one-cell-line-out development diagnostics; no tuning.
    cv = []
    for tumor in sorted(train_df.Tumor_Model.unique()):
        test_idx = train_df.index[train_df.Tumor_Model==tumor]
        candidates = train_df[train_df.Tumor_Model!=tumor]
        held = train_df.loc[test_idx].to_dict('records')
        fit_idx = [i for i,r in candidates.iterrows() if not any(overlap(r,h) for h in held)]
        for name,(model,cols) in models().items():
            if scored.loc[fit_idx,'Label'].nunique()!=2:
                continue
            model.fit(feats.loc[fit_idx,cols],scored.loc[fit_idx,'Label'])
            vals=model.predict_proba(feats.loc[test_idx,cols])[:,1]
            for i,v in zip(test_idx,vals):
                cv.append({'Peptide':scored.loc[i,'Peptide'],'MHC':scored.loc[i,'MHC'],
                           'Tumor_Model':tumor,'Label':int(scored.loc[i,'Label']),
                           'Model':name,'Score':float(v),'Fit_N':len(fit_idx)})
    pd.DataFrame(cv).to_csv(out/'development_out_of_fold.csv',index=False)
    for name,(model,cols) in fitted.items():
        scored[name]=model.predict_proba(feats[cols])[:,1]
    scored.to_csv(out/'scored_cohort.csv',index=False)
    print('Fixed models fitted; checkpoint saved before external evaluation.',flush=True)


if __name__=='__main__':
    train()
