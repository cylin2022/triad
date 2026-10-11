"""Read-only model/split verification, with an exported audit record."""
import importlib.metadata
import json
import platform
from datetime import datetime, timezone

import joblib
import numpy as np
import pandas as pd

from prepare import ROOT, overlap, sha
from train import BINDING, FEATURES, feature_frame, verify_lock
from assess import metrics


def audit():
    lock=verify_lock();out=ROOT/'run'
    manifest=json.loads((out/'model_manifest.json').read_text())
    started=json.loads((out/'evaluation_started.json').read_text())
    done=json.loads((out/'evaluation_completed.json').read_text())
    assert manifest['lock_sha256']==sha(ROOT/'locked/split_lock.json')
    assert manifest['train_code_sha256']==sha(ROOT/'train.py')
    assert started['assess_code_sha256']==sha(ROOT/'assess.py')
    assert started['model_manifest_sha256']==sha(out/'model_manifest.json')
    assert started['scored_cohort_sha256']==sha(out/'scored_cohort.csv')
    assert done['metrics_sha256']==sha(out/'external_metrics.csv')
    assert datetime.fromisoformat(lock['locked_utc']) < datetime.fromisoformat(manifest['fitted_utc']) < datetime.fromisoformat(started['opened_utc']) < datetime.fromisoformat(done['completed_utc'])
    df=pd.read_csv(out/'scored_cohort.csv',dtype={'PMID':str})
    f=feature_frame(df)
    original=pd.concat([pd.read_csv(ROOT/'locked'/(s+'.csv'),dtype={'PMID':str}) for s in ['train','validation','test']],ignore_index=True)
    keys=['Peptide','WT','MHC','Label','PMID','Tumor_Model','Split']
    pd.testing.assert_frame_equal(df[keys],original[keys])
    for name,h in manifest['model_files'].items():
        assert sha(out/name)==h
        model=joblib.load(out/name);cols=BINDING if name=='LR_binding.joblib' else FEATURES
        np.testing.assert_allclose(model.predict_proba(f[cols])[:,1],df[name.removesuffix('.joblib')],rtol=0,atol=1e-12)
        if name.startswith('LR_'):
            np.testing.assert_allclose(model.named_steps['standardscaler'].mean_,f.loc[df.Split=='train',cols].mean().to_numpy(),rtol=0,atol=1e-12)
    parts={s:df[df.Split==s] for s in ['train','validation','test']}
    for a,b in [('train','test'),('train','validation'),('validation','test')]:
        assert not any(overlap(x,y) for x in parts[a].to_dict('records') for y in parts[b].to_dict('records'))
    table=pd.read_csv(out/'external_metrics.csv')
    for row in table.itertuples():
        p=parts[row.Split]
        if row.Cohort.startswith('PMID_'):p=p[p.PMID==row.Cohort.removeprefix('PMID_')]
        elif row.Cohort!='pooled':p=p[p.MHC==row.Cohort]
        vals=metrics(p,row.Model)
        for k,v in vals.items():
            if v is None:assert pd.isna(getattr(row,k))
            else:assert np.isclose(getattr(row,k),v,atol=1e-12)
    report={
        'audited_utc':datetime.now(timezone.utc).isoformat(), 'rows':len(df),
        'split_and_source_hashes_valid':True,'lock_before_fit_before_test_evaluation':True,
        'serialized_predictions_match':True,'scalers_fit_training_only':True,
        'external_metrics_recomputed_without_refitting':True,'cross_split_sequence_source_overlap':0,
        'python':platform.python_version(),
        'packages':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scikit-learn','scipy','joblib','openpyxl','matplotlib']},
        'code_hashes':{p.name:sha(p) for p in ROOT.glob('*.py')},
    }
    (out/'artifact_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':audit()
