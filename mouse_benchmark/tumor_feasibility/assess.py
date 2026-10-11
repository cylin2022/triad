"""One locked external evaluation; no model selection or retraining."""
import json
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from prepare import ROOT, overlap, sha
from train import verify_lock
sys.path.insert(0,str(ROOT.parent))
from evaluate import top_k

MODELS=['Calis','EL_Rank_Score','BA_Rank_Score','Mouse_Base','Mouse_Composite',
        'LR_binding','LR_tumor_primary','GB_tumor_exploratory']


def metrics(part,model):
    y=part.Label.to_numpy();s=part[model].to_numpy()
    both=len(set(y))==2
    result={'N':len(part),'Positive':int(y.sum()),'Prevalence':float(y.mean()) if len(y) else None,
            'AP':float(average_precision_score(y,s)) if both else None,
            'ROC_AUC':float(roc_auc_score(y,s)) if both else None}
    for k in [5,10,20]:
        result['Precision_at_'+str(k)],result['Recall_at_'+str(k)]=top_k(y,s,k)
    return result


def similarity_groups(part):
    rows=part.to_dict('records');parent=list(range(len(rows)))
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i in range(len(rows)):
        for j in range(i):
            if overlap(rows[i],rows[j]):parent[root(i)]=root(j)
    return np.unique([root(i) for i in range(len(rows))],return_inverse=True)[1]


def uncertainty(part):
    # Stratified by study. No claim of new-study sampling coverage.
    group_indices=[]
    for _,sub in part.groupby('PMID'):
        groups=similarity_groups(sub)
        group_indices.append([part.index.get_indexer(sub.index[groups==g]) for g in np.unique(groups)])
    rng=np.random.default_rng(20261010)
    deltas={m:[] for m in ['Mouse_Base','BA_Rank_Score','Mouse_Composite','LR_binding']}
    y=part.Label.to_numpy()
    for _ in range(1000):
        w=np.zeros(len(part))
        for blocks in group_indices:
            for g in rng.integers(0,len(blocks),size=len(blocks)):w[blocks[g]]+=1
        if len(set(y[w>0]))!=2:continue
        ap=average_precision_score(y,part.LR_tumor_primary,sample_weight=w)
        for m in deltas:deltas[m].append(float(ap-average_precision_score(y,part[m],sample_weight=w)))
    return {'replicates':1000,'valid_replicates':len(deltas['Mouse_Base']),
            'groups_per_study':[len(x) for x in group_indices],
            'method':'paired sequence/source-group bootstrap, within fixed studies',
            'scope':'conditional candidate sampling only; not new-study uncertainty',
            'primary_minus_comparator_AP_95':{m:np.quantile(v,[.025,.975]).tolist() if v else None for m,v in deltas.items()}}


def assess():
    verify_lock()
    out=ROOT/'run'
    marker=out/'evaluation_started.json'
    if marker.exists():raise FileExistsError('External evaluation already opened; do not tune against it')
    manifest=json.loads((out/'model_manifest.json').read_text())
    if manifest['lock_sha256']!=sha(ROOT/'locked/split_lock.json') or manifest['train_code_sha256']!=sha(ROOT/'train.py'):
        raise ValueError('Training checkpoint provenance mismatch')
    for name,h in manifest['model_files'].items():
        if sha(out/name)!=h:raise ValueError('Model artifact changed')
    df=pd.read_csv(out/'scored_cohort.csv',dtype={'PMID':str})
    if not np.isfinite(df[MODELS].to_numpy()).all():raise ValueError('Missing prediction')
    marker.write_text(json.dumps({'opened_utc':datetime.now(timezone.utc).isoformat(),
                                 'model_manifest_sha256':sha(out/'model_manifest.json'),
                                 'scored_cohort_sha256':sha(out/'scored_cohort.csv'),
                                 'assess_code_sha256':sha(ROOT/'assess.py')},indent=2)+'\n')
    results=[]
    for split in ['validation','test']:
        part=df[df.Split==split]
        subsets={'pooled':part}
        subsets.update({a:p for a,p in part.groupby('MHC')})
        subsets.update({'PMID_'+str(pmid):p for pmid,p in part.groupby('PMID')})
        for name,p in subsets.items():
            for model in MODELS:results.append({'Split':split,'Cohort':name,'Model':model,**metrics(p,model)})
    table=pd.DataFrame(results);table.to_csv(out/'external_metrics.csv',index=False)
    # No training-set resubstitution accuracy is presented as validation.
    cv=pd.read_csv(out/'development_out_of_fold.csv')
    cv_metrics=[]
    for (tumor,model),p in cv.groupby(['Tumor_Model','Model']):
        cv_metrics.append({'Tumor_Model':tumor,'Model':model,**metrics(p.rename(columns={'Score':model}),model)})
    pd.DataFrame(cv_metrics).to_csv(out/'development_metrics.csv',index=False)
    test=df[df.Split=='test'].copy()
    intervals={'test_pooled':uncertainty(test)}
    for pmid,p in test.groupby('PMID'):intervals['PMID_'+pmid]=uncertainty(p)
    (out/'uncertainty.json').write_text(json.dumps(intervals,indent=2)+'\n')
    print(table[(table.Split=='test') & table.Cohort.isin(['pooled','PMID_36950115','PMID_40586324'])].to_string(index=False))
    print(json.dumps(intervals,indent=2))
    (out/'evaluation_completed.json').write_text(json.dumps({'completed_utc':datetime.now(timezone.utc).isoformat(),
                                                            'metrics_sha256':sha(out/'external_metrics.csv')},indent=2)+'\n')


if __name__=='__main__':assess()
