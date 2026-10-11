"""Export fixed-test comparisons; baseline prevalence shown for each cohort."""
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR', tempfile.mkdtemp(prefix='tumor_plot_'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from prepare import ROOT

df=pd.read_csv(ROOT/'run/external_metrics.csv')
names={'BA_Rank_Score':'netMHCpan BA rank','Mouse_Composite':'Current composite',
       'LR_binding':'Binding-only logistic','LR_tumor_primary':'Tumor logistic (primary)',
       'GB_tumor_exploratory':'Tumor gradient boosting'}
fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
for ax,cohort,title in zip(axes,['pooled','PMID_36950115','PMID_40586324'],
                         ['Pooled test','LLC1 / 2023','CMT93 / 2025']):
    part=df[(df.Split=='test') & (df.Cohort==cohort)].set_index('Model')
    values=[part.loc[m,'AP'] for m in names]
    bars=ax.barh(list(names.values()),values,color=['#377a9f','#7895a8','#729e79','#bc6d62','#aa8bb6'])
    ax.invert_yaxis()
    for bar,v in zip(bars,values):ax.text(v+.012,bar.get_y()+bar.get_height()/2,f'{v:.3f}',va='center',fontsize=9)
    prev=part.Prevalence.iloc[0]
    ax.axvline(prev,ls='--',color='#555555',lw=1,label=f'Prevalence {prev:.3f}')
    ax.set(xlim=(0,1.12),xlabel='Average precision (AP)',title=f'{title}\nN={int(part.N.iloc[0])}')
    ax.legend(loc='lower right',fontsize=8)
fig.suptitle('Tumor neoantigen feasibility: fixed research holdout, no production deployment',fontsize=12)
for ext in ['png','svg']:fig.savefig(ROOT/'run'/('test_comparison.'+ext),dpi=180)
plt.close(fig)
