# _pySCENIC_EC_new2.py
# 2022-05-12

## pySCENIC's AUC matrix retrieval
# Environment (scanpy_1.8.1)
import scanpy as sc
import loompy as lp
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import json
import zlib
import base64
from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.rss import regulon_specificity_scores
from pyscenic.cli.utils import load_signatures
#from pyscenic.plotting import plot_rss # 이거 그대로  function 쓰면 오류 생김 따라서 def() 로 함.
from pyscenic.binarization import binarize
from adjustText import adjust_text

%autoindent
%matplotlib
sns.set(font="Arial", font_scale=1, style='ticks')

lf = lp.connect(
    "/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/test3_endo_pyscenic_output.loom",
    mode='r+', validate=False)
lf.ca.keys()
# ['CellID', 'RegulonsAUC', 'nGene', 'nUMI']
lf.ra.keys()
# ['Gene', 'Regulons']
lf.attrs.keys()
# ['CreationDate', 'LOOM_SPEC_VERSION', 'MetaData', 'last_modified']

auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
auc_mtx.columns = auc_mtx.columns.str.replace('\(', '_(')

test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")
test3_endo = test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '1', '2', '3'])]

sig = load_signatures(
    '/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/test3_endo_reg.csv')
test3_endo = add_scenic_metadata(test3_endo, auc_mtx, sig)  # AUCell score가 test3_endo에 추가된다.

endo_leiden_to_celltype_dict = {'0': 'EC1',
'1': 'EC4',
'2': 'EC2',
'3': 'EC3'}

test3_endo.obs['EC_subclusters'] = test3_endo.obs['endo_leiden_r05'].map(lambda x: endo_leiden_to_celltype_dict[x]).astype('category')

reordered = ('EC1', 'EC2', 'EC3', 'EC4')
test3_endo.obs['EC_subclusters'] = test3_endo.obs['EC_subclusters'].cat.reorder_categories(list(reordered), ordered=True)

cellAnnot = test3_endo.obs[['batch', 'EC_subclusters']]
rss_EC_subclusters = regulon_specificity_scores(auc_mtx, cellAnnot['EC_subclusters'])

cats = sorted(list(set(cellAnnot['EC_subclusters'])))

data = rss_EC_subclusters.T['EC1'].sort_values(ascending=False)[0:rss_EC_subclusters.shape[1]]

# https://github.com/aertslab/pySCENIC/blob/master/src/pyscenic/plotting.py
from math import ceil, floor
def plot_rss(rss, cell_type, top_n=5, max_n=None, ax=None):
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(4, 4))
    if max_n is None:
        max_n = rss.shape[1]
    data = rss.T[cell_type].sort_values(ascending=False)[0:max_n]
    ax.plot(np.arange(len(data)), data, '.', color='darkblue')
    ax.set_ylim([floor(data.min() * 100.0) / 100.0, ceil(data.max() * 100.0) / 100.0])
    ax.set_ylabel('RSS')
    ax.set_xlabel('Regulon')
    ax.set_title(cell_type)
    ax.set_xticklabels([])
    font = {
        'color': 'darkred',
        'weight': 'normal',
        'size': 4,
    }
    for idx, (regulon_name, rss_val) in enumerate(zip(data[0:top_n].index, data[0:top_n].values)):
        ax.plot([idx, idx], [rss_val, rss_val], 'r.')
        ax.text(
            idx + (max_n / 25),
            rss_val,
            regulon_name,
            fontdict=font,
            horizontalalignment='left',
            verticalalignment='center',
        )



fig = plt.figure(figsize=(15, 8))

for c,num in zip(cats, range(1,len(cats)+1)):
    x = rss_EC_subclusters.T[c]
    ax = fig.add_subplot(1,4,num)
    plot_rss(rss_EC_subclusters, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )


sns.despine()

