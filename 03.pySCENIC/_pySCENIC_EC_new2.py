# _pySCENIC_EC_new2.py
# 2022-05-12

## pySCENIC's AUC matrix retrieval
# Environment (pySCENIC)
# ipython --profile=pyscenic

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
    ax.plot(np.arange(len(data)), data, '.', color='#104e8b')
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


# CELL TYPE SPECIFIC REGULATORS - Z-SCORE

signature_column_names = list(test3_endo.obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda x: x.startswith('Regulon('), signature_column_names))
test3_endo_scores = test3_endo.obs[signature_column_names + ['EC_subclusters']]
test3_endo_results = ((test3_endo_scores.groupby('EC_subclusters').mean() - test3_endo.obs[signature_column_names].mean()) / test3_endo.obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0: 'Z'})
test3_endo_results['regulon'] = list(map(lambda x: x[8:-1], test3_endo_results.regulon))

signature_column_names = list(test3_endo[~test3_endo.obs['EC_subclusters'].isin(['EC5', 'EC6'])].obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda x: x.startswith('Regulon('), signature_column_names))
test3_endo_scores = test3_endo[~test3_endo.obs['EC_subclusters'].isin(['EC5', 'EC6'])].obs[signature_column_names + ['EC_subclusters']]
test3_endo_results = ((test3_endo_scores.groupby('EC_subclusters').mean() - test3_endo[~test3_endo.obs['EC_subclusters'].isin(['EC5', 'EC6'])].obs[signature_column_names].mean()) / test3_endo[~test3_endo.obs['EC_subclusters'].isin(['EC5', 'EC6'])].obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0: 'Z'})

# 위를 Z-score가 아니라 일반 scale로 하는 것도 방법

test3_endo_results['regulon'] = list(map(lambda x: x[8:-1], test3_endo_results.regulon))


test3_endo_heatmap = pd.pivot_table(data=test3_endo_results[(test3_endo_results.Z >= 1.25)].sort_values('Z', ascending=False), index='EC_subclusters', columns='regulon', values='Z')

ec1_regulon_list = list(test3_endo_results[test3_endo_results['EC_subclusters'] == 'EC_1'].sort_values('Z', ascending=False).iloc[:10]['regulon'].values)
ec2_regulon_list = list(test3_endo_results[test3_endo_results['EC_subclusters'] == 'EC_2'].sort_values('Z', ascending=False).iloc[:10]['regulon'].values)
ec3_regulon_list = list(test3_endo_results[test3_endo_results['EC_subclusters'] == 'EC_3'].sort_values('Z', ascending=False).iloc[:10]['regulon'].values)
ec4_regulon_list = list(test3_endo_results[test3_endo_results['EC_subclusters'] == 'EC_4'].sort_values('Z', ascending=False).iloc[:10]['regulon'].values)

test3_endo_heatmap = pd.pivot_table(data=test3_endo_results[test3_endo_results['regulon'].isin(ec1_regulon_list + ec2_regulon_list + ec3_regulon_list + ec4_regulon_list)].sort_values('Z', ascending=False), index='EC_subclusters', columns='regulon', values='Z')

import matplotlib
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])

sns.clustermap(test3_endo_heatmap.sort_index(ascending=False),
               row_cluster=False,
               method='ward',
               metric='euclidean', z_score=None, standard_scale=None, cmap=cmap, xticklabels=True, yticklabels=True)
  


# Regulon to files
import pandas as pd
import loompy as lp
adjacencies = pd.read_csv("test3_endo_adj.csv", index_col=False)
lf = lp.connect(
    "/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/test3_endo_pyscenic_output.loom",
        mode='r+', validate=False)
exprMat = pd.DataFrame(lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = dict()
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)

from pyscenic.utils import modules_from_adjacencies
modules = list(modules_from_adjacencies(adjacencies, exprMat))

tf = 'Isl1'
tf_mods = [x for x in modules if x.transcription_factor == tf]
for i, mod in enumerate(tf_mods):
    print( f'{tf} module {str(i)}: {len(mod.genes)} genes' )

print( f'{tf} regulon: {len(regulons[tf+"(+)"])} genes' )

for i,mod in enumerate( tf_mods ):
    with open( tf+'_module_'+str(i)+'.txt', 'w') as f:
        for item in mod.genes:
            f.write("%s\n" % item)
            
with open( tf+'_regulon.txt', 'w') as f:
    
    for item in regulons[tf+'(+)']:
        f.write("%s\n" % item)

"""
type(modules) == list
len(modules) == 7,580

modules는 ctxcore.genesig.Regulon이라고하는 type들의 list로 구성되어 있다
ex) type(modules[0]) == ctxcore.genesig.Regulon

각각은 다음과 같이 구성되어 있다.
ex) modules[601]
Regulon(name='Regulon for Nfia', gene2weight=frozendict.frozendict({'Pam': 41.442711481206985, 'Gxylt2': 34.37402163050602, 'Nfix': 32.006428302211745,......)

이들 각각의 element들은 다음과 같이 접근이 가능하다.

modules[601].name == 'Regulon for Nfia'

각각의 element의 이름들은 아래와 같다 (각각이 의미하는 것들 나중에 적어놓을 것)
name
gene2weight
gene2occurrence
transcription_factor
context
score
nes
orthologous_identity
similarity_qvalue
annotation

"""

a, b, c, d, e, f = [], [], [], [], [], []