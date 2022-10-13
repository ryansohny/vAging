from anndata import AnnData
import anndata
from scipy import sparse, io
import scipy
import pandas as pd
import scipy.io
import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors
matplotlib.use('TkAgg')
import numpy as np
import seaborn as sb
import seaborn as sns
import math
import scanpy.external as sce
import scrublet as scr
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
sns.set(font="Arial", font_scale=1, style='ticks')
sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (6,6)
plt.rcParams['figure.dpi'] = 100
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#000000", "#8b0a50"])
batch_palette=['#689aff', '#fdbf6f', '#b15928']
%matplotlib
%autoindent

# Open an EC scRNA-seq file
test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")

test3_endo.obs['Age'] = test3_endo.obs['batch']
test3_endo.uns['Age_colors'] = ['#689aff', '#fdbf6f', '#b15928']

endo_leiden_to_celltype_dict = {'0': 'EC1',
'1': 'EC4',
'2': 'EC2',
'3': 'EC3',
'4': 'EC5',
'5': 'EC6'}
test3_endo.obs['Subpopulation of Endothelial Cells'] = test3_endo.obs['endo_leiden_r05'].map(lambda x: endo_leiden_to_celltype_dict[x]).astype('category')

# Figure A

fig, axes = plt.subplots(1, 2, figsize=(12,6))
sns.despine()
sc.pl.umap(test3_endo, color=['Age'], add_outline=False, legend_loc='right margin', color_map=cmap, show=False, ax=axes[0])
axes[0].legend(frameon=False, loc='lower center', bbox_to_anchor=(0.45, 0))
sc.pl.umap(test3_endo, color=['Subpopulation of Endothelial Cells'], add_outline=False, legend_loc='on data', palette='Set3', show=False, ax=axes[1])
axes[0].text(-7.3, 8, "A", size=20, weight='bold')
axes[1].text(-7.3, 8, "B", size=20, weight='bold')
plt.savefig('./figures/Figure_A', dpi=600)

# Figure B

lin = ('EC1', 'EC2', 'EC3', 'EC4', 'EC5', 'EC6')
#test3_endo.obs['leiden_r05']
test3_endo.obs['Subpopulation of Endothelial Cells'] = test3_endo.obs['Subpopulation of Endothelial Cells'].cat.reorder_categories(list(lin), ordered=True)

uni_ec = {'Universal EC markers': ["Pecam1", "Cdh5", 'Erg', 'Vwf', "Nos3", 'Procr', 'Icam1', 'Cd34', 'Cd36', ]}
ax_dict = sc.pl.matrixplot(test3_endo, uni_ec, layer='magic', groupby='Subpopulation of Endothelial Cells', dendrogram=False, cmap='viridis', standard_scale='var', colorbar_title='Scaled\nexpression', var_group_rotation=0, show=False)
ax_dict['mainplot_ax'].set_xticklabels(labels=list(uni_ec.values())[0], fontstyle='italic', rotation=45)

# Diff version
uni_ec = {'Universal EC markers': ["Pecam1", "Cdh5", "Kdr", 'Erg', 'Vwf', "Nos3", 'Procr', 'Icam1', 'Cd34', 'Cd36', ]}
ax_dict = sc.pl.matrixplot(test3_endo, uni_ec, layer='magic', groupby='Subpopulation of Endothelial Cells', dendrogram=False, cmap='viridis', standard_scale='var', colorbar_title='Scaled\nexpression', var_group_rotation=0, show=False)
ax_dict['mainplot_ax'].set_xticklabels(labels=list(uni_ec.values())[0], fontstyle='italic', rotation=45)


ec_others = {'Pericyte': ['Rgs5', 'Cspg4', 'Kcnj8', 'Des'], 'Lymphatic EC': ['Reln', 'Flt4'], 'Vasa Vasorum': ['Ackr1', 'Lrg1']}
ax_dict2 = sc.pl.matrixplot(test3_endo, ec_others, layer='magic', groupby='Subpopulation of Endothelial Cells', dendrogram=False, cmap='viridis', standard_scale='var', colorbar_title='Scaled\nexpression', var_group_rotation=0, show=False)
ax_dict2['mainplot_ax'].set_xticklabels(labels=list(x for xs in list(ec_others.values()) for x in xs), fontstyle='italic', rotation=45)
