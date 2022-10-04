#export PATH=/data/Projects/phenomata/99.Tools/anaconda3/bin:$PATH
#source activate scanpy_1.9.1
#ipython --profile=vaging (in cm03)

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
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'Arial'
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
batch_palette=['#689aff', '#fdbf6f', '#b15928']
%matplotlib
%autoindent

m01 = sc.read_10x_h5("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/01.Cell-Ranger/01month_filtered_feature_bc_matrix.h5")
m01.var_names_make_unique()

m10 = sc.read_10x_h5("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/01.Cell-Ranger/10months_filtered_feature_bc_matrix.h5")
m10.var_names_make_unique()

m20 = sc.read_10x_h5("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/01.Cell-Ranger/20months_filtered_feature_bc_matrix.h5")
m20.var_names_make_unique()

mito_genes = m01.var_names.str.startswith('mt-')
m01.obs['percent_mito'] = np.ravel(np.sum(m01[:, mito_genes].X, axis=1)) / np.ravel(np.sum(m01.X, axis=1))
m10.obs['percent_mito'] = np.ravel(np.sum(m10[:, mito_genes].X, axis=1)) / np.ravel(np.sum(m10.X, axis=1))
m20.obs['percent_mito'] = np.ravel(np.sum(m20[:, mito_genes].X, axis=1)) / np.ravel(np.sum(m20.X, axis=1))

for sample in [m01, m10, m20]:
    sce.pp.scrublet(sample, adata_sim=None, sim_doublet_ratio=2.0, expected_doublet_rate=0.05, stdev_doublet_rate=0.02, synthetic_doublet_umi_subsampling=1.0, knn_dist_metric='euclidean', n_prin_comps=30, verbose=True)

sc.pp.filter_cells(m01, min_counts=2000)
sc.pp.filter_cells(m01, min_genes=1500)

sc.pp.filter_cells(m10, min_counts=3000)
sc.pp.filter_cells(m10, min_genes=1500)

sc.pp.filter_cells(m20, min_counts=3000)
sc.pp.filter_cells(m20, min_genes=1500)

m01 = m01[m01.obs['percent_mito'] < 0.2]
m10 = m10[m10.obs['percent_mito'] < 0.2]
m20 = m20[m20.obs['percent_mito'] < 0.2]

integrated = AnnData.concatenate(m01, m10, m20, join='outer', batch_categories = ['m01', 'm10', 'm20'], index_unique = '-')
integrated.obs['Doublet'] = integrated.obs['predicted_doublet'].astype(str).astype('category')
integrated.obs[['Doublet', 'batch']].value_counts()
del integrated.obs['predicted_doublet']

sc.pp.filter_genes(integrated, min_cells=5) # 'n_cells' added in integrated.var 
integrated.layers["counts"] = integrated.X.copy()
integrated.raw = integrated

import rpy2.rinterface_lib.callbacks
import logging
from rpy2.robjects import pandas2ri
import anndata2ri
pandas2ri.activate()
anndata2ri.activate()
%load_ext rpy2.ipython
%%R
library(scran)
library(dplyr)



%config InlineBackend.figure_format = 'retina'

adata_pp = integrated.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp) # works on anndata.X
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = integrated.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))




del adata_pp
del data_mat

integrated.obs['size_factors'] = size_factors

integrated.X /= integrated.obs['size_factors'].values[:, None]
integrated.layers['scran'] = integrated.X # For cellphoneDB or CelChat maybe?
sc.pp.log1p(integrated) # works on anndata.X
integrated.X = scipy.sparse.csr_matrix(integrated.X)
integrated.raw = integrated ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

test3 = integrated.copy()
test3.raw = test3
test3.layers['scran_log1p'] = test3.X

#del integrated

sc.pp.highly_variable_genes(test3)
test3.var['highly_variable'].value_counts() # 2,410 ==> 2021-08-10 # 2,513 ==> 2022-09-26

sc.pp.scale(test3, max_value=10) # tabula muris senis default (2021-08-10) # mean and std on adata.var
#sc.pp.scale(test3, zero_center=True, max_value=10, copy=False, layer=None, obsm=None)

cell_cycle_genes=[x.strip()[0] + x.strip()[1:].lower() for x in open("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/regev_lab_cell_cycle_genes.txt")]
s_genes= cell_cycle_genes[:43]
g2m_genes= cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in test3.var_names]
sc.tl.score_genes_cell_cycle(test3, s_genes=s_genes, g2m_genes=g2m_genes)
"""
Used 'raw' attribute of adata (use_raw = True if .raw is present)
So, log-tranformed scran-normalized counts are put into score_genes_cell_cycle function
"""
sc.tl.pca(test3, n_comps=100, use_highly_variable=True, svd_solver='arpack')

sc.pl.pca_variance_ratio(test3, n_pcs=100, log=False)
sc.pl.pca(test3, color=['batch'], legend_loc='right margin', size=8, add_outline=False, color_map='CMRmap', components=['1,2'])

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3, batch_key='batch', n_pcs=20, neighbors_within_batch=5, trim=None)
sc.tl.umap(test3, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
test3.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']
sc.pl.umap(test3, color=['batch'], add_outline=False, legend_loc='right margin', size=20)

sc.tl.leiden(test3, resolution=0.5, key_added='leiden_r05') #### 0 ~ 13 ==> 2021-09-28
sc.tl.leiden(test3, resolution=1.0, key_added='leiden_r10')
sc.pl.umap(test3, color=['batch', 'leiden_r05', 'leiden_r10'], add_outline=False, legend_loc='right margin', size=20)

fig, axes = plt.subplots(1,3)
sc.pl.umap(test3, color=['batch'], add_outline=False, legend_loc='right margin', size=20, groups=['m01'], title='1 month', ax=axes[0])
sc.pl.umap(test3, color=['batch'], add_outline=False, legend_loc='right margin', size=20, groups=['m10'], title='10 months', ax=axes[1])
sc.pl.umap(test3, color=['batch'], add_outline=False, legend_loc='right margin', size=20, groups=['m20'], title='20 months', ax=axes[2])

sc.tl.rank_genes_groups(test3, 'leiden_r05', method='wilcoxon', corr_method='benjamini-hochberg', use_raw=True, pts=True) # key_added=''
sc.pl.rank_genes_groups(test3, n_genes=5, sharey=False)
sc.pl.rank_

########################################################################################
# ComBat batch correction test (2022-09-26)
test3_combat = test3.copy()
test3_combat.X = test3_combat.layers['scran_log1p'].copy()
sc.pp.combat(test3_combat, key='batch')
"""
sc.pp.combat(..., covariates=None)
위에서처럼, correction되는 것을 원치 않는 (e.g. biological signal) 것을 covariates에 넣어야 하는데 안 넣음
"""
sc.pp.highly_variable_genes(test3_combat)
sc.tl.pca(test3_combat, n_comps=100, use_highly_variable=True, svd_solver='arpack')
#sc.pl.pca_variance_ratio(test3_combat, n_pcs=100, log=False)
sc.pp.neighbors(test3_combat, n_neighbors=15, n_pcs=20, knn=True, method='umap', metric='euclidean')
sc.tl.umap(test3_combat, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
test3_combat.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']
sc.pl.umap(test3_combat, color=['batch'], add_outline=False, legend_loc='right margin', size=20)

########################################################################################




test3 = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3.h5ad")
test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")

