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
sns.set(font_scale=0.75, style='ticks')
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'Arial'
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
batch_palette=['#689aff', '#fdbf6f', '#b15928']
%matplotlib
%autoindent

test3 = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3.h5ad")
test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")

endo_leiden_to_celltype_dict = {'0': 'EC_1',
'1': 'EC_4',
'2': 'EC_2',
'3': 'EC_3',
'4': 'EC_5',
'5': 'EC_6'}
test3_endo.obs['Subpopulation of Endothelial Cells'] = test3_endo.obs['endo_leiden_r05'].map(lambda x: endo_leiden_to_celltype_dict[x]).astype('category')
sc.pl.umap(test3_endo, color=['Subpopulation of Endothelial Cells'], add_outline=False, legend_loc='right margin', color_map=cmap, palette='Accent')
test3_endo2 = test3_endo[~test3_endo.obs['Subpopulation of Endothelial Cells'].isin(['EC_5', 'EC_6'])].copy()
colormap = {'EC_1': '#8dd3c7',
            'EC_2': '#80b1d3',
            'EC_3': '#fccde5',
            'EC_4': '#bebada'}
sc.pl.umap(test3_endo2, color='Subpopulation of Endothelial Cells', palette=colormap) # update colormap
lin = ('EC_1', 'EC_2', 'EC_3', 'EC_4')
test3_endo2.obs['Subpopulation of Endothelial Cells'] = test3_endo2.obs['Subpopulation of Endothelial Cells'].cat.reorder_categories(list(lin), ordered=True)


def umap_all(LIST):
    if type(LIST) == list:sc.pl.umap(test3, layer='magic', color=LIST, add_outline=False, legend_loc='right margin', color_map=cmap, ncols=4)
    elif type(LIST) == str:sc.pl.umap(test3, layer='magic', color=LIST.split(), add_outline=False, legend_loc='right margin', color_map=cmap, ncols=4)


def umap_endo(LIST):
    if type(LIST) == list:sc.pl.umap(test3_endo, layer='magic', color=LIST, add_outline=False, legend_loc='right margin', color_map=cmap, ncols=4)
    elif type(LIST) == str:sc.pl.umap(test3_endo, layer='magic', color=LIST.split(), add_outline=False, legend_loc='right margin', color_map=cmap, ncols=4)




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

sce.pl.scrublet_score_distribution(m01)

integrated = AnnData.concatenate(m01, m10, m20, join='outer', batch_categories = ['m01', 'm10', 'm20'], index_unique = '-')
integrated.obs['Doublet'] = integrated.obs['predicted_doublet'].astype(str).astype('category')

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

sc.pp.filter_genes(integrated, min_cells=5) # integrated.var에 n_cells 추가
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
integrated.obs['size_factors'] = size_factors

integrated.X /= integrated.obs['size_factors'].values[:, None]
integrated.layers['scran'] = integrated.X # For cellphoneDB
sc.pp.log1p(integrated) # works on anndata.X
integrated.X = scipy.sparse.csr_matrix(integrated.X)
integrated.raw = integrated ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

test3 = integrated.copy()
test3.raw = test3
test3.layers['scran_log1p'] = test3.X

sc.pp.highly_variable_genes(test3)
test3.var['highly_variable'].value_counts() # 2,410 ==> 2021-08-10 # 2,513 ==>

sc.pp.scale(test3, max_value=10) # tabula muris senis default (2021-08-10) # mean and std on adata.var
#sc.pp.scale(test3, zero_center=True, max_value=10, copy=False, layer=None, obsm=None)
cell_cycle_genes=[x.strip()[0] + x.strip()[1:].lower() for x in open("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/regev_lab_cell_cycle_genes.txt")]
s_genes= cell_cycle_genes[:43]
g2m_genes= cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in test3.var_names]
sc.tl.score_genes_cell_cycle(test3, s_genes=s_genes, g2m_genes=g2m_genes)

sc.tl.pca(test3, n_comps=100, use_highly_variable=True, svd_solver='arpack')
matplotlib.use('TkAgg')
%matplotlib
sc.pl.pca_variance_ratio(test3, n_pcs=100, log=False)
#sc.pl.pca(test3, color=['batch'], legend_loc='right margin', size=8, add_outline=False, color_map='CMRmap', components=['1,2'])

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3, batch_key='batch', n_pcs=20, neighbors_within_batch=20, trim=None) #####
sc.tl.umap(test3, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
#test3.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
test3.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']
#sc.pl.umap(test3, color=['batch'], add_outline=False, legend_loc='right margin', size=20, color_map='CMRmap')

sc.tl.leiden(test3, resolution=0.5, key_added='leiden_r05') #### 0 ~ 13 ==> 2021-09-28
sc.tl.leiden(test3, resolution=1.0, key_added='leiden_r10')
sc.pl.umap(test3, color=['batch', 'leiden_r05', 'leiden_r10'], add_outline=False, legend_loc='right margin', size=20, color_map='CMRmap')

sc.tl.rank_genes_groups(test3, 'leiden_r05', method='wilcoxon', corr_method='benjamini-hochberg', use_raw=True, pts=True) # key_added=''
sc.pl.rank_genes_groups(test3, n_genes=5, sharey=False)
sc.pl.rank_genes_groups_heatmap(test3, n_genes=10, min_logfoldchange=2, cmap='cividis', show_gene_labels=True)

markers = ["Pecam1", "Cdh5", "Nos3", "Acta2", "Cnn1", "Tagln", "Rgs5", "Kcnj8", "Col1a1", "Col5a1", "Dpt", "Cd19", "Ighm", "Cd14", "Cd68", "Cd3d"] # Cd3g 없음
sc.pl.stacked_violin(test3, markers, groupby='batch')

result = test3.uns['rank_genes_groups']
groups = result['names'].dtype.names
deg_wilcoxon = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals_adj']})
deg_wilcoxon.to_csv("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/20210916_scanpy_deg.csv", mode='w')

# leiden_r0.5 ==> 20210916_scanpy_deg.csv ==> Fold change cutoff 1.0 / padj < 0.05 ==> canonical markers (2021-11-16)
leiden_to_celltype_dict = {'0': 'Smooth muscle cells',
'1': 'Smooth muscle cells',
'2': 'Smooth muscle cells',
'3': 'Fibroblasts',
'4': 'Smooth muscle cells',
'5': 'Endothelial cells',
'6': 'Fibroblasts',
'7': 'Endothelial cells',
'8': 'Smooth muscle cells',
'9': 'Fibroblasts',
'10': 'B cells',
'11': 'M\u03A6',
'12': 'T cells',
'13': 'Fibroblasts'}
test3.obs['celltype'] = test3.obs['leiden_r05'].map(lambda x: leiden_to_celltype_dict[x]).astype('category')

sc.pl.umap(test3, color=['batch', 'leiden_r05', 'celltype'], add_outline=False, legend_loc='right margin', size=50, color_map='CMRmap')

celltype_marker = {'EC': ['Pecam1', 'Cdh5', 'Vwf', 'Nos3'],
'SMC': ['Acta2', 'Tagln', 'Cnn1', 'Cnn2'],
'FB': ['Dpt', 'Col1a1', 'Col5a1', 'Pdgfra'],
'B cells': ['Ighm', 'Cd19'],
'MΦ':['Cd14', 'Cd68'],
'T cells':['Cd3d', 'Cd3g']}
reordered = ('Endothelial cells', 'Smooth muscle cells', 'Fibroblasts', 'B cells', 'M\u03A6', 'T cells')
#reordered = ('Endothelial cells', 'Smooth muscle cells', 'Fibroblasts', 'B cells', 'Macrophages', 'T cells')
test3.obs['celltype'] = test3.obs['celltype'].cat.reorder_categories(list(reordered), ordered=True)

sc.pl.umap(test3, color=['celltype'], add_outline=False, legend_loc='right margin', size=30, color_map=cmap, palette='Paired')

leiden_to_celltype2_dict = {'0': 'VSMC_1',
'1': 'VSMC_2',
'2': 'VSMC_3',
'3': 'FB_1',
'4': 'VSMC_4',
'5': 'EC_1',
'6': 'FB_2',
'7': 'EC_2',
'8': 'VSMC_5',
'9': 'FB_3',
'10': 'B-lympho',
'11': 'M\u03A6',
'12': 'T-lmpho',
'13': 'Neuronal'}
test3.obs['celltype2'] = test3.obs['leiden_r05'].map(lambda x: leiden_to_celltype2_dict[x]).astype('category')
sc.pl.umap(test3, color=['celltype2'], add_outline=False, legend_loc='on data', size=30, color_map=cmap, palette='tab20')

colormap = dict()
c = 0
for i in leiden_to_celltype2_dict.values():
    colormap[i] = test3.uns['celltype2_colors'][c]
    c += 1

df = pd.concat([test3.obs['batch'], test3.obs['celltype2']], axis=1)
ax = pd.crosstab(df['batch'], df['celltype2'], normalize=0).plot.bar(stacked=True, color=colormap)
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

# Imputed expression matrix using MAGIC
import magic
test3_MAGIC = test3.copy()
test3_MAGIC.X = test3.layers['scran_log1p']
test3_MAGIC = magic.MAGIC().fit_transform(test3_MAGIC)

test3.layers['magic'] = test3_MAGIC.X
sc.pl.umap(test3, layer='magic', color=['Pecam1'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')
sc.pl.umap(test3, color=['Pecam1'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')


dp = sc.pl.dotplot(test3, layer='magic', celltype_marker, groupby='celltype', return_fig=True)
dp.add_totals(size=1.5, color=['#a6cee3', '#b2df8a', '#fb9a99', '#ff7f00', '#6a3d9a', '#b15928']).legend(colorbar_title='log(SizeFactorNormlized+1)', width=1.5, show_size_legend=False, show_colorbar=False).style(cmap='winter', dot_edge_color='black', dot_edge_lw=1, size_exponent=1.5, grid=True, x_padding=0.4, y_padding=0.6).swap_axes().show()
# MAGIC imputed expression
dp = sc.pl.dotplot(test3, celltype_marker, layer='magic', groupby='celltype', return_fig=True)
dp.add_totals(size=1.5, color=['#a6cee3', '#b2df8a', '#fb9a99', '#ff7f00', '#6a3d9a', '#b15928']).legend(colorbar_title='log(SizeFactorNormlized+1)', width=1.5, show_size_legend=False, show_colorbar=False).style(cmap='winter', dot_edge_color='black', dot_edge_lw=1, size_exponent=1.5, grid=True, x_padding=0.4, y_padding=0.6).swap_axes().show()
sc.pl.StackedViolin(test3, celltype_marker, layer='magic', groupby='celltype')
#test3.write(filename="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3.h5ad")

sc.pl.matrixplot(test3, celltype_marker, layer='magic', groupby='leiden_r05', dendrogram=True, cmap=cmap, standard_scale='var', colorbar_title='column scaled\nexpression')
sc.pl.matrixplot(test3, celltype_marker, layer='magic', groupby='celltype', dendrogram=True, cmap=cmap, standard_scale='var', colorbar_title='column scaled\nexpression')

# Using Different celltype markers
celltype_marker2 = {'EC': ['Pecam1', 'Cdh5', 'Vwf', 'Nos3'],
'SMC': ['Acta2', 'Tagln', 'Cnn1', 'Cnn2'],
'FB': ['Dpt', 'Lum', 'Col1a1', 'Pdgfra'],
'B cells': ['Ighm', 'Cd19'],
'MΦ':['Cd14', 'Cd68'],
'T cells':['Cd3d', 'Cd3g']}
sc.pl.matrixplot(test3, celltype_marker2, layer='magic', groupby='celltype', dendrogram=False, cmap=cmap, standard_scale='var', colorbar_title='Scaled\nexpression', swap_axes=True)

# Only Endothelial cells
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

test3_endo = anndata.AnnData(X=test3[test3.obs['leiden_r05'].isin(['5', '7'])].layers['counts'], obs=test3[test3.obs['leiden_r05'].isin(['5', '7'])].obs, var=test3[test3.obs['leiden_r05'].isin(['5', '7'])].var)
test3_endo.layers["counts"] = test3_endo.X.copy()

# Doublet information
test3_endo.obs['Doublet'] = integrated.obs['Doublet'].loc[test3_endo.obs.index]

# Doublet removal
#test3_endo = test3_endo[test3_endo.obs['Doublet'] == 'False']

adata_pp = test3_endo.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = test3_endo.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))



del adata_pp
test3_endo.obs['size_factors'] = size_factors

test3_endo.X /= test3_endo.obs['size_factors'].values[:, None]
test3_endo.X = scipy.sparse.csr_matrix(test3_endo.X) #왜 이게 새로 들어가야될까????? # 아니면 ERRROR 남 (highly_variable_genes에서)

test3_endo.layers['scran'] = test3_endo.X

sc.pp.log1p(test3_endo) # works on anndata.X
#integrated.X = scipy.sparse.csr_matrix(integrated.X)
test3_endo.layers['scran_log1p'] = test3_endo.X

test3_endo.raw = test3_endo ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

sc.pp.highly_variable_genes(test3_endo)
test3_endo.var['highly_variable'].value_counts() # 2,612 ==> 2021-08-20, # 2,941 ==> 2021-09-28

sc.pp.filter_genes(test3_endo, min_cells=0) # integrated.var에 n_cells 추가 ==> test3에서 이루어졌던 n_cells UPDATE

sc.pp.scale(test3_endo, max_value=10) # ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
# adata.raw.X의 mean 과 std를 output함
sc.tl.pca(test3_endo, n_comps=100, use_highly_variable=True, svd_solver='arpack')

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3_endo, batch_key='batch', n_pcs=20, neighbors_within_batch=5, trim=None) #####
#sce.pp.bbknn(test3_endo, batch_key='batch', n_pcs=50, neighbors_within_batch=5, trim=None)
sc.tl.umap(test3_endo, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
#test3_endo.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
sc.tl.leiden(test3_endo, resolution=0.5, key_added='endo_leiden_r05')
sc.tl.leiden(test3_endo, resolution=1.0, key_added='endo_leiden_r10')

test3_endo.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']

sc.pl.umap(test3_endo, color=['endo_leiden_r05', 'endo_leiden_r10', 'leiden_r05'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')
sc.pl.umap(test3_endo, color=['batch', 'phase', 'percent_mito'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')
sc.pl.umap(test3_endo, color=['batch'], group_by='Month1', add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')

sc.tl.rank_genes_groups(test3_endo, 'endo_leiden_r05', method='wilcoxon', pts=True, key_added='endo_leiden_r05_rank_genes_groups')
#sc.pl.rank_genes_groups(test3_endo, n_genes=5, sharey=False)
sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=10, min_logfoldchange=2, cmap='cividis', show_gene_labels=True, key='endo_leiden_r05_rank_genes_groups')

sc.pl.umap(test3_endo, color='endo_leiden_r05', palette='Set3')

test3_endo.write(filename="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")

markers = ["Pecam1", "Cdh5", "Nos3", "Acta2", "Cnn1", "Tagln", "Rgs5", "Kcnj8", "Col1a1", "Col5a1", "Dpt", "Cd19", "Ighm", "Cd14", "Cd68", "Cd3d"] # Cd3g 없음
sc.pl.stacked_violin(test4, markers, groupby='batch')

result = test3_endo.uns['endo_leiden_r05_rank_genes_groups']
groups = result['names'].dtype.names
#pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals_adj']}).head(5)
#deg_wilcoxon = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pts', 'pts_rest', 'pvals_adj']})
deg_wilcoxon = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals_adj']})
deg_wilcoxon.to_csv("20210930_scanpy_deg_endo.csv", mode='w')

lin = ('0','1','2','3','4','5')
lin = ('1', '3', '5', '4', '2', '0')
#test3_endo.obs['leiden_r05']
test3_endo.obs['endo_leiden_r05'] = test3_endo.obs['endo_leiden_r05'].cat.reorder_categories(list(lin), ordered=True)
sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=10, key='rank_genes_groups', show_gene_labels=True, min_logfoldchange=1, dendrogram=False)

sc.tl.rank_genes_groups(test3_endo, 'endo_leiden_r05', method='wilcoxon', pts=True, key_added='endo_leiden_r05_rank_genes_groups')
sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=10, key='endo_leiden_r05_rank_genes_groups', show_gene_labels=True, min_logfoldchange=1, dendrogram=False, cmap='cividis')
# 아래와 같이 groups을 지정하면 이 group의 순서대로 column이 배치됨 !!!!!!!!!!!!!!!!!!!!!!!! ==> 아닌 것 같다 다시 한 번 확인.
sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=10, key='endo_leiden_r05_rank_genes_groups', groups=['1', '3', '5', '4', '2', '0'], show_gene_labels=True, min_logfoldchange=1, dendrogram=False, cmap='cividis')

#### endo_leiden_r05의 4,5는 아래와 같이 정리 ####
ec_others = {'Pericyte': ['Rgs5', 'Cspg4', 'Kcnj8', 'Des'], 'Lymphatic EC': ['Reln', 'Flt4'], 'Vasa Vasorum': ['Ackr1', 'Lrg1']}
sc.pl.matrixplot(test3_endo, ec_others, layer='magic', groupby='endo_leiden_r05', dendrogram=False, cmap=cmap, standard_scale='var', colorbar_title='Scaled\nexpression', return_fig=False, var_group_rotation=45)



ec_litvinukova = ['Rgcc', 'Car4', 'Sema3g', 'Gja5', 'Plvap']
sc.pl.heatmap(test3_endo, var_names=ec_litvinukova, groupby='leiden_r05', cmap='coolwarm')

artery_schupp = ['Ltbp4', 'Fbln5', 'Bmx', 'Gja4', 'Gja5', 'Efnb2', 'Sox17', 'Sema3g', 'Hey1']
sc.pl.heatmap(test3_endo, var_names=artery_schupp, groupby='leiden_r05', cmap='coolwarm')

capillary_schupp = ['Car4', 'Rgcc', 'Sgk1', 'Sparc', 'Prx']
sc.pl.heatmap(test3_endo, var_names=capillary_schupp, groupby='leiden_r05', cmap='coolwarm')

vein_schupp = ['Nr2f2', 'Selp', 'Vcam1']
sc.pl.heatmap(test3_endo, var_names=vein_schupp, groupby='leiden_r05', cmap='coolwarm')


a = list(test3_endo.obs['batch'].values)
b = list(test3_endo.obs['endo_leiden_r05'].values)
c = list(map(lambda x: a[x] + '_' + b[x], list(range(len(a)))))
test3_endo.obs['aging_leiden_r05'] = c
sc.tl.rank_genes_groups(test3_endo, 'aging_leiden_r05', method='wilcoxon', pts=True, key_added='rank_genes_groups_batch_aging')
# ... storing 'aging_leiden_r05' as categorical ==> 만들어 놓은 aging_leiden_r05를 categorical하게 바꿔야함

lin = ('m01_0', 'm10_0', 'm20_0', 'm01_1', 'm10_1', 'm20_1', 'm01_2', 'm10_2', 'm20_2', 'm01_3', 'm10_3', 'm20_3', 'm01_4', 'm10_4', 'm20_4', 'm01_5', 'm10_5', 'm20_5')
#lin = ('m01_0', 'm01_1', 'm01_2', 'm01_3', 'm01_4','m01_5', 'm10_0', 'm10_1', 'm10_2', 'm10_3', 'm10_4','m10_5', 'm20_0', 'm20_1', 'm20_2', 'm20_3', 'm20_4','m20_5')
test3_endo.obs['aging_leiden_r05']
test3_endo.obs['aging_leiden_r05'] = test3_endo.obs['aging_leiden_r05'].astype('category').cat.reorder_categories(list(lin), ordered=True)

sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=10, groups=['m01_0', 'm10_0', 'm20_0'], groupby='aging_leiden_r05', key='rank_genes_groups_batch_aging', show_gene_labels=True, min_logfoldchange=1, dendrogram=False, cmap='coolwarm')

sc.tl.rank_genes_groups(test3_endo, 'batch', method='wilcoxon', pts=True, key_added='rank_genes_groups_batch')
sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=15, key='rank_genes_groups_batch', show_gene_labels=True, min_logfoldchange=1)


# Diffusion pseudotime
sc.tl.diffmap(test3_endo)
sc.pl.diffmap(test3_endo, color=['batch', 'Pecam1', 'Cdh5'], add_outline=False, legend_loc='right margin', size=70, color_map='CMRmap')

sc.tl.draw_graph(test3, layout='fa', init_pos=None, neighbors_key=None) ## init_pos가 .obsm에 있는 pca, umap, paga 등이 될 수 있다.
sc.pl.draw_graph(test3, color=['batch', 'PECAM1', 'CDH5', 'phase'], add_outline=True, legend_loc='right margin', size=10, color_map='CMRmap')

start_cell = np.isin(test3_endo.obs['endo_leiden_r05'], '0') # boolean numpy array ==> array([False, False, False, ..., False, False, False])
#max_start_id = np.argmin(test3_endo.obsm['X_diffmap'][start_cell,1]) # 262
max_start_id = np.argmax(test3_endo.obsm['X_diffmap'][start_cell,1])
root_id = np.arange(len(start_cell))[start_cell][max_start_id] # 341
test3_endo.uns['iroot'] = root_id

sc.tl.dpt(test3_endo, n_branchings=1, n_dcs=10) # n_branchings를 0으로 하면 (recommended by Scanpy developer) dpt_groups가 생성 안 됨.
#computing Diffusion Pseudotime using n_dcs=10
sc.pl.dpt_groups_pseudotime(test3_endo) # 여기에서 pseudotime trajecgory 확인.

lin = ('2', '0', '3', '1') # DPT pseudotime group ordering에 맞게 배치
test3_endo.obs['dpt_groups'] = test3_endo.obs['dpt_groups'].cat.reorder_categories(list(lin), ordered=True)
sc.pl.dpt_groups_pseudotime(test3_endo) # 다시 ordering에 맞게 plotting
sc.pl.dpt_timeseries(test3_endo[:, test3_endo.var.highly_variable])


################## aEC 에서 비교 ##################
# endo_leiden_r05에서 '0'과 '2'의 비교? '2'는 m01과 m20이 enrich를 판단하기 어렵고, '0'은 확실히 m01이 많음.

# 0 vs 2 (upregulated in 0)
sc.tl.rank_genes_groups(test3_endo, 'endo_leiden_r05', method='wilcoxon', pts=True, key_added='endo_leiden_r05_0vs2_rank_genes_groups', groups=['0'], reference='2')
sc.pl.rank_genes_groups_heatmap(test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '2'])], n_genes=20, groups=['0'], key='endo_leiden_r05_0vs2_rank_genes_groups', show_gene_labels=True, min_logfoldchange=2, dendrogram=False, cmap='cividis')

# 2 vs 0 (upregulated in 2)
sc.tl.rank_genes_groups(test3_endo, 'endo_leiden_r05', method='wilcoxon', pts=True, key_added='endo_leiden_r05_2vs0_rank_genes_groups', groups=['2'], reference='0')
sc.pl.rank_genes_groups_heatmap(test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '2'])], n_genes=20, groups=['2'], key='endo_leiden_r05_2vs0_rank_genes_groups', show_gene_labels=True, min_logfoldchange=2, dendrogram=False, cmap='cividis')












# Only "Arterial" Endothelial cells
test3_aEC = anndata.AnnData(X=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '2'])].layers['counts'], obs=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '2'])].obs, var=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '2'])].var)
test3_aEC.layers["counts"] = test3_aEC.X.copy()

adata_pp = test3_aEC.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = test3_aEC.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))




del adata_pp
test3_aEC.obs['size_factors'] = size_factors

test3_aEC.X /= test3_aEC.obs['size_factors'].values[:, None]
test3_aEC.X = scipy.sparse.csr_matrix(test3_aEC.X) #왜 이게 새로 들어가야될까?????
test3_aEC.X

sc.pp.log1p(test3_aEC) # works on anndata.X
#integrated.X = scipy.sparse.csr_matrix(integrated.X)
test3_aEC.raw = test3_aEC ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

###################################### 2. BBKNN  ######################################
sc.pp.highly_variable_genes(test3_aEC)

test3_aEC.var['highly_variable'].value_counts() # 2,798 ==> 2021-10-14

sc.pp.scale(test3_aEC, max_value=10) # ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
sc.tl.pca(test3_aEC, n_comps=100, use_highly_variable=True, svd_solver='arpack')

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3_aEC, batch_key='batch', n_pcs=7, neighbors_within_batch=3, trim=None)
sc.tl.umap(test3_aEC, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
test3_aEC.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
sc.tl.leiden(test3_aEC, resolution=0.5, key_added='aEC_leiden_r05')
sc.pl.umap(test3_aEC, color=['batch', 'aEC_leiden_r05', 'endo_leiden_r05'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')

#test3_aEC.write(filename="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_aEC.h5ad")
#test3_aEC = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_aEC.h5ad")


# Diffusion pseudotime (########################### Fail ###########################)
sc.tl.diffmap(test3_aEC, n_comps=50) # n_comps=15 (default)
sc.pl.diffmap(test3_aEC, color=['batch', 'Pecam1', 'Cdh5'], add_outline=False, legend_loc='right margin', size=70, color_map='CMRmap')

start_cell = np.isin(test3_aEC.obs['aEC_leiden_r05'], '0') # boolean numpy array ==> array([False, False, False, ..., False, False, False])
max_start_id = np.argmax(test3_aEC.obsm['X_diffmap'][start_cell,2])
root_id = np.arange(len(start_cell))[start_cell][max_start_id] # 253
test3_aEC.uns['iroot'] = root_id

sc.tl.dpt(test3_aEC, n_branchings=1, n_dcs=50) # n_branchings를 0으로 하면 (recommended by Scanpy developer) dpt_groups가 생성 안 됨.
#computing Diffusion Pseudotime using n_dcs=10
sc.pl.dpt_groups_pseudotime(test3_aEC) # 여기에서 pseudotime trajecgory 확인.

lin = ('1', '0', '2') # DPT pseudotime group ordering에 맞게 배치
test3_aEC.obs['dpt_groups'] = test3_aEC.obs['dpt_groups'].cat.reorder_categories(list(lin), ordered=True)
sc.pl.dpt_groups_pseudotime(test3_aEC) # 다시 ordering에 맞게 plotting
sc.pl.dpt_timeseries(test3_aEC[:, test3_aEC.var.highly_variable])

# PAGA

sc.tl.paga(test3_aEC, groups='aEC_leiden_r05')
sc.pl.paga(test3_aEC, color=['aEC_leiden_r05'], threshold=0.2)
#--> added 'pos', the PAGA positions (adata.uns['paga'])
sc.tl.umap(test3_aEC, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='paga', method='umap')
sc.pl.umap(test3_aEC, color=['batch', 'aEC_leiden_r05', 'endo_leiden_r05'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')
sc.pl.paga_compare(test3_aEC, threshold=0.02, size=10, frameon=True, edges=True)











# age difference ####################################################
a = list(test3_endo.obs['batch'].values)
b = list(test3_endo.obs['endo_leiden_r05'].values)
c = list(map(lambda x: a[x] + '_' + b[x], list(range(len(a)))))
test3_endo.obs['aging_endo_leiden_r05'] = c
sc.tl.rank_genes_groups(test3_endo, 'aging_endo_leiden_r05', method='wilcoxon', groups=['m01_0', 'm10_0', 'm20_0'], pts=True, key_added='rank_genes_groups_aging0')
sc.tl.rank_genes_groups(test3_endo, 'aging_endo_leiden_r05', method='wilcoxon', pts=True, key_added='rank_genes_groups_aging2')


sc.pl.rank_genes_groups_heatmap(test3_endo[test3_endo.obs['aging_endo_leiden_r05'].isin(['m01_0', 'm10_0', 'm20_0'])], n_genes=20, groups=['m01_0', 'm10_0', 'm20_0'], key='rank_genes_groups_aging0', show_gene_labels=True, min_logfoldchange=2, dendrogram=False, cmap='cividis', use_raw=False)







# Removing vasa vasorum, lymphatic EC and SMC-like EC
test3_aEC = anndata.AnnData(X=test3_endo[~test3_endo.obs['endo_leiden_r05'].isin(['4', '5'])].layers['counts'], obs=test3_endo[~test3_endo.obs['endo_leiden_r05'].isin(['4', '5'])].obs, var=test3_endo[~test3_endo.obs['endo_leiden_r05'].isin(['4', '5'])].var)
test3_aEC.layers["counts"] = test3_aEC.X.copy()

adata_pp = test3_aEC.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = test3_aEC.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))




del adata_pp
test3_aEC.obs['size_factors'] = size_factors

test3_aEC.X /= test3_aEC.obs['size_factors'].values[:, None]
test3_aEC.X = scipy.sparse.csr_matrix(test3_aEC.X) #왜 이게 새로 들어가야될까?????
test3_aEC.X

sc.pp.log1p(test3_aEC) # works on anndata.X
#integrated.X = scipy.sparse.csr_matrix(integrated.X)
test3_aEC.raw = test3_aEC ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

###################################### 2. BBKNN  ######################################
sc.pp.highly_variable_genes(test3_aEC)

test3_aEC.var['highly_variable'].value_counts() # 2,798 ==> 2021-10-14

sc.pp.scale(test3_aEC, max_value=10) # ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
sc.tl.pca(test3_aEC, n_comps=100, use_highly_variable=True, svd_solver='arpack')

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3_aEC, batch_key='batch', n_pcs=7, neighbors_within_batch=3, trim=None)
sc.tl.umap(test3_aEC, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
test3_aEC.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
sc.tl.leiden(test3_aEC, resolution=0.5, key_added='aEC_leiden_r05')
sc.pl.umap(test3_aEC, color=['batch', 'aEC_leiden_r05', 'endo_leiden_r05'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')

#test3_aEC.write(filename="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_aEC.h5ad")
#test3_aEC = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_aEC.h5ad")







############# cellrank (CytoTrace) ######################

# on test3_endo
import scvelo as scv
from cellrank.tl.kernels import CytoTRACEKernel

# scVelo hack
test3_endo.layers['spliced'] = test3_endo.raw.X
test3_endo.layers['unspliced'] = test3_endo.raw.X

scv.pp.moments(test3_endo)
ctk_endo = CytoTRACEKernel(test3_endo)
sc.pl.umap(test3_endo, color=['ct_pseudotime'], add_outline=False, legend_loc='right margin', size=150)

ctk_endo.compute_transition_matrix(threshold_scheme="soft", nu=0.5)
ctk_endo.compute_projection(basis="umap")

scv.pl.velocity_embedding_stream(test3_endo, color="batch", vkey="T_fwd", basis="umap", legend_loc="right")

from cellrank.tl.estimators import GPCCA
g_fwd = GPCCA(ctk_endo)
g_fwd.compute_schur(n_components=20)
g_fwd.plot_spectrum(real_only=True)









############################################## cellrank ##############################################################

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
plt.rcParams['figure.figsize'] = (7,7)

import cellrank as cr
from cellrank.external.kernels import WOTKernel
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA

test3_endo = sc.read("test3_endo.h5ad")
timepoint = [1] * batches.count('m01') + [10] * batches.count('m10') + [20] * batches.count('m20')
test3_endo.obs['months'] = timepoint

wk = WOTKernel(test3_endo, time_key="months")

wk.compute_initial_growth_rates(organism="mouse", key_added="growth_rate_init")
#WARNING: genes are not in var_names and ignored: ['Hn1', 'Mlf1ip', 'Fam64a']
#WARNING: genes are not in var_names and ignored: ['Adck3', 'Nhlh2', 'Ikbkap', 'Gnb2l1', 'Krt17']

wk.compute_transition_matrix(growth_iters=3, growth_rate_key="growth_rate_init", last_time_point="connectivities")
# ==> ERROR: TypeError: compute_transport_map() got an unexpected keyword argument 'cost_matrix'









############################################## MAGIC ##############################################################

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

test3_endo = anndata.AnnData(X=test3[test3.obs['leiden_r05'].isin(['5', '7'])].layers['counts'], obs=test3[test3.obs['leiden_r05'].isin(['5', '7'])].obs, var=test3[test3.obs['leiden_r05'].isin(['5', '7'])].var)
test3_endo.layers["counts"] = test3_endo.X.copy()

adata_pp = test3_endo.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = test3_endo.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))



del adata_pp
test3_endo.obs['size_factors'] = size_factors
test3_endo.X /= test3_endo.obs['size_factors'].values[:, None]
test3_endo.X = scipy.sparse.csr_matrix(test3_endo.X) #왜 이게 새로 들어가야될까????? # 아니면 ERRROR 남 (highly_variable_genes에서)
test3_endo.X
sc.pp.log1p(test3_endo) # works on anndata.X
#integrated.X = scipy.sparse.csr_matrix(integrated.X)
test3_endo.raw = test3_endo ## ==> log transforamtion 된 것이 raw로 들어가게 됨.
sc.pp.highly_variable_genes(test3_endo)
test3_endo.var['highly_variable'].value_counts() # 2,612 ==> 2021-08-20, # 2,941 ==> 2021-09-28
sc.pp.filter_genes(test3_endo, min_cells=0) # integrated.var에 n_cells 추가 ==> test3에서 이루어졌던 n_cells UPDATE
sc.pp.scale(test3_endo, max_value=10) # ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
sc.tl.pca(test3_endo, n_comps=100, use_highly_variable=True, svd_solver='arpack')
#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3_endo, batch_key='batch', n_pcs=20, neighbors_within_batch=5, trim=None) #####
sc.tl.umap(test3_endo, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
test3_endo.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
sc.tl.leiden(test3_endo, resolution=0.5, key_added='endo_leiden_r05')
test3_endo.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']

import magic
test3_endo_MAGIC = test3_endo.copy()
test3_endo_MAGIC.X = test3_endo_MAGIC.layers['scran_log1p']
test3_endo_MAGIC = magic.MAGIC().fit_transform(test3_endo_MAGIC)
test3_endo.layers['magic'] = test3_endo_MAGIC.X
del test3_endo_MAGIC
sc.pl.umap(test3_endo, layer='magic', color=['Pecam1'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')



############################################## Slingshot ##############################################################
library(slingshot)
library(tidyverse)
library(tidymodels)

sds <- slingshot(Embeddings(seu.subset.endo, "umap"), clusterLabels=seu.subset.endo$seurat_clusters, start.clus=3, stretch=2)
# A Function (that assigns colors to each cell in base R graphics https://bustools.github.io/BUS_notebooks_R/slingshot.html)
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

library(scales)
cell_colors <- cell_pal(seu.subset.endo$orig.ident, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(seu.subset.endo$seurat_clusters, hue_pal())
plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 1)
#lines(sds, lwd = 2, type = 'lineages', col = 'black')
lines(sds, lwd = 2, col = 'black') ###

library(viridis)
nc <- 3
pt <- slingPseudotime(sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end=0.95)

##
par(mfrow = c(nr, nc))
for (i in nms) {
	colors <- pal[cut(pt[,i], breaks=100)]
	plot(reducedDim(sds), col=colors, pch=16, cex=1, main=i)
	#lines(sds, lwd=2, col='black', type='lineages')
	lines(sds, lwd=2, col='black')
}
##





top_hvg <- VariableFeatures(seu.subset.endo)
dat_use <- t(GetAssayData(seu.subset.endo, slot = "data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(sds)[,2], dat_use)
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])
dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)
model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>% set_engine("ranger", importance = "impurity", num.threads = 3) %>% fit(pseudotime ~ ., data = dat_train)
val_results <- dat_val %>% mutate(estimate = predict(model, .[,-1]) %>% pull()) %>% select(truth = pseudotime, estimate)
metrics(data = val_results, truth, estimate)





#### endo_leiden_r05 '4', '5' removal
test3_endoclean = anndata.AnnData(X=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '1', '2', '3'])].layers['counts'], obs=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '1', '2', '3'])].obs, var=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '1', '2', '3'])].var)

test3_endoclean.layers["counts"] = test3_endoclean.X.copy()

adata_pp = test3_endoclean.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']

data_mat = test3_endoclean.X.T

%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))



del adata_pp
test3_endoclean.obs['size_factors'] = size_factors

test3_endoclean.X /= test3_endoclean.obs['size_factors'].values[:, None]
test3_endoclean.X = scipy.sparse.csr_matrix(test3_endoclean.X) #왜 이게 새로 들어가야될까????? # 아니면 ERRROR 남 (highly_variable_genes에서)

test3_endoclean.layers['scran'] = test3_endoclean.X

sc.pp.log1p(test3_endoclean) # works on anndata.X
#integrated.X = scipy.sparse.csr_matrix(integrated.X)
test3_endoclean.layers['scran_log1p'] = test3_endoclean.X

test3_endoclean.raw = test3_endoclean ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

sc.pp.highly_variable_genes(test3_endoclean)

test3_endoclean.var['highly_variable'].value_counts() # 2,891 (2021-12-15)

sc.pp.filter_genes(test3_endoclean, min_cells=0) # integrated.var에 n_cells 추가 ==> test3에서 이루어졌던 n_cells UPDATE

sc.pp.scale(test3_endoclean, max_value=10) # ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
# adata.raw.X의 mean 과 std를 output함
sc.tl.pca(test3_endoclean, n_comps=100, use_highly_variable=True, svd_solver='arpack')

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3_endoclean, batch_key='batch', n_pcs=20, neighbors_within_batch=5, trim=None) #####
sc.tl.umap(test3_endoclean, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
#test3_endoclean.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
sc.tl.leiden(test3_endoclean, resolution=0.5, key_added='endoclean_leiden_r05')
sc.tl.leiden(test3_endoclean, resolution=1.0, key_added='endoclean_leiden_r10')

test3_endoclean.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']

sc.pl.umap(test3_endoclean, color=['endo_leiden_r05', 'endoclean_leiden_r05', 'leiden_r05'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')

# Imputed expression matrix using MAGIC
import magic
MAGIC = test3_endoclean.copy()
MAGIC.X = test3_endoclean.layers['scran_log1p']
MAGIC = magic.MAGIC().fit_transform(MAGIC)

test3_endoclean.layers['magic'] = MAGIC.X





for gene in test3_endo.var_names:
    sc.pl.umap(test3_endo, layer='magic', color=gene, color_map=cmap, legend_loc='on data', add_outline=True, outline_width=(0.01, 0.01), show=False, save='_' + gene + '_test3_endo_magic.png')

for gene in test3.var_names:
    sc.pl.umap(test3, layer='magic', color=gene, color_map=cmap, legend_loc='on data', add_outline=True, outline_width=(0.01, 0.01), show=False, save='_' + gene + '_test3_magic.png')


senescence = ['Il6', 'Il1a', 'Il1b', 'Timp1', 'Mmp3', 'Mmp12', 'Cxcl1', 'Cxcl2', 'Ccl8', 'Cdkn1a', 'Cdkn2a']
sc.tl.score_genes(test3_endo, senescence, score_name='senescence_score')

sb.histplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m01'], kde=True, stat='probability', bins=100, color='#689aff', element='bars')
sb.histplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m10'], kde=True, stat='probability', bins=100, color='#fdbf6f', element='bars')
sb.histplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m20'], kde=True, stat='probability', bins=100, color='#b15928', element='bars')

sb.histplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m01'], kde=True, stat='count', bins=100, color='#689aff', element='bars')
sb.histplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m10'], kde=True, stat='count', bins=100, color='#fdbf6f', element='bars')
sb.histplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m20'], kde=True, stat='count', bins=100, color='#b15928', element='bars')

sb.histplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m01'], kde=True, stat='density', bins=50, color='#689aff', element='step')
sb.histplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m10'], kde=True, stat='density', bins=50, color='#fdbf6f', element='step')
sb.histplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m20'], kde=True, stat='density', bins=50, color='#b15928', element='step')

sb.kdeplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m01'], color='#689aff', fill=True, cumulative=True)
sb.kdeplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m10'], color='#fdbf6f', fill=True, cumulative=True)
sb.kdeplot(test3_endo.obs['senescence_score'][test3_endo.obs['batch']=='m20'], color='#b15928', fill=True, cumulative=True)

test3_endokey = test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '1', '2', '3'])]
test3_endokey = test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['1', '3'])]

sb.kdeplot(test3_endokey.obs['senescence_score'][test3_endo.obs['batch']=='m01'], color='#689aff', fill=True)
sb.kdeplot(test3_endokey.obs['senescence_score'][test3_endo.obs['batch']=='m10'], color='#fdbf6f', fill=True)
sb.kdeplot(test3_endokey.obs['senescence_score'][test3_endo.obs['batch']=='m20'], color='#b15928', fill=True)


dbf = open("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/REACTOME_SASP_genes_mm.txt", 'r')
sasp_reactome = list(map(lambda x: x.strip('\n'), dbf.readlines()))
sc.tl.score_genes(test3_endo, sasp_reactome, score_name='senescence_score')

#### Plotting using scvelo
scv.set_figure_params(style='scvelo', figsize=[5.5,5], frameon=True, color_map=cmap)
scv.pl.scatter(test3_endo, color='batch', groups=[['m01'], ['m10'], ['m20']], ncols=3)
scv.pl.scatter(test3, color='batch', groups=[['m01'], ['m10'], ['m20']], ncols=3)

### ANOVA
import statsmodels.api as sm
from statsmodels.formula.api import ols
fuck = sc.get.obs_df(test3_endo, keys=['Cdkn1a', 'batch'], obsm_keys=(), layer=None, use_raw=True)
fuck_lm = ols('Cdkn1a ~ C(batch)', data=fuck).fit() # C() ==> categorical data (not necessary here because batch is already categorical)
print(sm.stats.anova_lm(fuck_lm, typ=2))

### Chi-squared test
from scipy.stats import chi2_contingency
df = pd.concat([test3.obs['batch'], test3.obs['celltype2']], axis=1)
df_pivot = pd.crosstab(df['batch'], df['celltype2'], normalize=False, margins=True)
chi2_contingency(np.array(list(map(lambda x: [df_pivot['VSMC_1'][x], df_pivot['All'][x] - df_pivot['VSMC_1'][x]], ['m01','m10','m20']))) )

from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

df = pd.concat([test3.obs['batch'], test3.obs['celltype2']], axis=1)
df_pivot = pd.crosstab(df['batch'], df['celltype2'], normalize=False, margins=True)

df_pivot = df_pivot[sorted('VSMC_1  VSMC_2  VSMC_3  FB_1  VSMC_4  EC_1  FB_2  EC_2  VSMC_5  FB_3  All'.split())]
# Immune cell  제거안 했을 때는 밑의 df_pivot.columns를 [:-1] 바꿔줘야함

Chi2s_df, Pvalues_df = list(), list()
months = df_pivot.index[:-1] # ['m01', 'm10', 'm20']
for month in months:
    chi2s, pvalues = list(), list()
    for celltype in df_pivot.columns[1:]:
        chi2, p, dof, ex = chi2_contingency( np.array(list(map(lambda x: [df_pivot[celltype][x], df_pivot['All'][x] - df_pivot[celltype][x]], months))) )
        chi2s.append(chi2)
        pvalues.append(p)
    Chi2s_df.append(chi2s)
    multipletests(pvals=pvalues, alpha=0.01, method='fdr_bh')[1] # B-H correction
    Pvalues_df.append(pvalues)

Chi2s = pd.DataFrame(Chi2s_df, index=df_pivot.index[:-1], columns=df_pivot.columns[1:]).loc['m01'].rename('Chi2')
Pvalues = - ( np.log(pd.DataFrame(Pvalues_df, index=df_pivot.index[:-1], columns=df_pivot.columns[1:]).loc['m01'].rename('-log10Padj')) / np.log(10) )

df_final = pd.concat([Chi2s, Pvalues], axis=1)
df_final

################# EC subclusters ################# Chi-squared test
df = pd.concat([test3_endo.obs['batch'], test3_endo.obs['endo_leiden_r05']], axis=1)
# same as : df = test3_endo.obs[['batch', 'endo_leiden_r05']]
df_pivot = pd.crosstab(df['batch'], df['endo_leiden_r05'], normalize=False, margins=True)

Chi2s_df, Pvalues_df = list(), list()
months = df_pivot.index[:-1] # ['m01', 'm10', 'm20']
for month in months:

    chi2s, pvalues = list(), list()

    for celltype in df_pivot.columns[:-1]:
        chi2, p, dof, ex = chi2_contingency( np.array(list(map(lambda x: [df_pivot[celltype][x], df_pivot['All'][x] - df_pivot[celltype][x]], months))) )
        chi2s.append(chi2)
        pvalues.append(p)

    Chi2s_df.append(chi2s)
    pvalues = multipletests(pvals=pvalues, alpha=0.01, method='fdr_bh')[1] # B-H correction
    Pvalues_df.append(pvalues)

Chi2s = pd.DataFrame(Chi2s_df, index=df_pivot.index[:-1], columns=df_pivot.columns[:-1]).loc['m01'].rename('Chi2')
Pvalues = - ( np.log(pd.DataFrame(Pvalues_df, index=df_pivot.index[:-1], columns=df_pivot.columns[:-1]).loc['m01'].rename('-log10Padj')) / np.log(10) )

df_final = pd.concat([Chi2s, Pvalues], axis=1)
df_final

### Odds ratio calculation (cell type enrichment) m01,m10,m20 간의 비교가 아니라 m01-rest, m10-rest, m20-rest 간의 비교

oddsratio_df, pvalue_df = list(), list()
for month in df_pivot.index[:-1]:
    oddsratio, pvalues = list(), list()
    for celltype in df_pivot.columns[:-1]:
#        table = np.array([ [df_pivot[celltype][month], df_pivot[celltype]['All'] - df_pivot[celltype][month] ], [df_pivot['All'][month], df_pivot['All']['All'] - df_pivot['All'][month]] ])
        table = np.array([ [df_pivot[celltype][month], df_pivot['All'][month] - df_pivot[celltype][month]], [df_pivot[celltype]['All'] - df_pivot[celltype][month], (df_pivot['All']['All'] - df_pivot['All'][month]) - (df_pivot[celltype]['All'] - df_pivot[celltype][month])] ])
        oddsr, p = fisher_exact(table, alternative='two-sided')
        oddsratio.append(oddsr)
        pvalues.append(p)
    oddsratio_df.append(oddsratio)
    pvalues = multipletests(pvals=pvalues, alpha=0.01, method='fdr_bh')[1]
    pvalue_df.append(pvalues)

Odds = pd.DataFrame(oddsratio_df, index=df_pivot.index[:-1], columns=df_pivot.columns[:-1])
Pvalues = pd.DataFrame(pvalue_df, index=df_pivot.index[:-1], columns=df_pivot.columns[:-1])

df_final = pd.concat([Odds, Pvalues], axis=0)
df_final

df_final.iloc[:3].T.plot.bar(color=batch_palette)

Odds.T.iloc[:-4,:].plot.bar(color=batch_palette)
Odds.T.iloc.plot.bar(color=batch_palette) # ALL



#celltype2              VSMC_1              VSMC_2              VSMC_3                FB_1  ...            B-lympho               MΦ             T-lmpho    #        Neuronal
#oddsratio  0.5961995451094636  1.3527062235318188  0.9149824472209976  0.8098978964599354  ...  15.737582161577503  1.8023622383825  2.1760607433292347  #0.5508153756552125

[1 rows x 14 columns]


pd.DataFrame({'celltype2':df_pivot.columns[:-1], 'oddsratio':['0.5961995451094636', '1.3527062235318188', '0.9149824472209976', '0.8098978964599354', '0.9441401406198028', '0.6008021113167676', '1.0882776512945411', '1.2827207378272074', '0.33504584893869377', '2.731315912339897', '15.737582161577503', '1.8023622383825', '2.1760607433292347', '0.5508153756552125']}).set_index('celltype2').T


#### Table generation 2021-09-16

a = list(test3_endo.obs['batch'].values)
b = list(test3_endo.obs['endo_leiden_r05'].values)
c = list(map(lambda x: a[x] + '_' + b[x], list(range(len(a)))))
test3_endo.obs['aging_endo_leiden_r05'] = c
lin = ('m01_0', 'm10_0', 'm20_0', 'm01_1', 'm10_1', 'm20_1', 'm01_2', 'm10_2', 'm20_2', 'm01_3', 'm10_3', 'm20_3', 'm01_4', 'm10_4', 'm20_4', 'm01_5', 'm10_5', 'm20_5')
#lin = ('m01_0', 'm01_1', 'm01_2', 'm01_3', 'm01_4','m01_5', 'm10_0', 'm10_1', 'm10_2', 'm10_3', 'm10_4','m10_5', 'm20_0', 'm20_1', 'm20_2', 'm20_3', 'm20_4','m20_5')
test3_endo.obs['aging_endo_leiden_r05']
test3_endo.obs['aging_endo_leiden_r05'] = test3_endo.obs['aging_endo_leiden_r05'].astype('category').cat.reorder_categories(list(lin), ordered=True)

tipcell_markers = ['Kdr', 'Flt4', 'Nrp1', 'Nrp2', 'Pdgfb', 'Dll4', 'Angpt2', 'Apln', 'Unc5b', 'Robo4', 'Plxnd1', 'Efnb2', 'Cxcr4']
sc.tl.score_genes(test3_endo, tipcell_markers, score_name='tipcell_score', use_raw=True)




















# '4', '5' 제거
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

test3_endo2 = anndata.AnnData(X=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '1', '2', '3'])].layers['counts'], obs=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '1', '2', '3'])].obs, var=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '1', '2', '3'])].var)

test3_endo2.layers["counts"] = test3_endo2.X.copy()

# Doublet information
#test3_endo.obs['Doublet'] = integrated.obs['Doublet'].loc[test3_endo.obs.index]

# Doublet removal
#test3_endo = test3_endo[test3_endo.obs['Doublet'] == 'False']

adata_pp = test3_endo2.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = test3_endo2.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))



del adata_pp
test3_endo2.obs['size_factors'] = size_factors

test3_endo2.X /= test3_endo2.obs['size_factors'].values[:, None]
test3_endo2.X = scipy.sparse.csr_matrix(test3_endo2.X) #왜 이게 새로 들어가야될까????? # 아니면 ERRROR 남 (highly_variable_genes에서)

test3_endo2.layers['scran'] = test3_endo2.X

sc.pp.log1p(test3_endo2) # works on anndata.X
#integrated.X = scipy.sparse.csr_matrix(integrated.X)
test3_endo2.layers['scran_log1p'] = test3_endo2.X

test3_endo2.raw = test3_endo2 ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

sc.pp.highly_variable_genes(test3_endo2)
test3_endo2.var['highly_variable'].value_counts() # 2,897 ==> 2022-01-26

sc.pp.filter_genes(test3_endo2, min_cells=0) # integrated.var에 n_cells 추가 ==> test3에서 이루어졌던 n_cells UPDATE

sc.pp.scale(test3_endo2, max_value=10) # ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
# adata.raw.X의 mean 과 std를 output함
sc.tl.pca(test3_endo2, n_comps=100, use_highly_variable=True, svd_solver='arpack')

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3_endo2, batch_key='batch', n_pcs=15, neighbors_within_batch=5, trim=None) #####
#sce.pp.bbknn(test3_endo, batch_key='batch', n_pcs=50, neighbors_within_batch=5, trim=None)
sc.tl.umap(test3_endo2, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')

sc.tl.leiden(test3_endo2, resolution=0.5, key_added='endo_leiden2_r05')
sc.tl.leiden(test3_endo2, resolution=1.0, key_added='endo_leiden2_r10')

test3_endo2.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']

sc.pl.umap(test3_endo2, color=['endo_leiden2_r05', 'endo_leiden2_r10', 'endo_leiden_r05'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')
sc.pl.umap(test3_endo2, color=['batch', 'phase', 'percent_mito'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')
sc.pl.umap(test3_endo2, color=['batch'], group_by='Month1', add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')

sc.tl.rank_genes_groups(test3_endo2, 'endo_leiden_r05', method='wilcoxon', pts=True, key_added='endo_leiden_r05_rank_genes_groups')
#sc.pl.rank_genes_groups(test3_endo, n_genes=5, sharey=False)
sc.pl.rank_genes_groups_heatmap(test3_endo2, n_genes=10, min_logfoldchange=2, cmap='cividis', show_gene_labels=True, key='endo_leiden_r05_rank_genes_groups')

#test3_endo.write(filename="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")



######################## EC endo_leiden_r05 이름 바꾸기 ########################
endo_leiden_to_celltype_dict = {'0': 'EC_1',
'1': 'EC_4',
'2': 'EC_2',
'3': 'EC_3',
'4': 'EC_5',
'5': 'EC_6'}
test3_endo.obs['EC_subclusters'] = test3_endo.obs['endo_leiden_r05'].map(lambda x: endo_leiden_to_celltype_dict[x]).astype('category')


reordered = ('EC_1', 'EC_2', 'EC_3', 'EC_4')
test3_endo2.obs['EC_subclusters'] = test3_endo2.obs['EC_subclusters'].cat.reorder_categories(list(reordered), ordered=True)
df = test3_endo2.obs[['EC_subclusters', 'phase']]
ax = pd.crosstab(df['EC_subclusters'], df['phase'], normalize='index', margins=True).plot.bar(stacked=True)

################# EC subclusters ################# Chi-squared test (w/o EC_5 and EC_6)
df = test3_endo2.obs[['batch', 'EC_subclusters']]
df_pivot = pd.crosstab(df['batch'], df['EC_subclusters'], normalize=False, margins=True)

Chi2s_df, Pvalues_df = list(), list()
months = df_pivot.index[:-1] # ['m01', 'm10', 'm20']
for month in months:

    chi2s, pvalues = list(), list()

    for celltype in df_pivot.columns[:-1]:
        chi2, p, dof, ex = chi2_contingency( np.array(list(map(lambda x: [df_pivot[celltype][x], df_pivot['All'][x] - df_pivot[celltype][x]], months))) )
        chi2s.append(chi2)
        pvalues.append(p)

    Chi2s_df.append(chi2s)
    pvalues = multipletests(pvals=pvalues, alpha=0.01, method='fdr_bh')[1] # B-H correction
    Pvalues_df.append(pvalues)

Chi2s = pd.DataFrame(Chi2s_df, index=df_pivot.index[:-1], columns=df_pivot.columns[:-1]).loc['m01'].rename('Chi2')
Pvalues = - ( np.log(pd.DataFrame(Pvalues_df, index=df_pivot.index[:-1], columns=df_pivot.columns[:-1]).loc['m01'].rename('-log10Padj')) / np.log(10) )

df_final = pd.concat([Chi2s, Pvalues], axis=1)
df_final

### Odds ratio calculation (cell type enrichment) m01,m10,m20 간의 비교가 아니라 m01-rest, m10-rest, m20-rest 간의 비교

oddsratio_df, pvalue_df = list(), list()
for month in df_pivot.index[:-1]:
    oddsratio, pvalues = list(), list()
    for celltype in df_pivot.columns[:-1]:
#        table = np.array([ [df_pivot[celltype][month], df_pivot[celltype]['All'] - df_pivot[celltype][month] ], [df_pivot['All'][month], df_pivot['All']['All'] - df_pivot['All'][month]] ])
        table = np.array([ [df_pivot[celltype][month], df_pivot['All'][month] - df_pivot[celltype][month]], [df_pivot[celltype]['All'] - df_pivot[celltype][month], (df_pivot['All']['All'] - df_pivot['All'][month]) - (df_pivot[celltype]['All'] - df_pivot[celltype][month])] ])
        oddsr, p = fisher_exact(table, alternative='two-sided')
        oddsratio.append(oddsr)
        pvalues.append(p)
    oddsratio_df.append(oddsratio)
    pvalues = multipletests(pvals=pvalues, alpha=0.01, method='fdr_bh')[1]
    pvalue_df.append(pvalues)

Odds = pd.DataFrame(oddsratio_df, index=df_pivot.index[:-1], columns=df_pivot.columns[:-1])
Pvalues = pd.DataFrame(pvalue_df, index=df_pivot.index[:-1], columns=df_pivot.columns[:-1])

df_final = pd.concat([Odds, Pvalues], axis=0)
df_final

df_final.iloc[:3].T.plot.bar(color=batch_palette)
plt.tight_layout()

Odds.T.iloc[:-4,:].plot.bar(color=batch_palette)
Odds.T.iloc.plot.bar(color=batch_palette) # ALL

############ MAST testing ############
# Create new Anndata object for use in MAST (scran normalized and log transformed data)
test3_endo_test = test3_endo.copy()
test3_endo_test.X = test3_endo.layers['scran_log1p']
test3_endo_test.obs['n_genes'] = (test3_endo_test.X > 0).sum(1)








