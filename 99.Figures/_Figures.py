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
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'Arial'
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
batch_palette=['#689aff', '#fdbf6f', '#b15928']
%matplotlib
%autoindent

test3 = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3.h5ad")
test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")

cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#000000", "#8b0a50"])
###### Figure 1 ######

# Figure 1B
test3.obs['Age'] = test3.obs['batch']
sc.pl.umap(test3, color=['Age'], add_outline=False, legend_loc='right margin', color_map=cmap, size=20)
test3.uns['Age_colors'] = ['#689aff', '#fdbf6f', '#b15928']
sc.pl.umap(test3, color=['Age'], add_outline=False, legend_loc='right margin', color_map=cmap, size=20)
sns.despine()

# Figure 1C
leiden_to_celltype_dict = {'0': 'Vascular smooth muscle cells',
'1': 'Vascular smooth muscle cells',
'2': 'Vascular smooth muscle cells',
'3': 'Fibroblasts',
'4': 'Vascular smooth muscle cells',
'5': 'Endothelial cells',
'6': 'Fibroblasts',
'7': 'Endothelial cells',
'8': 'Vascular smooth muscle cells',
'9': 'Fibroblasts',
'10': 'B cells',
'11': 'M\u03A6',
'12': 'T cells',
'13': 'Neuronal cells'}
test3.obs['celltype'] = test3.obs['leiden_r05'].map(lambda x: leiden_to_celltype_dict[x]).astype('category')

lin = ('Endothelial cells', 'Vascular smooth muscle cells', 'Fibroblasts', 'B cells', 'M\u03A6', 'T cells', 'Neuronal cells')
test3.obs['Annotated Cell Types'] = test3.obs['celltype'].cat.reorder_categories(list(lin), ordered=True)

sc.pl.umap(test3, color=['Annotated Cell Types'], add_outline=False, legend_loc='right margin', color_map=cmap, size=20, palette='tab20b')
sns.despine()

# Figure 1D
celltype_marker = {'EC markers': ['Pecam1', 'Cdh5', 'Vwf', 'Nos3'],
'VSMC markers': ['Acta2', 'Tagln', 'Cnn1', 'Cnn2'],
'FB markers': ['Dpt', 'Col1a1', 'Col5a1', 'Pdgfra'],
'Bc markers': ['Ighm', 'Cd19'],
'MΦ markers':['Cd14', 'Cd68'],
'Tc markers':['Cd3d', 'Cd3g'],
'Neuronal markers':['Mbp', 'Cnp']}
#reordered = ('Endothelial cells', 'Smooth muscle cells', 'Fibroblasts', 'B cells', 'M\u03A6', 'T cells', 'Neuronal')
#reordered = ('Endothelial cells', 'Smooth muscle cells', 'Fibroblasts', 'B cells', 'Macrophages', 'T cells')
#test3.obs['Annotated Cell Types'] = test3.obs['Annotated Cell Types'].cat.reorder_categories(list(reordered), ordered=True)
combined_celltype_marker = list(x for xs in list(celltype_marker.values()) for x in xs)


ax_dict = sc.pl.matrixplot(test3, celltype_marker, layer='magic', groupby='Annotated Cell Types', dendrogram=False, cmap=cmap, standard_scale='var', colorbar_title='Scaled\nexpression', var_group_rotation=45, show=False)
ax_dict['mainplot_ax'].set_xticklabels(labels=combined_celltype_marker, fontstyle='italic', rotation=45)
plt.xlabel("")
plt.tight_layout()

dp = sc.pl.dotplot(test3, celltype_marker, layer='magic', groupby='Annotated Cell Types', return_fig=True)
dp.add_totals(size=1.5, color=['#393b79', '#9c9ede', '#b5cf6b', '#e7ba52', '#ad494a', '#7b4173', '#de9ed6']).legend(colorbar_title='log(SizeFactorNormlized+1)', width=1.5, show_size_legend=False, show_colorbar=False).style(cmap='Reds', dot_edge_color='black', dot_edge_lw=1, size_exponent=1.5, grid=True, x_padding=0.4, y_padding=0.6).swap_axes().show()

# Figure 1E
#colormap = dict(zip(list(test3.obs['Annotated Cell Types'].unique()), list(test3.uns['Annotated Cell Types_colors'])))
colormap = {'Endothelial cells': '#a6cee3',
            'Smooth muscle cells': '#b2df8a',
            'Fibroblasts': '#fb9a99',
            'B cells': '#ff7f00',
            'MΦ': '#6a3d9a',
            'T cells': '#b15928'}

df = pd.concat([test3.obs['Age'], test3.obs['Annotated Cell Types']], axis=1)
ax = pd.crosstab(df['Age'], df['Annotated Cell Types'], normalize=0).sort_values(by='Age', ascending=False).plot.barh(stacked=True, color=colormap)
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.xlabel('Fraction of Annotated Cell Types')
plt.tight_layout()
sns.despine()

###### Figure 2 ######

# Fig2.A

endo_leiden_to_celltype_dict = {'0': 'EC1',
'1': 'EC4',
'2': 'EC2',
'3': 'EC3',
'4': 'EC5',
'5': 'EC6'}
test3_endo.obs['Subpopulation of Endothelial Cells'] = test3_endo.obs['endo_leiden_r05'].map(lambda x: endo_leiden_to_celltype_dict[x]).astype('category')
sc.pl.umap(test3_endo, color=['Subpopulation of Endothelial Cells'], add_outline=False, legend_loc='right margin', color_map=cmap, palette='Set3')

test3_endo.obs['Age'] = test3_endo.obs['batch']

sc.pl.umap(test3_endo, color=['Age'], add_outline=False, legend_loc='right margin', color_map=cmap)
test3_endo.uns['Age_colors'] = ['#689aff', '#fdbf6f', '#b15928']
sc.pl.umap(test3_endo, color=['Age'], add_outline=False, legend_loc='right margin', color_map=cmap)
sns.despine()

# Fig2.B,C
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

# Fig2.D
colormap = {'EC_1': '#8dd3c7',
            'EC_2': '#80b1d3',
            'EC_3': '#fccde5',
            'EC_4': '#bebada'}

df = pd.concat([ test3_endo[~test3_endo.obs['Subpopulation of Endothelial Cells'].isin(['EC_5', 'EC_6'])].obs['Age'], \
                test3_endo[~test3_endo.obs['Subpopulation of Endothelial Cells'].isin(['EC_5', 'EC_6'])].obs['Subpopulation of Endothelial Cells'] ], axis=1)
ax = pd.crosstab(df['Age'], df['Subpopulation of Endothelial Cells'], normalize=0).sort_values(by='Age', ascending=False).plot.barh(stacked=True, color=colormap)
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.xlabel('Fraction of EC Subpopulations')
plt.tight_layout()
sns.despine()



# Fig2.E
df = test3_endo[~test3_endo.obs['EC_subclusters'].isin(['EC_5', 'EC_6'])].obs[['batch', 'EC_subclusters']]
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

batch_palette= {'m01': '#689aff',
                'm10': '#fdbf6f',
                'm20': '#b15928'}

ax = df_final.iloc[:3].sort_values(axis=1, by='EC_subclusters').T.plot.bar(color=batch_palette, rot=45) # Odds ratio
ax.legend(loc='upper left', bbox_to_anchor=(0.75, 1.0))
ax.set_ylabel('Odds Ratio')
plt.tight_layout()
sns.despine()


# Fig2.F

test3_endo2 = test3_endo[~test3_endo.obs['Subpopulation of Endothelial Cells'].isin(['EC5', 'EC6'])].copy()
colormap = {'EC1': '#8dd3c7',
            'EC2': '#80b1d3',
            'EC3': '#fccde5',
            'EC4': '#bebada'}
sc.pl.umap(test3_endo2, color='Subpopulation of Endothelial Cells', palette=colormap) # update colormap
sns.despine()

lin = ('EC1', 'EC2', 'EC3', 'EC4')
test3_endo2.obs['Subpopulation of Endothelial Cells'] = test3_endo2.obs['Subpopulation of Endothelial Cells'].cat.reorder_categories(list(lin), ordered=True)

sc.tl.rank_genes_groups(test3_endo2, 'Subpopulation of Endothelial Cells', method='wilcoxon', pts=True, key_added='Subpopulation of Endothelial Cells_rank_genes_groups')

ax_dict = sc.pl.rank_genes_groups_heatmap(test3_endo2, n_genes=15, groupby='Subpopulation of Endothelial Cells', key='Subpopulation of Endothelial Cells_rank_genes_groups', groups=['EC1', 'EC2', 'EC3', 'EC4'], show_gene_labels=True, min_logfoldchange=1, dendrogram=False, cmap=cmap, use_raw=False, vmin=-1.5, vmax=1.5, swap_axes=True, show=False, var_group_rotation=90)
ax_dict['heatmap_ax'].set_yticklabels(labels=ax_dict['heatmap_ax'].get_yticklabels(), fontstyle='italic')

result = test3_endo2.uns['Subpopulation of Endothelial Cells_rank_genes_groups']
groups = result['names'].dtype.names
deg_wilcoxon = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals_adj']})


ec_1 = list(deg_wilcoxon[(deg_wilcoxon['EC_1_logfoldchanges'] > 1) & (deg_wilcoxon['EC_1_pvals_adj'] < 0.05)].iloc[0:15]['EC_1_names'].values)
ec_2 = list(deg_wilcoxon[(deg_wilcoxon['EC_2_logfoldchanges'] > 1) & (deg_wilcoxon['EC_2_pvals_adj'] < 0.05)].iloc[0:15]['EC_2_names'].values)
ec_3 = list(deg_wilcoxon[(deg_wilcoxon['EC_3_logfoldchanges'] > 1) & (deg_wilcoxon['EC_3_pvals_adj'] < 0.05)].iloc[0:15]['EC_3_names'].values)
ec_4 = list(deg_wilcoxon[(deg_wilcoxon['EC_4_logfoldchanges'] > 1) & (deg_wilcoxon['EC_4_pvals_adj'] < 0.05)].iloc[0:15]['EC_4_names'].values)

import gseapy as gp
gp.get_library_name(organism='Mouse')

go_EC_1 = gp.enrichr(gene_list=list(ec_1),
                     gene_sets=['GO_Biological_Process_2021'],
                     organism='Mouse',
                     description='GO Biological Process',
                     outdir='./GO_EC_subclsters',
                     cutoff=0.05,
                     no_plot=True)
gp.barplot(go_EC_1.res2d,
           title='GO Biological Process',
           column='Adjusted P-value',
           cutoff=0.05,
           top_term=5,
           cmap='viridis_r',
           ofname='./GO_EC_subclsters/Barplot_EC_1_GO_Biological_Process.pdf')
go_EC_2 = gp.enrichr(gene_list=list(ec_2),
			 gene_sets=['GO_Biological_Process_2021'],
			 organism='Mouse',
			 description='GO Biological Process',
			 outdir='./GO_EC_subclsters',
			 cutoff=0.05,
			 no_plot=True)
gp.barplot(go_EC_2.res2d,
           title='GO Biological Process',
           column='Adjusted P-value',
           cutoff=0.05,
           top_term=5,
           cmap='viridis_r',
           ofname='./GO_EC_subclsters/Barplot_EC_2_GO_Biological_Process.pdf')
go_EC_3 = gp.enrichr(gene_list=list(ec_3),
			 gene_sets=['GO_Biological_Process_2021'],
			 organism='Mouse',
			 description='GO Biological Process',
			 outdir='./GO_EC_subclsters',
			 cutoff=0.05,
			 no_plot=True)
gp.barplot(go_EC_3.res2d,
           title='GO Biological Process',
           column='Adjusted P-value',
           cutoff=0.05,
           top_term=5,
           cmap='viridis_r',
           ofname='./GO_EC_subclsters/Barplot_EC_3_GO_Biological_Process.pdf')
go_EC_4 = gp.enrichr(gene_list=list(ec_4),
			 gene_sets=['GO_Biological_Process_2021'],
			 organism='Mouse',
			 description='GO Biological Process',
			 outdir='./GO_EC_subclsters',
			 cutoff=0.05,
			 no_plot=True)
gp.barplot(go_EC_4.res2d,
           title='GO Biological Process',
           column='Adjusted P-value',
           cutoff=0.05,
           top_term=5,
           cmap='viridis_r',
           ofname='./GO_EC_subclsters/Barplot_EC_4_GO_Biological_Process.pdf')

# Figure 3.A
# Isl1, Twist1, Procr
genes = ['Isl1', 'Twist1', 'Procr']
ax = sc.pl.umap(test3_endo2, color=genes, color_map=cmap, layer='magic', show=False, ncols=len(genes))
for x in range( len(ax) ):
    ax[x].set_title(genes[x], style='italic')
sns.despine()

# Cd34, Prom1, Nkx2-5
genes = ['Cd34', 'Prom1', 'Nkx2-5']
ax = sc.pl.umap(test3_endo2, color=genes, color_map=cmap, layer='magic', show=False, ncols=len(genes))
for x in range( len(ax) ):
    ax[x].set_title(genes[x], style='italic')

sns.despine()

# Isl1 regulon
dfh = open("/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/Isl1_regulon.txt", 'r')
regulons = list()
for i in dfh:
    regulons.append(i.strip())
isl1 = sc.get.obs_df(test3_endo2, keys=['Subpopulation of Endothelial Cells', *regulons])
g = sns.clustermap(isl1.groupby('Subpopulation of Endothelial Cells').mean().T, 
cmap='viridis',
z_score=0, 
standard_scale=None, 
method='ward', 
metric='euclidean')
g.ax_heatmap.set_yticklabels(labels=g.ax_heatmap.get_yticklabels(), fontstyle='italic')
g.cax.set_visible(False)

# Twist1 regulon
dfh = open("/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/Twist1_regulon.txt", 'r')
regulons = list()
for i in dfh:
    regulons.append(i.strip())
twist1 = sc.get.obs_df(test3_endo2, keys=['Subpopulation of Endothelial Cells', *regulons])
g = sns.clustermap(twist1.groupby('Subpopulation of Endothelial Cells').mean().T, 
cmap='viridis',
z_score=0, 
standard_scale=None, 
method='ward', 
metric='euclidean')
g.ax_heatmap.set_yticklabels(labels=g.ax_heatmap.get_yticklabels(), fontstyle='italic')
g.cax.set_visible(False)

# Pparg regulon
dfh = open("/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/Pparg_regulon.txt", 'r')
regulons = list()
for i in dfh:
    regulons.append(i.strip())
pparg = sc.get.obs_df(test3_endo2, keys=['Subpopulation of Endothelial Cells', *regulons])
g = sns.clustermap(pparg.groupby('Subpopulation of Endothelial Cells').mean().T, 
cmap='viridis',
z_score=0, 
standard_scale=None, 
method='ward', 
metric='euclidean')
g.ax_heatmap.set_yticklabels(labels=g.ax_heatmap.get_yticklabels(), fontstyle='italic')
g.cax.set_visible(False)

# Figure 3.E
df = test3_endo2.obs[['EC_subclusters', 'phase']]
ax = pd.crosstab(df['EC_subclusters'], df['phase'], normalize='index', margins=True).plot.bar(stacked=True, rot=45)
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0))
ax.set_ylabel('Proportion of Cell Cycle Phase')
plt.tight_layout()
sns.despine()



# Research Letter Figures
colormap = {'EC1': '#8dd3c7',
            'EC2': '#80b1d3',
            'EC3': '#fccde5',
            'EC4': '#bebada'}
fig, axes = plt.subplots(1,3)
sns.despine()
sc.pl.umap(test3_endo2[test3_endo2.obs['batch'] == 'm01'], color='Subpopulation of Endothelial Cells', palette=colormap, show=False, legend_loc=None, size=200, title='1 month', ax=axes[0])
sc.pl.umap(test3_endo2[test3_endo2.obs['batch'] == 'm10'], color='Subpopulation of Endothelial Cells', palette=colormap, show=False, legend_loc=None, size=200, title='10 months', ax=axes[1])
sc.pl.umap(test3_endo2[test3_endo2.obs['batch'] == 'm20'], color='Subpopulation of Endothelial Cells', palette=colormap, show=False, legend_loc=None, size=200, title='20 months', ax=axes[2])



# Figure X
# 
genes = ['Atf3', 'Atf4', 'Ddit3', 'Ppp1r15a']
ax = sc.pl.umap(test3_endo2, color=genes, color_map=cmap, layer='magic', show=False, ncols=5)
for x in range( len(ax) ):
    ax[x].set_title(genes[x], style='italic')
sns.despine()

genes = ['Nfkb1', 'Nfkb2', 'Rela', 'Relb', 'Rel']
ax = sc.pl.umap(test3_endo2, color=genes, color_map=cmap, layer='magic', show=False, ncols=5)
for x in range( len(ax) ):
    ax[x].set_title(genes[x], style='italic')

sns.despine()

# Figure XX
plot = sc.pl.StackedViolin(test3_endo2, ec4_markers, groupby='Subpopulation of Endothelial Cells', layer='magic', cmap='viridis', figsize=(7,5), return_fig=True, show=False)
plot.swap_axes(swap_axes=True).show()
sns.despine()
plot.get_axes()['mainplot_ax'].set_xticklabels(labels=['EC1', 'EC2', 'EC3', 'EC4'], rotation=45)
plot.get_axes()['mainplot_ax'].set_yticklabels(labels=ec4_markers, fontstyle='italic')

plot = sc.pl.StackedViolin(test3_endo2, ec3_markers, groupby='Subpopulation of Endothelial Cells', layer='magic', cmap='viridis', figsize=(7,5), return_fig=True, show=False)
plot.swap_axes(swap_axes=True).show()
sns.despine()
plot.get_axes()['mainplot_ax'].set_xticklabels(labels=['EC1', 'EC2', 'EC3', 'EC4'], rotation=45)
plot.get_axes()['mainplot_ax'].set_yticklabels(labels=ec3_markers, fontstyle='italic')

ec1_markers = ['Twist1', 'Timp2', 'Mmp14', 'Mmp2', 'Timp4']
ec1_markers = ['Timp2', 'Mmp14', 'Mmp2', 'Timp4']
plot = sc.pl.StackedViolin(test3_endo2, ec1_markers, groupby='Subpopulation of Endothelial Cells', layer='magic', cmap='viridis', figsize=(7,5), return_fig=True, show=False)
plot.swap_axes(swap_axes=True).show()
sns.despine()
plot.get_axes()['mainplot_ax'].set_xticklabels(labels=['EC1', 'EC2', 'EC3', 'EC4'], rotation=45)
plot.get_axes()['mainplot_ax'].set_yticklabels(labels=ec1_markers, fontstyle='italic')
