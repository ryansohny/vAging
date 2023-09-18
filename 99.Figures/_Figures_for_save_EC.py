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
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#FFFF00", "#000000", "#0066CC"])
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
lin = ('EC1', 'EC2', 'EC3', 'EC4', 'EC5', 'EC6')
test3_endo.obs['Subpopulation of Endothelial Cells'] = test3_endo.obs['Subpopulation of Endothelial Cells'].cat.reorder_categories(list(lin), ordered=True)

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

uni_ec = {'Universal EC markers': ["Pecam1", "Cdh5", 'Erg', 'Vwf', "Nos3", 'Procr', 'Icam1', 'Cd34', 'Cd36']}
ax_dict = sc.pl.matrixplot(test3_endo, uni_ec, layer='magic', groupby='Subpopulation of Endothelial Cells', dendrogram=False, cmap='viridis', standard_scale='var', colorbar_title='Scaled\nexpression', var_group_rotation=0, show=False)
ax_dict['mainplot_ax'].set_xticklabels(labels=list(uni_ec.values())[0], fontstyle='italic', rotation=45)

# Diff version
uni_ec = {'Universal EC markers': ["Pecam1", "Cdh5", "Kdr", 'Erg', 'Vwf', "Nos3", 'Procr', 'Icam1', 'Cd34', 'Cd36', ]}
ec_others = {'Pericyte markers': ['Rgs5', 'Cspg4', 'Kcnj8', 'Des'], 'Lymphatic EC markers': ['Reln', 'Flt4'], 'Vasa Vasorum markers': ['Ackr1', 'Lrg1']}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sns.despine()
ax_dict1 = sc.pl.matrixplot(test3_endo, uni_ec, layer='magic', groupby='Subpopulation of Endothelial Cells', dendrogram=False, cmap='viridis', standard_scale='var', colorbar_title='Scaled\nexpression', var_group_rotation=0, show=False, ax=axes[0])
ax_dict1['mainplot_ax'].set_xticklabels(labels=list(uni_ec.values())[0], fontstyle='italic', rotation=45)
ax_dict1['color_legend_ax'].remove()

ax_dict2 = sc.pl.matrixplot(test3_endo, ec_others, layer='magic', groupby='Subpopulation of Endothelial Cells', dendrogram=False, cmap='viridis', standard_scale='var', colorbar_title='Scaled\nexpression', var_group_rotation=0, show=False, ax=axes[1])
ax_dict2['mainplot_ax'].set_xticklabels(labels=list(x for xs in list(ec_others.values()) for x in xs), fontstyle='italic', rotation=45)
plt.tight_layout()



# 아래 고치기
fig, axes = plt.subplots(1, 2, figsize=(17, 6))
sns.despine()
mp1 = sc.pl.matrixplot(test3_endo, uni_ec, layer='magic', groupby='Subpopulation of Endothelial Cells', dendrogram=False, cmap='viridis', standard_scale='var', colorbar_title='Scaled\nexpression', var_group_rotation=0, return_fig=True, ax=axes[0])
mp1.style(edge_color='none')
mp1.show()
mp1.get_axes()['mainplot_ax'].set_xticklabels(labels=list(uni_ec.values())[0], fontstyle='italic', rotation=45)
mp1.get_axes()['color_legend_ax'].remove()
axes[0].axhline(y=0.71, xmin=0.015, xmax=0.645, linewidth=2.0, c="black")

mp2 = sc.pl.matrixplot(test3_endo, ec_others, layer='magic', groupby='Subpopulation of Endothelial Cells', dendrogram=False, cmap='viridis', standard_scale='var', colorbar_title='Scaled\nexpression', var_group_rotation=0, return_fig=True, ax=axes[1])
mp2.var_group_rotation = 30
mp2.show()
mp2.get_axes()['mainplot_ax'].set_xticklabels(labels=list(x for xs in list(ec_others.values()) for x in xs), fontstyle='italic', rotation=45)
mp2.get_axes()
axes[1].axhline(y=0.71, xmin=0.016, xmax=0.28, linewidth=2.0, c="black")
axes[1].axhline(y=0.71, xmin=0.325, xmax=0.435, linewidth=2.0, c="black")
axes[1].axhline(y=0.71, xmin=0.485, xmax=0.59, linewidth=2.0, c="black")

axes[0].text(-0.1, 0.8, "A", size=20, weight='bold')
axes[1].text(-0.1, 0.8, "B", size=20, weight='bold')

plt.tight_layout()


test3_endo2 = test3_endo[~test3_endo.obs['Subpopulation of Endothelial Cells'].isin(['EC5', 'EC6'])].copy()
colormap = {'EC1': '#8dd3c7',
            'EC2': '#80b1d3',
            'EC3': '#fccde5',
            'EC4': '#bebada'}
sc.pp.filter_genes(test3_endo2, min_cells=0)

#sc.pl.umap(test3_endo2, color='Subpopulation of Endothelial Cells', palette=colormap) # update colormap
#sns.despine()
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#FFFF00", "#000000", "#197fe6"])

sc.tl.rank_genes_groups(test3_endo2, 'Subpopulation of Endothelial Cells', method='wilcoxon', groups=['EC1', 'EC2', 'EC3', 'EC4'], pts=True, key_added='rank_genes_groups_subEC')

#sc.pl.rank_genes_groups_heatmap(test3_endo2, n_genes=20, groups=['EC1', 'EC2', 'EC3', 'EC4'], key='rank_genes_groups_subEC', show_gene_labels=True, min_logfoldchange=1, dendrogram=False, cmap=cmap2, use_raw=False)

deg_ec_result = test3_endo2.uns['rank_genes_groups_subEC']
#deg_ec_result.keys() #dict_keys(['params', 'pts', 'pts_rest', 'names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges'])


deg_ec_groups = deg_ec_result['names'].dtype.names
deg_ec_wilcoxon = pd.DataFrame({group + '_' + key: deg_ec_result[key][group] for group in deg_ec_groups for key in ['names', 'logfoldchanges', 'scores', 'pvals_adj', 'pts', 'pts_rest']})

'''
a = list(test3_endo2.obs['batch'].values)
b = list(test3_endo2.obs['Subpopulation of Endothelial Cells'].values)
c = list(map(lambda x: a[x] + '_' + b[x], list(range(len(a)))))
test3_endo2.obs['aging_subec'] = c
sc.tl.rank_genes_groups(test3_endo2, 'aging_subec', method='wilcoxon', groups=['m01_EC4', 'm10_EC4', 'm20_EC4'], pts=True, key_added='rank_genes_groups_agingEC4')
sc.pl.rank_genes_groups_heatmap(test3_endo2[test3_endo2.obs['aging_subec'].isin(['m01_EC4', 'm10_EC4', 'm20_EC4'])], n_genes=10, groups=['m01_EC4', 'm10_EC4', 'm20_EC4'], key='rank_genes_groups_agingEC4', show_gene_labels=True, min_logfoldchange=1, dendrogram=False, cmap=cmap2, use_raw=False)
sc.pl.rank_genes_groups_heatmap(test3_endo2[test3_endo2.obs['aging_subec'].isin(['m01_EC4', 'm10_EC4', 'm20_EC4'])], n_genes=10, groups=['m01_EC4', 'm10_EC4', 'm20_EC4'], key='rank_genes_groups_agingEC4', show_gene_labels=True, min_logfoldchange=2, dendrogram=False, cmap=cmap2, use_raw=False)
'''

import gseapy as gp
gp.get_library_name(organism='Mouse')
total_no_ec = test3_endo2.shape[0]  # number of cells in EC cluster (test3_endo2) => for identifying expressed genes which later be used as backgroud gene set!!
background_genesets_ec = test3_endo2.var[(test3_endo2.var['n_cells'] > total_no_ec * 0.005)]  # amounts to 0.5% of total EC cells # sc.pp.filter_genes(test3_endo2, min_cells=0) this must be done prior to execution

# Test1
ec1_log2fc_up = list(deg_ec_wilcoxon[(deg_ec_wilcoxon['EC1_pvals_adj'] < 0.00001) & (deg_ec_wilcoxon['EC1_logfoldchanges'] >= 1.5)]['EC1_names'].values)
ec2_log2fc_up = list(deg_ec_wilcoxon[(deg_ec_wilcoxon['EC2_pvals_adj'] < 0.00001) & (deg_ec_wilcoxon['EC2_logfoldchanges'] >= 1.5)]['EC2_names'].values)
ec3_log2fc_up = list(deg_ec_wilcoxon[(deg_ec_wilcoxon['EC3_pvals_adj'] < 0.00001) & (deg_ec_wilcoxon['EC3_logfoldchanges'] >= 1.5)]['EC3_names'].values)
ec4_log2fc_up = list(deg_ec_wilcoxon[(deg_ec_wilcoxon['EC4_pvals_adj'] < 0.00001) & (deg_ec_wilcoxon['EC4_logfoldchanges'] >= 1.5)]['EC4_names'].values)

go_ec1_log2fc_up = gp.enrichr(gene_list=ec1_log2fc_up,
                              gene_sets=['GO_Biological_Process_2021'],
                              background=list(background_genesets_ec.index),
                              organism='Mouse',
                              outdir=None,
                              cutoff=0.05,
                              no_plot=True) 
go_ec1_log2fc_up.res2d.Term = go_ec1_log2fc_up.res2d.Term.str.split(" \(GO").str[0]
del go_ec1_log2fc_up.res2d['Combined Score']
go_ec1_log2fc_up.res2d['-log10(Adusted P-value)'] = -np.log10(go_ec1_log2fc_up.res2d['Adjusted P-value'])
ax = gp.dotplot(go_ec1_log2fc_up.res2d.sort_values(by='Odds Ratio', ascending=False)[:10],  
                figsize=(12,10),
                size=10, 
                top_term=10,
                fontsize=0.1,
                x='-log10(Adusted P-value)',
                show_ring=False,
                cmap='viridis_r',
                column = 'Odds Ratio',
                title="EC1 UP (log2FC>1.5)")
plt.tight_layout()
plt.savefig("GO_EC1_UP.pdf")

go_ec2_log2fc_up = gp.enrichr(gene_list=ec2_log2fc_up,
                              gene_sets=['GO_Biological_Process_2021'],
                              background=list(background_genesets_ec.index),
                              organism='Mouse',
                              outdir=None,
                              cutoff=0.05,
                              no_plot=True) 
go_ec2_log2fc_up.res2d.Term = go_ec2_log2fc_up.res2d.Term.str.split(" \(GO").str[0]
del go_ec2_log2fc_up.res2d['Combined Score']
go_ec2_log2fc_up.res2d['-log10(Adusted P-value)'] = -np.log10(go_ec2_log2fc_up.res2d['Adjusted P-value'])
ax = gp.dotplot(go_ec2_log2fc_up.res2d.sort_values(by='Odds Ratio', ascending=False)[:10],  
                figsize=(12,10),
                size=10, 
                top_term=10,
                fontsize=0.1,
                x='-log10(Adusted P-value)',
                show_ring=False,
                cmap='viridis_r',
                column = 'Odds Ratio',
                title="EC2 UP (log2FC>1.5)")
plt.tight_layout()
plt.savefig("GO_EC2_UP.pdf")

go_ec3_log2fc_up = gp.enrichr(gene_list=ec3_log2fc_up,
                              gene_sets=['GO_Biological_Process_2021'],
                              background=list(background_genesets_ec.index),
                              organism='Mouse',
                              outdir=None,
                              cutoff=0.05,
                              no_plot=True) 
go_ec3_log2fc_up.res2d.Term = go_ec3_log2fc_up.res2d.Term.str.split(" \(GO").str[0]
del go_ec3_log2fc_up.res2d['Combined Score']
go_ec3_log2fc_up.res2d['-log10(Adusted P-value)'] = -np.log10(go_ec3_log2fc_up.res2d['Adjusted P-value'])
ax = gp.dotplot(go_ec3_log2fc_up.res2d.sort_values(by='Odds Ratio', ascending=False)[:10],  
                figsize=(12,10),
                size=10, 
                top_term=10,
                fontsize=0.1,
                x='-log10(Adusted P-value)',
                show_ring=False,
                cmap='viridis_r',
                column = 'Odds Ratio',
                title="EC3 UP (log2FC>1.5)")
plt.tight_layout()
plt.savefig("GO_EC3_UP.pdf")

go_ec4_log2fc_up = gp.enrichr(gene_list=ec4_log2fc_up,
                              gene_sets=['GO_Biological_Process_2021'],
                              background=list(background_genesets_ec.index),
                              organism='Mouse',
                              outdir=None,
                              cutoff=0.05,
                              no_plot=True) 
go_ec4_log2fc_up.res2d.Term = go_ec4_log2fc_up.res2d.Term.str.split(" \(GO").str[0]
del go_ec4_log2fc_up.res2d['Combined Score']
go_ec4_log2fc_up.res2d['-log10(Adusted P-value)'] = -np.log10(go_ec4_log2fc_up.res2d['Adjusted P-value'])
ax = gp.dotplot(go_ec4_log2fc_up.res2d.sort_values(by='Odds Ratio', ascending=False)[:10],  
                figsize=(12,10),
                size=10, 
                top_term=10,
                fontsize=0.1,
                x='-log10(Adusted P-value)',
                show_ring=False,
                cmap='viridis_r',
                column = 'Odds Ratio',
                title="EC4 UP (log2FC>1.5)")
plt.tight_layout()
plt.savefig("GO_EC4_UP.pdf")





# New
go_ec1_log2fc_up = gp.enrichr(gene_list=ec1_log2fc_up,
                              gene_sets=['GO_Biological_Process_2021'],
                              background=list(background_genesets_ec.index),
                              organism='Mouse',
                              outdir=None,
                              cutoff=0.05,
                              no_plot=True) 
go_ec1_log2fc_up.res2d.Term = go_ec1_log2fc_up.res2d.Term.str.split(" \(GO").str[0]
go_ec1_log2fc_up.res2d['-log10(Adusted P-value)'] = -np.log10(go_ec1_log2fc_up.res2d['Adjusted P-value'])
ax = gp.dotplot(go_ec1_log2fc_up.res2d.sort_values(by='Combined Score', ascending=False)[:10],  
                figsize=(16,10),
                size=10, 
                top_term=10,
                fontsize=0.1,
                x='-log10(Adusted P-value)',
                show_ring=False,
                cmap='viridis_r',
                column = 'Combined Score',
                title="EC1 UP (log2FC>1.5)")
plt.tight_layout()


go_ec2_log2fc_up = gp.enrichr(gene_list=ec2_log2fc_up,
                              gene_sets=['GO_Biological_Process_2021'],
                              background=list(background_genesets_ec.index),
                              organism='Mouse',
                              outdir=None,
                              cutoff=0.05,
                              no_plot=True) 
go_ec2_log2fc_up.res2d.Term = go_ec2_log2fc_up.res2d.Term.str.split(" \(GO").str[0]
go_ec2_log2fc_up.res2d['-log10(Adusted P-value)'] = -np.log10(go_ec2_log2fc_up.res2d['Adjusted P-value'])
ax = gp.dotplot(go_ec2_log2fc_up.res2d.sort_values(by='Combined Score', ascending=False)[:10],  
                figsize=(16,10),
                size=10, 
                top_term=10,
                fontsize=0.1,
                x='-log10(Adusted P-value)',
                show_ring=False,
                cmap='viridis_r',
                column = 'Combined Score',
                title="EC2 UP (log2FC>1.5)")
plt.tight_layout()

go_ec3_log2fc_up = gp.enrichr(gene_list=ec3_log2fc_up,
                              gene_sets=['GO_Biological_Process_2021'],
                              background=list(background_genesets_ec.index),
                              organism='Mouse',
                              outdir=None,
                              cutoff=0.05,
                              no_plot=True) 
go_ec3_log2fc_up.res2d.Term = go_ec3_log2fc_up.res2d.Term.str.split(" \(GO").str[0]
go_ec3_log2fc_up.res2d['-log10(Adusted P-value)'] = -np.log10(go_ec3_log2fc_up.res2d['Adjusted P-value'])
ax = gp.dotplot(go_ec3_log2fc_up.res2d.sort_values(by='Combined Score', ascending=False)[:10],  
                figsize=(16,10),
                size=10, 
                top_term=10,
                fontsize=0.1,
                x='-log10(Adusted P-value)',
                show_ring=False,
                cmap='viridis_r',
                column = 'Combined Score',
                title="EC3 UP (log2FC>1.5)")
plt.tight_layout()


go_ec4_log2fc_up = gp.enrichr(gene_list=ec4_log2fc_up,
                              gene_sets=['GO_Biological_Process_2021'],
                              background=list(background_genesets_ec.index),
                              organism='Mouse',
                              outdir=None,
                              cutoff=0.05,
                              no_plot=True) 
go_ec4_log2fc_up.res2d.Term = go_ec4_log2fc_up.res2d.Term.str.split(" \(GO").str[0]
go_ec4_log2fc_up.res2d['-log10(Adusted P-value)'] = -np.log10(go_ec4_log2fc_up.res2d['Adjusted P-value'])
ax = gp.dotplot(go_ec4_log2fc_up.res2d.sort_values(by='Combined Score', ascending=False)[:10],  
                figsize=(16,10),
                size=10, 
                top_term=10,
                fontsize=0.1,
                x='-log10(Adusted P-value)',
                show_ring=False,
                cmap='viridis_r',
                column = 'Combined Score',
                title="EC4 UP (log2FC>1.5)")
plt.tight_layout()


















# concat results
go_vsmc_cluster0_log2fc_1up.res2d['UP_DW'] = "UP"
go_vsmc_cluster0_log2fc_1down.res2d['UP_DW'] = "DOWN"
go_vsmc_cluster0_res = pd.concat([go_vsmc_cluster0_log2fc_1up.res2d.head(), go_vsmc_cluster0_log2fc_1down.res2d.head()]) # 각각 5개씩만


ax = gp.barplot(go_vsmc_cluster0_res,
                group=["UP_DW"],
                title='vSMC Cluster_0 UP (log2FC>2) DOWN (log2FC<-2)',
                column='Adjusted P-value',
                color= ['b', 'r'],
                figsize=(15,10))
plt.tight_layout()

ax = gp.dotplot(go_vsmc_cluster0_log2fc_2up.res2d, 
                figsize=(10, 10), 
                size=100, 
                cmap='NbDr_r',
                show_ring=True,
                title="vSMC Cluster_0 UP (log2FC>2) DOWN (log2FC<-2)")

#ofname='./GO_EC_subclsters/Barplot_EC_1_GO_Biological_Process.pdf')





# Senescence-associated Secreted Protein (SASP)
dbf = open("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/REACTOME_SASP_genes_mm.txt", 'r')
sasp_reactome = list(map(lambda x: x.strip('\n'), dbf.readlines()))
sc.tl.score_genes(test3_endo2, sasp_reactome, score_name='senescence_score')

df_senescence_score = test3_endo2.obs[['senescence_score', 'Age']]

ax = sns.violinplot(df_senescence_score, x="Age", y="senescence_score", inner=None, palette=batch_palette)
sns.stripplot(df_senescence_score, x="Age", y="senescence_score", color='black', size=1, alpha=0.5, jitter=True, ax=ax)
ax.set_xlabel('')
ax.set_ylabel('Senescence Score')
