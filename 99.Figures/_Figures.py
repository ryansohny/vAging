cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#000000", "#8b0a50"])
###### Figure 1 ######

# Figure 1B
test3.obs['Age'] = test3.obs['batch']
sc.pl.umap(test3, color=['Age'], add_outline=False, legend_loc='right margin', color_map=cmap, size=20)
test3.uns['Age_colors'] = ['#689aff', '#fdbf6f', '#b15928']
sc.pl.umap(test3, color=['Age'], add_outline=False, legend_loc='right margin', color_map=cmap, size=20)
sns.despine()

# Figure 1C
test3.obs['Annotated Cell Types'] = test3.obs['celltype']
sc.pl.umap(test3, color=['Annotated Cell Types'], add_outline=False, legend_loc='right margin', color_map=cmap, size=20, palette='Paired')
sns.despine()

# Figure 1D
celltype_marker = {'EC markers': ['Pecam1', 'Cdh5', 'Vwf', 'Nos3'],
'SMC markers': ['Acta2', 'Tagln', 'Cnn1', 'Cnn2'],
'FB markers': ['Dpt', 'Col1a1', 'Col5a1', 'Pdgfra'],
'Bc markers': ['Ighm', 'Cd19'],
'MΦ markers':['Cd14', 'Cd68'],
'Tc markers':['Cd3d', 'Cd3g']}
reordered = ('Endothelial cells', 'Smooth muscle cells', 'Fibroblasts', 'B cells', 'M\u03A6', 'T cells')
#reordered = ('Endothelial cells', 'Smooth muscle cells', 'Fibroblasts', 'B cells', 'Macrophages', 'T cells')
test3.obs['Annotated Cell Types'] = test3.obs['Annotated Cell Types'].cat.reorder_categories(list(reordered), ordered=True)
sc.pl.matrixplot(test3, celltype_marker, layer='magic', groupby='Annotated Cell Types', dendrogram=False, cmap=cmap, standard_scale='var', colorbar_title='Scaled\nexpression', var_group_rotation=45)

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
endo_leiden_to_celltype_dict = {'0': 'EC_1',
'1': 'EC_4',
'2': 'EC_2',
'3': 'EC_3',
'4': 'EC_5',
'5': 'EC_6'}
test3_endo.obs['Subpopulation of Endothelial Cells'] = test3_endo.obs['endo_leiden_r05'].map(lambda x: endo_leiden_to_celltype_dict[x]).astype('category')
sc.pl.umap(test3_endo, color=['Subpopulation of Endothelial Cells'], add_outline=False, legend_loc='right margin', color_map=cmap, palette='Accent')

test3_endo.obs['Age'] = test3_endo.obs['batch']

sc.pl.umap(test3_endo, color=['Age'], add_outline=False, legend_loc='right margin', color_map=cmap)
test3_endo.uns['Age_colors'] = ['#689aff', '#fdbf6f', '#b15928']
sc.pl.umap(test3_endo, color=['Age'], add_outline=False, legend_loc='right margin', color_map=cmap)
sns.despine()

# Fig2.B,C
lin = ('EC_1', 'EC_2', 'EC_3', 'EC_4', 'EC_5', 'EC_6')
#test3_endo.obs['leiden_r05']
test3_endo.obs['Subpopulation of Endothelial Cells'] = test3_endo.obs['Subpopulation of Endothelial Cells'].cat.reorder_categories(list(lin), ordered=True)

uni_ec = {'Universal EC markers': ["Pecam1", "Cdh5", 'Erg', 'Vwf', "Nos3", 'Procr', 'Icam1', 'Cd34', 'Cd36', ]}
sc.pl.matrixplot(test3_endo, uni_ec, groupby='Subpopulation of Endothelial Cells', dendrogram=False, cmap='viridis', standard_scale='var', colorbar_title='Scaled\nexpression', return_fig=False, var_group_rotation=45)

ec_others = {'Pericyte': ['Rgs5', 'Cspg4', 'Kcnj8', 'Des'], 'Lymphatic EC': ['Reln', 'Flt4'], 'Vasa Vasorum': ['Ackr1', 'Lrg1']}
sc.pl.matrixplot(test3_endo, ec_others, groupby='Subpopulation of Endothelial Cells', dendrogram=False, cmap='viridis', standard_scale='var', colorbar_title='Scaled\nexpression', return_fig=False, var_group_rotation=45)

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

lin = ('0','1','2','3','4','5')
lin = ('EC_1', 'EC_2', 'EC_3', 'EC_4')
#test3_endo.obs['leiden_r05']


test3_endo2 = test3_endo[~test3_endo.obs['Subpopulation of Endothelial Cells'].isin(['EC_5', 'EC_6'])].copy()
colormap = {'EC_1': '#8dd3c7',
            'EC_2': '#80b1d3',
            'EC_3': '#fccde5',
            'EC_4': '#bebada'}
sc.pl.umap(test3_endo2, color='Subpopulation of Endothelial Cells', palette=colormap) # update colormap
lin = ('EC_1', 'EC_2', 'EC_3', 'EC_4')
test3_endo2.obs['Subpopulation of Endothelial Cells'] = test3_endo2.obs['Subpopulation of Endothelial Cells'].cat.reorder_categories(list(lin), ordered=True)

sc.tl.rank_genes_groups(test3_endo2, 'Subpopulation of Endothelial Cells', method='wilcoxon', pts=True, key_added='Subpopulation of Endothelial Cells_rank_genes_groups')

sc.pl.rank_genes_groups_heatmap(test3_endo2, n_genes=15, groupby='Subpopulation of Endothelial Cells', key='Subpopulation of Endothelial Cells_rank_genes_groups', groups=['EC_1', 'EC_2', 'EC_3', 'EC_4'], show_gene_labels=True, min_logfoldchange=1, dendrogram=False, cmap=cmap, use_raw=False, vmin=-1.5, vmax=1.5, swap_axes=True, show=False, var_group_rotation=0)
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









# Research Letter Figures
colormap = {'EC_1': '#8dd3c7',
            'EC_2': '#80b1d3',
            'EC_3': '#fccde5',
            'EC_4': '#bebada'}
fig, axes = plt.subplots(1,3)
sns.despine()
sc.pl.umap(test3_endo2[test3_endo2.obs['batch'] == 'm01'], color='Subpopulation of Endothelial Cells', palette=colormap, show=False, legend_loc=None, size=200, title='1 month', ax=axes[0])
sc.pl.umap(test3_endo2[test3_endo2.obs['batch'] == 'm10'], color='Subpopulation of Endothelial Cells', palette=colormap, show=False, legend_loc=None, size=200, title='10 months', ax=axes[1])
sc.pl.umap(test3_endo2[test3_endo2.obs['batch'] == 'm20'], color='Subpopulation of Endothelial Cells', palette=colormap, show=False, legend_loc=None, size=200, title='20 months', ax=axes[2])