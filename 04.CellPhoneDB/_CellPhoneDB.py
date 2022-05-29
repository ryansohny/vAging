# _CellPhoneDB.py

# source activate cpdb
# ipython

# Get mouse ortholog of human genes using biomaRt
library(biomaRt)

mart <- useMart('ENSEMBL_MART_ENSEMBL')
# > listMarts()
#                biomart                version
# 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 105
# 2   ENSEMBL_MART_MOUSE      Mouse strains 105
# 3     ENSEMBL_MART_SNP  Ensembl Variation 105
# 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 105
mart <- useDataset('hsapiens_gene_ensembl', mart)
# listDatasets(mart) ==> mart에 존재하는 여러 datasets 중에서 'hsapiens_gene_ensembl을 쓰는 것

annotLookup <- getBM(
    mart = mart,
    attributes = c('ensembl_gene_id',
                   'gene_biotype',
                   'external_gene_name',
                   'uniprot_gn_symbol',
                   'uniprot_gn_id',
                   'hgnc_symbol'),
    uniqueRows=TRUE)



# source activate cellphonedb
# /home/phenomata/cpdb/releases/v2.0.0/data
# /data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/CellPhoneDB

#####################################################################################################################################
#####################################################################################################################################

# From https://www.ensembl.org/biomart/martview/e76ba4ad9601a50ce208f4f2a1742513
# Dataset: Ensembl Genes 104 ==> Mouse genes (GRCm39)
# ==> Attributes: Gene stable ID, Gene stable ID version, Transcript stable ID, Gene name, Human gene stable ID, Human gene name and Human homology type
# tsv download and pick only "ortholog_one2one" in Human homology type field
# ==> /data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/CellPhoneDB/Ensembl_biomaRt_hs-mm_ortholog_one2one.txt

## only Ortholog gene_input generation
dbf = open("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/CellPhoneDB/Ensembl_biomaRt_hs-mm_ortholog_one2one.txt", 'r')
db = list(x.strip('\n').split('\t')[3] for x in dbf)
db.remove('Human gene stable ID')
dbf.close()

dfh = open("gene_input.csv", 'r')
rfh = open("gene_input_hs-mm_ortholog.csv", 'w')
rfh.write(dfh.readline())
for i in dfh:
    if i.strip('\n').split(',')[3] in db:
        rfh.write(i)

## Custom Database Generation (i)
cellphonedb database generate \
--user-protein protein_curated.csv \
--user-gene gene_input_hs-mm_ortholog.csv \
--user-complex complex_curated.csv \
--user-interactions interaction_curated.csv \
--user-interactions-only \
--result-path ./CustomDatabase_20211216_1 \
--log-file customdatabase_20211216_1.out

## Custom Database Generation (ii)
cellphonedb database generate \
--user-protein protein_curated.csv \
--user-gene gene_input_hs-mm_ortholog.csv \
--user-complex complex_curated.csv \
--user-interactions interaction_curated.csv \
--user-interactions-only \
--fetch \
--result-path ./CustomDatabase_20211216_2 \
--log-file customdatabase_20211216_2.out

#####################################################################################################################################
#####################################################################################################################################

cellphonedb method statistical_analysis \
--project-name test3_stat \
--output-path test3_stat \
--counts-data gene_name \
--threads 20 \
test3_meta.txt \
test3_forCellPhoneDB_uppercase.txt

cellphonedb method statistical_analysis \
--project-name test3_stat5 \
--output-path test3_stat5 \
--counts-data gene_name \
--threads 20 \
--threshold 0.1 \
--pvalue 0.005 \
test3_meta3.txt \
test3_forCellPhoneDB_uppercase.txt

cellphonedb method statistical_analysis \
--project-name test3_stat2 \
--output-path test3_stat2 \
--counts-data gene_name \
--threads 20 \
test3_meta2.txt \
test3_forCellPhoneDB_uppercase.txt

# threshold 10%
nohup cellphonedb method statistical_analysis \
--project-name test3_stat3 \
--output-path test3_stat3 \
--counts-data gene_name \
--threads 20 \
--threshold 0.1 \
test3_meta.txt \
test3_forCellPhoneDB_uppercase.txt > stat3.out &

# threshold 10% and pvalue < 0.005
nohup cellphonedb method statistical_analysis \
--project-name test3_stat4 \
--output-path test3_stat4 \
--counts-data gene_name \
--threads 20 \
--threshold 0.1 \
--pvalue 0.005 \
test3_meta.txt \
test3_forCellPhoneDB_uppercase.txt > stat4.out &


# threshold 10% and pvalue < 0.005 (Using Custom database)
/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/CellPhoneDB/CustomDatabase_20211216/cellphonedb_user_2021-12-16-17_07.db


## Plotting
# only using means.txt
cellphonedb plot dot_plot \
--means-path means.txt \
--pvalues-path pvalues.txt \
--output-path ./Figure \
--output-name dotplot.pdf \
--verbose

cellphonedb plot heatmap_plot \
../../test3_meta3.txt \
--pvalues-path pvalues.txt \
--output-path ./Figure \
--count-name heatmap_count.pdf \
--log-name heatmap_log_count.pdf \
--count-network-name test_count_network.txt \
--interaction-count-name test_interactions_count.txt


#####################################################################################################################################
#####################################################################################################################################

#export PATH=/data/Projects/phenomata/99.Tools/anaconda3/bin:$PATH
#source activate scanpy_1.8.1

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
import math
import scanpy.external as sce
import scrublet as scr

sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (7,7)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])

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

# Using the normalized, non-log transformed data
#df_expr_matrix = integrated.X # integrated.X ==> 지금 sparse matrix가 아니라 numpy.matrix로 되어있음.
#df_expr_matrix = df_expr_matrix.T
#df_expr_matrix = pd.DataFrame(df_expr_matrix)
#df_expr_matrix.columns = test3.obs.index
#df_expr_matrix.set_index(test3.var.index, inplace=True)

## Generate CellPhoneDB input files (counts)
anndata.AnnData(X=integrated.X, obs=test3.obs, var=test3.var).write(filename="test3_CellPhoneDB.h5ad")

## Generate CellPhoneDB input files (metadata)
# (i)
df_meta = pd.DataFrame(data={'Cell': list(test3.obs.index), 'cell_type': list(test3.obs['celltype'])}) # 나중에는 age 또는 fine clusters (leiden)도 집어넣어보자
# (ii)
df_meta = pd.DataFrame(data={'Cell': list(test3.obs.index), 'cell_type': list(test3.obs['leiden_r05'])})
# (iii)
df_meta = pd.DataFrame(data={'Cell': list(test3.obs.index), 'cell_type': list(test3.obs['celltype2'])})

# (i)
df_meta.set_index('Cell', inplace=True)
df_meta.to_csv("test3_meta.txt", sep='\t')
# (ii)
df_meta.set_index('Cell', inplace=True)
df_meta.to_csv("test3_meta2.txt", sep='\t')
# (iii)
df_meta.set_index('Cell', inplace=True)
df_meta.to_csv("test3_meta3.txt", sep='\t')