# CellChat_Preparation_EC.py
# source activate scanpy_1.9.1
# ipython --profile=cellchat

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

# Reading in EC h5ad
test3 = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3.h5ad")
test3.X = scipy.sparse.csr_matrix.toarray(test3.layers['scran_log1p'])

# CellChat Input

leiden_to_celltype_dict = {'0': 'vSMC1',
'1': 'vSMC2',
'2': 'vSMC3',
'3': 'FB1',
'4': 'vSMC4',
'5': 'EC1',
'6': 'FB2',
'7': 'EC2',
'8': 'vSMC5',
'9': 'FB3',
'10': 'Bc',
'11': 'Mpahge',
'12': 'Tc',
'13': 'Neuronal'}
test3.obs['celltype'] = test3.obs['leiden_r05'].map(lambda x: leiden_to_celltype_dict[x]).astype('category')

test3.write(filename="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/CellChat/test3_forCellChat.h5ad")