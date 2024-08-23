import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy.sparse as sps
from scipy.stats import ranksums, spearmanr, kendalltau
import pickle as pkl

import os
import sys

module_path = os.path.abspath(os.path.join('../..'))
if module_path not in sys.path:
    sys.path.append(module_path)

import tools.util_probe as up
import tools.util as ut
import tools.NB_est as nb
import tools.countsplit as cs
import tools.ClusterDE as cd

import importlib

import warnings

warnings.filterwarnings("ignore")

figure_path = "data_summary_figures"
if not os.path.exists(figure_path):
    os.makedirs(figure_path)

dataset_paths = [
    "../../data/probe_Bac_scRNAseq_Rosenthal/B subtilis minimal media/filtered_data_maxpool_processed.h5ad",
    "synthetic-empirical_corr-filtered_data_maxpool_processed.h5ad",
    "synthetic-schaefer_strimmer-filtered_data_maxpool_processed.h5ad"
]

dataset_names = [
    "Bsub_minmed_PB (original)",
    "Synthetic null with empirical correlation",
    "Synthetic null with Schaefer-Strimmer correlation"
]

datasets = [sc.read_h5ad(p) for p in dataset_paths]

summary_stats = {
    "n_cells": [data.X.shape[0] for data in datasets],
    "n_genes": [data.X.shape[1] for data in datasets],
    "zero_counts": [((np.prod(data.X.shape) - data.layers["counts"].getnnz()) / np.prod(data.X.shape)).round(3) for data in datasets],
    "count_max": [np.max(data.layers["counts"].toarray()) for data in datasets],
    "count_95%": [np.percentile(data.layers["counts"].toarray(), 95) for data in datasets],
    "count_99%": [np.percentile(data.layers["counts"].toarray(), 99) for data in datasets],
}

summary_df = pd.DataFrame(summary_stats, index=dataset_names)
summary_df.columns = ["Cells", "Genes", "Zero counts (percentage)", "Maximum count", "95% quantile", "99% quantile"]

summary_df.to_csv("table_e1-synthetic.csv", index = False)
with open("table_e1-synthetic-latex.txt", "w") as summary_df_latex:
    summary_df_latex.write(summary_df.to_latex())
