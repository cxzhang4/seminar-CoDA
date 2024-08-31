import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy
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
    # "../../data/probe_Bac_scRNAseq_Rosenthal/B subtilis minimal media/synthetic-empirical_corr-filtered_data_maxpool_processed.h5ad",
    # "../../data/probe_Bac_scRNAseq_Rosenthal/B subtilis minimal media/synthetic-schaefer_strimmer-filtered_data_maxpool_processed.h5ad"
    # "../../data/probe_Bac_scRNAseq_Rosenthal/B subtilis minimal media/null_data_opt_tutorial.h5ad"
    "notebook-synthetic-schaefer_strimmer-filtered_data_maxpool_processed-260824.h5ad",
    "notebook-synthetic-empirical_corr-filtered_data_maxpool_processed-260824.h5ad",
    "synthetic_null_data/schaefer_strimmer/synthetic-filtered_data_maxpool_processed.h5ad",
    "synthetic_null_data/empirical/synthetic-filtered_data_maxpool_processed-1.h5ad"
]

dataset_names = [
    "Bsub_minmed_PB (original)",
    "Synthetic null with Schaefer-Strimmer correlation",
    "Synthetic null with empirical correlation",
    "Schaefer-Strimmer, generated recently",
    "Empirical, generated recently"
]

datasets = [sc.read_h5ad(p) for p in dataset_paths]

def get_dense_array(maybe_sparse_matrix):
    if isinstance(maybe_sparse_matrix, scipy.sparse._csr.csr_matrix):
        return maybe_sparse_matrix.toarray()
    return maybe_sparse_matrix

summary_stats = {
    "n_cells": [data.X.shape[0] for data in datasets],
    "n_genes": [data.X.shape[1] for data in datasets],
    "min_seq_depth": [np.min(data.obs["total_counts"]) for data in datasets],
    "max_seq_depth": [np.max(data.obs["total_counts"]) for data in datasets],
    "median_seq_depth": [np.median(data.obs["total_counts"]) for data in datasets],
    "zero_counts": [((np.prod(data.layers["counts"].shape) - np.count_nonzero(get_dense_array(data.layers["counts"]))) / np.prod(data.layers["counts"].shape)).round(3) for data in datasets],
    "count_max": [np.max(get_dense_array(data.layers["counts"])) for data in datasets],
    "count_95%": [np.percentile(get_dense_array(data.layers["counts"]), 95) for data in datasets],
    "count_99%": [np.percentile(get_dense_array(data.layers["counts"]), 99) for data in datasets],
}
summary_df = pd.DataFrame(summary_stats, index=dataset_names)
summary_df.columns = ["Cells", "Genes", "Minimum seq. depth", "Maximum seq. depth", "Median seq. depth", "Zero counts (percentage)", "Maximum count", "95% quantile", "99% quantile"]

summary_df.insert(0, "Dataset", dataset_names)

summary_df.to_csv("table_e1-synthetic_data.csv", index = False)
with open("table_e1-synthetic_data.txt", "w") as summary_df_latex:
    summary_df_latex.write(summary_df.to_latex())
