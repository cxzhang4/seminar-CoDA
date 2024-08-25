import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.sparse as sps
from scipy.stats import ranksums, spearmanr, kendalltau

import tools.util_probe as up
import tools.util as ut
import tools.NB_est as nb
import tools.countsplit as cs
import tools.ClusterDE as cd

import multiprocessing

import warnings

def main():
    warnings.filterwarnings("ignore")

    data_dir = "../../"
    data_path = data_dir + "data/probe_Bac_scRNAseq_Rosenthal/B subtilis minimal media"
    data_file_name = "filtered_data_maxpool_processed.h5ad"
    data_gene = sc.read_h5ad(data_path + "/" + data_file_name)
    # data_gene = sc.read_10x_h5(data_path + "/" + data_file_name)

    warnings.filterwarnings("ignore")
    nb.estimate_overdisp_nb(data_gene, layer="counts", flavor="statsmod_auto")

    # find optimal correlation matrix scaling
    try:
        xmin, fval, R_est_noscale = cd.select_covariance_scaling(data_gene, cor_cutoff=0.1, min_scale=1, max_scale=2,
                                                                 maxiter=20, rng_seed=1234)
    except:
        _, R_est_noscale = cd.generate_nb_data_copula(data_gene, rng_seed=1234, nb_flavor="statsmod_auto",
                                                      auto_dist=True, correct_var=True, return_R=True, corr_factor=1,
                                                      R_est=None, check_pd=True)
        xmin = 1

    # generate scaling matrix
    cor_orig = cd.schaefer_strimmer(data_gene.layers["counts"].toarray(), use_corr=True)
    factor_cor = (np.abs(cor_orig) > 0.1)
    cf = factor_cor * xmin
    cf[cf == 0] = 1
    np.fill_diagonal(cf, 1)
    print(cf)

    n, p = data_gene.X.shape
    data_null_gen, R_est = cd.generate_nb_data_copula(data_gene, rng_seed=5678, nb_flavor="statsmod_auto",
                                                      auto_dist=True, correct_var=False, return_R=True,
                                                      new_data_shape=(2 * n, p),
                                                      corr_factor=cf, R_est=R_est_noscale, check_pd=False,
                                                      min_nonzero=2,
                                                      R_metric="corr")

    sc.pp.calculate_qc_metrics(data_null_gen)
    data_null_gen.var["var_counts"] = np.asarray(np.var(data_null_gen.X, axis=0)).squeeze()
    data_null_gen.var["mean_counts"] = np.asarray(np.mean(data_null_gen.X, axis=0)).squeeze()

    nb.estimate_overdisp_nb(data_null_gen, flavor="sctransform", seed=1234)
    data_null_gen.layers["counts"] = data_null_gen.X.copy()

    k_opt = data_gene.uns["BacSC_params"]["k_opt"]
    n_neighbors_opt = data_gene.uns["BacSC_params"]["n_neighbors_opt"]
    min_dist_opt = data_gene.uns["BacSC_params"]["min_dist_opt"]
    res_opt = data_gene.uns["BacSC_params"]["res_opt"]

    sc.pp.calculate_qc_metrics(data_null_gen, var_type="genes", percent_top=None, log1p=True, inplace=True)
    sc.pp.normalize_total(data_null_gen, target_sum=None, layer=None)
    data_null_gen.X = sps.csr_matrix(np.log(data_null_gen.X + np.array(data_null_gen.var["nb_overdisp"] / 4)))
    data_null_gen.layers["vst_counts"] = data_null_gen.X.copy()
    sc.pp.scale(data_null_gen, max_value=10, zero_center=True)
    data_null_gen.X[np.isnan(data_null_gen.X)] = 0
    sc.tl.pca(data_null_gen, svd_solver='arpack')
    sc.pp.neighbors(data_null_gen, n_neighbors=n_neighbors_opt, n_pcs=k_opt)
    sc.tl.umap(data_null_gen, neighbors_key="neighbors", min_dist=min_dist_opt, spread=1)

    sc.pl.umap(data_null_gen, color="total_counts", alpha=1, cmap="viridis", title="Null data")

    synthetic_data_file_name = "synthetic-empirical_corr-" + data_file_name
    data_null_gen.write(data_path + "/" + synthetic_data_file_name)

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
