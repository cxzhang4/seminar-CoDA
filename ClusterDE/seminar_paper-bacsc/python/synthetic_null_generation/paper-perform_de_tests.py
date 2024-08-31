# this script was generated from the BacSC tutorial, and modifications were made thereafter
# https://github.com/bio-datascience/BacSC/blob/main/tutorials/BacSC_cluster_analysis_tutorial.ipynb

#!/usr/bin/env python
# coding: utf-8

# # Analysis of clusters in BacSC
#
# To perform differential expression analysis, we use a modified version of the ClusterDE approach (Song et al., 2023). For this, we will first generate a synthetic null dataset, and then perform FDR control on the p-values of a Wilcoxon ranksum test.
#
# This tutorial covers testing each cluster against the union of all other clusters. For testing just two clusters against each other, subset the data to the clusters of interest immediately at the beginning and remove all genes that are not expressed in any of the two clusters. Then apply the rest of the pipeline as shown here

# In[1]:


# These imports and path modifications are only necessary for development
import importlib
#
# import os
# import sys
# module_path = os.path.abspath(os.path.join('..'))
# if module_path not in sys.path:
#     sys.path.append(module_path)


# In[2]:


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

import warnings
warnings.filterwarnings("ignore")

import multiprocessing
import argparse

def main(corr_type: str, n_datasets: int):

    # ## Read data from main pipeline
    #

    # In[3]:


    # Adjust the data path and filename based on your folder structure
    data_path = "../../data/probe_Bac_scRNAseq_Rosenthal/B subtilis minimal media"
    data_gene = sc.read_h5ad(data_path + "/filtered_data_maxpool_processed.h5ad")

    input_dir = "synthetic_null_data" + "/" + corr_type
    output_dir = "de_results" + "/" + corr_type


    # In[4]:


    # sc.pl.umap(data_gene, color="leiden_opt", palette="tab10", title="Sample 3")
    # plt.show()
    # plt.close()


    # ## Generate synthetic null data for ClusterDE
    #
    # ### Estimate distribution family and parameters for every gene
    #
    # Depending on the size of your data, this step might take a while, but is strictly necessary

    # In[5]:


    importlib.reload(nb)
    import warnings
    warnings.filterwarnings("ignore")
    nb.estimate_overdisp_nb(data_gene, layer="counts", flavor="statsmod_auto")


    # ### Find optimal correlation matrix scaling
    #
    # This step can be skipped if a slightly suboptimal null data is sufficient. To speed up the tutorial, the calculated scaling factor and correlation matrix are given below (commented out)

    # In[6]:


    importlib.reload(cd)
    try:
        xmin, fval, R_est_noscale = cd.select_covariance_scaling(data_gene, cor_cutoff=0.1, min_scale=1, max_scale=2, maxiter=20, rng_seed=1234)
    except:
        _, R_est_noscale = cd.generate_nb_data_copula(data_gene, rng_seed=1234, nb_flavor="statsmod_auto",
                                                      auto_dist=True, correct_var=True, return_R=True, corr_factor=1,
                                                      R_est=None, check_pd=True)
        xmin = 1
    print(xmin)


    # In[7]:


    # Use this if the previous code cell was skipped
    # xmin = 2.088637562328564
    # _, R_est_noscale = cd.generate_nb_data_copula(data_gene_01, rng_seed=3456, nb_flavor="statsmod_auto",
    #                                                   auto_dist=True, correct_var=True, return_R=True, corr_factor=1,
    #                                                   R_est=None, check_pd=True)


    # In[8]:


    importlib.reload(cd)

    use_schaefer_strimmer = (corr_type == "schaefer_strimmer")
    print(use_schaefer_strimmer)

    if use_schaefer_strimmer:
        # Generate scaling matrix
        cor_orig = cd.schaefer_strimmer(data_gene.layers["counts"].toarray(), use_corr=True)
        factor_cor = (np.abs(cor_orig) > 0.1)
        cf = factor_cor * xmin
        cf[cf == 0] = 1
        np.fill_diagonal(cf, 1)
    else:
        cf = 1

    print(cf)
    # ### Generate synthetic null data

    # In[9]:
    for i in range(n_datasets):


        importlib.reload(cd)

        data_null_gen = sc.read_h5ad(input_dir + "/synthetic-filtered_data_maxpool_processed-" + str(i + 1) + ".h5ad")
        # Filter out genes in the original data that were sampled as all zeros in the null data
        data_gene_nonzero = data_gene[:, data_null_gen.var_names].copy()


        # ## Differential expression testing
        #
        # Now we perform the actual DE testing.
        #
        # ### Find a good split for the null data
        #
        # This is a modification from ClusterDE - we randomly cluster the synthetic data in two halves multiple times and select the clustering that has the largest minimal p-value. This step again takes a moment, but is necessary

        # In[20]:


        n_splits = 1
        rng = np.random.default_rng(1234)
        seeds = rng.choice(1000, n_splits)
        null_pvals_dict = {}
        min_pvals_null = []
        c = 0
        res_start = 0.01
        res_step = 0.01

        for s in seeds:
            res2 = res_start
            twoclust = False
            was_greater = False
            max_res_1_cluster = 0.01

            while twoclust is False:
                sc.tl.leiden(data_null_gen, resolution=res2, key_added=f"leiden_{c}", random_state=s)
                nclust = len(data_null_gen.obs[f"leiden_{c}"].unique())
                print(f"resolution: {res2}, clusters: {nclust}")

                if nclust == 2:
                    twoclust = True
                    break
                elif nclust < 2:
                    if res2 > max_res_1_cluster:
                        max_res_1_cluster = res2
                    else:
                        res_step = res_step/2
                    if was_greater:
                        res2 += res_step
                    else:
                        res2 += 5*res_step
                else:
                    was_greater = True
                    min_res_2plus_clusters = res2
                    res2 -= res_step
                res2 = np.round(res2, 15)


            X_null_gen_0 = data_null_gen.X[data_null_gen.obs[f"leiden_{c}"] == "0"]
            X_null_gen_1 = data_null_gen.X[data_null_gen.obs[f"leiden_{c}"] != "0"]
            null_pvals = ranksums(X_null_gen_0, X_null_gen_1, alternative="two-sided").pvalue
            null_pvals_dict[c] = null_pvals

            clusters = data_gene.obs["leiden_opt"].unique()

            min_pvals_null.append(np.min(null_pvals))

            print(f"split {c+1}/{n_splits} - Resolution {res2}")

            res_start = res2 - 5*res_step
            c += 1

        best_split = np.where(min_pvals_null == np.max(min_pvals_null))[0][0]
        print(f"Best split: No. {best_split} - seed: {seeds[best_split]} - minimal p-value: {min_pvals_null[best_split]}")


        # In[21]:


        sc.pl.umap(data_null_gen, color=f"leiden_{best_split}", alpha=1, cmap="viridis", title="Null data")
        plt.show()
        plt.savefig(output_dir + "/synthetic-leiden_best-" + str(i + 1) + ".png")
        plt.close()

        # In[22]:


        data_null_gen.obs["leiden_best"] = data_null_gen.obs[f"leiden_{best_split}"]

        data_null_gen.write(output_dir + "/synthetic-opt-filtered_data_maxpool_processed-" + str(i + 1) + ".h5ad")


        # ### DE testing with FDR control
        #
        # Now we actually do the DE testing for every cluster against the rest of the population. You can adjust the FDR here or subset the respective tables in pvals_log_gen afterwards

        # In[23]:


        importlib.reload(cd)
        clusters = data_gene.obs["leiden_opt"].unique()
        DEs_log_gen = {}
        pvals_log_gen = {}
        fdr = 0.05
        rng = np.random.default_rng(1234)

        for c in clusters:

            X_data_0 = data_gene_nonzero.X[data_gene_nonzero.obs["leiden_opt"] == c].copy()
            X_data_1 = data_gene_nonzero.X[data_gene_nonzero.obs["leiden_opt"] != c].copy()

            n_cells_0 = X_data_0.shape[0]
            n_cells_1 = X_data_1.shape[0]

            X_null_gen_0 = data_null_gen.X[data_null_gen.obs[f"leiden_best"] != "0"]
            X_null_gen_0 = X_null_gen_0[rng.integers(X_null_gen_0.shape[0], size=n_cells_0),:]
            X_null_gen_1 = data_null_gen.X[data_null_gen.obs[f"leiden_best"] == "0"]
            X_null_gen_1 = X_null_gen_1[rng.integers(X_null_gen_1.shape[0], size=n_cells_1),:]
            null_pvals = ranksums(X_null_gen_0, X_null_gen_1, alternative="two-sided").pvalue

            pvals_data = ranksums(X_data_0, X_data_1, alternative="two-sided").pvalue
            p_data = pd.DataFrame({"pval_data": pvals_data}, index=data_gene_nonzero.var.index)
            pval_null_gen = pd.DataFrame({"pval_null": null_pvals}, index=data_null_gen.var.index)

            DE_TU, pval_TU = cd.call_de(p_data, pval_null_gen, FDR=fdr, correct=False, nlog=True)
            data_gene_nonzero.var[f"pval_cluster_{c}_gen"] = pval_TU["pval_data"]
            data_gene_nonzero.var[f"q_cluster_{c}_gen"] = pval_TU["q"]
            data_gene_nonzero.var[f"DE_cluster_{c}_gen"] = (data_gene_nonzero.var[f"q_cluster_{c}_gen"] < fdr)

            DEs_log_gen[c] = DE_TU
            pvals_log_gen[c] = pval_TU
            print(f"Cluster {c} - DE genes: {len(DEs_log_gen[c])}; Minimum q value: {np.min(pvals_log_gen[c]['q'])}")


        # Diagnostic plots - The distribution of contrast scores should have a peak around 0 and should be approx. symmetric

        # In[24]:


        c = "4"

        fig, ax = plt.subplots(2, 3, figsize=(12,6))
        sns.histplot(pvals_log_gen[c], x="pval_trafo_data", ax=ax[0,0], log_scale=True)
        ax[0,0].set_title("Target data p-values (log-transformed)")
        sns.histplot(pvals_log_gen[c], x="pval_trafo_null", ax=ax[0,1], log_scale=True)
        ax[0,1].set_title("Null data p-values (log-transformed)")

        sns.histplot(pvals_log_gen[c], x="cs", ax=ax[0,2])
        ax[0,2].set_title("Contrast scores (with log-transformation)")
        ax[0,2].set(xscale="symlog", ylim=(0, 50))

        sns.histplot(pvals_log_gen[c], x="pval_data", ax=ax[1,0], bins=100)
        ax[1,0].set_title("Target data p-values")
        sns.histplot(pvals_log_gen[c], x="pval_null", ax=ax[1,1], bins=100)
        ax[1,1].set_title("Null data p-values")

        sns.histplot(pvals_log_gen[c], x="cs", ax=ax[1,2], bins=100)
        ax[1,2].set_title("Contrast scores (with log-transformation)")


        plt.tight_layout()
        plt.show()


        # In pvals_log_gen, every cluster has a table that shows pvalues, q-values (which can be thresholded to get gene subsets at certain FDR levels), ... for each gene.

        # In[25]:


        pvals_log_gen["0"]


        # Attach DE results to the data object and write to disk

        # In[26]:


        data_gene_nonzero.uns["ClusterDE_results"] = pvals_log_gen

        data_gene_nonzero.uns["ClusterDE_results"]['params'] = {'groupby': 'leiden_opt',
          'reference': 'rest',
          'use_raw': False,
          'layer': None,
        }


        # In[27]:


        data_gene_nonzero.write(output_dir + "/de_results-filtered_data_maxpool_processed.h5ad")


        # ## DE plots, ...
        #
        # Here, we find some examples for plots that show the DE genes (or top-ranked genes) for every cluster

        # In[28]:


        # Plot 20 top-ranked genes (lowest q-values) for each cluster
        n_genes = 20
        for c in clusters:
            plot_genes = data_gene_nonzero.uns["ClusterDE_results"][c].iloc[:n_genes,:].index.tolist()
            with plt.rc_context({"figure.figsize": (10, 5)}):
                sc.pl.rank_genes_groups_violin(data_gene_nonzero, gene_names=plot_genes, key="ClusterDE_results", groups=c)


        # In[29]:


        # Show gene name and symbol for top genes per cluster
        all_marker_genes = []
        for c in clusters:
            plot_genes = data_gene_nonzero.uns["ClusterDE_results"][c].iloc[:n_genes,:].index.tolist()
            all_marker_genes += plot_genes
            print(f"Cluster {c}")
            print(data_gene.var.loc[plot_genes, ["Locus tag", "Name", "Symbol"]])


        # In[30]:


        # Heatmap of normalized expression per cluster
        sc.pl.heatmap(data_gene_nonzero, all_marker_genes, groupby='leiden_opt', swap_axes=True, show_gene_labels=True)


        # In[ ]:





if __name__ == '__main__':
    multiprocessing.freeze_support()
    parser = argparse.ArgumentParser(description="generate synthetic null data using the specified estimation procedure for the correlation matrix")
    parser.add_argument("corr_type", type=str, help="one of \"empirical\" or \"schaefer_strimmer\"")
    parser.add_argument("n_datasets", type=int, help="number of synthetic null datasets")
    args = parser.parse_args()
    main(args.corr_type, args.n_datasets)