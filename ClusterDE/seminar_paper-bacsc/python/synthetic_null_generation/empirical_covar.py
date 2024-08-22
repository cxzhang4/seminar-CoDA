# import prodg
from prodg.source import generator
import pandas as pd
import numpy as np
import scanpy as sc

# import sys
# import os
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))

dg = generator.DataGenerator()
print("instantiated data generator")

Bsub_minmed_PB = sc.read_10x_h5("../../data/probe_Bac_scRNAseq_Rosenthal/B subtilis minimal media/filtered_feature_bc_matrix.h5")
Bsub_minmed_PB_df = Bsub_minmed_PB.to_df()
Bsub_minmed_PB_df.index = range(0, len(Bsub_minmed_PB_df))

dg.fit(Bsub_minmed_PB_df)

print(Bsub_minmed_PB_df.marginal)
print(Bsub_minmed_PB_df.marginal["params"])
print(Bsub_minmed_PB_df.marginal["params"][0])

Bsub_minmed_PB_df_synthetic = dg.generate(Bsub_minmed_PB_df)

print(Bsub_minmed_PB_df_synthetic.head())
Bsub_minmed_PB_df_synthetic.to_csv("synthetic-empirical_covar-filtered_feature_bc_matrix.csv", index = False)