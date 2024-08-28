# import prodg
from prodg2.source import generator
import pandas as pd
import numpy as np
import scanpy as sc

# import sys
# import os
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# attempts to use ProDG, but fitting is slow and generation does not converge
dg = generator.DataGenerator()
print("instantiated data generator")

Bsub_minmed_PB = sc.read_h5ad("../../data/probe_Bac_scRNAseq_Rosenthal/B subtilis minimal media/filtered_data_maxpool_processed.h5ad")
Bsub_minmed_PB_df = Bsub_minmed_PB.to_df(layer = "counts").T.head()
Bsub_minmed_PB_df.index = range(0, len(Bsub_minmed_PB_df))

dg.fit(X = Bsub_minmed_PB_df)

Bsub_minmed_PB_df_synthetic = dg.generate(Bsub_minmed_PB_df)

print(Bsub_minmed_PB_df_synthetic.head())
Bsub_minmed_PB_df_synthetic.to_csv("prodg-synthetic-empirical_covar-filtered_feature_bc_matrix.csv", index = False)

