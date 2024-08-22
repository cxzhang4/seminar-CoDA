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

Bsub_minmed_PB = sc.read_h5ad("../../data/probe_Bac_scRNAseq_Rosenthal/B subtilis minimal media/filtered_data_maxpool_processed.h5ad")
Bsub_minmed_PB_df = Bsub_minmed_PB.to_df()
Bsub_minmed_PB_df.index = range(len(Bsub_minmed_PB_df))

dg.fit(Bsub_minmed_PB_df)

print(Bsub_minmed_PB_df.marginal)
print(Bsub_minmed_PB_df.marginal["params"])
print(Bsub_minmed_PB_df.marginal["params"][0])

Bsub_minmed_PB_df_synthetic = dg.generate(Bsub_minmed_PB_df)
#
print(Bsub_minmed_PB_df_synthetic.head())
Bsub_minmed_PB_df_synthetic.to_csv("synthetic_data_1_h5ad_maxpoolprocessed.csv", index = False)

#
# df = pd.DataFrame(np.random.randint(0,100,size=(100, 4)), columns=['Sample1', 'Sample2', 'Sample3', 'Sample4'])
#
# # call generator instance
# prodg = generator.DataGenerator()
#
# # Fit the models to the data
# prodg.fit(df)
#
# # Generate new data
# synthetic_data = prodg.generate(df)
#
# # Print the synthetic data
# print(synthetic_data)
#
# synthetic_data.to_csv("toy_synth_data.csv", index = False)