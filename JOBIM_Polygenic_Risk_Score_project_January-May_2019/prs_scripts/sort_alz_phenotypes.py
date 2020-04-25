#!/usr/bin/env python

import os
import pandas as pd 


UKB_RESULTS = "/home/tr256641/Documents/polygenic_risk_score_project/results/ukb"
HCP_RESULTS = "/home/tr256641/Documents/polygenic_risk_score_project/results/hcp"


### filter most significant phenotypes for ukb

ukb_alz = pd.read_csv(os.path.join(UKB_RESULTS, "ukb_alz_best.tsv"), sep="\t")

# sort by pval
ukb_alz = ukb_alz.sort_values("P")

# keep only 10 most significants
ukb_alz = ukb_alz[:10]

# save dataframe
ukb_alz.to_csv(
    os.path.join(UKB_RESULTS,"ukb_alz_10best.tsv"), sep="\t", header=True, index=False)


### filter most significant phenotypes for hcp

hcp_alz = pd.read_csv(os.path.join(HCP_RESULTS, "hcp_alz_best.tsv"), sep="\t")

# sort by pval
hcp_alz = hcp_alz.sort_values("P")

# keep only 10 most significants
hcp_alz_10best = hcp_alz[:10]
# save dataframe
hcp_alz_10best.to_csv(
    os.path.join(HCP_RESULTS,"hcp_alz_10best.tsv"), sep="\t", header=True, index=False)

# select tbbs phenotypes in hcp
hcp_alz_tbss = hcp_alz[(hcp_alz["Measure"] == "tbss") | (hcp_alz["Measure"] == "all_OPENING_scaled")]
# save dataframe
hcp_alz_tbss.to_csv(
    os.path.join(HCP_RESULTS,"hcp_alz_tbss.tsv"), sep="\t", header=True, index=False)

# concatenate 10best and tbss df
df_concat = pd.concat([hcp_alz_10best, hcp_alz_tbss])

# save dataframe
df_concat.to_csv(
    os.path.join(HCP_RESULTS,"hcp_alz_10best_tbss.tsv"), sep="\t", header=True, index=False)








