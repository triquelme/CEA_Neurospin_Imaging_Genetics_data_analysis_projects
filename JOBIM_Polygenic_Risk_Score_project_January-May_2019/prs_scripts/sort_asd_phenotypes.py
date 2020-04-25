#!/usr/bin/env python

import os
import pandas as pd 


UKB_RESULTS = "/home/tr256641/Documents/polygenic_risk_score_project/results/ukb"
HCP_RESULTS = "/home/tr256641/Documents/polygenic_risk_score_project/results/hcp"


### filter most significant phenotypes for ukb

ukb_asd = pd.read_csv(os.path.join(UKB_RESULTS, "ukb_asd_best.tsv"), sep="\t")

# sort by pval
ukb_asd = ukb_asd.sort_values("P")

# keep only 10 most significants
ukb_asd = ukb_asd[:10]

# save dataframe
ukb_asd.to_csv(
    os.path.join(UKB_RESULTS,"ukb_asd_10best.tsv"), sep="\t", header=True, index=False)


### filter most significant phenotypes for hcp

hcp_asd = pd.read_csv(os.path.join(HCP_RESULTS, "hcp_asd_best.tsv"), sep="\t")

# sort by pval
hcp_asd = hcp_asd.sort_values("P")

# keep only 10 most significants
hcp_asd_10best = hcp_asd[:10]
# save dataframe
hcp_asd_10best.to_csv(
    os.path.join(HCP_RESULTS,"hcp_asd_10best.tsv"), sep="\t", header=True, index=False)

# select tbbs phenotypes in hcp
hcp_asd_tbss = hcp_asd[(hcp_asd["Measure"] == "tbss") | (hcp_asd["Measure"] == "all_OPENING_scaled")]
# save dataframe
hcp_asd_tbss.to_csv(
    os.path.join(HCP_RESULTS,"hcp_asd_tbss.tsv"), sep="\t", header=True, index=False)

# concatenate 10best and tbss df
df_concat = pd.concat([hcp_asd_10best, hcp_asd_tbss])

# save dataframe
df_concat.to_csv(
    os.path.join(HCP_RESULTS,"hcp_asd_10best_tbss.tsv"), sep="\t", header=True, index=False)








