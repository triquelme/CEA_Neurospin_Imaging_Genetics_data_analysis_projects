#!/usr/bin/env python

import os
import pandas as pd 


UKB_RESULTS = "/home/tr256641/Documents/polygenic_risk_score_project/results/ukb"
HCP_RESULTS = "/home/tr256641/Documents/polygenic_risk_score_project/results/hcp"
OUTDIR = "/home/tr256641/Documents/polygenic_risk_score_project/results/common_phenotypes"


### search for common phenotypes for ASD

ukb_asd = pd.read_csv(os.path.join(UKB_RESULTS, "ukb_asd_best.tsv"), sep="\t")
hcp_asd = pd.read_csv(os.path.join(HCP_RESULTS, "hcp_asd_best.tsv"), sep="\t")

#select only phenotypes in common between HCP and UKB
hcp_asd = hcp_asd[hcp_asd["Measure"] == "tbss"]
#hcp_asd[(hcp_asd["Measure"] == "tbss") | (hcp_asd["Measure"] == "all_OPENING_scaled")]

# get phenotypes in common between HCP and UKB
asd_common_phenotypes = pd.merge(ukb_asd, hcp_asd, on=['Phenotype'], how='inner')

# save just phenotype names
asd_common_phenotypes['Phenotype'].to_csv(
    os.path.join(OUTDIR,"asd_common_phenotypes.txt"), sep="\t", header=True, index=False)

# save phenotypes names + values
asd_common_phenotypes.to_csv(
    os.path.join(OUTDIR,"asd_common_phenotypes.csv"), sep="\t", header=True, index=False)


### search for common phenotypes for ALZ

ukb_alz = pd.read_csv(os.path.join(UKB_RESULTS, "ukb_alz_best.tsv"), sep="\t")
hcp_alz = pd.read_csv(os.path.join(HCP_RESULTS, "hcp_alz_best.tsv"), sep="\t")

#select only phenotypes in common between HCP and UKB
hcp_alz = hcp_alz[hcp_alz["Measure"] == "tbss"]
#hcp_alz[(hcp_alz["Measure"] == "tbss") | (hcp_alz["Measure"] == "all_OPENING_scaled")]

# get phenotypes in common between HCP and UKB
alz_common_phenotypes = pd.merge(ukb_alz, hcp_alz, on=['Phenotype'], how='inner')

# save just phenotype names
alz_common_phenotypes['Phenotype'].to_csv(
    os.path.join(OUTDIR,"alz_common_phenotypes.txt"), sep="\t", header=True, index=False)

# save phenotypes names + values
alz_common_phenotypes.to_csv(
    os.path.join(OUTDIR,"alz_common_phenotypes.csv"), sep="\t", header=True, index=False)

