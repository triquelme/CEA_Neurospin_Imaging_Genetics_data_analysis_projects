import os
import glob
import pandas as pd
import csv

INDIR = "/neurospin/ukb/GENETIC_DATA_500k/IMPUTED_DATA/"
INDIR2 = "/neurospin/ukb/GENETIC_DATA_500k/IMPUTED_DATA/plink_extracted_25k_maf0.01/"
OUTDIR = "/home/tr256641/Documents/dMRI_project/data/"

fnames = []
df_result = pd.DataFrame()

with open(OUTDIR + "most_highly_associated_SNP-dMRI_phenotypes.rsid") as f:
	rsid = f.read().splitlines()

for ch in range(1,23):
    fname = INDIR + "ukb_mfi_chr" + str(ch) + "_v3.txt"
    fnames.append(fname)

for f in fnames:
	df = pd.read_csv(f, sep='\t', header=None)
	df2 = df[df[1].isin(rsid)][[1,5]]
	df_result = df_result.append(df2)

output = OUTDIR + "rsID_MAF"
df_result.to_csv(output, sep='\t', index=False, header=False)

