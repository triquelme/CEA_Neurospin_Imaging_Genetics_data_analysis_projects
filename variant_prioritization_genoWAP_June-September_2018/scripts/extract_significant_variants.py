#!/usr/bin/env python

import os
import pandas as pd 


input_dir = ("/neurospin/brainomics/2018_variant_prioritization/"
"prioritized_data/STAP/UKB_British/")

ouput_dir = ("/neurospin/brainomics/2018_variant_prioritization/"
"prioritized_data/STAP/UKB_British/Most_prioritized_variants_0.5_cutoff/")


def rank_variants_pval(gwas_file, ranked_variants_pval):
    """rank variant according to their p-values"""
    print("STEP11: Rank variants according to their their p-values...\n")
    df = pd.read_csv(gwas_file, sep="\t")
    df = df[df["P"] < 0.00005]
    df = df.sort_values("P")
    df.to_csv(ranked_variants_pval, header=True, index=False)


###################################  Main  ###################################

# check if output directory exists
if os.path.exists(output_dir) is False:
    # if not, creates it
    os.makedirs(output_dir)

#find file with p-val < 5x10-8
ranked_variants_pval = output_dir + gwas_name + ".ranked_variants_pval"

rank_variants_pval(args.i, ranked_variants_pval)

