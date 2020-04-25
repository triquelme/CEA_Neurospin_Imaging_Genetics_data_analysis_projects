#! /usr/bin/env python
# -*- coding: utf-8 -*

##########################################################################
# Copyright (C) CEA, 2017
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import glob
import argparse
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

def _load_plink_gwas(plink_gwas, keep=None):
    """
    Load the result of a Plink GWAS as a DataFrame.
    Any row is removed if p-value is NaN or TEST != "ADD".

    Parameters
    ----------
    plink_gwas: str
        Path to the Plink GWAS output.
    keep: list of str, default None
        Fields in the GWAS output to keep.
        By default ["CHR", "SNP", "TEST", "P"].
    """

    if keep is None:
        keep = ["CHR", "SNP", "TEST", "P"]

    # Load Plink GWAS result
    df = pd.read_csv(plink_gwas, delim_whitespace=True,
                         usecols=["CHR", "SNP", "TEST", "P"])

    # Remove rows where pvalue is NaN
    df.dropna(how="any", inplace=True, axis=0)

    # Drop rows where TEST != "ADD"
    df = df[df["TEST"] == "ADD"]

    return df


def qq_plot(plink_gwas, outdir=None, basename=None,
                          ext="_qqplot.png"):
    """
    Check distribution of p-values for a Plink GWAS and output a PNG.
    P-values are removed if NaN or TEST != "ADD".

    Parameters
    ----------
    plink_gwas: str
        Path to the Plink GWAS output.
    outdir: str, default None
        Directory where output, by default in the input directory.
    basename: str, default None
        Filename without extension of the output snashot.
        Output snapshot is <outdir>/<basename><ext>
    ext: str, default "_qqplot.png"
        Extension added to the basename.
    """
    if outdir is None:
        outdir = os.path.dirname(plink_gwas)

    if basename is None:
        basename = os.path.basename(plink_gwas)

    # Load Plink GWAS result
    df = _load_plink_gwas(plink_gwas=plink_gwas)

    # Compute logarithm of p-values
    observed = list(df["P"])
    median_test_chisq = np.median(stats.chi2.ppf(1-np.asarray(observed), 1))
    median_chisq = stats.chi2.ppf(0.5, 1)
    lambda_median = median_test_chisq/median_chisq
    phen = os.path.basename(plink_gwas[:-13])
    print "Lambda median: "+str(lambda_median)+", "+phen
    
    
    observed.sort()
    observed = [-np.log10(x) for x in observed]

    L = len(observed)
    expected = range(1,L+1)
    L = float(L)
    expected = [-np.log10(float(x)/L) for x in expected]

    plt.figure()
    plt.scatter(expected, observed, marker='o', color='k')
    plt.plot(expected, expected, linestyle='-', color='r')

    plt.xlabel('Expected (-logP)')
    plt.ylabel('Observed (-logP)')
    plt.xlim([0, max(expected)+1])
    plt.ylim([0, max(observed)+1])
    plt.title('Lambda median '+str(lambda_median))
    qqplot = os.path.join(outdir, basename + ext)
    plt.savefig(qqplot)
    plt.close()
    return qqplot


def manhattan_plot_of_plink_gwas(plink_gwas, outdir=None, basename=None,
                                 ext="_manhattan.png", critical_value=0.05,
                                 rs_ids=None):
    """
    Create a Manhattan plot for a Plink GWAS and output a PNG.
    Assuming one file and not one per chromosome.

    P-values are removed if NaN or TEST != "ADD".

    Parameters
    ----------
    plink_gwas: str
        Path to the Plink GWAS output.
    outdir: str, default None
        Path to directory where to output the PNG. By default input directory.
    basename: str, default None
        Output PNG is <outdir>/<basename><ext>.
    ext: str, default "_manhattan.png"
        Extension added to basename.
    critical_value: float, default 0.05
        A horizontal is plot with y=<critical_value>/<# p-values> (Bonferroni).
    rs_ids: list of str, default None
        List of RS IDs for which the marker will be '+' instead of '.'.
    """

    if outdir is None:
        outdir = os.path.dirname(plink_gwas)

    if basename is None:
        basename = os.path.basename(plink_gwas)

    rs_ids = set(rs_ids if rs_ids is not None else [])

    # Load Plink GWAS result
    df = _load_plink_gwas(plink_gwas=plink_gwas)

    # Compute -log10(pvalues)
    df["-log10(P)"] = -np.log10(df["P"])

    # Create Manhattan plot
    fig, ax = plt.subplots()

    # Bonferroni threshold
    bonferroni_threshold = 0.05/df.shape[0]
    minuslog10_threshold = -np.log10(bonferroni_threshold)

    # Set scales
    ax.set_xlim([0, df.shape[0]])
    y_max = max(bonferroni_threshold, df["-log10(P)"].max())
    ax.set_ylim([0, y_max * 1.1])  # Add 10% to the highest y-value

    # List available chromosomes
    chromosomes = df["CHR"].unique()

    # Use different colors for different chromosomes, circle with these
    #colors = itertools.cycle(['c', 'g', 'r', 'm', 'b', 'k'])
    colors = itertools.cycle(['0.25', '0.75'])
    # Abscissae where to put chromosome labels
    xticks = []

    # For each chromosome
    for chrom, color in zip(chromosomes, colors):

        # Extract chromosome data
        subset = df[df["CHR"] == chrom]

        # Abscissa of label
        middle_index = subset.index[subset.shape[0] / 2]
        xticks.append(middle_index)

        # Divide p-values: SNP in/out 'rs_ids' argument
        split = subset["SNP"].isin(rs_ids)
        df_special_marker = subset[split]
        df_normal_marker = subset[~split]

        # Plot p-values
        ax.plot(df_normal_marker.index, df_normal_marker["-log10(P)"],
                ls="", marker=".", c=color, markersize=1.2, markeredgewidth=.7)
        ax.plot(df_special_marker.index, df_special_marker["-log10(P)"],
                ls="", marker="+", c='r', markersize=3, markeredgewidth=1 )

    # Set labels
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(p-value)")
    ax.set_xticks(xticks)
    ax.set_xticklabels(chromosomes, fontsize="xx-small")

    # Plot Bonferroni threshold as an horizontal line
    plt.axhline(y=minuslog10_threshold, color='r', linestyle='-')

    # Save plot
    manhattan_plot = os.path.join(outdir, basename + ext)
    fig.savefig(manhattan_plot, dpi=300)
    plt.close()
    return manhattan_plot


if __name__ == '__main__':
    # Directory containing PLINK output to parse
    dirin = ''
    # Reference GWAS directory
    dir_ref = ''
    outdir = ''

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in',
                        help='path containing .assoc.linear to parse',
                        default=dirin)
    parser.add_argument('-r', '--ref',
                        help='directory containing previously identified SNPs',
                        default=dir_ref)
    parser.add_argument('-o', '--out',
                        help='output directory for PLINK parsed output',
                        default=outdir)
    args = vars(parser.parse_args())

    dirin = args['in']    
    dir_ref = args['ref']
    outdir = args['out']
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    """
    ref_gwas = os.path.join(dir_ref, ".tsv")
    df_ref = pd.read_csv(ref_gwas, sep="\t")
    rs_ids = list(df_ref.SNPS)
    """
    phenos = glob.glob(dirin+'/*.assoc.linear')
    set_parameters = []
    for plink_gwas in phenos:
        qq_plot(plink_gwas, outdir=outdir)
        manhattan_plot_of_plink_gwas(plink_gwas, outdir=outdir)
                                     #, rs_ids = rs_ids)
