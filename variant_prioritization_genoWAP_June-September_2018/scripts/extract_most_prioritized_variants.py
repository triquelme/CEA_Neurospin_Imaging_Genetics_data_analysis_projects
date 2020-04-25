#!/usr/bin/env python

import os
import pandas as pd
from shutil import copyfile
import glob

input_dir = (
    "/neurospin/brainomics/2018_variant_prioritization/prioritized_data/STAP/"
    "UKB_not_British/")

output_dir = input_dir + "most_prioritized_variants_0.5_cutoff/"


def get_most_prioritized_variants_files(input_dir):
    """Return the files contained in the input directory"""
    most_prioritized_variants = []
    # walk through all files in the input directory
    for root, dirs, files in os.walk(input_dir):
        for file_name in files:
            # add absolute path to the name
            file = os.path.join(os.path.abspath(root), file_name)
            basename, extension = os.path.splitext(file)
            if extension == ".most_prioritized_variants":
                most_prioritized_variants.append(file)
    return most_prioritized_variants


def get_ranked_variants_files(input_dir):
    """Return the files contained in the input directory"""
    ranked_variants = []
    # walk through all files in the input directory
    for root, dirs, files in os.walk(input_dir):
        for file_name in files:
            # add absolute path to the name
            file = os.path.join(os.path.abspath(root), file_name)
            basename, extension = os.path.splitext(file)
            if extension == ".ranked_variants" and os.stat(file).st_size > 0:
                ranked_variants.append(file)
    return ranked_variants


def extract_most_prioritized_variants_from_file(prioritized_scores_file,
                                                most_prioritized_variants, x):
    """Filter variants with prioritized scores over a stringent cutoff x, 
    and then sort the variants according to their scores in a descending 
    manner with greater values first, write filtered variants in new files
    """
    df = pd.read_csv(prioritized_scores_file, sep="\t", header=None)
    df2 = df[df[3] > x]
    if not df2.empty:
        df2 = df2.sort_values(3, ascending=False)
        df2.to_csv(most_prioritized_variants, sep='\t', header=False,
                   index=False)


def find_files_with_prioritized_scores_over_threshold(ranked_variants_files,
                                                      threshold):
    """Find files with prioritized scores > 0.5 in ranked variants files
    scan ranked files, filter for variant with prioritized_scored > 0.5,
    write filtered variants in new files
    """
    for file in ranked_variants_files:
        file_path, file_name = os.path.split(file)
        file_path, tissue_name = os.path.split(file_path)
        output_name = output_dir + tissue_name + "_" + file_name
        extract_most_prioritized_variants_from_file(
            file, output_name, threshold)


def copy_most_prioritized_variants_files_in_output_dir(
        most_prioritized_variants_files):
    """Find files with prioritized scores > 0.5 in most prioritized variants 
    files and copy them in output_dir
    """
    for file in most_prioritized_variants_files:
        if os.stat(file).st_size > 0:
            file_path, file_name = os.path.split(file)
            file_path, tissue_name = os.path.split(file_path)
            output_name = ouput_dir + tissue_name + "_" + file_name
            copyfile(file, output_name)


###################################  Main  ###################################

# STEP1: Creates output directory

# check if output directory exists
if os.path.exists(output_dir) is False:
    # if not, creates it
    os.makedirs(output_dir)


# STEP2: Find files with prioritized scores > 0.5 in ranked variants files
# write variants with prioritized scores > 0.5 in new files

ranked_variants_files = get_ranked_variants_files(input_dir)

print("Find files with prioritized scores > 0.5 in ranked variants files...\n")

find_files_with_prioritized_scores_over_threshold(ranked_variants_files, 0.5)


# STEP2 alternative:
# find files with prioritized scores > 0.5 in most_prioritized_variants_files
# and copy them in output_dir

# most_prioritized_variants_files = get_most_prioritized_variants_files(input_dir)

# copy_most_prioritized_variants_files_in_output_dir(
#     most_prioritized_variants_files)


# STEP3: list most prioritized variants across all files
# 1) open files' content into a common dataframe
# 2) sort descending

all_files = glob.glob(output_dir + "*.ranked_variants")

df = pd.concat((pd.read_csv(file, sep="\t", header=None)
                for file in all_files))

df = df.sort_values(3, ascending=False)

# most_prioritized_variants_summary = (
#     output_dir + "most_prioritized_variants_summary")

# df.to_csv(most_prioritized_variants_summary, sep="\t", header=False,
#            index=False)

# STEP4: merge duplicates and mean their scores
# second filtering with more stringent threshold
# in order not to lower the mean of highly prioritized variant 
# with their lower score in another tissue

df = df[df[3] > 0.7]
df = df.groupby([0,1,2]).mean().reset_index()
df = df.sort_values(3, ascending=False)

most_prioritized_variants_merged_duplicates = (
    output_dir + "most_prioritized_variants_merged_duplicates")

df.to_csv(most_prioritized_variants_merged_duplicates, sep="\t", header=False,
           index=False)


# TODO change threshold to 0.5 or higher
# do the same with pval: find file with p-val < 5x10-8

# list most prioritized variants across all files
# list most significant variants across all files

# find variants in common
# find new prioritized variants
