#!/usr/bin/env python

import sys
import os
import argparse
import textwrap
import pandas as pd
import subprocess
import shlex


# CONSTANTS

DOC = """
This command runs GenoWAP prioritization pipeline on a GWAS file based on
the chosen tissue-specific functional annotation file.

python run_GenoWAP_pipeline.py
    -i /neurospin/brainomics/2019_Riquelme/
       genoWAP_variant_prioritization_september_2018/GWAS_data/STAP/
       UKB_British/covar_SexAgeMRI_10PCSArray_19k.pits_asym_avg_depth_STAP.
       asym_avg_depth.pval
    -o /neurospin/brainomics/2019_Riquelme/
       genoWAP_variant_prioritization_september_2018/prioritized_data/
    -ts /neurospin/brainomics/2019_Riquelme/
        genoWAP_variant_prioritization_september_2018/GenoSuite/
        GenoSkyline_Plus_127_specific_tissues_BedGraph/
        Brain_Angular_Gyrus.bedGraph

"""
GENERAL_TISSUES = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/"
    "GenoSuite/GenoSkyline_9_general_tissues_BedGraph")

SPECIFIC_TISSUES = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/"
    "GenoSuite/GenoSkyline_Plus_127_specific_tissues_BedGraph")

EXTRACT_GENOCANYON_GENERAL_SCORES_SCRIPT = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/"
    "scripts/Script_GenoCanyon10K.R")

EXTRACT_GENOSKYLINE_TISSUE_SPECIFIC_SCORES_SCRIPT = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/"
    "scripts/extractScores.py")

GENOWAP_SCRIPT = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/"
    "GenoSuite/GenoWAP-V1.2/GenoWAP.py")

BEDGRAPH_TO_BIGWIG_SCRIPT = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/"
    "scripts/bedGraphToBigWig")

CHROM_SIZES = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/"
    "scripts/hg19.chrom.sizes")

# EXTRACT_GENOCANYON_GENERAL_SCORES_SCRIPT = (
#     "/neurospin/brainomics/bio_resources/genoWAP/"
#     "prioritization/scripts/Script_GenoCanyon10K.R")

# EXTRACT_GENOSKYLINE_TISSUE_SPECIFIC_SCORES_SCRIPT = (
#     "/neurospin/brainomics/bio_resources/genoWAP/"
#     "prioritization/scripts/extractScores.py")

# GENOWAP_SCRIPT = (
#     "/neurospin/brainomics/bio_resources/genoWAP/"
#     "prioritization/scripts/GenoWAP.py")

# BEDGRAPH_TO_BIGWIG_SCRIPT = (
#     "/neurospin/brainomics/bio_resources/genoWAP/"
#     "prioritization/scripts/bedGraphToBigWig")

# CHROM_SIZES = (
#     "/neurospin/brainomics/bio_resources/genoWAP/"
#     "prioritization/scripts/hg19.chrom.sizes")


def get_cmd_line_args():
    """Return the command-line arguments parser"""
    parser = argparse.ArgumentParser(
        prog="run_GenoWAP_pipeline.py",
        description=textwrap.dedent(DOC),
        formatter_class=argparse.RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i",
        metavar="<path>",
        help="input GWAS file")
    required.add_argument(
        "-o",
        metavar="<path>",
        help="output directory that will contain the prioritized data")
    required.add_argument(
        "-ts",
        metavar="<path>",
        help="path to tissue-specific functional annotation file")

    # Optional argument
    parser.add_argument(
        "-s", "--show",
        help="show list of tissue-specific functional annotation files "
             "and exit",
        action="store_true")
    args = parser.parse_args()
    return args


def show_tissue_annotation_files():
    """Print the names of the 9 general tissue annotation files
    and the 127 specific tissue annotation files and exit
    """
    print("\n9 general tissues:\n")
    for root, dirs, files in os.walk(GENERAL_TISSUES):
        for file_name in files:
            print(file_name)
    # print("\n")
    print("\n127 specific tissues:\n")
    for root, dirs, files in os.walk(SPECIFIC_TISSUES):
        for file_name in files:
            print(file_name)
    exit()


def make_output_directory(gwas_input_file, output_path, tissue):
    """Make output directory following the automated format:
    'input_path'/'input_dir_name'/'gwas_name'/tissue_name'
    and returns it with gwas_file, gwas_name and tissue_name for reuse
    """
    # split file path and file name
    gwas_dir, gwas_file = os.path.split(gwas_input_file)
    # split directory path and directory name
    dir_path, dir_name = os.path.split(gwas_dir)
    # remove extention from file name
    gwas_name, extension = os.path.splitext(gwas_file)
    # remove absolute path from tissue name
    tissue_file = os.path.basename(tissue)
    # remove extention from tissue name
    tissue_name, extension = os.path.splitext(tissue_file)
    # format output directory
    output_dir = (output_path + dir_name + "/" + gwas_name +
                  "/" + tissue_name + "/")
    # check if output directory exists
    if os.path.exists(output_dir) is False:
        # if not, make output directory that will store the results
        os.makedirs(output_dir)
    return gwas_file, gwas_name, tissue_name, output_dir


def convert_gwas_input_file_to_GenoWAP_format(gwas_input_file, formated_gwas_file):
    """Convert GWAS input file into GenoWAP format:
    3 columns: SNP chr | SNP genomic position | SNP pval 
    """
    print("Read GWAS input file...\n")
    # read csv file and save it as a pandas DataFrame
    df = pd.read_csv(gwas_input_file, sep="\t")

    print("STEP1: Convert GWAS input file to GenoWAP readable format...\n")
    # select the 3 columns of interest and save it as another dataframe
    df2 = df[["CHR", "BP", "P"]]

    # write the new dataframe with desired columns into a new csv file
    df2.to_csv(formated_gwas_file, sep='\t', header=False, index=False)


def extract_GenoCanyon_scores(formated_gwas_file, GenoCanyon_scores_out,
                              GenoCanyon_scores_log):
    """Extract GenoCanyon General functional scores for GWAS identified SNPs

    usage: R CMD BATCH --slave '--args input_file output_file' Script.R 
    Rcodeoutput.txt 
    """
    print("STEP2: Extract GenoCanyon General functional scores for GWAS "
          "identified SNPs...")
    print("~ about 10 minutes\n")
    # build command line
    command_line = ("R CMD BATCH --slave '--args " + formated_gwas_file + " "
                    + GenoCanyon_scores_out + "' " +
                    EXTRACT_GENOCANYON_GENERAL_SCORES_SCRIPT + " "
                    + GenoCanyon_scores_log)
    # split command line into a list of arguments using shell-like syntax to
    # determine the correct tokenization
    args = shlex.split(command_line)
    # run the unix command described by args
    subprocess.call(args)


def extract_GenoSkyline_scores(formated_input_file, tissue_scores,
                               GenoSkyline_scores_out, GenoSkyline_scores_log):
    """Extract GenoSkyline tissue specific functional scores for GWAS
    identified SNPs 

    usage: python extractScores.py -o outputName gwasSNPsPosition.csv
    tissueSpeAnnotBed
    """
    print("STEP3: Extract GenoSkyline tissue specific functional scores for "
          "GWAS identified SNPs...")
    print("~about 20 minutes\n")
    # construct command line
    command_line = ("python " +
                    EXTRACT_GENOSKYLINE_TISSUE_SPECIFIC_SCORES_SCRIPT +
                    " -o " + GenoSkyline_scores_out + " " + formated_input_file
                    + " " + tissue_scores)
    args = shlex.split(command_line)

    # run subprocess and redirect its stdout and stderr to a log file
    with open(GenoSkyline_scores_log, "w") as outfile:
        subprocess.call(args, stdout=outfile, stderr=outfile)


def prioritize(formated_gwas_file, GenoCanyon_scores_out,
               GenoSkyline_scores_out, prioritized_scores_out,
               prioritized_scores_log):
    """Compute SNPs prioritization 

    usage: python GenoWAP -o DESTINATION_PATH -a GENERAL_ANNOTATION_PATH
    -ts TISSUE_ANNOTATION_PATH GWAS_DATA_PATH 
    """
    print("STEP4: Compute prioritization...")
    print("~about 5 minutes\n")

    # build command line
    command_line = ("python " + GENOWAP_SCRIPT + " -o " + prioritized_scores_out +
                    " -a " + GenoCanyon_scores_out + " -ts " +
                    GenoSkyline_scores_out + " " + formated_gwas_file)
    args = shlex.split(command_line)

    # run subprocess and redirect its stdout and stderr to a log file
    with open(prioritized_scores_log, "w") as outfile:
        subprocess.call(args, stdout=outfile, stderr=outfile)


def convert_prioritized_scores_to_bedGraph(prioritized_scores_out, bedGraph):
    """Convert prioritized scores to bedGraph file format readable in 
    Integrative Genome Viewer (IGV)

    genoWAP output a result file with 3 columns:
    chr# | SNPposition | score
    (score between 0-1, the closer to 1, the stronger a SNP is prioritized)

    To obtain a bedGraph file readable in IGV: need to double SNP position 
    to obtain a start and stop positions
    bedGraph format:
    chr# | startPosition | endPosition | score 
    """

    # read prioritized_scores and save it as a Pandas dataframe
    df = pd.read_csv(prioritized_scores_out, sep="\t", header=None)
    print("STEP5: Convert prioritized scores output file to bedGraph...\n")
    # insert a third column containing second columns' SNPpositions values to
    # obtain a start and a stop position
    # with loc = insertion index, column = label of newly inserted column
    df.insert(loc=2, column=3, value=df[1])
    df.to_csv(bedGraph, sep='\t', header=False, index=False)


def convert_prioritized_scores_to_bed(prioritized_scores_out, bed):
    """Convert prioritized scores to bed file format readable in 
    Integrative Genome Viewer (IGV)

    genoWAP output a result file with 3 columns:
    chr# | SNPposition | score
    (score between 0-1, the closer to 1, the stronger a SNP is prioritized)

    To obtain a bed file readable in IGV: need to double SNPposition 
    to obtain a start and stop positions, and optionally add rsID
    bed format:
    chr# | startPosition | endPosition | rsID | score
    """

    # read prioritized_scores and save it as a Pandas dataframe
    df = pd.read_csv(prioritized_scores_out, sep="\t", header=None)
    print("STEP6: Convert prioritized scores output file to bed...\n")
    # 1) insert a third column containing second columns' SNPpositions values
    # to obtain a start and a stop position
    # with loc = insertion index, column = label of newly inserted column
    df.insert(loc=2, column=3, value=df[1])

    # 2) add rsID column (args.i['SNP'])
    df2 = pd.read_csv(args.i, sep="\t", usecols=['SNP'])
    df.insert(loc=3, column=4, value=df2['SNP'])
    df.to_csv(bed, sep='\t', header=False, index=False)


def change_chr_num_to_char(bedGraph):
    """Convert chr number (1) to chr character (chr1) in bedGraph first column
    """
    df = pd.read_csv(bedGraph, sep="\t", header=None)
    df[0] = "chr" + df[0].astype(str)
    df[0].replace("chr23", "chrX")
    df[0].replace("chr24", "chrY")
    df.to_csv(bedGraph, sep='\t', header=False, index=False)


def sort_bedGraph(bedGraph):
    """Sort bedGraph with unix sort command:
    sort -k1,1 -k2,2n unsorted.bedGraph > sorted.bedGraph
    """
    file_name, ext = os.path.splitext(bedGraph)
    sorted_bedGraph = file_name + ".sortedBedGraph"
    command_line = "sort -k1,1 -k2,2n " + bedGraph
    args = shlex.split(command_line)
    with open(sorted_bedGraph, "w") as outfile:
        subprocess.call(args, stdout=outfile)
    return sorted_bedGraph


def convert_bedGraph_to_bigWig(sorted_bedGraph, bigWig):
    """convert bedGraph to bigWig"""
    # usage: bedGraphToBigWig in.bedGraph chrom.sizes out.bw
    command_line = [BEDGRAPH_TO_BIGWIG_SCRIPT,
                    sorted_bedGraph, CHROM_SIZES, bigWig]
    subprocess.call(command_line)


def process_and_convert_bedGraph_to_bigWig(bedGraph, bigWig):
    """Do the required steps before converting bedGraph file into bigWig file:
    (1) Convert chr number (1) to chr character (chr1) in bedGraph first column
    (2) Sort bedGraph with the unix sort command:
       sort -k1,1 -k2,2n unsorted.bedGraph > sorted.bedGraph 
    after what bedGraph can be converted into bigWig (3)
    """
    print("STEP7: Convert bedGraph to bigWig...\n")
    # 1) change chr number to chr character
    change_chr_num_to_char(bedGraph)
    # 2) sort bedGraph
    sorted_bedGraph = sort_bedGraph(bedGraph)
    # 3) convert bedGraph to bigWig
    convert_bedGraph_to_bigWig(sorted_bedGraph, bigWig)
    # remove temporary sorted bedGraph
    os.remove(sorted_bedGraph)


def add_rsID_column_to_prioritized_scores(prioritized_scores_out):
    """Insert rsID column to prioritized scores file
    """
    df = pd.read_csv(prioritized_scores_out, sep="\t", header=None)
    df2 = pd.read_csv(args.i, sep="\t", usecols=['SNP'])

    print("STEP8: Insert rsID column to prioritized scores file...\n")
    df.insert(loc=2, column=3, value=df2['SNP'])

    df.to_csv(prioritized_scores_out, sep='\t', header=False, index=False)


def extract_most_prioritized_variants(prioritized_scores_out,
                                      most_prioritized_variants, x):
    """Filter variants with prioritized scores over a stringent cutoff x, 
    and then sort the variants according to their scores in a descending manner
    with greater values first
    """
    print("STEP9: Extract most prioritized variants...\n")
    df = pd.read_csv(prioritized_scores_out, sep="\t", header=None)
    df2 = df[df[3] > x]
    df2 = df2.sort_values(3, ascending=False)
    df2.to_csv(most_prioritized_variants, sep='\t', header=False, index=False)


def rank_variants(prioritized_scores_out, ranked_variants):
    """Filter variants with prioritized scores over a low threshold of 0.1
    and then sort the variants according to their scores in a descending manner
    with greater values first
    """
    print("STEP10: Rank variants according to their prioritization scores...\n")
    df = pd.read_csv(prioritized_scores_out, sep="\t", header=None)
    df2 = df[df[3] > 0.1]
    df2 = df2.sort_values(3, ascending=False)
    df2.to_csv(ranked_variants, sep='\t', header=False, index=False)



###################################  Main  ###################################


# Get user's input, output and tissue arguments
args = get_cmd_line_args()

# If user asks for tissue list, show it
if args.show:
    show_tissue_annotation_files()

# STEP0: Create output directory

gwas_file, gwas_name, tissue_name, output_dir = make_output_directory(
    args.i, args.o, args.ts)


# STEP 1,2,3: preliminary steps before variant prioritization
# need to convert GWAS file to GenoWAP readable format
# need to extract general functional scores in GenoCanyon
# and tissue specific functional scores in GenoSkyline for GWAS SNPs positions

# STEP1: Convert GWAS input file to GenoWAP format

formated_gwas_file = output_dir + "formated_" + gwas_file

convert_gwas_input_file_to_GenoWAP_format(args.i, formated_gwas_file)

# STEP2: Extract GenoCanyon General functional scores for GWAS identified SNPs

GenoCanyon_scores_out = (output_dir + "extracted_GenoCanyon_general_scores_"
                         + gwas_name + ".scores")
GenoCanyon_scores_log = (output_dir + "extracted_GenoCanyon_general_scores_"
                         + gwas_name + ".log")

extract_GenoCanyon_scores(formated_gwas_file, GenoCanyon_scores_out,
                          GenoCanyon_scores_log)

# STEP3: Extract GenoSkyline tissue specific functional scores for GWAS
# identified SNPs

GenoSkyline_scores_out = (output_dir + "extracted_GenoSkyline_" + tissue_name
                          + "_specific_scores_" + gwas_name + ".scores")
GenoSkyline_scores_log = (output_dir + "extracted_GenoSkyline_" + tissue_name
                          + "_specific_scores_" + gwas_name + ".log")

extract_GenoSkyline_scores(formated_gwas_file, args.ts, GenoSkyline_scores_out,
                           GenoSkyline_scores_log)

# STEP4: Compute variant prioritization

prioritized_scores_out = (output_dir + "prioritized_scores_" + gwas_name
                          + ".prioritized_scores")
prioritized_scores_log = (output_dir + "prioritized_scores_" + gwas_name
                          + ".log")

prioritize(formated_gwas_file, GenoCanyon_scores_out, GenoSkyline_scores_out,
           prioritized_scores_out, prioritized_scores_log)

# Remove preliminary files after prioritization
if os.path.exists(formated_gwas_file):
    os.remove(formated_gwas_file)
if os.path.exists(GenoCanyon_scores_out):
    os.remove(GenoCanyon_scores_out)
if os.path.exists(GenoSkyline_scores_out):
    os.remove(GenoSkyline_scores_out)

# STEP 5,6,7: convert prioritized scores file into Integrative Genome Viewer
# readable formats

# STEP5: Convert prioritized scores to bedGraph file readable in IGV
# Advantage of bedGraph file format: display variants as bars
# proportionnal to their prioritized scores in Genome Viewer
# Disadvantage: can not display rsID in Genome Viewer

bedGraph = (output_dir + "prioritized_scores_" + gwas_name + ".bedGraph")

convert_prioritized_scores_to_bedGraph(prioritized_scores_out, bedGraph)

# STEP6: Convert prioritized scores to bed file readable in IGV
# Advantage of bed file format: can display rsID
# Disadvantage: display variant as a fixed size bar in IGV,
# can not display a bar proportionnal to the variant's prioritized score
# so bed and bedgraph are complementary

bed = (output_dir + "prioritized_scores_" + gwas_name + ".bed")

convert_prioritized_scores_to_bed(prioritized_scores_out, bed)

# STEP7: Convert bedGraph to bigWig
# Convert bedGraph (flat file) to bigWif (binaray file)
# Advantage: faster to load in Integrative Genome Viewer (IGV)
# Disadvantage: only display variants when zooming enough in a genome region
# display nothing genome-wide or chromosome-wide (whereas bedGraph do)

file_name, ext = os.path.splitext(bedGraph)
bigWig = file_name + ".bigWig"

process_and_convert_bedGraph_to_bigWig(bedGraph, bigWig)


# STEP 8,9,10: extract information from prioritized scores

# STEP8: Add rsID column to prioritized_scores file

add_rsID_column_to_prioritized_scores(prioritized_scores_out)

# STEP9: Extract most prioritized variants

most_prioritized_variants = (
    output_dir + "prioritized_scores_" + gwas_name
    + ".most_prioritized_variants")

extract_most_prioritized_variants(
    prioritized_scores_out, most_prioritized_variants, 0.5)

# STEP10: Rank variants according to their prioritization scores

ranked_variants = (
    output_dir + "prioritized_scores_" + gwas_name + ".ranked_variants")

rank_variants(prioritized_scores_out, ranked_variants)


# TODO Rank variant according to their p-values

# TODO compare ranked p-val and ranked prioritized scores
