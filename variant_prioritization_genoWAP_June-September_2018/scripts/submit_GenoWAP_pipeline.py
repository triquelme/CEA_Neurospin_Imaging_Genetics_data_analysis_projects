#!/usr/bin/env python

import os
import subprocess
import argparse
import textwrap


# CONSTANTS

DOC = """
This command runs GenoWAP prioritization pipeline on the GWAS files listed in
the chosen input directory and based on the chosen tissue-specific functional
annotation file.

python submit_GenoWAP_pipeline.py
    -i /neurospin/brainomics/2019_Riquelme/
       genoWAP_variant_prioritization_september_2018/GWAS_data/
    -o /neurospin/brainomics/2019_Riquelme/
       genoWAP_variant_prioritization_september_2018/prioritized_data/
    -ts /neurospin/brainomics/2019_Riquelme/
        genoWAP_variant_prioritization_september_2018/GenoSuite/
        GenoSkyline_Plus_127_specific_tissues_BedGraph/
        Brain_Angular_Gyrus.bedGraph

"""
# path to 9 general tissue-specific functional annotation file
GENERAL_TISSUES = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/"
    "GenoSuite/GenoSkyline_9_general_tissues_BedGraph")

# path to 127 high resolution tissue-specific functional annotation file
SPECIFIC_TISSUES = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/"
    "GenoSuite/GenoSkyline_Plus_127_specific_tissues_BedGraph")

RUN_GENOWAP_PIPELINE = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/"
    "scripts/run_GenoWAP_pipeline.py")

# RUN_GENOWAP_PIPELINE = ("/neurospin/brainomics/bio_resources/"
#                         "genoWAP/prioritization/scripts/"
#                         "run_GenoWAP_pipeline.py")

# default parameters
INPUT_DIR = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/GWAS_data/")

OUTPUT_DIR = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/prioritized_data/")

TISSUE_ANNOTATION1 = (
    "/neurospin/brainomics/2019_Riquelme/"
    "genoWAP_variant_prioritization_september_2018/GenoSuite/"
    "GenoSkyline_Plus_127_specific_tissues_BedGraph/"
    "Brain_Angular_Gyrus.bedGraph")


def get_cmd_line_args():
    """Return the command-line arguments parser"""
    parser = argparse.ArgumentParser(
        prog="submit_GenoWAP_pipeline.py",
        description=textwrap.dedent(DOC),
        formatter_class=argparse.RawTextHelpFormatter)
    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i",
        metavar="<path>",
        help="input directory that contains the GWAS data",
        default=INPUT_DIR)
    required.add_argument(
        "-o",
        metavar="<path>",
        help="output directory that will contain the prioritized data",
        default=OUTPUT_DIR)
    required.add_argument(
        "-ts",
        metavar="<path>",
        help="path to tissue-specific functional annotation file",
        default=TISSUE_ANNOTATION1)
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


def get_input_files(input_dir):
    """Return the files contained in the input directory"""
    input_files = []
    # walk through all files in the input directory
    for root, dirs, files in os.walk(input_dir):
        for file_name in files:
            # add absolute path to the name
            file_name = os.path.join(os.path.abspath(root), file_name)
            input_files.append(file_name)
    return input_files


def run_GenoWAP_pipeline_on_multiple_files(input_files, output_dir, tissue):
    """Run GenoWAP pipeline for each file from the input directory
    as parallelized subprocesses
    """
    commands = []
    for file in input_files:
        # usage run_GenoWAP_pipeline.py: python run_GenoWAP_pipeline.py
        # -i gwas_input_file -o output_dir -ts tissue_annotation
        command_line = ["python", RUN_GENOWAP_PIPELINE, "-i", file,
                        "-o", output_dir, "-ts", tissue]
        commands.append(command_line)
        # print("Run GenoWAP pipeline for " + file + " and specific "
        # "tissue: " + tissue + " ...\n")
    procs = [subprocess.Popen(i) for i in commands]
    for p in procs:
        p.wait()


################################## Main ######################################

# Get user's input directory, output directory and tissue annotation file
# arguments
args = get_cmd_line_args()

# If user asks for tissue list, show it
if args.show:
    show_tissue_annotation_files()

# Browse files in the input directory
input_files = get_input_files(args.i)

# Run GenoWAP pipeline for each file
run_GenoWAP_pipeline_on_multiple_files(input_files, args.o, args.ts)
