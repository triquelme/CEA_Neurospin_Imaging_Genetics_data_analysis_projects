import os
import glob
import pandas as pd
import numpy as np

INDIR = '/neurospin/ukb/PLINK_OUTPUT/25k_diffusion_TBSS/imp_variants/concatenated/'
OUTDIR = "/home/tr256641/Documents/dMRI_project/data/"
BIM_DIR = "/neurospin/ukb/tmp/THOMAS/"
LD_DIR = "/neurospin/tmp/yleguen/formatted/"
LD_THRESHOLD = 0.1
BONFERRONI_THRESHOLD = 1.736e-10

# # 1) Format .sel files containing plink association results
# # filtered for SNPs with pval over genomic threshold corrected with Bonferroni
# # format to phenogram tool recognizable input file
# # concatenate them into one file

# fnames = glob.glob(INDIR + "*.sel")
# fnames2 = []
# for fname in fnames:
#     df = pd.read_csv(fname, sep="\t")

#     # select only necessary columns
#     df2 = df[["CHR", "SNP", "BP", "P"]]

#     # change col_name "BP" to "POS"
#     df2 = df2.rename(columns={"BP": "POS"})

#     # add column PHENOTYPE
#     basename = os.path.basename(fname)
#     basename, extension = os.path.splitext(basename)
#     basename = basename.replace("covar_SexAgeMRI_10PCSArray_25k.all_", "")
#     phenotype = basename.split(".", 1)[1]
#     df2["PHENOTYPE"] = phenotype

#     # add column group
#     # strip basename to get only measure type (FA, ICVF, MO...)
#     measure_type = basename.split("_", 1)[0]
#     df2["GROUP"] = measure_type

#     # save formatted df to csv
#     output = os.path.join(OUTDIR, os.path.basename(fname))
#     df2.to_csv(output, sep='\t', header=True, index=False)
#     fnames2.append(output)

# # concatenate all files in one
# df = pd.concat(pd.read_csv(fname, sep='\t') for fname in fnames2)
# # order snp by group, chr, position
# df = df.sort_values(by=['GROUP', 'CHR', 'POS'])
# output = (OUTDIR + "concatenated.sel")
# df.to_csv(output, sep='\t', header=True, index=False)


# # 2) Reduce numerous SNPs in a same region to one representative

# # 2-1) Select top SNPs and remove the ones in ld with them
# # So the goal here is to only keep top snps not in ld with others for each group and each tract
# # NB for now we assume the top SNP in a LD region will be the same for various phenotypes
# # if that's not the case we will adjust this maybe to merge the top SNPs
# # as ones

# df = pd.read_csv(output, sep="\t")
# df.index = df.SNP
# # DataFrame that will contain all the data kept
# df_result = pd.DataFrame()

# for ch in df['CHR'].unique():
#     # Work on one chromosome at the time
#     df0 = df.loc[df.CHR == ch]
#     input_bim = (BIM_DIR + "chr" + str(ch) + "_imp_filtered.bim")
#     bim = pd.read_csv(input_bim, sep='\t', header=None)

#     # Load LD matrix
#     input_ld = (LD_DIR + "ld_matrix_ch" + str(ch) + ".ld")
#     ld = pd.read_csv(input_ld, sep=' ', header=None)
#     ld = ld.drop(len(ld.columns) - 1, 1)
#     # Identified rows and columns in this LD matrix
#     ld.index = bim[1]
#     ld.columns = bim[1]

#     # Then let's work feature by feature
#     for gr in df0['GROUP'].unique():
#         df1 = df0.loc[df0.GROUP == gr]
#         # Then work tract by tract
#         for ph in df1['PHENOTYPE'].unique():
#             df2 = df1.loc[df1.PHENOTYPE == ph]
#             # Sort pvalues
#             df2 = df2.sort_values('P', ascending=True)
#             df2.index = df2.SNP
#             count = 0
#             # Then loop through snps and keep the one not in ld
#             while len(list(df2.SNP)) > count:
#                 # Select current top snp
#                 snp = list(df2.SNP)[count]
#                 for snp2 in list(df2.SNP)[count + 1:]:
#                     # Weirdly there are duplicates SNP in bim
#                     # This is probably due to the conversion
#                     # from .bgen to .bed/.bim/.fam format
#                     # for now we ignore this problem that should be dealt
#                     # with by looking at the reference that should be
#                     # different for each line
#                     """
#                     # if unique snp this lines are sufficient
#                     if ld.loc[snp][snp2] > LD_THRESHOLD:
#                         # Drop SNP in LD
#                         df2 = df2.drop(snp2)
#                     """
#                     # Normal case SNP is not duplicated
#                     if ld.loc[snp][snp2].size == 1:
#                         if ld.loc[snp][snp2] > LD_THRESHOLD:
#                             # Drop SNP in LD
#                             df2 = df2.drop(snp2)
#                     # Unusual case SNP is duplicated
#                     elif max(ld.loc[snp][snp2]) > LD_THRESHOLD:
#                         # print ld.loc[snp][snp2].size
#                         # print ld.loc[snp][snp2]
#                         # Check that we did not already removed it
#                         if snp2 in df2.index:
#                             # Drop SNP in LD
#                             df2 = df2.drop(snp2)
#                 count += 1
#             df_result = df_result.append(df2)

# df_result = df_result.sort_values(by=['CHR', 'POS'])
# df_result = df_result.loc[df_result.P < BONFERRONI_THRESHOLD]
# output = os.path.join(OUTDIR, 'ld_snp_removed_result.sel')
# df_result.to_csv(output, sep='\t', index=False, header=True)


# # 2-2) Second filter to reduce to one snp in a region of 200000 bp
# # reduce to min pval in a region of 200000 bp

# #df = pd.read_csv(LD_DIR + 'final_result_bonf_corrected.sel', sep="\t")
# #df = pd.read_csv(OUTDIR + 'ld_snp_removed_result.sel', sep="\t")

# df_result = pd.DataFrame()

# for ch in df.CHR.unique():
#     df2 = df.loc[df.CHR == ch]
#     for lab1, row1 in df2.iterrows():
#         for lab2, row2 in df2.iterrows():
#             if lab2 > lab1:
#                 if (row2.POS == row1.POS) or (row2.POS < (row1.POS + 200000)):
#                     if row2.P > row1.P:
#                         df2 = df2.drop(lab2)
#     df_result = df_result.append(df2)

# df_result = df_result.sort_values(by=['CHR', 'POS'])
# output = os.path.join(OUTDIR, 'ld_snp_removed_window_corrected_result.sel')
# df_result.to_csv(output, sep='\t', index=False, header=True)


# 3) Split dataframe into one with ICVF,ISOVF,MD and another with FA,MO,OD

df = pd.read_csv(
    OUTDIR + 'final_result_bonf_corrected_LD0.001_unique_200kBP_phenogram.csv', sep="\t")
# remove last column with NaN
df = df.drop('Unnamed: 6', 1)

df_IIO = df[np.logical_or(np.logical_or(df['GROUP'] == 'ICVF', df[
    'GROUP'] == 'ISOVF'), df['GROUP'] == 'OD')]

# df_FMM = df[np.logical_or(np.logical_or(df['GROUP'] == 'FA', df[
#     'GROUP'] == 'MO'), df['GROUP'] == 'MD')]

df_IIO.to_csv(
    OUTDIR + "final_result_bonf_corrected_LD0.001_unique_200kBP_phenogram_ICVF_ISOVF_OD.sel",
    sep='\t', header=True, index=False)

# df_FMM.to_csv(
#     OUTDIR + "final_result_bonf_corrected_LD0.001_unique_200kBP_phenogram_FA_MO_MD.sel",
#     sep='\t', header=True, index=False)


# 4) Split dataframe into ICVF, OD and FA

df = pd.read_csv(
    OUTDIR + 'final_result_bonf_corrected_LD0.001_unique_200kBP_phenogram.csv', sep="\t")
# remove last column with NaN
df = df.drop('Unnamed: 6', 1)

df_ICVF = df[df['GROUP'] == 'ICVF']
df_OD = df[df['GROUP'] == 'OD']
df_FA = df[df['GROUP'] == 'FA']
df_ISOVF = df[df['GROUP'] == 'ISOVF']

df_ICVF.to_csv(OUTDIR + "final_result_bonf_corrected_LD0.001_unique_200kBP_phenogram_ICVF.sel",
               sep='\t', header=True, index=False)
df_OD.to_csv(OUTDIR + "final_result_bonf_corrected_LD0.001_unique_200kBP_phenogram_OD.sel",
               sep='\t', header=True, index=False)
df_ISOVF.to_csv(OUTDIR + "final_result_bonf_corrected_LD0.001_unique_200kBP_phenogram_ISOVF.sel",
               sep='\t', header=True, index=False)
df_FA.to_csv(OUTDIR + "final_result_bonf_corrected_LD0.001_unique_200kBP_phenogram_FA.sel",
               sep='\t', header=True, index=False)

# #ideas

# df = pd.read_csv(output, sep="\t")
# groups = df['GROUP'].unique()
# chromosomes = df['CHR'].unique()
# chromosomes.sort()

# for g in groups:
#     for c in chromosomes:
#         positions = df['POS'][df['GROUP'] == g][df['CHR'] == c]
#         # print(g, c)
#         for p in positions:
#             limit = p - 300000
#             # if other snp exist before this snp in a window of 200,000 bp:
#             if len(df[df['GROUP'] == g][df['CHR'] == c][df['POS'] > limit][df['POS'] < p]) > 0:
#                 # delete snp at this position from the dataframe
#                 df = df[df['POS'] != p]

# def region_min_pval(df, position):
#     lower_limit = position - 100000
#     upper_limit = position + 100000
#     region_snp_pval = []
#     for snp in SNP:
#         if df[np.logical_and(df['POS'] > low_limit, df['POS'] < high_limit)]:
#             region_snp_pval.append(df[df['POS'] == position]['P'])
#     min = min(region_snp_pval)


# # for loop over group:
# for group in GROUP:
#     # for loop over chr:
#     for chromosome in CHR:
#         # for loop over positions inside a chr:
#         for position in POS:
#             # for 1st snp if other snp position inside window of 200,000 bp:
#             region_max(df, position)
#             if df[df['POS'] == position]['P'] != region_max:
#                 # keep max p-val, remove others
#                 df = df[df.POS != 'position']

# # "pandas.rolling" compute sliding window
# df['P'].rolling(60, center=True).max()
# df['POS'].rolling(30, center=True).mean()
# # TODO: replace mean() by custom function
# df['POS'].rolling(200000,).custom()
# # groupby
# df2 = df.groupby(['GROUP', 'CHR'])
# for name, group in df2:
#     print name
#     print group['POS'] - 200000
