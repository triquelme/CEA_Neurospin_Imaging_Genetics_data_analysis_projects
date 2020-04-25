import os
import glob
import pandas as pd

INDIR = ("/home/tr256641/Documents/dMRI_project/Data/")

#significant threshold after Bonferroni correction
#(=genomic threshold / nb of samples)
threshold = 5e-8/(48*6)

files = glob.glob(INDIR + "*.sel")

for file in files:
	df = pd.read_csv(file, sep="\t")
	# threshold P
	df2 = df[df["P"] < threshold]
	df2.to_csv(file, sep='\t', header=True, index=False)