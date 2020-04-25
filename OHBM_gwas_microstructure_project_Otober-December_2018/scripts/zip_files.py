import os
import glob

INDIR = ("/neurospin/ukb/PLINK_OUTPUT/25k_diffusion_TBSS/imp_variants/"
		"concatenated/")

files = glob.glob(INDIR + "*")

for f in files:
	# zip
	os.system("gzip " + f)
	# unzip
	# os.system("gzip -d " + f)
