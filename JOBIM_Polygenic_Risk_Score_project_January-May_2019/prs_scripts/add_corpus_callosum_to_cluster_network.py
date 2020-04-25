import os
import glob
import pandas as pd

clusters = glob.glob('cluster*.csv')

for cluster in clusters:
	df = pd.read_csv(cluster, sep=',')

	# select columns 
	df1 = df[['name','-log10(pval)']]

	# split 'name' col containg the connexions between two regions
	# into two cols containing one region each 
	df1['ctx_region_1'] = df['name'].str.split(' ').str[0]
	df1['ctx_region_2'] = df['name'].str.split(' ').str[-1]

	# delete 'name' col
	del df1['name']

	# add col 'left' and 'right'
	df1['orientation_1'] = df1['ctx_region_1'].str.split('-').str[0]
	df1['orientation_2'] = df1['ctx_region_2'].str.split('-').str[0]

	# add 'corpus callosum' commissure region between two connected inter-hemispheric regions (L and R)
	df2 = df1[df1.orientation_1 != df1.orientation_2]
	df3 = df1[df1.orientation_1 != df1.orientation_2]
	df4 = df1[df1.orientation_1 == df1.orientation_2]

	df2['ctx_region_2'] = 'CorpusCallosum'
	df3['ctx_region_1'] = 'CorpusCallosum'

	df5 = pd.concat([df2, df3, df4], ignore_index=True)

	# rearrange col order
	df5 = df5[['orientation_1', 'orientation_2', 'ctx_region_1', 'ctx_region_2', '-log10(pval)']]

	# save df as csv
	(basename, extension) = os.path.splitext(cluster)
	df5.to_csv(basename + "_corpus_callosum.tsv", sep='\t', index=False)
	