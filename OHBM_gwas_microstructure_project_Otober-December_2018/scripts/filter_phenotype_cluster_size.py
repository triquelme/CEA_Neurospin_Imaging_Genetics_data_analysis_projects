import os
import pandas as pd

INDIR = "/home/tr256641/Documents/dMRI_project/data/"

df = pd.read_csv(INDIR + "reduced_table2_most_highly_associated_SNP-dMRI_phenotypes_annotated.csv", sep='\t')

df1 = df[df['# phenotypes'] >= 6]
# results in 10 rows

df2 = df[df['# phenotypes'] >= 5]
# results in 13 rows

df3 = df[df['# phenotypes'] >= 4]
# results in 23 rows


df1.to_csv(INDIR + "reduced_table_10rows.csv", sep='\t', index=False, header=True)
df2.to_csv(INDIR + "reduced_table_13rows.csv", sep='\t', index=False, header=True)
df3.to_csv(INDIR + "reduced_table_23rows.csv", sep='\t', index=False, header=True)
