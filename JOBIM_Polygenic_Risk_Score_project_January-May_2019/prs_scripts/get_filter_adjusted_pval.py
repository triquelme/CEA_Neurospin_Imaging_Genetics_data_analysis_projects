import pandas as pd 
import statsmodels.stats.multitest as multi

df = pd.read_csv("ukb_alzheimer_all_OPENING_scaled.summary", sep='\t')

df = df[['Phenotype','P']]

df = df.sort_values('P')

df = df.reset_index(drop=True)

# adjust pval with Bonferroni
df['P_Bonferroni'] = multi.multipletests(df.P, alpha=0.05, method='bonferroni')[1]

# adjust pval with FDR
df['P_FDR'] = multi.multipletests(df.P, alpha=0.05, method='fdr_bh')[1]

df.to_csv("ukb_alzheimer_all_OPENING_scaled_p_adjusted.tsv", sep='\t', index=False)

# filter p_val < 0.01
df2 = df[df.P < 0.01]

df2.to_csv("ukb_alzheimer_all_OPENING_scaled_p_adjusted_p_filter_0.01.tsv", sep='\t', index=False)