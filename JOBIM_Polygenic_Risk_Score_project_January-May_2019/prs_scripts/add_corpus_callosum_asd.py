import pandas as pd

df = pd.read_csv("hcp_asd_ctx_connectivity_FDR_0.005.tsv", sep='\t')

df1 = df[df.orientation1 != df.orientation2]
df2 = df[df.orientation1 != df.orientation2]
df3 = df[df.orientation1 == df.orientation2]

df1.ctx_region_2 = 'CorpusCallosum'
df2.ctx_region_1 = 'CorpusCallosum'

df4 = pd.concat([df1, df2, df3], ignore_index=True)
df4.to_csv("hcp_asd_ctx_connectivity_FDR_0.005_corpus_callosum.tsv", sep='\t', index=False)