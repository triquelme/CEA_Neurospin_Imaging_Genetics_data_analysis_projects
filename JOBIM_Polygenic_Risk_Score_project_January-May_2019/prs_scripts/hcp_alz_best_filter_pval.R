setwd("/home/tr256641/Téléchargements/JOBIM_abstract/prs_results_hcp/cytoscape_data/")

# load tsv 
hcp_alz = read.csv("hcp_alz_best_cytoscape.csv", sep="\t")

# select columns of interest
hcp_alz = hcp_alz[,c("Phenotype","P","X.log.P.")]
names(hcp_alz) = c("ctx_connexion", "pval", "-log10(pval)")

# split column of connexions between 2 ctx regions into 2 columns of ctx_regions and save as a new df
ctx_regions <- data.frame(do.call('rbind', strsplit(as.character(hcp_alz$ctx_connexion),'_',fixed=TRUE)))
names(ctx_regions) <- c("ctx_region_1","ctx_region_2")
ctx_regions$ctx_region_1 <- gsub('ctx-rh', 'R', ctx_regions$ctx_region_1)
ctx_regions$ctx_region_1 <- gsub('ctx-lh', 'L', ctx_regions$ctx_region_1)
ctx_regions$ctx_region_1 <- gsub('Right', 'R', ctx_regions$ctx_region_1)
ctx_regions$ctx_region_1 <- gsub('Left', 'L', ctx_regions$ctx_region_1)
ctx_regions$ctx_region_2 <- gsub('ctx-rh', 'R', ctx_regions$ctx_region_2)
ctx_regions$ctx_region_2 <- gsub('ctx-lh', 'L', ctx_regions$ctx_region_2)
ctx_regions$ctx_region_2 <- gsub('Right', 'R', ctx_regions$ctx_region_2)
ctx_regions$ctx_region_2 <- gsub('Left', 'L', ctx_regions$ctx_region_2)

# bind the 2 dfs
hcp_alz = hcp_alz[,c("pval", "-log10(pval)")]
hcp_alz = cbind(ctx_regions, hcp_alz)

# create orientation Left-right column
orientation1 <- sapply(strsplit(hcp_alz$ctx_region_1,"-"), `[`, 1)
orientation2 <- sapply(strsplit(hcp_alz$ctx_region_2,"-"), `[`, 1)
hcp_alz = cbind(orientation1, orientation2, hcp_alz)

# sort by pval
hcp_alz = hcp_alz[order(hcp_alz$pval),]

# fdr corrected pval in a new column
hcp_alz_FDR = p.adjust(p=hcp_alz$pval, method = "fdr")

# bind this new column to the first table
hcp_alz = cbind(hcp_alz, FDR = hcp_alz_FDR)

# filter FDR under 0.05
hcp_alz_FDR_0.05 = hcp_alz[hcp_alz$FDR<0.05,]

# filter pval under 0.005
hcp_alz_pval_0.005 = hcp_alz[hcp_alz$pval<0.005,]

# filter FDR under 0.005
hcp_alz_FDR_0.005 = hcp_alz[hcp_alz$FDR<0.005,]

# filter pval under 0.001
hcp_alz_pval_0.001 = hcp_alz[hcp_alz$pval<0.001,]

# filter FDR under 0.001
hcp_alz_FDR_0.001 = hcp_alz[hcp_alz$FDR<0.001,]

# write complete table
write.table(hcp_alz, file="hcp_alz_ctx_connectivity.tsv", sep="\t", quote = F, row.names = F)

# write pval 0.005 filtered table
write.table(hcp_alz_pval_0.005, file="hcp_alz_ctx_connectivity_pval_0.005.tsv", sep="\t", quote = F, row.names = F)

# write pval 0.001 filtered table
write.table(hcp_alz_pval_0.001, file="hcp_alz_ctx_connectivity_pval_0.001.tsv", sep="\t", quote = F, row.names = F)

# write FDR 0.005 filtered table
write.table(hcp_alz_FDR_0.005, file="hcp_alz_ctx_connectivity_FDR_0.005.tsv", sep="\t", quote = F, row.names = F)
