setwd("/home/tr256641/Téléchargements/JOBIM_abstract/prs_results_hcp/cytoscape_data/")

# load tsv 
hcp_asd = read.csv("hcp_asd_best_cytoscape.csv", sep="\t")

# select columns of interest
hcp_asd = hcp_asd[,c("Phenotype","P","X.log.P.")]
names(hcp_asd) = c("ctx_connexion", "pval", "-log10(pval)")

# split column of connexions between 2 ctx regions into 2 columns of ctx_regions and save as a new df
ctx_regions <- data.frame(do.call('rbind', strsplit(as.character(hcp_asd$ctx_connexion),'_',fixed=TRUE)))
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
hcp_asd = hcp_asd[,c("pval", "-log10(pval)")]
hcp_asd = cbind(ctx_regions, hcp_asd)

# create orientation Left-right column
orientation1 <- sapply(strsplit(hcp_asd$ctx_region_1,"-"), `[`, 1)
orientation2 <- sapply(strsplit(hcp_asd$ctx_region_2,"-"), `[`, 1)
hcp_asd = cbind(orientation1, orientation2, hcp_asd)

# sort by pval
hcp_asd = hcp_asd[order(hcp_asd$pval),]

# fdr corrected pval in a new column
hcp_asd_FDR = p.adjust(p=hcp_asd$pval, method = "fdr")

# bind this new column to the first table
hcp_asd = cbind(hcp_asd, FDR = hcp_asd_FDR)

# filter FDR under 0.05
hcp_asd_FDR_0.05 = hcp_asd[hcp_asd$FDR<0.05,]

# filter pval under 0.005
hcp_asd_pval_0.005 = hcp_asd[hcp_asd$pval<0.005,]

# filter FDR under 0.005
hcp_asd_FDR_0.005 = hcp_asd[hcp_asd$FDR<0.005,]

# filter pval under 0.001
hcp_asd_pval_0.001 = hcp_asd[hcp_asd$pval<0.001,]

# filter FDR under 0.001
hcp_asd_FDR_0.001 = hcp_asd[hcp_asd$FDR<0.001,]

# write complete table
write.table(hcp_asd, file="hcp_asd_ctx_connectivity.tsv", sep="\t", quote = F, row.names = F)

# write pval 0.005 filtered table
write.table(hcp_asd_pval_0.005, file="hcp_asd_ctx_connectivity_pval_0.005.tsv", sep="\t", quote = F, row.names = F)

# write pval 0.001 filtered table
write.table(hcp_asd_pval_0.001, file="hcp_asd_ctx_connectivity_pval_0.001.tsv", sep="\t", quote = F, row.names = F)

# write FDR 0.005 filtered table
write.table(hcp_asd_FDR_0.005, file="hcp_asd_ctx_connectivity_FDR_0.005.tsv", sep="\t", quote = F, row.names = F)
