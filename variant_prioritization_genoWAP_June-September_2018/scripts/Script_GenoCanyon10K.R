# load required package "bigmemory"

library(bigmemory, lib.loc="/neurospin/brainomics/2019_Riquelme/genoWAP_variant_prioritization_september_2018/GenoSuite/GenoCanyon/lib_big_memory/")
    
GenoCanyon10K_dir <- "/neurospin/brainomics/2019_Riquelme/genoWAP_variant_prioritization_september_2018/GenoSuite/GenoCanyon/GenoCanyon_10K"

setwd(GenoCanyon10K_dir)

source("/neurospin/brainomics/2019_Riquelme/genoWAP_variant_prioritization_september_2018/GenoSuite/GenoWAP-V1.2/Script_Appendix.R")

# The user needs to define two vectors: Chr and Pos
# Chr is the chromosome information for all SNPs. Use 23 and 24 to denote chromosome X & Y, respectively. CharToNum function in the appendix can help convert "chr#" -> #.
# Pos is hg19 position for all SNPs.
# Examples (3 SNPs in total): Chr=c(1,2,3), Pos=c(1001, 40000, 12345)

args <- (commandArgs(TRUE))
input_file <- args[[1]]
output_file <- args[[2]]

print("Loading GWAS input file")
gwas_table <- read.table(file = input_file)

print("Reading SNP's chromosome number and genomic position")
Chr <- c(gwas_table[,1])
Pos <- c(gwas_table[,2])

# Given Chr and Pos, the following script can extract GenoCanyon10K scores for all SNPs
GenoCanyon = rep(NA, length(Chr))

pb = txtProgressBar(0, 45, style=3)

print("Extracting GenoCanyon10K general functional score for each SNP's position")

for(i in 1:45){
    Region = as.logical((Chr == Chr.GenoCanyon10K[i])*(Pos %in% (PosStart.GenoCanyon10K[i]:PosStop.GenoCanyon10K[i])))
    if(sum(Region) > 0){
        mat = attach.big.matrix(files.GenoCanyon10K[i])
        GenoCanyon[Region] = mat[Pos[Region] - PosStart.GenoCanyon10K[i] + 1]
    }
    setTxtProgressBar(pb, i)
}

# Add Chr and Pos information to GenoCanyon extracted scores to match the correct data format for GenoWAP
GenoCanyon <- cbind(gwas_table[,1:2],GenoCanyon)

# Output GenoCanyon10K scores
print(paste("Writing scores in file: ", output_file, sep = ""))
write.table(GenoCanyon, output_file, quote=F, row.names=F, col.names=F, sep = "\t")




