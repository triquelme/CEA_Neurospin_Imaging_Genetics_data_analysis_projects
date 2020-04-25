## Geno-tool suite explanations

Author : Thomas Riquelme
Date : September 5th 2018

### Geno-tool suite

The Geno-tool suite:
- GenoCanyon: compute general(=not tissue specific)-functionality scores for each locus
- GenoSkyline: compute tissue-specific functionality scores for each locus
- GenoWAP: GWAS signal prioritization

### GenoCanyon

GenoCanyon general scores for each locus of the hg19 human genome can be downloaded here:
http://genocanyon.med.yale.edu/GenoCanyon_Downloads.html
about 25Gb, binary files 
Script available to extract GenoCanyon10K scores.

### GenoSkyline

GenoSkyline tissue specific and GenoSkyline-plus cell specific scores can be downloaded here:
http://genocanyon.med.yale.edu/GenoSkyline

Genoskyline: 9 tissues bedGraph files, 300-500Mb each

Genoskyline-plus: 127 cells/tissues bedGraph files, 100Mb each (11Gb in total)

### GenoWAP

For SNP prioritization with GenoWAP, it is necessary to extract the functionality scores from GenoCanyon and GenoSkyline only for GWAS identified-SNPs loci.

To do so use:
- Script_GenoCanyon10K.R to extract general functionality scores from GenoCanyon10K binary files
(available here: http://genocanyon.med.yale.edu/GenoCanyon_Downloads.html)
- extractScores.py script to extract tissue specific scores (available with genoWAP)

Then you can call genoWAP on your GWAS data using the two custome files containing the extracted scores for GWAS-SNPs positions:

python GenoWAP [-a EXTRACTED_GENERAL_FUNCTIONAL_SCORES_DATA_PATH] [-ts EXTRACTED_TISSUE_FUNCTIONAL_SCORES_DATA_PATH] GWAS_DATA_PATH

