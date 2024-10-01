# regeQTL - eQTL analysis of Regulatory trios
regeQTL performs regulatory eQTL (expression Quantitative Trait Locus) analysis on TG-TF-SNV (Target Gene-Transcription Factor-Single Nucleotide Variant) trios, to identify potential regulatory relationships.

# Install
devtools::install_github('liliulab/reg-eQTL')

# Usage
library(regeQTL)<br />
process.regeqtl(expr.data, cov.data, trio.data, gt.data, out.dir)<br /> 
process.seqtl(expr.data, cov.data, pair.data, gt.data, out.dir)<br />


# Input Data Requirements
Expression Data: A matrix where rows are genes and columns are samples.<br />
Covariate Data: A matrix of covariates with a "sample_id" column to link with the expression data.<br />
TG-TF-SNV Trio Data: Data containing columns for Target Gene (gene), Transcription Factor (TF), and single nucleotide variant (SNP).<br />
Genotype Data: A matrix where rows are samples and columns are SNVs, and entries represent genotypes.

# Output
eQTL results will be saved to directory mentioned in out.dir.


# Example Input
A sample input is provided in the Examples folder. 

# Reference
Manscript under review.

# Contributors
This package was developed by Li Liu and Rekha Mudappathi. Please contact Li Liu at [liliu@asu.edu] for any questions or suggestions.

