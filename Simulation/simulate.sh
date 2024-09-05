#!/bin/bash

# GT.sim is a text file with first column as number of SNVs, MAF category, lower allele frequency range, upper allele frequency range, odds ratio for disease (heterozygote), odds ratio for disease (homzygote) 
# example GT.sim file
# 200  A     0.01 0.05  1.00 1.00
# 200  B     0.05 0.10  1.00 1.00
# 200  C     0.10 0.20  1.00 1.00
# 200  D     0.20 0.50  1.00 1.00

# PLINK command to simulate genotype data for different MAF categories 
plink --simulate GT.sim --make-bed --out simulated_genotypes --simulate-ncontrols 670 --simulate-ncases 0

