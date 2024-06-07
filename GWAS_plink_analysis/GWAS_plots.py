'''
Files used: 
lab3_gwas.phen: Normalized LDL values for each sample
lab3_gwas.vcf.gz: A VCF file containing the LD-pruned SNPs for a subset of the 1000 Genomes samples
lab3_gwas.vcf.gz.tbi: The VCF file index

Commands ran:
* Run plink on data to perform GWAS:

plink --vcf ~/public/lab3/lab3_gwas.vcf.gz --linear \
--pheno ~/public/lab3/lab3_gwas.phen --allow-no-sex --maf 0.05 --out lab3_gwas
NOTE: After plotting results, there seemed to be a confounding factor affecting the GWAS


* Run plink's PCA function on data to highlight main principle components affecting data:

plink --vcf ~/public/lab3/lab3_gwas.vcf.gz --pca 3 --out lab3_gwas


* Rerun plink's GWAS while accounting for the top 3 PCs from prior command stored in lab3_gwas.eigenvec:

plink --vcf ~/public/lab3/lab3_gwas.vcf.gz --pheno ~/public/lab3/lab3_gwas.phen \
--covar lab3_gwas.eigenvec --linear hide-covar --allow-no-sex --maf 0.05 --out lab3_gwas_covar


* Run plink's clump tool and bcftools norm tool to take the adjusted data and clump it to find unique genetic varation

bcftools norm --rm-dup all ~/public/lab3/lab3_gwas.vcf.gz > lab3_gwas_refined.vcf.gz
plink --vcf lab3_gwas_refined.vcf.gz --clump lab3_gwas_covar.assoc.linear
'''


"""code for Manhattan and QQ plot after correcting for non linear relationship"""
from qqman import qqman
import matplotlib.pyplot as plt
import pandas as pd
import zipfile

# Extract zip file
zip_filename = "GWAS_plink_analysis/lab3_gwas_covar.assoc.linear.zip"
with zipfile.ZipFile(zip_filename, "r") as zip_ref:
    zip_ref.extractall()

# Read the extracted CSV file
with zipfile.ZipFile(zip_filename, "r") as zip_ref:
    with zip_ref.open("lab3_gwas_covar.assoc.linear", "r") as csv_file:
        covdata = pd.read_csv(csv_file, delim_whitespace=True)

# With extracted data built plots
fig, (ax0, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
fig.set_size_inches((15, 5))
qqman.manhattan(covdata, ax=ax0)
qqman.qqplot(covdata, ax=ax1)

# Save plots in directory
fig.savefig("manhattan_and_qq_plot.png")