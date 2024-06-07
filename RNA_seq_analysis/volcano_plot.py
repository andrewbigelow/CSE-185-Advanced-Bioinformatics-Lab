'''
Working files and directories:
    *.genes.results: gene-level expression results output by RSEM for each condition/replicate.
    GRCm38.75.gene_names contains a mapping between ENSEMBL gene ids (ENSMUSG...) and gene names (e.g. Cdc45, Klf6)

After running an R script to perform Differential Analysis Expression using DESeq2 on the given RNA-seq data 
and outputing the results to "chow_vs_hfd_deseq2.csv", I used the ouput data to create a volcano plot highlighting
genes that are significantly differentially expressed in red and annotating the names of the top differentially 
expressed genes.

'''
import pandas as p
import matplotlib.pyplot as plt
import numpy as np

# Store CSV data
csv_data_untrimmed = p.read_csv('RNA_seq_analysis/chow_vs_hfd_deseq2.csv')

# Map gene names
gene_mapping = p.read_csv('RNA_seq_analysis/GRCm38.75.gene_names', sep='\t', header=None, index_col=0)
gene_mapping_dict = gene_mapping[1].to_dict()

# Store data as [log2FoldChange, -log10 pvalue]
csv_data = csv_data_untrimmed.iloc[:, [2, 5]]
csv_data.iloc[:, 1] = - np.log10(csv_data.iloc[:, 1])

# Calculate how many genes are significant
num_sig = 0
for index, row in csv_data.iterrows():
    if row[1] > 5:
        num_sig += 1

print('Total number of significant genes ', num_sig, '\n')

# Plot the data and highlight genes with -log10 pvalue > 5
plt.scatter(csv_data.iloc[:, 0], csv_data.iloc[:, 1], c=np.where(csv_data.iloc[:, 1] < 5, 'black', 'red'), alpha = 1, s = 1)


# Edit axis
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10 p value')

# Store 10 largest genes
top_genes = csv_data.nlargest(10, csv_data.columns[1])  
for i, row in top_genes.iterrows():
    # Annotate top genes with their complements from the other file
    gene_name = csv_data_untrimmed.iloc[i, 0]
    complement = gene_mapping_dict.get(gene_name, "Unknown")
    plt.annotate(complement, (row[csv_data.columns[0]], row[csv_data.columns[1]]), fontsize=8)
    print(f"Gene name: {complement}, log2 Fold Change: {row[csv_data.columns[0]]}, p-value: {10 ** -(row[csv_data.columns[1]])}")

# Save plot in directory
plt.savefig("volcano_plot.png")