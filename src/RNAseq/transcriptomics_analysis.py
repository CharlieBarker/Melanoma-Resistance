import scanpy as sc
import decoupler as dc

# Only needed for processing
import numpy as np
import pandas as pd
from anndata import AnnData
import os

# Set the desired directory path
directory_path = '/Users/charliebarker/Desktop/Melanoma_Resistance'

# Change the current working directory to the desired directory
os.chdir(directory_path)

# Read raw data and process it
file_path = './data/RNAseq/data/geneCounts_fixed.csv'
design_path = './data/RNAseq/Study_design.csv'

adata = pd.read_csv(file_path)
design = pd.read_csv(design_path)
column_mapping = dict(zip(adata.columns[1:], design.set_index('Study_ID')['New_Sample_name']))
adata.rename(columns=column_mapping, inplace=True)

# Retrieve gene symbols
annot = sc.queries.biomart_annotations("hsapiens",
        ["ensembl_gene_id", "external_gene_name"],
        use_cache=False
    ).set_index("ensembl_gene_id")

# Filter genes not in annotation
adata.set_index('ENSEMBL_ID', inplace=True)
adata = adata[adata.index.isin(annot.index)]
# Assign gene symbols
adata['gene_symbol'] = [annot.loc[ensembl_id,'external_gene_name'] for ensembl_id in adata.index]
adata = adata.reset_index().rename(columns={'index': 'ensembl_gene_id'}).set_index('gene_symbol')

# Remove rows with all zero values (empty rows)
adata = adata[(adata != 0).any(axis=1)]

# Rename columns with underscores to replace with spaces
adata.columns = adata.columns.str.replace('__', ' ')

# Remove the 'ENSEMBL_ID' column
adata = adata.drop(columns='ENSEMBL_ID')

# Transform to AnnData object
adata = adata.T.head()
# Remove columns with NaN column names
adata = adata.dropna(axis=1, how='all')

print(adata.head())
adata = AnnData(adata, dtype=np.float32)
adata.var_names_make_unique()
adata