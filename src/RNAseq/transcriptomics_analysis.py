#%%

import scanpy as sc
import decoupler as dc

# Only needed for processing https://decoupler-py.readthedocs.io/en/latest/notebooks/bulk.html
import numpy as np
import pandas as pd
from anndata import AnnData
import omnipath
import pybiomart
import pydeseq2
# Import DESeq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import os

#%%


# Retrieve CollecTRI gene regulatory network
collectri = dc.get_collectri(organism='human', split_complexes=False)
# Read raw data and process it
file_path = '/Users/charliebarker/Desktop/Melanoma_Resistance/data/RNAseq/data/geneCounts_fixed.csv'
design_path = '/Users/charliebarker/Desktop/Melanoma_Resistance/data/RNAseq/Study_design.csv'

adata = pd.read_csv(file_path)
design = pd.read_csv(design_path)
column_mapping = dict(zip(adata.columns[1:], design.set_index('Study_ID')['New_Sample_name']))
adata.rename(columns=column_mapping, inplace=True)

#%%

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
adata = adata.T
# Remove columns with NaN column names
adata = adata.loc[:, ~adata.columns.isna()]

adata = AnnData(adata, dtype=np.float32)
adata.var_names_make_unique()
#%%

#Inside an AnnData object, there is the .obs attribute where we can store the metadata of our samples.

# Process treatment information
adata.obs['sample-name'] = adata.obs.index.tolist()
adata.obs['trametinib'] = ['treatment' if 'Trametinib' in sample_id else 'control' for sample_id in adata.obs.index]
adata.obs['vemurafenib'] = ['treatment' if 'Vermurafenib' in sample_id else 'control' for sample_id in adata.obs.index]
adata.obs['combination'] = ['treatment' if 'and' in sample_id else 'control' for sample_id in adata.obs.index]
adata.obs['ARID1A_KO'] = ['treatment' if 'ARID1A_KO' in sample_id else 'control' for sample_id in adata.obs.index]

#%%
#Filter genes by expression

# Visualize metadata

min_count = 5
min_total_count = 10


dc.plot_filter_by_expr(adata, group=None, min_count=min_count, min_total_count=min_total_count, large_n=1, min_prop=.4)

# Obtain genes that pass the thresholds
genes = dc.filter_by_expr(adata, group=None, min_count=min_count, min_total_count=min_total_count, large_n=1, min_prop=.4)

# Filter by these genes
adata = adata[:, genes].copy()


# %%
# Run DESEQ2

# List of contrasts focusing on untreated vs treated within each genetic background
contrasts = [
    ("Untreated WT", "Untreated ARID1A-KO"),
    ("Untreated WT", "Vermurafenib-1uM WT"),
    ("Untreated WT", "Trametinib-10nM WT"),
    ("Untreated WT", "vemurafenib-and-trametinib WT"),
    ("Untreated ARID1A-KO", "Vermurafenib-1uM ARID1A-KO"),
    ("Untreated ARID1A-KO", "Trametinib-10nM ARID1A-KO"),
    ("Untreated ARID1A-KO", "vemurafenib-and-trametinib ARID1A-KO")
]

# Output directory
output_dir = '/Users/charliebarker/Desktop/Melanoma_Resistance/results/transcriptomics'

# Function to run DESeq2 analysis
def run_deseq2_analysis(adata, contrast, output_dir):
    condition1, condition2 = contrast
    # Set the correct reference level for the contrast
    ref_level = ['sample-name', condition1]
    
    # Build DESeq2 object with the specific reference level for the contrast
    dds = DeseqDataSet(
        adata=adata,
        design_factors='sample-name',
        ref_level=ref_level,
        refit_cooks=True
    )
    # Compute LFCs
    dds.deseq2()
    
    # Create the DeseqStats object with the specific contrast
    stat_res = DeseqStats(dds, contrast=["sample-name", condition2, condition1])
    
    # Compute Wald test
    stat_res.summary()
    
    # Extract results
    results_df = stat_res.results_df
    
    # Create a filename based on the contrast
    contrast_name = f"{condition1.replace(' ', '_').replace('-', '_')}_vs_{condition2.replace(' ', '_').replace('-', '_')}"
    file_path = f"{output_dir}/{contrast_name}_lfc.csv"
    
    # Save results to CSV
    results_df.to_csv(file_path)
    
    print(f"Results saved to {file_path}")

# Initialize DESeq2 and run for each contrast
for condition1, condition2 in contrasts:
    # Call the analysis function
    run_deseq2_analysis(adata, (condition1, condition2), output_dir)


#%%
#transcription factor inference 

# Define the directory where DESeq2 results are saved
results_dir = '/Users/charliebarker/Desktop/Melanoma_Resistance/results/transcriptomics'

# Define the directory for transcription factor activity results
tf_activity_dir = '/Users/charliebarker/Desktop/Melanoma_Resistance/results/tf_activity'

# Define the CollecTRI network
collectri = dc.get_collectri(organism='human', split_complexes=False)

# Function to run transcription factor inference from DESeq2 results
def run_tf_inference(results_dir, tf_activity_dir, collectri):
    # Get all result files in the results directory
    result_files = [f for f in os.listdir(results_dir) if f.endswith('_lfc.csv')]
    
    # Create the TF activity directory if it does not exist
    if not os.path.exists(tf_activity_dir):
        os.makedirs(tf_activity_dir)
    
    # Process each result file
    for result_file in result_files:
        # Construct experiment name from result file name
        exp_name = result_file.replace('_lfc.csv', '')
        
        # Full path to the result file
        result_path = os.path.join(results_dir, result_file)
        
        # Read the DESeq2 results
        results_df = pd.read_csv(result_path, index_col=0)
        
        # Retrieve CollecTRI gene regulatory network
        mat = results_df[['stat']].T.rename(index={'stat': exp_name})
        
        # Run ULM for transcription factor activity and p-values
        tf_acts, tf_pvals = dc.run_ulm(mat=mat, net=collectri, verbose=True)
        
        # Save TF activities and p-values to CSV files
        tf_acts_file = f"{tf_activity_dir}/{exp_name}_tf_acts.csv"
        tf_pvals_file = f"{tf_activity_dir}/{exp_name}_tf_pval.csv"
        tf_acts.T.to_csv(tf_acts_file)
        tf_pvals.T.to_csv(tf_pvals_file)
        
        print(f"TF activities saved to {tf_acts_file}")
        print(f"TF p-values saved to {tf_pvals_file}")


# Run the transcription factor inference
run_tf_inference(results_dir, tf_activity_dir, collectri)

# %%
