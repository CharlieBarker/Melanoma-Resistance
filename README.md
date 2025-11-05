# ARID1A-Mediated Resistance in Melanoma: Multi-Omics Analysis

In this repository is the code to explore **ARID1A-driven resistance mechanisms** in melanoma treated with BRAF and MAPK inhibitors. This repository contains the code integrating transcriptomics, proteomics, phosphoproteomics, and kinomics to identify key signaling pathways and resistance mechanisms using **network-based analysis**.

## Overview

This project examines the signaling and immune dynamics underlying drug resistance in melanoma in the context of ARID1A loss. Using **graph-theoretical and matrix factorization approaches**, the study uncovers actionable resistance pathways.

### Key Findings:
- **ARID1A-KO cells** sustain MAPK1/3 and JNK signaling post-treatment, bypassing feedback mechanisms.
- Resistance is mediated by **elevated RTK activity** (e.g., EGFR, ROS1), disruption of PKC dynamics, and JUN activation.
- ARID1A KO elicits immune evasion through reduced HLA expression and enriched extracellular matrix proteins.
- ARID1A KO function feeds through critical signalling nodes: **PRKD1, JUN, and NCK1**.

## Repository Structure

scripts to reproduce figures

figure 1
1D - mofa_supplementary.r

figure 2
2A - ~/Desktop/Melanoma_Resistance/src/phuego/supplementary_plot_networks.R set factor variable to factor1
2B - ~/Desktop/Melanoma_Resistance/src/phuego/enrichr_ggplot.R
2C - ~/Desktop/Melanoma_Resistance/src/phuego/ego_nets.R set factor variable to factor1
2D - ~/Desktop/Melanoma_Resistance/src/visualisation/figure2_heatmap.R
2E - ~/Desktop/Melanoma_Resistance/src/visualisation/figure2_heatmap.R

figure 3
3A - ~/Desktop/Melanoma_Resistance/src/phuego/supplementary_plot_networks.R set factor variable to factor2
3B - ~/Desktop/Melanoma_Resistance/src/phuego/ego_nets.R set factor variable to factor2
3C - ~/Desktop/Melanoma_Resistance/src/visualisation/combination_phosphosite.R
3D - ~/Desktop/Melanoma_Resistance/src/phuego/enrichr_ggplot.R

figure 4
4A - ~/Desktop/Melanoma_Resistance/src/ARID1A_heatdiffusion/ARID1A_heatmaps.R
4B - ~/Desktop/Melanoma_Resistance/src/ARID1A_heatdiffusion/ARID1A_heatmaps.R
4C - ~/Desktop/Melanoma_Resistance/src/ARID1A_heatdiffusion/random_walk_diffusion.R
4E - ~/Desktop/Melanoma_Resistance/src/ARID1A_heatdiffusion/ephrin_diffusion.R
4G - ~/Desktop/Melanoma_Resistance/src/kinomic_analysis/prkd1_story.R

figure 5
5A -
5B -
5C -
5D - 
5E -

