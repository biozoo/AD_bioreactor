**AD_bioreactor
**
This repository contains R scripts, data, and outputs supporting the main findings presented in Figures 1–5. The materials are organized into three sections:

_1. Code_
a) Scripts for Generating Figures:
Fig1_Demo_DEV_CH4_D47.R: Generates Figure 1.
Fig2_3_Demo_Heatmap_D47.R: Generates Figures 2 and 3 heatmaps.
Fig4_Demo_Warning_D47.R: Produces Figure 4.
Fig5_Demo_Cascade_D47.R: Generates Figure 5.
Additional metagenomics/metaproteomics analysis scripts (e.g., for Figure 2c) are included in the file:
Fig2C_Rscripts_data_for_metaomics.zip.

b) Analytical Scripts:
Demo_CCMSelect_OTU_metabolite_D47.R: Performs CCM causality analysis to select those OTUs and metabolites with signficant causal effects on CH4 production.
Demo_DEV_CRNOTU_D47.R: Conducts DEV early-warning signal analysis for OTUs, metabolites, and network properties.
Demo_MDR_metaboliteOTUs.R: Reconstructs the dynamic interactome between 234 microbial OTUs, 147 chemical fingerprints, and 7 VFAs using MDR S-map.
Demo_MDR_InfluenceOnCH4.R: Quantifies OTU and metabolite influences on CH4 production near the tipping point (Day 37).
Demo_NetworkProperty_computation.R: Computes network properties displayed in Figure 3 for each sampling day.
c) Supporting Functions:
Reusable functions for computations in (a) and (b):
Function_DEV_V2.R
Function_MDR.R
Function_MDR_inferM.R: MDR S-map for generating non-square interaction Jacobian matrices.
Function_network.R

_2. Output_
Derived datasets from the computations performed in the R scripts (Sections 1a and 1b).

_3. Data_
a) Supplementary dataset:
Data_S1.xlsx: Details of selected OTUs, metabolites, and network properties.
b) Reactor Datasets:
Reactor_data_envi.csv: Time series data for CH4 and relevant performance measurements, VFA, and environmental variables.
OLR_Perturbed_finger_D47.csv: Time series data for chemical fingerprint in the perturbed reactor.
OLR_Control_finger_D47.csv: Time series data for chemical fingerprint in the control reactor.
OLR_Perturbed_otu_D47_Species.csv: Time series data for microbial OTUs in the perturbed reactor.
OLR_Control_otu_D47_Species.csv: Time series data for microbial OTUs in the control reactor.
OLR_Perturbed_crnotu_verc_CH4_Species_D47.csv: Combined time series data for CH4, OTUs, and metabolites in the perturbed reactor.
OLR_Control_crnotu_verc_CH4_Species_D47.csv: Combined time series data for CH4, OTUs, and metabolites in the control reactor.
demo_otu_taxa_D47.csv: Taxonomic details of microbial OTUs
Demo_AD_small_summedLFQ_scaled.csv: Metaproteomics data involved in core biochemical pathways.
This structure ensures reproducibility of all main findings presented in the study.

