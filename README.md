**#AD_bioreactor
**
This repository contains R scripts, data, and outputs supporting the main findings presented in Figures 1–5. The materials are organized into three sections:

_1. Code_
a) Scripts for Generating Figures:
Fig1_Demo_DEV_CH4.R: Generates Figure 1.
Fig2_3_Demo_Heatmap.R: Generates Figures 2 and 3 heatmaps.
Fig4_Demo_Warning.R: Produces Figure 4.
Fig5_Demo_Cascade.R: Generates Figure 5.
Additional metagenomics/metaproteomics analysis scripts (e.g., for Figure 2c) are included in the file:
Fig2C_Rscripts_data_for_metaomics.zip.

b) Analytical Scripts:
Demo_CCM_Causality_Test.R: Performs CCM causality analysis to calculate the causal strength of OTUs and metabolites on CH4 production.
Demo_DEV_CRNOTU_ls_s.R: Conducts DEV early-warning signal analysis for OTUs, metabolites, and network properties.
Demo_MDR_metaboliteOTUs.R: Reconstructs the dynamic interactome between 250 OTUs and 131 metabolites using MDR S-map.
Demo_MDR_InfluenceOnCH4.R: Quantifies OTU and metabolite influences on CH4 production near the tipping point (Day 37).
Demo_NetworkProperty_computation.R: Computes network properties displayed in Figure 3 for each sampling day.
c) Supporting Functions:
Reusable functions for computations in (a) and (b):

Function_DEV.R
Function_MDR.R
Function_MDR_inferM.R: MDR S-map for generating non-square interaction Jacobian matrices.
Function_network.R

_2. Output_
Derived datasets from the computations performed in the R scripts (Sections 1a and 1b).

_3. Data_
a) Supplementary Tables:
Extended Data Tables 3–5 (Extended Data Table 3-5.xlsx): Details of selected OTUs, metabolites, and network properties.
b) Reactor Datasets:
reactor_ch4_ews.csv: Time series data for CH4 production and OTU richness.
OLR_crn_verc_ch4p_Species_tau132.csv: Time series data for metabolites.
OLR_otu_verc_ch4p_Species_tau132.csv: Time series data for microbial OTUs.
Perturbed_otucrn_ch4_verc_Species.csv: Combined time series data for CH4, OTUs, and metabolites in the perturbed reactor.
Control_otucrn_ch4_verc_Species.csv: Combined time series data for CH4, OTUs, and metabolites in the control reactor.
demo_otu_taxa.csv: Detailed taxonomic information for microbial OTUs.
This structure ensures reproducibility of all main findings presented in the study.


    
    
    


