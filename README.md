# AD_bioreactor


R script for producing the main findings presented in Fig. 1-5.
There are three documents include all the materials needed for showing the findings in main figures, including i) Code, ii) Data, iii) Output.
i) Code:
  a)For producing the main figures, we need to execute, 'Fig1_Demo_DEV_CH4.R', 'Fig2_3_Demo_Heatmap.R', 'Fig4_Demo_Warning.R', and 'Fig5_Demo_Cascade.R'.
  In addition, the R script for metagenomics/metaproteomics analysis (e.g., Fig. 2c) were included separately in the zip file, 'Fig2C_Rscripts_data_for_metaomics.zip'.
  b) The script for otaining the analytical results for producing figures were enclosed, including: 
  'Demo_CCM_Causality_Test.R': Performing CCM cuasal analysis to obtained the causal strength of OTUs and metabolites to CH4 production;    
  'Demo_DEV_CRNOTU_ls_s.R': Performing DEV early warning signal analysis for OTUs, metabolites, and network properties.
  'Demo_MDR_metaboliteOTUs.R': Reconstructing dynamic interactome that included the mutual interactions between 250 OTUs and 131 metabolites using MDR S-map.
  'Demo_MDR_InfluenceOnCH4.R':Quantifying the influences of OTUs and metabolites on CH4 production around the tipping point (Day 37)
  'Demo_NetworkProperty_computation.R': Calculating network properties shown in the figure 3 at each sampling day.
  c) The functions for supporting the computations in a) and b), including 'Function_DEV.R', 'Function_MDR.R', 'Function_MDR_inferM.R' (MDR S-map for producing non-square     
  interaciton Jacobian matrix), and 'Function_network.R'
ii) Output: The derived dataset from the computation of R-script espcially in code i-a and i-b.
iii) Data:
  a) Supplementary tables:  showing the details of the selected OTUs, metabolites, and network properties, Extended Data Table 3-5.xlsx
  b) Reactor dataset: 
    'reactor_ch4_ews.csv': time series data for CH4 production and OTU richness.
    'OLR_crn_verc_ch4p_Species_tau132.csv':  Time series data for metabolites.
    'OLR_otu_verc_ch4p_Species_tau132.csv':  Time series data for microbial OTUs
    'Perturbed_otucrn_ch4_verc_Species.csv': Compiled time series data for CH4, microbial OTUs, & metabolites in the perturbed reactor.
    'Control_otucrn_ch4_verc_Species.csv': Compiled time series data for CH4, microbial OTUs, & metabolites in the control reactor.
    'demo_otu_taxa.csv': detailed taxonomic information for microbial OTUS.
    
    
    


