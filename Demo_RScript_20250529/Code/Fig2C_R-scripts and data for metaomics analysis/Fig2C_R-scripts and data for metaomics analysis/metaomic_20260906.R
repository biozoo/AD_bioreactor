# =============================================================================
# PRODIGIO Functional Metaomics Analysis
# Author: Magnus Arntzen
# =============================================================================

# =============================================================================
# 1. Load libraries
# =============================================================================
options(stringsAsFactors = FALSE)

library(ampvis2)
library(ggplot2)
library(tibble)
library(tidyverse)
library(foreach)
library(gridExtra)
library(ggedit)
library(data.table)
library(patchwork)
library(data.table)
library(MetQy)      # v1.1.0 from Github
library(pheatmap)
library(ComplexHeatmap)
library(readxl)     # For reading Excel files
library(cowplot)
library(grid)
library(RColorBrewer)
library(circlize)
wd <- getwd()
setwd(dirname(file.choose()))
# =============================================================================
# 2. Load and preprocess data
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1 16S data
# -----------------------------------------------------------------------------
d <- amp_load(
  otutable = 'D2.2_Task2.1_45dayShock_otutable.txt',
  metadata = 'D2.2_Task2.1_45dayShock_metadata.xlsx',
)
d_rarefied <- amp_rarefy(d, 16693)
d_rarefied$metadata$samplingPoint[23] <- "R3"  # Error fix in annotation

# -----------------------------------------------------------------------------
# 2.2 Reactor chemical data
# -----------------------------------------------------------------------------
chem_control <- read_delim("D4_Table.EDM.Long.Shocks_NMBU_R1R2.csv", 
                           delim = ";", progress = FALSE, 
                           na = c("", "NaN"), comment = "#", 
                           locale = locale(encoding = "latin1", decimal_mark = ",")) %>%
  select(-Samples, -Days, -starts_with("OTU")) %>%
  mutate(Date = as.Date(`Data of sampling`, format = "%d.%m.%Y"))

chem_OLR <- read_delim("D4_Table.EDM.Long.Shocks_NMBU_R3R4.csv", 
                       delim = ";", progress = FALSE, 
                       na = c("", "NaN"), comment = "#", 
                       locale = locale(encoding = "latin1", decimal_mark = ",")) %>%
  select(-Samples, -Days, -starts_with("OTU")) %>%
  mutate(Date = as.Date(`Data of sampling`, format = "%d.%m.%Y"))

Chem <- rbindlist(list(chem_control, chem_OLR), use.names = TRUE, fill = TRUE) %>%
  select(-`Data of sampling`) %>%
  rename(Reactor = Group, Temp = 3) %>%
  relocate(Date)

Chem_long <- Chem %>%
  pivot_longer(cols = c(4:26))

# -----------------------------------------------------------------------------
# 2.3 Metagenomics data
# -----------------------------------------------------------------------------
MAG_annotations <- read_delim("D4_R3_OLR_310MAGs_withONT_annotated.txt", 
                              delim = "\t", progress = FALSE, 
                              na = c("", "NaN"), comment = "#")
MAG_annotations_df <- as.data.frame(MAG_annotations) %>%
  mutate(Bin = if_else(str_detect(Accn, "_sub"), 
                       paste0(str_split_i(Accn, "_", 1), "_sub"), 
                       str_split_i(Accn, "_", 1))) %>%
  mutate(Contig_id = str_remove(Accn, paste0(Bin, "_"))) %>%
  mutate(Protein_id = str_split_i(Contig_id, "_", -1)) %>%
  mutate(Contig_id = str_remove(Contig_id, paste0("_", Protein_id, "$"))) %>%
  select(Accn, Bin, Contig_id, Protein_id, CAZy, KoFamScan, Pfam, Pfam_desc, 
         InterPro, InterPro_desc)

# -----------------------------------------------------------------------------
# 2.4 Metaproteomics data
# -----------------------------------------------------------------------------
metaP <- read_delim("D4_18411proteins.txt", delim = "\t", 
                    progress = FALSE, na = c("", "NaN"), comment = "#")
metaP_df <- as.data.frame(metaP) %>%
  filter(is.na(`Indistinguishable Proteins`)) %>%
  select(-`Combined Total Peptides`, -`Protein ID`, -`Indistinguishable Proteins`) %>%
  left_join(MAG_annotations_df %>% select(Accn, Bin), by = c("Protein" = "Accn")) %>%
  mutate(across(1:12, ~ 2^.x))

sampling_dates <- c(`1:2022-02-17` = "1",
                    `2:2022-03-11` = "2",
                    `3:2022-03-24` = "3",
                    `4:2022-04-01` = "4",
                    `5:2022-04-09` = "5",
                    `6:2022-04-17` = "6",
                    `7:2022-04-25` = "7",
                    `8:2022-05-03` = "8",
                    `9:2022-05-13` = "9",
                    `10:2022-05-26` = "10",
                    `11:2022-06-17` = "11",
                    `12:2022-07-05` = "12")

metaP_df_for_export <- as.data.frame(metaP) %>%
  left_join(MAG_annotations_df %>% select(-Protein_id), by = c("Protein" = "Accn")) %>%
  rename(all_of(sampling_dates)) %>%
  mutate(across(1:12, ~ 2^.x))

# -----------------------------------------------------------------------------
# 2.5 Taxonomy data for MAGs
# -----------------------------------------------------------------------------
tax <- read_delim("D4_MAG_taxonomy.txt", delim = "\t", 
                  progress = FALSE, na = c("", "NaN"), comment = "#")
tax_df <- as.data.frame(tax)

# -----------------------------------------------------------------------------
# 2.6 Peptidase data for MAGs
# -----------------------------------------------------------------------------
MEROPS <- read_delim("D4_MEROPS_Interpro.txt", delim = "\t", 
                     progress = FALSE, na = c("", "NaN"), comment = "#")
MEROPS_df <- as.data.frame(MEROPS) %>%
  group_by(Interpro) %>%
  summarise(MEROPS = paste0(unique(MEROPS), collapse = ";"))

# =============================================================================
# 3. Stacked barplot for 16S community structure
# Basis for supplementary figure
# =============================================================================
d_Ctr_OLR_shock <- amp_subset_samples(d_rarefied, 
                                      samplingPoint %in% c("R1", "R2", "R3", "R4") &
                                        Date > as.Date("2022-03-21") &                
                                        Date < as.Date("2022-05-06"))

g_colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
              '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080',
              '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', 
              '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

ctr_olr_shock_genus <- amp_timeseries(d_Ctr_OLR_shock,
                                      time_variable = "Date",
                                      tax_aggregate = "Genus",
                                      tax_show = 20,
                                      group_by = "samplingPoint",
                                      tax_add = "Kingdom")

ctr_olr_shock_genus_data <- ctr_olr_shock_genus$data %>%
  select(Date, Genus, samplingPoint, Value) %>%
  mutate(Day = as.numeric(Date - as.Date("2022-03-21"))) %>%
  mutate(Condition = case_when(samplingPoint %in% c("R1", "R2") ~ "Control",
                               samplingPoint %in% c("R3", "R4") ~ "OLR",
                               .default = "Undetermined")) %>%
  left_join(Chem %>% select(Date, `CH4 (%)`, Reactor) %>% distinct(), 
            by = c("Date" = "Date", "samplingPoint" = "Reactor")) %>%
  select(-Date, -samplingPoint)

p_16S <- ggplot(ctr_olr_shock_genus_data, mapping = aes(x = Day, y = Value * 100)) +
  facet_wrap(vars(Condition), scales = "free_y", ncol = 2) + 
  geom_bar(aes(fill = Genus), stat = "identity", position = 'fill', 
           color = "black", size = 0.2) +
  scale_fill_manual(values = g_colors) +
  scale_y_continuous(expand = expansion(add = c(0, 0)), 
                     labels = scales::percent_format()) +
  labs(x = "Day", y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", linewidth = 0.5),
        strip.text = element_text(size = 14, face = "bold")) +
  geom_line(aes(x = Day, y = `CH4 (%)` / 100, group = 1), 
            linewidth = 2, color = "white", linetype = "solid")

ggsave("16S_community_structure.pdf", p_16S, width = 18, height = 7)

# =============================================================================
# 4. KEGG annotation and calculation of module completion fractions
# Purpose: Generate a lookup-object of links between KEGG K-numbers,
#          KEGG modules, and broader KEGG categories
# =============================================================================

# Using only top MAGs to simplify the computation
topMAGs <- c("S3.MetaBAT2000.30", "S7.MetaBAT2000.1", "COR3.MetaBAT2000.276", 
             "S8.MetaBAT2000.28", "S6.MetaBAT2000.28",
             "COR3.MetaBAT2000.103", "COR3.MetaBAT2000.65", "PROD-ONT2-BATall.6", 
             "S6.MaxBin2000.004", "PROD-ONT2-BATall.39",
             "COR3.MaxBin2000.062", "S3.MaxBin2000.052", "PROD-ONT2-BATall.59", 
             "PROD-ONT2-BATall.4", "S7.MetaBAT2000.6",
             "S6.MetaBAT2000.104", "S6.MetaBAT2000.121", "PROD-ONT2-MAXall.017", 
             "S6.MetaBAT2000.129", "S6.MetaBAT2000.90")

# Get all KEGG K-numbers per MAG
bin_vs_KO_allMAGs <- MAG_annotations_df %>%
  select(Bin, KoFamScan) %>%
  filter(!is.na(KoFamScan)) %>%
  separate_rows(KoFamScan, sep = ";") %>%
  group_by(Bin) %>%
  summarise(K = paste0(KoFamScan, collapse = ";"))

# Focus on top-MAGs
bin_vs_KO_selection <- bin_vs_KO_allMAGs %>%
  filter(Bin %in% topMAGs)

# MetQy analysis (slow, precomputed data provided)
# all_metqy <- query_genomes_to_modules(as.data.frame(bin_vs_KO_selection), 
#                                       GENOME_ID_COL = "Bin", 
#                                       GENES_COL = "K", 
#                                       META_OUT = TRUE, 
#                                       ADD_OUT = TRUE)
# saveRDS(all_metqy, file = "../2024 ONT/R_metqy.rds")

all_metqy <- readRDS(file = "R_metqy.rds")

# Create look-up object for later
K_M_lookup <- all_metqy$METADATA %>%
  select(MODULE_ID, DEFINITION) %>%
  separate_rows(DEFINITION, sep = "[\\|&() ]") %>%
  filter(DEFINITION != "") %>%
  group_by(DEFINITION) %>%
  summarise(Module = paste0(unique(MODULE_ID), collapse = ";"))

# =============================================================================
# 5. Metaproteomic abundance of curated AD functions
#    Irrespective of MAG origin
# =============================================================================

# Metaproteomics data in the relevant dates with functional annotation
MP <- metaP_df_for_export %>%
  select(2:10, CAZy, KoFamScan, InterPro, InterPro_desc) %>%
  pivot_longer(cols = 1:9, names_to = "Date", values_to = "LFQ") %>%
  mutate(Date = as.numeric(str_split_i(Date, ":", 1)))

# Main KEGG-metaproteomics object for later use
MP_kegg <- MP %>%
  select(KoFamScan, Date, LFQ) %>%
  filter(!is.na(KoFamScan)) %>%
  separate_rows(KoFamScan, sep = ";") %>%
  left_join(K_M_lookup, by = c("KoFamScan" = "DEFINITION")) %>%
  filter(!is.na(Module)) %>%
  select(-KoFamScan) %>%
  separate_rows(Module, sep = ";") %>%
  group_by(Module, Date) %>%
  summarise(summedLFQ = log2(sum(LFQ, na.rm = TRUE)),
            nExpressed = sum(!is.na(LFQ))) %>%
  mutate(summedLFQ = na_if(summedLFQ, -Inf))

# -----------------------------------------------------------------------------
# 5.1 Hydrolysis: Glycolysis, CAZymes, Lipases, Peptidases
# -----------------------------------------------------------------------------

# Glycolysis
MP_glycolysis <- MP_kegg %>%
  filter(Module == "M00002") %>%
  mutate(subclass = "Glycolysis, core module") %>%
  mutate(Class = "Glycolysis")

# CAZymes (selected)
MP_cazy <- MP %>%
  select(-KoFamScan, -InterPro, -InterPro_desc) %>%
  filter(!is.na(CAZy)) %>%
  mutate(subclass = case_when(str_detect(CAZy, "GH") ~ "Glycoside hydrolases",
                              str_detect(CAZy, "GT") ~ "Glycosyl transferases",
                              str_detect(CAZy, "PL") ~ "Polysaccharide lyases",
                              str_detect(CAZy, "AA") ~ "Auxiliary activities",
                              str_detect(CAZy, "CE") ~ "Carbohydrate esterases")) %>%
  filter(!is.na(subclass)) %>%
  group_by(subclass, Date) %>%
  summarise(summedLFQ = log2(sum(LFQ, na.rm = TRUE)),
            nExpressed = sum(!is.na(LFQ))) %>%
  mutate(summedLFQ = na_if(summedLFQ, -Inf)) %>%
  mutate(Class = "CAZyme")

# Lipases
search_lipases <- c("triacylglycerol lipase", "lipase", "butyrinase", "tributyrinase", 
                    "Tween hydrolase", "steapsin", "triacetinase", "tributyrin esterase", 
                    "Tweenase", "amno N-AP", "Takedo 1969-4-9", "Meito MY 30", 
                    "Tweenesterase", "GA 56", "capalase L", "triglyceride hydrolase", 
                    "triolein hydrolase", "tween-hydrolyzing esterase", "amano CE", 
                    "cacordase", "triglyceridase", "triacylglycerol ester hydrolase",
                    "amano P", "amano AP", "PPL", "glycerol-ester hydrolase", "GEH", 
                    "meito Sangyo OF lipase", "hepatic lipase", "lipazin", 
                    "post-heparin plasma protamine-resistant lipase",
                    "salt-resistant post-heparin lipase", "heparin releasable hepatic lipase", 
                    "amano CES", "amano B", "tributyrase", "triglyceride lipase", 
                    "liver lipase", "hepatic monoacylglycerol acyltransferase", 
                    "PNLIP", "LIPF")
search_lipases_re <- paste0("(?:", paste0(search_lipases, collapse = "|"), ")")

MP_lipases <- MP %>%
  select(-KoFamScan, -CAZy) %>%
  filter(!is.na(InterPro)) %>%
  filter(str_detect(InterPro_desc, search_lipases_re)) %>%
  mutate(Class = "Lipase") %>%
  group_by(Date, Class) %>%
  summarise(summedLFQ = log2(sum(LFQ, na.rm = TRUE)),
            nExpressed = sum(!is.na(LFQ))) %>%
  mutate(summedLFQ = na_if(summedLFQ, -Inf)) %>%
  mutate(subclass = "Lipase")

# Peptidases
MP_peptidases <- MP %>%
  select(-KoFamScan, -CAZy) %>%
  filter(!is.na(InterPro)) %>%
  left_join(MEROPS_df, by = c("InterPro" = "Interpro")) %>%
  filter(!is.na(MEROPS)) %>%
  mutate(subclass = case_when(str_starts(MEROPS, "A") ~ "Aspartic peptidases",
                              str_starts(MEROPS, "C") ~ "Cysteine peptidases",
                              str_starts(MEROPS, "G") ~ "Glutamic peptidases",
                              str_starts(MEROPS, "M") ~ "Metallo peptidases",
                              str_starts(MEROPS, "N") ~ "Asparagine peptidases",
                              str_starts(MEROPS, "P") ~ "Mixed peptidases",
                              str_starts(MEROPS, "S") ~ "Serine peptidases",
                              str_starts(MEROPS, "T") ~ "Threonine peptidases",
                              str_starts(MEROPS, "U") ~ "Unknown peptidases",
                              str_starts(MEROPS, "I") ~ "Inhibitor")) %>%
  filter(subclass != "Inhibitor") %>%
  group_by(subclass, Date) %>%
  summarise(summedLFQ = log2(sum(LFQ, na.rm = TRUE)),
            nExpressed = sum(!is.na(LFQ))) %>%
  mutate(summedLFQ = na_if(summedLFQ, -Inf)) %>%
  mutate(Class = "Peptidase")

# -----------------------------------------------------------------------------
# 5.2 Acidogenesis: Beta-oxidation and selected enzymes
# -----------------------------------------------------------------------------

# Beta-oxidation (https://www.genome.jp/pathway/map00071)
MP_b_oxidation <- MP_kegg %>%
  filter(Module %in% c("M00086", "M00087")) %>%
  left_join(all_metqy$METADATA %>% select(MODULE_ID, NAME_SHORT), 
            by = c("Module" = "MODULE_ID")) %>%
  rename(subclass = NAME_SHORT) %>%
  mutate(Class = "Beta-oxidation")

# Selected enzymes
other_enzymes_acidogenesis <- c(
  "K00016",  # L-lactate dehydrogenase (converts pyruvate to lactate)
  "K13951",  # alcohol dehydrogenase (converts acetaldehyde to ethanol)
  "K00656",  # pyruvate-formate lyase
  "K01568",  # pyruvate decarboxylase (converts pyruvate to acetaldehyde and CO2)
  "K01610",  # phosphoenolpyruvate carboxykinase
  "K01034",  # butyryl-CoA:acetate CoA-transferase - subunit alpha
  #"K01035",  # butyryl-CoA:acetate CoA-transferase - subunit beta
  "K19709",  # butyryl-CoA:acetate CoA-transferase
  "K00925",  # acetate kinase (converts acetyl phosphate to acetate)
  "K00929",  # butyrate kinase (involved in butyrate formation)
  "K00932",  # propionate kinase
  "K19697",  # propionate kinase
  "K01035",  # acetate CoA/acetoacetate CoA-transferase beta subunit
  "K01026"   # propionate CoA-transferase (involved in propionate formation)
)


MP_other_enzymes_acidogenesis <- MP %>%
  select(KoFamScan, Date, LFQ) %>%
  filter(!is.na(KoFamScan)) %>%
  separate_rows(KoFamScan, sep = ";") %>%
  filter(KoFamScan %in% other_enzymes_acidogenesis) %>%
  mutate(Class = "Pyruvate, lactate, butyrate, propionate metabolisms") %>%
  group_by(Date, Class) %>%
  summarise(summedLFQ = log2(sum(LFQ, na.rm = TRUE)),
            nExpressed = sum(!is.na(LFQ))) %>%
  mutate(summedLFQ = na_if(summedLFQ, -Inf)) %>%
  mutate(subclass = Class)


# -----------------------------------------------------------------------------
# 5.3 Acetogenesis: Wood-Ljungdahl pathway and selected enzymes
# -----------------------------------------------------------------------------

# Wood-Ljungdahl pathway (https://www.genome.jp/pathway/map00071)
MP_WLP <- MP_kegg %>%
  filter(Module %in% c("M00377")) %>%
  left_join(all_metqy$METADATA %>% select(MODULE_ID, NAME_SHORT), 
            by = c("Module" = "MODULE_ID")) %>%
  rename(subclass = NAME_SHORT) %>%
  mutate(Class = "Wood-Ljungdahl pathway")

# Selected enzymes (excluding pyruvate:formate lyase pathway as these are part of acidogenesis already)
other_enzymes_acetogenesis <- c(
  "K01895",  # Acetyl-CoA Synthetase (converts acetate into acetyl-CoA) also in M00357
  "K00122",  # Formate Dehydrogenase (oxidizes formate to CO2)
  "K00436",  # Hydrogenase (involved in production and consumption of hydrogen)
  "K18005",  # Hydrogenase (involved in production and consumption of hydrogen)
  "K18006",  # Hydrogenase (involved in production and consumption of hydrogen)
  "K18007",  # Hydrogenase (involved in production and consumption of hydrogen)
  "K06281",  # Hydrogenase (involved in production and consumption of hydrogen)
  "K00169"   # Pyruvate Ferredoxin Oxidoreductase (converts pyruvate to acetyl-CoA)
)

MP_other_enzymes_acetogenesis <- MP %>%
  select(KoFamScan, Date, LFQ) %>%
  filter(!is.na(KoFamScan)) %>%
  separate_rows(KoFamScan, sep = ";") %>%
  filter(KoFamScan %in% other_enzymes_acetogenesis) %>%
  mutate(Class = "(De-)hydrogenases and other oxidoreductases") %>%
  group_by(Date, Class) %>%
  summarise(summedLFQ = log2(sum(LFQ, na.rm = TRUE)),
            nExpressed = sum(!is.na(LFQ))) %>%
  mutate(summedLFQ = na_if(summedLFQ, -Inf)) %>%
  mutate(subclass = Class)

# -----------------------------------------------------------------------------
# 5.4 Methanogenesis: Four pathways
# -----------------------------------------------------------------------------
MP_methanogenesis <- MP_kegg %>%
  filter(Module %in% c("M00567", "M00357", "M00356", "M00563")) %>%
  left_join(all_metqy$METADATA %>% select(MODULE_ID, NAME_SHORT), 
            by = c("Module" = "MODULE_ID")) %>%
  rename(subclass = NAME_SHORT) %>%
  mutate(Class = "Methanogenesis")

exday <- c(-10, 3, 11, 19, 27, 35, 43)
# ===== MCR (methyl-CoM reductase) isolation=====
MCR_subunits <- c("K00399", "K00401", "K00402")

MP_MCR <- MP %>%
  select(KoFamScan, Date, LFQ) %>%
  filter(!is.na(KoFamScan)) %>%
  separate_rows(KoFamScan, sep = ";") %>%
  filter(KoFamScan %in% MCR_subunits & Date%in%c(2:8)) %>%
  mutate(Class = "MCR_subunits") %>%
  group_by(Date, Class) %>%
  summarise(summedLFQ = log2(sum(LFQ, na.rm = TRUE)),
            nExpressed = sum(!is.na(LFQ))) %>%
  mutate(summedLFQ = na_if(summedLFQ, -Inf)) %>%
  mutate(subclass = Class)


# =============================================================================
# 6. Plots summarizing the above categories into AD functions
# =============================================================================

# -----------------------------------------------------------------------------
# 6.1 Combine all of the above as one big and detailed plot
# -----------------------------------------------------------------------------

MP_all_AD <- rbindlist(list(MP_glycolysis,
                            MP_cazy,
                            MP_lipases,
                            MP_peptidases,
                            MP_other_enzymes_acidogenesis,
                            MP_b_oxidation,
                            MP_other_enzymes_acetogenesis,
                            MP_WLP,
                            MP_methanogenesis), use.names = TRUE, fill = TRUE)

# Function for normalizing a matrix to values between 0-1
normalize <- function(x) {
  x <- as.matrix(x)
  minAttr <- apply(x, 2, min)
  maxAttr <- apply(x, 2, max)
  x <- sweep(x, 2, minAttr, FUN = "-")
  x <- sweep(x, 2, maxAttr - minAttr, "/")
  attr(x, 'normalized:min') <- minAttr
  attr(x, 'normalized:max') <- maxAttr
  return(x)
}

MP_all_AD_wide <- MP_all_AD %>%
  select(-nExpressed) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(summedLFQ = ifelse(is.na(summedLFQ), rnorm(1, 10, 0.5), summedLFQ)) %>%
  pivot_wider(names_from = Date, values_from = summedLFQ) %>%
  mutate(AD = case_when(Class %in% c("Glycolysis", "CAZyme", "Lipase", "Peptidase") ~ "Hydrolysis",
                        Class %in% c("Beta-oxidation", "Pyruvate, lactate, butyrate, propionate metabolisms") ~ "Acidogenesis",
                        Class %in% c("Wood-Ljungdahl pathway", "(De-)hydrogenases and other oxidoreductases") ~ "Acetogenesis",
                        Class == "Methanogenesis" ~ "Methanogenesis")) %>%
  column_to_rownames(var = "subclass")

MP_all_AD_wide_scaled <- MP_all_AD_wide %>% 
  select(-Module, -Class, -AD) %>% 
  as.matrix()
MP_all_AD_wide_scaled <- t(normalize(t(MP_all_AD_wide_scaled)))[,1:7]
colnames(MP_all_AD_wide_scaled) <- paste('Day', exday)
# Annotation data
annot <- data.frame(EnzymeClass = factor(MP_all_AD_wide$Class, 
                                         levels = unique(MP_all_AD_wide$Class)),
                    AD = factor(MP_all_AD_wide$AD, 
                                levels = unique(MP_all_AD_wide$AD)))
rownames(annot) <- rownames(MP_all_AD_wide)

#annot_shock <- data.frame(Shock = c(0, 1, 1, 1, 1, 1, 2))
#rownames(annot_shock) <- c(rep(2:8, each = 1))

annot_colors <- list(
  EnzymeClass = setNames(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                           "#66A61E", "grey90", "#E6AB02", "#A6761D", "grey30"), 
                         levels(annot$EnzymeClass)),
  AD = setNames(c("#f2f0f7", "#cbc9e2", "#9e9ac8", "#6a51a3"), 
                levels(annot$AD)),
  Shock = c("white", "grey80", "grey40")
)

# Draw the big plot
pdf("AD_functions_heatmap_full.pdf", width = 8, height = 5)
ht_scaled <- ComplexHeatmap::pheatmap(MP_all_AD_wide_scaled,
                                      cluster_cols = FALSE,
                                      cluster_rows = FALSE,
                                      fontsize = 6.5,
                                      border_color = "grey40",
                                      annotation_row = annot,
                                      #annotation_col = annot_shock,
                                      annotation_colors = annot_colors,
                                      gaps_row = c(1, 6, 7, 14, 17, 19),
                                      heatmap_legend_param = list(
                                        title = "LFQ",
                                        labels_gp = gpar(fontsize = 7)
                                      ))
draw(ht_scaled)
dev.off()

#------------------------------------------------------------------------------
# Combined 
# ---- 0. row_split ----
n <- nrow(MP_all_AD_wide_scaled)
gaps_row <- c(1, 6, 7, 14, 17, 19)
row_split_vec <- cut(seq_len(n), breaks = c(0, gaps_row, n), labels = FALSE)
stopifnot(length(row_split_vec) == n)

# ---- 0b. pheatmap-style color palette ----
pheatmap_colors <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100)
lfq_range <- range(MP_all_AD_wide_scaled, na.rm = TRUE)
lfq_col_fun <- colorRamp2(
  seq(lfq_range[1], lfq_range[2], length.out = 100),
  pheatmap_colors
)

# ---- 1. Heatmap ----
ht <- Heatmap(
  MP_all_AD_wide_scaled,
  col             = lfq_col_fun,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  row_split       = row_split_vec,
  row_gap         = unit(1, "mm"),
  left_annotation = rowAnnotation(
    AD          = annot$AD,
    EnzymeClass = annot$EnzymeClass,
    col         = annot_colors
  ),
  row_names_gp    = gpar(fontsize = 6.5),
  column_names_gp = gpar(fontsize = 7),
  border          = TRUE,
  rect_gp         = gpar(col = "grey40", lwd = 0.5),   # 确认 col 和 lwd 在 gpar 内
  heatmap_legend_param = list(
    title          = "LFQ",
    direction      = "horizontal",
    title_position = "topcenter",
    legend_width   = unit(3, "cm"),
    labels_gp      = gpar(fontsize = 7),
    title_gp       = gpar(fontsize = 8)
  )
)

# ---- 2. Draw (修正 merge_legend -> merge_legends) ----
ht_grob <- grid.grabExpr(
  draw(ht,
       heatmap_legend_side    = "bottom",
       annotation_legend_side = "bottom",
       merge_legends           = TRUE)   # 注意是 merge_legends，复数
)

# ---- 2. Line plot ----
tipp  <- 38 + 1
tipp2 <- 44 + 1

p_line <- ggplot(MP_MCR, aes(x = exday, y = summedLFQ)) +
  geom_rect(aes(xmin = tipp - 1, xmax = tipp2 - 1,
                ymin = -Inf, ymax = Inf),
            fill = "orange", alpha = 0.8, inherit.aes = FALSE) +
  geom_line(linewidth = 1) +
  geom_point(size = 4, shape = 16) +
  labs(y = "Summed LFQ of methyl-CoM reductase", x = "Day") +
  theme_classic(base_size = 14) +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0))

# ---- 3. Stack ----
combined <- plot_grid(
  ht_grob, p_line,
  ncol        = 1,
  rel_heights = c(2, 1),
  align       = "v",
  axis        = "lr",
  labels      = c("A", "B"),
  label_size  = 16
)


# ---- 4. Save ----
ggsave(
  filename = "C:/data/Cooperation/Prodigio/Manuscript/Figure/D47/protein_level.pdf",
  plot     = combined,
  width    = 8,
  height   = 10,
  units    = "in",
  device   = cairo_pdf
)

dev.off()

windows()
combined






# -----------------------------------------------------------------------------
# 6.2 A smaller plot summarizing categories into main four AD functions
# Basis for heatmap-strip in Fig5
# -----------------------------------------------------------------------------

MP_all_AD_wide_small <- MP_all_AD %>%
  select(-nExpressed, -subclass, -Module) %>%
  mutate(AD = case_when(Class %in% c("Glycolysis", "CAZyme", "Lipase", "Peptidase") ~ "Hydrolysis",
                        Class %in% c("Beta-oxidation", "Pyruvate, lactate, butyrate, propionate metabolisms") ~ "Acidogenesis",
                        Class %in% c("Wood-Ljungdahl pathway", "(De-)hydrogenases and other oxidoreductases") ~ "Acetogenesis",
                        Class == "Methanogenesis" ~ "Methanogenesis")) %>%
  select(-Class) %>%
  mutate(summedLFQ = 2^summedLFQ) %>%
  ungroup() %>%
  group_by(AD, Date) %>%
  summarise(summedLFQ = sum(summedLFQ, na.rm = TRUE)) %>%
  mutate(summedLFQ = log2(summedLFQ)) %>%
  pivot_wider(names_from = Date, values_from = summedLFQ) %>%
  column_to_rownames(var = "AD") %>%
  slice(3, 2, 1, 4) %>%
  as.matrix()

MP_all_AD_wide_small_scaled <- t(normalize(t(MP_all_AD_wide_small)))[,1:7]
colnames(MP_all_AD_wide_small_scaled) <- paste('Day', exday)

pdf("AD_functions_heatmap_small.pdf", width = 3, height = 1)
ht_scaled_small <- ComplexHeatmap::pheatmap(MP_all_AD_wide_small_scaled,
                                            cluster_cols = FALSE,
                                            cluster_rows = FALSE,
                                            fontsize = 6.5,
                                            border_color = "grey40",
                                            #annotation_col = annot_shock,
                                            annotation_colors = annot_colors,
                                            heatmap_legend_param = list(
                                              title = "LFQ",
                                              labels_gp = gpar(fontsize = 7)
                                            ))
draw(ht_scaled_small)
dev.off()

# =============================================================================
# 7. Unsupervised analysis of all KEGG functions
# Not curated to known AD-functions as above
# Basis for Fig2c
# =============================================================================

MP_kegg_wide_export <- MP_kegg %>%
  select(-nExpressed) %>%
  ungroup() %>%
  pivot_wider(names_from = Date, values_from = summedLFQ) %>%
  left_join(all_metqy$METADATA %>% 
              select(MODULE_ID, MODULE_NAME, CLASS_I, CLASS_II, CLASS_III), 
            by = c("Module" = "MODULE_ID"))

# KEGG CLASS III annotation (broad categories)
MP_kegg_class_III <- MP_kegg_wide_export %>%
  select(-Module, -MODULE_NAME, -CLASS_I, -CLASS_II, -`9`, -`10`) %>%
  rename(`2` = `-10`, `3` = `3`, `4` = `11`, `5` = `19`, 
         `6` = `27`, `7` = `35`, `8` = `43`) %>%
  pivot_longer(cols = 1:7) %>%
  mutate(value = 2^value) %>%
  group_by(CLASS_III, name) %>%
  summarise(LFQ = sum(value, na.rm = TRUE)) %>%
  mutate(LFQ = log2(LFQ)) %>%
  group_by(CLASS_III) %>%
  filter(!all(is.infinite(LFQ))) %>%
  mutate(LFQ = ifelse(is.infinite(LFQ), rnorm(1, 14, 0.5), LFQ)) %>%
  pivot_wider(names_from = name, values_from = LFQ) %>%
  relocate(`3`, .before = `4`)

# Class II annotation of the Class III heatmap
annot_Class_II <- MP_kegg_class_III %>%
  select(CLASS_III) %>%
  left_join(MP_kegg_wide_export %>% select(CLASS_II, CLASS_III) %>% distinct(), 
            by = c("CLASS_III" = "CLASS_III"), multiple = "first") %>%
  pull(CLASS_II, CLASS_III)

annotation_colors <- c("Genetic information processing" = "#18d675",
                       "Environmental information processing" = "#005aff",
                       "Secondary metabolism" = "#8aef56",
                       "Gene set" = "#9fecca",
                       "Energy metabolism" = "#414141",
                       "Methane metabolism" = "#ff32ff",
                       "Nucleotide and amino acid metabolism" = "#ea6f56",
                       "Carbohydrate and lipid metabolism" = "#caab56",
                       "Metabolism" = "#8a6f56")

mat_MP_kegg_class_III <- MP_kegg_class_III %>%
  column_to_rownames("CLASS_III") %>%
  as.matrix()

# Plot
pdf("KEGG_unsupervised_heatmap.pdf", width = 6, height = 4)
ComplexHeatmap::pheatmap(t(normalize(t(mat_MP_kegg_class_III))),
                         cluster_cols = FALSE,
                         fontsize = 5,
                         cellwidth = 8,
                         border_color = "white",
                         treeheight_row = 20,
                         column_title = "Day",
                         column_title_side = "bottom",
                         column_title_gp = gpar(fontsize = 7),
                         heatmap_legend_param = list(
                           title = "LFQ",
                           title_gp = gpar(fontsize = 7, fontface = "bold"),
                           labels_gp = gpar(fontsize = 6)
                         ),
                         left_annotation = HeatmapAnnotation(which = "row",
                                                             Function = annot_Class_II,
                                                             col = list(Function = annotation_colors),
                                                             annotation_width = unit(0.1, "cm"),
                                                             annotation_name_gp = gpar(fontsize = 6),
                                                             annotation_legend_param = list(Function = list(
                                                               title = "Function",
                                                               title_gp = gpar(fontsize = 7, fontface = "bold"),
                                                               height = unit(1, "cm"),
                                                               labels_gp = gpar(fontsize = 6)
                                                             ))))
dev.off()

# =============================================================================
# 8. Unsupervised clustering at genus level, heatmap, two reactors
# Basis for supplementary figure
# =============================================================================

# Using object from 16S community structure above
d_OLR_Ctr_shock_plot <- amp_timeseries(d_Ctr_OLR_shock,
                                       time_variable = "Date",
                                       tax_aggregate = "Genus",
                                       tax_show = 886,  # There are 886 unique genera
                                       group_by = "samplingPoint"
)

d_OLR_Ctr_shock_plot_data <- d_OLR_Ctr_shock_plot$data %>%
  select(Genus, Value, Date, samplingPoint) %>%
  mutate(Day = as.numeric(Date - as.Date("2022-03-21"))) %>%
  mutate(Condition = case_when(samplingPoint %in% c("R1", "R2") ~ "Control",
                               samplingPoint %in% c("R3", "R4") ~ "OLR",
                               .default = "Undetermined")) %>%
  select(-Date, -samplingPoint)

# OLR
df_OLR_shock <- d_OLR_Ctr_shock_plot_data %>%
  filter(Condition == "OLR") %>%
  select(Genus, Value, Day) %>%
  group_by(Genus) %>%
  filter(!all(Value == 0)) %>%
  ungroup()

# CTR
df_CTR_shock <- d_OLR_Ctr_shock_plot_data %>%
  filter(Condition == "Control") %>%
  select(Genus, Value, Day) %>%
  group_by(Genus) %>%
  filter(!all(Value == 0)) %>%
  ungroup()

# Convert to matrix
mat_OLR_shock <- df_OLR_shock %>%
  pivot_wider(names_from = Day, values_from = Value, names_sort = TRUE) %>%
  column_to_rownames("Genus") %>%
  as.matrix()

mat_CTR_shock <- df_CTR_shock %>%
  pivot_wider(names_from = Day, values_from = Value, names_sort = TRUE) %>%
  column_to_rownames("Genus") %>%
  as.matrix()

# Normalize levels
mat_OLR_shock_norm <- t(normalize(t(mat_OLR_shock)))
mat_CTR_shock_norm <- t(normalize(t(mat_CTR_shock)))

# Draw heatmaps
ht_OLR <- ComplexHeatmap::pheatmap(mat_OLR_shock_norm,
                                   cluster_cols = FALSE,
                                   fontsize = 6,
                                   cellwidth = 3,
                                   cellheight = 2,
                                   treeheight_row = 20,
                                   column_title_gp = gpar(fontsize = 3),
                                   heatmap_legend_param = list(
                                     title = "LFQ",
                                     title_gp = gpar(fontsize = 8, fontface = "bold")),
                                   show_rownames = FALSE)

ht_CTR <- ComplexHeatmap::pheatmap(mat_CTR_shock_norm,
                                   cluster_cols = FALSE,
                                   fontsize = 6,
                                   cellwidth = 3,
                                   cellheight = 2,
                                   treeheight_row = 20,
                                   legend = FALSE,
                                   show_rownames = FALSE)

ht_combined <- ht_OLR + ht_CTR

pdf(file = paste0(wd, "./886genera_heatmaps.pdf"), width = 10, height = 25)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_combined, newpage = FALSE)
popViewport()

# Plot the second heatmap in the second column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht_CTR, newpage = FALSE)
popViewport()

popViewport()
dev.off()