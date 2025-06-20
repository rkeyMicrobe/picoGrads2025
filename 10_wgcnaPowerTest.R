# Clean workspace and unload loaded packages
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

##############################################################################################
# SECTION 6: TAXONOMIC COMPOSITION AND WGCNA MODULE ANALYSIS

# This script processes taxonomic data and performs a weighted correlation network analysis (WGCNA)
# to investigate the relationships between microbial communities and environmental variables.
# The analysis includes module-trait correlations and the taxonomic composition of various WGCNA modules.
##############################################################################################

# Change this per each power testing 
power = 8

# LOAD PACKAGES
library("tidyverse")
library("viridis")
library("cowplot")
library("feather")
library("bestNormalize")
library("WGCNA")
library("compositions")

# Set  Themes
theme_set(theme_cowplot())
facet_theme <- theme(strip.background = element_rect(fill = "#44475A"),
                     strip.text = element_text(color = "white", face = "bold"))

# Set Paths
dat_dir = "data_out/10_wgcnaPowerTest/dataframes/"
fig_dir = "data_out/10_wgcnaPowerTest/figures/"
tab_dir = "data_out/10_wgcnaPowerTest/tables/"
lmm_dir = "data_out/07_lmGmatrix/dataframes/"
in_dir = "data_out/02_qiime2_asv/dataframes/"

############################################################################################
# LOAD IN AMPLICONS
############################################################################################

# Fetch persistent phytoplankton amplicons used in MLMM models
amps <- read_feather(paste0(lmm_dir, "g123_phytos_clrData.feather"))

# Fetch 18S taxonomies for those known to photosynthesize 
getTaxonomy <- function(file = NULL) {
  x <- read_feather(file) %>% 
    select(featureID:species) %>% distinct %>% 
    rowwise() %>%
    filter(any(str_detect(c_across(domain:species), str_c(subset, collapse = "|")))) %>% 
    filter(!class %in% nones) %>% 
    mutate(specific = kingdom,
           specific = case_when(
             kingdom == "Alveolata" ~ "Dinoflagellata",
             superkingdom == "Archaeplastida" ~ "Archaeplastida",
             TRUE ~ specific
           ))
}
subset <- c("Archaeplastida", "Haptophyta", "Dinoflagellata", "Stramenopiles") 
nones <- c("Opalozoa", "Peronosporomycetes", "Bicoecea", "MOCH-2", "MOCH-4", "MOCH-5", "Developea")
tax_18s <- rbind(getTaxonomy(file = paste0(in_dir, "g1_18s_master_0-200m.feather")),
                 getTaxonomy(file = paste0(in_dir, "g2_18s_master_0-200m.feather")),
                 getTaxonomy(file = paste0(in_dir, "g3_18s_master_0-200m.feather"))) %>% 
  select(featureID, phylum:specific) %>% distinct %>% print

# Fetch 16s taxonomies
getTaxonomy <- function(file = NULL) {
  x <- read_feather(file) %>% 
    select(featureID:species) %>% distinct %>% 
    rowwise() %>%
    filter(any(str_detect(c_across(domain:species), str_c(subset, collapse = "|")))) %>% 
    mutate(specific = phylum, 
           specific = case_when(grepl("Prochlorococcus", genus) ~ "Prochlorococcus",
                                grepl("Synechococcus", genus) ~ "Synechococcus", 
                                grepl("UCYN", genus) ~ "UCYNA",
                                T ~ specific
           ))
}
subset <- c("Synechococcus", "Prochlorococcus", "UCYN")
tax_16s <- rbind(getTaxonomy(file = paste0(in_dir, "g1_16s_master_0-200m.feather")),
                 getTaxonomy(file = paste0(in_dir, "g2_16s_master_0-200m.feather")),
                 getTaxonomy(file = paste0(in_dir, "g3_16s_master_0-200m.feather"))) %>% 
  select(featureID, phylum:specific) %>% distinct %>% print

tax <- rbind(tax_18s, tax_16s)

############################################################################################
# LOAD IN PHYSIOCHEMICAL VARIABLES
variables = "Testing" # We are looking at all variables: NCP, POC, and PON
############################################################################################

# Fetch normally distributed NCP values used in MLMM models
smps_ncp <- read.csv(paste0(lmm_dir, "boxcox_NCP.csv")) %>% tibble
smps_poc <- read.csv(paste0(lmm_dir, "boxcox_POC.csv")) %>% tibble
smps_pon <- read.csv(paste0(lmm_dir, "boxcox_PON.csv")) %>% tibble

############################################################################################
# FILTER DOWN TO ONLY SAMPLES THAT CONTAIN NCP, POC, AND PON MEASUREMENT
############################################################################################

# Find samples shared between amps and the 3 environment variables
length(amps$sampleID)
temp_smps <- Reduce(intersect, list( # Only look to see how many NCP samples
  amps$sampleID,
  smps_ncp$sampleID
))
length(amps$sampleID)
length(temp_smps)

shared_smps <- Reduce(intersect, list( # Consider NCP, POC, and PON
  amps$sampleID,
  smps_ncp$sampleID,
  smps_poc$sampleID,
  smps_pon$sampleID
))
length(shared_smps)

# Handle Var datasets
cleanUp <- function(data = NULL, smp_list = NULL, column = NULL){
  data %>% filter(sampleID %in% smp_list) %>% 
    select(sampleID, {{column}})
}
vars <- smps_ncp %>% filter(sampleID %in% shared_smps) %>% 
  left_join(., cleanUp(data = smps_poc, smp_list = shared_smps, column = pc)) %>% 
  left_join(., cleanUp(data = smps_pon, smp_list = shared_smps, column = pn))

# Handle Amplicon dataset and chronologically order samples based on var dataset
x <- amps %>% filter(sampleID %in% shared_smps) %>% as.data.frame()
dim(x)
rownames(x) <- x$sampleID
rows_to_order <- intersect(vars$sampleID, rownames(x))
x <- x[rows_to_order, ]
amps <- x %>% tibble

# Bind the Var dataset with Amp datset
data <- cbind(vars %>% select(NCP, pc, pn), 
              amps) %>% tibble %>% 
  relocate(sampleID, .before = NCP) %>% 
  select(-cruise, -sampleID)

# transpose for samples as columns
x <- as.data.frame(t(data))
colnames(x) <- as.character(amps$sampleID)
x %>% tibble

############################################################################################
# FIND OPTIMAL POWER THRESHOLD FOR NETWORK CONSTRUCTION
############################################################################################

# Transpose
df = x
df <- t(df)  
dim(df)

# check for any NAs
good_asv <- goodSamplesGenes(df)
if(good_asv$allOK == FALSE){ stop("Quality Check for NAs: Failed!") 
} else{ print("Quality Check for NAs: There are no NAs present. YAY! Keep going!") }

############################################################################################
# MAKE THE CORRELATION NETWORK
############################################################################################

# Call your network parameter values
modsize = 10
samp = "allSamp" 
correlationMatrix <- bicor(as.matrix(df), use = 'pairwise.complete.obs')
TOM <- TOMsimilarity(correlationMatrix, TOMType = 'signed')
TOM_type <- "signed" # both positive/negative directions

tom <- paste0(dat_dir, "2_tom_phytos_", variables, "_p", power, "_m", modsize, "_", TOM_type)

# Begin network generation and 'result file saves' during the process
netwk = blockwiseModules(as.matrix(df), power = power, TOMType = TOM_type, 
                         TOMDenom = "min", TOM = TOM,
                         # Block and Tree Options
                         deepSplit = 4, pamRespectsDendro = F, 
                         minModuleSize = modsize,
                         # Module Adjustments
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         # Output Options
                         numericLabels = T, verbose = 3,
                         # Save TOM files
                         saveTOMs = T, saveTOMFileBase = tom)
length(unique(netwk$colors))
unique(labels2colors(netwk$colors))

moduleColors = labels2colors(netwk$colors)

# Custom Names for network clusters
custom_names <- c(
  "yellow"      = "Yellow",
  "red"         = "Red",
  "grey"        = "Gray",
  "blue"        = "Blue",
  "turquoise"   = "Cyan",
  "green"       = "Purple",
  "purple"      = "Green",
  "salmon"      = "Clay",
  "brown"       = "Brown",
  "tan"         = "Orange",
  "pink"        = "Pink",
  "black"       = "Black",
  "magenta"     = "Plum",
  "greenyellow" = "Lime"
) #Get WGCNA color names and rename them
module_colors <- labels2colors(netwk$colors)
module_names <- custom_names[module_colors]
netwk$moduleColors <- module_names

# Custom colors for network clusters
custom_colors <- c(
  "Yellow"      = "#F1FA8C",  # Yellow - dracula palette
  "Red"         = "#FF5555",  # Red - dracula palette
  "Gray"        = "#44475A",  # Gray - dracula palette
  "Blue"        = "blue",     # Blue - Standard
  "Cyan"        = "turquoise",# Turquoise- Standard
  "Purple"      = "#BD93F9",  # Purple - dracula palette
  "Green"       = "#50FA7B",  # Green - dracula palette
  "Clay"        = "#FF9E8A",  # Clay - Custom made hex
  "Brown"       = "#4B2E2B",    # Brown - Standard
  "Orange"      = "#FFB86C",  # Orange - dracula palette
  "Pink"        = "#FF9CCF",  # Pink - Custom made hex
  "Black"       = "black",    # Black - Standard
  "Plum"        = "#A44E9A",  # Plum - Custom made hex
  "Lime"        = "#C5F779"   # Lime - Custom made hex
) # Apply chosen color palette
netwk$moduleColorsHex <- custom_colors[netwk$moduleColors]

############################################################################################
# EXPORT TO CYTOSCAPE
############################################################################################

# Pull adjacency matrix, connect to modules, and get correlation values
adjacency = adjacency(df, power = power, type = TOM_type)
TOM = TOMsimilarity(adjacency, TOMType = TOM_type)

# Select modules
moduleColors =  netwk$moduleColorsHex
modules = unique(moduleColors)
inModule = is.finite(match(moduleColors, modules))

# Pull and Clean the amplicons
nodes = colnames(df)
nodes = nodes[inModule]
nodes <- gsub("p_", "", nodes)

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(nodes, nodes)

# Extract cytoscape network to get imports
cyt = exportNetworkToCytoscape(modTOM,
                               threshold = 0.1,
                               nodeNames = nodes,
                               nodeAttr = moduleColors[inModule]
)

# Make sure your new cluster names are in replace of old names or hex colors
hex_to_name <- setNames(names(custom_colors), custom_colors)

# Clean edge and node data
taxonomy <- cyt$nodeData %>% tibble %>% 
  rename(cluster = `nodeAttr[nodesPresent, ]`,
         featureID = nodeName) %>%
  mutate(cluster = ifelse(cluster %in% names(hex_to_name),
                          hex_to_name[cluster],
                          cluster)) %>% 
  select(-altName) %>% 
  left_join(., tax, by = "featureID") %>% 
  mutate(asv = str_sub(featureID, 1, 5),
         asv = paste0(asv, "_", specific, "_", species)) %>% 
  distinct
unique(taxonomy$specific)
unique(taxonomy$cluster)
write.csv(taxonomy, 
          paste0(dat_dir, "4_WGCNA_clusterTax_", variables, "_p", 
                 power, "_m", modsize, "_", TOM_type, ".csv"))

# Generate Cytoscape edge input 
edg <- cyt$edgeData %>% select(-fromAltName, -toAltName)
edges_result <- paste0(dat_dir, "4_EDGES_phytos_", variables, "_p", power, "_m", modsize, "_", TOM_type, ".txt")
write_delim(edg, file = edges_result, delim = "\t")

# Generate Cytoscape node input 
nod <- cyt$nodeData %>% select(-altName) %>% 
  rename(cluster = `nodeAttr[nodesPresent, ]`,
         featureID = nodeName) %>% 
  mutate(cluster = ifelse(cluster %in% names(hex_to_name),
                          hex_to_name[cluster],
                          cluster)) %>% 
  left_join(., taxonomy, by = "featureID") %>% 
  mutate_all(~replace_na(., "NA")) 
nodes_result <- paste0(dat_dir, "4_NODES_phytos_", variables, "_p", power, "_m", modsize, "_", TOM_type, ".txt")
write_delim(nod, file = nodes_result, delim = "\t")

############################################################################################
# Modules of Interest -- Arch placements? 
############################################################################################

ncp_module <- taxonomy %>% filter(featureID == "NCP") %>% pull(cluster)
poc_module <- taxonomy %>% filter(featureID == "POC") %>% pull(cluster)
pon_module <- taxonomy %>% filter(featureID == "PON") %>% pull(cluster)
bathy_module <- taxonomy %>% filter(species == "Bathycoccus_prasinos") %>% pull(cluster)
modules <- c(ncp_module, bathy_module)

taxonomy_summary <- taxonomy %>%
  filter(specific == "Archaeplastida") %>% 
  group_by(cluster, species) %>%
  summarise(count = n(), .groups = "drop") %>% 
  complete(cluster, species, fill = list(count = 0)) %>% 
  mutate(species = case_when(is.na(species) ~ "Unclassified", T ~ species)) %>% 
  arrange(species, cluster)
unique(taxonomy_summary$species)

subset_order <- rev(c(
  "Chloroparvula_pacifica",
  "Bathycoccus_prasinos",
  "Micromonas_commoda_A2",
  "Picozoa_XXXXX_sp.",
  "Chloroparvula_B2_sp.",
  "Chloroparvula_B3_sp.",
  "Cymbomonas_tetramitiformis",
  "Halosphaera_sp.",
  "Prasino-Clade-9_XXX_sp.",
  "Pterosperma_cristatum",
  "Unclassified"
))

p <- taxonomy_summary %>%
  filter(cluster %in% modules) %>%
  mutate(
    cluster = if (length(ncp_module) > 0) {
      case_when(
        cluster == ncp_module ~ "Yellow",
        cluster == bathy_module ~ "Purple",
        TRUE ~ cluster
      )
    } else {
      case_when(
        cluster == bathy_module ~ "Purple",
        TRUE ~ cluster
      )
    }
  ) %>% 
  mutate(count = as.numeric(count)) %>% 
  ggplot(aes(x = NA, 
             y = factor(species, levels = subset_order),
             fill = count)) + 
  geom_tile(color = "black") + 
  facet_wrap(~cluster, ncol = 1) +
  scale_fill_gradient2(low = "#000000", mid = "white", high = "#80B1D3", 
                       midpoint = median(taxonomy_summary$count), name = "Taxa Number") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  facet_theme
p
svg(paste0(fig_dir, "6_comps_archNCP_species_p", power, ".svg"), height = 6, width = 6); p; dev.off()
