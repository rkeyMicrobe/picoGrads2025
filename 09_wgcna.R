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
dat_dir = "data_out/09_wgcna/dataframes/"
fig_dir = "data_out/09_wgcna/figures/"
tab_dir = "data_out/09_wgcna/tables/"
lmm_dir = "data_out/07_lmGmatrix/dataframes/"
in_dir = "data_out/02_qiime2_asv/dataframes/"

############################################################################################
# LOAD IN DATA
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
variables = "allVars" # We are looking at all variables: NCP, POC, and PON
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

# ------------------------------------------------------------------------------
# FIND OPTIMAL POWER
## Define optimal power values
powers = c(1:10, seq(from = 12, to = 20, by = 2)) 
sft = pickSoftThreshold(df, powerVector = powers, verbose = 5) 

## Diagnostic Plots
p1 <- sft$fitIndices %>% 
  mutate(value = -sign(slope * SFT.R.sq)) %>% 
  select(Power, SFT.R.sq) %>% 
  ggplot(., aes(x = Power, y = SFT.R.sq)) +
  geom_smooth(method = "loess", se = FALSE, colour = "gray") +
  geom_hline(yintercept = 0.9, colour = "hotpink", size = 1) +
  geom_hline(yintercept = 0.8, colour = "black", size = 1) +
  geom_hline(yintercept = 1, colour = "black", size = 1) +
  geom_point(size = 2.5) +
  geom_text(aes(label = Power), vjust = -0.8, size = 4) + 
  labs(x = "Soft-Threshold Power", 
       y = "Scale-Free Model Fit (signed R^2)",
       title = "Scale Independence") +
  theme_cowplot() + theme(legend.position = "none")
p2 <- sft$fitIndices %>% 
  select(Power, mean.k.) %>% 
  ggplot(., aes(x = Power, y = mean.k.)) +
  geom_smooth(method = "loess", se = FALSE, colour = "gray") +
  geom_hline(yintercept = 0, colour = "black", size = 1) +
  geom_point(size = 2.5) +
  geom_text(aes(label = Power), vjust = -0.8, size = 4) + 
  labs(x = "Soft-Threshold Power", 
       y = "Expected Node Clusters",
       title = "Mean Connectivity") +
  theme_cowplot() + theme(legend.position = "none")  
p <- plot_grid(p1, p2); p

svg(paste0(fig_dir, "1_power_phytos_", variables, ".svg"), height = 5, width = 10); p; dev.off()

sft$fitIndices %>% 
  mutate(value = -sign(slope * SFT.R.sq)) %>% 
  select(Power, SFT.R.sq) %>% 
  filter( SFT.R.sq >= 0.8)

sft$fitIndices %>% 
  select(Power, mean.k.) %>% 
  filter(mean.k. >= .5)

############################################################################################
# MAKE THE CORRELATION NETWORK
############################################################################################

# Call your network parameter values
power = 6
modsize = 10
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

# Make Cluster Dendrogram
colors <- netwk$moduleColorsHex

svg(paste0(fig_dir, "3_dendro_phytos_", variables, "_p", power, "_", TOM_type, ".svg"), height = 5, width = 10)
plotDendroAndColors(
  dendro = netwk$dendrograms[[1]],   # The dendrogram from WGCNA
  colors = colors,                   # Apply the new colors
  groupLabels = "Module colors",     # Label for the color bar (e.g., Module colors)
  dendroLabels = FALSE,              # Do not show labels for individual genes/ASVs
  hang = 0.03,                       # Control how the leaves of the dendrogram hang
  addGuide = TRUE,                   # Add a guide to the height of the branches
  guideHang = 0.05,                  # Space between the branches and guide
  main = "Dendrogram with New Colors",  # Add a main title for the plot
  cex.colorLabels = 0.8,             # Adjust the text size of module labels
  cex.dendroLabels = 0.6,            # Adjust the text size of dendrogram labels
  marAll = c(1, 4, 3, 1),            # Adjust margins for better layout
  savePlot = FALSE                   # Set to TRUE if you want to save it
)
dev.off()

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
############################################################################################
############################ NETWORK CLUSTER ANALYSIS ######################################
############################################################################################
############################################################################################

# MODULE TRAIT CORRELATIONS 
# pull eigens
eigens <- moduleEigengenes(df, colors = netwk$moduleColors)$eigengenes
colnames(eigens) <- gsub("^ME", "", colnames(eigens))
length(colnames(eigens))
colnames(eigens)

samples <- rownames(eigens)
clusterNumber = netwk$colors
ClusterName  <- netwk$moduleColors 

modInfo <- data.frame(
  asv = names(netwk$colors),           
  number = netwk$colors,              
  module = netwk$moduleColors          
) %>%
  tibble()
modCol <- modInfo %>% 
  select(number, module) %>% distinct %>% 
  arrange(number)

# Bring in count and environment variable before input into WGCNA network
x <- data
dim(data) 
trait <- x[,c(1:3)]
dim(trait) 

# Perform Correlation between Module Eigens and NCP
var_corr <- cor(eigens, trait, use = "pairwise.complete.obs")
corr_data <- as.data.frame(var_corr) %>%
  rownames_to_column(var = "Module") %>%
  pivot_longer(cols = -Module, names_to = "Variable", values_to = "Correlation") %>% 
  filter(Module %in% unique(taxonomy$cluster)) 

cluster_order <- corr_data %>% filter(Variable == "NCP") %>% 
  arrange(desc(Correlation)) %>% pull(Module)

p <- ggplot(corr_data, aes(x = Variable, 
                           y = factor(Module, levels = cluster_order), fill = Correlation)) +
  geom_tile(color = "black") +  
  scale_fill_gradient2(low = "#FFA500", mid = "white", high = "#800080", midpoint = 0) +
  theme_cowplot() +  
  labs(title = "Correlation Matrix", x = "Traits", y = "Modules", fill = "Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right")
p
svg(paste0(fig_dir, "5_corr_phytos_", variables, ".svg"), height = 5, width = 3); p; dev.off()

############################

# Obtain p-values from pair-wise correlation test on values
results <- corr_data %>%
  mutate(p_value = map2_dbl(Module, Variable, function(r, c) {
    x <- eigens[[r]]
    y <- trait[[c]]
    if (is.numeric(x) && is.numeric(y)) {
      return(cor.test(x, y, method = "pearson")$p.value)
    } else {
      return(NA)
    }
  }))
results$Significant <- results$p_value < 0.05

x <- results %>% filter(Correlation < 0)
min(x$Correlation)
max(x$Correlation)

x <- results %>% filter(Correlation > 0) %>% filter(!Module %in% c("Yellow", "Purple"))
x


write_csv(results, paste0(dat_dir, "4_corrEigens_phytos_", 
                          variables, "_p", power, "_m", modsize, "_", TOM_type, ".txt"))

############################################################################################
# MODULE TAXONOMY COMPOSITION
############################################################################################
ncp_clust <- taxonomy %>% filter(featureID == "NCP") %>% pull(cluster) %>% print
pc_clust <- taxonomy %>% filter(featureID == "pc") %>% pull(cluster) %>% print
pn_clust <- taxonomy %>% filter(featureID == "pn") %>% pull(cluster) %>% print

taxonomy_summary <- taxonomy %>%
  filter(!featureID %in% c("NCP", "pc", "pn")) %>% 
  group_by(cluster, specific) %>%
  reframe(count = n()) %>% 
  complete(cluster, specific, fill = list(count = 0)) 
  

p <- taxonomy_summary %>% 
  mutate(count = as.numeric(count)) %>% 
  ggplot(aes(y = factor(cluster, levels = cluster_order), 
             x = count, 
             fill = "cluster")) + 
  geom_bar(stat = "identity", position = "stack", color = "#44475A") +  
  facet_wrap(~specific, scales = "free_x", nrow = 1) +
  labs(y = "Taxa Number", x = "Module", fill = "Taxa Group") +
  scale_fill_manual(values = "lightgray") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 16)) +
  facet_theme
p
svg(paste0(fig_dir, "6_comps_phytos_", variables, ".svg"), height = 5, width = 13); p; dev.off()

taxonomy_summary <- taxonomy %>%
  filter(specific == "Archaeplastida") %>% 
  group_by(cluster, species) %>%
  summarise(count = n(), .groups = "drop") %>% 
  complete(cluster, species, fill = list(count = 0)) %>% 
  mutate(species = case_when(is.na(species) ~ "Unclassified", T ~ species)) %>% 
  arrange(species, cluster)
 unique(taxonomy_summary$species)

species_order <- c("Bathycoccus_prasinos", 
           "Micromonas_commoda_A2",  
           "Chloroparvula_pacifica", 
           "Picozoa_XXXXX_sp.",
           "Chloroparvula_B2_sp.",
           "Chloroparvula_B3_sp.", 
           "Cymbomonas_tetramitiformis", 
           "Halosphaera_sp.", 
           "Pterosperma_cristatum",
           "Prasino-Clade-9_XXX_sp.",
           "Unclassified")    

p <- taxonomy_summary %>% 
  mutate(count = as.numeric(count)) %>% 
  ggplot(aes(y = factor(cluster, levels = cluster_order), 
             x = factor(species, levels = species_order),
             fill = count)) + 
  geom_tile(color = "black") +  
  scale_color_manual(values = colors) +
  scale_fill_gradient2(low = "#000000", mid = "white", high = "#80B1D3", 
                       midpoint = median(taxonomy_summary$count), name = "Taxa Number") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = paste0("Module Taxonomy of Prasinophytes in Network Power ", power))
p
svg(paste0(fig_dir, "6_comps_archae_", variables, ".svg"), height = 4.5, width = 4); p; dev.off()

############################################################################################
# SubNetwork Table
############################################################################################
library(flextable)
library(webshot2)

cluster_info <- taxonomy %>% filter(cluster %in% c("Yellow", "Purple")) %>% 
  select(-asv) %>% rename(Group = specific) %>% 
  rename(ASV_ID = featureID, Cluster = cluster,
         Phylum = phylum, Class = class, Order = order, Family = family, 
         Genus = genus, Species = species) %>% 
  select(Cluster, ASV_ID, Group, Phylum, Class, Order, Family, Genus, Species) %>% 
  arrange(Cluster, Phylum) 

# Save clusters of interest
cluster_yell <- cluster_info %>% filter(Cluster == "Yellow")
cluster_purp <- cluster_info %>% filter(Cluster == "Purple")
write.csv(cluster_yell, paste0(dat_dir, "wgcna_yellow_modASVs_p6.csv"), row.names = F)
write.csv(cluster_purp, paste0(dat_dir, "wgcna_purple_modASVs_p6.csv"), row.names = F)

# Clean up Names for Table
cluster_info <- cluster_info %>% 
  # Simple ASV names
  mutate(ASV_ID = paste0("ASV", substr(ASV_ID, 1, 5))) %>% 
  mutate(ASV_ID = case_when(str_detect(ASV_ID, "NCP") ~ "NCP", 
                            str_detect(ASV_ID, "pc") ~ "POC",
                            str_detect(ASV_ID, "pn") ~ "PON",
                            T ~ ASV_ID)
         )

cluster_info %>% filter(ASV_ID %in% c("NCP", "POC", "PON"))
cluster_colors <- c("Purple" = "#BC80BD", "Yellow" = "#FFED6F")

table <- cluster_info %>%
  flextable()  %>%
  color(i = ~Cluster == "Purple", j = "Cluster", color = "white") %>%
  color(i = ~Cluster == "Yellow", j = "Cluster", color = "black") %>%
  bg(i = ~Cluster == "Purple", j = "Cluster", bg = cluster_colors["Purple"]) %>%
  bg(i = ~Cluster == "Yellow", j = "Cluster", bg = cluster_colors["Yellow"]) %>%
  bg(i = ~Cluster == "Red", j = "Cluster", bg = cluster_colors["Red"]) %>%
  bold(i = ~Cluster %in% c("Purple", "Yellow", "Red"), j = "Cluster", bold = TRUE) %>%
  autofit()

save_as_image(table, path = paste0(tab_dir, "clusterInfoTable_purpleYellow.png"))

# End


