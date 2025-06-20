# Clean workspace and unload loaded packages
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

##############################################################################################
# 
##############################################################################################

library("tidyverse")
library("viridis")
library("cowplot")
library("feather")
library("scales")
library("bestNormalize")
library("compositions")
library("SpiecEasi")
library("igraph")
library("Matrix")
library("patchwork")

# Set  Themes
theme_set(theme_cowplot())
facet_theme <- theme(strip.background = element_rect(fill = "#44475A"),
                     strip.text = element_text(color = "white", face = "bold"))

# Set Paths
dat_dir = "data_out/11_speicEasi/dataframes/"
fig_dir = "data_out/11_speicEasi/figures/"
tab_dir = "data_out/11_speicEasi/tables/"
lmm_dir = "data_out/07_lmGmatrix/dataframes/"
per_dir = "data_out/06_persistent/dataframes/"
in_dir = "data_out/02_qiime2_asv/dataframes/"
wgcna_dir = "data_out/09_wgcna/dataframes/"
  
############################################################################################
# LOAD IN AMPLICONS
############################################################################################

# Fetch persistent phytoplankton amplicons used in MLMM models
amps <- read_feather(paste0(lmm_dir, "g123_phytos_clrData.feather"))

# Fetch persistent ids
persist_18S <- readLines(paste0(per_dir, "persistent_18S_featureIDs.txt"))
persist_16S <- readLines(paste0(per_dir, "persistent_16S_featureIDs.txt"))

# Load primary dataframes
cts_18s_clr <- read_feather(paste0(lmm_dir, "g123_18S_counts_CLR.feather"))
tax_18s <- read_feather(paste0(lmm_dir, "g123_18S_taxonomy.feather"))
cts_16s_clr <- read_feather(paste0(lmm_dir, "g123_16S_counts_CLR.feather"))

# Load WGCNA information
wgcna_yells <- read.csv(paste0(wgcna_dir, "wgcna_yellow_modASVs_p6.csv")) 

# Fetch persistent 18S photosynthesizers and taxonomies
euk_cts <- cts_18s_clr %>% 
  pivot_longer(., cols = -c(sampleID, cruise, domain), 
               names_to = "featureID", values_to = "clrVals") %>% 
  filter(featureID %in% persist_18S) 
tax_18s <- tax_18s %>% filter(featureID %in% unique(euk_cts$featureID))
write_feather(tax_18s, paste0(dat_dir, "g123_18S_taxonomy.feather"))

# Fetch persistent 16S photosynthesizers and taxonomies
pro_cts <- cts_16s_clr %>% 
  pivot_longer(., cols = -c(sampleID, cruise, domain), 
               names_to = "featureID", values_to = "clrVals") %>% 
  filter(featureID %in% persist_16S) 

read_clean <- function(file = NULL) {
  read_feather(file) %>%
    filter(order != "Chloroplast")  %>% 
    select(featureID:species) %>% distinct %>% 
    mutate(specific = case_when(
      str_detect(phylum, "SAR406_clade") ~ "SAR406",
      str_detect(phylum, "SAR324") ~ "SAR324",
      str_detect(class, "Alphaproteobacteria") ~ "Alphaproteobacteria",
      str_detect(class, "Gammaproteobacteria") ~ "Gammaproteobacteria",
      str_detect(genus, "Prochlorococcus") ~ "Prochlorococcus",
      str_detect(genus, "Synechococcus") ~ "Synechococcus",
      str_detect(genus, "UCYN-A") ~ "UCYN_A",
      str_detect(class, "Acidimicrobiia-A") ~ "Acidimicrobiia",
      str_detect(class, "Actinobacteria-A") ~ "Actinobacteria",
      TRUE ~ phylum)) %>% distinct
}
tax_16s <- rbind(read_clean(file = paste0(in_dir, "g1_16s_master_0-200m.feather")),
                 read_clean(file = paste0(in_dir, "g2_16s_master_0-200m.feather")), 
                 read_clean(file = paste0(in_dir, "g3_16s_master_0-200m.feather"))) %>% 
  select(featureID, phylum:specific) %>% distinct %>% 
  filter(featureID %in% unique(pro_cts$featureID))
write_feather(tax_16s, paste0(dat_dir, "g123_16S_taxonomy.feather"))

# find samples shared between domains and only consider those
year = "G3"
smps1 <-  euk_cts %>% filter(cruise == year) %>% pull(sampleID) %>% unique
smps2 <-  pro_cts %>% filter(cruise == year) %>% pull(sampleID) %>% unique
smps <- intersect(smps1, smps2)

euk_cts <- euk_cts %>% filter(sampleID %in% smps)
pro_cts <- pro_cts %>% filter(sampleID %in% smps)

# --------------------------------------------------------------------
# Grab WGCNA Archaeplastida ids and positive Braarudosphaeraceae control
modIDS <- wgcna_yells %>% 
  filter(!ASV_ID %in% c("NCP", "pc", "pn"),
         Group == "Archaeplastida") %>% 
  pull(ASV_ID) 

posC <- euk_cts %>% filter(featureID == "180a15328cf14383e0992ad4552cb976") %>% 
  pull(featureID) %>% unique 
cyanos <- pro_cts %>% filter(featureID == "772a3aa834519cf5afe559078b235b9b" |
                               featureID == "5f4ba6ed21ca4b46e699e30e401c5e45") %>% 
  pull(featureID)  %>% unique 

speicCandidates <- c(modIDS, posC, cyanos)

# Generate final dataframe structures for network analysis:
archs <- rbind(euk_cts, pro_cts) %>% filter(featureID %in% speicCandidates) %>% 
  select(-cruise, -domain) %>% 
  pivot_wider(., names_from = featureID, values_from = clrVals) %>% 
  relocate(sampleID) %>% 
  column_to_rownames(var = "sampleID") %>%
  as.matrix()
length(speicCandidates)
dim(archs)

proks <- pro_cts %>% filter(!featureID %in% speicCandidates) %>% 
  select(-cruise, -domain) %>% 
  pivot_wider(., names_from = featureID, values_from = clrVals) %>% 
  relocate(sampleID) %>% 
  column_to_rownames(var = "sampleID") %>%
  as.matrix()
dim(proks)
euks <- euk_cts %>% filter(!featureID %in% speicCandidates) %>% 
  pivot_wider(., names_from = featureID, values_from = clrVals) %>% 
  select(-cruise, -domain) %>% 
  relocate(sampleID) %>% 
  column_to_rownames(var = "sampleID") %>%
  as.matrix()
dim(euks)

all <- cbind(proks, euks)
dim(all)
taxonomy <- rbind(tax_18s %>% mutate(domain = "Euk"), 
                  tax_16s %>% mutate(domain = "Prok"))

# Loop for each column (ASV) in euks
dataframes <- list()
for (i in 1:ncol(archs)) {
  euk_column <- archs[, i, drop = FALSE]
  combined_data <- cbind(euk_column, all)
  dataframes[[colnames(archs)[i]]] <- combined_data
}

############################################################################################
# SPIEC-EASI NETWORK CONSTRUCTION
############################################################################################

gc()

get_ntwk <- function(data = NULL) {
  start_time <- Sys.time()
  
  pulsar <- list(rep.num = 50, ncores = 15)  # Define PULSAR/STAR parameters
  ntwk <- spiec.easi(data, method = "mb",    # Run SpiecEasi
                     pulsar.params = pulsar,
                     sel.criterion = 'stars',
                     lambda.min.ratio = 0.01,
                     nlambda = 20)
  
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "mins")
  cat("Model completed in", round(elapsed_time, 2), "minutes.\n")
  return(ntwk)
}
get_network_diagnostics <- function(network = NULL) {
  # Ensure object has components we need
  if (!("refit" %in% names(network)) || !("stars" %in% names(network$refit))) {
    cat("The network object does not contain the necessary 'refit$stars' data\n")
    return(NULL)
  }# Extract and print model estimation details
  if ("est" %in% names(network)) {
    cat("Model Estimation Details:\n")
    if (is.list(network$est)) {
      print(network$est)
    } else {
      cat(network$est, "\n")  
    }
  } else {
    cat("Model estimation details are not available.\n")
  }# Extract and print selection details
  if ("select" %in% names(network)) {
    cat("Selection Details:\n")
    print(network$select)
  } else {
    cat("Selection details are not available.\n")
  } # Extract edge count,and density
  edge_count <- sum(network$refit$stars != 0)
  num_nodes <- nrow(network$refit$stars)
  total_possible_edges <- num_nodes * (num_nodes - 1) / 2
  network_density <- edge_count / total_possible_edges
  cat("Edge Count:", edge_count, "\n")
  cat("Network Density:", network_density, "\n\n")
}

set.seed(1)
x = cbind(archs, all)
dim(archs)
dim(all)
dim(x)
model1 <- get_ntwk(data = x)
get_network_diagnostics(network = model1)

############################################################################################
# EXTRACT NETWORK FOR CYTOSCAPE IMPORT
############################################################################################

get_graph <- function(ntwk = NULL, ntw_input = NULL){
  bm <- symBeta(getOptBeta(ntwk), mode = "maxabs")
  diag(bm) <- 0
  weights <- Matrix::summary(t(bm))[,3]
  node_info <- colnames(ntw_input)
  graph <- adj2igraph(Matrix::drop0(getRefit(ntwk)),
                      edge.attr = list(weight = weights),
                      vertex.attr = list(name = node_info))
}
get_nodes <- function(edgeDF = NULL, taxonomyDF = taxonomy){
  nodes <- data.frame(featureID = edgeDF$from %>% unique) %>% tibble %>% 
    merge(., taxonomyDF, by = "featureID") 
  return(nodes)
}

# Get Graphs
graph <- get_graph(ntwk = model1, ntw_input = cbind(archs, all))
edg <- as_data_frame(graph, what = "edges")
nod <- get_nodes(edgeDF = edg, taxonomyDF = taxonomy)

name = "speicEasi"
write_delim(edg, file = paste0(dat_dir, "cyto_", name, "_MetComm_edges.txt"), delim = "\t")
write_delim(nod, file = paste0(dat_dir, "cyto_", name, "_MetComm_nodes.txt"), delim = "\t")

############################################################################################
# End
############################################################################################
