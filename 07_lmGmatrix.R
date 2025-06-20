# Detach loaded packages and clear environment
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls())
gc()
cat("\014")

##############################################################################################
# SECTION 7: PREPARATION AND ANALYSIS OF PHYTOPLANKTON AMPLICON AND PHYSIOCHEMICAL DATA
#
# This script processes amplicon data for eukaryotic (18S rRNA) and prokaryotic (16S rRNA) 
# communities, across multiple cruises (G1, G2, G3). It includes:
# 
# 1. Handling Amplicons: Reading, filtering, and preparing amplicon data for downstream analysis,
#    including identifying persistent taxa.
# 2. Handling Physicochemical Variables: Loading and normalizing environmental data variables 
#    such as temperature, salinity, and carbon content.
# 3. Preparing G-Matrices: Generating G-matrices from transformed amplicon data, integrating 
#    environmental metadata for multivariate mixed modeling (MLMM).
# 4. Saving Data Objects: Exporting G-matrices and environmental data objects for further analysis.
#
# The focus is on understanding the spatio-temporal distribution and persistence of key phytoplankton 
# communities and linking them to environmental drivers across the three annual surveys (2016, 2017, 2019).
##############################################################################################

# Load Libraries
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(feather)
library(compositions)
library(bestNormalize)

# Set paths
dat_dir = "data_out/07_lmGmatrix/dataframes/"
fig_dir = "data_out/07_lmGmatrix/figures/"
tab_dir = "data_out/07_lmGmatrix/tables/"
in_dir <- "data_out/02_qiime2_asv/dataframes/"

############################################################################################
#  HANDLING AMPLICONS
############################################################################################

# Fetch Amplicons
getData <- function(file = NULL, domain_label = NULL) {
  read_feather(file) %>% 
    select(sampleID, featureID, counts, cruise) %>%
    mutate(domain = domain_label) %>%
    filter(counts != 0)
}
# 18S Amplicons
g123_18s <- bind_rows(
  getData(file = paste0(in_dir, "g1_18s_master_0-200m.feather"), domain_label = "eukaryote"),
  getData(file = paste0(in_dir, "g2_18s_master_0-200m.feather"), domain_label = "eukaryote"),
  getData(file = paste0(in_dir, "g3_18s_master_0-200m.feather"), domain_label = "eukaryote")
)

# 16S Amplicons
g123_16s <- bind_rows(
  getData(file = paste0(in_dir, "g1_16s_master_0-200m.feather"), domain_label = "prokaryote"),
  getData(file = paste0(in_dir, "g2_16s_master_0-200m.feather"), domain_label = "prokaryote"),
  getData(file = paste0(in_dir, "g3_16s_master_0-200m.feather"), domain_label = "prokaryote")
)

# Fetch 18S taxonomies for those known to photosynthesize 
subset <- c("Archaeplastida", "Haptophyta", "Dinoflagellata", "Stramenopiles") 
nones <- c("Opalozoa", "Peronosporomycetes", "Bicoecea", "MOCH-2", "MOCH-4", "MOCH-5", "Developea")
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
tax_18s <- rbind(
  getTaxonomy(file = paste0(in_dir, "g1_18s_master_0-200m.feather")),
  getTaxonomy(file = paste0(in_dir, "g2_18s_master_0-200m.feather")),
  getTaxonomy(file = paste0(in_dir, "g3_18s_master_0-200m.feather"))
) %>% 
  select(featureID, phylum:specific) %>% distinct %>% print
ids_18S <- unique(tax_18s$featureID)
write_feather(tax_18s, paste0(dat_dir , "g123_18S_taxonomy.feather"))

# Fetch 16s taxonomies
subset <- c("Synechococcus", "Prochlorococcus")
getTaxonomy <- function(file = NULL) {
  x <- read_feather(file) %>% 
    select(featureID:species) %>% distinct %>% 
    rowwise() %>%
    filter(any(str_detect(c_across(domain:species), str_c(subset, collapse = "|")))) %>% 
    mutate(specific = phylum, 
           specific = case_when(grepl("Prochlorococcus", genus) ~ "Prochlorococcus",
                                grepl("Synechococcus", genus) ~ "Synechococcus", T ~ specific
           ))
}
tax_16s <- rbind(
  getTaxonomy(file = paste0(in_dir, "g1_16s_master_0-200m.feather")),
  getTaxonomy(file = paste0(in_dir, "g2_16s_master_0-200m.feather")),
  getTaxonomy(file = paste0(in_dir, "g3_16s_master_0-200m.feather"))
) %>% 
  select(featureID, phylum:specific) %>% distinct %>% print
ids_16s <- unique(tax_16s$featureID)
write_feather(tax_16s, paste0(dat_dir , "g123_16S_taxonomy.feather"))

#---------------------------------------------------------------------------------------
# Fetch persistent ids
fetchIds <- function(data = NULL){
  filtered_list <- list(
    G1 = data %>% filter(cruise == "G1") %>% pull(featureID) %>% unique(),
    G2 = data %>% filter(cruise == "G2") %>% pull(featureID) %>% unique(),
    G3 = data %>% filter(cruise == "G3") %>% pull(featureID) %>% unique()
  )
  Reduce(intersect, filtered_list)
}
persistents_IDs <- c(fetchIds(data = g123_18s), fetchIds(data = g123_16s))

#---------------------------------------------------------------------------------------
# Apply Center-Log Transformations and subset persistent ids
applyCLR <- function(data = NULL){
  x <- data %>% pivot_wider(., names_from = featureID, values_from = counts, values_fill = 0)
  x[, -(1:3)] <- x[, -(1:3)] + 1e-6
  x[,-(1:3)] <- clr(x[, -(1:3)])
  #
  return(x)
}
euks <- applyCLR(data = g123_18s)
proks <- applyCLR(data = g123_16s)

write_feather(euks,  paste0(dat_dir , "g123_18S_counts_CLR.feather"))
write_feather(proks, paste0(dat_dir , "g123_16S_counts_CLR.feather"))

applyCLR_pt2 <- function(data = NULL, idsDomain = NULL, idsPersitent = NULL){
  x <- data %>% 
    pivot_longer(., cols = -c(sampleID, cruise, domain), names_to = "featureID", values_to = "rpa_clr")
  x %>% select(-domain) %>% 
    filter(featureID %in% idsDomain) %>% 
    filter(featureID %in% idsPersitent)
}
euks <- applyCLR_pt2(data = euks, idsDomain = ids_18S, idsPersitent = persistents_IDs)
proks <- applyCLR_pt2(data = proks, idsDomain = ids_16s, idsPersitent = persistents_IDs)

# Bind the final dataframe
phytos <- rbind(euks, proks) %>% distinct 
min <- min(phytos$rpa_clr) 
phytos <- phytos %>% pivot_wider(., names_from = featureID, values_from = rpa_clr, values_fill = min)

write_feather(phytos, paste0(dat_dir , "g123_phytos_clrData.feather"))


g1 <- phytos %>% filter(cruise == "G1")
g2 <- phytos %>% filter(cruise == "G2")
g3 <- phytos %>% filter(cruise == "G3")
g123 <- rbind(g1, g2, g3)

cat(
  paste(
    paste0("Cruise 1 Sample Counts: ", length(unique(g1$sampleID))),
    paste0("Cruise 2 Sample Counts: ", length(unique(g2$sampleID))),
    paste0("Cruise 3 Sample Counts: ", length(unique(g3$sampleID))),
    sep = "\n"
  )
)

amplicons <- list(cruise1 = g1,
                  cruise2 = g2,
                  cruise3 = g3,
                  all = g123)

############################################################################################
#  HANDLING PHYSIOCHEMICAL VARIABLES
############################################################################################

# Load Variables and Select Relevant Columns
samp_info <- c("sampleID", "cruise", "latitude", "depth", "filter")
vars <- c("temp", "sal", "pn", "pc", "NCP")

 # Load and Combine Data
vars_df <- read.csv("data_out/03_physioChems/dataframes/g123_meta_physioChems.csv") %>% 
  select(all_of(samp_info), all_of(vars)) 

time_df <- read.csv("data_out/03_physioChems/dataframes/g123_smpTimes.csv") %>%
  select(sampleID, time_HST, month)

# Function to Get Sample Counts per Cruise
get_n_samples <- function(df) {
  cat(
    paste(
      paste0("Cruise 1 Sample Counts: ", length(unique(df %>% filter(cruise == "G1") %>% pull(sampleID)))),
      paste0("Cruise 2 Sample Counts: ", length(unique(df %>% filter(cruise == "G2") %>% pull(sampleID)))),
      paste0("Cruise 3 Sample Counts: ", length(unique(df %>% filter(cruise == "G3") %>% pull(sampleID)))),
      sep = "\n"
    ), "\n"
  )
}

# Update Data with Time Information and Get Sample Order
vars_df <- vars_df %>%
  left_join(time_df, by = "sampleID") %>%
  relocate(any_of(c("time_HST", "month")), .before = latitude)
get_n_samples(vars_df)

# Get Sample Order Based on Time
samp_order <- vars_df %>% arrange(time_HST) %>% pull(sampleID)

# Function to Normalize Variables
getSamps_and_normalize <- function(df, var) {
  samp_info <- c("sampleID", "cruise", "time_HST", "month", "latitude", "depth", "filter")
  df %>%
    select(all_of(samp_info), all_of(var)) %>%
    drop_na() %>%
    filter(cruise %in% cruises) %>%
    arrange(match(sampleID, samp_order)) %>%
    mutate(!!var := {
      bc_transform <- bestNormalize(!!sym(var))
      bc_transform$x.t
    })
}

# Define Cruises
cruises <- c("G1", "G2", "G3")

# Normalize Data for Each Variable
df_pn <- getSamps_and_normalize(vars_df, "pn") %>% select(-time_HST)
df_pc <- getSamps_and_normalize(vars_df, "pc") %>% select(-time_HST)
df_ncp <- getSamps_and_normalize(vars_df, "NCP") %>% select(-time_HST)


# Plot Histograms of Normalized Variables
par(mfrow = c(1, 3))
hist(df_pn$pn, main = "Total Nitrogen", xlab = "Values")
hist(df_pc$pc, main = "Total Carbon", xlab = "Values")
hist(df_ncp$NCP, main = "Net Comm Production", xlab = "Values")
par(mfrow = c(1, 1))

 # Create List of Normalized Variables
var_data <- list(
  carbon = df_pc,
  nitrogen = df_pn,
  netComm = df_ncp
)

# Function to Pull Unique Sample IDs
get_sampIDs <- function(df) {
  df %>% pull(sampleID) %>% unique()
}

# Get Sample IDs for Each Variable
var_smps <- lapply(var_data, get_sampIDs)

# Print Sample Counts for Each Variable
lapply(var_data, get_n_samples)

############################################################################################
#  PREPARING G-MATRICES
############################################################################################
# Clean R Environment
rm(list = setdiff(ls(), c(
  "directory", "tax_18s", "tax_16s", "amplicons", "ids_16s", "ids_18S",
  "samp_order", "var_data", "var_smps", "get_n_samples"
)))
gc()

 #-----------------------------------------------------------------------------------------#
# Refine Taxonomy
taxonomy <- bind_rows(tax_18s %>% select(featureID, specific),
                      tax_16s %>% select(featureID, specific)) %>% 
  distinct() %>% 
  split(.$specific)

# Prepare dataframes by (1) getting sample ids for each var condition
get_sampsIDs <- function(amplicon_df = NULL, var_df = NULL){
  clade_18S <- amplicon_df
  vars <- var_df %>% mutate(sampleID = as.character(sampleID))
  xt <- vars %>% full_join(clade_18S, by = "sampleID") %>% drop_na()
  ids <- xt %>% pull(sampleID) %>% unique
  cat(paste("Total Samples: ", length(ids)))
  return(ids)
}

# Get Sample IDs for Each Variable
smps <- list(
     carbon = get_sampsIDs(amplicon_df = amplicons$all, var_df = var_data$carbon),
     nitrogen = get_sampsIDs(amplicon_df = amplicons$all, var_df = var_data$nitrogen),
     netComm = get_sampsIDs(amplicon_df = amplicons$all, var_df = var_data$netComm)
     )

 getSamples <- function(smpList = NULL, data = amplicons$all) {
  ids <- smpList
  df <- data %>% filter(sampleID %in% ids) %>% distinct
  print(dim(df))
  return(df)
}
data_carb <- getSamples(smpList = smps$carbon, data = amplicons$all)
data_nitr <- getSamples(smpList = smps$nitrogen, data = amplicons$all)
data_netc <- getSamples(smpList = smps$netComm, data = amplicons$all)

# Split up dataframes into Taxonomy Sections
getData <- function(data = NULL, tax = taxonomy){
taxonomy = tax
df <- 
  list(Archaeplastida = data %>% select(sampleID, cruise, any_of(unique(taxonomy$Archaeplastida$featureID))),
     Dinoflagellata = data %>% select(sampleID, cruise, any_of(unique(taxonomy$Dinoflagellata$featureID))),
     Haptophyta = data %>% select(sampleID, cruise, any_of(unique(taxonomy$Haptophyta$featureID))),
     Stramenopiles = data %>% select(sampleID, cruise, any_of(unique(taxonomy$Stramenopiles$featureID))),
     Prochlorococcus = data %>% select(sampleID, cruise, any_of(unique(taxonomy$Prochlorococcus$featureID))),
     Synechococcus= data %>% select(sampleID, cruise, any_of(unique(taxonomy$Synechococcus$featureID)))
     )
return(df)
}
carbon <- getData(data = data_carb, tax = taxonomy)
nitrogen <- getData(data = data_nitr, tax = taxonomy)
netComm <- getData(data = data_netc, tax = taxonomy)

make_Gmatrix <- function(data = NULL, type = NULL){
  print(paste0(type, " --> Number of taxa = ", ncol(data) -1, ". Number of samples = ", nrow(data)))
  
  df = data
  taxa_matrix <- df[,-c(1:2)] %>% droplevels() %>% as.matrix()
  taxa_mult <- taxa_matrix %*% t(taxa_matrix)
  dist <- as.matrix(dist(taxa_matrix, method = "euclidean", diag = TRUE, upper = TRUE))
  dist <- 1-dist / max(dist, na.rm = TRUE) + diag(0.001, ncol(dist), ncol(dist))
  rownames(dist) = df$sampleID; colnames(dist) = df$sampleID
  return(dist)
}
pc_gmat <- lapply(names(carbon), function(df_name) make_Gmatrix(data = carbon[[df_name]], type = df_name))
pn_gmat <- lapply(names(nitrogen), function(df_name) make_Gmatrix(data = nitrogen[[df_name]], type = df_name))
ncp_gmat <- lapply(names(netComm), function(df_name) make_Gmatrix(data = netComm[[df_name]], type = df_name))

############################################################################################
############################################################################################
#  SAVING GMATRICE DATA OBJECT
############################################################################################
############################################################################################

gmats <- list(carbon = pc_gmat, 
              nitrogen = pn_gmat,
              netComm = ncp_gmat)
saveRDS(gmats, paste0(dat_dir, "/lmm_gMatrices_phyto_all_CLR.RDS"))

get_metaDF <- function(var_df = NULL, sample_list = NULL) {
  samples = sample_list
  var <- var_df %>% tibble() %>% 
    filter(sampleID %in% samples) %>% 
    arrange(match(sampleID, samp_order))
  
  # Ensure correct matrix structure
  var <- data.frame(
    sampleID = as.factor(var$sampleID),
    cruise = as.factor(var$cruise),
    month = as.factor(var$month),
    filter = as.factor(var$filter),
    depth = as.factor(var$depth),
    latitude = as.factor(var$latitude),
    var = var[,7]
  )
  str(var)
  return(var)
}
vars <- list(carbon = get_metaDF(var_df = var_data$carbon, sample_list = smps$carbon),
             nitrogen = get_metaDF(var_df = var_data$nitrogen , sample_list = smps$nitrogen),
             netComm = get_metaDF(var_df = var_data$netComm, sample_list = smps$netComm)
             )

saveRDS(vars, paste0(dat_dir, "/lmm_vars_phyto_all.RDS"))



