# Clean workspace and unload loaded packages
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

 ##############################################################################################
# LOAD UP PACKAGES AND NEEDED FILES
##############################################################################################

# Packages
library("tidyverse")
library("viridis")
library("cowplot")
library("feather")
library("flextable")

# Set  Themes
theme_set(theme_cowplot())
facet_theme <- theme(strip.background = element_rect(fill = "#44475A"),
                     strip.text = element_text(color = "white", face = "bold"))

# Set Paths
dat_dir = "data_out/12_speicEasiAnalysis/dataframes/"
fig_dir = "data_out/12_speicEasiAnalysis/figures/"
tab_dir = "data_out/12_speicEasiAnalysis/tables/"
per_dir = "data_out/06_persistent/dataframes/"
lmm_dir = "data_out/07_lmGmatrix/dataframes/"
wgcna_dir = "data_out/09_wgcna/dataframes/"
speic_dir = "data_out/11_speicEasi/dataframes/"

############################################################################################
# LOAD IN DATA
############################################################################################

# Fetch persistent ids
tax_16s <- read_feather(paste0(speic_dir, "g123_16S_taxonomy.feather"))
tax_18s <- read_feather(paste0(speic_dir, "g123_18S_taxonomy.feather"))
taxonomy <- rbind(tax_18s %>% mutate(domain = "Euk"), tax_16s %>% mutate(domain = "Prok"))

# Load up Network-related files from wgcna and speic
wgcna <- read.delim(paste0(wgcna_dir, "4_NODES_phytos_allVars_p6_m10_signed.txt")) %>% 
  tibble %>%
  rename(cluster = cluster.x) %>% select(-cluster.y) %>% 
  filter(!featureID %in% c("NCP", "pc", "pn")) %>% 
  select(featureID, cluster)

speic <-  read_delim(paste0(speic_dir, "cyto_speicEasi_MetComm_edges.txt")) %>% 
  tibble %>% 
  rename(candidate = from, associate = to) %>% 
  left_join(taxonomy, by = c("associate" = "featureID")) %>% 
  left_join(wgcna, by = c("associate" = "featureID")) %>% 
  mutate(associateID = paste0("asv", substr(associate, 1, 5), "_", class, "_", genus)) 

# Save the original interaction matrix for reference
tax_intx <- taxonomy %>% unite(full_string, phylum:species, sep = "_", remove = F, na.rm = F) %>% 
  select(featureID, full_string) %>% distinct

temp <-  read_delim(paste0(speic_dir, "cyto_speicEasi_MetComm_edges.txt")) %>% tibble %>% 
  rename(CandidateID = from, AssociateID = to) %>% 
  left_join(taxonomy, by = c("CandidateID" = "featureID")) %>% 
  left_join(wgcna, by = c("CandidateID" = "featureID")) %>% 
  unite(Candidate_Taxonomy, phylum:species, sep = "_", remove = F, na.rm = F) %>% 
  select(-c(phylum:species), -domain) %>% 
  rename(Candidate_Cluster = cluster,
         Candidate_Group = specific) %>% 
  relocate(., weight, .before = CandidateID) %>% 
  relocate(., AssociateID, .after = Candidate_Cluster) %>% 
  filter(weight >= 0) 

temp2 <- read_delim(paste0(speic_dir, "cyto_speicEasi_MetComm_edges.txt")) %>% tibble %>% 
  rename(candidate = from, associate = to) %>% 
  left_join(taxonomy, by = c("associate" = "featureID")) %>% 
  left_join(wgcna, by = c("associate" = "featureID")) %>% 
  select(-candidate, -domain) %>% 
  unite(Associate_taxonomy, phylum:species, sep = "_", remove = F, na.rm = F) %>% 
  rename(Associate_Cluster = cluster,
         Associate_Group = specific,
         asv = associate) %>% 
  select(-c(phylum:species), -weight) %>% 
  distinct

tops = 10
intx <- temp %>% 
  left_join(., temp2, by = c("AssociateID" = "asv")) %>% 
  distinct() %>% 
  group_by(CandidateID) %>%
  arrange(CandidateID, desc(weight)) %>% 
  slice_max(order_by = weight, n = tops, with_ties = FALSE) %>%  
  ungroup()
write.csv(intx,  paste0(dat_dir, "Association_Index_allVars.csv"), row.names = F)

############################################################################################
# EXPLORE first neighbor ASSOCIATES FOR YELLOW AND PURPLE MEMBERS FROM WGCNA
############################################################################################

edge_table <- speic
wgcna_index <- wgcna

# Extract positive and negative weights
pos_edges <- edge_table %>% filter(weight > 0)
neg_edges <- edge_table %>% filter(weight < 0)

# Extract Yellow and Purple cluster features, with control
yellow_features <- wgcna_index %>% filter(cluster == "Yellow") %>% pull(featureID) # n = 23
purple_features <- wgcna_index %>% filter(cluster == "Purple") %>% pull(featureID) # n = 19
control_feature <- wgcna_index %>% filter(featureID == "180a15328cf14383e0992ad4552cb976") %>% pull(featureID) 
selected_features <- unique(c(yellow_features, purple_features, control_feature))

# Create a SpiecEasi-derived index for Y-P WGCNA candidates
candidate_index <- edge_table %>% 
  filter(candidate %in% c(selected_features)) %>% 
  distinct(candidate) %>% 
  left_join(taxonomy, by = c("candidate" = "featureID")) %>% 
  left_join(wgcna_index, by = c("candidate" = "featureID")) %>% 
  mutate(candidateID = paste0("asv", substr(candidate, 1, 5), "_", class, "_", genus)) %>% 
  select(-domain) %>% 
  rename(C_Group = specific)
candidateIDs <- selected_features
Y_Ps <-  unique(c(yellow_features, purple_features))

# Create a SpiecEasi-derived index for neighbors (associates) based on candidate IDs
associate_index <- edge_table %>% 
  filter(candidate %in% c(selected_features)) %>% 
  distinct(candidate, associate, weight) %>% 
  left_join(taxonomy, by = c("associate" = "featureID")) %>% 
  left_join(wgcna_index, by = c("associate" = "featureID")) %>% 
  mutate(associateID = paste0("asv", substr(associate, 1, 5), "_", class, "_", genus)) %>% 
  rename(A_Group = specific)
associateIDs <- unique(associate_index$associate)

id_to_check <- "c0b80fec1d759345a2b73d866184e825"
list <- speic %>% pull(candidate) %>% unique %>% sort
id_to_check %in% list

# Pull out speic easi associates containing WGCNA assignments, only consider positive weighted ones
network_overlaps <- edge_table %>% filter(!is.na(cluster)) %>% 
  filter(candidate %in% Y_Ps) %>% 
  distinct(associate) %>% 
  left_join(wgcna_index, by = c("associate" = "featureID"))
unique(network_overlaps$cluster)
overlapIDs <- unique(network_overlaps$associate)

# Filter and summarise edges
edgeTable_clean <- edge_table %>% 
  filter(associate %in% associateIDs, candidate %in% candidateIDs) %>% 
  select(candidate, associate, weight, cluster) %>% 
  mutate(cluster = replace_na(cluster, "No Assignment")) %>% 
  mutate(status = case_when(weight > 0 ~ "Positive", 
                            weight < 0 ~ "Negative"))

#------------------------------------------------------------------------------
# PLOT OUT DISTRIBUTION OF POS/NEG WEIGHTS FOR CANDIDATES
temp <- edgeTable_clean %>% filter(associate %in% overlapIDs)
divider = 0
p <- ggplot(temp, aes(x = weight, fill = weight >= divider)) +
  geom_histogram(binwidth = 0.01, color = "black") +
  scale_fill_manual(values = c("FALSE" = "#F8F8F2", "TRUE" = "#44475A")) +
  labs(title = "Distribution of Associate Weights for Yellow and Purple Candidates", x = "Weight", y = "Frequency") +
  theme(legend.position = "none") 
p
svg(paste0(fig_dir, "Y-P_candidates_neighbor_weightDistributions.svg"), height = 3, width = 6); p; dev.off()

# Ranges
temp %>% filter(cluster %in% c("Yellow", "Purple")) %>% 
  summarise(
    min_weight = min(weight, na.rm = TRUE),
    max_weight = max(weight, na.rm = TRUE),
    mean_weight = mean(weight, na.rm = TRUE),
    positive_min = min(weight[weight > 0], na.rm = TRUE),
    positive_max = max(weight[weight > 0], na.rm = TRUE),
    positive_mean = mean(weight[weight > 0], na.rm = TRUE),
    negative_min = min(weight[weight < 0], na.rm = TRUE),
    negative_max = max(weight[weight < 0], na.rm = TRUE)
  ) 

# plot out Pos/Neg Weights per cluster
all_YP_associates <- edge_table %>% filter(candidate %in% Y_Ps) %>% 
  distinct(associate) %>% pull(associate)
temp <- edgeTable_clean %>% filter(associate %in% all_YP_associates) %>% 
  group_by(cluster, status) %>% 
  reframe(total_neighbors = n())
unique(temp$cluster)

cluster_order <- rev(c("No Assignment", "Red", "Cyan", "Blue", "Orange",  
                       "Brown", "Gray", "Black",  "Clay", "Purple",  "Yellow"))

p <- ggplot(temp, 
            aes(y = factor(cluster, levels = cluster_order), 
                x = total_neighbors, 
                fill = status)) +
  geom_col(position = "dodge", color = "black") +
  scale_fill_manual(values = c("Positive" = "#44475A", 
                               "Negative" = "#F8F8F2")) +
  labs(x = "Number of 1° Network Neighbors", 
       title = "WGCNA and Spiec-Easi Network Overlaps") +
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom")
p
svg(paste0(fig_dir, "Y-P_candidates_cluster_compositions.svg"), height = 3, width = 6); p; dev.off()

temp <- edgeTable_clean %>% filter(associate %in% all_YP_associates) %>% 
  select(-candidate) %>% 
  left_join(taxonomy, by = c("associate" = "featureID")) 
write.csv(temp,  paste0(dat_dir, "Y-P_relationship_index_all_w.Cluster.csv"), row.names = F)

############################################################################################
# ASV-ASV RELATIONSHIP TABLE
############################################################################################

# Make dataframe containing important info for plotting of 

 # Save Y-P connections to a file for easy access and record keeping
data <- associate_index %>% 
  mutate(Associate_taxonomy = paste(phylum, class, order, family, genus, species, sep = " | ")) %>% 
  select(-c(phylum:species)) %>%
  
  left_join(., candidate_index, by = "candidate") %>% 
  
  filter(candidate %in% candidateIDs, weight > 0) %>% 
  mutate(Candidate_taxonomy = paste(phylum, class, order, family, genus, species, sep = " | ")) %>% 
  mutate(Binomial_species_Candidate = ifelse(!is.na(species), "yes", "no")) %>% 
  arrange(Candidate_taxonomy, desc(weight)) %>% 
  select(-c(phylum:species)) %>% 
  rename(CandidateCluster = cluster.y, 
         AssociateCluster = cluster.x,
         AssociateDomain = domain,
         CandidateID_Full = candidate,
         AssociateID_Full = associate) 

data2 <- data %>% 
  select(weight, candidateID, associateID, 
         Binomial_species_Candidate, CandidateCluster, C_Group, Candidate_taxonomy, CandidateID_Full,
         AssociateCluster, AssociateDomain, A_Group, Associate_taxonomy, AssociateID_Full) %>% 
  group_by(candidateID, Binomial_species_Candidate) %>%
  ungroup %>% 
  group_by(candidateID, AssociateDomain) %>% 
  slice_head(n = 1)
  
write.csv(data2,  paste0(dat_dir, "Y-P_relationship_index_all.csv"), row.names = F)

#------------------------------------------------------------------------------
# Supplement Table

relationshipTable <- data2 %>%
  select(-CandidateID_Full, -AssociateID_Full) %>% 
  mutate(
    Associate_taxonomy = str_extract(Associate_taxonomy, "(\\|[^|]+){0,3}[^|]+$"),
    Candidate_taxonomy = str_extract(Candidate_taxonomy, "(\\|[^|]+){0,3}[^|]+$")
    ) %>% 
  mutate(
    Associate_taxonomy = substr(Associate_taxonomy, 3, nchar(Associate_taxonomy)),
    Candidate_taxonomy = substr(Candidate_taxonomy, 3, nchar(Candidate_taxonomy))
    ) %>% 
  select(-Binomial_species_Candidate) %>% 
  mutate(
    candidate_asvID = substr(candidateID, 1, 8),
    associate_asvID = substr(associateID, 1, 8),
    weight = round(weight, 2)
) %>% ungroup %>% 
  select(CandidateCluster, candidate_asvID, C_Group, Candidate_taxonomy, weight, 
         AssociateCluster, associate_asvID, AssociateDomain, A_Group, Associate_taxonomy) %>% 
  mutate(AssociateCluster = case_when(is.na(AssociateCluster) ~ "Unassigned", T ~ AssociateCluster)) %>% 
  arrange(CandidateCluster, Candidate_taxonomy)

unique(relationshipTable$CandidateCluster)
unique(relationshipTable$AssociateCluster)


flex_table <- relationshipTable %>% arrange(CandidateCluster, Candidate_taxonomy) %>% 
  flextable() %>%
  set_header_labels(
    CandidateCluster = "Candidate \n Cluster",
    candidate_asvID = "Candidate \n ASV ID",
    C_Group = "Candidate \n Group",
    Candidate_taxonomy = "Candidate Taxonomy \n (Family : Species)",
    weight = "Weight",
    AssociateCluster = "Associate \n Cluster",
    AssociateDomain = "Associate \n Domain",
    associate_asvID = "Associate \n ASV ID",
    A_Group = "Associate \n Group",
    Associate_taxonomy = "Associate Taxonomy \n (Family : Species)"
    ) %>%
  theme_vanilla() %>%
  autofit() %>% 
  width(j = "CandidateCluster", width = .25) %>%
  width(j = "candidate_asvID", width = .5) %>%
  width(j = "weight", width = .5) %>%
  width(j = "AssociateCluster", width = .5) %>%
  width(j = "associate_asvID", width = .5) %>%
  width(j = "AssociateDomain", width = .25) %>%
  bg(bg = "white", part = "body") %>%
  bg(i = ~ AssociateDomain == "Prok", bg = "#D9D9D9", part = "body"
     ) %>%
  # Candidate Colors
  bg(i = ~ CandidateCluster == "Yellow", j = 1, bg = "#FFED6F", part = "body") %>%
  bg(i = ~ CandidateCluster == "Cyan", j = 1, bg = "turquoise", part = "body") %>%
  bg(i = ~ CandidateCluster == "Purple", j = 1, bg = "#BC80BD", part = "body") %>%
  # Associate Colors
  bg(i = ~ AssociateCluster == "Yellow", j = 5, bg = "#FFED6F", part = "body") %>%
  bg(i = ~ AssociateCluster == "Cyan", j = 5, bg = "turquoise", part = "body") %>%
  bg(i = ~ AssociateCluster == "Purple", j = 5, bg = "#BC80BD", part = "body") %>%
  bg(i = ~ AssociateCluster == "Clay", j = 5, bg = "#FF9E8A", part = "body") %>%
  bg(i = ~ AssociateCluster == "Black", j = 5, bg = "black", part = "body") %>%
  color(i = ~ AssociateCluster == "Black", j = 5, color = "white", part = "body") %>% 
  bg(i = ~ AssociateCluster == "Brown", j = 5, bg = "brown", part = "body") %>%
  color(i = ~ AssociateCluster == "Brown", j = 5, color = "white", part = "body") %>% 
  color(part = "header", color = "black") %>% 
  autofit()

flex_table


 save_as_image(
  flex_table, 
  path = paste0(tab_dir, "relationship_Table.png"),
  webshot = "webshot2", 
  zoom = 2,  
  expand = 1
)
#------------------------------------------------------------------------------

# Only consider candidates with species classifications to refine down
tax <- taxonomy %>% filter(featureID %in% candidateIDs) %>% 
  mutate(species = case_when(str_detect(genus, "CC9902") ~ "Synechococcus_CC9902",
                             str_detect(species, "Braarudosphaeraceae") ~ "Braarudosphaeraceae_X",
                             T ~ species))
candidate_w_speciesName <- tax %>% filter(!is.na(species)) %>% unique %>% 
  filter(!str_detect(species, "_sp")) %>% 
  pull(featureID) %>% unique

# Identify cases of overlap where candidate serves as an associate with another candidate
dual_role_cases <- data2 %>% 
  filter(CandidateID_Full %in% associate_index$associate) %>% 
  arrange(associateID) %>% 
  filter(CandidateID_Full %in% candidate_w_speciesName) %>% 
  mutate(note = "dual_role") %>% 
  ungroup %>% 
  select(CandidateID_Full, note)

# Clean data object
relation_index <- data2 %>% filter(CandidateID_Full %in% candidate_w_speciesName) %>% 
  left_join(., dual_role_cases, by = "CandidateID_Full") %>% 
  rename(candidate = CandidateID_Full,
         associate = AssociateID_Full)

############################################################################################
# IN SITU DISTRIBUTION ANALYSIS
############################################################################################

# Load in meta, only consider G3 samples
meta <- read.csv("data_in/meta/g123_meta.csv") %>% tibble %>% 
  select(sampleID, cruise, latitude, region) %>% 
  mutate(sampleID = as.character(sampleID))
smps <- meta %>% filter(cruise == "G3") %>% 
  pull(sampleID) %>% unique 

# ????????????????? HOW DO I TREAT FILTER???

# Load in amplicon data
readInSitu <- function(countFile = NULL, taxonomyFile = NULL, cruiseSamples = NULL){
  x <- read_feather(countFile) %>% 
    select(-cruise, -filter, -group) %>% 
    left_join(., read_feather(taxonomyFile), by = "featureID") %>% 
    filter(sampleID %in% cruiseSamples) %>% 
    left_join(., meta, by = "sampleID")
  
  return(x)
}

all <- rbind(readInSitu(countFile = paste0(per_dir, "g123_16S_counts.feather"), 
                        taxonomyFile =  paste0(per_dir, "g123_16S_taxonomy.feather"),
                        cruiseSamples = smps),
             readInSitu(countFile =  paste0(per_dir, "g123_18S_counts.feather"), 
                        taxonomyFile = paste0(per_dir, "g123_18S_taxonomy.feather"),
                        cruiseSamples = smps)
)

# Subset those we care about
ids_of_interest <- c(unique(relation_index$candidate), 
                     unique(relation_index$associate))
insitu <- all %>% filter(featureID %in% ids_of_interest)


unique(insitu$sampleID)


# Continue cleaning..,
insitu <- insitu %>% 
  mutate(lat_rnd = round(latitude)) %>% 
  group_by(featureID, sampleID, lat_rnd) %>% 
  reframe(cnts_av = mean(counts[counts >0], na.rm = T),
          cnts_md = median(counts[counts >0], na.rm = T),
          domain, phylum, class, order, family, genus, species, specific, region) %>% 
  mutate(species = case_when(str_detect(genus, "CC9902") ~ "Synechococcus_CC9902",
                             str_detect(species, "Braarudosphaeraceae") ~ "Braarudosphaeraceae_X",
                             T ~ species)) %>% 
  mutate(microbeID = paste0("asv", substr(featureID, 1, 5), "_", class, "_", species)) %>% 
  select(-c(domain:species))

# Review inputs
transect_data = insitu
connections = relation_index %>% distinct %>% ungroup
candidateList = unique(connections$candidate)


length(unique(transect_data$sampleID))
for(asv in candidateList) {
  # Prepare candidate-specific data
  candidate_data <- connections %>% filter(candidate == asv)
  id_subset <- unique(c(candidate_data$candidate, 
                        candidate_data$associate))

  # Determine if the candidate ID has a dual relationship
  dual_situation_ids <- connections %>% filter(note == "dual_role") %>% pull(candidate) %>% unique
  asv_in_list <- asv %in% dual_situation_ids
  
  if (asv_in_list) {
    dual_candidateASV <- connections %>% filter(associate == asv) %>% pull(candidate)
    if (length(dual_candidateASV) > 0) {
      dual_data <- connections %>% filter(candidate == dual_candidateASV)
      ids <- unique(c(id_subset, dual_data$candidate, dual_data$associate))
    } else {
      # Nothing found → fallback
      ids <- id_subset
    }
  } else {
    ids <- id_subset
  }
  
   # Pull out all important featureIDs from transect data, classify Candidates vs. Associate
  countData <- transect_data %>% filter(featureID %in% ids) %>% 
    mutate(type = case_when(featureID %in% candidateList ~ "Candidate", TRUE ~ "Associate"))
  
  # Covert long featureIDs into the nice microbe IDs made earlier
  names <- 
    rbind(
      connections %>% select(candidate, candidateID) %>% rename(featureID = candidate, microbeID = candidateID),
      connections %>% select(associate, associateID) %>% rename(featureID = associate, microbeID = associateID)
  ) %>% unique
   
  clean_names <- names %>% filter(featureID %in% ids) %>% distinct 
  
  # Update CountData names
  countData <- countData %>% select(-microbeID) %>% 
    left_join(., clean_names, by = "featureID")
  
  if (asv_in_list) {
  asvs = c(asv, dual_candidateASV)
  asvs2 = setdiff(ids, asvs)
  } else {
    asvs = asv
    asvs2 = setdiff(ids, asvs)
  }

  name = clean_names %>% filter(featureID == asv) %>% pull(microbeID)
  
  # Generate the plot
  p <- ggplot(countData, aes(x = as.factor(lat_rnd), y = cnts_av, fill = type)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    stat_summary(fun = mean, geom = "line", 
                 aes(group = interaction(type, featureID)), 
                 color = "black", size = 0.5, linetype = "solid") +
    facet_wrap(~microbeID, scales = "free_y", ncol = 1) +
    
    scale_fill_manual(values = c("Candidate" = "#FF79C6", "Associate" = "#74838F", "Dual Role" = "#FFB86C") ) +
    scale_x_discrete(limits = as.character(26:42), labels = as.character(26:42)) +
    labs(title = paste0("Candidate: ", name),
         x = "Latitude (°N)",
         y = "RPA") +
    theme(legend.position = 'right',
          axis.text.x = element_text(angle = 45, hjust = 1)
          ) + facet_theme
  p
  
  # Save plot to SVG
  file_name <- paste0(fig_dir, "g3_speic_", name, "_inSitu.svg")
  svg(file = file_name, height = 5, width = 7.5)
  print(p)
  dev.off()
}

#-----------------------------------------------------------
# Save important files:
temp <- taxonomy %>% filter(featureID %in% ids_of_interest)
write.csv(temp,  paste0(dat_dir, "Y-P_taxonomy_candidates.csv"), row.names = F)
temp <- relation_index
write.csv(temp,  paste0(dat_dir, "Y-P_relationship_index.csv"), row.names = F)
#-----------------------------------------------------------
