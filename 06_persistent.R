lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

##############################################################################################
# SECTION 6: Spatio-Temporal Analysis of Persistent vs. Ephemeral Amplicons

# This script processes and visualizes amplicon data for photosynthetic microbial communities, 
# focusing on both eukaryotic (18S rRNA) and prokaryotic (16S rRNA) organisms across three annual 
# surveys conducted in 2016, 2017, and 2019. The goal is to examine the spatio-temporal dynamics 
# of these communities, specifically identifying persistent vs. ephemeral taxa within the 
# photosynthetic microbial groups

##############################################################################################

# LOAD PACKAGES, COLORS, AND THEMES

library("RColorBrewer")
library("scales")
library("cowplot")
library("feather")
library("ggforce")
library("ggvenn")
library("ggpubr")
library("vegan")
library("tidyverse")

# Set themes
theme_set(theme_cowplot())
facet_theme <- theme(strip.background = element_rect(fill = "#44475A"),
                     strip.text = element_text(color = "white", face = "bold"))

# Set paths
dat_dir = "data_out/06_persistent/dataframes/"
fig_dir = "data_out/06_persistent/figures/"
tab_dir = "data_out/06_persistent/tables/"
in_dir <- "data_out/02_qiime2_asv/dataframes/"

############################################################################################
# FETCH PHYTOPLANKTON GROUPS FROM EUKARYOTE AND PROKARYOTE DATAFRAMES
############################################################################################

# Fetch 18s amplicons
getData <- function(file = NULL){
  df <- read_feather(file) %>% 
    mutate(latitude = substr(latitude, start = 1, stop = 4),
           group = kingdom) %>% 
    select(sampleID, featureID, group, counts, cruise, filter) %>% 
    mutate(domain = "eukaryote") %>% 
    filter(counts != 0)
}
g123_18s <- rbind(getData(file = paste0(in_dir, "g1_18s_master_0-200m.feather")),
                  getData(file = paste0(in_dir, "g2_18s_master_0-200m.feather")),
                  getData(file = paste0(in_dir, "g3_18s_master_0-200m.feather"))
)
# Fetch 18S taxonomies for those known to photosynthesize 
subset <- c("Archaeplastida", "Haptophyta", "Dinoflagellata", "Stramenopiles") 
nones <- c("Opalozoa", "Peronosporomycetes", "Bicoecea", "MOCH-2", "MOCH-4", "MOCH-5", "Developea")
getTaxonomy <- function(file = NULL) {
  x <- read_feather(file) %>% 
    #mutate(featureID = substr(featureID, start = 1, stop = 10)) %>% 
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
  getTaxonomy(file = paste0(in_dir, "g3_18s_master_0-200m.feather"))) %>%
  select(featureID, phylum:specific) %>% 
  distinct %>% print

euks <- g123_18s %>% 
  left_join(., tax_18s %>% 
              select(featureID, specific) %>% 
              rename(group = specific), 
            by = "featureID") %>% 
  select(-group.x) %>% 
  rename(group = group.y, rpa = counts) %>% 
  filter(!is.na(group))

write_feather(g123_18s, paste0(dat_dir, "g123_18S_counts.feather"))
write_feather(tax_18s,  paste0(dat_dir, "g123_18S_taxonomy.feather"))

# Fetch 16s amplicons
getData <- function(file = NULL) {
  read_feather(file) %>% 
    mutate(group = phylum) %>% 
    select(sampleID, featureID, group, counts, cruise, filter) %>% 
    mutate(domain = "prokaryote") %>% 
    filter(counts != 0)
}
g123_16s <- rbind(
  getData(file = paste0(in_dir, "g1_16s_master_0-200m.feather")),
  getData(file = paste0(in_dir, "g2_16s_master_0-200m.feather")),
  getData(file = paste0(in_dir, "g3_16s_master_0-200m.feather"))
)
# Fetch 16s taxonomies
subset <- c("Synechococcus", "Prochlorococcus")
getTaxonomy16S <- function(file = NULL) {
  x <- read_feather(file) %>% 
    select(featureID:species) %>% distinct %>% 
    rowwise() %>%
    mutate(specific = phylum, 
           specific = case_when(grepl("Prochlorococcus", genus) ~ "Prochlorococcus",
                                grepl("Synechococcus", genus) ~ "Synechococcus", T ~ specific))
}
tax_16s <- rbind(
  getTaxonomy16S(file = paste0(in_dir, "g1_16s_master_0-200m.feather")),
  getTaxonomy16S(file = paste0(in_dir, "g2_16s_master_0-200m.feather")),
  getTaxonomy16S(file = paste0(in_dir, "g3_16s_master_0-200m.feather"))) %>%
  select(featureID, phylum:specific) %>% 
  distinct %>% print

proks <- g123_16s %>% 
  left_join(., tax_16s %>% 
              select(featureID, specific) %>% 
              rename(group = specific), 
            by = "featureID") %>% 
  select(-group.x) %>% 
  rename(group = group.y, rpa = counts) %>% 
  filter(!is.na(group))

write_feather(g123_16s, paste0(dat_dir, "g123_16S_counts.feather"))
write_feather(tax_16s,  paste0(dat_dir, "g123_16S_taxonomy.feather"))

############################################################################################
# Persistent AMPLICONS (amplicons that are found every year we sampled)
# Comparisons at richness (# of ASVs) and RPA totals
############################################################################################

# Fetch persistence within domains
grab_persistents <- function(data = NULL){
  g1 <- data %>% filter(cruise == "G1")
  g2 <- data %>% filter(cruise == "G2")
  g3 <- data %>% filter(cruise == "G3")
  filtered_list <- list(
    G1 = g1 %>% pull(featureID) %>% unique(),
    G2 = g2 %>% pull(featureID) %>% unique(),
    G3 = g3 %>% pull(featureID) %>% unique()
  )
  persistent_ids <- Reduce(intersect, filtered_list)
  return(persistent_ids)
}
euk_ids <- grab_persistents(data = euks)
prok_ids <- grab_persistents(data = proks)
writeLines(euk_ids, paste0(dat_dir, "persistent_18S_featureIDs.txt"))
writeLines(prok_ids, paste0(dat_dir, "persistent_16S_featureIDs.txt"))

# Modify proks to only consider cyanos
tax_16s <- tax_16s %>% 
  filter(any(str_detect(c_across(phylum:species), str_c(subset, collapse = "|"))))
ids <- tax_16s %>% pull(featureID) %>% unique
proks <- proks %>% filter(featureID %in% ids)
  
# Fetch persistent phytos
phytos <- rbind(euks, proks) %>% distinct
phyto_taxonomy <- rbind(tax_18s, tax_16s) %>% distinct

g1_phytos <- phytos %>% filter(cruise == "G1")
g2_phytos <- phytos %>% filter(cruise == "G2")
g3_phytos <- phytos %>% filter(cruise == "G3")

filtered_list <- list(
  G1 = g1_phytos %>% pull(featureID) %>% unique(),
  G2 = g2_phytos %>% pull(featureID) %>% unique(),
  G3 = g3_phytos %>% pull(featureID) %>% unique()
)
persistent_ids <- Reduce(intersect, filtered_list)

pe_theme <- 
  theme_cowplot() + facet_theme +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text =  element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

# Richness
fetchCategories <- function(data = NULL, survey = NULL){
  g1_status <- data %>%
    select(featureID) %>%
    mutate(status = case_when(
      featureID %in% persistent_ids ~ "persistent",
      TRUE ~ "ephemeral"),
      cruise = survey)
}
unique(g1_phytos$sampleID)
g1_status <- fetchCategories(data = g1_phytos, survey = "G1") %>% distinct
g2_status <- fetchCategories(data = g2_phytos, survey = "G2") %>% distinct
g3_status <- fetchCategories(data = g3_phytos, survey = "G3") %>% distinct
df <- rbind(g1_status, g2_status, g3_status) %>% 
  left_join(., phyto_taxonomy, by = "featureID")
df2 <- df %>% 
  group_by(specific, status, cruise) %>%
  summarise(count = n(), .groups = 'drop') %>% 
  pivot_wider(., names_from = "specific", values_from = "count")

ordered <- c( "Synechococcus","Prochlorococcus", "Archaeplastida", "Haptophyta", "Stramenopiles", "Dinoflagellata") 

color_palette <- c("ephemeral" = "#6272A4", "persistent" = "#8BE9FD")
p1 <- ggplot(df, aes(y = factor(specific, levels = rev(ordered)), fill = status)) +
  geom_bar(position = "stack", color = "black", width = .75) +  # Add black outline
  facet_wrap(~ cruise, ncol = 1) +
  labs(title = "Persistent and Ephemeral Richness",
       y = "Taxonomic Group",
       x = "Number of Unique Taxa",
       fill = "Status") +
  scale_fill_manual(values = color_palette) + 
  pe_theme

survey_colors <- c("persistent" = "#B4F4FF", "ephemeral" = "#F8F8F2")
table_plot <- ggtexttable(df2, rows = NULL, theme = ttheme("classic"))
for (i in 1:nrow(df2)) {
  for (j in 2:ncol(df2)) { 
    table_plot <- table_cell_bg(
      table_plot, row = i + 1, column = j,
      fill = survey_colors[df2$status[i]],
      color = "black" 
    )
  }
}; table_plot
ggsave(paste0(tab_dir, "g123_PE_richness_df.png"), table_plot, width = 7, height = 2, dpi = 300)

x <- df %>%
  group_by(specific, status, cruise) %>%
  summarise(count = n(), .groups = 'drop') %>% 
  group_by(specific, cruise) %>% 
  mutate(total = sum(count),
         percent = (count / total) *100)

ggplot(x, aes(x = specific, y = percent, fill = status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(percent, 1)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  facet_wrap(~cruise) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# RPA
fetchCategories <- function(data = NULL, survey = NULL) {
  data %>%
    mutate(status = case_when(
      featureID %in% persistent_ids ~ "persistent",
      TRUE ~ "ephemeral"),
      cruise = survey) %>%
    group_by(featureID, status, cruise) %>%
    summarise(total_RPA = sum(rpa, na.rm = TRUE)) %>%  # Assuming counts represents RPA
    ungroup()
}
g1_status <- fetchCategories(data = g1_phytos, survey = "G1")
g2_status <- fetchCategories(data = g2_phytos, survey = "G2")
g3_status <- fetchCategories(data = g3_phytos, survey = "G3")
df <- bind_rows(g1_status, g2_status, g3_status) %>%
  left_join(phyto_taxonomy, by = "featureID") %>% 
  group_by(specific, status, cruise) %>% 
  reframe(total_rpa = sum(total_RPA))

# Plot RPA
p2 <- ggplot(df, aes(y = factor(specific, levels = rev(ordered)), x = total_rpa, fill = status)) +
  geom_bar(position = "stack", stat = "identity", color = "black", width = 0.75) +
  facet_wrap(~ cruise, ncol = 1) +
  labs(title = "Persistent and Ephemeral RPA", y = "Taxonomic Group", x = "Relative Percent Abundance", fill = "Status") +
  scale_fill_manual(values = color_palette) +
  pe_theme + facet_theme
p2 

x <-  bind_rows(g1_status, g2_status, g3_status) %>%
  left_join(phyto_taxonomy, by = "featureID") %>% 
  group_by(specific, status, cruise) %>%
  summarise(count = sum(total_RPA), .groups = 'drop') %>% 
  group_by(specific, cruise) %>% 
  mutate(total = sum(count),
         percent = (count / total) *100)

ggplot(x, aes(x = specific, y = percent, fill = status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(percent, 1)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  facet_wrap(~cruise) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Create Table for RPA
df2 <- df %>% pivot_wider(., names_from = "specific", values_from = "total_rpa")

table_plot_rpa <- ggtexttable(df2, rows = NULL, theme = ttheme("classic"))
for (i in 1:nrow(df2)) {
  for (j in 2:ncol(df2)) {
    table_plot_rpa <- table_cell_bg(
      table_plot_rpa, row = i + 1, column = j,
      fill = survey_colors[df2$status[i]],
      color = "black"
    )
  }
}
ggsave(paste(tab_dir, "g123_PE_RPA_df.png"), table_plot_rpa, width = 10.5, height = 2, dpi = 300)

# Combine and Save Both Plots
p <- plot_grid(p1, p2, nrow = 1)
svg(paste0(fig_dir, "/g123_PE_Cruise.svg"), height = 7, width = 10); p; dev.off()

# Save persistent vs ephemeral statuses
status_meta <- rbind(fetchCategories(data = g1_phytos, survey = "G1"),
                     fetchCategories(data = g2_phytos, survey = "G2"),
                     fetchCategories(data = g3_phytos, survey = "G3")) %>% 
  select(featureID, status) %>% distinct
write.csv(status_meta, paste0(dat_dir, "g123_persistents_asvMeta.csv"), row.names = F)

############################################################################################
# PERSISTENT VS EPHEMERAL AMPLICONS -- Region Amplicon Compositions
############################################################################################

meta <- read.csv("data_in/meta/g123_meta.csv") %>% tibble %>% 
  select(sampleID, latitude, region) %>% 
  mutate(sampleID = as.character(sampleID))
status <- status_meta

# Function to Add Metadata and Filter by Cruise
mergeAndFilter <- function(data, cruise_label) {
  data %>%
    merge(., meta, by = "sampleID") %>%
    filter(cruise == cruise_label)
}

# Merging Metadata and Filtering for Cruises (Eukaryotes and Prokaryotes)
g_phytos <- list(
  euks = list(
    G1 = mergeAndFilter(euks, "G1"),
    G2 = mergeAndFilter(euks, "G2"),
    G3 = mergeAndFilter(euks, "G3")
  ),
  proks = list(
    G1 = mergeAndFilter(proks, "G1"),
    G2 = mergeAndFilter(proks, "G2"),
    G3 = mergeAndFilter(proks, "G3")
  )
)
# Persistent IDs for Prokaryotes
filtered_list_proks <- lapply(g_phytos$proks, function(x) unique(x$featureID))
persistent_ids_p <- Reduce(intersect, filtered_list_proks)

# Subset Data Function
fetchSubset <- function(data) {
  data %>%
    select(cruise, latitude, filter, sampleID, region, group, rpa, featureID) %>%
    drop_na(group) %>%
    filter(rpa != 0)
}

# Bind All Cruise Data and Classify Persistent and Ephemeral (Eukaryotes & Prokaryotes)
g123_data <- function(data_list, persistent_ids) {
  bind_rows(lapply(data_list, fetchSubset)) %>%
    mutate(
      status = ifelse(featureID %in% persistent_ids, "Persistent", "Ephemeral")
    ) %>%
    filter(filter == "0.2")
}
g123_e <- g123_data(g_phytos$euks, persistent_ids)
g123_p <- g123_data(g_phytos$proks, persistent_ids_p)

# Classify Oceanic Regions Based on Latitude
fetchTile <- function(data, group_list) {
  # Fronts Information
  fronts <- tibble(
    cruise = c("G1", "G2", "G3"),
    sal_front = c(32.15, 32.5, 32.45),
    chl_front = c(33, 36.2, 35)
  )
  
  # Join with Fronts and Classify Regions
  data <- data %>%
    left_join(fronts, by = "cruise") %>%
    mutate(
      region = case_when(
        latitude < sal_front ~ "NPSG",
        latitude >= sal_front & latitude <= chl_front ~ "STZ",
        latitude > chl_front ~ "NTZ"
      )
    ) %>%
    select(-sal_front, -chl_front)
  
  # Summarize RPA Data
  rpa_percentage <- data %>%
    group_by(region, cruise, status) %>%
    mutate(total_rpa_region = sum(rpa, na.rm = TRUE)) %>% # Total RPA per region
    group_by(group, region, cruise, status) %>%
    reframe(
      total_rpa = sum(rpa, na.rm = TRUE),
      percentage_rpa = total_rpa / total_rpa_region * 100
    )
  
  # Plotting
  ggplot(rpa_percentage, aes(x = region, y = factor(group, levels = group_list), fill = percentage_rpa)) +
    geom_tile(color = "black") +
    geom_text(aes(label = round(percentage_rpa, 1)), size = 4, color = "white") +
    facet_grid(status ~ cruise) +
    scale_fill_viridis_c(option = "D", direction = 1) +
    labs(
      title = "Percentage of Total RPA by Kingdom, Region, and Cruise",
      subtitle = "Filter Size = 0.2uM",
      fill = "RPA Percentage"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
      strip.background = element_rect(fill = "#282A36"),
      strip.text = element_text(color = "white", size = 15, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(size = 16)
    )
}
tile_18S <- fetchTile(data = g123_e, group_list = c("Haptophyta", "Stramenopiles", "Archaeplastida", "Dinoflagellata"))
svg(paste0(fig_dir, "/g123_18S_PE_Region.svg"), height = 6, width = 8); tile_18S; dev.off()

tile_16S <- fetchTile(data = g123_p, group_list = c("Synechococcus", "Prochlorococcus", "Other Taxa"))
svg(paste0(fig_dir, "/g123_16S_PE_Region.svg"), height = 6, width = 8); tile_16S; dev.off()

## Dominants
x <- g123_e %>% 
  filter(filter == "0.2", status == "Persistent",
         group == "Archaeplastida") %>% 
  group_by(featureID, group) %>% 
  reframe(rpa = sum(rpa)) %>% 
  arrange(desc(rpa)) %>% 
  left_join(., tax_18s, by = "featureID")

############################################################################################
# PERSISTENT VS EPHEMERAL AMPLICONS -- Regional Proportions
############################################################################################

x <- phytos %>% 
  left_join(., meta, by = "sampleID") %>% 
  left_join(., status, by = "featureID") %>% 
  filter(filter == "0.2") %>% 
  group_by(cruise, region, status, group) %>% 
  reframe(rpa = sum(rpa))

get_Regional_Donuts <- function(data = x, gradient = NULL){
  donut_data <- x %>%
    filter(cruise == gradient) %>% 
    group_by(region, group, status) %>%
    summarise(total_rpa = sum(rpa), .groups = 'drop') %>%
    group_by(region, group) %>%
    mutate(status, 
           percentage = total_rpa / sum(total_rpa) * 100,
           total_rpa_region = sum(total_rpa)) %>%
    ungroup() %>% 
    mutate(region = factor(region, levels = c("NTZ", "STZ", "NPSG")))
  
  p <- ggplot(donut_data, aes(x = 2, y = percentage, fill = status)) +
    geom_bar(stat = "identity", color = "black", width = 1) +
    coord_polar(theta = "y", start = 0) +
    facet_wrap(~ region + group, ncol = 6) +
    xlim(0.5, 2.5) +  # Adjust to create space in the center
    geom_text(aes(label = round(total_rpa_region, 1), x = 0.5), 
              position = position_stack(vjust = 0.5), color = "black", size = 4) +
    labs(title = paste0("Persistent and Ephemeral RPA: ", gradient), fill = "Status") +
    scale_fill_manual(values = c("#6272A4", "#8BE9FD")) +  
    theme_void() +
    theme(legend.position = "top",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14),
          strip.background = element_rect(fill = "#44475A"),
          strip.text = element_text(color = "white", face = "bold"))
  return(p)
}
p <- get_Regional_Donuts(data = x, gradient = "G1")
svg(paste0(fig_dir, "g1_PE_Region_Donut.svg"), height = 6, width = 8); p; dev.off()
p <- get_Regional_Donuts(data = x, gradient = "G2")
svg(paste0(fig_dir, "g2_PE_Region_Donut.svg"), height = 6, width = 8); p; dev.off()
p <- get_Regional_Donuts(data = x, gradient = "G3")
svg(paste0(fig_dir, "g3_PE_Region_Donut.svg"), height = 6, width = 8); p; dev.off()

############################################################################################
# MANTEL ANALYSIS
############################################################################################
g123_data <- function(data_list, persistent_ids) {
  bind_rows(lapply(data_list, fetchSubset)) %>%
    mutate(
      status = ifelse(featureID %in% persistent_ids, "Persistent", "Ephemeral")
    ) %>%
    filter(filter == "0.2")
}
g_euks <- g123_data(g_phytos$euks, persistent_ids)
g_proks <- g123_data(g_phytos$proks, persistent_ids) 

getdf <- function(data = NULL){
  rpa_summary <- data %>%
    left_join(., meta, by = "sampleID") %>% 
    select(-region.y) %>% rename(region = region.x) %>% 
    group_by(sampleID, region, cruise, group, status) %>%
    summarise(, total_rpa = sum(rpa, na.rm = TRUE), .groups = 'drop') %>% 
    group_by(sampleID, cruise, group) %>% 
    mutate(collective = sum(total_rpa, na.rm = TRUE)) %>% 
    pivot_wider(., names_from = "status", values_from = "total_rpa") %>% 
    pivot_longer(., cols = c(-sampleID, -cruise, -group, -region), names_to = "status", values_to = "rpa")
}
df_e <- getdf(data =  g_euks)
df_p <- getdf(data =  g_proks)

combined <- rbind(df_e, df_p)

ephemeral_data <- combined %>% filter(status == "Ephemeral")
persistent_data <- combined %>% filter(status == "Persistent")
collective_data <- combined %>% filter(status == "collective")

processDfs_samplewise <- function(data = NULL){
  pseudocount <- 0.0001
  data2 <- data %>%
    pivot_wider(names_from = group, values_from = rpa, values_fill = 0) %>%
    select(-status, -cruise, -region)
  data2[is.na(data2)] <- 0
  data2[,-c(1:2)] <- data2[,-c(1,2)] + pseudocount
  return(data2)
}
# Apply to create matrices with samples as rows and groups as columns
ephemeral_matrix_samplewise <- processDfs_samplewise(data = ephemeral_data)
persistent_matrix_samplewise <- processDfs_samplewise(data = persistent_data)
collective_matrix_samplewise <- processDfs_samplewise(data = collective_data)

# Calculate Bray-Curtis distances between ephemeral vs. collective for each sample
ephemeral_distances <- mapply(function(ephemeral, collective) {
  vegdist(rbind(ephemeral, collective), method = "bray")[1]
}, as.data.frame(t(ephemeral_matrix_samplewise[,-c(1,2)])), as.data.frame(t(collective_matrix_samplewise[,-c(1,2)])))

# Calculate Bray-Curtis distances between persistent vs. collective for each sample
persistent_distances <- mapply(function(persistent, collective) {
  vegdist(rbind(persistent, collective), method = "bray")[1]
}, as.data.frame(t(persistent_matrix_samplewise[,-c(1,2)])), as.data.frame(t(collective_matrix_samplewise[,-c(1,2)])))

# Create a dataframe with the distances
distance_df <- data.frame(
  SampleID = ephemeral_matrix_samplewise$sampleID,
  Ephemeral = ephemeral_distances,
  Persistent = persistent_distances
)

wilcox <- wilcox.test(distance_df$Ephemeral, distance_df$Persistent, paired = T)
formatted_p_value <- 2.2e-16


p <- pivot_longer(distance_df, cols = c(Ephemeral, Persistent), 
                  names_to = "Comparison", values_to = "Distance") %>% 
  ggplot(., aes(x = Comparison, y = Distance, group = SampleID)) +
  geom_line(aes(group = SampleID), color = "lightgray", alpha = 0.4) +
  geom_point(aes(color = "black"), size = 3.5, alpha = 0.5) +
  geom_point(aes(color = Comparison), size = 3, alpha = 0.5) +
  scale_color_manual(values = c("Ephemeral" = "#6272A4", "Persistent" = "#8BE9FD")) +
  labs(x = "Comparison", y = "Bray-Curtis Distance",
       subtitle = paste(" Wilcoxon P-value: ", formatted_p_value)) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank()) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        plot.subtitle = element_text(size = 16))

p <- pivot_longer(distance_df, cols = c(Ephemeral, Persistent), 
                  names_to = "Comparison", values_to = "Distance") %>% 
  filter(Comparison == "Persistent")
mean(p$Distance)
max(p$Distance)
min(p$Distance)
median(p$Distance)

p <- pivot_longer(distance_df, cols = c(Ephemeral, Persistent), 
                  names_to = "Comparison", values_to = "Distance") %>% 
  ggplot(aes(x = Comparison, y = Distance, fill = Comparison)) +
  geom_violin(alpha = 0.7) +
  scale_fill_manual(values = c("Ephemeral" = "#6272A4", "Persistent" = "#8BE9FD")) +
  labs(x = "Comparison", y = "Bray-Curtis Distance",
       subtitle = paste("Wilcoxon P-value: ", formatted_p_value)) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank()) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        plot.subtitle = element_text(size = 16))
p
svg(paste0(fig_dir, "/g123_distanceTest2.svg"), height = 5.5, width = 4); p; dev.off()
