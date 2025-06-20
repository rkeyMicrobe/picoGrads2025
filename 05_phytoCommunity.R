lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

##############################################################################################
# SECTION 5: SPATIO-TEMPORAL DISTRIBUTION OF PYTOPLANKTON COMMUNITIES

# This script processes and visualizes amplicon data for PHOTOSYNTHETIC eukaryotic (18S rRNA) and prokaryotic 
# (16S rRNA) communities across three annual surveys conducted in 2016, 2017, and 2019.
# The analysis focuses on relative proportion abundances (RPA), diversity 
# patterns, and the spatial distribution of major microbial groups.
##############################################################################################

# LOAD PACKAGES, COLORS, AND THEMES
library("tidyverse")
library("feather")
library("vegan")
library("cowplot")
library("viridis")

# Set themes
theme_set(theme_cowplot())
facet_theme <- theme(strip.background = element_rect(fill = "#44475A"),
                     strip.text = element_text(color = "white", face = "bold"))

# Set paths
dat_dir = "data_out/05_phytoCommunity/dataframes/"
fig_dir = "data_out/05_phytoCommunity/figures/"
tab_dir = "data_out/05_phytoCommunity/tables/"
in_dir <- "data_out/02_qiime2_asv/dataframes/"

############################################################################################
# FIND SAMPLES BORDERING THE TWO FRONTS FOR SPATIAL AMPLICON PLOTS
############################################################################################

meta <- read_csv("data_in/meta/g123_meta.csv") 
cruise_types <- c("G1", "G2", "G3")
sal_fronts <- c(32.15, 32.5, 32.45)
chl_fronts <- c(33, 36.2, 35)

# Extract samples IDs that border the each front
getSamp_byFront <- function(data, cruise_type, target_lat) {
  closest_samples <- data %>%
    filter(cruise == cruise_type) %>%
    mutate(diff = abs(latitude - target_lat)) %>%
    arrange(diff) %>%
    slice_head(n = 8)
  # Pull the upper and lower bounds
  samples_below <- closest_samples %>% filter(latitude == min(latitude)) %>% pull(sampleID)
  samples_above <- closest_samples %>% filter(latitude == max(latitude)) %>% pull(sampleID)
  # Merge them together
  c(samples_below, samples_above)
}

samp_sal <- c()
samp_chl <- c()
cruise_types <- c("G1", "G2", "G3")
for (i in seq_along(cruise_types)) {
  samp_sal <- c(samp_sal, getSamp_byFront(data = meta, cruise_type = cruise_types[i], target_lat = sal_fronts[i]))
  samp_chl <- c(samp_chl, getSamp_byFront(data = meta, cruise_type = cruise_types[i], target_lat = chl_fronts[i]))
}

# Function to get samples based on list and filter size
getSamps <- function(data, list, filterSize) {
  data %>% 
    filter(sampleID %in% list, filter == as.numeric(filterSize)) %>% 
    select(sampleID, latitude, cruise, filter) %>% 
    arrange(cruise)
}
# Define filters
small <- "0.2"
large <- "3.0"

# Fetch samples for both small and large filters
samp_sal_sm <- getSamps(data = meta, list = samp_sal, filterSize = small)
samp_chl_sm <- getSamps(data = meta, list = samp_chl, filterSize = small)
samp_sal_lg <- getSamps(data = meta, list = samp_sal, filterSize = large)
samp_chl_lg <- getSamps(data = meta, list = samp_chl, filterSize = large)

############################################################################################
# FETCH PHYTOPLANKTON GROUPS FROM EUKARYOTE AND PROKARYOTE DATAFRAMES
############################################################################################

# Get 18s Counts
getData <- function(file = NULL){
  df <- read_feather(file) %>% 
    mutate(featureID = substr(featureID, start = 1, stop = 10),
           latitude = substr(latitude, start = 1, stop = 4),
           group = kingdom) %>% 
    select(sampleID, featureID, group, counts, cruise, filter) %>% 
    mutate(domain = "eukaryote") %>% 
    filter(counts != 0)
}
g123_18s <- rbind(getData(file = paste0(in_dir, "g1_18s_master_0-200m.feather")),
                  getData(file = paste0(in_dir, "g2_18s_master_0-200m.feather")),
                  getData(file = paste0(in_dir, "g3_18s_master_0-200m.feather"))
)


# Get 18S taxonomies for photosynthetic groups
subset <- c("Archaeplastida", "Haptophyta", "Dinoflagellata", "Stramenopiles") 
nones <- c("Opalozoa", "Peronosporomycetes", "Bicoecea", "MOCH-2", "MOCH-4", "MOCH-5", "Developea")
getTaxonomy <- function(file = NULL) {
  x <- read_feather(file) %>% 
    mutate(featureID = substr(featureID, start = 1, stop = 10)) %>% 
    select(featureID:species) %>% distinct %>% 
    rowwise() %>%
    filter(any(str_detect(c_across(domain:species), str_c(subset, collapse = "|")))) %>% 
    filter(!class %in% nones) %>% 
    mutate(specific = kingdom,
           specific = case_when(
             kingdom == "Alveolata" ~ "Dinoflagellata",
             superkingdom == "Archaeplastida" ~ "Archaeplastida",
             TRUE ~ specific))
}
tax_18s <- rbind(
  getTaxonomy(file = paste0(in_dir, "g1_18s_master_0-200m.feather")),
  getTaxonomy(file = paste0(in_dir, "g2_18s_master_0-200m.feather")),
  getTaxonomy(file = paste0(in_dir, "g3_18s_master_0-200m.feather"))) %>%
  distinct() %>%
  print()
write_feather(tax_18s, paste0(dat_dir, "g123_18S_phytos_taxonomy.feather"))

# Merge count data and taxonomy data
tax <- tax_18s %>% select(featureID, specific) %>% rename(group = specific)
m <- meta %>% select(sampleID, latitude, depth, region)

euks <- g123_18s %>% 
  left_join(., tax, by = "featureID") %>% 
  select(-group.x) %>% rename(group = group.y, rpa = counts) %>% 
  drop_na(group) 

write_feather(euks, paste0(dat_dir, "g123_18S_phytos_0-200m.feather"))

euks <- euks %>% 
  select(sampleID, cruise, filter, domain, group, rpa) %>% 
  merge(., m, by = "sampleID")

#-------------------------------------------------------------------------------
# Get 16s counts
getData <- function(file = NULL) {
  read_feather(file) %>% 
    mutate(featureID = substr(featureID, start = 1, stop = 10),
           group = phylum) %>% 
    select(sampleID, featureID, group, counts, cruise, filter) %>% 
    mutate(domain = "prokaryote") %>% 
    filter(counts != 0)
}
g123_16s <- rbind(
  getData(file = paste0(in_dir, "g1_16s_master_0-200m.feather")),
  getData(file = paste0(in_dir, "g2_16s_master_0-200m.feather")),
  getData(file = paste0(in_dir, "g3_16s_master_0-200m.feather"))
)

# Get 16s taxonomies
subset <- c("Synechococcus", "Prochlorococcus")
getTaxonomy16S <- function(file = NULL) {
  x <- read_feather(file) %>% 
    mutate(featureID = substr(featureID, start = 1, stop = 10)) %>% 
    select(featureID:species) %>% distinct %>% 
    rowwise() %>%
    filter(any(str_detect(c_across(domain:species), str_c(subset, collapse = "|")))) %>% 
    mutate(specific = phylum, 
           specific = case_when(grepl("Prochlorococcus", genus) ~ "Prochlorococcus",
                                grepl("Synechococcus", genus) ~ "Synechococcus", T ~ specific))
}
tax_16s <- rbind(
  getTaxonomy16S(file = paste0(in_dir, "g1_16s_master_0-200m.feather")),
  getTaxonomy16S(file = paste0(in_dir, "g2_16s_master_0-200m.feather")),
  getTaxonomy16S(file = paste0(in_dir, "g3_16s_master_0-200m.feather"))) %>%
  distinct() %>%
  print()
write_feather(tax_16s, paste0(dat_dir, "g123_16S_phytos_taxonomy.feather"))

# Merge count data and taxonomy data
tax <- tax_16s %>% select(featureID, specific) %>% rename(group = specific)
m <- meta %>% select(sampleID, latitude, depth, region)

proks <- g123_16s %>% 
  left_join(., tax, by = "featureID") %>% 
  select(-group.x) %>% rename(group = group.y, rpa = counts) %>% 
  drop_na(group)

write_feather(proks, paste0(dat_dir, "g123_16S_phytos_0-200m.feather"))

proks <- proks %>% 
  select(sampleID, cruise, filter, domain, group, rpa)  %>% 
  merge(., m, by = "sampleID")

# Checking column name to see if they match for next section 
colnames(euks); colnames(proks)

############################################################################################
# SPATIAL MAPPING - HEAT MAP
############################################################################################

# Function to generate heatmap of spatial data
fetch_heatmap <- function(data = NULL, size = NULL,
                          salinity_ids = NULL, 
                          chlorophyll_ids = NULL){
  x <- data %>% 
    filter(filter == size) %>% 
    group_by(sampleID, filter, latitude, group, depth) %>% 
    summarise(rpa = sum(rpa))
  
  order <- x %>%
    group_by(sampleID, latitude, filter, group, depth) %>% 
    reframe(rpa = mean(rpa)) %>% 
    filter(filter ==  size) %>% 
    group_by(latitude) %>% 
    arrange(desc(latitude), desc(depth)) %>% 
    pull(sampleID) %>% unique
  
  plot <- x %>% 
    ggplot(.) + 
    geom_tile(aes(x = group, y = factor(sampleID, levels = rev(order)), 
                  fill = rpa), linewidth = 0.1) +
    scale_fill_viridis_c(option = "viridis") +
    scale_x_discrete(expand = c(0, 0)) +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          #axis.text.x = element_blank(),
          panel.background = element_rect("black")) +
    geom_hline(yintercept = match(salinity_ids, rev(order)), 
               color = "#FF7F00", linetype = "solid") +
    geom_hline(yintercept = match(chlorophyll_ids, rev(order)), 
               color = "#50FA7B", linetype = "solid")
  plot
}
fetch_frontIDs <- function(df, cruise_type = survey) {
  df <- df %>% filter(cruise == cruise_type) %>% group_by(latitude)
  ids <- c(df %>% head(1) %>% pull(sampleID), df %>% tail(1) %>% pull(sampleID))
  return(ids)
}
# Function to get spatial maps
fetch_spatials <- function(data1, data2, size = NULL, 
                           sal_front = NULL, chl_front = NULL, directory = "3.0_plots/") {
  surveys <- c("G1", "G2", "G3")
  fractions <- rep(size, length(surveys))
  
  # Perform Subsets
  g1_eu <- data1 %>% filter(cruise == "G1")
  g2_eu <- data1 %>% filter(cruise == "G2")
  g3_eu <- data1 %>% filter(cruise == "G3")
  g1_cy <- data2 %>% filter(cruise == "G1")
  g2_cy <- data2 %>% filter(cruise == "G2")
  g3_cy <- data2 %>% filter(cruise == "G3")
  
  euk_data_list <- list(g1_eu, g2_eu, g3_eu)
  cya_data_list <- list(g1_cy, g2_cy, g3_cy)
  
  for (i in seq_along(surveys)) {
    survey <- surveys[i]
    fraction <- fractions[i]
    euk_data <- euk_data_list[[i]]
    cya_data <- cya_data_list[[i]]
    
    ids <- intersect(euk_data$sampleID, cya_data$sampleID)
    euk_data <- euk_data %>% filter(sampleID %in% ids)
    cya_data <- cya_data %>% filter(sampleID %in% ids)
    
    sal_samps <- fetch_frontIDs(df = sal_front, cruise_type = survey)
    chl_samps <- fetch_frontIDs(df = chl_front, cruise_type = survey)
    
    plot_cr <- plot_grid(
      fetch_heatmap(data = cya_data, size = fraction, salinity_ids = sal_samps, chlorophyll_ids = chl_samps), 
      fetch_heatmap(data = euk_data, size = fraction, salinity_ids = sal_samps, chlorophyll_ids = chl_samps), 
      rel_widths = c(0.75, 1.25),
      labels = c('Cyanobacteria', 'Eukaryotes')
    )
    file_name <- paste0(fig_dir, survey, "_", fraction, "um_spatials_all.svg")
    svg(file_name, height = 6, width = 6)
    print(plot_cr)
    dev.off()
  }
}

# Make spatial heatmaps and save

euks_surface <- euks %>% filter(depth <= 15) 
proks_surface <- proks %>% filter(depth <= 15) 

size <- 3.0
sal_front <- samp_sal_lg
chl_front <- samp_chl_lg
fetch_spatials(data1 = euks_surface, data2 = proks_surface, size = size, 
               sal_front = sal_front, chl_front = chl_front)

x <- meta %>% select(sampleID, depth) %>% mutate(sampleID = as.character(sampleID))
euks %>% select(sampleID) %>% distinct %>% 
  left_join(., x, by = "sampleID")

size <- 0.2
sal_front <- samp_sal_sm
chl_front <- samp_chl_sm
fetch_spatials(data1 = euks_surface, data2 = proks_surface, size = size, 
               sal_front = sal_front, chl_front = chl_front)

############################################################################################
# SPATIAL MAPPING - LINES/POINTS
############################################################################################
cruise_types <- c("G1", "G2", "G3")
sal_fronts <- c(32.15, 32.5, 32.45)
chl_fronts <- c(33, 36.2, 35)
phytos <- rbind(euks, proks)

# Function to plot, RPA every 2 degrees latitude
get_plot <- function(data = NULL, gradient = NULL, palette = NULL) {
  z <- data %>%
    filter(cruise == gradient) %>%
    filter(filter == "0.2") %>% 
    mutate(lat = round(latitude / 2) * 2) %>%
    group_by(lat, group) %>%
    summarise(rpa = sum(rpa)) %>%
    arrange(desc(lat)) %>% 
    group_by(lat) %>%
    mutate(total = sum(rpa, na.rm = TRUE),
           perc = (rpa / total) * 100) %>% ungroup() %>%
    complete(lat = seq(20, 44, by = 2), 
             fill = list(group = "Unsampled Range", rpa = 0)) %>%
    ungroup() %>% 
    mutate(group = if_else(is.na(group), "Unsampled", group),
           perc = if_else(is.na(perc), 100, perc)) %>% 
    filter(group != "Unsampled Range")
 
   ggplot(z, aes(x = lat, y = rpa, color = group, group = group)) +
    geom_line(size = 2.5, aes(color = group), linetype = "solid") + 
     geom_hline(yintercept = 0, linetype = "dashed", color = "#44475A") +
    geom_point(size = 1.75, shape = 21, aes(fill = group), color = "#282A36", stroke = .5) + 
    scale_fill_manual(values = palette) + 
    scale_color_manual(values = palette) +
     xlim(20, 45) +
    labs(x = "Latitude °N", y = "RPA") + 
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18))
}
# Eukaryotic phytoplankton
euk_colors = c("#49C5B1", "#ff8686", "#DAA520", "#8BE9FD", "#44475A")
p1 <- get_plot(data = euks, gradient = "G1", palette = euk_colors); p1
p2 <- get_plot(data = euks, gradient = "G2", palette = euk_colors); p2
p3 <- get_plot(data = euks, gradient = "G3", palette = euk_colors); p3
p <- plot_grid(p1, p2, p3, ncol = 1)
p 
svg(paste0(fig_dir, "g123_18S_3_unfrac_lineSpatials_noLegend.svg"), height = 6, width = 4); print(p); dev.off()
#svg(paste0(plot_dir, "g123_18S_3_unfrac_dotSpatials_wLegend.svg"), height = 20, width = 6); print(p); dev.off()

# Prokaryotic phytoplankton
prok_colors = c("#bd93f9", "#ffb86c", "#44475A")
p1 <- get_plot(data = proks, gradient = "G1", palette = prok_colors)
p2 <- get_plot(data = proks, gradient = "G2", palette = prok_colors)
p3 <- get_plot(data = proks, gradient = "G3", palette = prok_colors)
p <- plot_grid(p1, p2, p3, ncol = 1); p 
svg(paste0(fig_dir, "g123_16S_3_unfrac_lineSpatials_noLegend.svg"),height = 6, width = 4); print(p); dev.off()

###########################################################################################3
# Box Plots
###########################################################################################3

phytos <- rbind(euks, proks)

getPlotStats <- function(data = NULL, gradient = NULL){
  df <- data %>% mutate(sampleID = as.character(sampleID))
  
  df <- df %>%
    filter(cruise == gradient) %>% 
    group_by(sampleID, group, filter) %>% 
    reframe(region, rpa = sum(rpa)) %>% 
    distinct
  
  order <-c("Archaeplastida", "Dinoflagellata", "Haptophyta", "Stramenopiles",
            "Prochlorococcus", "Synechococcus")
  
  
  df$region <- factor(df$region, levels = rev(c("NPSG", "STZ", "NTZ")))
  df$group <- factor(df$group, levels = order)
  color_palette <- c("3" = "#44475A", "0.2" = "#D49A3A")
  
  # Extract statistics for the boxplot (mean, max, min)
  boxplot_stats <- df %>% 
    group_by(region, group, filter) %>% 
    summarise(
      mean_rpa = mean(rpa) *100,
      max_rpa = max(rpa) *100,
      min_rpa = min(rpa) *100
    ) %>% 
    arrange(desc(mean_rpa))
  
  # Plot the data
  plot <- ggplot(df, aes(x = group, y = rpa, color = as.factor(filter))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3, width = .4, 
                 position = position_dodge(width = .9)) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.6) +  
    facet_wrap(~ region, ncol = 1) +
    ylim(0, 1) +
    scale_color_manual(values = color_palette) +
    labs(y = "Relative Percent Abundance", color = "Filter Size (um)",
         title = paste0("Survey Year: ", gradient)) +
    theme(axis.text.x = element_text(size = 15, angle = 45, hjust = .9, vjust = 1), 
          axis.text.y = element_text(size = 15), 
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          plot.title = element_text(size = 20),
          strip.text = element_text(size = 15),
          legend.position = "top") +
    facet_theme
  
  return(list(plot = plot, stats = boxplot_stats))
}

p1 <- getPlotStats(data = phytos, gradient = "G1")
p2 <- getPlotStats(data = phytos, gradient = "G2")
p3 <- getPlotStats(data = phytos, gradient = "G3")

p <- plot_grid(p1$plot, p2$plot, p3$plot, nrow = 1)
p
svg(paste0(fig_dir, "g123_Phyto_filterCompares.svg"), height = 6, width = 12); p; dev.off()

# Look at stats
write.csv(p1$stats, paste0(dat_dir, "g1_16s_boxPlotStats_metaComm.csv"), row.names = F)
write.csv(p2$stats, paste0(dat_dir, "g2_16s_boxPlotStats_metaComm.csv"), row.names = F)
write.csv(p3$stats, paste0(dat_dir, "g3_16s_boxPlotStats_metaComm.csv"), row.names = F)

fold_changes <- p1$stats %>%
  select(region, group, filter, mean_rpa) %>%
  pivot_wider(names_from = filter, values_from = mean_rpa, names_prefix = "filter_") %>%
  mutate(fold_change = filter_0.2 / filter_3) %>% 
  select(-filter_0.2, -filter_3) %>% 
  pivot_wider(., names_from = region, values_from = fold_change)

 ############################################################################################
# TOTAL RPA FILTER FOLD CHANGE
############################################################################################

phytos <- rbind(euks, proks)

data <- phytos %>% 
  filter(rpa != 0) %>%
  group_by(sampleID, group) %>% 
  reframe(cruise, filter, latitude, region,
          rpa = sum(rpa)) %>% 
  distinct %>% ungroup %>%  
  group_by(cruise, group, filter) %>%
  summarize(total_rpa = sum(rpa), .groups = 'drop') %>% 
  mutate(cruise_filter = paste0(cruise, "_", filter),
         group_filter = paste0(group, "_", filter))

order <- c("Archaeplastida", "Dinoflagellata", "Haptophyta", "Stramenopiles", 
           "Prochlorococcus", "Synechococcus")
data <- data %>% mutate(group = factor(group, levels = order))

# Acquire total rpa for each filter across years
all_colors <- c("Archaeplastida_0.2" = "#49C5B1", "Archaeplastida_3" = "#317A6A",
                "Haptophyta_0.2" = "#FFD700", "Haptophyta_3" = "#BFA500",
                "Dinoflagellata_0.2" = "#DA4949", "Dinoflagellata_3" = "#A83636",
                "Stramenopiles_0.2" = "#ADECF9", "Stramenopiles_3" = "#7DB2C9",
                "Prochlorococcus_0.2" = "#7570B3", "Prochlorococcus_3" = "#4D2C7A",
                "Synechococcus_0.2" = "#ff7f00", "Synechococcus_3" = "#B35900")

p <- data %>% 
  ggplot(., aes(y = cruise_filter, x = total_rpa, fill = group_filter)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  geom_text(aes(label = sprintf("%.2f", total_rpa)), 
            position = position_dodge(width = 0.8), 
            hjust = -0.5, size = 5, color = "black") +
  facet_wrap(~ group) +
  xlim(0, 40) +
  scale_fill_manual(values = all_colors) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = .5, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 15),
    strip.background = element_rect(fill = "#282A36"),
    strip.text = element_text(color = "white", size = 20, face = "bold"),
  )
p
svg(paste0(fig_dir, "g123_totalRPA_filterCompares.svg"), height = 6, width = 12); p; dev.off()





