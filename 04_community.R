lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

##############################################################################################
# SECTION 4: SPATIO-TEMPORAL DISTRIBUTION OF EUKARYOTIC AND PROKARYOTIC COMMUNITIES

# This script processes and visualizes amplicon data for eukaryotic (18S rRNA) and prokaryotic 
# (16S rRNA) communities across three annual surveys conducted in 2016, 2017, and 2019.
# The analysis focuses on relative proportion abundances (RPA), diversity 
# patterns, and the spatial distribution of major microbial groups
##############################################################################################

# LOAD PACKAGES
library("skimr")
library("tidyverse")
library("viridis")
library("cowplot")
library("feather")
library("scales")
library("vegan")
library("ggpubr")

# Set themes
theme_set(theme_cowplot())
facet_theme <- theme(strip.background = element_rect(fill = "#44475A"),
                     strip.text = element_text(color = "white", face = "bold"))

# Set paths
dat_dir = "data_out/04_community/dataframes/"
fig_dir = "data_out/04_community/figures/"
tab_dir = "data_out/04_community/tables/"

############################################################################################
# LOAD IN DATAFRAMES
############################################################################################

meta <- read_csv("data_in/meta/g123_meta.csv") 
in_dir <- "data_out/02_qiime2_asv/dataframes/"

# Fetch 18s counts
getData <- function(file = NULL){
  df <- read_feather(file) %>% 
    mutate(featureID = substr(featureID, start = 1, stop = 10),
           group = kingdom) %>% 
    select(sampleID, featureID, group, counts, cruise, filter) %>% 
    mutate(domain = "eukaryote") %>% 
    filter(counts != 0)
}
g123_18s <- rbind(getData(file = paste0(in_dir, "g1_18s_master_0-200m.feather")),
                  getData(file = paste0(in_dir, "g2_18s_master_0-200m.feather")),
                  getData(file = paste0(in_dir, "g3_18s_master_0-200m.feather"))
                  )

# Fetch 18s taxonomies
getTaxonomy <- function(file = NULL) {
  read_feather(file) %>% 
    mutate(featureID = substr(featureID, start = 1, stop = 10)) %>% 
    select(featureID:species) %>% distinct %>% 
    mutate(group = case_when(superkingdom == "Archaeplastida" ~ superkingdom,
                             is.na(kingdom) ~ "Unknown Eukaryote",
                             TRUE ~ kingdom))
}
tax_18s <- rbind(
  getTaxonomy(file = paste0(in_dir, "g1_18s_master_0-200m.feather")),
  getTaxonomy(file = paste0(in_dir, "g2_18s_master_0-200m.feather")),
  getTaxonomy(file = paste0(in_dir, "g3_18s_master_0-200m.feather"))
) %>%
  distinct() %>%
  print()

# Combine counts with taxonomy
tax <- tax_18s %>% select(featureID, group)
euks <- g123_18s %>% 
  left_join(., tax, by = "featureID") %>% 
  select(-group.x) %>% 
  rename(group = group.y, rpa = counts)

#-------------------------------------------------------------------------------
# Fetch 16s counts
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

# Fetch 16s taxonomies
getTaxonomy16S <- function(file = NULL) {
  read_feather(file) %>%
    mutate(featureID = substr(featureID, start = 1, stop = 10)) %>%
    select(featureID:species) %>% distinct() %>%
    mutate(group = phylum, 
           group = case_when(
             str_detect(group, "SAR406") ~ "SAR406",
             str_detect(group, "SAR324") ~ "SAR324",
             grepl("Prochlorococcus", genus) ~ "Prochlorococcus",
             grepl("Synechococcus", genus) ~ "Synechococcus",
             str_detect(group, "Cyanobacteria") ~ "Other Cyanobacteria",
             domain == "Archaea" ~ "Archaea",
             TRUE ~ group))
}
tax_16s <- rbind(
  getTaxonomy16S(file = paste0(in_dir, "g1_16s_master_0-200m.feather")),
  getTaxonomy16S(file = paste0(in_dir, "g2_16s_master_0-200m.feather")),
  getTaxonomy16S(file = paste0(in_dir, "g3_16s_master_0-200m.feather"))
) %>%
  distinct() %>%
  print()

# Combine counts with taxonomy
tax <- tax_16s %>% select(featureID, group)
proks <- g123_16s %>% 
  left_join(., tax, by = "featureID") %>% 
  select(-group.x) %>% 
  rename(group = group.y, rpa = counts)

#-------------------------------------------------------------------------------
# IDENTIFYING MAJOR TAXONOMIC GROUPS, THROW THE LOWS IN THE 'OTHER CATEGORY'

recatGroups <- function(data = NULL, groupUp = NULL){
  group_totals <- data %>%
    group_by(group) %>%
    summarise(total_rpa = sum(rpa)) %>%
    arrange(total_rpa)
  
  total_rpa_sum <- sum(group_totals$total_rpa)
  threshold <- total_rpa_sum * groupUp 
  
  all <- unique(group_totals$group)
  lows <- group_totals %>% filter(total_rpa < threshold) %>% pull(group) %>% print
  highs <- setdiff(all, lows)
  
  data <- data %>% 
    mutate(group_unaltered = group,
           group = if_else(group %in% highs, group, "Other"))
  return(data)
}

skim(euks)
data <- recatGroups(data = euks, groupUp = 0.01)
skim(data)
euks <- data

skim(proks)
data <- recatGroups(data = proks, groupUp = 0.01)
skim(data)
proks <- data

############################################################################################
# REGIONAL STATISTICS
############################################################################################

# Eukaryotes
getPlotStats <- function(data = NULL, smpMeta = NULL, gradient = NULL){
  df <- data %>% mutate(sampleID = as.character(sampleID))
  x <- smpMeta %>% select(sampleID, region) %>% mutate(sampleID = as.character(sampleID))
  
  df <- df %>% left_join(., x, by = "sampleID") %>% 
    filter(cruise == gradient) %>% 
    group_by(sampleID, group, filter) %>% 
    reframe(region, rpa = sum(rpa)) %>% 
    distinct
  
  order <- df %>% 
    group_by(group) %>% 
    reframe(x = sum(rpa)) %>% 
    arrange(desc(x)) %>% 
    pull(group)
  
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
          legend.position = "top")
  
  return(list(plot = plot, stats = boxplot_stats))
}
p1 <- getPlotStats(data = euks, smpMeta = meta, gradient = "G1")
p2 <- getPlotStats(data = euks, smpMeta = meta, gradient = "G2")
p3 <- getPlotStats(data = euks, smpMeta = meta, gradient = "G3")


p <- plot_grid(p1$plot, p2$plot, p3$plot, nrow = 1); p
svg(paste0(fig_dir, "g123_euks_smpRegionals.svg"), height = 8, width = 16); p; dev.off()
write.csv(p1$stats, paste0(dat_dir, "g1_18s_boxPlotStats_metaComm.csv"), row.names = F)
write.csv(p2$stats, paste0(dat_dir, "g2_18s_boxPlotStats_metaComm.csv"), row.names = F)
write.csv(p3$stats, paste0(dat_dir, "g3_18s_boxPlotStats_metaComm.csv"), row.names = F)

# Prokaryotes
p1 <- getPlotStats(data = proks, smpMeta = meta, gradient = "G1")
p2 <- getPlotStats(data = proks, smpMeta = meta, gradient = "G2")
p3 <- getPlotStats(data = proks, smpMeta = meta, gradient = "G3")

p <- plot_grid(p1$plot, p2$plot, p3$plot, nrow = 1); p
svg(paste0(fig_dir, "g123_proks_smpRegionals.svg"), height = 8, width = 16); p; dev.off()
write.csv(p1$stats, paste0(dat_dir, "g1_16s_boxPlotStats_metaComm.csv"), row.names = F)
write.csv(p2$stats, paste0(dat_dir, "g2_16s_boxPlotStats_metaComm.csv"), row.names = F)
write.csv(p3$stats, paste0(dat_dir, "g3_16s_boxPlotStats_metaComm.csv"), row.names = F)

#-------------------------------------------------------------------------------
# BOX PLOT AVERAGES
getBoxplotStats <- function(data = NULL, smpMeta = NULL, gradient = NULL){
  df <- data %>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    left_join(smpMeta %>% select(sampleID, region) %>% mutate(sampleID = as.character(sampleID)), by = "sampleID") %>% 
    filter(cruise == gradient) %>% 
    group_by(sampleID, group, filter) %>% 
    reframe(region, rpa = sum(rpa)) %>% 
    distinct()
  
  df$region <- factor(df$region, levels = rev(c("NPSG", "STZ", "NTZ")))
  order <- df %>% 
    group_by(group) %>% 
    reframe(x = sum(rpa)) %>% 
    arrange(desc(x)) %>% 
    pull(group)
  df$group <- factor(df$group, levels = order)
  
  boxplot_stats1 <- df %>% 
    group_by(group, filter, region) %>% 
    summarise(
      min_rpa = min(rpa, na.rm = TRUE),
      median_rpa = median(rpa, na.rm = TRUE),
      max_rpa = max(rpa, na.rm = TRUE),
      mean_rpa = mean(rpa, na.rm = TRUE),
      .groups = 'drop'
    ) %>% 
    mutate(cruise = gradient)
  
  boxplot_stats2 <- df %>% 
    group_by(group, filter) %>% 
    summarise(
      min_rpa = min(rpa, na.rm = TRUE) *100,
      median_rpa = median(rpa, na.rm = TRUE)*100,
      max_rpa = max(rpa, na.rm = TRUE)*100,
      mean_rpa = mean(rpa, na.rm = TRUE)*100,
      .groups = 'drop'
    ) %>% 
    mutate(cruise = gradient)
  
  return(list(region = boxplot_stats1, cruise = boxplot_stats2))
}
getPlot <- function(data = NULL, gyre = NULL){
  data %>% 
    filter(region == gyre) %>% 
    ggplot(aes(x = group, y = as.character(filter), fill = mean_rpa)) +
    geom_tile(color = "white") + 
    geom_text(aes(label = round(mean_rpa*100, 1)), color = "white", size = 4) + 
    labs(
      fill = "Mean RPA",
      y = "Cruise",
      x = "Group"
    ) +
    facet_wrap(~ cruise) +
    facet_theme + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 14),
      legend.position = "right"
    ) +
    coord_flip() 
}

# Eukaryotes
p1 <- getBoxplotStats(data = euks, smpMeta = meta, gradient = "G1") 
p2 <- getBoxplotStats(data = euks, smpMeta = meta, gradient = "G2") 
p3 <- getBoxplotStats(data = euks, smpMeta = meta, gradient = "G3")

stats <- rbind(p1$cruise, p2$cruise, p3$cruise)
write.csv(stats, paste0(dat_dir, "g1_18s_boxPlotStats_Cruise.csv"), row.names = F)
stats <- rbind(p1$region, p2$region, p3$region)
write.csv(stats, paste0(dat_dir, "g1_18s_boxPlotStats_Region.csv"), row.names = F)

order <- stats %>% group_by(group) %>% 
  reframe(x = sum(mean_rpa)) %>% arrange(desc(x)) %>% 
  pull(group)

stats$group <- factor(stats$group, levels = rev(order))
p <- plot_grid(getPlot(data = stats, gyre = "NPSG"), 
               getPlot(data = stats, gyre = "STZ"), 
               getPlot(data = stats, gyre = "NTZ"), 
               ncol = 1); p
svg(paste0(fig_dir, "g123_euks_smpRegionals_average.svg"), height = 8, width = 10); p; dev.off()

# Prokaryotes
p1 <- getBoxplotStats(data = proks, smpMeta = meta, gradient = "G1") 
p2 <- getBoxplotStats(data = proks, smpMeta = meta, gradient = "G2") 
p3 <- getBoxplotStats(data = proks, smpMeta = meta, gradient = "G3")

stats <- rbind(p1$cruise, p2$cruise, p3$cruise)
write.csv(stats, paste0(dat_dir, "g1_16s_boxPlotStats_Cruise.csv"), row.names = F)
stats <- rbind(p1$region, p2$region, p3$region)
write.csv(stats, paste0(dat_dir, "g1_16s_boxPlotStats_Region.csv"), row.names = F)

order <- stats %>% group_by(group) %>% 
  reframe(x = sum(mean_rpa)) %>% arrange(desc(x)) %>% 
  pull(group)

stats$group <- factor(stats$group, levels = rev(order))
p <- plot_grid(getPlot(data = stats, gyre = "NPSG"), 
               getPlot(data = stats, gyre = "STZ"), 
               getPlot(data = stats, gyre = "NTZ"), 
               ncol = 1); p
svg(paste0(fig_dir, "g123_proks_smpRegionals_average.svg"), height = 8, width = 10); p; dev.off()

############################################################################################
# LATITUDINAL GRADIENT DISTRIBUTIONS
############################################################################################

data <- euks
# Identify low and high RPA groups
group_totals <- data %>%
  group_by(group) %>%
  summarise(total_rpa = sum(rpa)) %>%
  arrange(total_rpa)

# Acquire threshold and apply it
total_rpa_sum <- sum(group_totals$total_rpa)
threshold <- total_rpa_sum * 0.01 # groups contributing less than 1% of the total RPA

all <- unique(group_totals$group)
lows <- group_totals %>% filter(total_rpa < threshold) %>% pull(group) %>% print
highs <- setdiff(all, lows)

data <- data %>% 
  mutate(group_unaltered = group,
         group = if_else(group %in% highs, group, "Other"))
euks_archive <- data %>% select(-group_unaltered)

# Function to get sample order based on latitude
getOrder <- function(data, gradient, size) {
  data %>%
    filter(cruise == gradient & filter == size) %>%
    arrange(desc(latitude)) %>%
    pull(sampleID)
}
# Function to generate ggplot of total group RPA across years
getPlot <- function(data, size, gradient, smp_order, colors) {
  data %>%
    mutate(group = factor(group, levels = names(colors))) %>%
    filter(filter == size & cruise == gradient, sampleID %in% smp_order) %>%
    group_by(sampleID, cruise, filter, group) %>%
    summarise(total_rpa = sum(rpa)) %>%
    mutate(sampleID = factor(sampleID, levels = smp_order)) %>%
    ggplot(aes(x = total_rpa * 100, y = sampleID, fill = group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = colors)
}
# Function to merge plots for both filter sizes
getMerge <- function(df, cruise_type, colors) {
  p <- plot_grid(
    getPlot(df, size = "0.2", gradient = cruise_type, smp_order = getOrder(meta, cruise_type, "0.2"), colors = colors),
    getPlot(df, size = "3", gradient = cruise_type, smp_order = getOrder(meta, cruise_type, "3"), colors = colors),
    nrow = 1
  )
  return(p)
}

euk_colors <- c(
  "Archaeplastida" = "#49C5B1",
  "Alveolata" = "#DA4949",
  "Haptophyta" = "#FFD700",
  "Stramenopiles" = "#ADECF9",
  "Opisthokonta" = "#B19CD9", 
  "Rhizaria" = "#32CD32",
  "Other" = "#44475A",
  "Unknown Eukaryote" = "#FFB6C1"
)
# Generate spatial plots for eukaryote datasets
p1 <- getMerge(df = euks_archive, cruise_type = "G1", colors = euk_colors)
p2 <- getMerge(df = euks_archive, cruise_type = "G2", colors = euk_colors)
p3 <- getMerge(df = euks_archive, cruise_type = "G3", colors = euk_colors)

svg(paste0(fig_dir, "G1_euks_totalRPA.svg"), height = 8, width = 15); p1; dev.off()
svg(paste0(fig_dir, "G2_euks_totalRPA.svg"), height = 8, width = 15); p2; dev.off()
svg(paste0(fig_dir, "G3_euks_totalRPA.svg"), height = 8, width = 15); p3; dev.off()


data <- proks
# Identify low and high RPA groups
group_totals <- data %>%
  group_by(group) %>%
  summarise(total_rpa = sum(rpa)) %>%
  arrange(total_rpa)

# Acquire threshold and apply it
total_rpa_sum <- sum(group_totals$total_rpa)
threshold <- total_rpa_sum * 0.01 # groups contributing less than 1% of the total RPA

all <- unique(group_totals$group)
lows <- group_totals %>% filter(total_rpa < threshold) %>% pull(group) %>% print
highs <- setdiff(all, lows)

data <- data %>% 
  mutate(group_unaltered = group,
         group = if_else(group %in% highs, group, "Other"))
proks_archive <- data %>% select(-group_unaltered)

# Functions to generate ggplot of total group RPA across years
getPlot <- function(data, size, gradient, smp_order, colors){
  x_noChloroplast <- data %>% 
    filter(filter == size & cruise == gradient) %>%
    filter(sampleID %in% smp_order) %>% 
    select(-featureID, -filter, -cruise)
  
  x_Chloroplast <- data %>%
    group_by(sampleID) %>%
    summarize(missing_rpa = 1 - sum(rpa)) %>%
    filter(missing_rpa > 0) %>%
    mutate(group = "Chloroplast",
           rpa = missing_rpa,
           total_rpa = missing_rpa * 100,
           domain = "prokaryote") %>% 
    select(sampleID, rpa, domain, group)
  rbind(x_noChloroplast, x_Chloroplast) %>% 
    filter(sampleID %in% smp_order) %>% 
    mutate(group = factor(group, levels = names(colors))) %>% 
    group_by(sampleID, group) %>% 
    reframe(total_rpa = sum(rpa)) %>% 
    mutate(sampleID = factor(sampleID, levels = smp_order)) %>% 
    arrange(sampleID) %>% 
    
    ggplot(., aes(x = total_rpa * 100, y = sampleID, fill = group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = colors)
}

prok_colors <- c(
  "Chloroplast" = "#D9D9D9",
  "Other" = "#44475A",
  "Synechococcus" = "#FDB462",       
  "Prochlorococcus" = "#BC80BD",     
  "Proteobacteria" = "#49C5B1",
  "Archaea" = "#FFD700",    
  "SAR406" = "#8BE9FD",              
  "SAR324" = "#B3DE69",
  "Bacteroidota" = "#6272A4", 
  "Actinobacteriota" =   "#FF79C6",   
  "Planctomycetota" = "#FFFFB3", 
  "Verrucomicrobiota" = "#80B1D3"
)

# GENERATE SPATIAL PLOTS FOR PROKARYOTE DATASETS
p1 <- getMerge(df = proks_archive, cruise_type = "G1", colors = prok_colors)
p2 <- getMerge(df = proks_archive, cruise_type = "G2", colors = prok_colors)
p3 <- getMerge(df = proks_archive, cruise_type = "G3", colors = prok_colors)

svg(paste0(fig_dir, "G1_prok_totalRPA.svg"), height = 8, width = 15); p1; dev.off()
svg(paste0(fig_dir, "G2_prok_totalRPA.svg"), height = 8, width = 15); p2; dev.off()
svg(paste0(fig_dir, "G3_prok_totalRPA.svg"), height = 8, width = 15); p3; dev.off()

############################################################################################
# NON-METRIC MULTIDIMENSIONAL SCALING
############################################################################################

euks <- euks_archive
proks <- proks_archive

data = euks
gradient = "G3"
getNMDS <- function(data = NULL, gradient = NULL, plot_title = NULL){
  temp <- meta %>% mutate(sampleID = as.character(sampleID)) %>% 
    select(sampleID, depth)
  x <- data %>%
    filter(cruise == gradient) %>%
    group_by(sampleID, group, filter) %>%
    summarise(total_rpa = sum(rpa)) %>%
    pivot_wider(names_from = group, values_from = total_rpa, values_fill = list(total_rpa = 0)) %>% 
    left_join(., temp, by = "sampleID") %>% 
    relocate(depth, .before = filter)
  mat <- x[, -c(1:3)]
  set.seed(123)
  nmds <- metaMDS(mat, distance = "bray", k = 2, trymax = 100)
  nmds_sites <- as.data.frame(scores(nmds)$sites)
  nmds_sites$sampleID <- x$sampleID 
  
  data <- nmds_sites %>% 
    merge(., meta %>% select(sampleID, region, filter, depth), by = "sampleID") %>% 
    relocate(., c("region", "filter", "depth"), .after = sampleID)
  
  correlations <- as.data.frame(cor(x[,-(1:3)], nmds_sites[,1:2]))
  colnames(correlations) <- c("NMDS1", "NMDS2")
  stress_value <- round(nmds$stress, 4)
  
  ggplot(data, aes(x = NMDS1, y = NMDS2, color = region)) +
    geom_point(aes(shape = as.character(filter)), size = 4.5, color = "black") +
    geom_point(aes(shape = as.character(filter)), size = 3) +
    scale_color_manual(values = c("#FF79C6", "#8BE9FD", "#FFD700")) +
    geom_segment(data = correlations, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 arrow = arrow(length = unit(.35, "cm")), 
                 color = "black") +
    geom_text(data = correlations, aes(x = NMDS1, y = NMDS2, label = rownames(correlations)),
              color = "black", size = 3, hjust = "inward", vjust = "inward") + 
    labs(x = "NMDS1", y = "NMDS2", title = plot_title) + 
    annotate("text", x = Inf, y = -Inf, label = paste("Stress =", stress_value), 
             hjust = 1.1, vjust = -0.5, size = 5, color = "black") +
    theme(legend.position = "top")
}
# Generate NMDS plots for cruises
nmds_plot_list <- list()
for (cruise_type in c("G1", "G2", "G3")) {
  p <- plot_grid(
    getNMDS(data = proks, gradient = cruise_type, plot_title = "Prokaryote"),
    getNMDS(data = euks, gradient = cruise_type, plot_title = "Eukaryotes"),
    nrow = 1
  )
  nmds_plot_list[[cruise_type]] <- p
}

svg(file = paste0(fig_dir, "G1_NMDS.svg"), height = 5, width = 12)
print(nmds_plot_list[["G1"]])
dev.off()

svg(file = paste0(fig_dir, "G1_NMDS.svg"), height = 5, width = 12)
print(nmds_plot_list[["G2"]])
dev.off()

svg(file = paste0(fig_dir, "G1_NMDS.svg"), height = 5, width = 12)
print(nmds_plot_list[["G3"]])
dev.off()


############################################################################################
# ENVFIT ANALYSIS
############################################################################################

# Function for Factor Contributions
getEnvFit_factors <- function(data = NULL, meta = NULL, survey = NULL){
  df <- data %>% filter(cruise == survey)
  meta_g1 <- meta %>% filter(cruise == survey)
  smps1 <- unique(df$sampleID)
  smps2 <- unique(meta_g1$sampleID)
  smps <- intersect(smps1, smps2)
  meta_g1 <- meta_g1 %>% filter(sampleID %in% smps) %>% 
    mutate(sampleID = factor(sampleID, levels = smps)) %>% 
    arrange(sampleID)
  temp <- meta_g1 %>% select(sampleID, depth)
  df <- df %>% filter(sampleID %in% smps) %>%
    left_join(., temp, by = "sampleID") %>% 
    mutate(sampleID = factor(sampleID, levels = smps)) %>% 
    arrange(sampleID)
  com1 <- df %>% select(sampleID, featureID, rpa) %>%
    spread(featureID, rpa, fill = 0)
  com_matrix <- as.matrix(com1 %>% select(-sampleID))
  rownames(com_matrix) <- com1$sampleID
  bc <- vegdist(com_matrix, method = "bray")
  nmds_results <- metaMDS(bc, k = 2, trymax = 500)
  meta_g1$filter <- as.factor(meta_g1$filter)
  meta_g1$depth <- as.factor(meta_g1$depth)
  envfit_meta <- envfit(nmds_results, meta_g1[, c("filter", "region", "depth")], permutations = 999)
  p <- envfit_meta
  return(p)
}

data = 
# Function for Taxa Contributions
getEnvFit_taxa <- function(data = NULL, meta = NULL, survey = NULL){
  df <- data %>% filter(cruise == survey)
  meta_g1 <- meta %>% filter(cruise == survey)
  smps1 <- unique(df$sampleID)
  smps2 <- unique(meta_g1$sampleID)
  smps <- intersect(smps1, smps2)
  meta_g1 <- meta_g1 %>% filter(sampleID %in% smps) %>% 
    mutate(sampleID = factor(sampleID, levels = smps)) %>% 
    arrange(sampleID)
  df <- df %>% filter(sampleID %in% smps) %>% 
    mutate(sampleID = factor(sampleID, levels = smps)) %>% 
    arrange(sampleID)
  com1 <- df %>% select(sampleID, featureID, rpa) %>%
    spread(featureID, rpa, fill = 0)
  group_means <- com1 %>%
    pivot_longer(., cols = -sampleID, names_to = "featureID", values_to = "rpa") %>% distinct %>% 
    left_join(df %>% select(featureID, group), by = "featureID") %>% distinct %>% 
    group_by(sampleID, group) %>%
    summarise(mean_rpa = mean(rpa, na.rm = TRUE)) %>%
    spread(group, mean_rpa, fill = 0) %>%
    as.matrix()
  com_matrix <- as.matrix(com1 %>% select(-sampleID))
  rownames(com_matrix) <- com1$sampleID
  bc <- vegdist(com_matrix, method = "bray")
  nmds_results <- metaMDS(bc, k = 2, trymax = 100)
  envfit_group <- envfit(nmds_results, group_means, permutations = 999)
  p <- envfit_group
  return(p)
}
# Function to make a table
create_summary_table <- function(g1, g2, g3) {
  g1_summary <- data.frame(
    Survey = "G1",
    Factor = c("Filter", "Region", "Depth"),
    r2 = c(g1$fit$factors$r[1], g1$fit$factors$r[2], g1$fit$factors$r[3]),
    p_value = c(g1$fit$factors$pvals[1], g1$fit$factors$pvals[2], g1$fit$factors$pvals[3])
  )
  g2_summary <- data.frame(
    Survey = "G2",
    Factor = c("Filter", "Region", "Depth"),
    r2 = c(g2$fit$factors$r[1], g2$fit$factors$r[2], g2$fit$factors$r[3]),
    p_value = c(g2$fit$factors$pvals[1], g2$fit$factors$pvals[2], g1$fit$factors$pvals[3])
  )
  g3_summary <- data.frame(
    Survey = "G3",
    Factor = c("Filter", "Region", "Depth"),
    r2 = c(g3$fit$factors$r[1], g3$fit$factors$r[2], g3$fit$factors$r[3]),
    p_value = c(g3$fit$factors$pvals[1], g3$fit$factors$pvals[2], g1$fit$factors$pvals[3])
  )
  summary_table <- rbind(g1_summary, g2_summary, g3_summary)
  return(summary_table)
}
# Function to color the table
color_table_rows <- function(table_plot, table, survey_colors) {
  for (i in 1:nrow(table)) {
    table_plot <- table_cell_bg(
      table_plot, row = i + 1, column = 1, # row = i + 1 to account for the header row
      fill = survey_colors[table$Survey[i]],
      color = "white"
    )
  }
  return(table_plot)
}

#-------------------------------------------------------------------------------
# Eukaryotes
g1 <- getEnvFit_factors(data = euks, meta = meta, survey = "G1")
g2 <- getEnvFit_factors(data = euks, meta = meta, survey = "G2")
g3 <- getEnvFit_factors(data = euks, meta = meta, survey = "G3")
table <- data.frame(
  Survey = c("G1", "G1", "G1", "G2", "G2", "G2", "G3", "G3", "G3"),
  Factor = c("Filter", "Region", "Depth", "Filter", "Region", "Depth", "Filter", "Region", "Depth"),
  r2 = c(0.2823, 0.4678, 0.0712, 0.1684, 0.3832, 0.0800, 0.3723, 0.4042, 0.1324),
  p_value = c(0.001, 0.001, 0.843, 0.001, 0.001, 0.007, 0.001, 0.001, 0.001)
)
survey_colors <- c("G1" = "#50FA7B", "G2" = "#FFB86C", "G3" = "#BD93F9")
table_plot <- ggtexttable(table, rows = NULL, theme = ttheme("classic"))
table_plot <- color_table_rows(table_plot, table, survey_colors)
ggsave(paste0(tab_dir, "envFit_factors.png"), table_plot, width = 5, height = 3, dpi = 300)

g1 <- getEnvFit_taxa(data = euks, meta = meta, survey = "G1")
g2 <- getEnvFit_taxa(data = euks, meta = meta, survey = "G2")
g3 <- getEnvFit_taxa(data = euks, meta = meta, survey = "G3")
table <- data.frame(
  Survey = rep(c("G1", "G2", "G3"), each = 9),
  Taxa = rep(c("Alveolata", "Archaeplastida", "Haptophyta", "Opisthokonta", "Other", "Rhizaria", "Stramenopiles", "Unknown Eukaryote", "sampleID"), times = 3),
  r2 = c(0.7638, 0.3428, 0.3312, 0.8372, 0.2943, 0.1834, 0.5232, 0.2900, 0.1392, # G1
         0.2596, 0.7301, 0.2257, 0.5000, 0.0154, 0.3012, 0.0709, 0.2856, 0.0318, # G2
         0.7113, 0.5673, 0.3922, 0.7304, 0.1602, 0.1035, 0.6061, 0.1947, 0.4940), # G3
  p_value = c(0.001, 0.001, 0.001, 0.001, 0.004, 0.019, 0.001, 0.004, 0.083,   # G1
              0.001, 0.001, 0.001, 0.001, 0.646, 0.001, 0.100, 0.001, 0.380,  # G2
              0.001, 0.001, 0.001, 0.001, 0.001, 0.003, 0.001, 0.001, 0.001)  # G3
) %>% 
  filter(Taxa != "sampleID")

table_plot <- ggtexttable(table, rows = NULL, theme = ttheme("classic"))
table_plot <- color_table_rows(table_plot, table, survey_colors)
ggsave(paste0(tab_dir, "envFit_taxa.png"), table_plot, width = 5, height = 7, dpi = 300)

#-------------------------------------------------------------------------------
# Prokaryotes

proks <- proks_archive

g1 <- getEnvFit_factors(data = proks, meta = meta, survey = "G1")
g2 <- getEnvFit_factors(data = proks, meta = meta, survey = "G2")
g3 <- getEnvFit_factors(data = proks, meta = meta, survey = "G3")
table <- data.frame(
  Survey = c("G1", "G1", "G1", "G2", "G2", "G2", "G3", "G3", "G3"),
  Factor = c("Filter", "Region", "Depth", "Filter", "Region", "Depth", "Filter", "Region", "Depth"),
  r2 = c(0.1235, 0.5019, 0.4137, 0.0406, 0.7589, 0.0038, 0.0936, 0.4639, 0.1808),
  p_value = c(0.009, 0.001, 0.033, 0.098, 0.001, 0.767, 0.001, 0.001, 0.001)
)
table_plot <- ggtexttable(table, rows = NULL, theme = ttheme("classic"))
table_plot <- color_table_rows(table_plot, table, survey_colors)
ggsave(paste0(tab_dir, "envFit_factors_proks.png"), table_plot, width = 5, height = 3, dpi = 300)

g1 <- getEnvFit_taxa(data = proks, meta = meta, survey = "G1")
g2 <- getEnvFit_taxa(data = proks, meta = meta, survey = "G2")
g3 <- getEnvFit_taxa(data = proks, meta = meta, survey = "G3")
table <- data.frame(
  Survey = rep(c("G1", "G2", "G3"), each = 11),
  Taxa = rep(c("sampleID", "Actinobacteriota", "Archaea", "Bacteroidota", "Other", "Planctomycetota", 
               "Prochlorococcus", "Proteobacteria", "SAR406", "Synechococcus", "Verrucomicrobiota"), times = 3),
  r2 = c(0.1688, 0.4017, 0.0100, 0.3729, 0.3214, 0.3689, 0.5345, 0.4819, 0.4289, 0.4469, 0.6447,  
         0.0104, 0.4473, 0.3344, 0.4571, 0.0675, 0.3859, 0.5557, 0.2218, 0.2987, 0.4298, 0.1903,  
         0.5240, 0.0316, 0.1512, 0.2533, 0.3419, 0.1436, 0.2730, 0.1258, 0.0804, 0.2287, 0.1074), 
  p_value = c(0.059, 0.001, 0.840, 0.002, 0.005, 0.001, 0.001, 0.002, 0.001, 0.001, 0.001,
              0.737, 0.001, 0.001, 0.001, 0.121, 0.001, 0.001, 0.001, 0.001, 0.001, 0.003,
              0.001, 0.176, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.005, 0.001, 0.007)
) %>% 
  filter(Taxa != "sampleID")
table_plot <- ggtexttable(table, rows = NULL, theme = ttheme("classic"))
table_plot <- color_table_rows(table_plot, table, survey_colors)
ggsave(paste0(tab_dir, "envFit_taxa_proks.png"), table_plot, width = 5, height = 9, dpi = 300)


############################################################################################
# ALPHA DIVERSITY SPATIALS -- PHYTO AND NON-PHYTO
############################################################################################
library("phyloseq")

meta2 <- meta %>% select(sampleID, latitude, region, depth) %>% 
  mutate(sampleID = as.character(sampleID))

getDivs <- function(phyloPath, domain, gradient, sal, chl) {
  # Load your phyloseq object
  object <- readRDS(phyloPath)
  
  # Extract OTU table and convert to a matrix
  otu_table <- as(otu_table(object), "matrix")
  otu_table <- t(otu_table)
  
  # Calculate diversity metrics using vegan
  richness <- specnumber(otu_table) # Observed richness
  shannon <- diversity(otu_table, index = "shannon") # Shannon diversity
  invsimpson <- diversity(otu_table, index = "invsimpson") # Simpson diversity
  
  # Combine the diversity metrics into a dataframe
  richness_df <- data.frame(
    sampleID = rownames(otu_table),
    Observed = richness,
    Shannon = shannon,
    invSimpson = invsimpson
  )
  
  # Add sample information
  data <- meta2
  
  colnames(data); colnames(richness_df)
  
  # Combine the diversity metrics with sample metadata
  divs <- data %>% tibble %>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    left_join(., richness_df, by = "sampleID") %>% 
    select(sampleID, latitude, depth, Observed, Shannon, invSimpson) %>% 
    filter(depth <= 15) %>% 
    pivot_longer(cols = c("Observed", "Shannon", "invSimpson"),
                 names_to = "DivMetric",
                 values_to = "Value")
  
  # Summarize the diversity metrics by latitude
  divs_summary <- divs %>%
    mutate(latitude = round(latitude)) %>%
    group_by(latitude, DivMetric) %>%
    summarize(MeanValue = mean(Value, na.rm = TRUE),
              SE = sd(Value, na.rm = TRUE) / sqrt(n())) %>%
    ungroup()
  
  # Plot the diversity metrics
  p <- ggplot(divs_summary, aes(x = latitude, y = MeanValue, 
                                color = latitude, group = DivMetric)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(size = 1) +
    xlim(21, 43) +
    geom_errorbar(aes(ymin = MeanValue - SE, ymax = MeanValue + SE), width = 0.2) +
    facet_wrap(~DivMetric, ncol = 1, scales = "free_y") +
    geom_vline(xintercept = sal, size = 1.5) +
    geom_vline(xintercept = chl, size = 1.5, linetype = "dotted") +
    labs(title = paste0("Diversity Metrics - ", domain, " ", gradient),
         subtitle = "Solid line = Salinity Front \nDashed Line = Chl(a) Front",
         x = "Latitude (°N)",
         y = "Mean Diversity Value ± SE") +
    theme(legend.position = "right",
          strip.background = element_rect(fill = "black"),
          strip.text = element_text(color = "white", face = "bold", size = 15)
    )
  
  return(p)
}

# Cruise 3
p <- getDivs(phyloPath = "data_out/02_qiime2_asv/dataframes/g3_18s_phylo_r70k.RDS", 
             domain = "Eukaryotes", gradient = "G3", sal = 32.45, chl = 35)
svg(paste0(fig_dir, "g3_18s_divMetrics.svg"), height = 5, width = 6); p; dev.off() 
p <- getDivs(phyloPath = "data_out/02_qiime2_asv/dataframes/g3_16s_phylo_r50k.RDS", 
             domain = "Prokaryotes", gradient = "G3", sal = 32.45, chl = 35)
svg(paste0(fig_dir, "g3_16s_divMetrics.svg"), height = 5, width = 6); p; dev.off() 

############################################################################################
# END
############################################################################################
