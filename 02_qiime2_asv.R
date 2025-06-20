lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

##############################################################################################
# SECTION 2:
##############################################################################################

# Load in packages
library("tidyverse")
library("RColorBrewer")
library("vegan")
library("maps")
library("mapproj")
library("cowplot")
library("phyloseq")

# Set theme and colors
drac <- c("#50FA7B", "#FFB86C", "#BD93F9", "#FF79C6", "#FF5555", "#F1FA8C", 
          "#6272A4", "#8BE9FD", "#282A36", "#44475A", "#F8F8F2")

RSK_THEME <-  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 20, margin = margin(t = 0, r = 0, b = 10, l = 0)),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(face = "bold", size = 15, margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(face = "bold", size = 15, margin = margin(t = 0, r = 20, b = 0, l = 0)),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.position = "top", 
    legend.key.width = unit(1, 'cm'),
    legend.text = element_text(size = 20),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    panel.grid = element_blank()
  )

# Set paths
dat_dir = "data_out/02_qiime2_asv/dataframes/"
fig_dir = "data_out/02_qiime2_asv/figures/"
tab_dir = "data_out/02_qiime2_asv/tables/"

##############################################################################################
# Load in dataframes
##############################################################################################

# Fetch sample metaData
meta <- read.csv("data_in/meta/g123_meta.csv")

# LOAD AMPLICON COUNT TABLES
read_count_table <- function(cruise_name, count_table_path) {
  samps <- read_csv("data_in/meta/g123_meta.csv") %>% 
    filter(cruise == cruise_name) %>% 
    select(sampleID) %>% 
    pull()
  asv <- read.csv(count_table_path, sep = "\t")
  colnames(asv) <- str_remove(colnames(asv), pattern = "X")
  asv <- asv %>%
    pivot_longer(., cols = 2:(ncol(.)), names_to = "sampleID", values_to = "counts") %>%
    filter(sampleID %in% samps)
  return(asv)
}

## Prokaryotes
asv1.pro <- read_count_table("G1", "data_in/g1/g1_16s_countTable.csv")
asv2.pro <- read_count_table("G2", "data_in/g2/g2_16s_countTable.csv")
asv3.pro <- read_count_table("G3", "data_in/g3/g3_16s_countTable.csv")

## Eukaryotes
asv1.euk <- read_count_table("G1", "data_in/g1/g1_18s_countTable.csv")
asv2.euk <- read_count_table("G2", "data_in/g2/g2_18s_countTable.csv")
asv3.euk <- read_count_table("G3", "data_in/g3/g3_18s_countTable.csv")

# MERGE GRADIENTS TOGETHER
get_mergedDF <- function(g1 = NULL, g2 = NULL, g3 = NULL){
  asv <- rbind(g1, g2, g3) %>% pivot_wider(., names_from = "sampleID", values_from = "counts") %>% replace(is.na(.), 0)
  asv[,-1] <- sapply(asv[,-1], as.numeric)
  return(asv)
}

asv_pro <- get_mergedDF(g1 = asv1.pro, g2 = asv2.pro, g3 = asv3.pro)
asv_euk <- get_mergedDF(g1 = asv1.euk, g2 = asv2.euk, g3 = asv3.euk)

##############################################################################################
# Exploratory Analysis
##############################################################################################

## Histogram of READ DISTRIBUTIONS
get_histo <- function(data = NULL){
  data <- data.frame(sampleID = colnames(data[,-1]), reads = colSums(data[,-1]), meta)
  p <- ggplot(data, aes(x = reads, fill = cruise)) +
    geom_histogram(position="identity", bins = 50, alpha = 0.5) +
    scale_fill_manual(values = drac) +
    labs(fill = "Cruise", x = "Mapped ASVs", y = "Frequency") +
    RSK_THEME
  print(p); return(p)
}

p <- get_histo(data = asv_pro); p 
svg(paste0(fig_dir, "hist_16S_mappedASVs.svg"), height = 5, width = 10); p; dev.off()
p <- get_histo(data = asv_euk); p
svg(paste0(fig_dir, "hist_18S_mappedASVs.svg"), height = 5, width = 10); p; dev.off()

# RARECURVES
quickRareCurve <- function(x, step = 1, sample, 
                           xlab = "Sample Size", ylab = "Species",
                           label = TRUE, col, lty, max.cores = TRUE, nCores = 1, ...) {
  require(parallel)
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) {
    stop("function accepts only integers (counts)")}
  if (missing(col)) {
    col <- par("col")}
  if (missing(lty)) {
    lty <- par("lty")}
  tot <- rowSums(x) # calculates library sizes
  S <- specnumber(x) # calculates n species for each sample
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]} # removes any empty rows
  nr <- nrow(x) # number of samples
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  # parallel mclapply, set number of cores
  mc <- getOption("mc.cores", ifelse(max.cores, detectCores(), nCores))
  message(paste("Using", mc, "cores"))
  out <- mclapply(seq_len(nr), mc.cores = mc, function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])}
    drop(rarefy(x[i, ], n))})
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),
                                           y = z, xout = sample, 
                                           rule = 1)$y)
    abline(h = rare, lwd = 0.5)}
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)}
  if (label) {
    ordilabel(cbind(tot, S), labels = rownames(x), ...)}
  invisible(out)
}

get_rareCurve <- function(data = NULL, g1_cutoff = NULL, g2_cutoff = NULL, g3_cutoff = NULL){
  x <- quickRareCurve(t(data[,-1]), nCores = 32, step = 1000)
  
  rareGraph <- data.frame()
  for(i in 1:length(x)){
    id <- names(attr(x[[i]], "Subsample"))[length(names(attr(x[[i]], "Subsample")))]
    tmp <- data.frame(sampleID = id, SampleSize = attr(x[[i]], "Subsample"), Species = as.vector(x[[i]]))
    rareGraph <- rbind(rareGraph, tmp)
  }
  rareGraph$sampleID <- as.integer(as.character(rareGraph$sampleID)) 
  rareGraph <- merge(rareGraph, meta, by.all = "sampleID")
  
  p <- ggplot(rareGraph, aes(x = SampleSize, y = Species, group = sampleID, color = cruise)) +
    geom_line(alpha = 1.0) + labs(col = "Cruise Type") +
    scale_color_manual(values = drac) + 
    facet_wrap(~cruise, ncol = 1) +
    RSK_THEME +
    geom_vline(data = filter(rareGraph, cruise == "G1"), aes(xintercept = g1_cutoff), colour="black", size = 1) + 
    geom_vline(data = filter(rareGraph, cruise == "G2"), aes(xintercept = g2_cutoff), colour="black", size = 1) + 
    geom_vline(data = filter(rareGraph, cruise == "G3"), aes(xintercept = g3_cutoff), colour="black", size = 1) +
    theme(legend.position = "none", 
          panel.border = element_rect(color = drac[11], fill = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = drac[11]))
  print(p); return(p)
}

pp <- get_rareCurve(data = asv_pro, g1_cutoff = 90000, g2_cutoff = 20000, g3_cutoff = 50000)
pe <- get_rareCurve(data = asv_euk, g1_cutoff = 35000, g2_cutoff = 10000, g3_cutoff = 70000)

svg(paste0(fig_dir, "rCurve_16S.svg"), height = 8, width = 6); pp; dev.off()
svg(paste0(fig_dir, "rCurve_18S.svg"), height = 8, width = 6); pe; dev.off()

# SAMPLE DEPTH CUT-OFFS
generate_plot <- function(df, cruise_type, threshold, color) {
  df <- df %>%
    pivot_wider(., names_from = "sampleID", values_from = "counts") %>%
    replace(is.na(.), 0)
  meta_data <- meta %>% dplyr::filter(cruise == cruise_type)
  asvCount <- data.frame(sampleID = colnames(df[,-1]), reads = colSums(df[,-1]), meta_data)
  asvCount$fill <- "x"
  
  ggplot(asvCount, aes(x = reads, y = fill, color = color)) +
    geom_point(aes(size = 2)) + 
    geom_vline(xintercept = threshold, color = "black", size = 2, linetype = "solid") +
    labs(x = "Sample Size", y = "Dots are Samples") +
    scale_color_manual(values = color) +
    RSK_THEME +
    theme(axis.text.y = element_blank(),
          legend.key.width = unit(1, 'cm'),
          legend.text = element_text(size = 20),
          legend.position = "none")
}

pp <- plot_grid(
  generate_plot(asv1.pro, "G1", 90000, drac[1]), 
  generate_plot(asv2.pro, "G2", 20000, drac[2]), 
  generate_plot(asv3.pro, "G3", 50000, drac[3]), 
  ncol = 1
)
svg(paste0(fig_dir, "smpDepthCuts_16S.svg"), height = 8, width = 6); pp; dev.off()

pe <- plot_grid(
  generate_plot(asv1.euk, "G1", 35000, drac[1]), 
  generate_plot(asv2.euk, "G2", 10000, drac[2]), 
  generate_plot(asv3.euk, "G3", 70000, drac[3]), 
  ncol = 1
)
svg(paste0(fig_dir, "smpDepthCuts_18S.svg"), height = 8, width = 6); pe; dev.off()


# Make tables for reporting 
rarefy_thresholds <- data.frame(
  Cruise = c("2016_Gradients1", "2017_Gradients2", "2019_Gradients3"),
  'Prokaryotes_16S' = c(90000, 20000, 50000),
  'Eukaryotes_18S' = c(35000, 10000, 70000)
); library(flextable)

rarefy_table <- flextable(rarefy_thresholds) %>%
  theme_vanilla() %>%
  bg(j = 1, i = 1, bg = "#50FA7B") %>%  # Green for G1
  bg(j = 1, i = 2, bg = "#FFB86C") %>%  # Orange for G2
  bg(j = 1, i = 3, bg = "#BD93F9") %>%  # Purple for G3
  autofit()
rarefy_table

save_as_image(
  rarefy_table,
  path = paste0(tab_dir, "RarefyThresholds.png"),
  webshot = "webshot2",
  zoom = 2,
  expand = 1
)

############################################################################################
# GENERATE OBJECTS FOR PHYLOSEQ
############################################################################################

# GET SAMPLE IDS PER GRADIENT CRUISE
get_cruise_samples <- function(cruise_name) {
  metadata <- read_csv("data_in/meta/g123_meta.csv")
  samples <- metadata %>% filter(cruise == cruise_name) %>% pull(sampleID)
  return(samples)
}
g1samps <- get_cruise_samples("G1"); length(g1samps)
g2samps <- get_cruise_samples("G2"); length(g2samps)
g3samps <- get_cruise_samples("G3"); length(g3samps)

# ----------------------------------------------------------------------------------------
# PROKARYOTES 

process_16s_data <- function(cruise_dir, cruise_name, sample_filter) {
  # Load asvs
  counts <- read.csv(paste0(cruise_dir, "/g", cruise_name, "_16s_countTable.csv"), sep = "\t")
  colnames(counts) <- str_remove(colnames(counts), pattern = "X")
  rownames(counts) <- counts$featureID
  counts[,-1] <- sapply(counts[,-1], as.numeric)
  counts <- counts[,-1]
  ASV <- otu_table(as.matrix(counts), taxa_are_rows = TRUE)
  # Load taxonomy
  tax <- read.csv(paste0(cruise_dir, "/g", cruise_name, "_16s_taxonomy.csv"), sep = "\t") %>% select(featureID, taxa)
  rownames(tax) <- tax$featureID
  tax$taxa <- gsub(":plas", "", as.character(tax$taxa))
  tax$taxa <- gsub(" ", "", as.character(tax$taxa))
  tax$taxa <- gsub("s__", "", as.character(tax$taxa))
  df <- str_split_fixed(tax$taxa, ";", 7)
  df[df == ""] <- NA
  df2 <- cbind(tax, df)
  df2 <- df2 %>% 
    dplyr::select(-taxa) %>% 
    dplyr::rename("domain" = '1', 
                  "phylum" = '2',
                  "class" = '3', 
                  "order" = '4',
                  "family" = '5', 
                  "genus" = '6',
                  "species" = '7')
  TAX <- tax_table(as.matrix(df2))
  # Load meta
  meta <- read.csv("data_in/meta/g123_meta.csv") %>%
    filter(cruise == paste0("G", cruise_name))
  rownames(meta) <- meta$sampleID
  SAMPS <- sample_data(meta)
  # Load seqs
  seqs <- read.csv(paste0(cruise_dir, "/g", cruise_name, "_16s_dnaSeqs.csv"), sep = ",")
  rownames(seqs) <- seqs$featureID
  seqs <- seqs %>% select(-1)
  SEQS <- Biostrings::DNAStringSet(as.matrix(seqs))
  names(SEQS) <- taxa_names(TAX)
  # Create phyloseq object
  phy <- phyloseq(ASV, TAX, SAMPS)
  phy_merged <- merge_phyloseq(phy, SEQS)
  return(phy_merged)
}
phy_g1_16s <- process_16s_data(cruise_dir = "data_in/g1", cruise_name = "1", sample_filter = g1samps)
phy_g2_16s <- process_16s_data(cruise_dir = "data_in/g2", cruise_name = "2", sample_filter = g2samps)
phy_g3_16s <- process_16s_data(cruise_dir = "data_in/g3", cruise_name = "3", sample_filter = g3samps)

saveRDS(phy_g1_16s, paste0(dat_dir, "g1_16s_phylo.RDS"))
saveRDS(phy_g2_16s, paste0(dat_dir, "g2_16s_phylo.RDS"))
saveRDS(phy_g3_16s, paste0(dat_dir, "g3_16s_phylo.RDS"))

# ----------------------------------------------------------------------------------------
# EUKARYOTES

process_18s_data <- function(cruise_dir, cruise_name, sample_filter) {
  # Load asvs
  counts <- read.csv(paste0(cruise_dir, "/g", cruise_name, "_18s_countTable.csv"), sep = "\t")
  colnames(counts) <- str_remove(colnames(counts), pattern = "X")
  rownames(counts) <- counts$featureID
  counts[,-1] <- sapply(counts[,-1], as.numeric)
  counts <- counts[,-1]
  ASV <- otu_table(as.matrix(counts), taxa_are_rows = TRUE)
  # Load taxonomy
  tax <- read.csv(paste0(cruise_dir, "/g", cruise_name, "_18s_taxonomy_pr2v5.csv"), sep = "\t") 
  tax$Taxon <- gsub(":plas", "", as.character(tax$Taxon))
  df <- str_split_fixed(tax$Taxon, ";", 9)
  df[df == ""] <- NA
  df2 <- cbind(tax, df)
  df2 <- df2 %>% 
    dplyr::select(-Taxon, -Confidence) %>% 
    dplyr::rename("domain" = "1",
                  "superkingdom" = "2",
                  "kingdom" = "3", 
                  "phylum" = "4",
                  "class" = "5", 
                  "order" = "6",
                  "family" = "7", 
                  "genus" = "8",
                  "species" = "9",
                  "featureID" = "Feature.ID")
  df2$kingdom <- gsub("_X", "", df2$kingdom)
  df2$species <- gsub("[\\;,]", "", df2$species)
  df2 <- df2 %>% column_to_rownames(var = "featureID")
  TAX <- tax_table(as.matrix(df2))
  # Load meta
  meta <- read.csv("data_in/meta/g123_meta.csv") %>%
    filter(cruise == paste0("G", cruise_name))
  rownames(meta) <- meta$sampleID
  SAMPS <- sample_data(meta)
  # Load seqs
  seqs <- read.csv(paste0(cruise_dir, "/g", cruise_name, "_18s_dnaSeqs.csv"), sep = ",")
  rownames(seqs) <- seqs$featureID
  seqs <- seqs %>% tibble() %>% select(-1)
  SEQS <- Biostrings::DNAStringSet(as.matrix(seqs))
  names(SEQS) <- taxa_names(TAX)
  # Create phyloseq object
  phy <- phyloseq(ASV, TAX, SAMPS)
  data.frame(as.matrix(tax_table(phy)), stringsAsFactors = FALSE)
  phy_merged <- merge_phyloseq(phy, SEQS)
  return(phy_merged)
}

phy_g1_18s <- process_18s_data(cruise_dir = "data_in/g1", cruise_name = "1", sample_filter = g1samps)
phy_g2_18s <- process_18s_data(cruise_dir = "data_in/g2", cruise_name = "2", sample_filter = g2samps)
phy_g3_18s <- process_18s_data(cruise_dir = "data_in/g3", cruise_name = "3", sample_filter = g3samps)

saveRDS(phy_g1_18s, paste0(dat_dir, "g1_18s_phylo.RDS"))
saveRDS(phy_g2_18s, paste0(dat_dir, "g2_18s_phylo.RDS"))
saveRDS(phy_g3_18s, paste0(dat_dir, "g3_18s_phylo.RDS"))

############################################################################################
# RAREFY OBJECTS
############################################################################################

set.seed(1); phyR_g1_16s <- rarefy_even_depth(phy_g1_16s, sample.size = 90000, rngseed = 1, replace = F) 
set.seed(1); phyR_g2_16s <- rarefy_even_depth(phy_g2_16s, sample.size = 20000, rngseed = 1, replace = F) 
set.seed(1); phyR_g3_16s <- rarefy_even_depth(phy_g3_16s, sample.size = 50000, rngseed = 1, replace = F) 
set.seed(1); phyR_g1_18s <- rarefy_even_depth(phy_g1_18s, sample.size = 35000, rngseed = 1, replace = F) 
set.seed(1); phyR_g2_18s <- rarefy_even_depth(phy_g2_18s, sample.size = 10000, rngseed = 1, replace = F) 
set.seed(1); phyR_g3_18s <- rarefy_even_depth(phy_g3_18s, sample.size = 70000, rngseed = 1, replace = F) 

saveRDS(phyR_g1_16s, paste0(dat_dir, "g1_16s_phylo_r90k.RDS"))
saveRDS(phyR_g2_16s, paste0(dat_dir, "g2_16s_phylo_r20k.RDS"))
saveRDS(phyR_g3_16s, paste0(dat_dir, "g3_16s_phylo_r50k.RDS"))
        
saveRDS(phyR_g1_18s, paste0(dat_dir, "g1_18s_phylo_r35k.RDS"))
saveRDS(phyR_g2_18s, paste0(dat_dir, "g2_18s_phylo_r10k.RDS"))
saveRDS(phyR_g3_18s, paste0(dat_dir, "g3_18s_phylo_r70k.RDS"))


############################################################################################
# GENERATE MASTER DATAFRAMES
############################################################################################
library("feather")

# ----------------------------------------------------------------------------------------
# PROKARYOTES

process_data <- function(df, rarefyDepth) {
  # pull asvs
  count <- data.frame(as.matrix(otu_table(df)), stringsAsFactors = FALSE) %>% t()
  count <- count / rarefyDepth
  count <- t(count) %>% as.data.frame()
  count$featureID <- rownames(count)
  colnames(count) <- str_remove(colnames(count), pattern = "X")
  count <- count %>%
    pivot_longer(., cols = 1:ncol(count) - 1, names_to = "sampleID", values_to = "counts")
  # pull taxonomies
  tax <- data.frame(as.matrix(tax_table(df)), stringsAsFactors = FALSE)
  tax <- tax %>% filter(order != "Chloroplast" & family != "Mitochondria")
  # pull sequences
  seqs <- data.frame(refseq(df))
  seqs <- data.frame(featureID = rownames(seqs), seqs) %>% tibble() %>%
    rename("taxa" = 2)
  # pull in basic meta
  meta <- read.csv("data_in/meta/g123_meta.csv") %>% 
    select(sampleID, cruise, station, latitude, longitude, depth, filter)
  # Make master
  tc <- merge(tax, count, by = "featureID")
  tcs <- merge(tc, seqs, by = "featureID")
  tcsm <- tcs %>% select(featureID, sampleID, domain:species, taxa, counts) %>%
    merge(., meta, by = "sampleID") %>% as_tibble() %>% 
    # filter(depth <= 15) %>% 
    mutate(species = case_when(
      grepl("uncultured_", species) ~ "uncultured_sp", TRUE ~ species,
      grepl("unidentified_", species) ~ "unidentified_sp", TRUE ~ species,
      grepl("metagenome", species) ~ "unidentified_sp", TRUE ~ species,
      grepl("Candidatus_", species) ~ "Candidatus_sp", TRUE ~ species)
    )
  return(tcsm)
}
g1_16s <- process_data(df = phyR_g1_16s, rarefyDepth = 90000)
g2_16s <- process_data(df = phyR_g2_16s, rarefyDepth = 20000)
g3_16s <- process_data(df = phyR_g3_16s, rarefyDepth = 50000)

write_feather(g1_16s, paste0(dat_dir, "g1_16s_master_0-200m.feather"))
write_feather(g2_16s, paste0(dat_dir, "g2_16s_master_0-200m.feather"))
write_feather(g3_16s, paste0(dat_dir, "g3_16s_master_0-200m.feather"))

# ----------------------------------------------------------------------------------------
# EUKARYOTES

process_data <- function(df, rarefyDepth) {
  # pull asvs
  count <- data.frame(as.matrix(otu_table(df)), stringsAsFactors = FALSE) %>% t()
  count <- count / rarefyDepth
  count <- t(count) %>% as.data.frame()
  count$featureID <- rownames(count)
  colnames(count) <- str_remove(colnames(count), pattern = "X")
  count <- count %>%
    pivot_longer(., cols = 1:ncol(count) - 1, names_to = "sampleID", values_to = "counts")
  # pull taxonomies
  tax <- data.frame(as.matrix(tax_table(df)), stringsAsFactors = FALSE)
  tax <- tax %>% mutate(domain = gsub(" ", "", as.character(domain)),
                        superkingdom = gsub(" ", "", as.character(superkingdom)),
                        kingdom = gsub(" ", "", as.character(kingdom)),
                        phylum = gsub(" ", "", as.character(phylum)),
                        phylum = gsub("_X", "", as.character(phylum)),
                        class = gsub(" ", "", as.character(class)),
                        order = gsub(" ", "", as.character(order)),
                        family = gsub(" ", "", as.character(family)),
                        genus = gsub(" ", "", as.character(genus)),
                        species = gsub(" ", "", as.character(species)))
  tax <- tax %>% filter(domain != "Bacteria" & domain != "Archaea" & domain != "Unassigned") 
  tax$featureID <- rownames(tax)
  # pull sequences
  seqs <- data.frame(refseq(df))
  seqs <- data.frame(featureID = rownames(seqs), seqs) %>% tibble() %>% rename("taxa" = 2)
  # pull in basic meta
  meta <- read.csv("data_in/meta/g123_meta.csv") %>% 
    select(sampleID, cruise, station, latitude, longitude, depth, filter)
  # Make master
  tc <- merge(tax, count, by = "featureID")
  tcs <- merge(tc, seqs, by = "featureID")
  tcsm <- tcs %>% select(featureID, sampleID, domain:species, taxa, counts) %>%
    merge(., meta, by = "sampleID") %>% as_tibble() %>% 
    #filter(depth <= 15) 
    return(tcsm)
}
g1_18s <- process_data(df = phyR_g1_18s, rarefyDepth = 35000)
g2_18s <- process_data(df = phyR_g2_18s, rarefyDepth = 10000)
g3_18s <- process_data(df = phyR_g3_18s, rarefyDepth = 70000)

write_feather(g1_18s, paste0(dat_dir, "g1_18s_master_0-200m.feather"))
write_feather(g2_18s, paste0(dat_dir, "g2_18s_master_0-200m.feather"))
write_feather(g3_18s, paste0(dat_dir, "g3_18s_master_0-200m.feather"))

############################################################################################
# PRINCIPLE COORDINATE ANALYSIS
############################################################################################

# LOAD OBJECTS
phy_pro <- merge_phyloseq(readRDS(paste0(dat_dir, "g1_16s_phylo_r90k.RDS")), 
                          readRDS(paste0(dat_dir, "g2_16s_phylo_r20k.RDS")), 
                          readRDS(paste0(dat_dir, "g3_16s_phylo_r50k.RDS"))
                          )
phy_euk <- merge_phyloseq(readRDS(paste0(dat_dir, "g1_18s_phylo_r35k.RDS")), 
                          readRDS(paste0(dat_dir, "g2_18s_phylo_r10k.RDS")), 
                          readRDS(paste0(dat_dir, "g3_18s_phylo_r70k.RDS"))
                          )

# LOAD COLOR TEMPLATE FOR LATITUDE
gradient <- colorRampPalette(drac[c(1,6,2,3)])
scale_grad <- scale_color_gradientn(colors = gradient(100))

# FUNCTION 
plot_PCoA <- function(ps.ord, color_var, custom_scale = NULL) {
  temp <- data.frame(sampleID = rownames(ps.ord$vectors),
                     axis1 = ps.ord$vectors[,1],
                     axis2 = ps.ord$vectors[,2])
  temp2 <- meta %>% select(sampleID, latitude, cruise, filter)
  temp <- merge(temp, temp2)
  temp$filter <- as.character(temp$filter)
  
  p <- ggplot(temp) +
    geom_point(aes(x = axis1, y = axis2, 
                   color = get(color_var), 
                   size = factor(filter))) +
    scale_size_manual(values = c(3,7), guide = "none") +
    coord_flip() +
    labs(x = "PCoA.1", y = "PCoA.2", title = "PCoA Plot") +
    RSK_THEME +
    theme(legend.position = "right")
  if (!is.null(custom_scale)) {p <- p + custom_scale}
  print(p)
}

# PROKARYOTES
set.seed(1); ps.ord.pro <- ordinate(phy_pro, "PCoA", "bray")
p1 <- plot_PCoA(ps.ord.pro, 'latitude', scale_grad)
p2 <- plot_PCoA(ps.ord.pro, "cruise", scale_color_manual(values = drac))
p3 <- plot_PCoA(ps.ord.pro, "filter", scale_color_manual(values = c("#FF79C6", "#8BE9FD")))

svg(paste0(fig_dir, "pcoa_16S_latitudeR.svg"), height = 5, width = 8); p1; dev.off()  
svg(paste0(fig_dir, "pcoa_16S_cruiseR.svg"), height = 5, width = 8); p2; dev.off() 
svg(paste0(fig_dir, "pcoa_16S_filterR.svg"), height = 5, width = 8); p3; dev.off() 

# EUKARYOTES
set.seed(1); ps.ord.euk <- ordinate(phy_euk, "PCoA", "bray")
p1 <- plot_PCoA(ps.ord.euk, 'latitude', scale_grad)
p2 <- plot_PCoA(ps.ord.euk, "cruise", scale_color_manual(values = drac))
p3 <- plot_PCoA(ps.ord.euk, "filter", scale_color_manual(values = c("#FF79C6", "#8BE9FD")))

svg(paste0(fig_dir, "pcoa_18S_latitudeR.svg"), height = 5, width = 8); p1; dev.off()  
svg(paste0(fig_dir, "pcoa_18S_cruiseR.svg"), height = 5, width = 8); p2; dev.off() 
svg(paste0(fig_dir, "pcoa_18S_filterR.svg"), height = 5, width = 8); p3; dev.off()

############################################################################################
# PERMANOVA 
############################################################################################

getAdonis2 <- function(g1, g2, g3){
  g123 <- merge_phyloseq(g1, g2, g3)
  m <- as(sample_data(g123), "data.frame") 
  x <- adonis2(distance(g123, method = "bray") ~ year + latitude + filter, data = m) #, strata = m$filter)
  
  x <- data.frame(
    Term = rownames(x),
    R2 = x$R2 * 100,
    p_value = x$"Pr(>F)"
  )
  
  p <- ggplot(x[1:nrow(x)-1, ], aes(x = Term, y = R2, fill = Term)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label = sprintf("p = %.3f", p_value)), 
              vjust = -0.5, 
              position = position_dodge(width = 0.9)) +
    geom_text(aes(label = paste0(sprintf("%.3f", R2), "%")), 
              vjust = 1.5, 
              position = position_dodge(width = 0.9), 
              color = "black") +
    labs(title = "R2 and p-values from PERMANOVA Analysis",
         x = "Term",
         y = "R2 * 100") +
    theme_cowplot()
  return(p)
}
p <- getAdonis2(
  readRDS(paste0(dat_dir, "g1_16s_phylo_r90k.RDS")),
  readRDS(paste0(dat_dir, "g2_16s_phylo_r20k.RDS")),
  readRDS(paste0(dat_dir, "g3_16s_phylo_r50k.RDS"))
)
svg(paste0(fig_dir, "permanova_16s.svg"), height = 5, width = 8); p; dev.off()

p <- getAdonis2(
  readRDS(paste0(dat_dir, "g1_18s_phylo_r35k.RDS")),
  readRDS(paste0(dat_dir, "g2_18s_phylo_r10k.RDS")),
  readRDS(paste0(dat_dir, "g3_18s_phylo_r70k.RDS"))
)
svg(paste0(fig_dir, "permanova_18s.svg"), height = 5, width = 8); p; dev.off()
