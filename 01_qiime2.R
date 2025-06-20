lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

##############################################################################################
# SECTION 1: 
##############################################################################################

# Load in Packages
library("tidyverse")
library("ggpubr")
library("cowplot")

# Set theme and colors
drac <- c("#50FA7B", "#FFB86C", "#BD93F9", "#FF79C6", "#FF5555", 
          "#F1FA8C", "#6272A4", "#8BE9FD", "#282A36", "#44475A", 
          "#F8F8F2")

RSK_THEME <- theme_bw() +
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

# Set Directories
dat_dir = "data_out/01_qiime2/dataframes/"
fig_dir = "data_out/01_qiime2/figures/"
tab_dir = "data_out/01_qiime2/tables/"


###########################################################################################
# PREPROCESS, CLEAN QIIME2 OUTPUTS
###########################################################################################

# function to read and process CSV files
read_process_csv <- function(path, drop_first_row = FALSE) {
  data <- read_csv(path)
  if (drop_first_row)
    data <- data[-1, ]
  data <- data %>% mutate_if(is.character, as.numeric)
  return(data)
}

# function to rename columns
rename_columns <- function(data, forward_name, reverse_name, sample_name, type = "raw") {
  prefix <- ifelse(type == "raw", "raw", "trim")
  data <- data %>% mutate_if(is.character, as.numeric) %>% 
    rename(
      !!paste0(prefix, "_F") := {{ forward_name }},
      !!paste0(prefix, "_R") := {{ reverse_name }},
      "samp" := {{ sample_name }}
    )
  return(data)
}

# Gradients 1, 2016 Cruise
base_path = "data_in/g1/qiimeVisuals/"
## 16S files, Prokaryotes
d_g1_16s <- read_process_csv(paste0(base_path, "denoiseStats_16s.csv"), TRUE)
r_g1_16s <- read_csv(paste0(base_path, "rawSeqs_counts_16s.csv")) %>% 
  rename_columns("forward sequence count", "reverse sequence count", 'sample ID', type = "raw")
t_g1_16s <- read_csv(paste0(base_path, "rawSeqs_counts_trimmed_16s.csv")) %>% 
  rename_columns("forward sequence count", "reverse sequence count", 'sample ID', type = "trim")
## 16S files, Eukaruotes
d_g1_18s <- read_process_csv(paste0(base_path, "denoiseStats_16s.csv"), TRUE)
r_g1_18s <- read_csv(paste0(base_path, "rawSeqs_counts_18s.csv")) %>% 
  mutate_if(is.character, as.numeric) %>%
  rename(samp = 'sample ID', raw_F = 'forward sequence count')
t_g1_18s <- read_csv(paste0(base_path, "rawSeqs_counts_trimmed_18s.csv")) %>% 
  mutate_if(is.character, as.numeric) %>%
  rename(samp = 'sample ID', trim_F = 'forward sequence count')

# -------------------------------------------------------------------------------------------------
# Gradients 2, 2017 Cruise
base_path = "data_in/g2/qiimeVisuals/"
## 16S files, Prokaryotes
d_g2_16s <- read_process_csv(paste0(base_path, "denoiseStats_16s.csv"), TRUE)
r_g2_16s <- read_csv(paste0(base_path, "rawSeqs_counts_16s.csv")) %>% 
  rename_columns("forward sequence count", "reverse sequence count", 'sample ID', type = "raw")
t_g2_16s <- read_csv(paste0(base_path, "rawSeqs_counts_trimmed_16s.csv")) %>% 
  rename_columns("forward sequence count", "reverse sequence count", 'sample ID', type = "trim")
## 16S files, Eukaruotes
d_g2_18s <- read_process_csv(paste0(base_path, "denoiseStats_16s.csv"), TRUE)
r_g2_18s <- read_csv(paste0(base_path, "rawSeqs_counts_18s.csv")) %>% 
  mutate_if(is.character, as.numeric) %>%
  rename(samp = `sample ID`, raw_F = `forward sequence count`)
t_g2_18s <- read_csv(paste0(base_path, "rawSeqs_counts_trimmed_18s.csv")) %>% 
  mutate_if(is.character, as.numeric) %>%
  rename(samp = `sample ID`, trim_F = `forward sequence count`)

# -------------------------------------------------------------------------------------------------
# Gradients 3, 2019 Cruise
base_path = "data_in/g3/qiimeVisuals/"
## 16S files, Prokaryotes
d_g3_16s <- read_process_csv(paste0(base_path, "denoiseStats_16s.csv"), TRUE)
r_g3_16s <- read_csv(paste0(base_path, "rawSeqs_counts_16s.csv")) %>% 
  rename_columns("forward sequence count", "reverse sequence count", 'sample ID', type = "raw")
t_g3_16s <- read_csv(paste0(base_path, "rawSeqs_counts_trimmed_16s.csv")) %>% 
  rename_columns("forward sequence count", "reverse sequence count", 'sample ID', type = "trim")
# 16S files, Eukaruotes
d_g3_18s <- read_process_csv(paste0(base_path, "denoiseStats_16s.csv"), TRUE)
r_g3_18s <- read_csv(paste0(base_path, "rawSeqs_counts_18s.csv")) %>% 
  mutate_if(is.character, as.numeric) %>%
  rename(samp = `sample ID`, raw_F = `forward sequence count`)
t_g3_18s <- read_csv(paste0(base_path, "rawSeqs_counts_trimmed_18s.csv")) %>% 
  mutate_if(is.character, as.numeric) %>%
  rename(samp = `sample ID`, trim_F = `forward sequence count`)

###########################################################################################
# MAKE FIGURES OF QIIME2 RESULTS
###########################################################################################

## Helper function to generate the plot
generate_histo_error <- function(df = NULL, color = NULL) {
  df %>% ggplot(., aes(x = `percentage of input passed filter`)) +
    geom_histogram(binwidth = 0.5, alpha = 0.5, fill = color) +
    scale_x_continuous(breaks = round(seq(min(df$`percentage of input passed filter`), 
                                          max(df$`percentage of input passed filter`), 
                                          by = 0.5), 1)) +
    labs(title = "Trimmed Seqs that Passed Max Error Filter (%)", 
         x = "Percent Passed", 
         y = "Frequency Across Samples") + 
    RSK_THEME
}

generate_histo_merged <- function(df = NULL, color = NULL) {
  df %>% ggplot(., aes(x = `percentage of input merged`), fill = color) +
    geom_histogram(binwidth = 1, alpha = 0.5, fill = color) +
    scale_x_continuous(breaks = round(seq(min(df$`percentage of input merged`), 
                                          max(df$`percentage of input merged`), 
                                          by = 5), 5)) +
    labs(title = "Trimmed Seqs Succesfully Merged (%)", 
         x = "Percent Passed", 
         y = "Frequency Across Samples") + 
    RSK_THEME
}

## 16S - G1
p <- generate_histo_error(df = d_g1_16s, color = drac[1]); p
svg(paste0(fig_dir, "G1_16sDenoise_errorPass.svg"), height = 5, width = 10); p; dev.off()
p <- generate_histo_merged(df = d_g1_16s, color = drac[1])
svg(paste0(fig_dir, "G1_16sDenoise_merged.svg"), height = 5, width = 10); p; dev.off()
## 16S - G2
p <- generate_histo_error(df = d_g2_16s, color = drac[2])
svg(paste0(fig_dir, "G2_16sDenoise_errorPass.svg"), height = 5, width = 10); p; dev.off()
p <- generate_histo_merged(df = d_g2_16s, color = drac[2])
svg(paste0(fig_dir, "G2_16sDenoise_merged.svg"), height = 5, width = 10); p; dev.off()
## 16S - G3
p <- generate_histo_error(df = d_g3_16s, color = drac[3])
svg(paste0(fig_dir, "G3_16sDenoise_errorPass.svg"), height = 5, width = 10); p; dev.off()
p <- generate_histo_merged(df = d_g3_16s, color = drac[3])
svg(paste0(fig_dir, "G3_16sDenoise_merged.svg"), height = 5, width = 10); p; dev.off()


## 18S - G1
p <- generate_histo_error(df = d_g1_18s, color = drac[1])
svg(paste0(fig_dir, "G1_18sDenoise_errorPass.svg"), height = 5, width = 10); p; dev.off()
p <- generate_histo_merged(df = d_g1_18s, color = drac[1])
svg(paste0(fig_dir, "G1_18sDenoise_merged.svg"), height = 5, width = 10); p; dev.off()
## 18S - G2
p <- generate_histo_error(df = d_g2_18s, color = drac[2])
svg(paste0(fig_dir, "G2_18sDenoise_errorPass.svg"), height = 5, width = 10); p; dev.off()
p <- generate_histo_merged(df = d_g2_18s, color = drac[2])
svg(paste0(fig_dir, "G2_18sDenoise_merged.svg"), height = 5, width = 10); p; dev.off()
## 18S - G3
p <- generate_histo_error(df = d_g3_18s, color = drac[3])
svg(paste0(fig_dir, "G3_18sDenoise_errorPass.svg"), height = 5, width = 10); p; dev.off()
p <- generate_histo_merged(df = d_g3_18s, color = drac[3])
svg(paste0(fig_dir, "G3_18sDenoise_merged.svg"), height = 5, width = 10); p; dev.off()

# --------------------------------------------------------------------------------------
# Function to plot sequence count loss after trimming
plot_sequence_counts <- function(data, title) {
  levels <- merge(data$r, data$t, by.all = "samp") %>%
    select(samp, raw_F, trim_F) %>%
    arrange(desc(raw_F)) %>%
    pull(samp) %>% unique

merge(data$r, data$t, by.all = "samp") %>%
    select(samp, raw_F, trim_F) %>%
    pivot_longer(., cols = 2:(ncol(.)), names_to = "type", values_to = "value") %>%
    ggplot(., aes(y = value, x = factor(samp, level = levels), fill = type)) +
    geom_col(position = "dodge") +
    labs(title = title, x = "Samples", y = "Number of Sequences") +
    RSK_THEME + theme(axis.text.x = element_blank())
}

# Prokaryote 16S
x <- plot_sequence_counts(data = list(r = r_g1_16s, t = t_g1_16s), title = "Sequence Counts Lost after Trimming (G1 16S)")
svg(paste0(fig_dir, "G1_16sCounts_changeAfterTrim.svg"), height = 5, width = 20); x; dev.off()
x <- plot_sequence_counts(data = list(r = r_g2_16s, t = t_g2_16s), title = "Sequence Counts Lost after Trimming (G2 16S)")
svg(paste0(fig_dir, "G2_16sCounts_changeAfterTrim.svg"), height = 5, width = 20); x; dev.off()
x <- plot_sequence_counts(data = list(r = r_g3_16s, t = t_g3_16s), title = "Sequence Counts Lost after Trimming (G3 16S)")
svg(paste0(fig_dir, "G3_16sCounts_changeAfterTrim.svg"), height = 5, width = 20); x; dev.off()

# Eukaryote 18S
x <- plot_sequence_counts(data = list(r = r_g1_18s, t = t_g1_18s), title = "Sequence Counts Lost after Trimming (G1 18S)")
svg(paste0(fig_dir, "G1_18sCounts_changeAfterTrim.svg"), height = 5, width = 20); x; dev.off()
x <- plot_sequence_counts(data = list(r = r_g2_18s, t = t_g2_18s), title = "Sequence Counts Lost after Trimming (G2 18S)")
svg(paste0(fig_dir, "G2_18sCounts_changeAfterTrim.svg"), height = 5, width = 20); x; dev.off()
x <- plot_sequence_counts(data = list(r = r_g3_18s, t = t_g3_18s), title = "Sequence Counts Lost after Trimming (G3 18S)")
svg(paste0(fig_dir, "G3_18sCounts_changeAfterTrim.svg"), height = 5, width = 20); x; dev.off()

