# Clean workspace and unload loaded packages
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

##############################################################################################
# SECTION 3: ENVIRONMENTAL CONDITIONS

# This script processes and visualizes physiochemical sample data from three annual surveys 
# conducted in 2016, 2017, and 2019 in the North Pacific Ocean.
##############################################################################################

# LOAD PACKAGES
library("tidyverse")
library("cowplot")
library("lubridate")
library("ggpubr")
library("scales")
library("gt")
library("ggridges")

# Set themes
theme_set(theme_cowplot())
facet_theme <- theme(strip.background = element_rect(fill = "#44475A"),
                     strip.text = element_text(color = "white", face = "bold"))

# Set paths
dat_dir = "data_out/03_physioChems/dataframes/"
fig_dir = "data_out/03_physioChems/figures/"
tab_dir = "data_out/03_physioChems/tables/"

##############################################################################################
# ACQUIRE PHYSIOCHEM VARIABLES OF INTEREST
##############################################################################################

# Load in meta 
smpMeta <- read_csv("data_in/meta/g123_meta.csv")

# Function to handle time observations
fetchTime_format <- function(df) {
  df %>%
    mutate(
      timeDate = ymd_hms(time),
      year = year(timeDate),
      month = month(timeDate),
      day = day(timeDate),
      time = format(timeDate, "%H:%M:%S"),
      dayStatus = case_when(
        hour(timeDate) >= 3 & hour(timeDate) < 15 ~ "day",
        TRUE ~ "night")
    )
}
# Function to process dataframes after load-in
processUW_data <- function(file_path, file_name, cols_to_select) {
  file <- read.csv(paste0(file_path, file_name))  %>%
    select(all_of(cols_to_select)) %>%
    group_by(time, lat) %>%
    reframe(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    fetchTime_format() %>%
    mutate(timeID = paste0(year, "_", month, "_", day, "_", dayStatus)) %>%
    distinct() %>%
    drop_na() %>%
    arrange(lat)
}

#-------------------------------------------------------------------------------
# Load Productivity datasets
path <- "data_in/cmap/productivity/"
columns <- c("time", "lat", "temp", "sal", "O2_Ar_sat", "NCP")
NCPRdata <- list(
  g1 = processUW_data(path, "KOK1606_Gradients1_Surface_O2Ar_NCP.csv", columns),
  g2 = processUW_data(path, "MGL1704_Gradients2_Surface_O2Ar_NCP.csv", columns),
  g3 = processUW_data(path, "KM1906_Gradients3_Surface_O2Ar_NCP.csv", columns)
)

# Load Particulate datasets
path <- "data_in/cmap/particulates/"
columns <- c("time", "lat", "pc", "pn")
PARTdata <- list(
  g1 = processUW_data(path, "Gradients1_KOK1606_PPPCPN_UW.csv", columns),
  g2 = processUW_data(path, "Gradients2_MGL1704_PPPCPN_UW.csv", columns),
  g3 = processUW_data(path, "Gradients3_KM1906_PCPN_UW.csv", columns)
)

#-------------------------------------------------------------------------------
# Convert sample time stamps to 
# Hawaiian Standard Time (HST) to Coordinated Universal Time (UTC) 
fetchTime_zone <- function(df, timeZone = NULL) {
  df %>%
    mutate(
      time = ymd_hms(time_HST),
      time = if (!is.null(timeZone)) {
        force_tz(time, tz = timeZone) + hours(10) # HST to UTC
      } else {
        time
      },
      year = year(time),
      month = month(time),
      day = day(time),
      time = format(time, "%H:%M:%S")
    )
}

# Load ASV time stamp information and apply fetchTime_zone()
samp_g12_times <- read.csv("data_in/meta/g12_meta_times.csv") %>%
  select(Sample_ID, Datetime_ISO8601) %>%
  mutate(
    time = substr(Datetime_ISO8601, 1, nchar(Datetime_ISO8601) - 6),
    time_HST = ymd_hms(time)
  ) %>%
  drop_na() %>%
  rename(sampleID = Sample_ID) %>%
  fetchTime_zone(timeZone = "HST") %>%
  select(-Datetime_ISO8601)

samp_g3_times <- read.csv("data_in/meta/g3_meta_times.csv", sep = "\t") %>%
  select(sampleID, date, time) %>%
  mutate(
    time = ifelse(nchar(time) == 4, paste0("0", time), time),
    time = paste0(time, ":00"),
    time_HST = mdy_hms(paste(date, time))
  ) %>%
  select(sampleID, time_HST) %>%
  fetchTime_zone(timeZone = "HST") %>%
  drop_na()

# Combine sample times for all 3 cruises into one dataframe and save
samp_times <- rbind(samp_g12_times, samp_g3_times) %>%
  mutate(
    time_UTC = ymd_hms(paste(year, month, day, time)),
    dayStatus = case_when(
      hour(time_HST) >= 3 & hour(time_HST) < 15 ~ "day",
      TRUE ~ "night"
    )
  )
write.csv(samp_times, paste0(dat_dir, "g123_smpTimes.csv"))

##############################################################################################
# PERFORM SAMPLE TO FLOWDATA SPACE/TIME MATCHING
##############################################################################################
library(tidyverse)
# This function will find variable observations that align the closest with amplicon collection times.
# If collection times are >12 hours or latitudes >0.5 degrees, they will not be considered in later analyses.

getMatches <- function(mydata = NULL, col_selection = NULL) {
  
  samp_times <- read.csv(paste0(dat_dir, "g123_smpTimes.csv"))
  dim(samp_times)
  
  # Load and prepare sample metadata
  g123 <- smpMeta %>% tibble() %>%
    select(-year) %>%
    merge(., samp_times, by = "sampleID") %>%
    select(-station, -longitude, -time_HST) %>%
    mutate(timeID = paste0(year, "_", month, "_", day), hr_min_sec = time) %>%
    select(-time, -year, -month, -day, -time_UTC)
  
  # Prepare environmental data
  part <- mydata %>%
    mutate(timeID = paste0(year, "_", month, "_", day), hr_min_sec = time) %>%
    select(-time, -timeDate, -year, -month, -day, -dayStatus) %>%
    relocate(timeID:hr_min_sec, .before = "lat")
  
  # Perform the match
  find_best_match <- function(sample, part_df) {
    if (!is.data.frame(sample) || !"latitude" %in% names(sample)) {
      stop("The 'sample' is not a dataframe or 'latitude' column is missing.")
    }
    # Filter part_df for matching timeID
    matched_part <- filter(part_df, timeID == sample$timeID)
    # If no matches, return NA
    if (nrow(matched_part) == 0) {
      return(data.frame(matrix(NA, ncol = ncol(part_df), nrow = 1, dimnames = list(NULL, names(part_df)))))
    }
    # Calculate the closest latitude
    lat_diffs <- abs(sample$latitude - matched_part$lat)
    closest_lat_index <- which.min(lat_diffs)
    closest_lat_row <- matched_part[closest_lat_index, ]
    # Calculate the closest time
    sample_time <- ymd_hms(paste(sample$timeID, sample$hr_min_sec))
    part_times <- ymd_hms(paste(closest_lat_row$timeID, closest_lat_row$hr_min_sec))
    time_diffs <- abs(difftime(sample_time, part_times, units = "mins"))
    closest_time_index <- which.min(time_diffs)
    # Return result dataframe
    return(closest_lat_row[closest_time_index, ])
  }
  
  g123_with_best_match <- g123 %>%
    rowwise() %>%
    mutate(best_match = list(find_best_match(cur_data(), part))) %>% 
    unnest(best_match, names_sep = "_") %>%
    rename(
      dateID_A = timeID, time_A = hr_min_sec,
      dateID_P = best_match_timeID, time_P = best_match_hr_min_sec
    ) %>%
    rename_with(~ gsub("best_match_", "", .x), starts_with("best_match_")) %>%
    mutate(
      datetime_A = as.POSIXct(paste(dateID_A, time_A), format = "%Y_%m_%d %H:%M:%S"),
      datetime_P = as.POSIXct(paste(dateID_P, time_P), format = "%Y_%m_%d %H:%M:%S"),
      dif_t = as.numeric(abs(difftime(datetime_A, datetime_P, units = "mins")) / 60),
      dif_l = abs(latitude - lat)
    )
  
  # Apply time (within 12hr window) and latitude (within 0.5 degree area) cutoffs
  passed <- g123_with_best_match %>%
    select(-latitude, -depth, -filter, -X, -dayStatus) %>% 
    filter(dif_t <= 12 & dif_l <= 0.5) %>%
    relocate(c(datetime_A, datetime_P, dif_t, dif_l), .before = lat) %>%
    mutate(status = "passed",
           failure_reason = "None")
  
  failed <- g123_with_best_match %>%
    select(-latitude, -depth, -filter, -X, -dayStatus) %>%
    filter(dif_t > 12 | dif_l > 0.5 | is.na(dif_t) | is.na(dif_l)) %>%
    relocate(c(datetime_A, datetime_P, dif_t, dif_l), .before = lat) %>%
    mutate(
      status = "failed",
      failure_reason = case_when(
        dif_t > 12 ~ "Time_Fail",
        dif_l > 0.5 ~ "Lat_Fail",
        is.na(dif_t) | is.na(dif_l) ~ "No_Match",
        TRUE ~ "Unknown Fail"  
      )
    )
  
  # Prepare result dataframes
  result <- passed %>%
    select(-status, -failure_reason) %>%
    select(sampleID, last_col():(last_col() - col_selection))
  checks <- bind_rows(passed, failed) %>%
    select(sampleID, cruise, datetime_A, lat, datetime_P, dif_t, dif_l, status, failure_reason)
  # Combine them to print
  list(result = result, checks = checks)
}

# Productivity Matches
dfs <- rbind(NCPRdata$g1, NCPRdata$g2, NCPRdata$g3)
ndf <- getMatches(mydata = dfs, col_selection = 3)
# Particulate Matches
dfs <- rbind(PARTdata$g1, PARTdata$g2, PARTdata$g3)
pdf <- getMatches(mydata = dfs, col_selection = 1)

# Make an Environmental Meta df for amplicon samples and save
final <- smpMeta %>% select(-year, -longitude, -station) %>% 
  merge(., pdf$result %>% select(sampleID, pn, pc), all.x = T) %>% 
  merge(., ndf$result %>% select(sampleID, NCP, O2_Ar_sat, temp, sal), all.x = T) %>% 
  mutate(across(where(is.numeric), ~replace(., is.nan(.), NA))) 
write.csv(final, paste0(dat_dir, "g123_meta_physioChems.csv"), row.names = F)

# Make an Environmental Meta dfs containing pass/fail information for each data type
write.csv(pdf$checks, paste0(dat_dir, "match_particulates.csv"), row.names = F)
write.csv(ndf$checks, paste0(dat_dir, "match_productivity.csv"), row.names = F)


times <- smpMeta %>% 
  select(sampleID, month, day, time) %>% 
  rename(UTC_time = time) # Universal time

final <- smpMeta %>% 
  left_join(., times, by = "sampleID") %>% 
  relocate(year, .before = month) %>% 
  left_join(., pdf$result, by = "sampleID") %>% 
  left_join(., ndf$result, by = "sampleID") %>% 
  rename(PON = pn, POC = pc) %>% 
  select(-O2_Ar_sat, -sal, -temp) %>% 
  arrange(cruise, latitude, depth, filter) %>% 
  mutate(across(where(is.numeric), ~replace(., is.nan(.), NA))) %>% 
  tibble %>% distinct

write.csv(final, paste0(dat_dir, "g123_meta_final.csv"), row.names = F)

x <- final %>%
  filter(!is.na(UTC_time)) %>% 
  mutate(depth = case_when(station == "underway" ~ "Underway",
                           station == "tow-fish" ~ "Underway",
                           T ~ as.character(depth))) %>% 
  mutate(collectionTime = paste(year, month, day, UTC_time, sep = "_")) %>% 
  group_by(cruise, collectionTime, depth, latitude, filter) %>%
  summarise(n_reps = n(), .groups = "drop") %>% 
  arrange(cruise)

write.csv(x, paste0(dat_dir, "g123_meta_sampleSites.csv"), row.names = F)

#-------------------------------------------------------------------------------
# Visualize pass and fail proportions

samps <- final %>% pull(sampleID)

prod = pdf$checks
part = ndf$checks

# Proportion of pass/fail samples
getFraction <- function(data = NULL, plot_title = NULL){
  data %>%
    filter(sampleID %in% samps) %>% 
    group_by(cruise, status) %>%
    summarize(count = n()) %>%
    mutate(percentage = count / sum(count)) %>% 
    
    ggplot(., aes(x = "", y = percentage, fill = status)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    facet_wrap(~cruise) +
    geom_text(aes(label = paste(count, "(", round(100 * percentage, 1), "%)", sep = "")), 
              position = position_stack(vjust = 0.5)) +
    labs(fill = "Status", title = plot_title) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 20),
          plot.title = element_text(size = 30))
}
plotFailureReasons <- function(data = NULL, plot_title = NULL){
  data %>%
    filter(status == "failed") %>%
    group_by(failure_reason) %>%
    summarize(count = n()) %>%
    mutate(percentage = count / sum(count)) %>%
    
    ggplot(., aes(x = failure_reason, y = percentage, fill = failure_reason)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = paste(count, "(", round(100 * percentage, 1), "%)", sep = "")), 
              position = position_stack(vjust = 0.5), size = 5) +
    labs(fill = "Failure Reason", title = plot_title, x = "Failure Reason", y = "Percentage") +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 15, face = "bold"))
}
# Make plot and save
p <- plot_grid(
  getFraction(data = prod, plot_title = "Proportion of Matches - Productivity"),
  getFraction(data = part, plot_title = "Proportion of Matches - Particulates"),
  plotFailureReasons(prod, "Failure Reasons - Productivity Matches"),
  plotFailureReasons(part, "Failure Reasons - Particulate Matches"),
  ncol = 2
)
svg(paste0(fig_dir, "g123_smp0-200_matchResults.svg"), height = 12, width = 20); p; dev.off()

##############################################################################################
# VISUALIZE THE CHLOROPHYLL AND ISOHALINE SALINITY FRONTS
##############################################################################################

## Pull CTD salinity and CTD Chlorophyll data
processWaterCol <- function(file_path, gradients) {
  read.csv(file_path) %>%
    as_tibble() %>%
    filter(depth <= 15) %>%
    group_by(lat, depth) %>%
    summarize(
      temp = mean(CTD_Temperature, na.rm = TRUE),
      sal = mean(CTD_Salinity, na.rm = TRUE),
      chl = mean(CTD_Chloropigment, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    distinct() %>%
    arrange(lat) %>%
    mutate(cruise = gradients)
}
gradient_files <- list(
  G1 = "KOK1606_Gradients1.csv",
  G2 = "MGL1704_Gradients2.csv",
  G3 = "KM1906_Gradients3.csv"
)
water_data <- lapply(names(gradient_files), function(grad) {
  processWaterCol(file_path = paste0("data_in/cmap/waterColumn/", gradient_files[[grad]]), gradients = grad)
})

getPlot <- function(data = NULL, sal_front = NULL, chl_front = NULL){
  scale_factor <- (max(data$chl) - min(data$chl)) / (max(data$sal) - min(data$sal))
  data$scaled_sal <- (data$sal - min(data$sal)) * scale_factor + min(data$chl)
  
  data %>% 
    pivot_longer(., cols = c(chl, scaled_sal), names_to = "measurement", values_to = "value") %>% 
    ggplot(., aes(x = lat, y = value, color = "measurement")) +
    geom_point(aes(shape = measurement, color = measurement), size = 3, alpha = .4) +
    geom_smooth(aes(linetype = measurement), size = 2, method = "loess", se = FALSE, color = "black") +
    ylim(min(data$chl), max(data$chl)) +
    xlim(19, 43) +
    theme_minimal() +
    scale_color_manual(values = c("chl" = "#2b6439", "scaled_sal" = "#ad5700")) +
    scale_shape_manual(values = c("chl" = 1, "scaled_sal" = 1)) +
    scale_linetype_manual(values = c("chl" = "11", "scaled_sal" = "32")) + 
    geom_vline(aes(xintercept = chl_front), color = "#50FA7B", size = 2) +
    geom_vline(aes(xintercept = sal_front), color = "#FF7F00", size = 2) +
    coord_flip() +
    labs(x = "Latitude (°N)") +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold", size = 25, margin = margin(t = 0, r = 20, b = 0, l = 0)),
    )
}
p1 <- getPlot(data = water_data[[1]], sal_front = 32.15, chl_front = 33)
p2 <- getPlot(data = water_data[[2]], sal_front = 32.5, chl_front = 36.2)
p3 <- getPlot(data = water_data[[3]], sal_front = 32.45, chl_front = 35)

p <- plot_grid(p1, p2, p3, ncol = 1) %>% print
svg(paste0(fig_dir, "g123_fronts.svg"), height = 15, width = 4.5); p; dev.off()

############################################################################################################################################################################################
# VISUALIZE PHYSIOCHEM SPATIALS
############################################################################################################################################################################################

# Load metadata and filter based on depth
meta <- read.csv(paste0(dat_dir, "g123_meta_physioChems.csv")) %>%
  as_tibble() %>%
  filter(depth <= 15) %>%
  rename(lat = latitude) %>%
  mutate(region = case_when(lat < 30 ~ "NSG_STZ", lat > 30 ~ "NTZ", TRUE ~ "NA"))

# Define salinity and chlorophyll fronts for different cruises
fronts <- data.frame(
  cruise = c("G1", "G2", "G3"),
  salF = c(32.15, 32.5, 32.45),
  chlF = c(33, 36.2, 35)
)

# Function to visualize physiochemical variables of interest
getPlot <- function(flowData, sampData, variable, fronts = NULL, flowAlpha = 0.2) {
  
  # Select relevant columns from sample data
  samp <- sampData %>% select(lat, cruise, depth, !!sym(variable))
  # Create regions for flow data
  flow <- flowData %>% mutate(region = case_when(lat < 30 ~ "NSG_STZ", lat > 30 ~ "NTZ", TRUE ~ "NA"))
  
  # Plot
  p <- ggplot() +
    geom_point(data = flow, aes_string(y = variable, x = "lat", color = "cruise"), alpha = flowAlpha) +
    geom_smooth(data = flow, aes_string(y = variable, x = "lat"), 
                color = "#282A36", size = 0.5, method = "loess", se = FALSE, span = 0.1) +  
    scale_color_manual(values = c("G1" = "#50FA7B", "G2" = "#FFB86C", "G3" = "#BD93F9")) +
    facet_wrap(~cruise, nrow = 1, scales = "free_x") +
    xlim(20, 44) +
    theme_cowplot() +
    labs(y = variable, x = "Latitude °N") +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 18))
  
  # Add front lines if provided
  if (!is.null(fronts)) {
    for (i in 1:nrow(fronts)) {
      p <- p +
        geom_vline(data = filter(fronts, cruise == fronts$cruise[i]), aes(xintercept = salF), linetype = "solid", color = "black") +
        geom_vline(data = filter(fronts, cruise == fronts$cruise[i]), aes(xintercept = chlF), linetype = "dashed", color = "black")
    }
  }
  
  return(p)
}
# Function to save plots
savePlots <- function(object, dir, variable) {
  svg(paste0(dir, "g123_spatials_", variable, ".svg"), height = 3, width = 10)
  print(object)
  dev.off()
}

#-------------------------------------------------------------------------------
# Process Productivity data
path <- "data_in/cmap/productivity/"
prData <- list(
  g1 = read.csv(paste0(path, "KOK1606_Gradients1_Surface_O2Ar_NCP.csv")),
  g2 = read.csv(paste0(path, "MGL1704_Gradients2_Surface_O2Ar_NCP.csv")),
  g3 = read.csv(paste0(path, "KM1906_Gradients3_Surface_O2Ar_NCP.csv"))
)
cols <- c("time", "lat", "temp", "sal", "O2_Ar_sat", "NCP")
pr <- bind_rows(
  prData$g1 %>% select(all_of(cols)) %>% mutate(cruise = "G1"),
  prData$g2 %>% select(all_of(cols)) %>% mutate(cruise = "G2"),
  prData$g3 %>% select(all_of(cols)) %>% mutate(cruise = "G3")
) %>%
  filter(NCP >= 0) %>%
  mutate(lat = round(lat, 2)) %>%
  group_by(lat, cruise) %>%
  summarize(across(c(temp, sal, NCP), mean, na.rm = TRUE), .groups = "drop")

# Generate and save temperature plot
p <- getPlot(flowData = pr, sampData = meta, variable = "temp", fronts = fronts, flowAlpha = 0.2) %>% print
savePlots(object = p, dir = fig_dir, variable = "Temperature")
# NCP plot
p <- getPlot(flowData = pr, sampData = meta, variable = "NCP", fronts = fronts, flowAlpha = 0.2) %>% print
savePlots(object = p, dir = fig_dir, variable = "NetCommPro")
# Salinity plot
p <- getPlot(flowData = pr, sampData = meta, variable = "sal", fronts = fronts, flowAlpha = 0.2) %>% print
savePlots(object = p, dir = fig_dir, variable = "salinity")

#-------------------------------------------------------------------------------
path <- "data_in/cmap/particulates/"
paData <- list(
  g1 = read.csv(paste0(path, "Gradients1_KOK1606_PPPCPN_UW.csv")), 
  g2 = read.csv(paste0(path, "Gradients2_MGL1704_PPPCPN_UW.csv")),
  g3 = read.csv(paste0(path, "Gradients3_KM1906_PCPN_UW.csv"))
)
cols <- c("time", "lat", "depth", "pc", "pn")  
pa <- bind_rows(
  paData$g1 %>% select(all_of(cols)) %>% mutate(cruise = "G1"),
  paData$g2 %>% select(all_of(cols)) %>% mutate(cruise = "G2"),
  paData$g3 %>% select(all_of(cols)) %>% mutate(cruise = "G3")
) %>%
  filter(depth <= 15) %>%
  select(-depth) %>%
  group_by(lat, cruise) %>%
  summarize(across(c(pc, pn), mean, na.rm = TRUE), .groups = "drop")

# Generate and save particulate carbon plot
p <- getPlot(flowData = pa, sampData = meta, variable = "pc", fronts = fronts, flowAlpha = 0.6)
savePlots(object = p, dir = fig_dir, variable = "TotalCarbon")
# Particulate nitrogen plot
p <- getPlot(flowData = pa, sampData = meta, variable = "pn", fronts = fronts, flowAlpha = 0.6)
savePlots(object = p, dir = fig_dir, variable = "TotalNitrogen")

############################################################################################
# TRANSECT MAPS
############################################################################################

# Load metadata and filter by depth
meta <- smpMeta %>% filter(depth <= 15)

# Define color gradient for latitude
gradient <- colorRampPalette(c("#50FA7B", "#F1FA8C", "#FFB86C", "#BD93F9"))

# Create ggplot object and extract latitude color mapping
latCol.tmp <- ggplot(meta, aes(latitude, longitude, color = latitude)) +
  geom_point() +
  scale_color_gradientn(colors = gradient(100)) %>%
  ggplot_build() %>%
  .$data[[1]] %>%
  select(colour, x) %>%
  distinct()

# Merge color mapping with metadata
latCol <- merge(meta, latCol.tmp, by.x = "latitude", by.y = "x")

# Function to generate map plot
generate_map_plot <- function(data, cruise.type, sal_front = NULL, chla_front = NULL) {
  data %>%
    filter(cruise == cruise.type) %>%
    ggplot(aes(x = round(-longitude), y = latitude)) +
    borders("world", colour = "black", fill = "#663524", xlim = c(-200, -50), ylim = c(0, 200)) +
    geom_hline(yintercept = sal_front, linetype = "solid", color = "#FF7F00", linewidth = 2) +
    geom_hline(yintercept = chla_front, linetype = "solid", color = "#50FA7B", linewidth = 2) +
    geom_point(color = "black", size = 5.8) +
    geom_point(size = 4.5, aes(color = I(colour))) +
    coord_map(xlim = c(-161, -153), ylim = c(19, 43)) +
    labs(y = "Latitude °N") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.text.y = element_text(size = 20),
      panel.grid.major.x = element_blank()
    )
}

# Generate and save the plots
directory <- "3.0_"
cruise_types <- c("G1", "G2", "G3")  # Add cruise types
sal_fronts <- c(32.15, 32.5, 32.45)  # Define salinity fronts
chl_fronts <- c(33, 36.2, 35)        # Define chlorophyll fronts

p <- list()

for (i in seq_along(cruise_types)) {
  p[[cruise_types[i]]] <- generate_map_plot(latCol, cruise_types[i], sal_fronts[i], chl_fronts[i])
}

# Save the plots
for (i in seq_along(p)) {
  file_name <- paste0("results/", directory, "plots/", cruise_types[i], "_map_white.svg")
  svg(file_name, height = 8, width = 8)
  print(p[[cruise_types[i]]])
  dev.off()
}

 ############################################################################################
# UPDATE SAMPLE METAFILE WITH NORTH PACIFIC REGION CATEGORIES
############################################################################################

# Define cruise types and corresponding fronts
cruise_types <- c("G1", "G2", "G3")
sal_fronts <- c(32.15, 32.5, 32.45)
chl_fronts <- c(33, 36.2, 35)

# Function to determine region based on latitude and cruise
determine_region <- function(latitude, cruise) {
  cruise_index <- match(cruise, cruise_types)
  
  sal_front <- sal_fronts[cruise_index]
  chl_front <- chl_fronts[cruise_index]
  
  if (latitude < sal_front) {
    return("NPSG")
  } else if (latitude < chl_front) {
    return("STZ")
  } else {
    return("NTZ")
  }
}

# Apply region determination and update metadata
meta2 <- meta %>% mutate(region = mapply(determine_region, latitude, cruise))

# Check and save
dim(smpMeta); dim(meta2)
write.csv(meta2, "data_in/meta/g123_meta.csv", row.names = FALSE)

##############################################################################################
# Supplement Table on Samples
##############################################################################################

meta <- read_csv(paste0(dat_dir, "g123_meta_physioChems.csv"))
x <- smpMeta %>% select(sampleID, longitude, region)

x <- meta %>%
  left_join(., x, by = "sampleID") %>% 
  mutate(latitude = round(latitude),
         longitude = round(longitude)) %>% 
  group_by(cruise, latitude, depth, filter) %>%
  reframe(pn, pc, NCP, longitude, region, n_reps = n()) %>% 
  distinct %>% 
  mutate(filter = paste0(filter, "um")) %>%
  pivot_wider(
    names_from = filter,
    values_from = n_reps,
    names_prefix = "nReps_"
  ) %>% 
  distinct %>% 
  select(cruise, depth, latitude, longitude, region, everything()) %>% 
  arrange(cruise, latitude, depth)
x
write_csv(x, paste0(dat_dir, "g123_meta_suppTable"))

##############################################################################################
# Supplement Table on Matches
##############################################################################################

meta <- read_csv(paste0(dat_dir, "g123_meta_physioChems.csv"))
library("flextable")
summary <- meta %>%
  group_by(cruise) %>%
  summarise(
    total_samples = n(),
    PN = sum(!is.na(pn)),
    PC = sum(!is.na(pc)),
    NCP = sum(!is.na(NCP))
  ) %>% 
  bind_rows(
    summarise(., 
              cruise = "Total",
              total_samples = sum(total_samples),
              PN = sum(PN),
              PC = sum(PC),
              NCP = sum(NCP)
    )
  )

table <- summary %>%
  flextable() %>%
  autofit() %>%
  set_caption("Summary of Samples and Matches") 
table

save_as_image(
  table, 
  path = paste0(tab_dir, "match_Table.png"),
  webshot = "webshot2", 
  zoom = 2,  
  expand = 1
)
#

##############################################################################################
# Supplement Figure on Depth Profiles
##############################################################################################

station_info <- smpMeta %>% 
  select(sampleID, depth, station) %>% 
  mutate(collection = case_when(station == "underway" ~ "UW", # SIS: Surface-intake system
                                station == "tow-fish" ~ "TOW", 
                                T ~ "CTD")) %>% 
  mutate(collection = case_when(depth == "0" ~ "UW", T ~ collection)) %>% 
  select(-station, -depth)
physio <- read_csv(paste0(dat_dir, "g123_meta_physioChems.csv")) %>% 
  select(sampleID, NCP, pn, pc) %>% 
  rename(PON = pn, POC = pc)


dat <- meta %>%
  left_join(., station_info, by = "sampleID") %>% 
  select(sampleID, cruise, latitude, depth, filter, collection, region) %>%
  distinct()

times <- read.csv("data_in/meta/g123_smpTimes.csv") %>% 
  select(sampleID, time_UTC)

dat <- dat %>% left_join(., times, by = "sampleID") %>% 
  relocate(time_UTC, .before = sampleID) %>% 
  left_join(., physio, by = "sampleID")

write.csv(dat, paste0(dat_dir, "g123_ASV_SampleMeta.csv"), row.names = F)

p <- dat %>% 
  ggplot(aes(x = latitude, y = depth, fill = depth)) +
  geom_point(alpha = 0.7, size = 2, color = "black", shape = 21, stroke = 0.3) + 
  scale_y_reverse() +
  facet_wrap(~ cruise, ncol = 1) +
  scale_fill_gradient(
    name = "Depth (m)",
    low = "#a1d99b", 
    high = "#225ea8" 
  ) +
  labs(x = "Latitude (°N)", y = "Depth (m)",
       title = "Sample Site Depth Profiles") +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "gray90"),
        legend.position = "right") +
  facet_theme
p
svg(paste0(fig_dir, "sample_depthProfiles.svg"), height = 6, width = 4); p; dev.off()