# Detach loaded packages and clear environment
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc()
cat("\014")

##############################################################################################
# SECTION 8: PREPARATION AND ANALYSIS OF PHYTOPLANKTON AMPLICON AND PHYSIOCHEMICAL DATA
# Continued....
#
# This script processes amplicon data for eukaryotic (18S rRNA) and prokaryotic (16S rRNA) 
# communities, across multiple cruises (G1, G2, G3). It includes:
# 
# 4. Linear Mixed Modeling (LMM) of Environmental Variables: Using linear mixed modeling to 
#    understand how environmental variables (temperature, salinity, particulate carbon, particulate 
#    nitrogen, and net community production) influence phytoplankton community structure.
# * Steps 1-3 can be found in prior script where gmatrices were prepared
#
# The focus is on understanding the spatio-temporal distribution and persistence of key phytoplankton 
# communities and linking them to environmental drivers across the three annual surveys (2016, 2017, 2019).
# Linear mixed models are used to estimate variance components, fixed effects, and assess the relationship 
# between phytoplankton and physicochemical variables.
##############################################################################################

# LOAD PACKAGES
library("tidyverse")
library("RColorBrewer")
library("cowplot")
library("sommer")
library("ggrepel")
library("reshape2")

# Set Paths
dat_dir = "data_out/08_lmModelRun/dataframes/"
fig_dir = "data_out/08_lmModelRun/figures/"
tab_dir = "data_out/08_lmModelRun/tables/"
in_dir = "data_out/07_lmGmatrix/dataframes/"

# Load in Data
gmats <- readRDS(paste0(in_dir, "lmm_gMatrices_phyto_all_CLR.RDS"))
vars <- readRDS(paste0(in_dir, "lmm_vars_phyto_all.RDS"))

# Apply Naming convention
applyNames <- function(object = NULL){
  names(object) <- c("Archaeplastida", "Dinoflagellata", "Haptophyta", "Stramenopiles",  
                     "Prochlorococcus", "Synechococcus")
  return(object)
}

# LOAD FUNCTIONS
get_varComponents <- function(result_output = NULL){
  # Pull information
  model_summary <- summary(result_output)
  stat <- data.frame(predictor = rownames(summary(result_output)$varcomp),
                     phyto = c("Archaeplastida", "Dinoflagellata", "Haptophyta", "Stramenopiles", 
                               "Prochlorococcus", "Synechococcus", "Unexplained"),
                     varComponent = abs(summary(result_output)$varcomp$VarComp) / 
                       sum(abs(summary(result_output)$varcomp$VarComp))) %>% 
    mutate(percent = abs(varComponent * 100),
           percent_total = sum(percent)) %>% 
    dplyr::select(-predictor)
  
  stat$bic <- model_summary$logo$BIC 
  stat$aic <- model_summary$logo$AIC
  stat$Converge <- model_summary$logo$Converge
  
  # Make variance proportion plot
  hole <- 2.5
  df <- stat %>%
    mutate(csum = rev(cumsum(rev(percent))),
           pos = percent / 2 + lead(csum, 1),
           pos = if_else(is.na(pos), percent / 2, pos))
  p <- ggplot(df, aes(x = hole, y = percent, fill = phyto)) +
    geom_col(color = "#282A36") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c(
                                 "#49C5B1", # Archs #"#FF69B4", # Cyano
                                 "#DA4949", # Dino
                                 "#F1FA8C", # Haptos
                                 "#7570B3",
                                 "#ADECF9", # diatoms
                                 "#ff7f00",
                                 "#44475A")) + # ERROR
    xlim(c(1.15, hole + 1)) +
    labs(fill = paste0("Response: ", var)) +
    theme_void()
  print(p)
  # Pull stats and visual
  df <-  list(plot = p, stats = stat)
  return(df)
}

getBetas <- function(result_output = NULL){
  result_output <- model
  coefficients <- data.frame(
    physioChem = result_output$Beta$Trait,  
    effect = result_output$Beta$Effect,
    estimate = result_output$Beta$Estimate,
    stdError =  sqrt(diag(result_output$VarBeta))) %>% 
    mutate(sigs = (estimate - 1.96 * stdError > 0 & estimate + 1.96 * stdError > 0) |
             (estimate - 1.96 * stdError < 0 & estimate + 1.96 * stdError < 0),
           Color = ifelse(sigs & estimate > 0, "blue", ifelse(sigs & estimate < 0, "red", "black")))
  
  p <- ggplot(coefficients, aes(x = effect, y = estimate)) +
    geom_hline(yintercept = 0, color = "darkgray", size = 1) +
    geom_point(size = 4, aes(color = Color)) +
    scale_color_identity() +
    geom_errorbar(aes(ymin = estimate - 1.96 * stdError, 
                      ymax = estimate + 1.96 * stdError), width = 0.5) +
    labs(title = paste0("Fixed Effect Coefficients: ", var), 
         y = "Effect Size", 
         x = "Fixed Variables") +
    coord_flip() +
    theme_cowplot()
  
  print(p); return(p)
}

getResidualPlots <- function(filename = NULL, result_output = model, 
                             residuals = NULL, fitted_values = NULL){
  model_summary <- summary(result_output)
  
  png(filename, width = 8, height = 8, units = "in", bg = "white", res = 300)
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
  # Residuals vs. Fitted Values
  plot(residuals ~ fitted_values, main = "Residuals vs. Fitted Values")
  abline(h = 0, col = "red")
  mtext(paste0("Model Convergence: ", model_summary$logo$Converge), side = 1, 
        line = 4, cex = 0.8, adj = 0.5)
  # Generate histogram to see residual distribution, Add normality test
  hist(residuals, main = "Residuals Distribution", xlab = "Residuals")
  test_message <- sprintf("Shapiro-Wilk normality test:\nW = %f; p-value = %.7f",
                          shapiro.test(residuals)$statistic, shapiro.test(residuals)$p.value)
  mtext(test_message, side = 1, line = 5, cex = 0.8, adj = 0.5)
  # Generates QQ plot to again assess if residuals follow a normal distribution
  qqnorm(residuals, main = "QQ Plot of Residuals")
  qqline(residuals, col = "red")
  # Close png device
  dev.off()
}

getQQ <- function(residuals = NULL){
  res <- residuals
  qq_data <- qqnorm(res, plot.it = FALSE)
  qq_df <- data.frame(Theoretical = qq_data$x, Sample = qq_data$y)
  
  text_message <- sprintf("Shapiro-Wilk normality test:\nW = %f; p-value = %.7f",
                          shapiro.test(res)$statistic, shapiro.test(res)$p.value)
  
  
  p <- ggplot(qq_df, aes(x = Theoretical, y = Sample)) +
    geom_point(alpha = .5, size = 1.5) +
    geom_abline(intercept = mean(res), slope = sd(res), color = "#FF79C6", size = 1, alpha = 1) +
    labs(title = "Residual Quantile-Quantile", x = "Predicted Values", y = "Observed Values") +
    theme_cowplot() +
    theme(axis.title.x = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 20),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20, face = "bold"))
  p <- p + annotate("text", x = Inf, y = Inf, label = text_message, 
                    hjust = 1.0, vjust = 2, size = 5, color = "black", 
                    position = position_nudge(y = 0.5)) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"))  
  print(p); return(p)  
}

############################################################################################
############################################################################################
# PERFORM LINEAR MIXED MODELING OF PHYSIOCHEM (PC) VARIABLES
############################################################################################
############################################################################################
pc <- applyNames(gmats$carbon)
pn <- applyNames(gmats$nitrogen)
net <- applyNames(gmats$netCom)

top = 0
# ------------------------------------------------------------------------------------------
# PARTICULATE NITROGEN
# ------------------------------------------------------------------------------------------
set.seed(1)
var = "Nitrogen"
meta = vars$nitrogen
predictor = pn

meta$depth <- as.numeric(as.character(meta$depth))
meta$depth_category <- cut(meta$depth, breaks = c(-Inf, 15, 75, Inf), 
                           labels = c("0-15m", "45-75m", "90-125m"),
                           right = FALSE)

meta$cruise <- relevel(meta$cruise, ref = "G3")
meta$time_intrxn <- interaction(meta$cruise, meta$month, drop = TRUE)
meta$time_intrxn <- relevel(meta$time_intrxn, ref = "G3.4")

library("sommer")
model_nitrogen <- mmer(
  data = meta, 
  pn ~ filter + time_intrxn + depth_category, 
  random = ~ predictor[[1]] + 
    predictor[[2]] + 
    predictor[[3]] + 
    predictor[[4]] +
    predictor[[5]] +
    predictor[[6]],
  getPEV = T, rcov = ~units)
summary(model_nitrogen)

model = model_nitrogen
# RESIDUAL ANALYSIS
getResidualPlots(filename = paste0(fig_dir, "/pclmm_", var, "_res_t", top, ".png"),
                 residuals = model$residuals, fitted_values = model$fitted)
# EXTRACT QQ
p <- getQQ(residuals = model$residuals); p
svg(paste0("results/", directory, "/", var, "/pclmm_", var, "_QQ_t", top, ".svg"), 
    height = 5, width = 6); p; dev.off()

# VARIANCE COMPONENTS
p <- get_varComponents(result_output = model); p
pn_stats <- p$stats
svg(paste0(fig_dir, "pclmm_", var, "_varComp_t", top, ".svg"), height = 8, width = 8); p$plot; dev.off()

# FIXED VARIABLE COEFFICIENTS
p <- getBetas(result_output = model); p
p; svg(paste0(fig_dir, "pclmm_", var, "_Coeff_t", top, ".svg"), height = 5, width = 8); p; dev.off()

# ------------------------------------------------------------------------------------------
# PARTICULATE Carbon
# ------------------------------------------------------------------------------------------
set.seed(1)
var = "Carbon"
meta = vars$carbon
predictor = pc

meta$depth <- as.numeric(as.character(meta$depth))
meta$depth_category <- cut(meta$depth, breaks = c(-Inf, 15, 75, Inf), 
                           labels = c("0-15m", "45-75m", "90-125m"),
                           right = FALSE)

meta$cruise <- relevel(meta$cruise, ref = "G3")
meta$time_intrxn <- interaction(meta$cruise, meta$month, drop = TRUE)
meta$time_intrxn <- relevel(meta$time_intrxn, ref = "G3.4")

model_carbon <- mmer(
  data = meta, 
  pc ~ filter + time_intrxn + depth_category, 
  random = ~ predictor[[1]] + 
    predictor[[2]] + 
    predictor[[3]] + 
    predictor[[4]] +
    predictor[[5]] +
    predictor[[6]],
  getPEV = T, rcov = ~units)
summary(model_carbon)

model = model_carbon
# RESIDUAL ANALYSIS
getResidualPlots(filename = paste0(fig_dir, "/pclmm_", var, "_res_t", top, ".png"),
                 residuals = model$residuals, fitted_values = model$fitted)
# EXTRACT QQ
p <- getQQ(residuals = model$residuals); p
svg(paste0("results/", directory, "/", var, "/pclmm_", var, "_QQ_t", top, ".svg"), 
    height = 5, width = 6); p; dev.off()

# VARIANCE COMPONENTS
p <- get_varComponents(result_output = model); p
pn_stats <- p$stats
svg(paste0(fig_dir, "pclmm_", var, "_varComp_t", top, ".svg"), height = 8, width = 8); p$plot; dev.off()

# FIXED VARIABLE COEFFICIENTS
p <- getBetas(result_output = model); p
p; svg(paste0(fig_dir, "pclmm_", var, "_Coeff_t", top, ".svg"), height = 5, width = 8); p; dev.off()


# ------------------------------------------------------------------------------------------
# NetComm
# ------------------------------------------------------------------------------------------
set.seed(1)
var = "NetComm"
meta = vars$netComm
predictor = net

meta$depth <- as.numeric(as.character(meta$depth))
meta$depth_category <- cut(meta$depth, breaks = c(-Inf, 15, 75, Inf), 
                           labels = c("0-15m", "45-75m", "90-125m"),
                           right = FALSE)

meta$cruise <- relevel(meta$cruise, ref = "G3")
meta$time_intrxn <- interaction(meta$cruise, meta$month, drop = TRUE)
meta$time_intrxn <- relevel(meta$time_intrxn, ref = "G3.4")

model_netComm <- mmer(
  data = meta, 
  NCP ~ filter + time_intrxn + depth_category, 
  random = ~ predictor[[1]] + 
    predictor[[2]] + 
    predictor[[3]] + 
    predictor[[4]] +
    predictor[[5]] +
    predictor[[6]],
  getPEV = T, rcov = ~units)
summary(model_netComm)

model = model_netComm

# RESIDUAL ANALYSIS
getResidualPlots(filename = paste0(fig_dir, "/pclmm_", var, "_res_t", top, ".png"),
                 residuals = model$residuals, fitted_values = model$fitted)
# EXTRACT QQ
p <- getQQ(residuals = model$residuals); p
svg(paste0("results/", directory, "/", var, "/pclmm_", var, "_QQ_t", top, ".svg"), 
    height = 5, width = 6); p; dev.off()

# VARIANCE COMPONENTS
p <- get_varComponents(result_output = model); p
pn_stats <- p$stats
svg(paste0(fig_dir, "pclmm_", var, "_varComp_t", top, ".svg"), height = 8, width = 8); p$plot; dev.off()

# FIXED VARIABLE COEFFICIENTS
p <- getBetas(result_output = model); p
p; svg(paste0(fig_dir, "pclmm_", var, "_Coeff_t", top, ".svg"), height = 5, width = 8); p; dev.off()


############################################################################################
# PIT Residual Analysis -- IN DEVELOPMENT!! 
############################################################################################
# Detach loaded packages and clear environment
lapply(names(sessionInfo()$otherPkgs), function(pkg) {
  detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
});  cat("\014")

# Load up your libraries
library("tidyverse")
library("RColorBrewer")
library("cowplot")
library("topmodels") # https://topmodels.r-forge.r-project.org/topmodels/

model = model_nitrogen
residuals <- model$residuals
hist(residuals, main = "Residual Histogram", breaks = 20)
qqnorm(residuals); qqline(residuals, col = "red")
shapiro.test(residuals)

# Compute PIT values
residuals <- model$residuals
pit_values <- pnorm(residuals, mean = mean(residuals), sd = sd(residuals))
hist(pit_values, main = "PIT Histogram", breaks = 20)
qqnorm(pit_values); qqline(pit_values, col = "blue")
shapiro.test(pit_values)

#-------------------------------------------------------------------------------
model = model_carbon
residuals <- model$residuals
hist(residuals, main = "Residual Histogram", breaks = 20)
qqnorm(residuals); qqline(residuals, col = "red")
shapiro.test(residuals)

# Compute PIT values
residuals <- model$residuals
pit_values <- pnorm(residuals, mean = mean(residuals), sd = sd(residuals))
hist(pit_values, main = "PIT Histogram", breaks = 20)
qqnorm(pit_values); qqline(pit_values, col = "blue")
shapiro.test(pit_values)

#-------------------------------------------------------------------------------
model = model_netComm
residuals <- model$residuals
hist(residuals, main = "Residual Histogram", breaks = 20)
qqnorm(residuals); qqline(residuals, col = "red")
shapiro.test(residuals)

# Compute PIT values
residuals <- model$residuals
pit_values <- pnorm(residuals, mean = mean(residuals), sd = sd(residuals))
hist(pit_values, main = "PIT Histogram", breaks = 20)
qqnorm(pit_values); qqline(pit_values, col = "blue")
shapiro.test(pit_values)
