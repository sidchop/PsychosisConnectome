# Load necessary library for network analysis
library(NetworkToolbox)

# Load brain connectivity data
source("~/Dropbox/Sid/R_files/STAGES_difussion/scripts/functions/get_STAGES_sc_data.R")
data_all <- get_STAGES_sc_data(vectorised = FALSE, cca_covars = TRUE)

# Extract covariates and structural connectivity matrices
covars <- data_all[[2]] # covariates
con_mats <- data_all[[1]] # structural connectivity matrices

# Convert structural connectivity matrices into a 3D array
con_mats_3d <- abind::abind(con_mats, along = 3)

# Load function to get NBS connectome
source("~/Dropbox/Sid/R_files/STAGES_difussion/scripts/functions/get_nbs_connectome.R")

#mask by connections indcluded in primary group analysis
mask <- get_nbs_connectome(type = 'no_thr', zero_thr = 30)
mask[mask != 0] <- 1
mask[is.na(mask)] <- 1

# Apply mask to connectivity matrices
for (m in 1:dim(con_mats_3d)[3]) {
  con_mats_3d[,,m] <- con_mats_3d[,,m] * mask
}

# Load clinical data (change from baseline to 3 months)
source("~/Dropbox/Sid/R_files/STAGES_difussion/scripts/functions/get_STAGES_clin_data.R")
bl.clin <- get_STAGES_clin_data(ses = 1)
id_list <- bl.clin$BIDS_ID
m3.clin <- get_STAGES_clin_data(ses = 2)
bl.clin <- bl.clin[bl.clin$BIDS_ID %in% m3.clin$BIDS_ID,]

# Extract and clean up medication dose data
olz_dose <- m3.clin$Olz_cum_dose
olz_dose <- olz_dose[m3.clin$BIDS_ID != "sub-049"]
bl.clin <- bl.clin[bl.clin$BIDS_ID != "sub-049",]
m3.clin <- m3.clin[m3.clin$BIDS_ID != "sub-049",]

# Select relevant clinical score columns and merge datasets
bl.clin <- bl.clin[, c("BIDS_ID", "bprs_score", "sofas_score", "sans_score", "hamd_score", "hama_score", "whoq_score", "Age", "sex", "BIDS_ID")]
m3.clin <- m3.clin[, c("BIDS_ID", "bprs_score", "sofas_score", "sans_score", "hamd_score", "hama_score", "whoq_score")]
colnames(m3.clin) <- paste0(colnames(m3.clin), "_2")
clin <- cbind(m3.clin[,-1], bl.clin[,-1])

# Calculate change scores
cs_bl3m <- matrix(nrow = nrow(m3.clin), ncol = 6)
for (i in 1:6) {
  cs_bl3m[,i] <- round(((c(clin[,i]) - c(clin[,i+6])) / c(clin[,i])), 3)
}
cs_bl3m <- as.data.frame(cs_bl3m)
colnames(resid_cs_bl3m) <- colnames(cs_bl3m) <- c("bprs_score", "sofas_score", "sans_score", "hamd_score", "hama_score", "whoq_score")

# Visualize distribution of change scores
# Load required libraries
library(ggplot2)         # For creating plots
library(scales)          # Provides additional formatting options for plots
library(patchwork)       # For combining plots

# Plot distribution of change scores
# Plot 1: Distribution of SOFAS score changes over 3 months
p1 <- sofas_dist_3m <- ggplot(cs_bl3m, aes(x = sofas_score)) +
  geom_histogram(bins = 8, fill = "#4E79A7", color = "black", alpha = 0.7, position = position_dodge(width = 0.5)) +
  theme_classic() +
  labs(title = "", x = "SOFAS (Δ3)", y = "Frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) # Adjust x-axis breaks

# Plot 2: Distribution of BPRS score changes over 3 months
p2 <- bprs_dist_3m <- ggplot(cs_bl3m, aes(x = bprs_score)) +
  geom_histogram(bins = 8, fill = "#4E79A7", color = "black", alpha = 0.7, position = position_dodge(width = 0.5)) +
  theme_classic() +
  labs(title = "", x = "BPRS (Δ3)", y = "Frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) # Adjust x-axis breaks

# Combine and display the plots
p1 | p2

# Select patients for analysis
con_mats_3d.patients.3m <- con_mats_3d[,,c(covars$BIDS_ID %in% m3.clin$BIDS_ID)]

# Keep only essential objects in memory
rm(list = setdiff(ls(), c("resid_cs_bl3m", "cs_bl3m", "con_mats_3d", "con_mats_3d.patients.3m", "covars","data_all")))

# Clinical data processing (baseline to 12 months change)
# Load required functions
source("~/Dropbox/Sid/R_files/STAGES_difussion/scripts/functions/get_STAGES_sc_data.R")
source("~/Dropbox/Sid/R_files/STAGES_difussion/scripts/functions/get_STAGES_clin_data.R")

# Fetch baseline clinical data
bl.clin <- get_STAGES_clin_data(ses = 1)
id_list <- bl.clin$BIDS_ID

# Fetch 12-month clinical data
m12.clin <- get_STAGES_clin_data(ses = 3)

# Filter to include only matching IDs in both time points
bl.clin <- bl.clin[bl.clin$BIDS_ID %in% m12.clin$BIDS_ID,]

# Extract medication dose data
olz_dose <- m12.clin$Olz_cum_dose

# Select relevant columns from baseline and 12-month clinical data
bl.clin <- bl.clin[, c("BIDS_ID", "bprs_score", "sofas_score", "sans_score", "hamd_score", "hama_score", "whoq_score", "Age", "sex", "BIDS_ID")]
m12.clin <- m12.clin[, c("BIDS_ID", "bprs_score", "sofas_score", "sans_score", "hamd_score", "hama_score", "whoq_score")]

# Handle specific values in sans_score
m12.clin[which(m12.clin$sans_score == 0), 4] <- 1 # Replace 0 for proportional change scores

# Rename columns for clarity
colnames(m12.clin) <- paste0(colnames(m12.clin), "_2")

# Combine baseline and 12-month clinical data
clin <- cbind(m12.clin[,-1], bl.clin[,-1])

# Initialize matrices for change scores
cs_bl12m <- matrix(nrow = nrow(m12.clin), ncol = ncol(m12.clin) - 1)
colnames(cs_bl12m) <- colnames(m12.clin)[-1]

# Calculate residual and change scores
for (i in 1:ncol(clin)) {
  lm.fit <- lm(clin[,i] ~ bl.clin$Age + as.factor(bl.clin$sex) + bl.clin$BIDS_ID)
  resid_cs_bl12m[,i] <- resid(lm.fit)
  cs_bl12m[,i] <- clin[,i] - bl.clin[,i+1]
}

# Save 12m change scores
save(resid_cs_bl12m, cs_bl12m, file = "~/Dropbox/Sid/R_files/STAGES_difussion/scripts/results/cs_bl12m.RData")


#Connectome-based Predictive Modeling (CPM) w/ nested K-fold CV
# Load required libraries
library(doParallel)
library(foreach)

# Setting up parallel computing cluster with 9 cores
cl <- makeCluster(9)
registerDoParallel(cl)

# Initialize list to store CPM results
cpm_results <- list()

# Timing the execution of the parallel computation
system.time(
  cpm_results  <- foreach(o = 1:100, .packages = c("NetworkToolbox")) %dopar% {
    x <- list()
    for (s in 1:ncol(cs_bl3m)) {
      x[[s]] <-  cpmIV(neuralarray = con_mats_3d.patients.3m,
                       bstat = cs_bl3m[, s],
                       kfolds = 4,
                       nEdges = 10,
                       method = "sum",
                       connections = "separate",
                       corr = "spearman",
                       thresh = 0.05,
                       plots = FALSE,
                       standardize = TRUE,
                       cores = 1)
    }
    return(x)
  }
)

# Processing the results to extract positive and negative correlations
kf_r_pos <- matrix(NA, 100,6)
kf_r_neg <- matrix(NA, 100,6)
for (p in 1:length(cpm_results)) {
  for (s in 1:6) {
    kf_r_pos[p,s] <- cpm_results[[p]][[s]]$results[1,1]
    kf_r_neg[p,s] <- cpm_results[[p]][[s]]$results[2,1]
  }
}

# Naming the matrices for clarity
names(kf_r_pos) <- names(kf_r_neg) <- paste0(names(cs_bl3m), "_krr_sum_0.05_spearman_zscore_NBS_masked")

# Stopping the parallel cluster after use
stopCluster(cl)

# Storing the final results
cpm_results <- list(kf_r_pos, kf_r_neg)
names(cpm_results) <- c("kf_r_pos", "kf_r_neg")
saveRDS(cpm_results, "~/Dropbox/Sid/R_files/STAGES_difussion/scripts/NBS_paper/CPM_results/long_bl3m_raw_cs_kf_sum_0.05_spearman_zscore_NBS_masked.RDS")

# CPM (K-Fold) from Baseline to 12-month prediction
# --------------------------------------------------------------------

# Re-initializing parallel computing cluster
cl <- makeCluster(9)
registerDoParallel(cl)

# Resetting the results list
cpm_results <- list()

# Executing the parallel computation for 12-month prediction
system.time(
  cpm_results  <- foreach(o = 1:100, .packages = c("NetworkToolbox")) %dopar% {
    x <- list()
    for (s in 1:ncol(cs_bl3m)) {
      x[[s]] <-  cpmIV(neuralarray = con_mats_3d.patients.12m,
                       bstat = cs_bl12m[, s],
                       kfolds = 4,
                       nEdges = 10,
                       method = "sum",
                       connections = "separate",
                       corr = "spearman",
                       thresh = 0.05,
                       plots = FALSE,
                       standardize = TRUE,
                       cores = 1)
    }
    return(x)
  }
)

# Processing the 12-month results
kf_r_pos <- matrix(NA, 100,6)
kf_r_neg <- matrix(NA, 100,6)
for (p in 1:length(cpm_results)) {
  for (s in 1:6) {
    kf_r_pos[p,s] <- cpm_results[[p]][[s]]$results[1,1]
    kf_r_neg[p,s] <- cpm_results[[p]][[s]]$results[2,1]
  }
}

# Naming and storing the results
names(kf_r_pos) <- names(kf_r_neg) <- paste0(names(cs_bl3m), "_krr_sum_0.05_spearman_zscore_NBS_masked")
stopCluster(cl)

cpm_results <- list(kf_r_pos, kf_r_neg)
names(cpm_results) <- c("kf_r_pos", "kf_r_neg")
saveRDS(cpm_results, "~/Dropbox/Sid/R_files/STAGES_difussion/scripts/NBS_paper/CPM_results/long_bl12m_raw_cs_kf_sum_0.05_spearman_zscore_NBS_masked.RDS")
