# NBS with ziGamma models

#### Load Libraries ####
library(parallel)
library(pbapply)
library(glmmTMB)
library(tidyverse)
library(boot)
library(parameters)
library(DHARMa)
library(sjstats)
library(multcomp)
library(car)
library(Rfast2)
library(fastglm)
library(Rfast)
source("scripts/functions/stages_dti_get_nbs_covars.R")

# Set Working Directory
setwd("~/Dropbox/Sid/R_files/STAGES_difussion/")

# Read in connectivity data for all subjects
connectomes <- list.files(path = "~/Dropbox/Sid/R_files/STAGES_difussion/data/connectomes/",
                          pattern = "*0.001*", full.names = TRUE, recursive = TRUE)

# Exclude subjects with specified IDs
exclude.t1 <- c("011", "025") # IDs for exclusion at T1
exclude.fov <- c("048","034c", "015", "060") # IDs for exclusion due to field of view issues; potentially add "063"

# Filter out excluded subjects
connectomes.qc <- grep(connectomes, pattern='*011_*|*025_*|*048_*|*034c_*|*015_*|*060_',
                       invert = TRUE, value = TRUE)

# Read connectivity matrices into a list
con_mat <- list()
for (i in seq_along(connectomes.qc)) {
  con_mat[[i]] <- read.table(connectomes.qc[[i]])
}

# Exclude ROIs with poor field of view coverage
removed_roi <- c(50, 51, 52, 53, 86, 90, 197, 199, 200, 201, 242, 243, 244)
con_mat_fov <- lapply(con_mat, function(mat) {
  temp <- as.data.frame(mat)
  temp[-removed_roi, -removed_roi]
})

# Calculate the number of subjects and ROIs
n_subj <- length(con_mat_fov) # Number of subjects
nroi <- nrow(con_mat_fov[[1]]) # Number of ROIs

# Prepare a matrix to vectorize each subject's connectivity matrix
con_mat_vec <- matrix(nrow = n_subj, ncol = nroi * (nroi - 1) / 2) # Use formula for the number of unique edges

# Vectorize lower triangle of each connectivity matrix
for (i in seq_len(n_subj)) {
  con_mat_vec[i, ] <- con_mat_fov[[i]][lower.tri(con_mat_fov[[i]], diag = FALSE)]
}

# Load covariates
covars <- stages_dti_get_nbs_covars(log_hm = TRUE, center_age = TRUE)

#### Zero-Inflated Gamma Model Definition ####

fast_ziGamma <- function(y, tol = 1e-07, maxiters = 100, full = TRUE) {
  x=data
  y1 <- y
  id <- which(y > 0)
  y1[id] <- 1
  prob <- Rfast::glm_logistic(x, y1, full = full, tol = tol,
                              maxiters = maxiters)
  p1 <- prob$info[2,4]
  z1 <- prob$info[2,3]

  if(p1==1) {
    y=y+10^-10     #add noise
  }
  x <- model.matrix(y ~ ., data = as.data.frame(x))
  mod <- summary(fastglm(x = x[id, ,drop = FALSE],y =  y[id], family = Gamma(link = "log"), tol = tol,maxit = maxiters))

  p2 <- try(mod$coefficients[2,4],silent = T)
  z2 <- try(mod$coefficients[2,3],silent = T)
  if(class(p2) == 'try-error') {p2 <- 1; z2 <- NA_integer_}
  return(list(min(p1,p2), absmax(c(z1,z2))))
}

# Run the fast version of the model on all edges using parallel computing
cl <- makeCluster(detectCores())
clusterExport(cl, list("data", "fastglm", "covars"))
fast_results <- pbapply(con_mat_vec, 2, fast_ziGamma)
fast_results_matrix <- do.call(rbind, fast_results)
colnames(fast_results_matrix) <- c("min_p", "max_z")
saveRDS(fast_results_matrix, "data/ziGamma_nbs/glmmTMB_ziGamma_fast.rds")

# Run a slower version on all edges which takes approximately 1 hour locally
# NOTE: The slow version underestimates standard errors as per https://github.com/glmmTMB/glmmTMB/issues/662
# However, it is still useful for model diagnostics, such as fits and outliers. We will use the fast version for the final results.
cl <- makeCluster(detectCores())
clusterExport(cl, list("glmmTMB", "data", "ziGamma", "testResiduals", "simulateResiduals"))

# Applying the model on each edge and handling different types of errors and convergence issues
results <- pbapply(con_mat_vec, 2, function(x) {
  try(mod <- glmmTMB(formula = x ~ group_bl + Age + sex + headmotion,
                     family = ziGamma(link = "log"),
                     data = data,
                     ziformula = ~ ., se = TRUE), silent = TRUE)

  if (class(mod) != 'try-error' && mod$fit$message == "singular convergence (7)") {
    try(mod <- glmmTMB(formula = (x + 10^-6) ~ group_bl + Age + sex + headmotion,
                       family = ziGamma(link = "log"),
                       data = data,
                       ziformula = ~ 0, se = TRUE), silent = TRUE)
  }
  if (class(mod) == 'try-error') {
    p <- m <- k <- d <- o <- NA
  } else {
    if (mod$fit$message == "iteration limit reached without convergence (10)" ||
        mod$fit$message == "false convergence (8)" ||
        mod$sdr$pdHess == FALSE) {
      m <- mod$fit$message
      p <- k <- d <- o <- NA
      return(c(p, k, d, o, m))
    } else {
      p <- min(parameters::p_value(mod)[2, 2], parameters::p_value(mod)[2, 7])
      m <- mod$fit$message
      x <- try(testResiduals(simulateResiduals(mod), plot = FALSE), silent = TRUE)
      if (class(x) == 'try-error') {
        k <- d <- o <- NA
        return(c(p, k, d, o, m))
      } else {
        k <- x$uniformity$p.value
        d <- x$dispersion$p.value
        o <- x$outliers$p.value
        return(c(p, k, d, o, m))
      }
    }
  }
}, cl = cl)

# Saving and reading results
# saveRDS(results, "data/temp/glmmTMB_ziGamma.rds")
results <- readRDS("data/ziGamma_nbs/glmmTMB_ziGamma.rds")
results_fast <- readRDS("data/ziGamma_nbs/glmmTMB_ziGamma_fast.rds")

# Create a vector with a certain percentage of non-zero values
pctNonZero <- function(x) { sum(x != 0) / length(x) }
pctNonZero_vec <- apply(con_mat_vec, 2, pctNonZero)

# Processing results for slow and fast functions and excluding certain edges based on criteria
results_p <- as.numeric(results[1, ]) # Slow function
results_p_fast <- as.numeric(results_fast[, 1]) # Fast function

pctZero <- function(x) { sum(x != 0) / length(x) }
pctZero_vec <- apply(con_mat_vec, 2, pctZero)

# Exclude edges with less than a certain percentage of non-zeros
results_p[pctZero_vec < 0.30] <- 1
results_p_fast[pctZero_vec < 0.30] <- 1

# Exclude NA and true 0 values
results_p[is.na(results_p)] <- 1
results_p_fast[is.na(results_p_fast)] <- 1

# Exclude significant Kolmogorov-Smirnov statistics (for QQ plot fit)
results_k <- as.numeric(results[2, ])
results_k[is.na(results_k)] <- 1
results_p[results_k < 0.05] <- 1
results_p_fast[results_k < 0.05] <- 1

# Exclude edges with significant over-dispersion
results_d <- as.numeric(results[3, ])
results_d[is.na(results_d)] <- 1
results_p[results_d < 0.05] <- 1
results_p_fast[results_d < 0.05] <- 1

# Exclude edges with outliers
results_o <- as.numeric(results[4, ])
results_o[is.na(results_o)] <- 1
results_p[results_o < 0.05] <- 1
results_p_fast[results_o < 0.05] <- 1

#Compute NBS

# Reading the maximum component list data from nulls NBS models run on HCP
max_comp_list <- readRDS("data/ziGamma_nbs/max_comp_nulls_zi30_ks_d_o.RDS")

# Sourcing functions from personal Dropbox folders
source("~/Dropbox/Sid/R_files/functions/vec_2_mat.R") # Converts a vector to a matrix
source("~/Dropbox/Sid/R_files/functions/get_components.R") # Retrieves components from a matrix

# Setting threshold values for primary threshold
p <- c(0.001, 0.01, 0.05)

# Looping over each primary threshold value
for (c in 1:3) {
  # Plotting histogram of maximum component sizes for each threshold
  hist(as.numeric(max_comp_list[[c]]),
       main = paste("Primary threshold:", p[c]),
       breaks = 30,
       xlab = "Maximum component size")

  # Adding a red line at the 95th percentile
  abline(v = quantile(as.numeric(max_comp_list[[c]]), .95), col = "red")

  # Processing the fast results to create a matrix and apply threshold
  P_mat <- vec_2_mat(results_p_fast, length = 319, diag = 0, lower.tri = TRUE, supress = TRUE)
  P_mat[P_mat >= p[c]] <- 0
  P_mat[P_mat != 0] <- 1

  # Finding the maximum component size in the zeta-inverse gamma model
  max_comp_ziGamma <- get_components(adj = P_mat, return_comp = FALSE, return_max = TRUE)

  # Adding a blue line at the 95th percentile for the ziGamma model
  abline(v = quantile(max_comp_ziGamma, .95), col = "blue")

  # Calculating and displaying the proportion of components
  message(paste("Proportion: ", sum(max_comp_list[[c]] >= max_comp_ziGamma) / length(max_comp_list[[c]])))
  message("Total size: ", sum(max_comp_ziGamma))

  # Compute percentages of positive and negative edges (pos z = hc > fep)
  max_comp_ziGamma <- get_components(adj = P_mat, return_comp = TRUE, return_max = FALSE)
  Z_mat_dir <- vec_2_mat(round(unlist(results_fast[, 2]), 3), length = 319, diag = 0, lower.tri = TRUE)
  Z_mat_dir <- max_comp_ziGamma * Z_mat_dir

  # Calculating percentages for specific conditions
  feplthc <- round((length((which(Z_mat_dir[upper.tri(Z_mat_dir)] > 0))) / sum(max_comp_ziGamma[upper.tri(max_comp_ziGamma)])) * 100, 2)
  fepgthc <- round((length((which(Z_mat_dir[upper.tri(Z_mat_dir)] < 0))) / sum(max_comp_ziGamma[upper.tri(max_comp_ziGamma)])) * 100, 2)

  # Displaying the percentages
  message(paste0(feplthc, "% of edges where fep < hc; and ", fepgthc, "% of edges where fep > hc"))
}
