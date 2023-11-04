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

# Read data for all subjects (choice between raw or preprocessed data)
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

