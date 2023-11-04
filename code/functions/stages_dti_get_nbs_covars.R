stages_dti_get_nbs_covars <- function(center_age = T, log_hm = T) {
  data <- readxl::read_xlsx("~/Dropbox/Sid/PhD/Thesis/STAGES_data/Stages_master_long.xlsx", sheet = "dti_covars")
  data <- data[order(data$BIDS_ID),]
  data <- data[data$Included_in_analysis==1,]
  
  #get headmotion
  source("~/Dropbox/Sid/R_files/STAGES_difussion/scripts/functions/get_headmotion.R")
  hm <- get_headmotion(type = "absolute")
  colnames(hm) <- c("BIDS_ID", "headmotion")
  data <- merge(data, hm, by = "BIDS_ID")
  

    select <- c("group_bl", "Age", "sex", "headmotion")
    data <- data[,c(select)]
    if (center_age == T) {data$Age <- scale(data$Age, center = T, scale = F)}
    data$sex <- as.factor(data$sex-1)
    data$group_bl <- data$group_bl-1
    data$group_bl <- as.factor(data$group_bl)
    if (log_hm == T) {data$headmotion <- log(data$headmotion+0.001)}

  return(data)
}


