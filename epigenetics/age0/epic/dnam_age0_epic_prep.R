# Load librariries to read spss files, parallel computing and M value transformation
library(foreign)
library(ENmix)

# Load original files from lab
load("data/original/beta_all_QN.RData")

# Exclude participants set as missing by data management and transpose, so that CpG are columns, participants are row
sample.data <- read.spss("data/original/Selection_GENR_MethylEPIC_release1_birth_20230512.sav", to.data.frame = T)
samples_to_include <- sample.data[sample.data$EwasChildMethylEpic != 888, "SampleID"]
GENR_EPICv1METH_Norm_Betas_birth_ALL.data <- as.data.frame(t(beta_all_QC[,(colnames(beta_all_QC) %in% samples_to_include)]))
save(GENR_EPICv1METH_Norm_Betas_birth_ALL.data, file = "data/final/GENR_EPICv1METH_Norm_Betas_birth_ALL.RData")
beta_all_QC <- NULL; gc()

# Convert to M-values
GENR_EPICv1METH_Norm_Mvalues_birth.data <- B2M(GENR_EPICv1METH_Norm_Betas_birth_ALL.data)
save(GENR_EPICv1METH_Norm_Mvalues_birth.data, file = "data/final/GENR_EPICv1METH_Norm_Mvalues_birth.RData")

# Winsorize and Trim

##### N.B. Script adapted from: https://github.com/matthieugomez/statar/blob/master/R/winsorize.R 

#' Winsorize a numeric vector
# @param x A vector of values
# @param cutpoints Cutpoints under and above which are defined outliers. 
#  Default is (median - five times interquartile range, median + five times interquartile range) - N.b. I changed this to 3*IQR.
#  Compared to bottom and top percentile, this takes into account the whole distribution of the vector.
# @param probs A vector of probabilities that can be used instead of cutpoints. Quantiles are computed as the inverse of the empirical distribution function (type = 1)
# @param replace Values by which outliers are replaced. Default to cutpoints. A frequent alternative is NA.
# @param verbose Boolean. Should the percentage of replaced values printed?
# @examples                          
# v <- c(1:4, 99)
#' winsorize(v)
#' winsorize(v, replace = NA)
#' winsorize(v, probs = c(0.01, 0.99))
#' winsorize(v, cutpoints = c(1, 50))
# @export


winsorize <- function(x, probs = NULL, cutpoints = NULL , replace = c(cutpoints[1], cutpoints[2]), verbose = FALSE){
  dummy = is.integer(x)
  if (!is.null(probs)){
    stopifnot(is.null(cutpoints))
    stopifnot(length(probs)==2)
    cutpoints <- quantile(x, probs, type = 1, na.rm = TRUE)
  } else if (is.null(cutpoints)){
    l <- quantile(x, c(0.25, 0.50, 0.75), type = 1, na.rm = TRUE) 
    cutpoints <- c(l[1]-3*(l[3]-l[1]), l[3]+3*(l[3]-l[1]))  ### Default was Median+-5*IQR but has been changed to +-1*IQR+-3*IQR
  } else{
    stopifnot(length(cutpoints)==2)
  }
  if (is.integer(x)) cutpoints <- round(cutpoints)
  bottom <- x < cutpoints[1]
  top <- x > cutpoints[2]
  if (verbose){
    length <- length(x)
    message(paste(100*sum(bottom, na.rm = TRUE)/length,"% observations replaced at the bottom"))
    message(paste(100*sum(top, na.rm = TRUE)/length,"% observations replaced at the top"))
  }
  x[bottom] <- replace[1]
  x[top] <- replace[2]
  if (dummy){
    x <- as.integer(x)
  }
  x
}

# Winsorize at 3IQR
GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.list <- lapply(GENR_EPICv1METH_Norm_Betas_birth_ALL.data,winsorize) 
GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.data <- as.data.frame(do.call(cbind, GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.list))
row.names(GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.data) <- row.names(GENR_EPICv1METH_Norm_Betas_birth_ALL.data)
save(GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.data, file = "data/final/GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.RData")

# Trim at 3IQR
dnam_age0_epic_beta_3IQR_NA.list <- lapply(GENR_EPICv1METH_Norm_Betas_birth_ALL.data, function(cpg) {
  winsorize(cpg, replace = NA)
}) 
GENR_EPICv1METH_Norm_Betas_birth_3IQRNA.data <- as.data.frame(do.call(cbind, GENR_EPICv1METH_Norm_Betas_birth_3IQRNA.list))
row.names(GENR_EPICv1METH_Norm_Betas_birth_3IQRNA.data) <- row.names(GENR_EPICv1METH_Norm_Betas_birth_ALL.data)
save(GENR_EPICv1METH_Norm_Betas_birth_3IQRNA.data, file = "data/final/GENR_EPICv1METH_Norm_Betas_birth_3IQRNA.RData")

# Check whether distributions make sense
summary(GENR_EPICv1METH_Norm_Betas_birth_ALL.data[100000:100003])
summary(GENR_EPICv1METH_Norm_Mvalues_birth.data[100000:100003])
summary(GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.data[100000:100003])
summary(GENR_EPICv1METH_Norm_Betas_birth_3IQRNA.data[100000:100003])

# Check whether correlations between different transformations make sense
cor(data.frame(GENR_EPICv1METH_Norm_Betas_birth_ALL.data[100000], GENR_EPICv1METH_Norm_Mvalues_birth.data[100000], GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.data[100000], GENR_EPICv1METH_Norm_Betas_birth_3IQRNA.data[100000]), use = "pairwise.complete.obs")
cor(data.frame(GENR_EPICv1METH_Norm_Betas_birth_ALL.data[100001], GENR_EPICv1METH_Norm_Mvalues_birth.data[100001], GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.data[100001], GENR_EPICv1METH_Norm_Betas_birth_3IQRNA.data[100001]), use = "pairwise.complete.obs")
cor(data.frame(GENR_EPICv1METH_Norm_Betas_birth_ALL.data[100002], GENR_EPICv1METH_Norm_Mvalues_birth.data[100002], GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.data[100002], GENR_EPICv1METH_Norm_Betas_birth_3IQRNA.data[100002]), use = "pairwise.complete.obs")
cor(data.frame(GENR_EPICv1METH_Norm_Betas_birth_ALL.data[100003], GENR_EPICv1METH_Norm_Mvalues_birth.data[100003], GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.data[100003], GENR_EPICv1METH_Norm_Betas_birth_3IQRNA.data[100003]), use = "pairwise.complete.obs")

# Descriptives
descriptives <- function(cpg_name){
  cpg <- dnam_age0_epic_beta.data[,cpg_name]
  cpg_mean <- mean(cpg, na.rm = T)
  cpg_sd <- sd(cpg, na.rm = T)
  distribution <- quantile(cpg, c(0,0.25,0.50,0.75,1), na.rm = T)
  cpg_min <- distribution[1]
  cpg_quartile_1 <- distribution[2]
  cpg_median <- distribution[3]
  cpg_quartile_3 <- distribution[4]
  cpg_max <- distribution[5]
  cpg_iqr <- IQR(cpg)
  cpg_bottom_3iqr <- cpg_quartile_1 - 3*cpg_iqr
  cpg_top_3iqr <- cpg_quartile_3 + 3*cpg_iqr
  cpg_outlier_bottom_n <- table(cpg < cpg_bottom_3iqr)["TRUE"]
  cpg_outlier_top_n <- table(cpg > cpg_top_3iqr)["TRUE"]
  cpg_na <- table(is.na(cpg))["TRUE"]
  data.frame(cpg_name,cpg_mean,cpg_sd,cpg_min,cpg_quartile_1,cpg_median,cpg_quartile_3,cpg_max,cpg_iqr,cpg_bottom_3iqr,cpg_top_3iqr,cpg_outlier_bottom_n,cpg_outlier_top_n,cpg_na)
}

descriptives.list <- lapply(names(GENR_EPICv1METH_Norm_Betas_birth_ALL.data), descriptives)
GENR_EPICv1METH_birth_descriptives.data <- do.call(rbind, descriptives.list)
save(GENR_EPICv1METH_birth_descriptives.data, file = "data/final/GENR_EPICv1METH_birth_descriptives.RData")

