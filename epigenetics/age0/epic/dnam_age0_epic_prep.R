# Load librariries to read spss files, parallel computing and M value transformation
library(foreign)
library(ENmix)

# Load original files from lab
load("data/original/beta_all_QN.RData")

# Exclude participants set as missing by data management and transpose, so that CpG are columns, participants are row
sample.data <- read.spss("data/original/Selection_GENR_MethylEPIC_release1_birth_20230512.sav", to.data.frame = T)
samples_to_include <- sample.data[sample.data$EwasChildMethylEpic != 888, "SampleID"]
dnam_age0_epic_beta.data <- as.data.frame(t(beta_all_QC[,(colnames(beta_all_QC) %in% samples_to_include)]))
save(dnam_age0_epic_beta.data, file = "data/final/dnam_age0_epic_beta.Rdata")
beta_all_QC <- NULL; gc()

# Convert to M-values
dnam_age0_epic_M.data <- B2M(dnam_age0_epic_beta.data)
save(dnam_age0_epic_M.data, file = "data/final/dnam_age0_epic_M.Rdata")

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
dnam_age0_epic_beta_3IQR_winsorized.list <- lapply(dnam_age0_epic_beta.data,winsorize) 
dnam_age0_epic_beta_3IQR_winsorized.data <- as.data.frame(do.call(cbind, dnam_age0_epic_beta_3IQR_winsorized.list))
save(dnam_age0_epic_beta_3IQR_winsorized.data, file = "data/final/dnam_age0_epic_beta_3IQR_winsorized.Rdata")

# Trim at 3IQR
dnam_age0_epic_beta_3IQR_NA.list <- lapply(dnam_age0_epic_beta.data, function(cpg) {
  winsorize(cpg, replace = NA)
}) 
dnam_age0_epic_beta_3IQR_NA.data <- as.data.frame(do.call(cbind, dnam_age0_epic_beta_3IQR_NA.list))
save(dnam_age0_epic_beta_3IQR_NA.data, file = "data/final/dnam_age0_epic_beta_3IQR_NA.Rdata")

# Check whether distributions make sense
summary(dnam_age0_epic_beta.data[100000:100003])
summary(dnam_age0_epic_M.data[100000:100003])
summary(dnam_age0_epic_beta_3IQR_winsorized.data[100000:100003])
summary(dnam_age0_epic_beta_3IQR_NA.data[100000:100003])

# Check whether correlations between different transformations make sense
cor(data.frame(dnam_age0_epic_beta.data[100000], dnam_age0_epic_M.data[100000], dnam_age0_epic_beta_3IQR_winsorized.data[100000], dnam_age0_epic_beta_3IQR_NA.data[100000]), use = "pairwise.complete.obs")
cor(data.frame(dnam_age0_epic_beta.data[100001], dnam_age0_epic_M.data[100001], dnam_age0_epic_beta_3IQR_winsorized.data[100001], dnam_age0_epic_beta_3IQR_NA.data[100001]), use = "pairwise.complete.obs")
cor(data.frame(dnam_age0_epic_beta.data[100002], dnam_age0_epic_M.data[100002], dnam_age0_epic_beta_3IQR_winsorized.data[100002], dnam_age0_epic_beta_3IQR_NA.data[100002]), use = "pairwise.complete.obs")
cor(data.frame(dnam_age0_epic_beta.data[100003], dnam_age0_epic_M.data[100003], dnam_age0_epic_beta_3IQR_winsorized.data[100003], dnam_age0_epic_beta_3IQR_NA.data[100003]), use = "pairwise.complete.obs")

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

descriptives.list <- lapply(names(dnam_age0_epic_beta.data), descriptives)


