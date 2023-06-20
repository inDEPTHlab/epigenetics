# Load packages
library(psych) #For descriptives
library(foreign) #To read spss files
library(missMDA) #For PCA based imputation
library(pbmcapply)
library(bacon)

# Load general data
general.data <- read.spss("~/perinatal_mediation/data/CHILD-ALLGENERALDATA_24102022.sav", to.data.frame=T)

# Load maternal smoking
smoking.data <- read.spss("~/perinatal_mediation/data/130815_ GEDRAGSGROEP Maternal smoking pregnancy.sav", to.data.frame=T)
general_smoking.data <- merge(general.data, smoking.data, by = "IDM", all = T)

# EPIC samples
samples_birth_epic.data <- read.spss("/home/599032/epic_qc/data/original/Selection_GENR_MethylEPIC_release1_birth_20230512.sav", to.data.frame = T)
samples_birth_epic.data <- samples_birth_epic.data[samples_birth_epic.data$EwasChildMethylEpic != 888, c("SampleID", "IDC", "Sample_Plate")]
general_smoking_samples.data <- merge(general_smoking.data, samples_birth_epic.data, by = "IDC")

# Merge cell count data
load("/home/599032/cell_proportion/wbc_birth_epic_combined.RData")
general_smoking_samples_counts.data <- merge(general_smoking_samples.data, wbc_birth_epic_combined.data, by = "SampleID")

# Reduce data to only relevant variables and complete maternal smoking
phenotype.data <- general_smoking_samples_counts.data[!is.na(general_smoking_samples_counts.data$msmoke),c("IDC", "MOTHER.x", "SampleID", "GENDER", "GESTBIR", "msmoke", "EDUCM", "CD8T", "NK", "CD4T", "Bcell", "Gran", "Mono", "nRBC", "Sample_Plate")]

# Keep randomly only one child per family
phenotype.data <- phenotype.data[sample(1:nrow(phenotype.data)),]
phenotype.data <- phenotype.data[!duplicated(phenotype.data$MOTHER.x), ]

# Convert education to continuous
phenotype.data$EDUCM <- as.numeric(phenotype.data$EDUCM)

# Make sure reference level for smoking is to never smoked
phenotype.data$msmoke <- relevel(phenotype.data$msmoke, ref = "never smoked")

# Set sample plate as factor
phenotype.data$Sample_Plate <- as.factor(phenotype.data$Sample_Plate)

# Check data
str(phenotype.data)
describe(phenotype.data)

# Variables to be included in imputation
variable_names <- c("GENDER", "GESTBIR", "msmoke", "EDUCM", "CD8T", "NK", "CD4T", "Bcell", "Gran", "Mono", "nRBC", "Sample_Plate")

# Maximum number of PCs, which can be fitted (number of variables - 1)
ncp.max <- length(variable_names)-1

# Use 10 fold cross-validation to determine optimum number of componencts
#set.seed(20230615)
#nb <- estim_ncpFAMD(phenotype.data[,variable_names], ncp.max = ncp.max, nbsim = 10)
#nb$ncp

# Save elbow plot
#png(file="figures/smoking_elbow_plot.png")
#plot(0:ncp.max, nb$criterion, xlab = "nb dim", ylab = "MSEP")
#dev.off()

# Perform PCA based imputation using the optimum number of components determined by CV
set.seed(20230423)
phenotype.imputePCA <- imputeFAMD(phenotype.data[,variable_names], ncp = 4)

# Recover the IDs, which were not included in imputation
phenotype_imputed.data <- as.data.frame(phenotype.imputePCA$completeObs)
phenotype_imputed.data <- cbind(phenotype.data[c("IDC", "SampleID")],phenotype_imputed.data)

# Make sure reference level for smoking is to never smoked
phenotype_imputed.data$msmoke <- relevel(phenotype_imputed.data$msmoke, ref = "never smoked")

# Load EPIC methylation data
load("data/final/dnam_age0_epic_beta.Rdata")

# Match methylation data participants to phenotype data
dnam_age0_epic_beta.data <- dnam_age0_epic_beta.data[as.character(phenotype_imputed.data$SampleID), ]
table(row.names(dnam_age0_epic_beta.data) == phenotype_imputed.data$SampleID)

regress <- function(cpg) {
  cpg_phenotype.data <- cbind(cpg, phenotype_imputed.data)
  coef(summary(lm(cpg ~ msmoke + GENDER + GESTBIR + CD8T + NK + CD4T + Bcell + Gran + Mono + nRBC + EDUCM + Sample_Plate, data = cpg_phenotype.data)))[2,]
}

smoking_ewas.list <- pbmclapply(dnam_age0_epic_beta.data, regress, mc.cores = 12)
smoking_ewas.data <- as.data.frame(do.call(rbind, smoking_ewas.list))
smoking_ewas.data$cpg <- names(dnam_age0_epic_beta.data)
save(smoking_ewas.data, file = "results/smoking_ewas.Rdata")

# Bacon
bc <- bacon(smoking_ewas.data["t value"])
bc
estimates(bc)
inflation(bc)
bias(bc)

png(file="figures/smoking_ewas_hist.png", width = 3000, height = 4000, pointsize = 100)
plot(bc, type="hist")
dev.off

png(file="figures/smoking_ewas_qq.png", width = 3000, height = 4000, pointsize = 100)
plot(bc, type="qq")
dev.off()
