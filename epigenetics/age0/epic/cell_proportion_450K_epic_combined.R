library(meffil)
library(psych)

### Birth 450K
load("/data/GENR3/Methylation/GENR_450KMETH_Norm_Release3/GENR_450KMETH_Release3_Betas_ALL_birth_20190813.RData")

# Estimate cell counts using combined reference panel
wbc_birth_450K_combined.data <-  as.data.frame(meffil.estimate.cell.counts.from.betas(x, "combined cord blood", verbose = T))

# Add participant ID
wbc_birth_450K_combined.data$SampleID <- colnames(x)

# Save cell count data
save(wbc_birth_450K_combined.data, file = "wbc_birth_450K_combined.RData")


### Birth EPIC
load("/home/599032/epic_qc/data/final/GENR_EPICv1METH_Norm_Betas_birth_ALL.RData")

# Estimate cell counts using combined reference panel
GENR_EPICv1METH_birth_CellTypes_combined.data <-  as.data.frame(meffil.estimate.cell.counts.from.betas(t(GENR_EPICv1METH_Norm_Betas_birth_ALL.data), "combined cord blood", verbose = T))

# Add participant ID
GENR_EPICv1METH_birth_CellTypes_combined.data$SampleID <- rownames(GENR_EPICv1METH_birth_CellTypes_combined.data)

# Save cell count data
save(GENR_EPICv1METH_birth_CellTypes_combined.data, file = "GENR_EPICv1METH_birth_CellTypes_combined.RData")

### Merge and compare 450K and EPIC
describe(wbc_birth_450K_combined.data)
describe(GENR_EPICv1METH_birth_CellTypes_combined.data)

wbc_birth_450K_epic_combined.data <- rbind(wbc_birth_450K_combined.data, GENR_EPICv1METH_birth_CellTypes_combined.data)
save(wbc_birth_450K_epic_combined.data, file = "wbc_birth_450K_epic_combined.RData")

### Birth 450K functional normalization with ALSPAC
load("~/alspac/GenR_METH_JointedQC_with_ALSPAC/Data/ch3_genr_180321.R")

# Estimate cell counts using combined reference panel
wbc_birth_450KwithALSPAC_combined.data <-  as.data.frame(meffil.estimate.cell.counts.from.betas(as.matrix(ch3_genr), "combined cord blood", verbose = T))

# Add participant ID
wbc_birth_450KwithALSPAC_combined.data$SampleID <- colnames(ch3_genr)

# Save cell count data
save(wbc_birth_450KwithALSPAC_combined.data, file = "wbc_birth_450KwithALSPAC_combined.RData")

### Merge and compare 450KwithALSPAC and EPIC
describe(wbc_birth_450KwithALSPAC_combined.data)
describe(wbc_birth_epic_combined.data)

wbc_birth_450KwithALSPAC_epic_combined.data <- rbind(wbc_birth_450KwithALSPAC_combined.data, wbc_birth_epic_combined.data)
save(wbc_birth_450KwithALSPAC_epic_combined.data, file = "wbc_birth_450KwithALSPAC_epic_combined.RData")