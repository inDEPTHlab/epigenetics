library(meffil)

load("data/final/dnam_age0_epic_beta.Rdata")

# Estimate cell counts using combined reference panel
wbc_birth_combined_450k.data <-  as.data.frame(meffil.estimate.cell.counts.from.betas(as.matrix(dnam_age0_epic_beta.data), "combined cord blood", verbose = T))
