# DNA methylation at birth EPIC array README
## Data description
### Methylation data
Methylation at birth as assessed by the EPIC array is available in four .RData files. Each file can be read into R with load("GENR_EPICv1METH_Norm_\*.Rdata"). It is then available in the R environment as a data.frame with the format GENR_EPICv1METH_Norm_\*.data. The data requires about 7GB of memory, but in practice at least 20GB
should be reserved. The data.frame contains 1115 participants (rows) and 808,183 CpG sites (columns). All CpGs passed quality control and both autosomal and allosomal CpGs are included. Row.names are SampleIDs.

- **GENR_EPICv1METH_Norm_Betas_birth_ALL.RData**: Matrix of BETA methylation values. Normalized with quantile normalization and values can range from 0-1.
- **GENR_EPICv1METH_Norm_Betas_birth_3IQR3IQRwinsorized.RData**: BETA methylation values, but outliers exceeding 3rd Quartile+3*IQR or below 1st Quarile-3\*IQR were replaced with highest/lowest non-outlying value.
- **GENR_EPICv1METH_Norm_Betas_birth_3IQRwinsorized.RData**: BETA methylation values, but outliers exceeding 3rd Quartile+3*IQR or below 1st Quarile-3\*IQR were replaced with NA.
- **GENR_EPICv1METH_Norm_Mvalues_birth.RData**: M-value methylation values. Transformed Beta values to achieve better statistical properties.

See also [dnam_age0_epic_prep.R](https://github.com/inDEPTHlab/epigenetics/blob/main/epigenetics/age0/epic/dnam_age0_epic_prep.R) for preparation script.

### Auxiliary files
- **Selection_GENR_MethylEPIC_release1_birth_20230627.sav**: Selection file in SPSS format (can be read into R with foreign::read.spss). Contains both SampleIDs, which match row.names in GENR_EPICv1METH_Norm_\*.data, but also IDC to merge with the rest of the Generation R data. Includes also batch information, most importantly Sample_Plate.
- **GENR_EPICv1METH_birth_CellTypes_combined.RData**: Cell proportions in .RData format. Can be read in with load("wbc_birth_450K_epic_combined.RData") and matched/merged using SampleID. Reference panel based on FlowSorted.CordBloodCombined.450k. Created using [cell_proportion_450K_epic_combined.R](https://github.com/inDEPTHlab/epigenetics/blob/main/epigenetics/age0/epic/cell_proportion_450K_epic_combined.R).
- **GENR_EPICv1METH_birth_Descriptives.RData**: A data.frame with descriptive information for each CpG.

## General Analysis Workflow
1. Read in methylation data
2. Read in selection file
3. Read in cell proportions
4. Read in phenotype and covariate file
5. Merge selection file, cell proportion and phenotype/covariate information
6. Remove randomly sibling/or be sure to account for it in analysis (e.g. random family effect)
7. Match DNAm rows to those in the phenotype files (or merge data)
8. Optional: keep autosomal CpGs only using an annotation file (e.g. https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html)
9. Run EWAS: Use mclapply to loop through CpG sites in parallel and associate them with the phenotype of interest. Most common covariates to consider are sex, sample_plate, cell proportions, phenotype and/or gestational age, maternal smoking and an indicator of SES (e.g. maternal education). 

See also [smoking_ewas.R](https://github.com/inDEPTHlab/epigenetics/blob/main/epigenetics/age0/epic/smoking_ewas_epic.R) for an example EWAS on maternal smoking.

#### Filter autosomal probes
```
# Load annotation data
annotation.data <- ENmix::readmanifest("data/annotation/infinium-methylationepic-v-1-0-b5-manifest-file.csv")
autosomal_probes <- annotation.data$assay[annotation.data$assay$chr %in% c(1:22), "Name"]
autosomal_probes <- autosomal_probes[autosomal_probes %in% names(GENR_EPICv1METH_Norm_Betas_birth_ALL.data)]
GENR_EPICv1METH_Norm_Betas_birth_ALL_aut.data <- GENR_EPICv1METH_Norm_Betas_birth_ALL.data[autosomal_probes]; GENR_EPICv1METH_Norm_Betas_birth_ALL.data <- NULL; annotation.data <- NULL; gc()
```

