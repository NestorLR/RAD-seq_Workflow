# Data and code from: Population genomics of a cosmopolitan weed provides insights into its local adaptation and recent demographic history

Dataset DOI: [10.5061/dryad.8w9ghx40k](10.5061/dryad.8w9ghx40k)

## Description of the data and file structure

This dataset accompanies the manuscript **“**&#x50;opulation genomics of a cosmopolitan weed provides insights into its local adaptation and recent demographic histor&#x79;**”**, published in *Molecular Ecology*. It includes genomic data, associated environmental variables, and R scripts used for data processing and analysis.

## 1. Genotype and Metadata Files

* Erigeron.canadensis.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf\
  Final SNP dataset in VCF format, containing 11,501 filtered SNPs across 640 individuals from 280 populations. SNPs were filtered for quality, minor allele frequency, missing data, and linkage (one SNP per contig).
* **`Erigeron_canadensis_metadata_cleaned.csv`** Metadata file providing population and environmental information (coordinates, range, region, climate variables, anthropogenic index, and performance data such as biomass and number of capitula). Variable definitions, units, and data sources are documented in the accompanying `Metadata_description.md` file.
* **`Metadata_description.md`** Data dictionary for `Erigeron_canadensis_metadata_cleaned.csv`. Describes each column, including the full meaning of abbreviations, units of measurement, and references for the climate, remote-sensing, and anthropogenic data sources.

## 2. R Scripts

These scripts reproduce the analyses from the manuscript, including population structure, genetic diversity, environmental associations, and genomic offset. Each script is commented and organized for reproducibility.

* `Admixture_Erigeron_canadensis.R` – Runs ADMIXTURE for ancestry inference
* `AMOVA_Erigeron_canadensis.R` – Performs analysis of molecular variance
* `Offset_field_Erigeron_canadensis.R` – Tests genomic offset against field performance
* `RDA_Erigeron_canadensis.R` – Performs dbRDA to explore genotype–environment associations
* `divMigrate_Erigeron_canadensis.R` – Estimates directional gene flow
* `divMigrate-27366003.log` – Output log from gene flow analysis
* `FEEMs_Erigeron_canadensis.R` – Fast Estimation of Effective Migration Surfaces
* `PCA_bio1_to_19_Erigeron_canadensis.R` – PCA of 19 BioClim variables
* `Erigeron.dist.out.nex` – Nexus-format genetic distance matrix
* `Trees_Erigeron_canadensis.R` – Constructs neighbor-joining trees
* `Mantel_Test_IBD_Erigeron_canadensis.R` – Partial Mantel tests for IBD/IBE
* `GeneticDiversity_Erigeron_canadensis.R` – Calculates Ho, He, Ar, FIS
* `RDA_Offset_Erigeron_canadensis.R` – Genomic offset estimation
* `cpDNA_loci_Erigeron_canadensis.R` – Placeholder for cpDNA analysis
* `PCA_NeiDist_Fst_Erigeron_canadensis.R` – PCA of Nei distance and Fst values

## 3. Software

All scripts were developed in R (version ≥ 4.0). The following packages are required: `adegenet`, `dartR`, `poppr`, `vegan`, `hierfstat`, `vcfR`, `phangorn`, `ape`, `ggtree`, `ggplot2`, and others as specified in each script.

### Files and variables

### File: `Erigeron.canadensis.SNPs.g66mac3minDP3bir15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf`

**Description:**\
This is the final filtered genotype dataset in Variant Call Format (VCF), containing 11,501 biallelic SNPs from 640 individuals across 280 populations. SNPs were filtered for quality, missingness, linkage (1 SNP per contig), and minor allele frequency. Standard VCF fields are included (CHROM, POS, REF, ALT, QUAL, FILTER, INFO), along with genotype data for each sample.

**Variables:**\
Standard VCF fields, plus sample genotypes in GT format.\
Missing data are represented as `./.`

---

### File: `Erigeron_canadensis_metadata_cleaned.csv`

**Description:**\
Tabular metadata for all sampled individuals and populations. Includes geographic coordinates, environmental variables, population classification (native/non-native), and field performance metrics (biomass, capitula count). Variables are described in detail in the associated manuscript. Missing values are denoted as `NA`.

**Key Variables:**

* `Population_ID`: Unique population identifier
* `Sample_ID`: Individual sample identifier
* `Latitude`, `Longitude`: Geographic location (decimal degrees)
* `Range`, `Region`, `Country`: Geopolitical and ecological groupings
* `CWD`, `Annual_Precip`, `Temp_Seasonality`, etc.: Climate-related variables
* `Human_Footprint`: Anthropogenic pressure index
* `Biomass`, `Capitula`: Field performance traits (log-transformed)

## Code

### File: `Admixture_Erigeron_canadensis.R`

**Description:**\
Runs ADMIXTURE analysis using the filtered SNP dataset. Includes code for cross-validation, optimal K selection, and plotting of ancestry proportions.

---

### File: `AMOVA_Erigeron_canadensis.R`

**Description:**\
Performs analysis of molecular variance (AMOVA) to partition genetic variation among predefined hierarchical groups (e.g., regions, ranges).

---

### File: `Offset_field_Erigeron_canadensis.R`

**Description:**\
Tests for a correlation between genomic offset values and observed fitness-related traits (biomass, capitula).

---

### File: `RDA_Erigeron_canadensis.R`

**Description:**\
Runs redundancy analysis (RDA) to assess genotype–environment associations using climate and human footprint variables.

---

### File: `divMigrate_Erigeron_canadensis.R`

**Description:**\
Estimates directional gene flow among regions using genetic distances. Outputs migration matrix for visualization.

---

### File: `divMigrate-27366003.log`

**Description:**\
Log file containing output from the divMigrate analysis run in the previous script.

---

### File: `FEEMs_Erigeron_canadensis.R`

**Description:**\
Calculates Fast Estimation of Effective Migration Surfaces

---

### File: `PCA_bio1_to_19_Erigeron_canadensis.R`

**Description:**\
Conducts PCA on 19 BioClim variables to reduce dimensionality before environmental association analysis.

---

### File: `Erigeron.dist.out.nex`

**Description:**\
Nexus-format genetic distance matrix for phylogenetic or clustering analyses (e.g., neighbor-joining trees).

---

### File: `Trees_Erigeron_canadensis.R`

**Description:**\
Builds and visualizes phylogenetic trees based on Nei distances among populations.

---

### File: `Mantel_Test_IBD_Erigeron_canadensis.R`

**Description:**\
Performs partial Mantel tests to evaluate isolation by distance and environment.

---

### File: `GeneticDiversity_Erigeron_canadensis.R`

**Description:**\
Calculates population-level genetic diversity statistics: observed/expected heterozygosity, allelic richness, and inbreeding coefficient.

---

### File: `RDA_Offset_Erigeron_canadensis.R`

**Description:**\
Estimates genomic offset using RDA-based projections of genotypes in non-native environments.

---

### File: `cpDNA_loci_Erigeron_canadensis.R`

**Description:**\
Placeholder script for analysis of chloroplast loci (no cpDNA sequences are included in this version).

---

### File: `PCA_NeiDist_Fst_Erigeron_canadensis.R`

**Description:**\
Performs PCA using pairwise Nei genetic distances and FST values between populations or regions.

#### Software

Any program that will open a spreadsheet, such as Excel is recommended.

 

RStudio (Software capable of running R) is used for statistical analysis and data visualization.
