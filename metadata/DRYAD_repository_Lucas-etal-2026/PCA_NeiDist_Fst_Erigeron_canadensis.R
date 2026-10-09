################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# exploratory PCA, FST and Nei Distance  Analysis
################################################################################

#####----------------------------------------------------------------------#####

# Load required packages
packages <- c("adegenet", "vcfR", "dartR", "adegraphics", "ggmap", "scatterpie", 
              "ggplot2", "gridExtra", "gplots", "qvalue", "LEA", "remotes", "sf", 
              "tidyverse", "StAMPP", "geodist", "segmented", "pegas", "hierfstat", 
              "devtools", "rvalues")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

#####----------------------------------------------------------------------#####

vcf <- "Erigeron.canadensis.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf"
New_genind_conyza <- read.vcfR(vcf_path) %>% 
  vcfR2genlight() %>% 
  gl.compliance.check()

indNames(New_genind_conyza) <- gsub("^.*\\/", "", indNames(New_genind_conyza))
pop(New_genind_conyza) <- substr(indNames(New_genind_conyza), 4, 11)

metadata <- "Metadata.txt"
d <- read.table(metadata_path, sep = "\t", header = TRUE, fileEncoding = "latin1")
d <- d[!names(d) %in% c("Plate")]
d <- d[!duplicated(d$sample.ID2), ]

s <- data.frame(sample.ID2 = substr(indNames(New_genind_conyza), 4, 11))
d.meta <- merge(s, d, by = "sample.ID2", sort = FALSE, all.x = TRUE)


gl.pc <- adegenet::glPca(New_genind_conyza, nf = 150, useC = TRUE)
write.table(gl.pc["scores"])
            
my.grp <- as.factor(substr(indNames(New_genind_conyza), 34, 36))
do_1_PCA_plot(filename = NA, gl.pc = gl.pc, pca.scores = gl.pc$scores[, 1:2], 
              fac = my.grp, col.grp = funky(length(levels(my.grp))), alpha = 1, 
              labels = levels(my.grp), plot = TRUE)

my.grp <- as.factor(d.meta$Range)
do_1_PCA_plot(filename = NA, gl.pc = gl.pc, pca.scores = gl.pc$scores[, 1:2], 
              fac = my.grp, col.grp = c("blue", "red3"), alpha = 1, 
              labels = levels(my.grp), main = list("Native vs. Non-native", cex = 0.8), 
              plot = TRUE)

my.grp <- as.factor(d.meta$Region_stat)
do_1_PCA_plot(filename = NA, gl.pc = gl.pc, pca.scores = gl.pc$scores[, 1:2], 
              fac = my.grp, col.grp = funky(length(levels(my.grp))), alpha = 1, 
              labels = levels(my.grp), main = list("Regions", cex = 1), plot = TRUE)

my.grp <- as.factor(d.meta$country)
do_1_PCA_plot(filename = NA, gl.pc = gl.pc, pca.scores = gl.pc$scores[, 1:2], 
              fac = my.grp, col.grp = funky(length(levels(my.grp))), alpha = 1, 
              labels = levels(my.grp), main = list("Countries", cex = 1.1), plot = TRUE)

my.grp <- as.factor(d.meta$region.in.country)
do_1_PCA_plot(filename = NA, gl.pc = gl.pc, pca.scores = gl.pc$scores[, 1:2], 
              fac = my.grp, col.grp = funky(length(levels(my.grp))), alpha = 1, 
              labels = levels(my.grp), main = list("Region.in.Country", cex = 1.5), plot = TRUE)

NeisD <- stamppNeisD(New_genind_conyza, pop = TRUE)
NeisD_df <- as.data.frame(NeisD)
nj.tree <- nj(as.matrix(NeisD_df))

range <- as.numeric(substr(nj.tree$tip.label, 1, 3))
range <- ifelse(range < 107, "blue", "red3")
plot.phylo(nj.tree, type = "phylogram", cex = 0.5, rotate.tree = 300, tip.color = range)

pop(New_genind_conyza) <- substr(indNames(New_genind_conyza), 4, 9)
Fst <- stamppFst(New_genind_conyza, nboots = 0, percent = 95, nclusters = 1)
save(Fst, file = "Fst.RData")

phylogenies <- njs(as.dist(Fst))
plot.phylo(phylogenies, "phylogram", edge.width = 2, font = 3, tip.color = "blue", cex = 1.7)
