################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# IBD and IBE Mantel tests with environmental and spatial data
################################################################################

#####----------------------------------------------------------------------#####

# Load required packages
packages <- c("vegan", "ggplot2", "dplyr", "ggsci", "ade4", "reshape2", "tibble", "gridExtra", "ggpubr")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

#####----------------------------------------------------------------------#####

Clim <- read.table("meta_data.txt", header = TRUE, row.names = 1)

pca <- rda(Clim, scale = TRUE, center = TRUE)
summary(pca)
axes_pca <- scores(pca)
plot(pca, xlab = "PCA axis 1 (63%)", ylab = "PCA axis 2 (25.6%)")
biplot(pca, xlab = "PCA axis 1 (63%)", ylab = "PCA axis 2 (25.6%)")

climdist <- vegdist(axes_pca$sites, "euclidean")
climdistNames <- attr(climdist, "Labels")

Space <- read.table("DistMatConyza.txt", header = TRUE, row.names = 1)
Space <- as.matrix(Space) + 1
diag(Space) <- 1
Space <- log10(Space)
Space <- as.data.frame(Space)

Fst <- read.table("Mirror_matrix_0_NoNA.txt", header = TRUE, row.names = 1)
Fst <- as.matrix(Fst)
Fst_trans <- Fst / (1 - Fst)

geodist <- as.dist(Space)
fstdist_trans <- as.dist(Fst_trans)

keep <- climdistNames[climdistNames %in% attr(geodist, "Labels")]
climdist <- as.dist(as.matrix(climdist)[keep, keep])
geodist <- as.dist(as.matrix(geodist)[keep, keep])
fstdist_trans <- as.dist(as.matrix(fstdist_trans)[keep, keep])

mantel.randtest(geodist, fstdist_trans, nrep = 999)
plot(geodist, fstdist_trans, xlab = "geographic distance (km)", ylab = "genetic distance (Fst/(1-Fst))", pch = 19)

mantel.randtest(climdist, fstdist_trans, nrep = 999)
plot(climdist, fstdist_trans, xlab = "geographic distance (km)", ylab = "genetic distance (Fst/(1-Fst))", pch = 19)

mantel.partial(climdist, fstdist_trans, geodist, method = "pearson", permutations = 999)
mantel.partial(geodist, fstdist_trans, climdist, method = "pearson", permutations = 999)

x1 <- melt(as.data.frame(as.matrix(climdist)) %>% rownames_to_column()) %>% rename(CLIM = value)
x2 <- melt(as.data.frame(as.matrix(geodist)) %>% rownames_to_column()) %>% rename(GEO = value)
x3 <- melt(as.data.frame(as.matrix(fstdist_trans)) %>% rownames_to_column()) %>% rename(FST = value)

x <- merge(x1, x2, by = c("rowname", "variable")) %>% merge(x3, by = c("rowname", "variable")) %>% rename(orig = rowname, dest = variable)

md <- read.delim("Conyza_Populations_294.csv", sep = ",")
y <- md %>% select(X3.digit_code, Range) %>% rename(Population_ID = X3.digit_code) %>% mutate(Population_ID = gsub("-", "", Population_ID)) %>% unique()

z <- merge(x, y, by.x = "orig", by.y = "Population_ID") %>% merge(y, by.x = "dest", by.y = "Population_ID")

P1 <- z %>% filter(Range.x == Range.y) %>% ggplot(aes(x = GEO, y = CLIM, color = FST)) + geom_point(alpha = 0.5) + geom_smooth(method = "lm", color = "red") + scale_color_viridis_c(option = "turbo") + theme_bw() + facet_wrap(Range.x ~ Range.y)
P2 <- z %>% filter(Range.x == Range.y) %>% ggplot(aes(x = GEO, y = FST, color = CLIM)) + geom_point(alpha = 0.5) + geom_smooth(method = "lm", color = "red") + scale_color_viridis_c(option = "turbo") + theme_bw() + facet_wrap(Range.x ~ Range.y)
P3 <- z %>% filter(Range.x == Range.y) %>% ggplot(aes(x = CLIM, y = FST, color = 10^GEO)) + geom_point(alpha = 0.5) + geom_smooth(method = "lm", color = "red") + scale_color_viridis_c(option = "turbo", name = "GEO") + theme_bw() + facet_wrap(Range.x ~ Range.y)

ggpubr::ggarrange(P1, P2, P3, nrow = 3)

#####----------------------------------------------------------------------#####

geodist_mat <- as.matrix(geodist)
climdist_mat <- as.matrix(climdist)
fstdist_mat  <- as.matrix(fstdist_trans)

rownames(geodist_mat) <- colnames(geodist_mat) <- toupper(trimws(gsub("-", "", rownames(geodist_mat))))
rownames(climdist_mat) <- colnames(climdist_mat) <- toupper(trimws(gsub("-", "", rownames(climdist_mat))))
rownames(fstdist_mat)  <- colnames(fstdist_mat)  <- toupper(trimws(gsub("-", "", rownames(fstdist_mat))))

native_ids <- md %>% filter(tolower(Range) == "native") %>% pull(X3.digit_code) %>% gsub("-", "", .) %>% toupper() %>% trimws()
non_native_ids <- md %>% filter(tolower(Range) == "non-native") %>% pull(X3.digit_code) %>% gsub("-", "", .) %>% toupper() %>% trimws()

common_native_ids <- native_ids[native_ids %in% rownames(geodist_mat) & native_ids %in% rownames(climdist_mat) & native_ids %in% rownames(fstdist_mat)]
common_non_native_ids <- non_native_ids[non_native_ids %in% rownames(geodist_mat) & non_native_ids %in% rownames(climdist_mat) & non_native_ids %in% rownames(fstdist_mat)]

geodist_native <- geodist_mat[common_native_ids, common_native_ids]
climdist_native <- climdist_mat[common_native_ids, common_native_ids]
fstdist_native  <- fstdist_mat[common_native_ids, common_native_ids]

geodist_non_native <- geodist_mat[common_non_native_ids, common_non_native_ids]
climdist_non_native <- climdist_mat[common_non_native_ids, common_non_native_ids]
fstdist_non_native  <- fstdist_mat[common_non_native_ids, common_non_native_ids]

mantpart_native_ibd <- mantel.partial(as.dist(climdist_native), as.dist(fstdist_native), as.dist(geodist_native), method = "pearson", permutations = 999)
mantpart_native_ibe <- mantel.partial(as.dist(geodist_native), as.dist(fstdist_native), as.dist(climdist_native), method = "pearson", permutations = 999)

mantpart_non_native_ibd <- mantel.partial(as.dist(climdist_non_native), as.dist(fstdist_non_native), as.dist(geodist_non_native), method = "pearson", permutations = 999)
mantpart_non_native_ibe <- mantel.partial(as.dist(geodist_non_native), as.dist(fstdist_non_native), as.dist(climdist_non_native), method = "pearson", permutations = 999)

cat("Native IBD (control climate): r =", round(mantpart_native_ibd$statistic, 3), ", p =", mantpart_native_ibd$signif, "\n")
cat("Native IBE (control geography): r =", round(mantpart_native_ibe$statistic, 3), ", p =", mantpart_native_ibe$signif, "\n")
cat("Non-Native IBD (control climate): r =", round(mantpart_non_native_ibd$statistic, 3), ", p =", mantpart_non_native_ibd$signif, "\n")
cat("Non-Native IBE (control geography): r =", round(mantpart_non_native_ibe$statistic, 3), ", p =", mantpart_non_native_ibe$signif, "\n")

#####----------------------------------------------------------------------#####