
# ============================================================
# Parte II. PCA genética de Erigeron canadensis
# Adaptado del análisis de Lucas et al. (2026)
# ============================================================

# 1. Cargar paquetes ------------------------------------------

library(vcfR)
library(adegenet)
library(ggplot2)

# 2. Definir rutas --------------------------------------------

vcf_file <- paste0(
  "data/published/",
  "Erigeron.canadensis.SNPs.g66mac3minDP3birl15rm75.",
  "F4.meanDP10maxDP1000maf05thin.recode.vcf"
)

metadata_file <- "data/published/Erigeron_canadensis_metadata_students.csv"

output_dir <- "results/genetic_PCA"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Importar metadatos ---------------------------------------

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 4. Importar VCF y convertir a genlight -----------------------

vcf <- read.vcfR(vcf_file, verbose = TRUE)

genotypes <- vcfR2genlight(vcf)

# Revisar dimensiones
cat("Individuos:", nInd(genotypes), "\n")
cat("Loci:", nLoc(genotypes), "\n")

# 5. Comprobar y vincular los metadatos -----------------------

genotype_ids <- indNames(genotypes)
metadata_ids <- metadata$Ind_name_genind

# Verificar identificadores únicos
stopifnot(!anyDuplicated(genotype_ids))
stopifnot(!anyDuplicated(metadata_ids))

# Reordenar metadatos según el orden de los genotipos
metadata_index <- match(genotype_ids, metadata_ids)

if (anyNA(metadata_index)) {
  stop("Hay individuos del VCF sin correspondencia en los metadatos.")
}

metadata_pca <- metadata[metadata_index, , drop = FALSE]

# Verificar que el orden coincida exactamente
stopifnot(
  identical(genotype_ids, metadata_pca$Ind_name_genind)
)

# Comprobar variables de agrupación
stopifnot(
  all(c("Range", "Region_stat") %in% names(metadata_pca))
)

cat("La correspondencia entre genotipos y metadatos es correcta.\n")

# 6. Calcular la PCA genética ---------------------------------

# Retenemos 10 ejes para explorar los primeros componentes.
# La gráfica utilizará PC1 y PC2.
set.seed(123)

pca_genetic <- adegenet::glPca(
  genotypes,
  nf = 10,
  useC = TRUE
)


# Porcentaje de variación explicado por los componentes
eig <- pca_genetic$eig

variance_explained <- 100 * eig / sum(eig)

data.frame(
  Component = paste0("PC", seq_along(eig)),
  Eigenvalue = eig,
  Variance_explained = variance_explained,
  Cumulative_variance = cumsum(variance_explained)
)[1:10, ]

# Etiquetas de los ejes
x_label <- paste0(
  "PC1 (", round(variance_explained[1], 2), "%)"
)

y_label <- paste0(
  "PC2 (", round(variance_explained[2], 2), "%)"
)

print(x_label)
print(y_label)


# Extraer los scores de los dos primeros componentes
pca_scores <- as.data.frame(pca_genetic$scores[, 1:2])

colnames(pca_scores) <- c("PC1", "PC2")



# Agregar metadatos respetando el orden de las muestras
pca_scores$Ind_name_genind <- genotype_ids
pca_scores$Range <- metadata_pca$Range
pca_scores$Region_stat <- metadata_pca$Region_stat

# Exportar los scores del PCA
write.csv(
  pca_scores,
  file.path(output_dir, "genetic_PCA_scores.csv"),
  row.names = FALSE
)

# Guardar el objeto completo para análisis posteriores
saveRDS(
  pca_genetic,
  file.path(output_dir, "genetic_PCA_object.rds")
)

# 7. Gráfica 1: PCA por rango nativo/no nativo ---------------

p_range <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = Range)
) +
  geom_point(alpha = 0.75, size = 2) +
  theme_classic(base_size = 12) +
  labs(
    title = "PCA genética de Erigeron canadensis",
    subtitle = "Individuos del rango nativo y no nativo",
    x = x_label,
    y = y_label,
    color = "Rango"
  )

ggsave(
  file.path(output_dir, "PCA_by_range.png"),
  plot = p_range,
  width = 8,
  height = 6,
  dpi = 300
)

# 8. Gráfica 2: PCA por región geográfica ---------------------

p_region <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = Region_stat)
) +
  geom_point(alpha = 0.75, size = 2) +
  theme_classic(base_size = 12) +
  labs(
    title = "PCA genética de Erigeron canadensis",
    subtitle = "Individuos coloreados por región geográfica",
    x = x_label,
    y = y_label,
    color = "Región"
  )

ggsave(
  file.path(output_dir, "PCA_by_region.png"),
  plot = p_region,
  width = 10,
  height = 7,
  dpi = 300
)

# 9. Mostrar gráficas -----------------------------------------

print(p_range)
print(p_region)

cat("\nAnálisis terminado.\n")
cat("Resultados guardados en:", output_dir, "\n")
