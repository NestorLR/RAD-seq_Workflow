
# ============================================================
# Parte II. Inspeccionar matriz de  genotipos
# Erigeron canadensis
# ============================================================

# 1. Cargar paquetes ---------------------------------

library(vcfR)

# 2. Rutas -----------------------------------------------------

vcf_file <- paste0(
  "data/published/",
  "Erigeron.canadensis.SNPs.g66mac3minDP3birl15rm75.",
  "F4.meanDP10maxDP1000maf05thin.recode.vcf"
)

output_dir <- "data/vcf_stats/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# 3. Importar el VCF -------------------------------------------

vcf <- read.vcfR(vcf_file, verbose = TRUE)

# En vcfR, la primera columna de @gt contiene FORMAT.
# Las restantes corresponden a los individuos.
genotype_ids <- colnames(vcf@gt)[-1]

cat("\n--- MATRIZ GENÓMICA ---\n")
cat("Individuos:", length(genotype_ids), "\n")
cat("Variantes:", nrow(vcf@fix), "\n")

cat(
  "Identificadores duplicados:",
  sum(duplicated(genotype_ids)), "\n"
)


# 4. Revisar datos faltantes en los genotipos ------------------

# Extraer únicamente el campo GT (genotipo), no DP, GQ, etc.
gt <- extract.gt(
  vcf,
  element = "GT",
  as.numeric = FALSE
)

# Identificar llamadas faltantes:
# ./., .|., 0/., ./1, etc.
missing_gt <- is.na(gt) | grepl("\\.", gt)

# Porcentaje de genotipos faltantes por individuo
missing_by_individual <- data.frame(
  Ind_name_genind = colnames(gt),
  missing_genotypes = colSums(missing_gt),
  total_sites = nrow(gt),
  missing_percent = 100 * colMeans(missing_gt)
)

# Porcentaje de genotipos faltantes por variante

missing_by_variant <- data.frame(
  CHROM = vcf@fix[, "CHROM"],
  POS = vcf@fix[, "POS"],
  missing_genotypes = rowSums(missing_gt),
  total_individuals = ncol(gt),
  missing_percent = 100 * rowMeans(missing_gt)
)

# Crear un identificador de variante a partir de cromosoma y posición
missing_by_variant$variant <- paste(
  missing_by_variant$CHROM,
  missing_by_variant$POS,
  sep = "_"
)

# Reordenar columnas
missing_by_variant <- missing_by_variant[, c(
  "variant",
  "CHROM",
  "POS",
  "missing_genotypes",
  "total_individuals",
  "missing_percent"
)]

# Revisar resultados
head(missing_by_variant)
summary(missing_by_variant$missing_percent)



# Guardar resultados de control de calidad
write.csv(
  missing_by_individual,
  file.path(output_dir, "vcf_missing_by_individual.csv"),
  row.names = FALSE
)

write.csv(
  missing_by_variant,
  file.path(output_dir, "vcf_missing_by_variant.csv"),
  row.names = FALSE
)

