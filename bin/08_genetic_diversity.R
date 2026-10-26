################################################################################
# Population genomics of a cosmopolitan weed provides insights into its
# local adaptation and recent demographic history
# Genetic diversity analysis
################################################################################

library(vcfR)
library(hierfstat)
library(adegenet)
library(dplyr)
library(stringr)

# -----------------------------------------------------------------------------#
# 1. Leer VCF y metadata
# -----------------------------------------------------------------------------#

vcf <- read.vcfR(
  "data/published/Erigeron.canadensis.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf"
)

genind_data <- vcfR2genind(vcf)

metadata <- read.csv(
  "data/published/Erigeron_canadensis_metadata_students.csv"
)

# Extraer el código de población del identificador de cada muestra
pop_data <- metadata %>%
  mutate(
    Population = str_extract(
      Ind_name_genind,
      "(?<=^\\d{3})[A-Z]+(?=_)"
    )
  )

# Comprobar identificadores y poblaciones
stopifnot(!anyDuplicated(pop_data$Ind_name_genind))
stopifnot(!anyNA(pop_data$Population))
stopifnot(!anyNA(pop_data$Region_stat))

# Reordenar metadata para que coincida exactamente con el VCF
idx <- match(indNames(genind_data), pop_data$Ind_name_genind)

if (anyNA(idx)) {
  stop("Hay individuos del VCF que no aparecen en el metadata.")
}

pop_data <- pop_data[idx, , drop = FALSE]

stopifnot(
  identical(indNames(genind_data), pop_data$Ind_name_genind)
)

# Asignar poblaciones a los individuos
pop(genind_data) <- factor(pop_data$Population)

# Guardar región como metadato
genind_data@other$Region_stat <- factor(pop_data$Region_stat)

# Comprobaciones
cat("Individuos:", nInd(genind_data), "\n")
cat("Poblaciones:", nPop(genind_data), "\n")
cat("Regiones:", n_distinct(pop_data$Region_stat), "\n")

# Verificar que cada población pertenezca a una sola región
population_region <- pop_data %>%
  distinct(Population, Region_stat)

if (anyDuplicated(population_region$Population)) {
  stop("Una o más poblaciones están asignadas a varias regiones.")
}

# -----------------------------------------------------------------------------#
# 2. Calcular índices de diversidad por población
# -----------------------------------------------------------------------------#

basic_stats <- hierfstat::basic.stats(genind_data)

allelic_richness <- hierfstat::allelic.richness(genind_data)

# Verificar que las columnas correspondan a las mismas poblaciones
stopifnot(
  setequal(colnames(basic_stats$Ho), colnames(basic_stats$Hs)),
  setequal(colnames(basic_stats$Ho), colnames(basic_stats$Fis)),
  setequal(colnames(basic_stats$Ho), colnames(allelic_richness$Ar))
)

# Reordenar las matrices para garantizar el mismo orden
population_order <- colnames(basic_stats$Ho)

basic_stats$Hs <- basic_stats$Hs[
  , population_order, drop = FALSE
]

basic_stats$Fis <- basic_stats$Fis[
  , population_order, drop = FALSE
]

allelic_richness$Ar <- allelic_richness$Ar[
  , population_order, drop = FALSE
]

# Promediar los valores por locus dentro de cada población
Ho_pop <- colMeans(basic_stats$Ho, na.rm = TRUE)
He_pop <- colMeans(basic_stats$Hs, na.rm = TRUE)
Fis_pop <- colMeans(basic_stats$Fis, na.rm = TRUE)
Ar_pop <- colMeans(allelic_richness$Ar, na.rm = TRUE)

# Tabla por población
diversity_pop <- data.frame(
  Population = population_order,
  Ho = unname(Ho_pop),
  He = unname(He_pop),
  Fis = unname(Fis_pop),
  Ar = unname(Ar_pop)
)

# Añadir región y número de muestras por población
population_metadata <- pop_data %>%
  group_by(Population, Region_stat) %>%
  summarise(
    n_samples = n_distinct(Ind_name_genind),
    .groups = "drop"
  )

diversity_pop <- diversity_pop %>%
  left_join(population_metadata, by = "Population")

stopifnot(!anyNA(diversity_pop$Region_stat))

# -----------------------------------------------------------------------------#
# 3. Diagnóstico adicional de Fis
# -----------------------------------------------------------------------------#

# Alternativa descriptiva: calcular Fis a partir de las medias Ho y He
# No es idéntica a la media de Fis por locus.
diversity_pop <- diversity_pop %>%
  mutate(
    Fis_from_mean_Ho_He = ifelse(
      is.finite(He) & He > 0,
      (He - Ho) / He,
      NA_real_
    )
  )

# Comparar ambos cálculos por población
print(
  diversity_pop %>%
    select(
      Population, Region_stat, Ho, He,
      Fis, Fis_from_mean_Ho_He
    )
)

# -----------------------------------------------------------------------------#
# 4. Resumir por región: media entre poblaciones +/- SD
# -----------------------------------------------------------------------------#

# Número total de individuos únicos por región
samples_region <- pop_data %>%
  group_by(Region_stat) %>%
  summarise(
    n_samples = n_distinct(Ind_name_genind),
    .groups = "drop"
  )

# Resumen de índices entre poblaciones
diversity_region <- diversity_pop %>%
  group_by(Region_stat) %>%
  summarise(
    n_populations = n_distinct(Population),
    
    Ho_mean = mean(Ho, na.rm = TRUE),
    Ho_sd = sd(Ho, na.rm = TRUE),
    
    He_mean = mean(He, na.rm = TRUE),
    He_sd = sd(He, na.rm = TRUE),
    
    Fis_mean = mean(Fis[is.finite(Fis)], na.rm = TRUE),
    Fis_sd = sd(Fis[is.finite(Fis)], na.rm = TRUE),
    
    Ar_mean = mean(Ar, na.rm = TRUE),
    Ar_sd = sd(Ar, na.rm = TRUE),
    
    # Diagnóstico: estimador alternativo de Fis
    Fis_from_mean_Ho_He_mean = mean(
      Fis_from_mean_Ho_He[is.finite(Fis_from_mean_Ho_He)],
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  left_join(samples_region, by = "Region_stat") %>%
  select(
    Region_stat,
    n_samples,
    n_populations,
    Ho_mean, Ho_sd,
    He_mean, He_sd,
    Fis_mean, Fis_sd,
    Fis_from_mean_Ho_He_mean,
    Ar_mean, Ar_sd
  ) %>%
  arrange(Region_stat)

print(diversity_region)

# -----------------------------------------------------------------------------#
# 5. Exportar resultados
# -----------------------------------------------------------------------------#

dir.create(
  "results/genetic_diversity",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  diversity_pop,
  "results/genetic_diversity/diversity_by_population.csv",
  row.names = FALSE
)

write.csv(
  diversity_region,
  "results/genetic_diversity/diversity_by_region.csv",
  row.names = FALSE
)

# -----------------------------------------------------------------------------#
# 6. Verificar totales
# -----------------------------------------------------------------------------#

cat(
  "Total de individuos únicos:",
  n_distinct(pop_data$Ind_name_genind), "\n"
)

cat(
  "Suma de muestras entre regiones:",
  sum(diversity_region$n_samples), "\n"
)

cat(
  "Total de poblaciones únicas:",
  n_distinct(pop_data$Population), "\n"
)

cat(
  "Suma de poblaciones entre regiones:",
  sum(diversity_region$n_populations), "\n"
)
