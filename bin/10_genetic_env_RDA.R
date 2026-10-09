################################################################################
# Erigeron canadensis
# Integración de datos genéticos y ambientales mediante dbRDA
################################################################################

# 1. Paquetes ------------------------------------------------------------------

library(dplyr)
library(vcfR)
library(adegenet)
library(vegan)
library(scatterpie)
library(ggplot2)

# 2. Importar datos ------------------------------------------------------------

vcf <- read.vcfR(
  "data/published/Erigeron.canadensis.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf"
)

metadata <- read.csv(
  "data/published/Erigeron_canadensis_metadata_env.csv"
)

# Convertir el VCF a genind
genind_data <- vcfR2genind(vcf)

# Preparar identificadores
rownames(metadata) <- metadata$Ind_name_genind

ids_genetic <- indNames(genind_data)
ids_metadata <- metadata$Ind_name_genind


# Conservar muestras comunes, en el orden de la matriz genética
common_ids <- ids_genetic[ids_genetic %in% ids_metadata]

genind_data <- genind_data[common_ids, ]
env <- metadata[match(common_ids, metadata$Ind_name_genind), ]


cat("Muestras genéticas:", length(ids_genetic), "\n")
cat("Muestras con metadata:", length(ids_metadata), "\n")
cat("Muestras coincidentes:", length(common_ids), "\n")


# 4. Preparar matriz genética --------------------------------------------------

# Representación de alelos por individuo
freq <- as.data.frame(genind_data)

# Mantener las columnas sin valores faltantes, como en el script original
freq <- freq[, vapply(freq, function(x) !anyNA(x), logical(1)), drop = FALSE]

cat("Columnas genéticas conservadas:", ncol(freq), "\n")

# Asegurar orden idéntico
freq <- freq[match(env$Ind_name_genind, rownames(freq)), , drop = FALSE]


# 5. Seleccionar muestras completas para dbRDA --------------------------------

environmental_predictors <- c(
  "climatic_water_deficit",
  "Latitude",
  "New_Climatic_PC1",
  "human_footprint"
)

# Comprobar que las variables explicativas sean numéricas
if (!all(vapply(env[, environmental_predictors, drop = FALSE],
                is.numeric, logical(1)))) {
  stop("Las variables explicativas deben ser numéricas.")
}

complete_env <- complete.cases(
  env[, environmental_predictors, drop = FALSE]
)

env_dbrda <- env[complete_env, , drop = FALSE]
freq_dbrda <- freq[complete_env, , drop = FALSE]


cat("Muestras utilizadas en dbRDA:", nrow(freq_dbrda), "\n")


# 6. Ejecutar dbRDA ------------------------------------------------------------

dbrda_all <- vegan::capscale(
  freq_dbrda ~ climatic_water_deficit +
    Latitude +
    New_Climatic_PC1 +
    human_footprint,
  data = env_dbrda,
  distance = "bray"
)

# Resumen del análisis
print(dbrda_all)
anova(dbrda_all, permutations = 999)
anova(dbrda_all, by = "term", permutations = 999)

# Proporción de variación explicada por el modelo
print(RsquareAdj(dbrda_all))


# 7. Obtener coordenadas y flechas ---------------------------------------------

site_scores <- as.data.frame(
  vegan::scores(dbrda_all, display = "sites")
)

env_dbrda$CAP1 <- site_scores[, 1]
env_dbrda$CAP2 <- site_scores[, 2]

bi_dat <- as.data.frame(dbrda_all$CCA$biplot)


# 8. Varianza de los ejes ------------------------------------------------------

eig_constrained <- dbrda_all$CCA$eig
total_constrained <- sum(eig_constrained)
total_inertia <- dbrda_all$tot.chi

cap1_fitted <- 100 * eig_constrained[1] / total_constrained
cap1_total <- 100 * eig_constrained[1] / total_inertia

cap2_fitted <- 100 * eig_constrained[2] / total_constrained
cap2_total <- 100 * eig_constrained[2] / total_inertia


# 9. Graficar individuos y ancestría ------------------------------------------

admixture_cols <- paste0("Admixture", 1:4)

p_dbrda <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scatterpie::geom_scatterpie(
    data = env_dbrda,
    aes(
      x = CAP1,
      y = CAP2,
      group = Ind_name_genind,
      r = 0.04
    ),
    cols = admixture_cols
  ) +
  geom_segment(
    data = bi_dat,
    aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
    arrow = arrow(length = grid::unit(10, "pt")),
    linewidth = 0.8
  ) +
  geom_label(
    data = bi_dat,
    aes(x = CAP1, y = CAP2, label = Variable),
    fontface = "bold",
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      Admixture1 = "#dc3912",
      Admixture2 = "#ff9900",
      Admixture3 = "#3366cc",
      Admixture4 = "#109618"
    )
  ) +
  facet_grid(~Range) +
  labs(
    title = "dbRDA: variación genética y ambiente",
    x = paste0(
      "CAP1 (", round(cap1_fitted, 1),
      "% restringido / ", round(cap1_total, 1), "% total)"
    ),
    y = paste0(
      "CAP2 (", round(cap2_fitted, 1),
      "% restringido / ", round(cap2_total, 1), "% total)"
    ),
    fill = "Ancestría"
  ) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_dbrda)

# Guardar resultados
dir.create("results/genotype_environment", recursive = TRUE,
           showWarnings = FALSE)

ggsave(
  "results/genotype_environment/dbRDA_genetics_environment.png",
  p_dbrda,
  width = 11,
  height = 7,
  dpi = 300
)

write.csv(
  env_dbrda,
  "results/genotype_environment/dbRDA_sample_scores.csv",
  row.names = FALSE
)

write.csv(
  bi_dat,
  "results/genotype_environment/dbRDA_environmental_vectors.csv",
  row.names = FALSE
)

