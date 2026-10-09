################################################################################
# Population genomics of a cosmopolitan weed provides insights into its
# local adaptation and recent demographic history
#
# PCA ambiental: variables bioclimáticas Bio01-Bio19
# Datos: metadata de Erigeron canadensis
################################################################################

# 1. Cargar paquetes -----------------------------------------------
library(tidyverse)
library(factoextra)
library(ggplot2)

# 2. Importar metadata ---------------------------------------------------------

metadata <- read.csv(
  "data/published/Erigeron_canadensis_metadata_students.csv"
)

# Identificar las variables bioclimáticas
bio_columns <- paste0("Bio", sprintf("%02d", 1:19))

# Extraer únicamente las variables ambientales
bio_data <- metadata %>%
  select(all_of(bio_columns))

# Comprobar que las variables sean numéricas
if (!all(vapply(bio_data, is.numeric, logical(1)))) {
  stop("Una o más variables bioclimáticas no son numéricas.")
}

# Comprobar valores faltantes
if (anyNA(bio_data)) {
  stop(
    "Hay valores faltantes en las variables ambientales. ",
    "Revísalos antes de ejecutar el PCA."
  )
}


# PCA usando únicamente filas completas

complete_rows <- complete.cases(bio_data)

bio_data_complete <- bio_data[complete_rows, ]

metadata_pca <- metadata[complete_rows, ]

# Verificar cuántas muestras se conservaron
cat("Muestras originales:", nrow(metadata), "\n")
cat("Muestras con datos completos:", nrow(metadata_pca), "\n")
cat(
  "Muestras excluidas:",
  nrow(metadata) - nrow(metadata_pca), "\n"
)


# 3. Ejecutar PCA ambiental ----------------------------------------------------

# Estandarizar variables para que tengan varianza comparable
pca_environmental <- prcomp(
  bio_data_complete,
  center = TRUE,
  scale. = TRUE
)

# Porcentaje de varianza explicado por cada componente
variance_explained <- 100 *
  pca_environmental$sdev^2 /
  sum(pca_environmental$sdev^2)

# Resumen del PCA
summary(pca_environmental)

# Porcentaje explicado por los dos primeros componentes
cat(
  "PC1:", round(variance_explained[1], 2), "%\n",
  "PC2:", round(variance_explained[2], 2), "%\n",
  "PC1 + PC2:", round(sum(variance_explained[1:2]), 2), "%\n"
)


# 4.  Asignar PC1 a las muestras que participaron en el PCA.
# Las demás conservan NA
metadata$New_Climatic_PC1 <- NA_real_
metadata$New_Climatic_PC1[complete_rows] <- pca_environmental$x[, 1]




# 5. Gráfico de individuos -----------------------------------------------------

# Incorporar información de procedencia

pca_scores <- as.data.frame(pca_environmental$x) %>%
  mutate(
    Ind_name_genind = metadata_pca$Ind_name_genind,
    Range = metadata_pca$Range,
    Region_stat = metadata_pca$Region_stat
  )


pca_ind_plot <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = Range)
) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "PCA ambiental de Erigeron canadensis",
    subtitle = "Variables bioclimáticas Bio01-Bio19",
    x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
    color = "Rango"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(pca_ind_plot)


# 6. Gráfico de variables ambientales -----------------------------------------

pca_var_plot <- fviz_pca_var(
  pca_environmental,
  col.var = "contrib",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE
) +
  labs(
    title = "Contribución de las variables bioclimáticas",
    x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 2), "%)")
  )

print(pca_var_plot)


# 7. Guardar resultados --------------------------------------------------------

dir.create("results/environmental_PCA", recursive = TRUE, showWarnings = FALSE)


ggsave(
  "results/environmental_PCA/environmental_PCA_individuals.png",
  pca_ind_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "results/environmental_PCA/environmental_PCA_variables.png",
  pca_var_plot,
  width = 8,
  height = 6,
  dpi = 300
)

write.csv(
  pca_scores,
  "results/environmental_PCA/environmental_PCA_scores.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    Component = paste0("PC", seq_along(variance_explained)),
    Variance_explained_percent = variance_explained
  ),
  "results/environmental_PCA/environmental_PCA_variance.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(pca_environmental$rotation),
  "results/environmental_PCA/environmental_PCA_loadings.csv",
  row.names = TRUE
)

write.csv(
  metadata,
  "data/published/Erigeron_canadensis_metadata_env.csv",
  row.names = FALSE
)

