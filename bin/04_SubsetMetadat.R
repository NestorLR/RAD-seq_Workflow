# ============================================================
# Parte II. Preparación de metadatos para el análisis
# Estudio: Erigeron canadensis
# Fuente: Lucas et al. (2026), archivos de Dryad
# ============================================================

# 1. Rutas ----------------------------------------------------

input_file <- "data/published/Erigeron_canadensis_metadata_cleaned.csv"
output_dir <- "data/published"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Leer metadatos originales -------------------------------

metadata <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensiones de la tabla original:\n")
cat("Individuos:", nrow(metadata), "\n")
cat("Columnas:", ncol(metadata), "\n\n")
print(names(metadata))

# 3. Definir las columnas que conservaremos -------------------

# Identificación y distribución geográfica
id_columns <- c(
  "Ind_name_genind",
  "Range",
  "Region_stat",
  "Latitude",
  "Longitude"
)

# Las 19 variables bioclimáticas
bio_columns <- paste0("Bio", sprintf("%02d", 1:19))

# Otras variables ambientales
environment_columns <- c(
  "climatic_water_deficit",
  "human_footprint"
)

# Variables genéticas publicadas
genetic_columns <- c(
  "He",
  "Ho",
  "Fis",
  "PC1",
  "PC2",
  "Admixture1",
  "Admixture2",
  "Admixture3",
  "Admixture4"
)

# 4. Verificar que todas las columnas existan -----------------

selected_columns <- c(
  id_columns,
  bio_columns,
  environment_columns,
  genetic_columns
)

missing_columns <- setdiff(selected_columns, names(metadata))

# 5. Crear la tabla reducida ----------------------------------

metadata_students <- metadata[, selected_columns, drop = FALSE]

# 6. Exportar -------------------------------------------------

output_file <- file.path(
  output_dir,
  "Erigeron_canadensis_metadata_students.csv"
)

write.csv(
  metadata_students,
  output_file,
  row.names = FALSE,
  na = ""
)

# Guardar el listado de columnas para documentar el proceso
writeLines(
  selected_columns,
  file.path(output_dir, "metadata_columns_retained.txt")
)

# 7. Resumen de comprobación ----------------------------------

cat("\nPreparación terminada.\n")
cat("Individuos conservados:", nrow(metadata_students), "\n")
cat("Columnas conservadas:", ncol(metadata_students), "\n")
cat("Archivo generado:", output_file, "\n")

# Comprobar las primeras filas y los nombres
print(names(metadata_students))