
# ============================================================
# Parte II. Ancestría genética de Erigeron canadensis (K = 4)
# Adaptado de Lucas et al. (2026)
#
# Figuras:
# 1. Ancestría por rango nativo / no nativo
# 2. Ancestría por región geográfica
#
# Los datos de ancestría proceden de los metadatos publicados.
# Este script no vuelve a estimar los componentes de ancestría.
# ============================================================


# 1. Paquetes -------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scatterpie)
library(maps)
library(stringr)

# 2. Rutas ----------------------------------------------------

metadata_file <- paste0(
  "data/published/",
  "Erigeron_canadensis_metadata_students.csv"
)

output_dir <- "results/admixture"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Importar metadatos ---------------------------------------

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ancestry_columns <- paste0("Admixture", 1:4)

required_columns <- c(
  "Ind_name_genind",
  "Range",
  "Region_stat",
  "Latitude",
  "Longitude",
  ancestry_columns
)


# 4. Preparar y validar los datos -----------------------------

# Conservar las columnas necesarias para ambas figuras
admixture_data <- metadata %>%
  select(
    Ind_name_genind,
    Range,
    Region_stat,
    all_of(ancestry_columns)
  )


# Comprobar que las proporciones sumen aproximadamente 1
admixture_data <- admixture_data %>%
  mutate(
    ancestry_sum = rowSums(across(all_of(ancestry_columns)))
  )

cat("\nResumen de la suma de proporciones por individuo:\n")
print(summary(admixture_data$ancestry_sum))

# No normalizamos automáticamente las proporciones:
# se conservan los valores publicados.


# 5. Definir colores para K = 4 -------------------------------

ancestry_colors <- c(
  Admixture1 = "#FF8C00",
  Admixture2 = "#FFD700",
  Admixture3 = "#3366CC",
  Admixture4 = "#109618"
)

ancestry_labels <- c(
  Admixture1 = "Grupo 1",
  Admixture2 = "Grupo 2",
  Admixture3 = "Grupo 3",
  Admixture4 = "Grupo 4"
)

# Formato del eje Y
proportion_labels <- function(x) {
  paste0(round(100 * x), "%")
}


# ============================================================
# FIGURA 1. ANCESTRÍA POR RANGO NATIVO / NO NATIVO
# ============================================================

# Ordenar primero por rango, después por región e individuo
range_order <- c("Native", "Non-native")

observed_ranges <- unique(admixture_data$Range)

range_levels <- c(
  range_order[range_order %in% observed_ranges],
  sort(setdiff(observed_ranges, range_order))
)

admixture_range <- admixture_data %>%
  mutate(
    Range = factor(Range, levels = range_levels)
  ) %>%
  arrange(Range, Region_stat, Ind_name_genind) %>%
  mutate(
    individual_order = row_number()
  )

# Formato largo: una fila por individuo y componente
range_long <- admixture_range %>%
  select(
    Ind_name_genind,
    Range,
    individual_order,
    all_of(ancestry_columns)
  ) %>%
  pivot_longer(
    cols = all_of(ancestry_columns),
    names_to = "Ancestry",
    values_to = "Proportion"
  )

p_range <- ggplot(
  range_long,
  aes(
    x = individual_order,
    y = Proportion,
    fill = Ancestry
  )
) +
  geom_col(width = 1) +
  facet_grid(
    ~ Range,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = ancestry_colors,
    labels = ancestry_labels,
    name = "Ancestría"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = proportion_labels,
    expand = c(0, 0)
  ) +
  labs(
    title = "Ancestría genética de Erigeron canadensis",
    subtitle = "K = 4; individuos del rango nativo y no nativo",
    x = "Individuos ordenados por región",
    y = "Coeficiente de ancestría"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(
      fill = "grey90",
      color = NA
    ),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p_range)

ggsave(
  filename = file.path(
    output_dir,
    "admixture_K4_native_vs_non_native.png"
  ),
  plot = p_range,
  width = 12,
  height = 6,
  dpi = 300
)


# ============================================================
# FIGURA 2. ANCESTRÍA POR REGIÓN GEOGRÁFICA
# ============================================================

# Orden alfabético de las regiones
# Orden geográfico de las regiones según la figura publicada
region_levels <- c(
  "South-Western North America",
  "Pacific Coast",
  "Intermountain US",
  "Central US",
  "North US Prairies",
  "Great Lakes",
  "South-Eastern US",
  "East Coast",
  "Western Mediterranean",
  "Western Europe",
  "Central Europe",
  "Baltic Europe",
  "Eastern Mediterranean",
  "Eastern Europe",
  "Armenia-Iran",
  "Middle Asia",
  "Siberia",
  "Far East Asia",
  "China"
)


admixture_region <- admixture_data %>%
  mutate(
    Region_stat = factor(
      Region_stat,
      levels = region_levels
    )
  ) %>%
  arrange(Region_stat, Ind_name_genind) %>%
  mutate(
    individual_order = row_number()
  )

# Posición central de cada región para las etiquetas del eje X
region_positions <- admixture_region %>%
  group_by(Region_stat) %>%
  summarise(
    midpoint = mean(individual_order),
    n_individuals = n(),
    .groups = "drop"
  )

# Límites entre regiones
region_boundaries <- region_positions %>%
  mutate(
    boundary = cumsum(n_individuals) + 0.5
  ) %>%
  filter(row_number() < n())

# Convertir a formato largo
region_long <- admixture_region %>%
  select(
    Ind_name_genind,
    Region_stat,
    individual_order,
    all_of(ancestry_columns)
  ) %>%
  pivot_longer(
    cols = all_of(ancestry_columns),
    names_to = "Ancestry",
    values_to = "Proportion"
  )

p_region <- ggplot(
  region_long,
  aes(
    x = individual_order,
    y = Proportion,
    fill = Ancestry
  )
) +
  geom_col(width = 1) +
  
  # Separadores blancos entre regiones
  geom_vline(
    data = region_boundaries,
    aes(xintercept = boundary),
    inherit.aes = FALSE,
    color = "white",
    linewidth = 0.5
  ) +
  
  scale_fill_manual(
    values = ancestry_colors,
    labels = ancestry_labels,
    name = "Ancestría"
  ) +
  
  scale_x_continuous(
    breaks = region_positions$midpoint,
    labels = as.character(region_positions$Region_stat),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = proportion_labels,
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Ancestría genética de Erigeron canadensis",
    subtitle = "K = 4; individuos ordenados por región geográfica",
    x = "Región geográfica",
    y = "Coeficiente de ancestría"
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8
    ),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

print(p_region)

ggsave(
  filename = file.path(
    output_dir,
    "admixture_K4_by_region.png"
  ),
  plot = p_region,
  width = 14,
  height = 7,
  dpi = 300
)


# ============================================================
# 6. Exportar tabla de individuos ordenados -------------------
# ============================================================

write.csv(
  admixture_range,
  file.path(
    output_dir,
    "admixture_K4_individuals_by_range.csv"
  ),
  row.names = FALSE
)

write.csv(
  admixture_region,
  file.path(
    output_dir,
    "admixture_K4_individuals_by_region.csv"
  ),
  row.names = FALSE
)

cat("\nAnálisis terminado.\n")
cat("Individuos incluidos:", nrow(admixture_data), "\n")
cat("Regiones representadas:", n_distinct(admixture_data$Region_stat), "\n")
cat("Resultados guardados en:", output_dir, "\n")

# ============================================================
# FIGURA 3. MAPA MUNDIAL DE ANCESTRÍA POR POBLACIÓN
# Adaptado de Lucas et al.
# ============================================================

# 1. Preparar datos geográficos y extraer población
map_data_individual <- metadata %>%
  mutate(
    Population = str_extract(
      Ind_name_genind,
      "(?<=^\\d{3})[A-Z]+(?=_)"
    )
  ) %>%
  select(
    Population,
    Latitude,
    Longitude,
    all_of(ancestry_columns)
  )

# Validar coordenadas y proporciones
if (anyNA(map_data_individual$Population)) {
  stop("No se pudieron extraer algunos códigos de población.")
}

if (anyNA(map_data_individual$Latitude) ||
    anyNA(map_data_individual$Longitude)) {
  warning("Hay muestras sin coordenadas; se excluirán de los mapas.")
}

map_data_individual <- map_data_individual %>%
  filter(
    is.finite(Latitude),
    is.finite(Longitude),
    if_all(all_of(ancestry_columns), is.finite)
  )

# 2. Calcular centro geográfico y ancestría media por población
data_population <- map_data_individual %>%
  group_by(Population) %>%
  summarise(
    latitude = mean(Latitude, na.rm = TRUE),
    longitude = mean(Longitude, na.rm = TRUE),
    across(
      all_of(ancestry_columns),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_individuals = n(),
    .groups = "drop"
  )

# Comprobar que las proporciones medias sumen aproximadamente 1
data_population <- data_population %>%
  mutate(
    ancestry_sum = rowSums(
      across(all_of(ancestry_columns))
    )
  )

cat("\nResumen de la suma de ancestría por población:\n")
print(summary(data_population$ancestry_sum))

# 3. Mapa base mundial
world <- map_data("world")

p_world_base <- ggplot() +
  geom_map(
    data = world %>% filter(lat > -60),
    map = world,
    aes(map_id = region),
    color = "#999999",
    fill = "antiquewhite",
    linewidth = 0.2
  ) +
  scale_x_continuous(
    limits = c(-155, 155),
    breaks = seq(-160, 160, 40),
    labels = c(
      "160°W", "120°W", "80°W", "40°W",
      "0°", "40°E", "80°E", "120°E", "160°E"
    )
  ) +
  scale_y_continuous(
    limits = c(10, 80),
    breaks = seq(10, 80, 10),
    labels = paste0(seq(10, 80, 10), "°N")
  ) +
  coord_equal(expand = FALSE) +
  labs(
    title = "Ancestría genética de Erigeron canadensis",
    subtitle = "Composición promedio por población (K = 4)",
    x = "Longitud",
    y = "Latitud",
    fill = "Ancestría"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# 4. Añadir gráficos circulares
p_map_world <- p_world_base +
  scatterpie::geom_scatterpie(
    data = data_population,
    aes(
      x = longitude,
      y = latitude,
      group = Population,
      r = 1.1
    ),
    cols = ancestry_columns,
    color = NA,
    alpha = 1
  ) +
  scale_fill_manual(
    values = ancestry_colors,
    labels = ancestry_labels,
    breaks = ancestry_columns,
    name = "Ancestría"
  )

print(p_map_world)

ggsave(
  file.path(output_dir, "admixture_K4_world_map.png"),
  p_map_world,
  width = 13,
  height = 7,
  dpi = 300
)


# ============================================================
# FIGURA 4. MAPA DE NORTEAMÉRICA
# ============================================================

p_north_america_base <- ggplot() +
  geom_map(
    data = world,
    map = world,
    aes(map_id = region),
    color = "#999999",
    fill = "antiquewhite",
    linewidth = 0.2
  ) +
  scale_x_continuous(
    limits = c(-160, -40),
    breaks = seq(-160, -40, 20),
    labels = c(
      "160°W", "140°W", "120°W", "100°W",
      "80°W", "60°W", "40°W"
    )
  ) +
  scale_y_continuous(
    limits = c(5, 80),
    breaks = seq(10, 80, 10),
    labels = paste0(seq(10, 80, 10), "°N")
  ) +
  coord_equal(expand = FALSE) +
  labs(
    title = "Ancestría genética en Norteamérica",
    subtitle = "Composición promedio por población (K = 4)",
    x = "Longitud",
    y = "Latitud",
    fill = "Ancestría"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_map_north_america <- p_north_america_base +
  scatterpie::geom_scatterpie(
    data = data_population,
    aes(
      x = longitude,
      y = latitude,
      group = Population,
      r = 1.15
    ),
    cols = ancestry_columns,
    color = NA,
    alpha = 1
  ) +
  scale_fill_manual(
    values = ancestry_colors,
    labels = ancestry_labels,
    breaks = ancestry_columns,
    name = "Ancestría"
  )

print(p_map_north_america)

ggsave(
  file.path(output_dir, "admixture_K4_north_america_map.png"),
  p_map_north_america,
  width = 10,
  height = 7,
  dpi = 300
)


# ============================================================
# MAPA DE ANCESTRÍA GENÉTICA: RANGO NO NATIVO
# ============================================================

non_native_base <- ggplot() +
  geom_map(
    data = world,
    map = world,
    aes(map_id = region),
    color = "#999999",
    fill = "antiquewhite",
    linewidth = 0.2
  ) +
  scale_x_continuous(
    limits = c(-15, 150),
    breaks = seq(0, 150, 30),
    labels = c("0°", "30°E", "60°E", "90°E", "120°E", "150°E")
  ) +
  scale_y_continuous(
    limits = c(20, 75),
    breaks = seq(20, 70, 10),
    labels = paste0(seq(20, 70, 10), "°N")
  ) +
  coord_equal(expand = FALSE) +
  labs(
    title = "Ancestría genética en el rango no nativo",
    subtitle = "Composición promedio por población (K = 4)",
    x = "Longitud",
    y = "Latitud",
    fill = "Ancestría"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_map_non_native <- non_native_base +
  scatterpie::geom_scatterpie(
    data = data_population,
    aes(
      x = longitude,
      y = latitude,
      group = Population,
      r = 1.15
    ),
    cols = ancestry_columns,
    color = NA,
    alpha = 1
  ) +
  scale_fill_manual(
    values = ancestry_colors,
    labels = ancestry_labels,
    breaks = ancestry_columns,
    name = "Ancestría"
  )

print(p_map_non_native)

ggsave(
  file.path(output_dir, "admixture_K4_non_native_map.png"),
  p_map_non_native,
  width = 10,
  height = 7,
  dpi = 300
)









