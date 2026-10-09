################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# divMigrate analysis
################################################################################

#####----------------------------------------------------------------------#####

packages <- c("htmltools", "vcfR", "adegenet", "diveRsity", "devtools", 
              "VariantAnnotation", "pegas", "dartR", "divMigrate", "ggplot2", 
              "dplyr", "maps")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

sessionInfo()

#####----------------------------------------------------------------------#####

vcf <- read.vcfR("Erigeron.canadensis.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf")  
genind_conyza <- vcfR2genind(vcf)

pop_data <- read.csv("your_population_file.csv")  

pop_data_ordered <- pop_data[match(indNames(genind_conyza), pop_data$IndividualID), ]
populations <- as.factor(pop_data_ordered$Population)
pop(genind_conyza) <- populations

genpop_conyza <- genind2genpop(genind_conyza)




#####----------------------------------------------------------------------#####

log_file <- readLines("your_log_file.log")  

start <- grep("\\$nmRelMig", log_file)[1] + 1
end <- grep("\\$nmRelMigSig", log_file)[1] - 1

nmRelMig_lines <- log_file[start:end]
nmRelMig_lines <- gsub("^\\s+", "", nmRelMig_lines)

migration_values <- as.numeric(unlist(strsplit(paste(nmRelMig_lines, collapse = " "), "\\s+")))
migration_matrix <- matrix(migration_values, nrow = 19, ncol = 19, byrow = TRUE)

#####----------------------------------------------------------------------#####

edges <- data.frame(
  from = rep(populations, each = length(populations)),
  to = rep(populations, times = length(populations)),
  value = as.vector(migration_matrix)
) %>%
  filter(from != to, !is.na(value), value > 0.5)

edges_coords <- edges %>%
  left_join(coords, by = c("from" = "Population")) %>%
  rename(from_long = Longitude, from_lat = Latitude) %>%
  left_join(coords, by = c("to" = "Population")) %>%
  rename(to_long = Longitude, to_lat = Latitude)

#####----------------------------------------------------------------------#####

world <- map_data("world")

pop_map_plot <- ggplot() +
  geom_rect(aes(xmin = -155, xmax = 155, ymin = 10, ymax = 80), fill = "aliceblue", color = NA) +
  geom_map(data = world[world$lat > -60, ],
           map = world, aes(map_id = region),
           color = "#999999", fill = "antiquewhite", size = 0.1) +
  scale_y_continuous(limits = c(10, 80),
                     breaks = c(0, 20, 40, 60, 80),
                     labels = c("0°N", "20°N", "40°N", "60°N", "80°N")) +
  scale_x_continuous(limits = c(-155, 155),
                     breaks = seq(-160, 160, 40),
                     labels = c("-160°W", "-120°W", "-80°W", "-40°W", "0°",
                                "40°E", "80°E", "120°E", "160°E")) +
  theme(
    panel.background = element_rect(fill = "aliceblue", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(color = "white", size = 0.1),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  geom_curve(data = edges_coords, aes(x = from_long, y = from_lat, 
                                      xend = to_long, yend = to_lat, color = value),
             curvature = 0.25, 
             arrow = arrow(length = unit(0.1, "inches")), size = 0.7) +
  geom_point(data = coords, aes(x = Longitude, y = Latitude),
             color = "black", fill = "lightgray", shape = 21, size = 8, alpha = 0.3) +
  scale_color_gradient(low = "#59C173", high = "#662D8C", name = "Migration\nIntensity",
                       guide = guide_colorbar(barwidth = 0.4, barheight = 6)) +
  labs(x = "Longitude", y = "Latitude") +
  coord_fixed()

pop_map_plot


#####---------------------------------------------------------------------------#####
