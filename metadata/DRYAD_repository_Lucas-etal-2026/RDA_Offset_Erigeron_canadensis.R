################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# RDA-based Genomic Offset Analysis
################################################################################

#####----------------------------------------------------------------------#####

# Load required packages
packages <- c("SNPRelate", "vegan", "vcfR", "gdsfmt", "adegenet", "data.table", 
              "dplyr", "ggplot2", "ggsci", "cowplot", "tibble", "scatterpie", "maps")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

#####----------------------------------------------------------------------#####

data <- read.csv("metadata.csv", header = TRUE)

native <- data[data$Range == "Native", ]
nonnative <- data[data$Range == "Non-native", ]

load(file = "GenindObject.RData")
freq <- as.data.frame(New_genind_conyza)

freq <- freq %>% select_if(~ !any(is.na(.)))

abiotic <- data[, c("climatic_water_deficit", "New_Climatic_PC1", "Latitude", "Human_Foot_Print", "Range")]
abiotic_imputed <- abiotic %>% mutate_if(is.numeric, ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))

ord <- rda(freq[2:237, ] ~ ., data = abiotic_imputed[2:237, ], scale = TRUE)
predzO <- fitted(ord)

diz <- sapply(1:nrow(predzO), function(i) sqrt(sum((predzO[i, ] - freq[2:237, ][i, ])^2)))
write.csv(diz, file = "diz_nativeRange.csv", row.names = FALSE)

abiotic_imputed_subset <- abiotic_imputed[238:641, ] %>% mutate_if(is.numeric, ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))
predz <- predict(ord, newdata = abiotic_imputed_subset, type = 'response')
dizTest <- sapply(1:nrow(predz), function(i) sqrt(sum((predz[i, ] - freq[238:641, ][i, ])^2)))
write.csv(dizTest, file = "diz_nonnativeRange.csv", row.names = FALSE)

#####----------------------------------------------------------------------#####

world <- map_data("world")

map <- ggplot() +
  geom_map(data = world[world$lat > -60, ], map = world,
           aes(map_id = region), color = "#999999", fill = "antiquewhite", size = 0.1) +
  scale_y_continuous(limits = c(10, 80), breaks = c(0, 20, 40, 60, 80), labels = c("0°N", "20°N", "40°N", "60°N", "80°N")) +
  scale_x_continuous(limits = c(-155, 155), breaks = seq(-160, 160, 40),
                     labels = c("-160°W", "-120°W", "-80°W", "-40°W", "0°", "40°E", "80°E", "120°E", "160°E")) +
  theme(panel.background = element_rect(fill = "aliceblue", colour = "#000000"),
        plot.background = element_rect(fill = "#ffffff", colour = "#ffffff"),
        legend.position = "right",
        axis.title.y = element_text(angle = 90),
        axis.ticks = element_line()) +
  coord_equal(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  geom_point(data = data, aes(x = Longitude, y = Latitude, fill = Offset, color = Offset),
             size = 2, alpha = 0.8, shape = 21, stroke = 0.5) +
  scale_fill_gradient(low = "white", high = "black", name = "Offset",
                      guide = guide_colorbar(barwidth = 0.5, barheight = 10)) +
  scale_color_gradient(low = "white", high = "black", guide = "none")

print(map)

#####----------------------------------------------------------------------#####

set.seed(xxx)
permute_offsets <- function(data, n_perm = 999) {
  replicate(n_perm, mean(sample(data$Offset), na.rm = TRUE))
}

perm_results <- permute_offsets(data)

hist(perm_results, col = "gray", main = "Permutation Test for Offset", 
     xlab = "Offset Values", breaks = 30)
abline(v = mean(data$Offset, na.rm = TRUE), col = "red", lwd = 2)

p_value <- mean(perm_results >= mean(data$Offset, na.rm = TRUE))
cat("Permutation test p-value:", p_value, "\n")

#####----------------------------------------------------------------------#####

