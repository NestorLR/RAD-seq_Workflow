################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# Admixture analysis
################################################################################

#####----------------------------------------------------------------------#####


# Load required packages
packages <- c("dplyr", "ggplot2", "tidyr", "scatterpie")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


#####----------------------------------------------------------------------#####

# Load population-level admixture and geographic data
data <- read.csv("your_data_file.csv")  
str(data)

data_population <- data %>%
  group_by(Pop) %>%
  summarise(
    latitude = mean(latitude, na.rm = TRUE),
    longitude = mean(longitude, na.rm = TRUE),
    Admixture1 = mean(Admixture1, na.rm = TRUE),
    Admixture2 = mean(Admixture2, na.rm = TRUE),
    Admixture3 = mean(Admixture3, na.rm = TRUE),
    Admixture4 = mean(Admixture4, na.rm = TRUE)
  ) %>%
  ungroup()

str(data_population)

#####----------------------------------------------------------------------#####

world <- map_data("world")

pop_map_plot <- ggplot() +
  geom_map(data = world[world$lat > -60, ],
           map = world,
           aes(map_id = region),
           color = "#999999",
           fill = "antiquewhite",
           size = 0.1) +
  scale_y_continuous(limits = c(10, 80),
                     breaks = c(0, 20, 40, 60, 80),
                     labels = c("0°N", "20°N", "40°N", "60°N", "80°N")) +
  scale_x_continuous(limits = c(-155, 155),
                     breaks = seq(-160, 160, 40),
                     labels = c("-160°W", "-120°W", "-80°W", "-40°W", "0°",
                                "40°E", "80°E", "120°E", "160°E")) +
  theme(panel.background = element_rect(fill = "aliceblue", colour = "#000000"),
        plot.background = element_rect(fill = "#ffffff", colour = "#ffffff"),
        legend.position = "bottom",
        axis.title.y = element_text(angle = 90),
        axis.ticks = element_line()) +
  geom_hline(aes(yintercept = 10), size = 0.2) +
  geom_hline(aes(yintercept = 80), size = 0.2) +
  geom_vline(aes(xintercept = -155), size = 0.2) +
  geom_vline(aes(xintercept = 155), size = 0.2) +
  coord_equal(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude")

pop_map_plot

#####----------------------------------------------------------------------#####

pie_map_plot <- pop_map_plot +
  geom_scatterpie(aes(x = longitude, y = latitude, group = Pop, r = 1.1),
                  cols = c("Admixture1", "Admixture2", "Admixture3", "Admixture4"),
                  color = NA,
                  alpha = 1,
                  data = data_population) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("darkorange", "gold", "#3366cc", "#109618", "#990099", "#0099c6"))

pie_map_plot

#####----------------------------------------------------------------------#####


data_population <- data %>%
  group_by(Pop) %>%
  summarise(
    latitude = mean(latitude, na.rm = TRUE),
    longitude = mean(longitude, na.rm = TRUE),
    Admixture1 = mean(Admixture1, na.rm = TRUE),
    Admixture2 = mean(Admixture2, na.rm = TRUE),
    Admixture3 = mean(Admixture3, na.rm = TRUE),
    Admixture4 = mean(Admixture4, na.rm = TRUE)
  ) %>%
  ungroup()

str(data_population)

#####----------------------------------------------------------------------#####


world <- map_data("world")

pop_map_plot <- ggplot() +
  geom_map(data = world[world$lat > -60, ],
           map = world,
           aes(map_id = region),
           color = "#999999",
           fill = "antiquewhite",
           size = 0.1) +
  scale_y_continuous(limits = c(5, 80),
                     breaks = c(10, 20, 30, 40, 50, 60, 70, 80),
                     labels = c("10°N", "20°N", "30°N", "40°N", "50°N", "60°N", "70°N", "80°N")) +
  scale_x_continuous(limits = c(-160, -40),
                     breaks = seq(-160, -40, 20),
                     labels = c("160°W", "140°W", "120°W", "100°W", "80°W", "60°W", "40°W")) +
  theme(panel.background = element_rect(fill = "aliceblue", colour = "#000000"),
        plot.background = element_rect(fill = "#ffffff", colour = "#ffffff"),
        legend.position = "bottom",
        axis.title.y = element_text(angle = 90),
        axis.ticks = element_line()) +
  geom_hline(aes(yintercept = 0), size = 0.2) +
  geom_hline(aes(yintercept = 80), size = 0.2) +
  geom_vline(aes(xintercept = 0), size = 0.2) +
  geom_vline(aes(xintercept = 155), size = 0.2) +
  coord_equal(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude")

pop_map_plot

#####----------------------------------------------------------------------#####

pie_map_plot <- pop_map_plot +
  geom_scatterpie(aes(x = longitude, y = latitude, group = Pop, r = 1.15),
                  cols = c("Admixture1", "Admixture2", "Admixture3", "Admixture4"),
                  color = NA,
                  alpha = 1,
                  data = data_population) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("darkorange", "gold", "#3366cc", "#109618", "#990099", "#0099c6"))

pie_map_plot

#####----------------------------------------------------------------------#####
