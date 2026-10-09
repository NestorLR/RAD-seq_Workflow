################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# Offset analysis
################################################################################

#####----------------------------------------------------------------------#####

packages <- c("readxl", "tidyverse", "lme4", "lmerTest")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


#####----------------------------------------------------------------------#####

# setwd("your/path/here")

# Load offset data

model <- lm(log(mean_shoot_biomass) ~ mean_offset, data = data)
summary(model)

ggplot(data, aes(mean_offset, log(mean_shoot_biomass))) +
  geom_point(cex = 3, alpha = 0.3) +
  geom_smooth(method = "lm",
              se = TRUE, 
              color = "black",
              linetype = "solid",
              linewidth = 1.5) +
  labs(
    x = "Genomic offset (mean)",
    y = "Log(mean shoot biomass)",
    title = "Shoot biomass vs genomic offset"
  ) +
  theme_minimal()


#####----------------------------------------------------------------------#####

hist(data$mean_relative_field_performance)
hist(log(data$mean_relative_field_performance))

model <- lm(log(mean_relative_field_performance) ~ mean_offset, data = data)
summary(model)

ggplot(data, aes(mean_offset, log(mean_relative_field_performance))) +
  geom_point(cex = 3, alpha = 0.3) +
  geom_smooth(method = "lm",
              se = TRUE, 
              color = "black",
              linetype = "dashed",
              linewidth = 1.5) +
  labs(
    x = "Genomic offset (mean)",
    y = "Log(mean relative field performance)",
    title = "Relative performance vs genomic offset"
  ) +
  theme_minimal()


#####---------------------------------------------------------------------------#####
