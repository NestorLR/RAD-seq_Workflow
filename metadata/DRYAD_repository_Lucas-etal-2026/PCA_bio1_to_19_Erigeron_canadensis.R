################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# PCA analysis (Bio01 to Bio19 variables)
################################################################################

#####----------------------------------------------------------------------#####

packages <- c("tidyverse", "factoextra", "ggplot2")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

sessionInfo()

#####----------------------------------------------------------------------#####

data <- read.table("bioclim_file.txt", header = TRUE, sep = "\t") 

bio_data <- data %>%
  select(Bio01:Bio19)

pca_result <- prcomp(bio_data, scale. = TRUE)


#####----------------------------------------------------------------------#####

scree_plot <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))

pca_ind_plot <- fviz_pca_ind(pca_result,
                             geom.ind = "point",
                             col.ind = "cos2",
                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                             repel = TRUE) +
  labs(x = "PC1", y = "PC2")

pca_var_plot <- fviz_pca_var(pca_result,
                             col.var = "contrib",
                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                             repel = TRUE) +
  labs(x = "PC1", y = "PC2")

#####---------------------------------------------------------------------------#####
