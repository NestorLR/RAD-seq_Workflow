################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# RDA
################################################################################

#####----------------------------------------------------------------------#####


packages <- c("vegan", "dplyr", "ggplot2", "adegenet", "ggsci", "cowplot",
              "tibble", "scatterpie", "ggforce")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


#####----------------------------------------------------------------------#####

load("your_genind_file.RData") 

freq <- as.data.frame(New_genind_conyza)
freq <- freq %>% select_if(~ !any(is.na(.)))  

env <- read.csv("your_metadata_file.csv")  

env$climatic_water_deficit[is.na(env$climatic_water_deficit)] <- mean(env$climatic_water_deficit, na.rm = TRUE)
env$New_Climatic_PC1[is.na(env$New_Climatic_PC1)] <- mean(env$New_Climatic_PC1, na.rm = TRUE)

rownames(freq) <- tolower(trimws(rownames(freq)))
env$Ind_name_genind <- tolower(trimws(env$Ind_name_genind))

common_ids <- intersect(rownames(freq), env$Ind_name_genind)
freq <- freq[common_ids, ]
env <- env[match(common_ids, env$Ind_name_genind), ]

#####----------------------------------------------------------------------#####

dbrda_all <- capscale(freq ~ climatic_water_deficit + Latitude + New_Climatic_PC1 + human_footprint,
                      data = env, distance = "bray")

env$CAP1 <- scores(dbrda_all, display = "sites")[, 1]
env$CAP2 <- scores(dbrda_all, display = "sites")[, 2]

bi_dat <- as.data.frame(dbrda_all$CCA$biplot)
colnames(bi_dat)[1:2] <- c("CAP1", "CAP2")
bi_dat$Variable <- rownames(bi_dat)

bi_dat <- bi_dat %>%
  mutate(Variable = case_when(
    Variable == "climatic_water_deficit" ~ "Climatic water deficit",
    Variable == "Latitude" ~ "Latitude",
    Variable == "New_Climatic_PC1" ~ "PC1_BioClim",
    Variable == "human_footprint" ~ "Human Footprint",
    TRUE ~ Variable
  ))

eig_constr <- dbrda_all$CCA$eig
total_constr <- sum(eig_constr)
total_inertia <- dbrda_all$tot.chi

cap1_fitted <- eig_constr[1] / total_constr * 100
cap1_total <- eig_constr[1] / total_inertia * 100
cap2_fitted <- eig_constr[2] / total_constr * 100
cap2_total <- eig_constr[2] / total_inertia * 100

#####----------------------------------------------------------------------#####

admixture_cols <- c("Admixture1", "Admixture2", "Admixture3", "Admixture4")

plot_faceted <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_scatterpie(data = env,
                  aes(x = CAP1, y = CAP2, group = Ind_name_genind, r = 0.04),
                  cols = admixture_cols) +
  geom_segment(data = bi_dat,
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(10, "pt")),
               linewidth = 1) +
  geom_label(data = bi_dat,
             aes(x = CAP1, y = CAP2, label = Variable),
             fontface = "bold",
             size = 3) +
  scale_fill_manual(values = c("#dc3912", "#ff9900", "#3366cc", "#109618")) +
  facet_grid(~Range, labeller = labeller(Range = c("Native" = "Native", "Non-native" = "Non-native"))) +
  labs(
    title = "Redundancy Analysis (dbRDA) by Range",
    x = paste0("CAP1 (", round(cap1_fitted, 1), "% fitted / ", round(cap1_total, 1), "% total)"),
    y = paste0("CAP2 (", round(cap2_fitted, 1), "% fitted / ", round(cap2_total, 1), "% total)")
  ) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "gray25"),
    text = element_text(size = 10),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )
print(plot_faceted)

#####---------------------------------------------------------------------------#####
