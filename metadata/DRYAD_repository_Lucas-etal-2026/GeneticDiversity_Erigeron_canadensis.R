################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# Genetic diversity analysis
################################################################################

#####----------------------------------------------------------------------#####

packages <- c("vcfR", "hierfstat", "adegenet", "poppr", "tidyr", "dplyr")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

#####----------------------------------------------------------------------#####

vcf <- read.vcfR("Erigeron.canadensis.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf")

genind_data <- vcfR2genind(vcf)

pop_data <- read.csv("metadata.csv")
pop_data <- pop_data[match(indNames(genind_data), pop_data$Individual), ]

genind_data@pop <- as.factor(pop_data$Population)
genind_data@other$Region_stat <- as.factor(pop_data$Region_stat)

#####----------------------------------------------------------------------#####

basic_stats <- basic.stats(genind_data)

Ho <- rowMeans(basic_stats$Ho, na.rm = TRUE)
He <- rowMeans(basic_stats$Hs, na.rm = TRUE)
Fis <- rowMeans(basic_stats$Fis, na.rm = TRUE)

allelic_richness <- allelic.richness(genind_data)
Ar <- rowMeans(allelic_richness$Ar, na.rm = TRUE)

diversity_df <- data.frame(
  Individual = indNames(genind_data),
  Region_stat = pop_data$Region_stat,
  Ho = Ho,
  He = He,
  Fis = Fis,
  Ar = Ar
)


#####----------------------------------------------------------------------#####

perform_ttest <- function(trait) {
  combinations <- combn(unique(diversity_df$Region_stat), 2, simplify = FALSE)
  
  ttest_results <- lapply(combinations, function(pair) {
    group1 <- diversity_df %>% filter(Region_stat == pair[1]) %>% pull(!!sym(trait))
    group2 <- diversity_df %>% filter(Region_stat == pair[2]) %>% pull(!!sym(trait))
    test <- t.test(group1, group2)
    
    data.frame(
      Trait = trait,
      Group1 = pair[1],
      Group2 = pair[2],
      p.value = test$p.value
    )
  })
  
  ttest_results <- do.call(rbind, ttest_results)
  ttest_results$p.adjusted <- p.adjust(ttest_results$p.value, method = "bonferroni")
  return(ttest_results)
}

ttest_Ho <- perform_ttest("Ho")
ttest_He <- perform_ttest("He")
ttest_Fis <- perform_ttest("Fis")
ttest_Ar <- perform_ttest("Ar")

all_ttests <- bind_rows(ttest_Ho, ttest_He, ttest_Fis, ttest_Ar)

#####----------------------------------------------------------------------#####
