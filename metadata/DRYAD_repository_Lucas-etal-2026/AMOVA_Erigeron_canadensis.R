################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# AMOVA analysis
################################################################################

#####----------------------------------------------------------------------#####

# Load required packages
packages <- c("dartR", "hierfstat", "gaston", "vcfR", "reshape2", 
              "pegas", "cowplot", "adegenet", "PopGenReport", "poppr")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

sessionInfo()

#####----------------------------------------------------------------------#####

# setwd("your/path/here")

# Load VCF file
vcf <- read.vcfR("Erigeron.canadensis.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf")

# Convert to genlight object
new_Conyza_gl <- vcfR2genlight(vcf)
new_Conyza_gl

#####---------------------- Hierarchical AMOVA ----------------------------#####

pop.nr <- as.numeric(substr(indNames(new_Conyza_gl), 1, 3))
pop(new_Conyza_gl) <- as.factor(pop.nr)

status <- ifelse(pop.nr < 109, "native", "non-native")

hierarchy <- data.frame(
  status = status,
  pop = as.character(pop(new_Conyza_gl)),
  stringsAsFactors = FALSE
)

tapply(hierarchy$pop, hierarchy$status, function(x) length(unique(x)))

new_Conyza_gl@other$hierarchy <- ~status/pop
amova_result <- poppr.amova(new_Conyza_gl, ~status/pop)
print(amova_result)

amova_test <- randtest(amova_result, nrepet = 999)
print(amova_test)

new_Conyza_gl@other$hierarchy <- ~status
amova_result_simple <- poppr.amova(new_Conyza_gl, ~status)
print(amova_result_simple)

#####---------------------- Native / Non-native ---------------------------#####

pop.nr <- as.numeric(substr(indNames(new_Conyza_gl), 1, 3))
status <- ifelse(pop.nr < 109, "native", "non-native")
pop(new_Conyza_gl) <- as.factor(pop.nr)

hierarchy <- data.frame(
  status = status,
  pop = as.character(pop(new_Conyza_gl)),
  stringsAsFactors = FALSE
)
strata(new_Conyza_gl) <- hierarchy

natives_gl <- new_Conyza_gl[strata(new_Conyza_gl)$status == "native", ]
pop(natives_gl) <- as.factor(strata(natives_gl)$pop)

amova.nat <- poppr.amova(
  x = natives_gl,
  ~pop,
  missing = "loci",
  cutoff = 0.05,
  correction = "quasieuclid",
  method = "ade4"
)
amova.nat

nonnative_gl <- new_Conyza_gl[strata(new_Conyza_gl)$status == "non-native", ]
pop(nonnative_gl) <- as.factor(strata(nonnative_gl)$pop)

amova.inv <- poppr.amova(
  x = nonnative_gl,
  ~pop,
  missing = "loci",
  cutoff = 0.05,
  correction = "quasieuclid",
  method = "ade4"
)
amova.inv

#####----------------------------- Regional AMOVA ------------------------------#####


pop.nr <- as.numeric(substr(indNames(new_Conyza_gl), 1, 3))

category <- rep(NA, length(pop.nr))

category[pop.nr %in% 2:6 | pop.nr %in% 13:15 | pop.nr %in% 22:27] <- "South-Western North America"
category[pop.nr %in% 7:12 | pop.nr == 93] <- "Pacific Coast"
category[pop.nr %in% 16:20 | pop.nr %in% 28:32 | pop.nr %in% 90:92] <- "Intermountain US"
category[pop.nr == 21 | pop.nr %in% 33:41 | pop.nr %in% 45:55 | pop.nr %in% 62:64] <- "Central US"
category[pop.nr %in% 43:44 | pop.nr %in% 94:98] <- "North US Prairies"
category[pop.nr %in% 56:58 | pop.nr == 65 | pop.nr %in% 69:72 | pop.nr %in% 78:83 | pop.nr %in% 99:101] <- "Great Lakes"
category[pop.nr %in% 59:61 | pop.nr %in% 66:68 | pop.nr %in% 74:77] <- "South-Eastern US"
category[pop.nr %in% 78:80 | pop.nr %in% 84:89 | pop.nr %in% 102:109] <- "East Coast"
category[pop.nr %in% 111:118 | pop.nr %in% 128:130 | pop.nr %in% 136:139 | pop.nr %in% 148:150] <- "Western Mediterranean"
category[pop.nr %in% 119:123] <- "Western Europe"
category[pop.nr %in% 124:126 | pop.nr %in% 133:135 | pop.nr %in% 142:147 | pop.nr %in% 155:163 | pop.nr %in% 175:184] <- "Central Europe"
category[pop.nr %in% 131:132 | pop.nr %in% 140:141 | pop.nr %in% 151:154 | pop.nr %in% 167:174] <- "Baltic Europe"
category[pop.nr %in% 164:166 | pop.nr %in% 187:202] <- "Eastern Mediterranean"
category[pop.nr %in% 185:186 | pop.nr %in% 215:236] <- "Eastern Europe"
category[pop.nr %in% 203:208] <- "Armenia-Iran"
category[pop.nr %in% 209:214 | pop.nr %in% 253:264] <- "Middle Asia"
category[pop.nr %in% 237:245] <- "Siberia"
category[pop.nr %in% 246:252 | pop.nr %in% 277:294] <- "Far East Asia"
category[pop.nr %in% 265:276] <- "China"

category <- as.factor(category)

hierarchy <- data.frame(pop = category)
strata(new_Conyza_gl) <- hierarchy
new_Conyza_gl@other$hierarchy <- ~pop

amova_result <- poppr.amova(new_Conyza_gl, ~pop)
print(amova_result)

#####---------------------------------------------------------------------------#####


