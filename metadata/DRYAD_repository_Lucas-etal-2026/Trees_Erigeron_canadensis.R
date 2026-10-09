################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# SplitsTree analysis using Nei distances
################################################################################

#####----------------------------------------------------------------------#####

# Load required packages
packages <- c("vcfR", "adegenet", "poppr", "dplyr", "phangorn", "ggplot2")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


#####----------------------------------------------------------------------#####

load("genind_file.RData")  
my.gl <- New_genind_conyza

meta.dat <- read.csv("metadata_file.csv")

identical(my.gl@ind.names, meta.dat$Ind_name_genind)
pop(my.gl) <- meta.dat$Ind_Popcode

my.dist <- dartR::gl.dist.pop(my.gl, method = "nei", plot.out = FALSE, verbose = 5)

phangorn::writeDist(my.dist, file = "Conyza.dist.nex", format = "nexus")

#####----------------------------------------------------------------------#####

netout <- phangorn::read.nexus.splits("Conyza.dist.out.nex")
netout <- phangorn::as.networx(netout)

#####----------------------------------------------------------------------#####

admixture_colors <- c("darkorange", "gold", "#3366cc", "#109618")
netout$edge.length <- netout$edge.length * 0.7 

plot(netout, type = "equal angle", use.edge.length = TRUE,
     show.tip.label = FALSE, show.edge.label = FALSE,
     show.node.label = FALSE, show.nodes = FALSE,
     edge.width = 0.7)

title("SplitsTree on Nei-distances")

net.cols <- meta.dat$Range_col[match(netout$tip.label, meta.dat$Ind_Popcode)]

ape::tiplabels(pie = meta.dat[match(netout$tip.label, meta.dat$Ind_Popcode),
                              c("Admixture1", "Admixture2", "Admixture3", "Admixture4")], 
               piecol = admixture_colors, cex = 0.7)

ape::tiplabels(pch = 21, bg = net.cols, cex = 1)

legend("topright",
       title = "Range",
       legend = c("Native", "Non-Native"),
       pch = 21,
       pt.bg = c("dark blue", "red3"),
       pt.cex = 1.8,
       lty = 0,
       cex = 0.7)

legend("bottomright",
       title = "Admixture Cluster",
       legend = c("Admixture1", "Admixture2", "Admixture3", "Admixture4"),
       pch = 21,
       pt.bg = admixture_colors,
       pt.cex = 1.8,
       lty = 0,
       cex = 0.7)

#####---------------------------------------------------------------------------#####
