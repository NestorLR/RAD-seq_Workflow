#################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2025
# Check cpDNA Locus Effect
################################################################################

#####----------------------------------------------------------------------#####

packages <- c("adegenet", "vcfR", "dartR", "adegraphics", "ggmap", "scatterpie",
              "ggplot2", "gridExtra", "gplots", "qvalue", "LEA", "remotes", "sf",
              "tidyverse", "StAMPP", "geodist", "segmented", "readr", "expss")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


#####----------------------------------------------------------------------#####

add.LocInfo2gl <- function(gl, Megablast.file) {
  if (is.null(gl@chromosome)) {
    chr.names <- locNames(gl)
    for (i in 1:length(chr.names)) {
      chr.names[i] <- substr(chr.names[i], 1, unlist(gregexpr("_", chr.names[i]))[3] - 1)
    }
    gl@chromosome <- as.factor(chr.names)
  }
  
  if (file.exists(Megablast.file)) {
    d <- as.data.frame(read_tsv(Megablast.file, col_names =
                                  c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                                    "staxids", "sscinames", "sblastnames", "sskingdoms", "salltitles")))
    
    Megablast <- ifelse(
      grepl("chloroplast", d$salltitles, fixed = TRUE) | 
        grepl("plastid", d$salltitles, fixed = TRUE),
      "C", "U"
    )
    
    Megablast <- vlookup(lookup_value = gl@chromosome,
                         dict = as.data.frame(cbind(d$qseqid, Megablast)), result_column = 2)
    
    Megablast[is.na(Megablast)] <- "U"
    gl@other$LocInfo.Megablast <- Megablast
    gl@other$LocInfo.Genome <- Megablast
  } else {
    message(Megablast.file, " not found!")
  }
  
  return(gl)
}

#####----------------------------------------------------------------------#####


load(file = paste0(dataPath, "Genind.RData"))
d.gl <- New_genind_conyza

d.gl <- add.LocInfo2gl(d.gl, Megablast.file = paste0(dataPath, "C.megablast.out"))


d.gl.pc <- glPca(d.gl, nf = 100)
plot(d.gl.pc$eig)
plot(d.gl.pc$scores[, 2] ~ d.gl.pc$scores[, 1])
SNP.clusters <- find.clusters(d.gl, n.pca = 100)

d.gl.cp <- gl.drop.loc(d.gl, loc.list = locNames(d.gl)[d.gl@other$LocInfo.Genome != "C"])
d.gl.ocp <- gl.drop.loc(d.gl, loc.list = locNames(d.gl)[d.gl@other$LocInfo.Genome == "C"])

d.gl.cp.pc <- glPca(d.gl.cp, nf = 100)
d.gl.ocp.pc <- glPca(d.gl.ocp, nf = 100)

plot(d.gl.ocp.pc$eig)
plot(d.gl.ocp.pc$scores[, 2] ~ d.gl.ocp.pc$scores[, 1], pch = 16, col = as.numeric(SNP.clusters$grp))

plot(d.gl.cp.pc$eig)
plot(d.gl.cp.pc$scores[, 2] ~ d.gl.cp.pc$scores[, 1], pch = 16, col = as.numeric(SNP.clusters$grp))

plot(d.gl.cp.pc$scores[, 1] ~ d.gl.pc$scores[, 1])

png(filename = paste0(dataPath, "cpDNA_SNPs_PCA_kmeans4.png"), width = 300, height = 300, units = 'mm', res = 600)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 0.5))
plot(d.gl.pc$scores[, 2] ~ d.gl.pc$scores[, 1], main = paste("All loci, N=", nLoc(d.gl)), pch = 16, col = as.numeric(SNP.clusters$grp))
plot(-1 * d.gl.cp.pc$scores[, 2] ~ d.gl.cp.pc$scores[, 1], main = paste("cpDNA loci, N=", nLoc(d.gl.cp)), pch = 16, col = as.numeric(SNP.clusters$grp))
plot(d.gl.ocp.pc$scores[, 2] ~ d.gl.ocp.pc$scores[, 1], pch = 16, col = as.numeric(SNP.clusters$grp), main = paste("without cpDNA loci, N=", nLoc(d.gl.ocp)))
plot(d.gl.cp.pc$scores[, 1] ~ d.gl.ocp.pc$scores[, 1], main = "PCA1 cpDNA vs without cpDNA")

#####----------------------------------------------------------------------#####

