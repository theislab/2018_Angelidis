# R code
# Miko Simon
# Wed Sept 19th 2018
# This code will generate a scatter plot for aging foldchanges

# Load R libraries ####
library(readxl)
library(mcr)

# Load protein data ####
prot <- read_excel("../data/Table S2_DGE_bulk.xlsx", sheet = 'Proteome')
prot_fc <- prot$`Student's T-test Difference_old/young [log2]`
names(prot_fc) <- prot$`Gene names`

# Load whole lung bulk data ####
bulk <- read_excel("../data/Table S2_DGE_bulk.xlsx", sheet = 'Bulk transcriptome')
bulk_fc <- bulk$`log2FoldChange [old/young]`
names(bulk_fc) <- bulk$`gene name`

# Load in silico bulk data ####
insilico <- read_excel("../data/Table S2_DGE_bulk.xlsx", sheet = 'In silico transcriptome')
insilico_fc <- insilico$`log2FoldChange [old/young]`
names(insilico_fc) <- insilico$`Gene name`

# generate scatter plots ####
genScatterPlot <- function(fc1, fc2, lab1, lab2){
  ok <- intersect(names(fc1), names(fc2))
  tmp <- cbind(fc1[ok], fc2[ok])
  tmp <- na.omit(tmp)
  correl <- cor.test(tmp[,1], tmp[,2], use = "pairwise.complete.obs", method = "spearman")
  smoothScatter(tmp[,1:2], main = "Foldchange old/young [log2]", ylab = lab2, xlab = lab1)
  abline(v = 0, h = 0, lty = 2)
  legend("topright", paste("Spearman", "\nRho:", signif(correl$estimate, 2), "\nP:", signif(correl$p.value, 2)), bty = "n")
  legend("topleft", paste(nrow(tmp), "genes"), bty = "n")
  dem.reg <- mcreg(tmp[,1], tmp[,2], method.reg = "Deming")
  abline(dem.reg@para[1:2], col = "blue")  
}

#pdf(useDingbats = F, "/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/figs/FoldchangeScatterplots.pdf", width = 11, height = 4)

par(mfrow = c(1, 3))
genScatterPlot(insilico_fc, bulk_fc, lab1 = 'scRNA-seq', lab2 = 'RNA-seq')
genScatterPlot(insilico_fc, prot_fc, lab1 = 'scRNA-seq', lab2 = 'Mass spec')
genScatterPlot(bulk_fc, prot_fc, lab1 = 'RNA-seq', lab2 = 'Mass spec')

#dev.off()