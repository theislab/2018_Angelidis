# R 
# Lukas Simon
# Wed Nov 21st 2018
# This code will generate figures for the cell cycle analysis

set.seed(15)

# Load R libraries ####
library(Seurat)
library(Matrix)
library(wesanderson)

# Load Seurat object ####
load("../data/SeuratObject.RData")

# Load cell cycle genes ####
cc.genes <- scan(what = character(), "../data/regev_lab_cell_cycle_genes.txt", sep = "\n")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

# Convert to mouse symbols ####
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
s.genes <- convertHumanGeneList(s.genes)
g2m.genes <- convertHumanGeneList(g2m.genes)

# Generate Featureplot of entire data set ####
g2m.genes <- intersect(g2m.genes, rownames(seu.ica@data))
s.genes <- intersect(s.genes, rownames(seu.ica@data))
g2m_sig <- Matrix::colMeans(seu.ica@data[g2m.genes,])
s_sig <- Matrix::colMeans(seu.ica@data[s.genes,])
seu.ica@meta.data$g2m_sig <- g2m_sig
seu.ica@meta.data$s_sig <- s_sig
FeaturePlot(seu.ica, features.plot = c("g2m_sig", "s_sig"), cols.use = c("lightgrey", "blue"), pt.size = 1.5, vector.friendly = T)

# Restrict to cell cycle cluster ####
subset <- SubsetData(seu.ica, ident.use = 18)
subset <- FindVariableGenes(subset, do.plot = F)
subset <- ScaleData(subset, vars.to.regress = "nUMI")
subset <- RunPCA(subset, pc.genes = c(s.genes, g2m.genes), do.print = F)
subset <- CellCycleScoring(subset, g2m.genes = g2m.genes, s.genes = s.genes, set.ident = TRUE)
p <- PCAPlot(subset, do.return = T)
ggplot(p$data, aes(PC1, PC2, color = ident))+geom_point()+stat_ellipse()

# Regress out cell cycle and cluster ####
subset <- ScaleData(object = subset, vars.to.regress = c("S.Score", "G2M.Score", "nUMI"), display.progress = FALSE)
subset <- RunPCA(subset, pc.genes = c(s.genes, g2m.genes), do.print = F)
p <- PCAPlot(subset, do.return = T)
ggplot(p$data, aes(PC1, PC2, color = ident))+geom_point()+stat_ellipse()

# Load celltype markers ####
marker_table <- read.delim('../data/Table S1_AllMarkersCelltypes.txt')
markers <- unique(marker_table$gene[which(marker_table$p_val_adj < 0.05 & marker_table$avg_logFC > 0)])

# Rerun clustering ####
subset <- RunPCA(subset, pc.genes = markers, do.print = F)
subset <- RunTSNE(subset, dims.use = 1:4, reduction.use = 'pca')
subset <- FindClusters(subset, print.output = F, force.recalc = T, dims.use = 1:4)
TSNEPlot(subset, pt.size = 2, no.axes = T, no.legend = T, colors.use = wes_palette("GrandBudapest1", 3), title = "Louvain clustering")

# Generate Featureplot ####
FeaturePlot(subset, features.plot = c("Ear2", "Trbc2", "Sftpd"), cols.use = c("lightgrey", "blue"), pt.size = 1.5, no.axes = T)

# Calculate age association for the rate of falling into the proliferating cluster #### 
counts <- cbind(table(subset@meta.data$identifier),table(seu.ica@meta.data$identifier))
age <- as.character(seu.ica@meta.data$grouping[match(rownames(counts), seu.ica@meta.data$identifier)])
fit <- summary(glm(counts ~ age, family = binomial(link="logit")))
boxplot(split(counts[,1]/rowSums(counts), age), ylab = "Fraction of proliferating cells", col = c("red", "blue"),
        main = paste("General linearized binomial P", signif(coefficients(fit)["age3m", 4], 2)))
