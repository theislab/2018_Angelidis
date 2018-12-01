# R code
# Lukas Simon
# This code will generate Fig S1

# Load data ####
library(Seurat)
library(Matrix)
library(cluster)
load("../data/SeuratObject.RData")

# Calculate Silhouette coefficient for mouse clustering ####
dist_euc <- dist(seu.ica@dr$ica@cell.embeddings, method = "euclidean")
identifier <- as.numeric(as.factor(seu.ica@meta.data$identifier))
sil <- silhouette(identifier, dist_euc)
summary(sil)[["avg.width"]]

# Generate Fig S1f ####
DimPlot(seu.ica, reduction.use = "tsne", group.by = "grouping", cols.use = c("red", "blue"))

# Generate Fig S1e ####
DimPlot(seu.ica, reduction.use = "tsne", group.by = "identifier")
p <- DimPlot(seu.ica, reduction.use = "tsne", group.by = "identifier", do.return = T)
p
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]]
farben <- unlist(lapply(split(pdata$colour, seu.ica@meta.data$identifier), unique))

# Generate Fig S1 a and b ####
VlnPlot(seu.ica, c("nGene", "nUMI"), group.by = "identifier", x.lab.rot = T, point.size.use = 0)

# Generate Fig S1d ####
par(mar = c(15, 5, 2, 2), xpd=TRUE)
matr <- unclass(table(seu.ica@meta.data$celltype, seu.ica@meta.data$identifier))
tmp <- t(matr / rowSums(matr))
bad <- c("Gamma-Delta_T_cells", "Plasma_cells", "low_quality_cells")
tmp <- tmp[, setdiff(colnames(tmp), bad)]
barplot(tmp, col = farben, ylab = "Fraction of cells", las = 2)
legend(34, 1, names(farben), col = farben, pch = 15)

# Generate Fig S1c ####
logdata <- read.csv("../data/MappingLogdata.csv", row.names = 1)
par(mfrow = c(1, 3), mar = c(5,5,5,5))
barplot(log10(logdata$totalreads), horiz = T, col = farben, names = rownames(logdata), las = 2, xlab = "Total reads (log10)")
barplot(logdata$percentmapped, horiz = T, col = farben, names = rownames(logdata), las = 2, xlab = "Percent mapped to mm10")
barplot(logdata$readlength, horiz = T, col = farben, names = rownames(logdata), las = 2, xlab = "Average read length")
