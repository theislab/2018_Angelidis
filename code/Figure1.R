# R code
# Lukas Simon
# This code will generate Fig1 panels

# Load packages and data ####
library(Seurat)
library(Matrix)
library(dendextend)
library(dplyr)
library(tidyr)
load("../data/SeuratObject.RData")
cells <- rownames(seu.ica@meta.data)[which(!is.na(seu.ica@meta.data$celltype))]
seu.ica <- SubsetData(seu.ica, cells.use = cells)

# Load marker table ####
markers <- read.delim("../data/Table S1_AllMarkersCelltypes.txt", header = T)
markers$cluster <- as.character(markers$cluster)

# Define meta celltypes ####
epithelium <- c("Type_2_pneumocytes", "Type1_pneumocytes", "Ciliated_cells", "Club_cells", "Goblet_cells")
apcs <- c("Ccl17+/Cd103-/Cd11b-_dendritic_cells","Cd103+/Cd11b-_dendritic_cells", "CD209+/Cd11b+_dendritic_cells")
lympho <- c("Natural_Killer_cells","B_cells", "Cd4+_T_cells", "CD8+_T_cells", "Plasma_cells")
mesenchyme <- c("Mesothelial_cells", "Smooth_muscle_cells", "Interstitial_Fibroblast", "Lipofibroblast")
endothelium <- c("lymphatic_endothelial_cells", "Capillary_endothelial_cells", "vascular_endothelial_cells", "Vcam1+_endothelial_cells")
otherLeuko <- c("Neutrophils", "Eosinophils", "Alveolar_macrophage", "classical_monocyte_(Ly6c2+)", "non-classical_monocyte_(Ly6c2-)", "Interstitial_macrophages", "Fn1+_macrophage")
hemos <- c('Megakaryocytes')
weird <- c("Mki67+_proliferating_cells", "low_quality_cells", "Gamma-Delta_T_cells", 'red_blood_cells')
metacells <- rep(NA, nrow(seu.ica@meta.data))
metacells[which(seu.ica@meta.data$celltype %in% epithelium)] <- "epithelium"
metacells[which(seu.ica@meta.data$celltype %in% apcs)] <- "apcs"
metacells[which(seu.ica@meta.data$celltype %in% lympho)] <- "lympho"
metacells[which(seu.ica@meta.data$celltype %in% mesenchyme)] <- "mesenchyme"
metacells[which(seu.ica@meta.data$celltype %in% endothelium)] <- "endothelium"
metacells[which(seu.ica@meta.data$celltype %in% otherLeuko)] <- "otherLeuko"
metacells[which(seu.ica@meta.data$celltype %in% hemos)] <- "hemos"
metacells[which(seu.ica@meta.data$celltype %in% weird)] <- "weird"
seu.ica@meta.data$metacelltype <- metacells

# Define coloring ####
colorRampAlpha <- function(..., n, alpha) {
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
}

myColors1 <- colorRampAlpha(c("red", "darkred"), alpha = 0.8, n = length(epithelium))
names(myColors1) <- epithelium
myColors2 <- colorRampAlpha(c("cyan", "darkcyan"), alpha = 0.8, n = length(lympho))
names(myColors2) <- lympho
myColors3 <- colorRampAlpha(c("lightblue", "darkblue"), alpha = 0.8, n = length(otherLeuko))
names(myColors3) <- otherLeuko
myColors4 <- colorRampAlpha(c("orange", "brown"), alpha = 0.8, n = length(mesenchyme))
names(myColors4) <- mesenchyme
myColors5 <- colorRampAlpha(c("green", "darkgreen"), alpha = 0.8, n = length(apcs))
names(myColors5) <- apcs
myColors6 <- colorRampAlpha(c("orchid", "darkorchid"), alpha = 0.8, n = length(endothelium))
names(myColors6) <- endothelium
myColors7 <- colorRampAlpha(c("pink", "salmon"), alpha = 0.8, n = length(hemos))
names(myColors7) <- hemos
myColors8 <- rep("grey", length(weird))
names(myColors8) <- weird
myColors9 <- rgb(0, 0, 0, 0.8)
names(myColors9) <- "Mki67+_proliferating_cells"

myColors <- c(myColors1, myColors2, myColors3, myColors4, myColors5, myColors6, myColors7, myColors8, myColors9)

farben <- myColors[seu.ica@meta.data$celltype]

asplit <- split(1:nrow(seu.ica@meta.data), seu.ica@meta.data$celltype)
textPos <- do.call(rbind, lapply(asplit, function(x) apply(seu.ica@dr$tsne@cell.embeddings[x, 1:2], 2, median)))

# Generate hierarchial tree ####
genClusterTree <- function(des){
  asplit <- split(1:nrow(des), as.character(des[,"cluster"]))
  allgenes <- unique(as.character(des$gene))
  fcs <- do.call(cbind, lapply(asplit, function(x){
    tmp <- des[x,]
    as.numeric(as.character(tmp[match(allgenes, as.character(tmp$gene)),"avg_logFC"]))
  }))
  rownames(fcs) <- allgenes
  fcs[which(is.na(fcs))] <- 0
  correl <- cor(fcs, use = "pairwise.complete")
  d <- as.dist(1 - correl)
  h <- hclust(d)
  h
}
tree <- genClusterTree(markers)
celltypes <- tree$labels[tree$order]

# Fig 1c - part 1 ####
tree <- as.dendrogram(tree)
labels_colors(tree) <- myColors[celltypes]
neworder <- c("Type_2_pneumocytes", "Ciliated_cells", "Club_cells", "Goblet_cells", "Type1_pneumocytes",
              "Mesothelial_cells", "Smooth_muscle_cells", "Interstitial_Fibroblast", "Lipofibroblast",
              "lymphatic_endothelial_cells", "Vcam1+_endothelial_cells", "Capillary_endothelial_cells", "vascular_endothelial_cells",
              "Megakaryocytes", "Fn1+_macrophage", "Eosinophils", "Neutrophils", "Alveolar_macrophage",
              "Mki67+_proliferating_cells", "Natural_Killer_cells", "Plasma_cells", "B_cells", "Cd4+_T_cells", "CD8+_T_cells",
              "Ccl17+/Cd103-/Cd11b-_dendritic_cells", "Cd103+/Cd11b-_dendritic_cells", "CD209+/Cd11b+_dendritic_cells",
              "classical_monocyte_(Ly6c2+)","non-classical_monocyte_(Ly6c2-)", "Interstitial_macrophages")
tree  <- reorder(tree, wts = order(match(neworder, celltypes)), agglo.FUN = "mean")
celltypes <- labels(tree)
labels(tree) <- paste(1:length(labels(tree)), labels(tree))
plot(tree)


# Generate tSNE colored by celltype - Fig 1b ####
numCelltype <- rep("NA", nrow(seu.ica@meta.data))
tmp <- cbind(celltypes, 1:length(celltypes))
tmp <- rbind(tmp, cbind(weird[-1], as.character(300 + 2:length(weird))))
#tmp <- rbind(tmp, c("Mki67+_proliferating_cells", 400))
lapply(1:nrow(tmp), function(x) numCelltype[which(seu.ica@meta.data$celltype == tmp[x, 1])] <<- tmp[x, 2])
seu.ica@meta.data$numCelltype <- numCelltype
myColors_num <- myColors
names(myColors_num) <- tmp[match(names(myColors), tmp[,1]),2]
DimPlot(seu.ica, reduction.use = "tsne", cols.use = myColors_num, group.by = "numCelltype", do.label = T, no.legend = T, vector.friendly = T)

# Generate Dotplot for marker genes ####
# Define marker genes
genes <- c("Vwf", "Tmem100", "Rtkn2", "Sftpd", "Acta2", "Mzb1", "Ngp", "Gzma",
           "Msln", "Ppbp", "Mmrn1", "Inmt", "C1qb", "Col1a2", "Dcn", "Bpifb1",
           "Prg4", "Cxcr2", "Scgb1a1", "Plac8", "Ly6c2", "Foxj1", "Cd8b1",
           "Trbc2", "Itgax", "Itgam", "Cd209a", "Itgae", "Ccl17", "Ednrb",
           "Cd79b", "Ear2")

# Remove bad cells
cells <- rownames(seu.ica@meta.data)[-which(seu.ica@meta.data$celltype %in% c("low_quality_cells", "red_blood_cells", "Gamma-Delta_T_cells"))]
tmp <- SubsetData(seu.ica, cells.use = cells)

PercentAbove <- function (x, threshold) length(x = x[x > threshold])/length(x = x)

data.to.plot <- data.frame(FetchData(object = tmp, vars.all = genes))
data.to.plot$cell <- rownames(x = data.to.plot)
data.to.plot$id <- factor(tmp@meta.data$celltype, levels = rev(celltypes))
data.to.plot <- data.to.plot %>% gather(key = genes.plot, value = expression, -c(cell, id))
data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression, threshold = 0))
data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale, max = 2.5, min = -2.5))
data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, levels = rev(x = sub(pattern = "-", replacement = ".", x = genes)))
data.to.plot$pct.exp[data.to.plot$pct.exp < 0] <- NA

# Fig 1c - part 2 ####
p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
  geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) + 
  scale_radius(range = c(0, 6)) + scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p