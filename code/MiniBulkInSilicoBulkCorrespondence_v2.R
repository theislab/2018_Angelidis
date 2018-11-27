# Load R libraries ####
library(Seurat)
library(DESeq2)
library(biomaRt)
library(wesanderson)
library(pheatmap)
library(gridExtra)
library(ggrepel)
library(limma)
library(edgeR)
library(mcr)

# Load Aging Seurat object ####
load("/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/SeuratObject.RData")

# Load count table ####
fc <- read.table("/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/fcCounts", header = T)
counts <- fc[, -c(1:6)]
colnames(counts) <- gsub("STARalignmentOutputs.", "", fixed = T, colnames(counts))
colnames(counts) <- gsub("Aligned.sortedByCoord.out.bam", "", fixed = T, colnames(counts))
rownames(counts) <- fc[,1]

# Define treatment vectors ####
treat <- read.csv(file = "/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/MiniBulk_sampleinfo.csv")
tmp <- do.call(rbind, lapply(as.character(treat[,2]), function(x) strsplit(x, "_", fixed = T)[[1]]))
tmp[,1] <- unlist(lapply(tmp[,1], function(x) substr(x, 0, 1)))
colnames(tmp) <- c("age", "celltype")
rownames(tmp) <- treat[,1]
treat <- data.frame(tmp)

# Merge treatment and counts ####
ok <- intersect(rownames(treat), colnames(counts))
counts <- counts[, ok]
treat <- treat[ok,]

# Remove samples without enough reads ####
bad <- which(colSums(counts) < 1e6)
counts <- counts[,-bad]
treat <- treat[-bad,]

# Get gene symbols ####
ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
t2g <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
t2g <- t2g[which(t2g$ensembl_gene_id %in% rownames(counts)),]

# Calculate in silico bulks per mouse ####
clusters <- list(c("Type_2_pneumocytes"), c("Alveolar_macrophage"))
pseudobulks <- lapply(clusters, function(x){
  cells <- rownames(seu.ica@meta.data)[which(seu.ica@meta.data$celltype %in% x)]
  subset <- SubsetData(seu.ica, cells.use = cells)
  asplit <- split(1:nrow(subset@meta.data), subset@meta.data$identifier)
  #asplit <- lapply(asplit, function(x) sample(x, 30, replace = T))
  means <- do.call(cbind, lapply(asplit, function(x) Matrix::rowSums(subset@data[, x])))
  means
})
tmp <- cbind(pseudobulks[[1]], pseudobulks[[2]])
colnames(tmp) <- paste(colnames(tmp), c(rep("EP", 15), rep("MAC", 15)))
pseudobulks <- tmp

# Load scRNA-seq marker genes ####
markers <- read.table("/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/Table S1_AllMarkersCelltypes.txt", sep = "\t", header = T)
markers <- markers[which(markers$cluster %in% c("Type_2_pneumocytes","Alveolar_macrophage")),]
markers <- unique(as.character(markers$gene[which(markers$p_val_adj < 0.1 & markers$avg_logFC > 0)]))

# Match pseuodbulk and minibulk on scRNA-seq marker genes ####
ok <- intersect(markers, t2g$external_gene_name)
pseudobulks_ok <- pseudobulks[ok,]
ok <- t2g$ensembl_gene_id[match(ok, t2g$external_gene_name)]
counts_ok <- counts[ok,]
rownames(counts_ok) <- rownames(pseudobulks_ok)

# Run PCA on pseudobulk and project mini bulk into it ####
norm2 <- function(matr){
  matr <- matr[which(rowSums(matr) > 0),]
  matr <- DGEList(matr)
  matr <- calcNormFactors(matr, method = "TMM")
  tmp <- voom(matr)$E
  t(apply(tmp, 1, function(x) (x - mean(x))/sd(x)))
}
pca <- prcomp(t(norm(pseudobulks_ok)))
preds <- predict(pca, t(norm(counts_ok)))

# Plot PCA ####
farben <- c("#F1BB7B", "#FD6467")
names(farben) <- c("EP", "MAC")

treat_pb <- c(rep("EP", 15), rep("MAC", 15))
treat_pb <- data.frame(celltype = treat_pb, age = as.character(seu.ica@meta.data$grouping[match(unlist(lapply(colnames(pseudobulks), function(x) strsplit(x, " ", fixed = T)[[1]][1])), seu.ica@meta.data$identifier)]))
treat_mb <- as.character(treat$celltype)
treat_mb[which(treat_mb != "EP")] <- "MAC"

tmp <- data.frame(pca$x[, 1:2], celltype = treat_pb$celltype)
p <- ggplot(tmp, aes(PC1, PC2, color = celltype)) + geom_point(shape=3, size = 3) + scale_color_manual(values=farben)
tmp <- data.frame(preds[, 1:2], celltype = treat_mb)
pca_plot <- p + geom_point(data = tmp, aes(x=PC1, PC2, color = celltype), shape=2, size = 3) + ggtitle("PCA")

tmp <- data.frame(pca$rotation[,1:2], symbol = rownames(pca$rotation))
p <- ggplot(tmp, aes(PC1, PC2, color = PC1)) + geom_point() + scale_colour_gradient(low = farben[1], high = farben[2])
ok <- c("Sftpc", "Ear2", "Ccl6", "Mrc1", "Cd68", "Itgb2", "Itgax", "Scd1", "Sftpd")
loadings_plot <- p + geom_point(data = tmp[ok, ], aes(x=PC1, PC2), color="red") + geom_text_repel(data = tmp[ok, ], aes(label = symbol), color = "black") + ggtitle("Loadings")

p_final <- grid.arrange(pca_plot, loadings_plot, ncol = 2)
ggsave("/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/figs/MinibulkPCA.pdf", plot = p_final, height = 5, width = 11)

# Plot correlation heatmap ####
correl <- cor(pseudobulks_ok, counts_ok, method = "spearman")
anno_row <- data.frame(celltype = treat_pb$celltype)
rownames(anno_row) <- rownames(correl)
anno_col <- data.frame(celltype = treat_mb)
rownames(anno_col) <- colnames(correl)
ann_colors <- list(celltype = farben)
pheatmap(correl, annotation_col = anno_col, annotation_row = anno_row, annotation_colors = ann_colors, show_rownames = F, show_colnames = F)

# Match bulk and sc all genes ####
ok <- intersect(rownames(pseudobulks), t2g$external_gene_name)
pseudobulks_ok <- pseudobulks[ok,]
ok <- t2g$ensembl_gene_id[match(ok, t2g$external_gene_name)]
counts_ok <- counts[ok,]
rownames(counts_ok) <- rownames(pseudobulks_ok)

# Calculate differential expression ####
runRegression <- function(matr, treat){
  matr <- matr[-which(rowSums(matr) == 0),]
  matr <- DGEList(round(matr))
  matr <- calcNormFactors(matr, method = "TMM")
  v <- voom(matr)
  design <- -model.matrix(~ treat)
  fit <- lmFit(v, design = design)
  fit <- eBayes(fit)
  topTable(fit, adjust.method = 'BH', number = Inf, sort = "none")
}

ok <- which(treat_mb == "EP")
res_mb_ep <- runRegression(counts_ok[, ok], treat$age[ok])
ok <- which(treat_mb == "MAC")
res_mb_mac <- runRegression(counts_ok[, ok], treat$age[ok])

write.csv(res_mb_ep, file = "/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/Minibulk_EP.csv")
write.csv(res_mb_mac, file = "/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/Minibulk_MAC.csv")

# Compare to cell type resolved fold changes ####
sc_de <- read.delim("/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/AllDEtable.txt")

pdf(useDingbats = F, "/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/figs/MinibulkScatterplot.pdf", width = 12.6, height = 7)
par(mfrow = c(1, 2))
ok <- grep("Alveolar_macrophage", colnames(sc_de))
tmp <- as.data.frame(sc_de[,ok])
tmp <- tmp[which(!is.na(tmp$Alveolar_macrophage.avg_logFC)),]
tmp <- tmp[which(tmp$Alveolar_macrophage.p_val < 0.1),]
ok <- intersect(rownames(tmp), rownames(res_mb_mac))
aframe_mac <- data.frame(tmp[ok,], res_mb_mac[ok,])
smoothScatter(aframe_mac[, "logFC"], aframe_mac[, "Alveolar_macrophage.avg_logFC"], main = "Alveolar macrophages", xlab = "Mini bulk", ylab = "scRNA-seq", col = rgb(0, 0, 0, 0.5), pch = 19)
dem.reg <- mcreg(aframe_mac[, "logFC"], aframe_mac[, "Alveolar_macrophage.avg_logFC"], method.reg = "Deming")
#abline(dem.reg@para[1:2], col = "blue")
abline(v = 0, h = 0)
#fish <- cor.test(aframe_mac[, "Alveolar_macrophage.avg_logFC"], aframe_mac[, 'logFC'], method = "spearman")
fish <- fisher.test(table(aframe_mac[, "Alveolar_macrophage.avg_logFC"] > 0, aframe_mac[, 'logFC'] > 0))
legend("bottomright", paste("Fisher P:", signif(fish$p.value, 2),"OR:", signif(fish$estimate, 2)), bty = "n")
genes <- c("Marco", "Awat1", "Cebpb", "Fabp1", "Atp1b1", "mt-Cytb")
text(aframe_mac[genes, "logFC"], aframe_mac[genes, "Alveolar_macrophage.avg_logFC"], genes)
or_mac <- c(fish$estimate, fish$conf.int[1:2])

ok <- grep("Type_2_pneumocytes", colnames(sc_de))
tmp <- as.data.frame(sc_de[,ok])
tmp <- tmp[which(!is.na(tmp$Type_2_pneumocytes.avg_logFC)),]
tmp <- tmp[which(tmp$Type_2_pneumocytes.p_val < 0.1),]
ok <- intersect(rownames(tmp), rownames(res_mb_ep))
aframe_ep <- data.frame(tmp[ok,], res_mb_ep[ok,])
smoothScatter(aframe_ep[, "logFC"], aframe_ep[, "Type_2_pneumocytes.avg_logFC"], main = "Type 2 pneumocytes", xlab = "Mini bulk", ylab = "scRNA-seq", col = rgb(0, 0, 0, 0.5), pch = 19)
dem.reg <- mcreg(aframe_ep[, "logFC"], aframe_ep[, "Type_2_pneumocytes.avg_logFC"], method.reg = "Deming")
#abline(dem.reg@para[1:2], col = "blue")
abline(v = 0, h = 0)
#fish <- cor.test(aframe_ep[, "logFC"], aframe_ep[, "Type_2_pneumocytes.avg_logFC"], method = "spearman")
fish <- fisher.test(table(aframe_ep[, "logFC"] > 0, aframe_ep[, "Type_2_pneumocytes.avg_logFC"] > 0))
legend("bottomright", paste("Fisher P:", signif(fish$p.value, 2),"OR:", signif(fish$estimate, 2)), bty = "n")
genes <- c("Lyz1", "B2m", "Scd1", "H2-Q7", "mt-Nd3", "Sparc", "Scd1", "H2-K1")
text(aframe_ep[genes, "logFC"], aframe_ep[genes, "Type_2_pneumocytes.avg_logFC"], genes)
or_ep <- c(fish$estimate, fish$conf.int[1:2])
dev.off()

tmp <- data.frame(rbind(or_mac, or_ep))
tmp$cell <- c("mac", "ep")
colnames(tmp) <- c("OR", "low", "high", "cell")
ggplot(tmp, aes(x = OR, y = cell)) + geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = high, xmin = low), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") +
  theme_bw()+
  theme(panel.grid.minor = element_blank()) +
  xlab("Odds ratio")
ggsave("/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/figs/OddsRatioMiniBulkDE.pdf", height = 5, width = 6)


write.csv(aframe_mac, file = "/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/Minibulk_MAC_wSC.csv",quote = F)
write.csv(aframe_ep, file = "/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/Minibulk_EP_wSC.csv",quote = F)

# Highlight some example genes ####
ok <- which(treat_mb == "EP")
treat_tmp <- treat$age[ok]
matr <- counts_ok[, ok]
matr <- matr[-which(rowSums(matr) == 0),]
matr <- DGEList(round(matr))
matr <- calcNormFactors(matr, method = "TMM")
v_ep <- voom(matr)
pdf(useDingbats = F, "/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/figs/Scd1H2K1_minibulk_boxplot.pdf", height = 5, width = 7)
par(mfrow = c(1, 2))
boxplot(split(v_ep$E["Scd1",], treat_tmp)[c("Y", "A")], ylab = "Normalized Expression", main = "Scd1", col = c("blue", "red"))
boxplot(split(v_ep$E["H2-K1",], treat_tmp)[c("Y", "A")], ylab = "Normalized Expression", main = "H2-K1", col = c("blue", "red"))
dev.off()

seu.ica@meta.data$grouping <- factor(seu.ica@meta.data$grouping, levels = c("3m", "24m"))
VlnPlot(SubsetData(seu.ica, subset.name = "celltype", accept.value = "Type_2_pneumocytes"), features.plot = c("Scd1"),group.by = "grouping", cols.use = c("blue", "red"))
ggsave("/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/figs/Scd1_vlnPlot.pdf")
VlnPlot(SubsetData(seu.ica, subset.name = "celltype", accept.value = "Type_2_pneumocytes"), features.plot = c("H2-K1"),group.by = "grouping", cols.use = c("blue", "red"))
ggsave("/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/figs/H2K1_vlnPlot.pdf")

ok <- which(treat_mb == "MAC")
treat_tmp <- treat$age[ok]
matr <- counts_ok[, ok]
zeros <- which(rowSums(matr) == 0)
if(length(zeros) > 0) matr <- matr[-zeros,]
matr <- DGEList(round(matr))
matr <- calcNormFactors(matr, method = "TMM")
v_mac <- voom(matr)
#boxplot(split(v_mac$E["Ear2",], treat_tmp), ylab = "Expression", main = "Ear2", col = c("red", "blue"))
#boxplot(split(v_mac$E["Marco",], treat_tmp)[c("Y", "A")], ylab = "Expression", main = "Marco", col = c("blue", "red"))
#boxplot(split(v_mac$E["Cd63",], treat_tmp), ylab = "Expression", main = "Cd63", col = c("red", "blue"))
pdf(useDingbats = F, "/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/figs/Atp1b1_mtCytb_minibulk_boxplot.pdf", height = 5, width = 7)
par(mfrow = c(1, 2))
boxplot(split(v_mac$E["Atp1b1",], treat_tmp)[c("Y", "A")], ylab = "Normalized Expression", main = "Atp1b1", col = c("blue", "red"))
boxplot(split(v_mac$E["mt-Cytb",], treat_tmp)[c("Y", "A")], ylab = "Normalized Expression", main = "mt-Cytb", col = c("blue", "red"))
dev.off()

VlnPlot(SubsetData(seu.ica, subset.name = "celltype", accept.value = "Alveolar_macrophage"), features.plot = c("Atp1b1"),group.by = "grouping", cols.use = c("blue", "red"))
ggsave("/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/figs/Atp1b1_vlnPlot.pdf")
VlnPlot(SubsetData(seu.ica, subset.name = "celltype", accept.value = "Alveolar_macrophage"), features.plot = c("mt-Cytb"),group.by = "grouping", cols.use = c("blue", "red"))
ggsave("/Users/lukas.simon/OneDrive/Miko/Helmholtz/Schiller/outputs/Aging Project/figs/mt-Cytb_vlnPlot.pdf")

runRegression2 <- function(matr, treat){
  des <- DESeqDataSetFromMatrix(countData = round(matr), colData = data.frame(treat), design = ~ treat)
  des <- DESeq(des)
  res <- results(des, contrast = c("treat", "A", "Y"))
  res
}
ok <- which(treat_mb == "EP")
res_mb_ep <- runRegression2(counts_ok[, ok], treat$age[ok])
ok <- which(treat_mb == "MAC")
res_mb_mac <- runRegression2(counts_ok[, ok], treat$age[ok])


