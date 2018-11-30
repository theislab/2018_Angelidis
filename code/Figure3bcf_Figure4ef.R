# R code
# Lukas Simon
# This code will generate panels for Figure 3

# Load R libraries ####
library(biomaRt)
library(limma)
library(Seurat)
library(pheatmap)
library(preprocessCore)
library(corrplot)

# Load whole lung bulk expression ####
fc <- read.table("../data/raw_counts_whole_lung_bulk.txt", header = T)
counts <- fc[, -c(1:6)]
rownames(counts) <- fc[,1]
counts <- counts[which(rowSums(counts) > 0),]
bulk <- counts
age_bulk <- c(rep("old", 3), rep("young", 3))

# Perform DE analysis on whole lung bulk data ###
des <- DESeqDataSetFromMatrix(countData = counts, colData = data.frame(age_bulk), design = ~ age_bulk)
des <- DESeq(des)
res <- results(des, contrast = c('age_bulk', 'old', 'young')) 

# Get gene symbols from biomart ####
ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
t2g <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
t2g <- t2g[which(t2g$ensembl_gene_id %in% rownames(bulk)),]

# Save DE result table ####
res$symbol <- t2g$external_gene_name[match(rownames(res), t2g$ensembl_gene_id)]
final <- data.frame(res, counts(des, normalized = T))

# Perform bulk deconvolution on whole lung results ####
markers <- read.delim('../data/Table S1_AllMarkersCelltypes.txt')
tmp <- markers[which(markers$avg_logFC > 1 & markers$p_val_adj < 0.1),]
schiller_genes <- split(as.character(tmp$gene), tmp$cluster)

tmp <- do.call(rbind, lapply(schiller_genes, function(x){
  ok <- intersect(x, res$symbol)
  set <- res$log2FoldChange[match(ok, res$symbol)]
  rest <- res$log2FoldChange[-match(ok, res$symbol)]
  pval <- ks.test(set, rest)$p.value
  coef <- mean(set) - mean(rest)
  c(fc_diff = coef, pval = pval)
}))
tmp[which(tmp[,2] < 1e-50), 2] <- 1e-50
tmp <- data.frame(tmp)
tmp$celltype <- rownames(tmp)

# Generate Fig 4e ####
ggplot(tmp, aes(fc_diff, -log10(pval))) + geom_point() +
  geom_text_repel(data = tmp, aes(label = celltype), color = "black", size = 3) +
  geom_vline(xintercept = 0)

# Generate Fig 4f ####
genECDFplot2 <- function(celltype){
  ok <- intersect(schiller_genes[[celltype]], res$symbol)
  type <- rep('Rest', nrow(res))
  type[match(ok, res$symbol)] <- celltype 
  df <- data.frame(foldchange = res$log2FoldChange, type)
  ggplot(df, aes(foldchange, colour = type)) + scale_color_manual(values=c("red", "black")) + stat_ecdf() + scale_x_continuous(limits = c(-10, 10)) 
}
genECDFplot2('Ciliated_cells')

# Load and generate in silico bulk data ####
load("../data/SeuratObject.RData")
asplit <- split(1:nrow(seu.ica@meta.data), seu.ica@meta.data$identifier)
tmp <- do.call(cbind, lapply(asplit, function(x) Matrix::rowSums(seu.ica@raw.data[,x])))
tmp <- tmp[which(rowSums(tmp) > 0),]
insilico <- tmp
age_insilico <- as.character(seu.ica@meta.data$grouping[match(colnames(insilico), seu.ica@meta.data$identifier)])
age_insilico <- gsub("3m", "young", age_insilico)
age_insilico <- gsub("24m", "old", age_insilico)

# Show correspondence between bulk and insilico ####
ok <- intersect(rownames(insilico), t2g$external_gene_name)
tmp <- bulk[t2g$ensembl_gene_id[match(ok, t2g$external_gene_name)], ]

# Generate Fig 3b - part 1 ####
tmp2 <- cbind(tmp, insilico[ok, ])
correl <- cor(tmp2, method = 'spearman')
corrplot(correl, method = 'ellipse', type = 'upper', cl.lim = c(0,1))

# Generate Fig 3b - part 2 ####
insilico_means <- rowMeans(insilico[ok,])
bulk_means <- rowMeans(tmp)
plot(log(insilico_means), log(bulk_means), col = rgb(0, 0, 0, 0.5), pch = 19, xlab = 'Average in silico bulk [log]', ylab = 'Average whole lung bulk [log]')
cor.test(log(insilico_means), log(bulk_means), method = 'spearman')
dem.reg <- mcreg(log(insilico_means), log(bulk_means), method.reg = "Deming")
abline(dem.reg@para[1:2], col = "red", lwd = 2)
abline(0, 1)

# Load protein expression data ####
prot <- read.delim("../data/Table S2_TissueProteome_imputed.txt", header = T)
expr <- data.matrix(prot[,1:8])
genes <- as.character(prot$Gene.names)
protein <- do.call(rbind, lapply(1:length(genes), function(x){
  syms <- strsplit(genes[x], ";", fixed = T)[[1]]
  subm <- do.call(rbind, lapply(1:length(syms), function(y) expr[x,]))
  rownames(subm) <- syms
  subm
}))
age_prot <- c(rep("old", 4), rep("young", 4))

# Match genes for all three matrices ####
ok <- intersect(rownames(protein), rownames(insilico))
ok <- intersect(ok, t2g$external_gene_name)
protein_ok <- protein[ok,]
bulk_ok <- bulk[t2g$ensembl_gene_id[match(ok, t2g$external_gene_name)], ]
insilico_ok <- insilico[ok, ]
rownames(bulk_ok) <- rownames(protein_ok) <- rownames(insilico_ok) <- ok

# Normalize and merge all three matrices ####
bulk_ok <- voom(bulk_ok)$E
insilico_ok <- voom(insilico_ok)$E
merged <- cbind(bulk_ok, insilico_ok, protein_ok)
merged <- normalize.quantiles(merged)
rownames(merged) <- rownames(bulk_ok)

# Define sample attributes ####
batch <- c(rep("bulk", ncol(bulk_ok)), rep("insilico", ncol(insilico_ok)), rep("protein", ncol(protein_ok)))
age <- c(age_bulk, age_insilico, age_prot)

farben <- c("blue", "red")
names(farben) <- c('young', 'old')

shape <- c(1, 2, 3)
names(shape) <- c("bulk", "insilico", "protein")

# Calculate PCA ####
pca <- prcomp(t(merged))

# Generate Fig 3c ####
par(mfrow = c(1, 2))
plot(pca$x[,1:2], col = farben[age], pch = shape[batch])
legend("bottomright", c(names(farben), names(shape)), pch = c(16, 16, shape), bty = "n", col = c("blue", "red", "black", "black", "black"))

plot(pca$x[,2:3], col = farben[age], pch = shape[batch])
legend("bottomleft", c(names(farben), names(shape)), pch = c(16, 16, shape), bty = "n", col = c("blue", "red", "black", "black", "black"))
abline(h = 0)

# Generate Fig 3f ####
col4a_genes <- c('Col4a1', 'Col4a2', 'Col4a3', 'Col4a4', 'Col4a5', 'Col4a6')
norm1 <- function(matr, treat){
  t(apply(matr, 1, function(x) (x - mean(x))/sd(x)))
}

batch <- c(rep("insilico", ncol(insilico_ok)), rep("bulk", ncol(bulk_ok)),rep("protein", ncol(protein_ok)))
age <- c(age_insilico, age_bulk, age_prot)

final <- cbind(norm1(insilico_ok[col4a_genes,]),
               norm1(bulk_ok[col4a_genes,]),
               norm1(protein_ok[col4a_genes,]))

tmp <- do.call(c, apply(final, 1, function(x) split(x, paste(batch, age))))
ord <- c(seq(from = 1, to = 36, by = 6),
         seq(from = 2, to = 36, by = 6),
         seq(from = 3, to = 36, by = 6),
         seq(from = 4, to = 36, by = 6),
         seq(from = 5, to = 36, by = 6),
         seq(from = 6, to = 36, by = 6))
tmp <- tmp[ord]
order4plot <- c(paste(col4a_genes, 'insilico young', sep = '.'),
                paste(col4a_genes, 'insilico old', sep = '.'),
                paste(col4a_genes, 'bulk young', sep = '.'),
                paste(col4a_genes, 'bulk old', sep = '.'),
                paste(col4a_genes, 'protein young', sep = '.'),
                paste(col4a_genes, 'protein old', sep = '.'))
                
boxplot(tmp[order4plot], las = 2, ylab = 'Expression z-score', outline = F)
abline(h = 0)

