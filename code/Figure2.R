# R code
# Lukas Simon
# This code will generate Figure 2 panels

# Load R libraries ####
library(Seurat)
library(Matrix)
library(pheatmap)
library(sfsmisc)
library(MASS)

# Load Aging Seurat object ####
load("../data/SeuratObject.RData")

# Define the celltypes ####
celltypes <- unique(seu.ica@meta.data$celltype)
celltypes <- celltypes[which(!is.na(celltypes))]
celltypes <- setdiff(celltypes, c("red_blood_cells", "low_quality_cells", "Gamma-Delta_T_cells"))

# Define function to calculate euclidean distances accounting for cell number and nUMI ####
getEuclideanDistance <- function(celltype, lowcv = T){
  print(paste("Working on", celltype))
  library(hopach)
  tmp <- SubsetData(seu.ica, subset.name = "celltype", accept.value = celltype, subset.raw = T)
  expr <- tmp@raw.data
  
  zeros <- which(Matrix::rowSums(expr) == 0)
  expr <- data.matrix(expr[-zeros,])
  
  Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)
  
  nsample <- min(table(tmp@meta.data$grouping)[c("3m", "24m")])
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  } 
  old_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "24m")], nsample)
  young_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "3m")], nsample)
  ds_expr_r <- ds_expr[, c(young_r, old_r)]
  
  if(lowcv){
    getLowCVgenes <- function(matr){
      means <- Matrix::rowMeans(matr)
      bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)))
      mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
      asplit <- split(names(means), mean_bin)
      genes <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
        coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
        bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
      })))
      genes
    }
    genes <- getLowCVgenes(ds_expr_r)
  }
  else{
    genes <- rownames(ds_expr_r)
  }
    
  calcEuclDist <- function(matr, young, old){
    tmp <- data.matrix(sqrt(matr[genes, young]))
    mean <- rowMeans(sqrt(matr[genes, young]))
    d_young <- distancevector(t(tmp), mean , d="euclid")
    names(d_young) <- young
    
    tmp <- data.matrix(sqrt(matr[genes, old]))
    mean <- rowMeans(sqrt(matr[genes, old]))
    d_old <- distancevector(t(tmp), mean , d="euclid")
    names(d_old) <- old
    
    list(young = d_young, old = d_old)
  }
  ds <- calcEuclDist(matr = ds_expr_r, old = old_r, young = young_r)
  ds
}

# Run for all celltypes ####
res <- lapply(celltypes, function(x) getEuclideanDistance(x, lowcv = F))
names(res) <- celltypes
nulls <- which(unlist(lapply(res, class)) == "NULL")
celltypes[nulls]

# Remove celltypes without enough cells ####
res <- res[-nulls]
res_original <- res

# Calculate mean differences and p-values ####
diffs <- unlist(lapply(res_original, function(x) log2(mean(x[[2]]) / mean(x[[1]]))))
pvals <- unlist(lapply(res_original, function(x) wilcox.test(x[[1]], x[[2]])$p.value))
adj_pvals <- p.adjust(pvals, method = "BH")
sizes <- (-log10(adj_pvals))
sizes[which(sizes < 1)] <- 1
sizes[which(sizes > 4)] <- 4
sizes <- sizes * 0.75
farben <- rep("grey", length(adj_pvals))
farben[which(adj_pvals < 0.05)] <- "purple"

# Generate Fig 2a ####
ord <- rev(order(diffs))
par(mar = c(15,5,2,5))
boxplot(do.call(c, res[ord]), las = 2, outline = F, col = c("blue", "red"), ylab = "Transcriptional noise", xaxt = 'n')
axis(1, at = seq(from = 1.5, to = 56.5, by = 2), names(res)[ord], las = 2)

# Calculate mouse means ####
showMouseLevel <- function(celltype){
  tmp <- res_original[[celltype]]
  lapply(tmp, function(x){
    nom <- names(x)
    ids <- unlist(lapply(nom, function(y) strsplit(y, ':', fixed = T)[[1]][1]))
    unlist(lapply(split(x, ids), median))
  })
}
plotMouseLevel <- function(celltype){
  noise <- showMouseLevel(celltype)
  pval <- wilcox.test(noise[[1]], noise[[2]])$p.value
  boxplot(noise, main = paste(celltype, paste('Wilcox P',signif(pval, 2)), sep = '\n'), col = c('blue', 'red'), names = NA)
}

# Correlate mean differences ####
diffs <- unlist(lapply(res_original, function(x) log2(mean(x[[2]]) / mean(x[[1]]))))
diffs_mouselevel <- unlist(lapply(names(res_original), function(x){
  tmp <- showMouseLevel(x)
  log2(mean(tmp[[2]]) / mean(tmp[[1]]))
}))
names(diffs_mouselevel) <- names(diffs)

# Generate Fig 2b ####
highlight <- c("Type_2_pneumocytes", "Neutrophils", "Club_cells", "Plasma_cells")
plot(diffs_mouselevel, diffs, main = 'Transcriptional noise ratio [old/young log2]', ylab = 'Cell level', xlab = 'Mouse level', pch = 19, cex = sizes, col = rgb(0, 0, 0, 0.7))
abline(h = 0, v = 0, lty = 2)
text(diffs_mouselevel[highlight], diffs[highlight], highlight, cex = 0.8)
f <- f.robftest(rlm(diffs ~ diffs_mouselevel), var = "diffs_mouselevel")
abline(coefficients(rlm(diffs ~ diffs_mouselevel)), col = "red")
legend('bottomright', paste("Robust F test P", signif(f$p.value, 2)), bty = 'n')

# Define function to get Spearman correlations ####
getSpearmanCorrelations <- function(celltype){
  print(paste("Working on", celltype))
  
  tmp <- SubsetData(seu.ica, subset.name = "celltype", accept.value = celltype, subset.raw = T)
  expr <- tmp@raw.data
  zeros <- which(Matrix::rowSums(expr) == 0)
  expr <- data.matrix(expr[-zeros,])
  Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)
  nsample <- min(table(tmp@meta.data$grouping)[c("3m", "24m")])
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  }
  old_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "24m")], nsample)
  young_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "3m")], nsample)
  ds_expr_r <- ds_expr[, c(young_r, old_r)]
  cor_young <- cor(method = "spearman", ds_expr_r[, young_r])
  cor_young <- unlist(cor_young[upper.tri(cor_young)])
  cor_old <- cor(method = "spearman", ds_expr_r[, old_r])
  cor_old <- unlist(cor_old[upper.tri(cor_old)])
  list(cor_old, cor_young)
}

# Run for all celltypes ####
res <- lapply(celltypes, getSpearmanCorrelations)
names(res) <- celltypes
nulls <- which(unlist(lapply(res, class)) == "NULL")
celltypes[nulls]

# Remove celltypes without enough cells ####
res_cor <- res[-nulls]

# Generate Fig 2c ####
diffs_cor <- unlist(lapply(res_cor, function(x) log2(mean(1 - x[[1]]) / mean(1 - x[[2]]))))
plot(diffs_cor, diffs, ylab = "Euclidean distance [old/young log2]", xlab = "1 - Spearman correlation [old/young log2]", main = "Transcriptional noise", pch = 19, cex = sizes, col = rgb(0, 0, 0, 0.7))
abline(v = 0, h = 0, lty = 2)
text(diffs_cor[highlight], diffs[highlight], highlight)
f <- f.robftest(rlm(diffs ~ diffs_cor), var = "diffs_cor")
abline(coefficients(rlm(diffs ~ diffs_cor)), col = "red")
legend("topleft", paste("Robust F test P", signif(f$p.value, 2)), bty = "n")

# Generate Fig 2d ####
plot(density(1- res_cor[["Type_2_pneumocytes"]][[2]]), col = "blue", ylab = "Density", xlab = "1 - Spearman Correlation coefficient", main = "Type 2 pneumocytes", lwd = 2)
lines(density(1- res_cor[["Type_2_pneumocytes"]][[1]]), col = "red", lwd = 2)
correl <- ks.test(res_cor[["Type_2_pneumocytes"]][[2]], res_cor[["Type_2_pneumocytes"]][[1]])
legend("topright", paste("KS P", signif(correl$p.value, 2)), bty = "n")
legend("topleft", c("young", "old"), col = c("blue", "red"), bty = "n", pch = 19)

