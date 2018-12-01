# R code
# Lukas Simon
# This code will generate Figure S2

# Load R packages ####
library(readxl)
library(matchSCore)

# Define some functions ####
matchSCore2 <- function(gene_cl.ref,gene_cl.obs,tissue,ylab,xlab){
  
  score=0
  lab.ref=seq(1:length(gene_cl.ref))
  lab.obs=seq(1:length(gene_cl.obs))
  anno.lab=vector()
  max_ji=vector()
  ji_mat=vector()
  
  for(i in 1:length(lab.ref)){
    len1=length(gene_cl.ref[[i]])
    JI=vector()
    for(ind in 1:length(lab.obs)){
      len2=length(gene_cl.obs[[ind]])
      I=length(intersect(gene_cl.ref[[i]],gene_cl.obs[[ind]]))
      J=I/(len1+len2-I)
      JI=append(JI,J,length(JI))
    }
    score=sum(score,max(JI))
    anno.lab = append(anno.lab,which(JI==max(JI)),after = length(anno.lab)) 
    max_ji = append(max_ji,max(JI),after = length(max_ji))
    ji_mat =rbind(ji_mat,JI)
  }
  score=score/length(gene_cl.ref)
  colnames(ji_mat)=names(gene_cl.obs)
  rownames(ji_mat)=names(gene_cl.ref)
  file_name=paste("summary_",tissue,".pdf",sep="")
  gg<-summary_ggplot(data=ji_mat,file_name,ylab,xlab)
  return(list(matchScore=score,labels=anno.lab,max_JI=max_ji,JI.mat=ji_mat,ggplot=gg))
  
}
extendLabels <- function(tmp){
  tmp <- as.data.frame(tmp)
  col <- grep('cell', ignore.case = T, colnames(tmp))
  celltypes <- unique(tmp[, col])
  celltypes <- celltypes[which(!is.na(celltypes))]
  lapply(celltypes, function(x){
    ind <- match(celltypes, tmp[, col]) 
    hit <- match(x, celltypes)
    from <- ind[hit]
    to <- ind[hit+1] - 1
    if((hit + 1) > length(celltypes)) to <- nrow(tmp)
    tmp[from:to, col] <<- x
  })
  tmp  
}
summary_ggplot <- function(data,name_file,ylab,xlab){
  
  library(ggplot2)
  library(reshape2)      
  library(grid)          
  
  my_df <- data.frame(t(data),check.names = F,check.rows = F)
  my_df.melt <-  melt(cbind(x=1:nrow(my_df),my_df),id ="x")
  
  gg <- ggplot(my_df.melt, aes(x=factor(x),y=variable,fill=value)) + labs(x=xlab,y=ylab)+
    geom_tile(aes(fill = value)) + scale_x_discrete(lab=rownames(my_df))+ 
    theme(axis.text.x=element_text(angle=30,hjust = 1,size=16),axis.text.y = element_text(size=16),axis.title = element_text(size=16))+
    geom_text(aes(label = round(value, 2))) + 
    scale_fill_gradient(low = "white", high = "red",name="matchSCore\n") 
  
  
  return(gg)
}

# Load MCA data ####
path <- '../data/MCACelltypeMarkers.xlsx'
sheets <- excel_sheets(path = path)
mca <- as.data.frame(read_excel(path, sheet = "Lung"))
mca <- extendLabels(mca)
mca <- mca[which(mca$avg_logFC > 1 & mca$p_val_adj < 0.1),]
mca_genes_lung <- split(mca$gene, mca$`Cell type`)

mca <- as.data.frame(read_excel(path, sheet = "Peripheral blood"))
mca <- extendLabels(mca)
mca <- mca[which(mca$avg_logFC > 1 & mca$p_val < 3e-5),]
mca_genes_blood <- split(mca$gene, mca$`Cell Type`)

# Load Schiller data ####
markers <- read.delim('../data/Table S1_AllMarkersCelltypes.txt')
tmp <- markers[which(markers$avg_logFC > 1 & markers$p_val_adj < 0.1),]
schiller_genes <- split(as.character(tmp$gene), tmp$cluster)

# Generate Fig S2a ####
out_aging_vs_lung <- matchSCore2(gene_cl.ref = mca_genes_lung, gene_cl.obs = schiller_genes, tissue = "Lung", xlab = "Schiller", ylab = "MCA Lung")

# Generate Fig S2b ####
out_aging_vs_blood <- matchSCore2(gene_cl.ref = mca_genes_blood, gene_cl.obs = schiller_genes, tissue = "Blood", xlab = "Schiller", ylab = "MCA Blood")

# Generate Fig S2c ####
out_blood_vs_lung <- matchSCore2(gene_cl.ref = mca_genes_blood, gene_cl.obs = mca_genes_lung, tissue = "MCA Blood", xlab = "MCA Blood", ylab = "MCA Lung")

# Generate Fig S2d ####
tmp <- list(apply(out_aging_vs_blood$JI.mat, 2, max), apply(out_aging_vs_lung$JI.mat, 2, max), apply(out_blood_vs_lung$JI.mat, 1, max))
par(mar = c(10, 5, 5, 5))
boxplot(tmp, ylab = "Max matchSCore", names = c("Aging vs MCA blood", "Aging vs MCA lung", "MCA lung vs MCA blood"), las = 2)
