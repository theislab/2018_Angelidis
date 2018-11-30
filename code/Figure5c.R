# Load R libraries ####
library(Seurat)

# Load Aging Seurat object ####
load("../data/SeuratObject.RData")

# Generate Fig 5c ####
tmp <- seu.ica
tmp <- SetIdent(tmp, ident.use = seu.ica@meta.data$celltype) 
tmp <- SubsetData(tmp, ident.remove = c('Gamma-Delta_T_cells', 'low_quality_cells', 'red_blood_cells'))
DotPlot(tmp, c('Col14a1', 'Dcn'), group.by = 'celltype', x.lab.rot = T, dot.scale = 15, plot.legend = T)
