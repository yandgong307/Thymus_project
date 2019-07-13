#######################################################################################################################
#####################split epithelium cells from mesenchymal cells ############################################
#######################################################################################################################

options(stringsAsFactors = F)
rm(list=ls())
load('part1.rds')

#
mes_cell <- rownames(all_anno)[all_anno$cluster=='Mesenchyma']
mes_rawdata <- rawdata[,mes_cell]
mes_anno <- all_anno[mes_cell,'time',drop=F]

#
library(Seurat)
mes_anno$nUMI <- colSums(mes_rawdata)
mes_anno$nGene <- colSums(mes_rawdata>0)

mes <- CreateSeuratObject(mes_rawdata, meta.data = mes_anno)
mes <- NormalizeData(mes, scale.factor = 1e4)
mes <- ScaleData(mes, vars.to.regress = c("nUMI",'nGene','time'))
mes <- FindVariableGenes(mes, x.low.cutoff = 0.1, y.cutoff = 0.5)
mes <- RunPCA(mes, pcs.compute = 50, pc.genes = mes@var.genes)
PCElbowPlot(mes, num.pc = 50)

mes <- RunUMAP(mes, dims.use = 1:30,min_dist = 0.3, spread=0.5)
DimPlot(mes, group.by = 'time', reduction.use = 'umap')

mes <- FindClusters(mes, dims.use = 1:30, resolution = seq(0.2,2,0.2))
DimPlot(mes, group.by = 'res.0.6', reduction.use = 'umap',do.label = T)
FeaturePlot(mes, reduction.use = 'umap', features.plot = c('EPCAM','KRT8','PDGFRA','KRT18'))

#identify epi
epi_cell <- rownames(mes@meta.data)[mes@meta.data$res.0.6=='10']#cluster 10 express EPCAM/KRT8 indicated Epithelium identity
