################################################################################################################################################
########################################################################TEC analysis ###########################################################
################################################################################################################################################
options(stringsAsFactors = F)
rm(list = ls())
load('human_thy_epi_data.Rdata')


### Detect and identify all cell types
library(Seurat)
tec <- CreateSeuratObject(tec_rawdata, meta.data = tec_annot)
tec <- NormalizeData(tec, scale.factor = 1e5)
tec <- ScaleData(tec,vars.to.regress = 'nUMI')
tec <- FindVariableGenes(tec, x.low.cutoff = 0.5, y.cutoff =0.5 )
tec <- RunPCA(tec, pc.genes = tec@var.genes, pcs.compute = 50)
PCElbowPlot(tec,num.pc = 50)

tec <- RunUMAP(tec, dims.use = 1:20, n_neighbors = 15)
tec <- RunTSNE(tec, dims.use = 1:20)
DimPlot(tec,reduction.use = 'umap',group.by = 'stage')

FeaturePlot(tec, reduction.use = 'umap', features.plot = c('PTPRC','EPCAM','CLDN4','SIX1','FOXN1','PDGFRB','PDGFRA','GATA1','HBA1'),nCol = 3)
FeaturePlot(tec, reduction.use = 'umap', features.plot = c('PTPRC','EPCAM','CLDN4','SIX1','FOXN1','LY75','CCL25','CCL21'),nCol = 3)
FeaturePlot(tec, reduction.use = 'umap', features.plot = c('PTPRC','CD3D','LEF1','TCF7','KLF1','HBB','NKG7','CD79A','MAFF'),nCol = 3)

tec <- FindClusters(tec, dims.use = 1:20, resolution = seq(0.2,2,0.2),k.param = 10)
a <- DimPlot(tec,reduction.use ='tsne',group.by = 'res.2',pt.size = 2, do.label = T,do.return = T)+scale_color_manual(values = mclust::mclust.options()$classPlotColors)
b <- DimPlot(tec,reduction.use ='umap',group.by = 'res.2',pt.size = 2, do.label = T,do.return = T)+scale_color_manual(values = mclust::mclust.options()$classPlotColors)
c <- DimPlot(tec,reduction.use ='tsne',group.by = 'stage',pt.size = 2, do.label = T)
plot_grid(a,c,ncol = 2)

tec@meta.data$show <- ifelse(tec@meta.data$res.2 %in% c('1','2','3','4','5','6','7','9'),'tec','others')
DimPlot(tec,reduction.use ='tsne',group.by = 'show',pt.size = 2, do.label = F,do.return = T)+scale_color_manual(values = c('grey','red'))

FeaturePlot(tec, reduction.use = 'tsne', features.plot = c('PTPRC','EPCAM','HLA-DRA','CD74','HOXA3','PBX1','PAX9','PMP2','OLIG1'),nCol = 3)
FeaturePlot(tec, reduction.use = 'tsne', features.plot = c('PTPRC','EPCAM','CLDN4','SIX1','FOXN1','PDGFRB','PDGFRA','GATA1','HBA1'),nCol = 3)
FeaturePlot(tec, reduction.use = 'tsne', features.plot = c('PTPRC','EPCAM','CLDN4','SIX1','FOXN1','LY75','CCL25','CCL21'),nCol = 3)
FeaturePlot(tec, reduction.use = 'tsne', features.plot = c('PTPRC','CD3D','LEF1','TCF7','KLF1','C1QA','NKG7','CD79A','CXCL9'),nCol = 3)
FeaturePlot(tec, reduction.use = 'tsne', features.plot = c('PTPRC','CD3D','LEF1','KLF1','PDGFRA','PDGFRB'),nCol = 3)

### filter out erythroid cells, mesenchymal cells and neron cells
ana_cell <- rownames(tec@meta.data)[tec@meta.data$res.2 %in% c('1','2','3','4','5','6','7','9')]
all_annot <- tec@meta.data[,c('res.2','stage')]
names(all_annot)[1] <- 'cluster'
tec_annot <- all_annot[ana_cell,]
tec_rawdata <- rawdata[,rownames(tec_annot)]


### epithelial cells were subset and reanalyzed
gene <- list(c('HBB','HBA1','HBA2','KLF1','HBG2','GYPA','GYPC','RHOB','RHAG'))
ana <- CreateSeuratObject(tec_rawdata, meta.data = tec_annot)
ana <- NormalizeData(ana, scale.factor = 1e5)
ana <- AddModuleScore(ana, genes.list = gene,enrich.name = 'Erp')
ana <- ScaleData(ana,vars.to.regress =c('nGene','nUMI','Erp1'))

ana <- FindVariableGenes(ana, x.low.cutoff = 0.5, y.cutoff =0.5 )
ana <- RunPCA(ana, pc.genes = ana@var.genes, pcs.compute = 50)
PCElbowPlot(ana,num.pc = 50)



ana <- RunUMAP(ana, dims.use = 1:20,n_neighbors = 20)
DimPlot(ana,reduction.use = 'umap',group.by = 'stage',pt.size = 2)
FeaturePlot(ana, reduction.use = 'umap', features.plot = c('PTPRC','SIX1','HLA-DRA','CD74','CCL25','CLDN4'),nCol = 3)
ana <- FindClusters(ana, dims.use = 1:20, resolution = seq(0.1,2,0.1),k.param = 15)

p1 <- DimPlot(ana,reduction.use = 'umap',group.by = 'res.1.2',pt.size = 2,do.return = T)
p2 <- DimPlot(ana,reduction.use = 'umap',group.by = 'stage',pt.size = 2, do.return = T)
ggpubr::ggarrange(p1,p2)

save(file = 'thy_project_epi_data.Rdata',tec_anno, rawdata, all_annot)

