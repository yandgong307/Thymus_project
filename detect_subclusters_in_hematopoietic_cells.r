#######################################################################################################################
############################### Detect subclusters in hematopoietic cells #############################################
#######################################################################################################################

rm(list=ls())
load("data.Rdata")

###Detect subclusters use Seurat (version 2.3.4)
library(Seurat)
tmp <- CreateSeuratObject(rawdata, meta.data = anno)
tmp <- NormalizeData(tmp)
tmp <- ScaleData(tmp, vars.to.regress = c('time'))
tmp <- FindVariableGenes(tmp, x.low.cutoff = 0.0125, y.cutoff = 1)

library(org.Hs.eg.db)
gogenes <- select(org.Hs.eg.db, keys = c("GO:0007049"), columns = c("SYMBOL"), keytype = "GOALL")
tmp@var.genes <- setdiff(tmp@var.genes,unique(gogenes$SYMBOL))


tmp <- RunPCA(tmp, pcs.compute = 50, pc.genes = tmp@var.genes )
PCElbowPlot(tmp, num.pc = 50)
tmp <- FindClusters(tmp, dims.use = 1:30, resolution =seq(0.1,2,0.1),
                    temp.file.location = "/data2/gyd/", k.param = 15)

tmp@meta.data$cluster <- tmp@meta.data$res.0.3 #res.0.3 is selected
tmp@meta.data$cluster[tmp@meta.data$cluster=="0"] <- '10'
tmp <- RunTSNE(tmp, dims.use = 1:30, seed.use = 2)


col <- setNames(mclust::mclust.options()$classPlotColors[1:10],unique(tmp@meta.data$cluster))
col_time <- setNames(c('#1f77b4', '#ff7f0e', '#2ca02c'),c('week8','week 9',"week 11"))


cluster <- DimPlot(tmp, reduction.use = 'tsne', do.return = T,group.by = 'cluster', do.label = T,pt.size = 1)+
  scale_color_manual(values=col)+theme_nothing()

time <- DimPlot(tmp, reduction.use = 'tsne', do.return = T,group.by = 'time', do.label = T,pt.size = 1)+theme_nothing()+
  scale_color_manual('time',values = col_time)

plot_grid(cluster, time)

FeatureHeatmap(tmp, features.plot = c('CD34','CD7','SELL'), pt.size = 0.8,
               group.by = 'time', do.return = T, plot.horiz = T)+
               scale_color_gradientn(colours = colorRampPalette(c("lightgrey",'blue','red'))(100))

FeaturePlot(tmp, features.plot = c('CD34','CD7','SELL') ,cols.use = c('lightgrey','blue'),nCol=3, no.axes = T)
FeaturePlot(tmp, features.plot = c('SPINK2','PRSS57','IGLL1') ,cols.use = c('lightgrey','blue'),nCol=3, no.axes = T)
FeaturePlot(tmp, features.plot = c('CD63','CD74','CD68') ,cols.use = c('lightgrey','blue'),nCol=3, no.axes = T)


###Detect subclusters in ILC

ilc <- SubsetData(tmp ,cells.use = rownames(anno)[anno$cluster=="4"])

ilc <- FindVariableGenes(ilc, x.low.cutoff = 1, y.cutoff = 1)

ilc <- RunPCA(ilc, pc.genes = ilc@var.genes)

ilc <- FindClusters(ilc, dims.use = 1:10, temp.file.location = '/data2/gyd/', k.param = 10, resolution = seq(0.2,2,0.2),force.recalc = T)

DimPlot(ilc,do.return = T, reduction.use = 'pca',group.by = 'res.0.2' ,dim.1 = 1, dim.2 = 2,pt.size = 3)+
  scale_color_manual(values = col_ilc)

ilc <- SetAllIdent(ilc, id="res.0.2")
library(plyr)

ilc@meta.data$res.0.2 <- mapvalues(ilc@meta.data$res.0.2,c("0","1","2"),c('c1',"c2","c3"))

deg3 <- FindAllMarkers(ilc,logfc.threshold = 0.5, test.use = 'wilcox', only.pos = T)
deg3 <- deg3[deg3$p_val_adj<0.05,]

