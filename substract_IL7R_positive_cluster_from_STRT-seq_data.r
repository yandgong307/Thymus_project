########################################################################################################################################################################################
####################################################Detect and substract IL7R+ cluster from STRT-seq data in AGM, Liver and Blood ######################################################
########################################################################################################################################################################################

rm (list = ls() )
options( stringsAsFactors = F)


load('data.rds')

###qc: nGene 2,000-10,000, nUMI 10000-1,000,000
annot <- annot[(annot$nGene >2000 & annot$nGene <10000 & annot$nUMI>10000 & annot$nUMI<1e6 & annot$spike_in<0.5),]
rawdata <- rawdata [,rownames(annot)]

###remove ERCC
nonERCC <- grep('^ERCC-', rownames(rawdata), value = T, invert = T)
rawdata <- rawdata[nonERCC, ]


library(Seurat)
main <- CreateSeuratObject(raw.data = rawdata, meta.data = annot)
main <- NormalizeData(main, scale.factor = 1e5)
main <- ScaleData(main)
main <- FindVariableGenes(main, x.low.cutoff = 0.5, y.cutoff = 0.5)
main <- RunPCA(main, pc.genes = main@var.genes, pcs.compute = 50)
PCElbowPlot(main,num.pc = 50)

main@meta.data$site[main@meta.data$site=='FL']<- 'Liver'
main@meta.data$site[main@meta.data$site %in% c('CH','DA','AGM')]<- 'AGM'

#umap analysis and find clusters
main <- RunUMAP(main, dims.use = 1:20, min_dist = 0.5 , spread=1)
a <- DimPlot(main, reduction.use = 'umap', group.by = 'site',do.return = T)
b <- DimPlot(main, reduction.use = 'umap', group.by = 'stage',do.return = T)

main <- FindClusters(main, dims.use = 1:30, resolution = seq(0.1,2,0.1),force.recalc = T,k.param = 15)
c <- DimPlot(main, reduction.use = 'umap', group.by = 'res.0.8',do.return = T,pt.size = 1.8)+scale_color_manual(values = mclust::mclust.options()$classPlotColors)
d <- FeaturePlot(main, reduction.use = 'umap', features.plot = c('IL7R','CD7','CD3E','CSF1R','GATA1','GYPA','KLF1','PF4','CPA3','CD14','MPO','CEACAM8','JCHAIN','SELL','CD34','IL1R1'),do.return = T)


plot_grid(a,b,c,d$CD7,d$CD3E,d$IL7R,d$JCHAIN,d$SELL,d$CSF1R,ncol = 3)
plot_grid(a,b,c,d$IL7R,ncol = 2)
plot_grid(a,b,ncol = 2)

#feature genes on umap plot
FeaturePlot(main, reduction.use = 'umap', features.plot = c('IL7R','CD7','RORC','LTA','CD34','SPINK2'),cols.use = c('grey','red'),nCol = 3)
FeaturePlot(main, reduction.use = 'umap', features.plot = c('IL7R','CD7','CD3E','RORC','LTA','MS4A1','CD19','CD79A',"CD79B"))
DimPlot(main, reduction.use = 'umap',do.label = T, group.by = 'res.1',do.return = T,pt.size = 1.6)+scale_color_manual(values = mclust::mclust.options()$classPlotColors)

#####
ftable(table(main@meta.data$res.0.8, main@meta.data$site,main@meta.data$stage))

##select 0 from res.0.8
select_cell <- rownames(main@meta.data)[main@meta.data$res.0.8=='0']#cluster 0 express IL7R ,CD7 and CD34 and will be substracted later
select_rawdata <- rawdata[,select_cell]
select_anno <- annot[select_cell,]

###save data
all_anno <- data.frame(main@meta.data[,c('site','stage','group','res.0.8')],
                    main@dr$umap@cell.embeddings[])
names(all_anno)[4] <- 'cluster'
all_anno$cluster <- plyr::mapvalues(all_anno$cluster, c(0:5),
                         paste('c',1:6,sep=""))

main@meta.data$res.0.8 <- plyr::mapvalues(main@meta.data$res.0.8, c(0:5),
                         paste('c',1:6,sep=""))

library(ggplot2)
library(ggsci)
ggplot(all_anno, aes(x=UMAP1, y=UMAP2, col=cluster))+geom_point(size=1.5)+theme_bw()+theme(legend.position = c(0.9,0.82))+scale_color_aaas()


write.csv(rawdata, file='F:/analysis/project/human_thy_project/part3_progenitor_comparasion/Fig4_v5/get_cd127_data/rawdata_qc_controled.csv')
write.csv(all_anno, file='F:/analysis/project/human_thy_project/part3_progenitor_comparasion/Fig4_v5/get_cd127_data/anno_dm.csv')


save(file='F:/analysis/project/human_thy_project/part3_progenitor_comparasion/Fig4_v5/get_cd127_data/ss_cd127.rds',select_anno, select_rawdata)
save(file='F:/analysis/project/human_thy_project/part3_progenitor_comparasion/Fig4_v5/get_cd127_data/data_to_select_ss_cd127.rds',rawdata, annot, main)




