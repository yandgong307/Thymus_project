#######################################################################################################################################################################################
################################################### joint analysis of ETP in thymus and ILR+ clusters in AGM, Liver, Blood############################################################
#######################################################################################################################################################################################

setwd('F:/analysis/project/human_thy_project/part3_progenitor_comparasion/')
rm(list= ls())
options(stringsAsFactors = F)

###load data
load('FL_progenitor/fl_progenitor.rds') #IL7R+ cluster in Liver of CS20 and CS23
load('THY_progenitor/thy_etp.rds') #ETP and CTP in thymus
load('F:/analysis/project/human_thy_project/part3_progenitor_comparasion/Fig4_v5/get_cd127_data/ss_cd127.rds')# IL7R+ cluster in AGM, Liver, and Blood


###rename ITP as committed T cell precursors
etp_anno$cluster <- as.character(etp_anno$cluster)  
etp_anno$cluster[etp_anno$cluster=='intermedial thymus progenitor'] <- 'committed T cell precusor'
names(etp_anno)[2] <- 'site'
etp_anno <- etp_anno[,c('stage','site')]
etp_rawdata <- etp_rawdata[,rownames(etp_anno)]


### map AGM or DA to AGM
select_anno$site[select_anno$site %in% c('AGM','DA')] <- 'AGM'
select_anno$site[select_anno$site =="FL"] <- 'Liver'

###
fl_progenitor_anno$site <- 'Liver'
fl_progenitor_anno <- fl_progenitor_anno[,c('stage',"site")]
fl_progenitor_rawdata <- fl_progenitor_rawdata[,rownames(fl_progenitor_anno)]

###combine data
co_gene <-  intersect(intersect(rownames(select_rawdata), rownames(etp_rawdata)),
                      rownames(fl_progenitor_rawdata)) #26318 gene intersection


rawdata <- cbind(select_rawdata[co_gene,],
                 etp_rawdata[co_gene,],
                 as.matrix(fl_progenitor_rawdata)[co_gene,])

anno <- data.frame(rbind(select_anno[,c('stage',"site")],
                         fl_progenitor_anno,
                         etp_anno) )


###change the annotation
library(Seurat)
anno$stage <- tolower(anno$stage)
anno$cluster <- paste(anno$stage,anno$site, sep="_")
anno$cluster <- plyr::mapvalues(anno$cluster,
                                c("w8_fetal ETP","w8_embryonic ETP",'w9_fetal ETP','w10_fetal ETP',"w9_committed T cell precusor","w10_committed T cell precusor"),
                                c('fetal ETP','embryonic ETP','fetal ETP','fetal ETP','committed T cell precursor','committed T cell precursor'))

anno$stage[anno$stage=='w8'] <- 'w8(56d)'
anno$site[anno$site %in% c('embryonic ETP','fetal ETP','committed T cell precusor')]<-'Thymus'


anno$cluster[anno$stage=='cs12'] <- paste(anno$cluster[anno$stage=='cs12'],' (27d)',sep="")
anno$cluster[anno$stage=='cs13'] <- paste(anno$cluster[anno$stage=='cs13'],' (28d)',sep="")
anno$cluster[anno$stage=='cs14'] <- paste(anno$cluster[anno$stage=='cs14'],' (32d)',sep="")
anno$cluster[anno$stage=='cs15'] <- paste(anno$cluster[anno$stage=='cs15'],' (35d)',sep="")
anno$cluster[anno$stage=='cs17'] <- paste(anno$cluster[anno$stage=='cs17'],' (41d)',sep="")
anno$cluster[anno$stage=='cs20'] <- paste(anno$cluster[anno$stage=='cs20'],' (50d)',sep="")
anno$cluster[anno$stage=='cs23'] <- paste(anno$cluster[anno$stage=='cs23'],' (57d)',sep="")

###compute similarity between clusters
get_similarity_heatmap <- function(anno, data, method,label, method_cluster, cell_order){
  
  tmp_data <- data.frame(t(data[,rownames(anno)]),cluster=anno$cluster)
  aggr_data <- t(aggregate(.~cluster, tmp_data, mean))
  cluster_names<- aggr_data[1,]
  aggr_data <- apply(aggr_data[2:nrow(aggr_data),],2,as.numeric)
  colnames(aggr_data) <- cluster_names
  cor_data <- cor(aggr_data, method = method)
  if(!is.null(cell_order)){
    p <- pheatmap::pheatmap(cor_data[cell_order,cell_order],
                            color = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"))(100),display_numbers = T,number_color = 'red3',cluster_rows = F,cluster_cols = F)
  }else{
    p <- pheatmap::pheatmap(cor_data, treeheight_row = 10, treeheight_col = 10,
                            color = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"))(100),labels_row = label,labels_col = label,display_numbers = T,number_color = 'red3',clustering_method = method_cluster)
  }
  
  return(p)
}

###remove cluster with less than 3 cells
tmp <- names(table(anno$cluster))[table(anno$cluster)>=3]
cell <- rownames(anno)[anno$cluster %in% tmp]



rawdata1 <- rawdata[,rownames(anno)]
main <- CreateSeuratObject(rawdata1[,cell], meta.data = anno[cell,])
main <- NormalizeData(main, scale.factor = 1e4)
main <- ScaleData(main, vars.to.regress = c('nUMI'))
main <- FindVariableGenes(main, x.low.cutoff = 0.1, y.cutoff = 0.5)


table(anno$site,anno$stage)
dend1 <- get_similarity_heatmap(anno =anno , data = as.matrix(main@data[rowMeans(as.matrix(main@data))>0.5,]),method = 'spearman',label = NULL,method_cluster = 'ward.D',cell_order = NULL)


### cell numbers in each stage and cluster
tmp_table <- table(anno[cell,,drop=F]$stage,anno[cell,,drop=F]$cluster)
tmp_table <- tmp_table[c('cs14','cs15','cs17','cs20','cs23',"w8(56d)",'w9','w10'),]

library(ggpubr)
ggpubr::ggtexttable(tmp_table[,c(colnames(tmp_table)[-1],colnames(tmp_table)[1])], theme = ttheme('mBlue'))
write.csv(tmp_table, file='F:/analysis/project/human_thy_project/part3_progenitor_comparasion/Fig4_v5/conjoint_analysis/cell_num.csv')

###joint analysis without CTP
anno$cluster[anno$cluster %in% c('embryonic ETP','fetal ETP')]<-'ETP'
anno <- anno[!anno$cluster=='committed T cell precursor',]
main <- CreateSeuratObject(rawdata1[,rownames(anno)], meta.data = anno)
main <- NormalizeData(main, scale.factor = 1e4)

###compute cell cycle score
g1sGene <- toupper(c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"))
g2mGene <- toupper(c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa"))
main <- ScaleData(main,vars.to.regress = c('nUMI',"nGene"))
main <-CellCycleScoring(main, g2m.genes = g2mGene, s.genes = g1sGene, set.ident = T)

###define colors to use
library(ggsci)
col_site <- setNames(pal_nejm("default")(5),c('YS','AGM','Blood','Liver','Thymus'))
col_stage <- setNames(c("#5050FFFF",'magenta','purple3',"#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#1f77b4","#ff7f0e","#2ca02c"),
                      c('cs12','cs13','cs14','cs15','cs17','cs20','cs23','w8(56d)','w9','w10'))



###Find hvgs and perform pca 
main <- FindVariableGenes(main, x.low.cutoff = 0.0125, y.cutoff =0.5 )
main <- RunPCA(main, pc.genes = main@var.genes, pcs.compute = 50)
PCElbowPlot(main,num.pc = 50)
DimPlot(main, reduction.use = 'pca', group.by = 'stage',do.return = T)+
  scale_color_manual(values = mclust::mclust.options()$classPlotColors)


FeaturePlot(main, reduction.use = 'pca', features.plot = c('IL7R','SELL','CD34','IGLL1'),no.legend = F,do.return = T)+scale_color_gradientn("",colours = colorRampPalette(c('white','black'))(100))

seurat_anno <- data.frame(main@meta.data, main@dr$pca@cell.embeddings[,1:2])
seurat_anno <- data.frame(cbind(seurat_anno, t(as.matrix(main@data)[c('IL7R','SELL','CD34','CD7','CD3E'),rownames(seurat_anno)])))

###feature gene expression
col_cluster <- setNames(pal_d3("category20")(14),unique(main@meta.data$cluster))
library(ggpubr)

a <- ggplot(seurat_anno, aes(x=cluster,y=IL7R,fill=factor(cluster,levels = c('cs12_YS (27d)','cs15_AGM (35d)','cs15_Liver (35d)','cs17_YS (41d)','cs17_Blood (41d)','cs17_Liver (41d)','cs20_Liver (50d)','cs23_Liver (57d)','ETP','ITP')
)))+geom_boxplot(scale = 'width')+scale_fill_manual(values = col_cluster)+theme(axis.text.x = element_text(angle = 60,hjust = 1),legend.position = 'none')

b <- ggplot(seurat_anno, aes(x=cluster,y=SELL,fill=factor(cluster,levels = c('cs12_YS (27d)','cs15_AGM (35d)','cs15_Liver (35d)','cs17_YS (41d)','cs17_Blood (41d)','cs17_Liver (41d)','cs20_Liver (50d)','cs23_Liver (57d)','ETP','ITP')
)))+geom_boxplot(scale = 'width')+scale_fill_manual(values = col_cluster)+theme(axis.text.x = element_text(angle = 60,hjust = 1),legend.position = 'none')

c <- ggplot(seurat_anno, aes(x=cluster,y=CD34,fill=factor(cluster,levels = c('cs12_YS (27d)','cs15_AGM (35d)','cs15_Liver (35d)','cs17_YS (41d)','cs17_Blood (41d)','cs17_Liver (41d)','cs20_Liver (50d)','cs23_Liver (57d)','ETP','ITP')
)))+geom_boxplot(scale = 'width')+scale_fill_manual(values = col_cluster)+theme(axis.text.x = element_text(angle = 60,hjust = 1),legend.position = 'none')

d <- ggplot(seurat_anno, aes(x=cluster,y=CD7,fill=factor(cluster,levels = c('cs12_YS (27d)','cs15_AGM (35d)','cs15_Liver (35d)','cs17_YS (41d)','cs17_Blood (41d)','cs17_Liver (41d)','cs20_Liver (50d)','cs23_Liver (57d)','ETP','ITP')
)))+geom_boxplot(scale = 'width')+scale_fill_manual(values = col_cluster)+theme(axis.text.x = element_text(angle = 60,hjust = 1),legend.position = 'none')

ggarrange(a,b,c,d,ncol = 2,nrow = 2)


### project stage, site, and cluster information onto pca plot
p1 <- ggplot(seurat_anno, aes(x=PC1, y=PC2, col=factor(site,levels = c('YS','AGM','Blood','Liver','Thymus'))))+geom_point(size=2)+scale_color_manual("", values =col_site )+theme_bw()+theme(panel.grid = element_blank(),legend.position = c(0.1,0.8))
p2 <- ggplot(seurat_anno, aes(x=PC1, y=PC2, col=factor(stage,levels = c('cs12','cs13','cs14','cs15','cs17','cs20','cs23','w8(56d)','w9','w10'))))+geom_point(size=2)+scale_color_manual("",values = col_stage)+theme_bw()+theme(panel.grid = element_blank(),legend.position = c(0.1,0.8))
p3 <- ggplot(seurat_anno, aes(x=PC1, y=PC2, col=cluster))+geom_point(size=2)+scale_color_manual(values = col_cluster)+theme_bw()+theme(panel.grid = element_blank(),legend.position = c(0.15,0.7))

ggpubr::ggarrange(p1,p2,p3,ncol = 3,nrow = 1)

###Find clusters
main <- FindClusters(main, dims.use = 1:20, resolution = seq(0.1,2,0.1),force.recalc = T)

DimPlot(main,reduction.use ='pca',group.by = 'res.0.3',pt.size = 1.5,do.return = T,label.size = 8)+scale_color_manual(values = mclust::mclust.options()$classPlotColors)+theme(legend.position = c(0.1,0.8))
merge_anno <- data.frame(main@meta.data[,c('stage','site','res.0.3','Phase')], main@dr$pca@cell.embeddings[,1:2])
names(merge_anno)[3] <- 'cluster'
merge_anno$cluster <- plyr::mapvalues(merge_anno$cluster, as.character(0:4), 
                                      c('c1','c2','c3','c4','c5'))

####ETP were devided into 2 compartment, one part merge with IL7R+ clusters in Liver, and the other only cluster with ETP (c4)
table(merge_anno$stage, merge_anno$cluster)
etp_cell <- rownames(etp_anno)[as.character(etp_anno$site) %in% c('embryonic ETP','fetal ETP')]

merge_anno$group <- 'others'
merge_anno$group[!merge_anno$cluster=='c4' & rownames(merge_anno) %in% etp_cell] <- 'TSP-like ETP'
merge_anno$group[merge_anno$cluster=='c4' & rownames(merge_anno) %in% etp_cell] <- 'proliferating ETP'

merge_anno$Phase[!merge_anno$site=='Thymus']<-'others'
merge_anno$Phase <- factor(merge_anno$Phase, levels = c('G1','S','G2M','others'))


col_etp <- setNames(c('blue','red',scales::alpha('grey',alpha = 0.1)),c('TSP-like ETP','proliferating ETP','others'))
merge_anno$group <- factor(merge_anno$group, levels = c('TSP-like ETP','proliferating ETP','others'))
col_phase <- c('G1'='blue3','S'='green3','G2M'='red3','others'=scales::alpha('grey',alpha = 0.1))

f <- merge_anno %>% ggplot(aes(PC1,PC2,col=Phase))+geom_point(size=1.5)+scale_color_manual(values = col_phase)+theme_bw()+theme(legend.position = c(0.15,0.8),panel.grid = element_blank())
g <- merge_anno %>% ggplot(aes(PC1,PC2,col=group))+geom_point(size=1.5)+scale_color_manual(values = col_etp)+theme_bw()+theme(legend.position = c(0.2,0.8),panel.grid = element_blank())

ggpubr::ggarrange(c,d,f,g,nrow = 2, ncol = 2)

###save data
write.csv(file = 'F:/analysis/project/human_thy_project/part3_progenitor_comparasion/Fig4_v5/conjoint_analysis/pca_data_include_ctp.csv',seurat_anno)
write.csv(file = 'F:/analysis/project/human_thy_project/part3_progenitor_comparasion/Fig4_v5/conjoint_analysis/pca_data_exclude_ctp.csv',merge_anno)

###compare TSP-like ETP and proliferating ETP
etp <- CreateSeuratObject(rawdata[,etp_cell], meta.data = merge_anno[etp_cell,])
etp <- NormalizeData(etp, scale.factor = 1e4)
etp <- ScaleData(etp)
etp <- SetAllIdent(etp, id="group")

deg <- FindAllMarkers(etp, logfc.threshold = 0.25, test.use = 'wilcox', only.pos = T)
deg <- deg[deg$p_val_adj <0.05,]

write.csv(file = 'F:/analysis/project/human_thy_project/part3_progenitor_comparasion/Fig4_v5/deg_different_etp.csv',deg)

###pseudotime analysis combined IL7R+ cluster in Liver of CS20 and CS23 with ETP and CTP in thymus
monocle_cell <- rownames(seurat_anno)[seurat_anno$cluster %in% c('cs20_Liver (50d)','cs23_Liver (57d)','ETP','committed T cell precursor')]
monocle_rawdata <- rawdata[,monocle_cell]
monocle_anno <- data.frame(seurat_anno[monocle_cell,],t(as.matrix(main@data[,monocle_cell])))

library(monocle)#version 2.8.0
monocle_anno$cluster[rownames(monocle_anno) %in% rownames(merge_anno)[as.character(merge_anno$group)=='TSP-like ETP']]<- 'TSP-like ETP'
monocle_anno$cluster[rownames(monocle_anno) %in% rownames(merge_anno)[as.character(merge_anno$group)=='proliferating ETP']]<- 'proliferating ETP'

### select genes with large PCA loading
tmp <- SubsetData(main, cells.use = monocle_cell)
tmp <- ScaleData(tmp, vars.to.regress = c('nGene','nUMI'))
tmp <- FindVariableGenes(tmp, x.low.cutoff = 0.0125,y.cutoff = 0.5)
tmp <- RunPCA(tmp, pc.genes = tmp@var.genes)

PCHeatmap(tmp, pc.use = 1:5,num.genes = 10)
pca_data <- data.frame(tmp@dr$pca@gene.loadings[,1:3])

gene1 <- c(rownames(pca_data)[order(pca_data$PC1,decreasing =T)][1:50],
           rownames(pca_data)[order(pca_data$PC1,decreasing =F)][1:50],
           rownames(pca_data)[order(pca_data$PC3,decreasing =T)][1:50],
           rownames(pca_data)[order(pca_data$PC3,decreasing =F)][1:50])#PC2 mainly comprised of cell cycle genes, which is exclude 


pd <- new("AnnotatedDataFrame", data =monocle_anno)
HSMM <- newCellDataSet(as.matrix(monocle_rawdata), phenoData = pd,expressionFamily=negbinomial())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
#
disp_table <- dispersionTable(HSMM)
# ordering_genes <- subset(disp_table,
#                              mean_expression >= 0.5 &
#                                dispersion_empirical >= 2* dispersion_fit)$gene_id

HSMM<- setOrderingFilter(HSMM, gene1)
HSMM <- reduceDimension(HSMM, max_components=2,param.gamma=60)
HSMM <- orderCells(HSMM)
#HSMM <- orderCells(HSMM,root_state = '5')
plot_cell_trajectory(HSMM,color_by = "cluster",cell_size = 1.5,show_branch_points = F)+scale_color_manual(values=c(col_cluster,col_etp))
plot_cell_trajectory(HSMM,color_by = "site",cell_size = 1.5,show_branch_points = F)+scale_color_manual(values=c(col_site))
plot_cell_trajectory(HSMM,color_by = "CD34",cell_size = 1.5)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
plot_cell_trajectory(HSMM,color_by = "MEF2C",cell_size = 1.5)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
plot_cell_trajectory(HSMM,color_by = "CD74",cell_size = 1.5)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
plot_cell_trajectory(HSMM,color_by = "LY86",cell_size = 1.5)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
plot_cell_trajectory(HSMM,color_by = "IL7R",cell_size = 1.5)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
plot_cell_trajectory(HSMM,color_by = "SPINK2",cell_size = 1.5)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
plot_cell_trajectory(HSMM,color_by = "IL6R",cell_size = 1.5)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
plot_cell_trajectory(HSMM,color_by = "PROM1",cell_size = 1.5)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
plot_cell_trajectory(HSMM,color_by = "IRF8",cell_size = 1.5)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
plot_cell_trajectory(HSMM,color_by = "SPINK2",cell_size = 1.5)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
plot_cell_trajectory(HSMM,color_by = "PRSS57",cell_size = 1.5)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))


###
plot_cell_trajectory(HSMM,color_by = "Pseudotime",cell_size = 1.5,show_tree = T,show_branch_points = F)+scale_color_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(10))
plot_cell_trajectory(HSMM,color_by = "State",cell_size = 1.5,show_tree = T,show_branch_points = F)+scale_color_manual(values = mclust::mclust.options()$classPlotColors)

###IL7R+ clusters in Liver were devided into two part(state 1 and 3),cells merge with TSP-like ETP in thymus were identified as TSP, then DEGs were calculated
monocle_tmp <- pData(HSMM)
ana_cell <- rownames(monocle_tmp)[as.character(monocle_tmp$State) %in% c('1','3')]

pro <- CreateSeuratObject(monocle_rawdata[,ana_cell], meta.data = monocle_tmp[ana_cell,])
pro <- NormalizeData(pro, scale.factor = 1e4)
pro <- ScaleData(pro)

pro <- SetAllIdent(pro, id='State')
deg_ana <- FindAllMarkers(pro, logfc.threshold = 0.25,only.pos = T,min.pct = 0.2)
deg_ana <- deg_ana[deg_ana$p_val_adj<0.05,]


a <- plot_cell_trajectory(HSMM,color_by = "TYROBP",cell_size = 1.5,show_branch_points = F,show_tree = T)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
b <- plot_cell_trajectory(HSMM,color_by = "IL3RA",cell_size = 1.5,show_branch_points = F,show_tree = T)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
c <- plot_cell_trajectory(HSMM,color_by = "IRF8",cell_size = 1.5,show_branch_points = F,show_tree = T)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))

ggpubr::ggarrange(a,b,c,ncol = 3,nrow = 1)

e <- plot_cell_trajectory(HSMM,color_by = "CD79A",show_tree = T,cell_size = 1.5,show_branch_points = F)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
f <- plot_cell_trajectory(HSMM,color_by = "CD79B",show_tree = T,cell_size = 1.5,show_branch_points = F)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))
g <- plot_cell_trajectory(HSMM,color_by = "VPREB3",show_tree = T,cell_size = 1.5,show_branch_points = F)+scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(10))

ggpubr::ggarrange(a,b,c,e,f,g,ncol = 3,nrow = 2)

write.csv(file = 'F:/analysis/project/human_thy_project/part3_progenitor_comparasion/Fig4_v5/conjoint_analysis/deg_of_branch1_and_brach2.csv',deg_ana)

##################################################################################################################################################
################################DEG, IL7R+ cluster comparasion in YS, AGM, FL, and Thymus#########################################################
##################################################################################################################################################
###IL7R+ cluster in Liver were devide into 2 part for comparsion, one id TSP in Liver of CS20 and CS23, the others is IL7R+ cluster in Liver of CS15 and CS17

tsp <- rownames(monocle_anno)[as.character(monocle_anno$State)=='3' & monocle_anno$site=='Liver']
liver_cell2 <- rownames(merge_anno)[merge_anno$site=='Liver' & merge_anno$stage %in% c('cs15','cs17')]
other_progenitor <- rownames(merge_anno)[!merge_anno$site=='Liver']
anno <- merge_anno[c(tsp,liver_cell2, other_progenitor),]


anno$Site <- anno$site
anno$Site[rownames(anno) %in% tsp] <- 'Liver (cs20,cs23)'
anno$Site[rownames(anno) %in% liver_cell2] <- 'Liver (cs15,cs17)'

###calculate DEGs
expr <- rawdata [,rownames(anno)] 
pro <- CreateSeuratObject(expr, meta.data = anno)
pro <- NormalizeData(pro)
pro <- ScaleData(pro)
pro <- FindVariableGenes(pro, x.low.cutoff = 0.1, y.cutoff = 0.5)
#deg of site
pro <- SetAllIdent(pro, id='Site')
deg_site <- FindAllMarkers(pro, logfc.threshold = 0.5, only.pos = T ,test.use = 'wilcox',min.pct = 0.2)
deg_site <- deg_site[deg_site$p_val_adj< 0.05,]

write.csv(file = 'Fig4_v5/conjoint_analysis/deg_progenitor_in_different_site.csv',deg_site)


### combined analysis including cells, ETP(thymus), TSP and IL7R+ cluster in AGM ,Liver(cs15,cs17), Blood (cs15 and cs17)
tsp <- rownames(monocle_anno)[as.character(monocle_anno$State)=='3' & monocle_anno$site=='Liver']
etp <- rownames(anno)[anno$cluster=='ETP']
other_progenitor <- rownames(anno)[anno$cluster %in% unique(anno$cluster)[1:6]]

rawdata.use <- rawdata[,c(tsp,etp,other_progenitor)]
anno.use <- anno[colnames(rawdata.use),]

###
main <- CreateSeuratObject(rawdata.use, meta.data = anno.use)
main <- NormalizeData(main, scale.factor = 1e4)
main <-CellCycleScoring(main, g2m.genes = g2mGene, s.genes = g1sGene, set.ident = T)
main <- ScaleData(main, vars.to.regress = c('nUMI','nGene','S.Score','G2M.Score'))
main <- FindVariableGenes(main, x.low.cutoff = 0.5, y.cutoff = 0.5)
main <- RunPCA(main, pc.genes = main@var.genes, pcs.compute = 50)
PCElbowPlot(main, num.pc = 50)

main <- RunUMAP(main, dims.use = 1:10)
DimPlot(main, reduction.use = 'umap',group.by = 'site')
DimPlot(main, reduction.use = 'umap',group.by = 'stage')

###filter out outliers
tmpcell <- c("TTCTTAGAGAAGGTTT-2","ACCCACTTCTGCGTAA-3",
             "TGCGTGGAGACAAAGG-3","GGGCACTGTAGCGATG-1",
             "ACATCAGCACAAGCCC-2","CACAAACGTCCCTACT-2",
             "GACCTGGGTCAGATAA-2","GTGGGTCTCTATCCCG-2")
main <- SubsetData(main, cells.use = setdiff(rownames(anno.use),tmpcell))
main <- FindVariableGenes(main, x.low.cutoff = 0.0125, y.cutoff = 0.5)
main <- RunPCA(main, pc.genes = main@var.genes, pcs.compute = 50)
PCElbowPlot(main, num.pc = 50)

#
main <- RunUMAP(main, dims.use = 1:20)
DimPlot(main, reduction.use = 'umap',group.by = 'site',do.return = T,pt.size = 2)+scale_color_manual(values = col_site)

p1 <- DimPlot(main,pt.size = 2, reduction.use = 'umap',group.by = 'site',do.return = T)+scale_color_manual(values = col_site)+theme(legend.position = c(0.8,0.2))
p2 <- DimPlot(main,pt.size = 2, reduction.use = 'umap',group.by = 'stage',do.return = T)+scale_color_manual(values = col_stage)+theme(legend.position = c(0.8,0.2))
plot_grid(p1,p2)

FeaturePlot(main, reduction.use = 'umap',features.plot = c('CD34','IL7R','CD7','LEF1'),nCol = 2,cols.use = c('grey','red'),pt.size = 1.5)
FeaturePlot(main, reduction.use = 'umap',features.plot = c('CD34','IL7R','CD7'),cols.use = c('grey','red'),pt.size = 1.5,nCol = 3)

###
agm_anno <- data.frame(main@meta.data, main@dr$umap@cell.embeddings[])

projection_on_all <- function(select, anno, col,dm){
  
  tmp_anno <- anno
  tmp_anno$site[!tmp_anno$stage==select] <- 'others'
  tmp_anno <- data.frame(tmp_anno, dm)
  
  p <- ggplot(tmp_anno, aes(x=UMAP1, y=UMAP2, col=site))+geom_point()+scale_color_manual(values = c(col, 'others'='grey'))+theme(legend.position = 'none')
  return(p)
}



p1 <- projection_on_all(select = 'cs14', anno = agm_anno, col = col_site,dm = agm_anno[,c('UMAP1','UMAP2')])
p2 <- projection_on_all(select = 'cs15', anno = agm_anno, col = col_site,dm = agm_anno[,c('UMAP1','UMAP2')])
p3 <- projection_on_all(select = 'cs17', anno = agm_anno, col = col_site,dm = agm_anno[,c('UMAP1','UMAP2')])
plot_grid(p1,p2,p3,ncol = 3)

save(file = 'F:/analysis/project/human_thy_project/part3_progenitor_comparasion/Fig4_v5/agm_orign/agm_orgin_il7R_only.rds',
     rawdata.use,anno.use,agm_anno,col_stage,col_site)


