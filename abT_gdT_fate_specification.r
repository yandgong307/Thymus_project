########################################################################################################################################################################
#######################################################################abT and gdT specification #######################################################################
########################################################################################################################################################################
rm(list=ls())
load('human_thy_blood_data.Rdata')
blood_anno <- blood_anno[blood_anno$cluster %in% as.character(2:6),]
blood_rawdata <- rawdata[,rownames(blood_anno)]

###select gene with high PCA loadings for Monocle analysis
library(Seurat)
tmp <- CreateSeuratObject(blood_rawdata, meta.data = blood_anno)
tmp <- NormalizeData(tmp)
tmp <- ScaleData(tmp, vars.to.regress = c('time'))
tmp <- FindVariableGenes(tmp, x.low.cutoff = 0.0125, y.cutoff = 1)
tmp <- RunPCA(tmp, pcs.compute = 100, pc.genes = tmp@var.genes )
tmp <- ProjectPCA(tmp)


PCElbowPlot(tmp, num.pc = 50)
PCHeatmap(tmp, num.genes = 20,pc.use = 1:5,use.full = F)


library(org.Hs.eg.db)
gogenes <- select(org.Hs.eg.db, keys = c("GO:0007049"), columns = c("SYMBOL"), keytype = "GOALL")
pca_data <- data.frame(tmp@dr$pca@gene.loadings[,1:3])

gene1 <- c(rownames(pca_data)[order(pca_data$PC2,decreasing =T)][1:50],
          rownames(pca_data)[order(pca_data$PC2,decreasing =F)][1:50],
          rownames(pca_data)[order(pca_data$PC3,decreasing =T)][1:50],
          rownames(pca_data)[order(pca_data$PC3,decreasing =F)][1:50])


###perform Monocle analysis
monocle_data <- as.matrix(tmp@data)[,rownames(blood_anno)]
library('monocle')
cell_cycle = unique(gogenes$SYMBOL)
rawdata = blood_rawdata
anno = data.frame(blood_anno,t(monocle_data))

pd <- new("AnnotatedDataFrame", data = anno)
monocle1 <- newCellDataSet(as.matrix(rawdata), phenoData = pd,expressionFamily=negbinomial())
monocle1 <- estimateSizeFactors(monocle1)
monocle1 <- estimateDispersions(monocle1)
monocle1<- setOrderingFilter(monocle1, setdiff(gene1, cell_cycle))
monocle1 <- reduceDimension(monocle1, max_components=2,param.gamma=60)
monocle1 <- orderCells(monocle1)

###define ETP as root of pseudotime 
monocle1 <- reduceDimension(monocle1, max_components=2,param.gamma=100)
monocle1 <- orderCells(monocle1, reverse=FALSE,root_state = '7')


###BEAM analysis to detect genes playing pivital role in abT and gdT fate decision
BEAM_res <- BEAM(monocle1, branch_point=1, cores = 10)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("pval", "qval")]

dran <- plot_genes_branched_heatmap(monocle1[row.names(subset(BEAM_res, qval < 1e-7)),],
                                 branch_point = 1, num_clusters = 6,
                                 cores = 10, use_gene_short_name = T, show_rownames = T,return_heatmap = T)



###plot genes along pseudotime
plot_genes_branched_pseudotime(monocle1[c('SELL','JCHAIN','MKI67','HMGB2','MCM7','CENPM','IGLL1','IFITM1','TRAC','CD1A'),],
                               branch_point=1,
                               color_by="cluster",
                               ncol=2)+scale_color_manual(values = col)

