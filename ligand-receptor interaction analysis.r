#######################################################################################################################################################################################
########################################################################ligand-receptor interaction analysis ##########################################################################
#######################################################################################################################################################################################
get_enrich_matrix <- function(data1, data2, anno1, anno2, org){
  #input data should be normalize
  #celltype combination
  ctp1 <- unique(anno1[,1])
  ctp2 <- unique(anno2[,1])
  all_combination <- c()
  
  for(i in 1:length(ctp2)){
    a <- paste(ctp1, ctp2[i],sep="/")
    b <- paste(ctp2[i],ctp1,sep="/")
    all_combination <- c(all_combination,a,b)
  }
  
  # all l-r pair
  co_gene <- intersect(rownames(data1), rownames(data2))
  l_r_pair<- read.csv("/data2/gyd/data/human_ligand_receptor_paire_melted.csv", row.names = 1)$Pair
  
  com <- c()
  for(i in 1:length(all_combination)){
    tmp <- paste(l_r_pair,all_combination[i],sep=".")
    com <- c(com, tmp)
  }
  
  library(reshape2)
  com_data <- colsplit(com,"\\.",names = c('pair','celltype_combination'))
  
  lr_data <- colsplit(com_data$pair,"_",names = c('ligand','receptor'))
  celltype_data <- colsplit(com_data$celltype_combination,"/",names = c('celltype1','celltype2'))
  com_data <- data.frame(cbind(com_data, lr_data, celltype_data))
  com_data <- com_data[com_data$ligand %in% co_gene &com_data$receptor %in% co_gene,]
  # compute the interaction strength and enrichment pvalue 
  do_enrich <- function(m){
    
    get_enrich_pvalue <- function(data, anno, gene, cluster){
      
      #
      expr1 <- data[gene, rownames(anno)[anno[,1]==cluster]]
      
      expr2 <- c()
      
      for(i in 1:100){
        
        group1 <- rownames(anno)[anno[,1]==cluster]
        
        expr2[i] <- round(mean(data[gene, sample(colnames(data), length(group1))]),2)
      }
      
      #calculate the pvalue, use wilcox 
      
      pvalue <- wilcox.test(expr1, expr2, alternative = 'g')$p.value
      
      log2FC <- log2(mean(expr1)/(mean(expr2)+1e-10))
      
      if(log2FC < 0.5){
        
        pvalue<-NA
      }
      
      return(pvalue) 
    }  
    
    l_gene <- m[3]
    r_gene <- m[4]
    
    celltype1 <- m[5]
    celltype2 <- m[6]
    
    if(celltype1 %in% ctp1){
      
      data1.use <- data1;anno1.use <- anno1
      data2.use <- data2;anno2.use <- anno2
    }else if (celltype1 %in% ctp2){
      
      data1.use <- data2;anno1.use <- anno2
      data2.use <- data1;anno2.use <- anno1
    }else{
      print('expected celltype is imported')
    }  
    
    pvalue1 <- get_enrich_pvalue(data = data1.use, anno = anno1.use, gene = l_gene,cluster = celltype1)
    pvalue2 <- get_enrich_pvalue(data = data2.use, anno = anno2.use, gene = r_gene,cluster = celltype2)
    enrich_pvalue <- pvalue1*pvalue2
    
    pro1 <- round(rowSums(data1.use[l_gene,rownames(anno1.use)[anno1.use[,1]==celltype1],drop=F]>1)/length(rownames(anno1.use)[anno1.use[,1]==celltype1]),2)
    mean1 <- round(rowMeans(data1.use[l_gene,rownames(anno1.use)[anno1.use[,1]==celltype1],drop=F]),2) 
    score1 <- round(pro1*mean1,2)
    
    pro2 <- round(rowSums(data2.use[r_gene,rownames(anno2.use)[anno2.use[,1]==celltype2],drop=F]>1)/length(rownames(anno2.use)[anno2.use[,1]==celltype2]),2)
    mean2 <- round(rowMeans(data2.use[r_gene,rownames(anno2.use)[anno2.use[,1]==celltype2],drop=F]),2) 
    score2 <- round(pro2*mean2,2)
    
    result <- c(m[1:6],pvalue1, pvalue2, enrich_pvalue, pro1, mean1, score1, pro2, mean2, score2)      
    return(result)
  }
  
  enrich_matrix <-apply(com_data,1,function(x) do_enrich(x))
  return(enrich_matrix)
}
###间质上皮与造血的互作
load('/data3/gyd/tmp_file/blood.rds')
load('/data3/gyd/tmp_file/tec.rds')
load('/data3/gyd/tmp_file/mes.rds')
load('/data3/gyd/tmp_file/endo.rds')

###
table(blood_anno$cluster,blood_anno$time)
blood_anno <- blood_anno[,'cluster',drop=F]
blood_data <- blood_data[,rownames(blood_anno)]

tec_anno <- tec_annot[tec_annot$cluster %in% c('TEC1','TEC2','TEC3'),'cluster',drop=F]
tec_anno <- tec_anno[,'cluster',drop=F]
tec_data <- tec_data[,rownames(tec_anno)]

#mes_anno <- mes_anno[!as.character(mes_anno$time) =='w8','cluster',drop=F]
mes_anno$cluster <- 'Mesenchyme'
mes_anno <- mes_anno[,'cluster',drop=F]
mes_data <- mes_data[,rownames(mes_anno)]

#endo_anno <- endo_anno[!as.character(endo_anno$time)=='w8',]
endo_data <- apply(endo_rawdata,2,function(x) log2(1e4*x/sum(x)+1))
endo_anno <- endo_anno[,'cluster',drop=F]
endo_data <- endo_data[,rownames(endo_anno)]

#blood and mes
lr_result1 <- get_enrich_matrix(data1 = blood_data, data2 = mes_data,
                                anno1=blood_anno, anno2 = mes_anno,org = 'hsa')

#blood and epi
lr_result2 <- get_enrich_matrix(data1 = blood_data, data2 = tec_data,
                                anno1=blood_anno, anno2 = tec_anno,org = 'hsa')

#mes and epi
lr_result3 <- get_enrich_matrix(data1 = mes_data, data2 = tec_data,
                                anno1=mes_anno, anno2 = tec_anno,org = 'hsa')

#endo and blood
lr_result4 <- get_enrich_matrix(data1 = blood_data, data2 = endo_data,
                                anno1=blood_anno, anno2 = endo_anno,org = 'hsa')

#endo and epi
lr_result5 <- get_enrich_matrix(data1 = tec_data, data2 = endo_data,
                                anno1=tec_anno, anno2 = endo_anno,org = 'hsa')

#endo and mes
lr_result6 <- get_enrich_matrix(data1 = mes_data, data2 = endo_data,
                                anno1=mes_anno, anno2 = endo_anno,org = 'hsa')




