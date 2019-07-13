
####load packages

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy.api as sc
import os
from scipy.sparse import *

sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100)  # low dpi (dots per inch) yields small inline figures
sc.logging.print_versions()

pca_file = './write/pca.h5ad'
umap_file = './write/umap.h5ad'
tsne_file = './write/tsne.h5ad'
spring_file = './write/spring.h5ad'
dpt_file = './write/dpt.h5ad'


#read data

path = 'filter_data/'
adata = sc.read(path + 'matrix.mtx', cache=True).T  # transpose the data
adata.var_names = pd.read_csv(path + 'genes.tsv', header=None, sep='\t')[1]
adata.obs_names = pd.read_csv(path + 'barcodes.tsv', header=None)[0]


#load annotation
anno=pd.read_csv('anno.txt',index_col=0,header=0,sep=",")
adata.obs['time']= anno['time']

#basci filter
sc.pp.filter_cells(adata, min_genes=1000)#will compute the ngenes
sc.pp.filter_genes(adata, min_cells=3)


# compute n_genes, n_counts
mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse to transform to a dense array after summing
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
#add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

sc.pl.scatter(adata, x='n_counts', y='percent_mito')
sc.pl.scatter(adata, x='n_counts', y='n_genes')


#keep rawdata
adata.raw = sc.pp.log1p(adata, copy=True)

#normalizedata
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
adata_normal=adata.copy()

#get hvgs
filter_result = sc.pp.filter_genes_dispersion(
    adata.X, min_mean=0.0125, min_disp=0.5,max_mean=10)
sc.pl.filter_genes_dispersion(filter_result,save="hvg.png")
np.sum(filter_result.gene_subset)


#regress out n_counts and scale data
sc.pp.log1p(adata)
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)

adata = adata[:,filter_result.gene_subset]


#perfrom pca
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True,save="pca.screeplot.png")
sc.pl.pca_loadings(adata,save="pca.loadings.png")
adata.write(pca_file)

#find neighbors and RunUmap
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata,min_dist=0.6,spread=1,maxiter=50)

#Find clusters using louvain
sc.tl.louvain(adata,resolution=0.02)


#visulization on umap plot
sc.pl.umap(adata, color=['time','louvain'],edges=False,size=30)
sc.pl.umap(adata, color=['louvain','PTPRC','CD7','CDH5','SOX7'],edges=False,size=30,color_map='jet')
sc.pl.umap(adata, color=['PDGFRA','PDGFRB','EPCAM','KRT8'],edges=False,size=30,color_map='jet')
sc.pl.umap(adata, color=['EPCAM','KRT5','KRT8','HOXA3'],edges=False,size=30,color_map='bwr')
sc.pl.umap(adata, color=['CD34','CD38','CD7','SELL','CD1A'],edges=False,color_map='bwr',size=30)
sc.pl.umap(adata, color=['TRAC','TRBC2','TRDC','TRGC2'],edges=False,color_map='bwr',size=30)
sc.pl.umap(adata, color=['CD3E','CD4','CD8A','CD8B','MME'],edges=False,color_map='bwr',size=30)

#save data
adata.write(umap_file)
np.savetxt("./write/connection_for_umap.txt", adata.uns['neighbors']['connectivities'].todense())#稀疏矩阵
adata.obsm.to_df()[['X_umap1','X_umap2']].to_csv( './write/umap.csv')
adata.obs[['time','louvain']].to_csv( './scanpy.anno.csv')

