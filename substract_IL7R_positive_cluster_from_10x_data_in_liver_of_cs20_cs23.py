################################################################################################################################################################
#####################################Detect and substract IL7R+ clusters from 10x data in Liver of CS20 and CS23################################################
################################################################################################################################################################
###scanpy version 1.3.7

#load packages needed
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy.api as sc
import os
from scipy.sparse import *


sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=200)  # low dpi (dots per inch) yields small inline figures
sc.settings.n_jobs=1
sc.logging.print_versions()

#read data
path= '/data2/gyd/aboy_new_world/anlysis/human_thy/human_thy_v3/FL_10x/cs20_cs23_combine/'
os.chdir('/data2/gyd/aboy_new_world/anlysis/human_thy/human_thy_v3/FL_10x/cs20_cs23_combine/')
adata = sc.read(path + 'matrix.mtx').T  # transpose the data
adata.var_names = pd.read_csv(path + 'genes.tsv', header=None, sep='\t')[1]
adata.obs_names = pd.read_csv(path + 'barcodes.tsv', header=None)[0]

# read cell annotation
anno =pd.read_csv('anno.txt',header=0, index_col=0)
adata.obs['time']= anno.time

#plot gene-UMI and qc 
adata.obs['n_genes'] = (adata.X>0).sum(axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
sc.pl.scatter(adata, x='n_genes',y='n_counts',color='PTPRC',)

#filter cells
sc.pp.filter_cells(adata, min_genes=1000)
sc.pp.filter_cells(adata, max_genes=6000)
sc.pp.filter_genes(adata, min_cells=5)

#normalize data
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)


#Find hvgs
sc.pp.highly_variable_genes(adata,min_mean=0.0125, max_mean=10,min_disp=0.5)
sc.pl.highly_variable_genes(adata)
np.sum(adata.var.highly_variable)
adata.raw =adata.copy()

#regress out n_count
sc.pp.regress_out(adata, ['n_counts','n_genes'])
sc.pp.scale(adata, max_value=10)

#pca
sc.tl.pca(adata,use_highly_variable=True)
sc.pl.pca_variance_ratio(adata, log=True)
sc.pl.pca_loadings(adata)


#Find neighbors and Run UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata)


#Detect clusters use louvain
sc.tl.louvain(adata,resolution=0.4) #cluster 6 high express IL7R, CD34, and CD34, will be substracted later


#Define colors
from matplotlib import colors
colors2 = pl.cm.Reds(np.linspace(0, 1, 128))
colors3 = pl.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)


#visualization on umap plot
sc.pl.umap(adata, color=['louvain','time'],legend_loc='on data')
sc.pl.umap(adata, color=['PTPRC','SPINK2','PRSS57','ADGRG1','PROCR','IL7R','IGLL1','CD34','CD7','CD38','TNFRSF18','CD3E','TRAC','TRDC','KIT'],ncols=5,color_map=mymap,use_raw=True)


#save data
adata.obs.to_csv('cs20_cs23_anno.csv')
adata.write('cs20_cs23_umap.h5ad')
pd.DataFrame(adata.obsm['X_umap']).to_csv('cs20_cs23_umap.csv')

