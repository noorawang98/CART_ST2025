import scanpy as sc
import numpy as np
import pandas as pd 
import bbknn
import scrublet as scr
import os
import scvi
sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor="white")
'''
os.chdir('/home/wr/LJX/code/celescope/shell_rna')
os.listdir('./')
batchname= ['Control','NR','R']
filename = [batch+"/outs/raw" for batch in batchname]
adata_dic = [sc.read_10x_mtx(file) for file in filename]
#-----def functions ----------
def preprocess_step1(adata,batch):
    adata.raw = adata
    adata.layer['counts']=adata.X.copy()
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^Hb[^(P)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=[50], log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=100)
    #sc.pp.filter_genes(adata, min_cells=2)
    adata = adata[adata.obs.pct_counts_mt<5,:]
    sc.external.pp.scrublet(adata)
    adata = adata[adata.obs.predicted_doublet!=True,:]
    sc.pl.violin(adata,["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                 jitter=0.4,multi_panel=True,save=batch+'_qc.pdf')
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",save=batch+'_qc.pdf')
    adata.obs =adata.obs.fillna(0)
    sc.pp.normalize_total(adata,target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.pca(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    adata.write(batch+'.h5ad')
    return(adata)

def preprocess_step2(adata):
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.pca(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    return(adata)
#-----end def functions ----------


#----markers--------
marker_dic={
    'Exogene':['EGFP','mCherry-CAR'],
    'T':['Cd3e','Cd3d','Cd4','Cd8a','Cd8b1','Trdc','Trgc1'],
    'NK':['Klrb1c','Klrk1'],
    'B':['Cd19','Cd79a','Cd79b','Sdc1'],
    'Mono':['Ly6g','Cd14','Csf3r','Adgre1'],
    'Macro':['C1qa','Cd68','Mrc1','Cd163','Cd86','Cd80','Il1b'],
    'Neu':['Csf1r','S100a8','S100a9'],
    'DC':['Cst3','Cd74','Itgax','Siglech'],
    'Fibro':['Col1a1','Dcn','Acta2','Cxcl1','Cxcl10'],
    'EC':['Pecam1','Eng','Vwf'],
    'EpC':['Epcam','Wfdc2']
}


#----end markers------
adata_dic = [sc.read_h5ad(batch+'.h5ad') for batch in batchname]
adata_dic=[preprocess_step1(adata_dic[i],batch=batchname[i]) for i in range(3)]

'''

'''
adatalist = [preprocess(adata_dic[i],batch=batchname[i]) for i in range(3)]
def preprocess_step2(adata,batch):
    #remove the duplicated cells
    #scrub = scr.Scrublet(adata.X)
    #doublet_scores, predicted_doublets = scrub.scrub_doublets()
    #adata.obs['doublet_scores']=doublet_scores
    #adata.obs['predicted_doublet']=predicted_doublet
    sc.external.pp.scrublet(adata)
    adata = adata[adata.obs.predicted_doublet!=True,:]
    #adata = adata0
    adata.obs =adata.obs.fillna(0)
    # Saving count data
    #adata.layers["counts"] = adata.X.copy()
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    adata.write(batch+'.h5ad')
    return(adata)
'''

'''
#merge the data
#adata_dic=[preprocess_step1(adata_dic[i],batch=batchname[i]) for i in range(3)]

sc.tl.ingest(adata_dic[1],adata_dic[0], obs="leiden")
sc.tl.ingest(adata_dic[2],adata_dic[0], obs="leiden")

adata_concat = adata_dic[0].concatenate(adata_dic[1],adata_dic[2], batch_categories=batchname)

##%%time remove batch effect
sc.external.pp.bbknn(adata_concat, batch_key="batch")  # running bbknn 1.3.6
'''

'''
sc.pp.normalize_total(adata_concat,target_sum=1e4)
sc.pp.highly_variable_genes(adata_concat, n_top_genes=2000)
#basic clustering for concentrated adata
sc.tl.pca(adata_concat)
sc.tl.umap(adata_concat)
sc.tl.leiden(adata_concat,key_added='cluster')
sc.tl.dendrogram(adata_concat, groupby='cluster')
sc.pl.dotplot(adata_concat,marker_dic,groupby='cluster',dendrogram=True,save='concat',use_raw=False)
adata_concat.write('adata_concat.h5ad')
'''
#-------use scvi concact data and transfer label-------------------
import os
import tempfile

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
import torch
from scvi.model.utils import mde

from sklearn.neighbors import KNeighborsTransformer
transformer = KNeighborsTransformer(n_neighbors=30, metric='manhattan', algorithm='kd_tree')

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

os.chdir('/home/wr/LJX/sc_in/adata_h5ad/filter_adata_20240420')

adata = sc.read_h5ad('adata_concat2.h5ad')
adata = sc.read_h5ad('adata_concat1.h5ad')

adata.obs['immune']=np.array(adata.raw.X[:,adata.var.index =="Ptprc"].todense())
adata_immune = adata[adata.obs.immune >0,:]
adata_stromal = adata[adata.obs.immune==0,:]
adata = adata_immune

def Refmodel(adata,BATCH='batch'):
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=BATCH)
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    #Now we train scVI. This should take a couple of minutes on a Colab session
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()

    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata)
    adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"])
    
    '''
    sc.pl.embedding(
    adata,
    basis="X_mde",
    color=["batch", "cell_type"],
    frameon=False,
    ncols=1,save='embeding')

    sc.pl.embedding(adata, basis="X_mde", color=["leiden"], save='SVI_Test',frameon=False, ncols=1)
    lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="cell_type",
    unlabeled_category="Unknown")
    '''

#label transfer with labeled labels 
adata.var['gene_name']=adata.var.index

adata.var.index = adata.var.gene_ids

#train ref
ref = sc.read_h5ad('/home/wr/database/mouse_sc_ref/nature_2020.h5ad')
ref_all = sc.read_h5ad('/home/wr/database/mouse_sc_ref/nature_2020.h5ad')
ref_all.layers['counts'] = ref_all.raw.X
#ref_all = sc.read_h5ad('/home/wr/database/mouse_sc_ref/nature_2020.h5ad')
ref = ref_all[ref_all.obs.tissue=='bone marrow',:]
var_names = ref.var_names.intersection(adata.var_names)
ref = ref[ :,var_names]
#adata.var['gene_symbols'] = adata.var.index
#adata.var.index = adata.var.gene_ids    

#query = adata[ :,var_names]
#ref_all = ref

query = adata
ref = ref.copy()
scvi.model.SCVI.setup_anndata(ref, batch_key="batch", layer="counts")

scvi_ref = scvi.model.SCVI(
    ref,
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
)
scvi_ref.train()
SCVI_LATENT_KEY = "X_scVI"
SCVI_PREDICTIONS_KEY='scVI_prediction'
ref.obsm[SCVI_LATENT_KEY] = scvi_ref.get_latent_representation()

scvi_ref_path='nature_2020_scvi_ref_bm'

scvi_ref.save(scvi_ref_path, overwrite=True)

sc.pp.neighbors(adata, transformer=transformer,use_rep=SCVI_LATENT_KEY)
#sc.pp.neighbors(ref, use_rep=SCVI_LATENT_KEY)
#sc.tl.leiden(pancreas_ref)
sc.tl.umap(ref)


#train query
scvi.model.SCVI.prepare_query_anndata(query, scvi_ref)

scvi_query = scvi.model.SCVI.load_query_data(query,scvi_ref)

scvi_query.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})

query.obsm[SCVI_LATENT_KEY] = scvi_query.get_latent_representation()

#query.obs[SCVI_PREDICTIONS_KEY] = scvi_query.predict()

sc.pp.neighbors(query, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(query)
sc.tl.umap(query)
'''
full = anndata.concat([query, ref])
full.obsm[SCVI_LATENT_KEY] = scvi_query.get_latent_representation(
    full
)

sc.pp.neighbors(full, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(full)
sc.tl.umap(full)
'''


SCANVI_LABELS_KEY = "labels_scanvi"

ref.obs[SCANVI_LABELS_KEY] = ref.obs["cell_type"].values

scanvi_ref = scvi.model.SCANVI.from_scvi_model(
    scvi_ref,
    unlabeled_category="Unknown",
    labels_key=SCANVI_LABELS_KEY,
)

scanvi_ref.train(max_epochs=20, n_samples_per_label=100)

SCANVI_LATENT_KEY = "X_scANVI"

ref.obsm[SCANVI_LATENT_KEY] = scanvi_ref.get_latent_representation()
sc.pp.neighbors(ref, use_rep=SCANVI_LATENT_KEY)
sc.tl.leiden(ref)
sc.tl.umap(ref)

scanvi_ref_path = "scanvi_ref_nature_2020_bm"

scanvi_ref.save(scanvi_ref_path, overwrite=True)


# again a no-op in this tutorial, but good practice to use
scvi.model.SCANVI.prepare_query_anndata(query, scanvi_ref)
scvi.model.SCANVI.prepare_query_anndata(query, scanvi_ref_path)

scanvi_query = scvi.model.SCANVI.load_query_data(query, scanvi_ref)
scanvi_query = scvi.model.SCANVI.load_query_data(query, scanvi_ref_path)

scanvi_query.train(
    max_epochs=100,
    plan_kwargs={"weight_decay": 0.0},
    check_val_every_n_epoch=10,
)


SCANVI_PREDICTIONS_KEY = "predictions_scanvi"

query.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
query.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()

sc.pp.neighbors(query, transformer=transformer,use_rep=SCANVI_LATENT_KEY)
sc.pp.neighbors(query, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(query)

sc.pl.umap(query,color=['leiden',SCANVI_PREDICTIONS_KEY,'batch'],legend_loc= 'on data',save='query_scanvi_nature2020')

query.write('query.h5ad')

df = (
    query.obs.groupby(["leiden", SCANVI_PREDICTIONS_KEY])
    .size()
    .unstack(fill_value=0)
)
norm_df = df / df.sum(axis=0)

plt.figure(figsize=(8, 8))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("Predicted")
plt.ylabel("Observed")
plt.save('heatmap.pdf')

#lvae.train(max_epochs=20, n_samples_per_label=100)

adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

adata.obsm["X_mde_scanvi"] = scvi.model.utils.mde(adata.obsm["X_scANVI"])

sc.pl.embedding(
    adata, basis="X_mde_scanvi", color=["cell_type"], ncols=1, frameon=False
)




# in order to make colors matchup
adata.obs.C_scANVI = pd.Categorical(
    adata.obs.C_scANVI.values, categories=adata.obs.cell_ontology_class.cat.categories
)


SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTIONS_KEY = "C_scANVI"
scanvi_model = lvae
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
adata.obs[SCANVI_PREDICTIONS_KEY] = scanvi_model.predict(adata)

sc.pl.embedding(
    adata,
    basis=SCANVI_MDE_KEY,
    color=["cell_ontology_class", SCANVI_PREDICTION_KEY],
    ncols=1,
    frameon=False,
    palette=adata.uns["cell_ontology_class_colors"],
    save='cell_oncology_svi'
)

adata.write('adata_SVI_LT.h5ad')
#estimation
bm = Benchmarker(
    adata,
    batch_key="batch",
    label_key="cell_type",
    embedding_obsm_keys=["X_pca", "X_scVI", "X_scANVI"],
    n_jobs=-1,
)
bm.benchmark()

df = bm.get_results(min_max_scale=False)
print(df)

#---------use celltypelist -----------
import celltypist


import re
import os
import scvi
import scipy
import pickle
import anndata
import logging
import warnings
import matplotlib
import celltypist
import scipy.stats
import scanpy as sc
import scrublet as scr
import scanpy.external as sce

import pandas as pd
import numpy as np
import seaborn as sns

from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.legend import Legend
import matplotlib.gridspec as gridspec

os.chdir('~/LJX/sc_in/adata_h5ad/filter_adata_20240420')

sce.pp.harmony_integrate(adata, 'batch')