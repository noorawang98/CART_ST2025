#-*- encoding: UTF-8 -*-
#sys.setdefaultencoding('utf-8')

#import destvi_utils
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvi
from scvi.model import CondSCVI, DestVI
import os
import torch

# 设置GPU显存限制
#YOUR_LIMIT=
#torch.cuda.set_limit(torch.cuda.current_device(), limit=YOUR_LIMIT)

path='/home/wr/LJX/destVI_test/in'
os.chdir(path)

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

sc_adata = sc.read_h5ad('sc.h5ad')
sc_adata =sc.read_h5ad('/home/wr/LJX/sc_in/adata_h5ad/filter_adata_20240420/query2.h5ad')
sc_adata.layers["counts"] = sc_adata.X.copy()


st_adata = sc.read_h5ad('r_l2.h5ad')
st_adata.layers["counts"] = st_adata.X.copy()

sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)
st_adata.raw = st_adata

intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()
G = len(intersect)

sc.pp.highly_variable_genes(
    sc_adata, n_top_genes=G, subset=True, layer="counts", flavor="seurat_v3"
)

sc.pp.normalize_total(sc_adata, target_sum=10e4)
sc.pp.log1p(sc_adata)
sc_adata.raw = sc_adata

#fit scLVM
CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key="subtype")
sc_model = CondSCVI(sc_adata, weight_obs=False)
#sc_model.view_anndata_setup()
sc_model.train()
sc_model.history["elbo_train"].iloc[5:].plot()
plt.savefig('Elbo_sc_300epoch_r.pdf')

#R stLVM
#Deconvolution with stLVM

DestVI.setup_anndata(st_adata, layer="counts")

st_model = DestVI.from_rna_model(st_adata, sc_model)
#st_model.view_anndata_setup()

st_model.train(max_epochs=2500)
st_model.history["elbo_train"].iloc[10:].plot()
plt.savefig('Elbo_st_r_2500epoch.pdf')

#cell type proportions
st_adata.obsm["proportions"] = st_model.get_proportions()
st_adata.obsm["proportions"].head(5)

ct_list = st_adata.obsm['proportions'].keys()
ct = []
for ct0 in ct_list:
    ct.append(ct0)

ct_list = ct
#ct_list = ["B cells", "CD8 T cells", "Monocytes"]
for ct in ct_list:
    data = st_adata.obsm["proportions"][ct].values
    st_adata.obs[ct] = np.clip(data, 0, np.quantile(data, 0.99))

sc.pl.embedding(st_adata, basis="spatial", color=ct_list, cmap="viridis",s=80, save='sc2st_R.pdf')
ct_thresholds = destvi_utils.automatic_proportion_threshold(
    st_adata, ct_list=ct_list, kind_threshold="secondary"
)


#NR stLVM
#Deconvolution with stLVM
st_adata = sc.read_h5ad('nr_l2.h5ad')
st_adata.layers["counts"] = st_adata.X.copy()

sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)
st_adata.raw = st_adata

intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()
G = len(intersect)
sc.pp.highly_variable_genes(
    sc_adata, n_top_genes=G, subset=True, layer="counts", flavor="seurat_v3"
)

sc.pp.normalize_total(sc_adata, target_sum=10e4)
sc.pp.log1p(sc_adata)
sc_adata.raw = sc_adata

sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)
st_adata.raw = st_adata

intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()
G = len(intersect)


#fit scLVM
CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key="celltype")
sc_model = CondSCVI(sc_adata, weight_obs=False)
#sc_model.view_anndata_setup()
sc_model.train()
sc_model.history["elbo_train"].iloc[5:].plot()
plt.savefig('Elbo_sc_300epoch_nr.pdf')

DestVI.setup_anndata(st_adata, layer="counts")

st_model = DestVI.from_rna_model(st_adata, sc_model)
#st_model.view_anndata_setup()

st_model.train(max_epochs=2500)
st_model.history["elbo_train"].iloc[10:].plot()
plt.savefig('Elbo_st_nr_2500epoch.pdf')

#cell type proportions
st_adata.obsm["proportions"] = st_model.get_proportions()
st_adata.obsm["proportions"].head(5)

ct_list = st_adata.obsm['proportions'].keys()
ct = []
for ct0 in ct_list:
    ct.append(ct0)

ct_list = ct
#ct_list = ["B cells", "CD8 T cells", "Monocytes"]
for ct in ct_list:
    data = st_adata.obsm["proportions"][ct].values
    st_adata.obs[ct] = np.clip(data, 0, np.quantile(data, 0.99))


st_adata.write_h5ad('nr_l2_ano.h5ad')

#move data from cpu to gpu
#data = data.to("cpu")
#model = model.to("cpu")





