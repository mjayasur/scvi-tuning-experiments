import scvi
import scanpy
from metrics import silhouette
from ray.tune.schedulers import ASHAScheduler
from ray.tune import loguniform
from scvi.autotune import tune_scvi


adata = scvi.data.read_h5ad("pancreas.h5ad")
scanpy.pp.highly_variable_genes(adata, batch_key="tech", n_top_genes=2000, subset=True)
scvi.data.setup_anndata(adata, batch_key="tech", layer="counts")


tuned_model, analysis = tune_scvi(adata, 400)

adata.obsm["X_scVI"] = tuned_model.get_latent_representation()

bdata = scvi.data.read_h5ad("pancreas.h5ad")
scanpy.pp.highly_variable_genes(bdata, batch_key="tech", n_top_genes=2000, subset=True)
scvi.data.setup_anndata(bdata, batch_key="tech", layer="counts")
default_model = scvi.model.SCVI(bdata)

default_model.train()

bdata.obsm["X_scVI"] = default_model.get_latent_representation()

results = {
    "tuned": silhouette(adata, "celltype", "X_scvi"),
    "default": silhouette(bdata, "celltype", "X_scvi"),
}

out_file = open("results.json", "w")

json.dump(results, out_file, indent=6)

json.dump(analysis.get_best_config(), open("config.json", "w"), indent=6)
out_file.close()
