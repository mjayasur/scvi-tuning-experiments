def silhouette(adata, group_key, embed, metric="euclidean", scale=True):
    """
    wrapper for sklearn silhouette function values range from [-1, 1] with 1 being an ideal fit, 0 indicating
    overlapping clusters and -1 indicating misclassified cells
    :param group_key: key in adata.obs of cell labels
    :param embed: embedding key in adata.obsm, default: 'X_pca'
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f"{embed} not in obsm")
    asw = sklearn.metrics.silhouette_score(
        X=adata.obsm[embed], labels=adata.obs[group_key], metric=metric
    )
    if scale:
        asw = (asw + 1) / 2
    return asw


# LISI graph function (analoguous to lisi function)
def lisi_graph(
    adata,
    batch_key=None,
    label_key=None,
    k0=90,
    type_=None,
    subsample=None,
    scale=True,
    multiprocessing=None,
    nodes=None,
    verbose=False,
):
    """
    Compute lisi score (after integration)
    params:
        adata: adata object to calculate on
        batch_key: variable to compute iLISI on
        label_key: variable to compute cLISI on
        k0: number of nearest neighbors to compute lisi score
            Please note that the initial neighborhood size that is
            used to compute shortest paths is 15.
        type_: type of data integration, either knn, full or embed
        subsample: Percentage of observations (integer between 0 and 100)
                   to which lisi scoring should be subsampled
        scale: scale output values (True/False)
        multiprocessing: parallel computation of LISI scores, if None, no parallisation
                         via multiprocessing is performed
        nodes: number of nodes (i.e. CPUs to use for multiprocessing); ignored, if
               multiprocessing is set to None
    return:
        pd.DataFrame with median cLISI and median iLISI scores
        (following the harmony paper)
    """

    # recompute neighbours
    if type_ == "embed":
        adata_tmp = sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_emb", copy=True)
    if type_ == "full":
        if "X_pca" not in adata.obsm.keys():
            sc.pp.pca(adata, svd_solver="arpack")
        adata_tmp = sc.pp.neighbors(adata, n_neighbors=15, copy=True)
    else:
        adata_tmp = adata.copy()
    # if knn - do not compute a new neighbourhood graph (it exists already)

    # compute LISI score
    ilisi_score = lisi_graph_py(
        adata=adata_tmp,
        batch_key=batch_key,
        n_neighbors=k0,
        perplexity=None,
        subsample=subsample,
        multiprocessing=multiprocessing,
        nodes=nodes,
        verbose=verbose,
    )

    clisi_score = lisi_graph_py(
        adata=adata_tmp,
        batch_key=label_key,
        n_neighbors=k0,
        perplexity=None,
        subsample=subsample,
        multiprocessing=multiprocessing,
        nodes=nodes,
        verbose=verbose,
    )

    # iLISI: nbatches good, 1 bad
    ilisi_score = np.nanmedian(ilisi_score)
    # cLISI: 1 good, nbatches bad
    clisi_score = np.nanmedian(clisi_score)

    if scale:
        # get number of batches
        nbatches = len(np.unique(adata.obs[batch_key]))
        ilisi_score, clisi_score = scale_lisi(ilisi_score, clisi_score, nbatches)

    return ilisi_score, clisi_score


def sum_metric(model):
    """
    ilsi_score + silhouette_score
    """
    model.is_trained_ = True
    latent = model.get_latent_representation()
    model.is_trained_ = False
    adata = model.adata
    adata.obsm["X_scVI"] = latent
    silhouette_score = silhouette(adata, "celltype", "X_scVI")
    ilisi_score = lisi_graph(adata, batch_key="tech", type_="embed")[0]
    print("hi")
    return ilisi_score + silhouette_score