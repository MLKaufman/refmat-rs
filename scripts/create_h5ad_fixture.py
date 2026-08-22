import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


counts = np.array(
    [
        [1, 0, 3],
        [0, 2, 0],
        [5, 0, 1],
        [0, 4, 0],
    ],
    dtype=np.float64,
)
obs = pd.DataFrame(
    {
        "cell_type": pd.Categorical(["T cell", "B cell", "T cell", "B cell"]),
        "batch": ["one", "one", "two", "two"],
        "score": [0.1, 0.2, 0.3, 0.4],
        "flag": [True, False, True, False],
    },
    index=[f"Cell{i}" for i in range(1, 5)],
)
var = pd.DataFrame(index=["GeneA", "GeneB", "GeneC"])
adata = ad.AnnData(X=np.log1p(counts), obs=obs, var=var)
adata.layers["counts"] = sparse.csr_matrix(counts)
adata.write_h5ad("testdata/anndata.h5ad", compression="gzip")

expected = []
for cell_type in obs["cell_type"].cat.categories:
    expected.append(np.log1p(counts[obs["cell_type"] == cell_type].mean(axis=0)))
pd.DataFrame(
    np.stack(expected, axis=1),
    index=var.index,
    columns=obs["cell_type"].cat.categories,
).rename_axis("gene").to_csv("testdata/anndata.expected.tsv", sep="\t")
