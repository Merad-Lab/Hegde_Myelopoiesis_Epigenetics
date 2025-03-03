### Imports
import scanpy as sc
from pathlib import Path
import anndata as ad
import matplotlib.pyplot as plt

### Replace with path to object
ad_int = sc.read_h5ad(Path("[/path/to/your/object]"))

### Subset to celltypes of interest
keep_celltypes = [
    "HSPC",
    "LMPP",
    "CD14+ monocyte"
]

### Subsample equal number of cells per celltype of interest
n_obs = ad_int.obs["Level 2"].value_counts()[keep_celltypes].min()

ad_subset = []

for celltype in keep_celltypes:
    ad_celltype = sc.pp.subsample(
        ad_int[ad_int.obs["Level 2"] == celltype, :].copy(),
        n_obs = n_obs,
        random_state = 123,
        copy = True
    )

    ad_subset.append(ad_celltype)

ad_subset = ad.concat(ad_subset)

ad_subset.obs["Level 2"] = ad_subset.obs["Level 2"].cat.reorder_categories(
    keep_celltypes,
    ordered = True
)

### Gene list 
genes = {
    "HSPC / LMPP" : [
        "ZNF385D",
        "KIT",
        "PRSS57",
        "CRHBP",
        "MED12L",
        "CDK6",
        "MEIS1",
        "NKAIN2",
        "SPINK2",
        "CD34",
        "STMN1",
        "FAM30A",
        "ADAM28",
        "AFF3",
        "TCF4",
        "RUNX2",
        "DOCK1",
        "HOXA9",
        "COBLL1",
        "MME",
        "IGHM",
    ],
    "CD14+ monocyte" : [
        "NEAT1",
        "IRAK3",
        "CSF3R",
        "NLRP3",
        "VCAN",
        "CLEC7A",
        "CD36",
        "TLR2",
        "MYO1F",
        "TREM1",
        "LYZ",
        "CD163",
        "S100A8",
        "CD14",
    ]
}

### Plot heatmap
g = sc.pl.heatmap(
    ad_subset,
    cmap = "gist_heat_r",
    groupby = "Level 2",
    var_names = genes,
    standard_scale = "var",
    vmax = 0.5,
    show = False,
    figsize = (5, 2.5)
)

labels = g["heatmap_ax"].get_xticklabels()

for lbl in labels:
    lbl.set_style('italic')

g["groupby_ax"].set_ylabel("")

plt.savefig("[/path/to/your/output].pdf")

