### Imports
import pandas as pd
import seaborn as sb
import numpy as np

### Replace with path to ImmGen CUT&RUN DORC matrix
dorc_matrix = pd.read_csv(
    "[/path/to/your/object]",
    sep = "\t"
)

### Genes
genes = [
    "S100a9",
    "S100a8",
    "Plaur",
    "Il6ra",
    "Adrb2",
    "Hif1a",
    "Lcn2",
    "Chil3",
    "Plaur",
    "Vsir",
    "Bst1",
    "Cd177",
    "Fn1",
    "Cd14",
    "Irf8",
    "Ciita",
    "H2-Ab1",
    "H2-Eb1",
    "H2-Aa"
    "H2-Ob",
    "Ms4a4c",
    "Ifitm3",
    "Oasl1",
    "Mx1",
    "Usp18",
    "Gbp4",
]

### Subset and plot heatmap
samples = dorc_matrix.columns

celltypes = ["GMP"]
marks = ["H3K4me1", "H3K27ac", "H3K27me3"]

data = {}

kp_mask = samples.str.contains("KP")

for c in celltypes:
    celltype_mask = samples.str.contains(c)

    for m in marks:
        mark_mask = samples.str.contains(m)

        naive_data = (
            dorc_matrix.loc[
                :,
                celltype_mask & mark_mask & ~kp_mask
            ].mean(axis = 1)
        )

        kp_data = (
            dorc_matrix.loc[
                :,
                celltype_mask & mark_mask & kp_mask
            ].mean(axis = 1)
        )

        data[f"{c}__{m}"] = np.log2((kp_data + 0.001) / (naive_data + 0.001))

df_plot = pd.DataFrame.from_dict(data)

df_plot = df_plot.loc[flat_genes, :].copy()

g = sb.clustermap(
    data = df_plot,
    cmap = "vlag",
    vmax = 1.5,
    vmin = -1.5,    
    center = 0,
    col_cluster = False,
    row_cluster = False,
    figsize = (1.75, 12),
    cbar_pos = (-0.3, 0.8, 0.05, 0.18),
    linewidth = 0.1,
    linecolor = "black",
    yticklabels = True,
    dendrogram_ratio = 0.01
)

labels = g.ax_heatmap.get_yticklabels()

for lbl in labels:
    lbl.set_style("italic")

g.ax_heatmap.tick_params(right = False, bottom = False)

g.savefig("[path/to/your/output].pdf")