from scipy.stats import zscore
import pandas as pd
import matplotlib.pyplot as plt

chromVAR_TF_H3K4me1 = pd.read_csv("[/path/to/H3K4me1/output/]", index_col = 6)

chromVAR_TF_H3K4me3 = pd.read_csv("[/path/to/H3K4me3/output/]", index_col = 6)

chromVAR_TF_H3K27ac = pd.read_csv("[/path/to/H3K27ac/output/]", index_col = 6)

df_plot = pd.concat(
    [
        chromVAR_TF_H3K4me1["mean_diff"].rename("H3K4me1_mean_diff"),
        chromVAR_TF_H3K4me3["mean_diff"].rename("H3K4me3_mean_diff"),
        chromVAR_TF_H3K27ac["mean_diff"].rename("H3K27ac_mean_diff"),
    
    ],
    axis = 1
)

df_plot = df_plot.dropna()

df_plot["H3K4me3_mean_diff"] = df_plot[["H3K4me3_mean_diff"]].apply(zscore)
df_plot["H3K4me1_mean_diff"] = df_plot[["H3K4me1_mean_diff"]].apply(zscore)
df_plot["H3K27ac_mean_diff"] = df_plot[["H3K27ac_mean_diff"]].apply(zscore)

df_plot = df_plot[
    ~df_plot.index.str.contains("Zfp|Zscan|Zkscan|ENSMUS|Zbtb|28047")
].copy()

x_axis = "H3K4me3_mean_diff"
y_axis = "H3K27ac_mean_diff"

fig, ax = plt.subplots(figsize = (3, 3))

g = sb.scatterplot(
    data = df_plot,
    x = x_axis,
    y = y_axis,
    alpha = 0.75,
    s = 7,
    linewidth = 0,
    ax = ax
)

g.set_xlabel("H3K4me3 (KP - Naive)")
g.set_ylabel("H3K4me1 (KP - Naive)")

sb.despine()

fig.savefig("[/path/to/your/output].svg")

x_axis = "H3K4me3_mean_diff"
y_axis = "H3K27ac_mean_diff"

fig, ax = plt.subplots(figsize = (3, 3))

g = sb.scatterplot(
    data = df_plot,
    x = x_axis,
    y = y_axis,
    alpha = 0.75,
    s = 7,
    linewidth = 0,
    ax = ax
)

g.set_xlabel("H3K4me3 (KP - Naive)")
g.set_ylabel("H3K27ac (KP - Naive)")

sb.despine()

fig.savefig("[/path/to/your/output].svg")