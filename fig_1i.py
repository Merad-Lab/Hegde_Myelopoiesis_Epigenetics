import pycisTopic
from pathlib import Path
import pickle
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

with open(Path(OUT_PATH, "[/path/to/object/].pkl"), "rb") as f:
  cistopic_obj = pickle.load(f)

df_plot = cistopic_obj.selected_model.cell_topic.T[["Topic5"]]

df_plot["disease"] = cistopic_obj.cell_data.loc[df_plot.index, "disease"]
df_plot["Level 2"] = cistopic_obj.cell_data.loc[df_plot.index, "Level 2"].astype("str")

celltypes = ["HSPC", "LMPP", "CD14+ monocyte"]

df_plot = df_plot[df_plot["Level 2"].isin(celltypes)].copy()

df_plot.loc[df_plot["disease"].isna(), "disease"] = "HD"

fig, ax = plt.subplots(figsize = (4, 3))

g = sb.violinplot(
    data = df_plot,
    order = celltypes,
    hue_order = ["HD", "LUAD"],
    x = "Level 2",
    y = "Topic5",
    hue = "disease",
    split = True,
    inner = None,
    palette = {"LUAD": "red", "HD": "lightgray"},
    ax = ax
)

g.set_ylim(-0.0035, 0.03)
g.set_xlabel(None)
g.set_ylabel("Score")
g.set_title("Topic 5 score, per cell")
g.legend_.set_title(None)

sb.despine()
sb.move_legend(g, bbox_to_anchor = (1, 0.5), loc = "center left")

fig.savefig(
    Path(PLOTS_PATH, "[/path/to/output/].pdf"),
    bbox_inches = "tight"
)
