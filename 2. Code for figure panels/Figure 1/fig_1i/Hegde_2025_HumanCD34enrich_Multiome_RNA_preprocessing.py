### Imports

from glob import glob
import scanpy as sc
import pandas as pd
from pathlib import Path
import os

### Defining file paths and sample locations. As a note, the example below demonstrates how the CD34-enriched PBMC data was integrated with the non-CD34-enriched Multiome PBMC samples, to help with the number of cells and annotations. However, where indicated in the manuscript, only the CD34-enriched PBMC samples were used in single scoring analysis.

CELLRANGER_PATH = "/path/to/cellranger/inputs"
CELLBENDER_PATH = "/path/to/cellbender/intermediate/files"
PUMATAC_PATH = "/path/to/pumatac/outputs"

samples = {
    "Lu927": {
        "file_path": Path(
            CELLRANGER_PATH,
            "hs_multiome_PBMC_Lu927_0_v1/raw_feature_bc_matrix.h5",
        ),
        "learning_rate": 0.0001,
    },
    "Lu941": {
        "file_path": Path(
            CELLRANGER_PATH,
            "hs_multiome_PBMC_Lu941_0_v1/raw_feature_bc_matrix.h5",
        ),
        "learning_rate": 0.0001,
    },
    "Lu954": {
        "file_path": Path(
            CELLRANGER_PATH,
            "hs_multiome_PBMC_Lu954_0_v1/raw_feature_bc_matrix.h5",
        ),
        "learning_rate": 0.0001,
    },
    "Lu979": {
        "file_path": Path(
            CELLRANGER_PATH,
            "hs_multiome_PBMC_Lu979_0_v1/raw_feature_bc_matrix.h5",
        ),
        "learning_rate": 0.0001,
    },
    "Lu1024": {
        "file_path": Path(
            CELLRANGER_PATH,
            "hs_multiome_PBMC_Lu1024_0_v1/raw_feature_bc_matrix.h5",
        ),
        "learning_rate": 0.0001,
    },
    "Lu1027": {
        "file_path": Path(
            CELLRANGER_PATH,
            "hs_multiome_PBMC_Lu1027_0_v1/raw_feature_bc_matrix.h5",
        ),
        "learning_rate": 0.0001,
    },
    "Lu1032": {
        "file_path": Path(
            CELLRANGER_PATH,
            "hs_multiome_PBMC_Lu1032_0_v1/raw_feature_bc_matrix.h5",
        ),
        "learning_rate": 0.0001,
    },
    "HD4": {
        "file_path": Path(
            CELLRANGER_PATH,
            "hs_multiome_PBMC_HD4_0_v1/raw_feature_bc_matrix.h5",
        ),
        "learning_rate": 0.0001,
    },
    "HD5": {
        "file_path": Path(
            CELLRANGER_PATH,
            "hs_multiome_PBMC_HD5_0_v1/raw_feature_bc_matrix.h5",
        ),
        "learning_rate": 0.0001,
    },
}

### Cellbender was used to remove ambient RNA from the samples. Institution-specific LSF parameters have been removed, from this script.

submission_script = """
bsub << EOF

#!/bin/bash
cd {sample_folder}

cellbender remove-background \
--cuda \
--input {file_path} \
--output {sample_id}_CellBender_remove-background.h5 \
--learning-rate {learning_rate}

EOF
"""

for sample_id, sample_dict in samples.items():
    sample_folder = Path(CELLBENDER_PATH, sample_id)

    file_path = sample_dict["file_path"]
    learning_rate = sample_dict["learning_rate"]

    sample_folder.mkdir(exist_ok=True)

    os.system(
        submission_script.format(
            sample_folder=str(sample_folder),
            file_path=file_path,
            sample_id=sample_id,
            learning_rate=learning_rate,
        )
    )

files = {
    sample: os.path.join(
        CELLBENDER_PATH,
        sample,
        f"{sample}_CellBender_remove-background_filtered.h5",
    )
    for sample in samples
}

### Since these are Multiome samples, we can additionally filter based on ATAC quality control outputs, from PUMATAC.

pumatac_barcodes_txt = glob(
    os.path.join(
        PUMATAC_PATH,
        "*_bc_passing_filters_otsu.txt",
    )
)

df_pumatac_barcodes = pd.concat(
    [pd.read_csv(file, header=None) for file in pumatac_barcodes_txt]
).reset_index(drop=True)

pumatac_barcodes_txt = [
    file
    for file in pumatac_barcodes_txt
    if any([sample in file for sample in samples])
]

### Beginning import of Cellbender outputs and filtering based on PUMATAC outputs.

ad_int = []

for sample_id, file in files.items():
    print(sample_id)
    print(file)

    ad_sample = sc.read_10x_h5(file)

    ad_sample.var_names_make_unique()

    ad_sample.obs_names = ad_sample.obs_names.astype("str") + "_" + sample_id

    ad_sample.obs["sample_id"] = sample_id

    ad_int.append(ad_sample)

ad_int = sc.concat(ad_int)

ad_int.raw = ad_int

ad_int.layers["counts"] = ad_int.X

ad_int.obs["passed_pumatac_filter"] = ad_int.obs_names.isin(
    df_pumatac_barcodes[0].values
)

ad_int = ad_int[ad_int.obs["passed_pumatac_filter"] == True, :].copy()

### Beginning standard RNA-based quality control.

ad_int.var["mt"] = ad_int.var_names.str.startswith("MT-")
ad_int.var["ribo"] = ad_int.var_names.str.startswith(("RPS", "RPL"))
ad_int.var["hb"] = ad_int.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    ad_int,
    qc_vars=["mt", "ribo", "hb"],
    percent_top=None,
    log1p=None,
    inplace=True,
)

ad_int = ad_int[ad_int.obs["pct_counts_mt"] < 40, :]
ad_int = ad_int[ad_int.obs["pct_counts_ribo"] < 30, :]
ad_int = ad_int[ad_int.obs["pct_counts_hb"] < 1, :]

sc.pp.filter_cells(ad_int, max_genes=4500)
sc.pp.filter_genes(ad_int, min_cells=3)

sc.pp.scrublet(ad_int, batch_key="sample_id", verbose=True)

ad_int = ad_int[ad_int.obs["predicted_doublet"] == False, :].copy()

### Beginning standard RNA-based normalization and integration using Harmony.

sc.pp.normalize_total(ad_int)
sc.pp.log1p(ad_int)

sc.pp.highly_variable_genes(ad_int, batch_key="sample_id")

sc.tl.pca(ad_int)

sc.external.pp.harmony_integrate(ad_int, "sample_id")

sc.pp.neighbors(ad_int, use_rep="X_pca_harmony", key_added="pca_harmony")

sc.tl.umap(ad_int, neighbors_key="pca_harmony")

### Beginning clustering and combining with manual annotations. 

sc.tl.leiden(
    ad_int,
    key_added=f"leiden_res_2.00",
    resolution=2.00,
    flavor="igraph",
    neighbors_key="pca_harmony",
)

annotations = pd.read_csv("/path/to/annotations", sep="\t")

annotations["Cluster Number"] = annotations["Cluster Number"].astype("str")

ad_int.obs = (
    ad_int.obs.reset_index()
    .merge(annotations, left_on="leiden_res_2.00", right_on="Cluster Number")
    .set_index("index")
    .copy()
)