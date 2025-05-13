### Imports

import os
import pickle
import pypumatac as pum  # From the PUMATAC pipeline
from pycisTopic.qc import *
import ray
import pprint as pp
from glob import glob
import numpy as np

import pyranges as pr
import requests

import scrublet as scr

import anndata as ad
from pathlib import Path
import pickle

from pycisTopic.pseudobulk_peak_calling import *
from pycisTopic.iterative_peak_calling import *

from pycisTopic.cistopic_class import run_cgs_models

from pycisTopic.topic_binarization import binarize_topics

### The beginning portion of this script utilizes the PUMATAC pipeline to filter down to a list of good barcodes before counts matrix generation (see https://www.nature.com/articles/s41587-023-01881-x and https://github.com/aertslab/PUMATAC). All credits belong to the authors of the PUMATAC pipeline, and please see their Github repository for more information on installation and usage.

CELLRANGER_PATH = "/path/to/cellranger/inputs"
PUMATAC_PATH = "/path/to/pumatac/outputs"
PUMATAC_DEPENDENCIES_PATH = "/path/to/pumatac/dependencies"

fragments_path_dict = {
    "Lu927": Path(CELLRANGER_PATH, "hs_multiome_PBMC_Lu927_fragments.tsv.gz"),
    "Lu941": Path(CELLRANGER_PATH, "hs_multiome_PBMC_Lu941_fragments.tsv.gz"),
    "Lu954": Path(CELLRANGER_PATH, "hs_multiome_PBMC_Lu954_fragments.tsv.gz"),
    "Lu979": Path(CELLRANGER_PATH, "hs_multiome_PBMC_Lu979_fragments.tsv.gz"),
    "Lu1024": Path(CELLRANGER_PATH, "hs_multiome_PBMC_Lu1024_fragments.tsv.gz"),
    "Lu1027": Path(CELLRANGER_PATH, "hs_multiome_PBMC_Lu1027_fragments.tsv.gz"),
    "Lu1032": Path(CELLRANGER_PATH, "hs_multiome_PBMC_Lu1032_fragments.tsv.gz"),
    "HD4": Path(CELLRANGER_PATH, "hs_multiome_PBMC_HD4_fragments.tsv.gz"),
    "HD5": Path(CELLRANGER_PATH, "hs_multiome_PBMC_HD5_fragments.tsv.gz"),
}

os.chdir(PUMATAC_PATH)

output_dir = "qc_outputs"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

cistopic_qc_out = "cistopic_qc_outputs"
if not os.path.exists(cistopic_qc_out):
    os.makedirs(cistopic_qc_out)

selected_barcodes = "selected_barcodes"
if not os.path.exists(selected_barcodes):
    os.mkdir(selected_barcodes)

qc_plots = "qc_plots"
if not os.path.exists(qc_plots):
    os.mkdir(qc_plots)

saturation_stats_path = "saturation_stats"
if not os.path.exists(saturation_stats_path):
    os.mkdir(saturation_stats_path)


pipeline_dict = {
    sample_id: "cellranger_arc" for sample_id in fragments_path_dict.keys()
}

genome_dict = {sample_id: "hg38" for sample_id in fragments_path_dict.keys()}

alias_dict = {sample_id: sample_id for sample_id in fragments_path_dict.keys()}

inverse_genome_dict = {}
for sample, genome in genome_dict.items():
    if genome not in inverse_genome_dict:
        inverse_genome_dict[genome] = []
    inverse_genome_dict[genome].append(sample)

annotation_dict = pum.download_genome_annotation(inverse_genome_dict)

regions_path_dict = {
    "hg38": Path(
        PUMATAC_DEPENDENCIES_PATH, "V2.hg38-rDHS-Unfiltered.blacklisted.bed"
    ),
}

samples = fragments_path_dict.keys()

samples_todo = []
for sample in samples:
    metadata_file = os.path.join(cistopic_qc_out, sample + "_metadata_bc.pkl")
    print(f"Checking if {metadata_file} exist...")
    if os.path.exists(metadata_file):
        print("\tMetadata exists! Skipping...")
    else:
        samples_todo.append(sample)
        print("\tMetadata does not exist, adding to subdict to generate")

samples_sizes = [
    (sample, os.path.getsize(fragments_path_dict[sample]))
    for sample in samples_todo
]

samples_sizes_sorted = sorted(samples_sizes, key=lambda x: x[1], reverse=True)
samples_todo = [sample for sample, _ in samples_sizes_sorted]

n_cores = int(
    os.environ.get("LSB_DJOB_NUMPROC", 4)
)  # This is for an LSF submission service, adjust as needed.

ray.shutdown()

for genome, inverse_genome_dict_samples in inverse_genome_dict.items():
    if samples_todo != []:
        samples_sub = list(
            set(samples_todo).intersection(inverse_genome_dict_samples)
        )
        blocks = [
            samples_sub[i : i + n_cores]
            for i in range(0, len(samples_sub), n_cores)
        ]
        pp.pprint(blocks)
        for samples_torun_in_block in blocks:
            fragments_sub_dict_block = {
                key: fragments_path_dict[key] for key in samples_torun_in_block
            }
            regions_sub_dict_block = {
                key: regions_path_dict[genome_dict[key]]
                for key in samples_torun_in_block
            }
            print(
                f"Running samples {samples_torun_in_block} for genome {genome}"
            )

            metadata_bc_dict, profile_data_dict = compute_qc_stats(
                fragments_dict=fragments_sub_dict_block,
                tss_annotation=annotation_dict[genome],
                stats=[
                    "barcode_rank_plot",
                    "duplicate_rate",
                    "insert_size_distribution",
                    "profile_tss",
                    "frip",
                ],
                label_list=None,
                path_to_regions=regions_sub_dict_block,
                n_cpu=n_cores,
                valid_bc=None,
                n_frag=10,
                n_bc=None,
                tss_flank_window=2000,
                tss_window=50,
                tss_minimum_signal_window=100,
                tss_rolling_window=10,
                # min_norm = 0.2,
                remove_duplicates=True,
            )

            ray.shutdown()
            print(f"Done, writing files to {cistopic_qc_out}...")
            for sample in sorted(metadata_bc_dict.keys()):
                metadata_bc_dict[sample]["sample_id"] = sample
                metadata_bc_dict[sample].index = [
                    x + "_" + sample
                    for x in list(metadata_bc_dict[sample].index)
                ]
                with open(
                    os.path.join(cistopic_qc_out, f"{sample}_metadata_bc.pkl"),
                    "wb",
                ) as f:
                    pickle.dump(metadata_bc_dict[sample], f, protocol=4)

                with open(
                    os.path.join(cistopic_qc_out, f"{sample}_profile_data.pkl"),
                    "wb",
                ) as f:
                    pickle.dump(profile_data_dict[sample], f, protocol=4)
    else:
        print(f"All samples already processed  for genome {genome}.")


metadata_bc_pkl_path_dict = {
    os.path.basename(x).split("_metadata_bc.pkl")[0]: x
    for x in sorted(glob(f"{cistopic_qc_out}/*metadata_bc.pkl"))
}

standard_min_x_val = 100
standard_min_y_val = 1

min_otsu_frags_dict = {}
min_otsu_tss_dict = {}

for sample in metadata_bc_pkl_path_dict.keys():
    min_otsu_frags_dict[sample] = standard_min_x_val
    min_otsu_tss_dict[sample] = standard_min_y_val

x_threshold_dict = {}
y_threshold_dict = {}
bc_dict = {}
n_bc_dict = {}

for sample in metadata_bc_pkl_path_dict.keys():
    with open(metadata_bc_pkl_path_dict[sample], "rb") as fh:
        metadata_bc_df = pickle.load(fh)

    x_arr = np.log10(metadata_bc_df["Unique_nr_frag_in_regions"])
    x_threshold_log = pum.threshold_otsu(
        x_arr, nbins=5000, min_value=np.log10(min_otsu_frags_dict[sample])
    )
    x_threshold = 10**x_threshold_log
    x_threshold_dict[sample] = x_threshold

    y_arr = metadata_bc_df["TSS_enrichment"]
    y_threshold = pum.threshold_otsu(
        y_arr, nbins=5000, min_value=min_otsu_tss_dict[sample]
    )
    y_threshold_dict[sample] = y_threshold

    metadata_bc_df_passing_filters = metadata_bc_df.loc[
        (metadata_bc_df.Unique_nr_frag_in_regions > x_threshold)
        & (metadata_bc_df.TSS_enrichment > y_threshold)
    ]
    bc_passing_filters = metadata_bc_df_passing_filters.index
    bc_dict[sample] = bc_passing_filters
    n_bc_dict[sample] = len(bc_passing_filters)

    print(f"\tSaving...")
    with open(
        f"selected_barcodes/{sample}_bc_passing_filters_otsu.pkl", "wb"
    ) as fh:
        pickle.dump(bc_passing_filters, fh)
    fh.close()

    fh = open(f"selected_barcodes/{sample}_bc_passing_filters_otsu.txt", "w")
    for bc in list(bc_passing_filters):
        fh.write(bc + "\n")
    fh.close()

    metadata_bc_df.loc[bc_passing_filters].to_csv(
        f"selected_barcodes/{sample}_metadata_bc_df.tsv", sep="\t"
    )

### Begin cisTopic quality control and topic modeling. This section of the code is adapted from the pycisTopic pipeline (https://github.com/aertslab/pycisTopic) and (https://www.nature.com/articles/s41592-023-01938-4). All credits belong to the authors of the pycisTopic pipeline, and please see their Github repository for more information on installation and usage.

CISTOPIC_PATH = "/path/to/cistopic/outputs"
TMP_DIR = "/path/to/scratch/dir/with/lots/of/space"
RNA_OBJ = "/path/to/rna/obj/with/metadata"

# Output directory
OUT_PATH = Path(CISTOPIC_PATH, "output/")
OUT_PATH.mkdir(exist_ok=True)

PEAK_PATH = Path(OUT_PATH, "consensus_peak_calling/")
PEAK_PATH.mkdir(exist_ok=True)

BED_PATH = Path(PEAK_PATH, "pseudobulk_bed_files/")
BED_PATH.mkdir(exist_ok=True)

BW_PATH = Path(PEAK_PATH, "pseudobulk_bw_files")
BW_PATH.mkdir(exist_ok=True)

MACS2_PATH = Path(PEAK_PATH, "MACS2")
MACS2_PATH.mkdir(exist_ok=True)

MODELS_PATH = Path(OUT_PATH, "models")
MODELS_PATH.mkdir(exist_ok=True)

PLOTS_PATH = Path(OUT_PATH, "plots")
PLOTS_PATH.mkdir(exist_ok=True)

DAR_PATH = Path(OUT_PATH, "DAR")
DAR_PATH.mkdir(exist_ok=True)

ad_rna = ad.read_h5ad(RNA_OBJ)

df_rna_metadata = ad_rna.obs.copy()

df_rna_metadata["barcode"] = df_rna_metadata.index.str.split(
    "_", expand=True
).get_level_values(0)

df_rna_metadata["Level 2"] = df_rna_metadata[
    "Level 2"
].astype("str")

### Get chromosome sizes

target_url = (
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
)

chromsizes = pd.read_csv(target_url, sep="\t", header=None)
chromsizes.columns = ["Chromosome", "End"]
chromsizes["Start"] = [0] * chromsizes.shape[0]
chromsizes = chromsizes.loc[:, ["Chromosome", "Start", "End"]]

chromsizes["Chromosome"] = [
    chromsizes["Chromosome"][x].replace("v", ".")
    for x in range(len(chromsizes["Chromosome"]))
]
chromsizes["Chromosome"] = [
    (
        chromsizes["Chromosome"][x].split("_")[1]
        if len(chromsizes["Chromosome"][x].split("_")) > 1
        else chromsizes["Chromosome"][x]
    )
    for x in range(len(chromsizes["Chromosome"]))
]

chromsizes = pr.PyRanges(chromsizes)

### Export pseudobulks based on RNA annotations

ray.shutdown()

bw_paths, bed_paths = export_pseudobulk(
    input_data=df_rna_metadata,
    variable="broad_annotations",
    sample_id_col="sample_id",
    chromsizes=chromsizes,
    bed_path=str(BED_PATH),
    bigwig_path=str(BW_PATH),
    path_to_fragments=fragments_path_dict,
    n_cpu=8,
    normalize_bigwig=True,
    remove_duplicates=True,
    _temp_dir=str(Path(TMP_DIR, "ray_spill")),
    split_pattern="_",
)

with open(Path(BED_PATH, "bed_paths.pkl"), "wb") as f:
    pickle.dump(bed_paths, f)

with open(Path(BW_PATH, "bw_paths.pkl"), "wb") as f:
    pickle.dump(bw_paths, f)


# Call peaks

with open(Path(BED_PATH, "bed_paths.pkl"), "rb") as f:
    bed_paths = pickle.load(f)

bed_paths_quote = {k: '"' + str(v) + '"' for k, v in bed_paths.items()}

ray.shutdown()

MACS2_BIN_PATH = "/path/to/macs2"

narrow_peaks_dict = peak_calling(
    MACS2_BIN_PATH,
    bed_paths_quote,
    str(MACS2_PATH),
    genome_size="hs",
    n_cpu=12,
    input_format="BEDPE",
    shift=73,
    ext_size=146,
    keep_dup="all",
    q_value=0.05,
    _temp_dir=str(Path(TMP_DIR, "ray_spill")),
)

with open(Path(MACS2_PATH, "narrow_peaks_dict.pkl"), "wb") as f:
    pickle.dump(narrow_peaks_dict, f)


### Iterative peak reduction to get consensus peaks

with open(Path(MACS2_PATH, "narrow_peaks_dict.pkl"), "rb") as f:
    narrow_peaks_dict = pickle.load(f)

peak_half_width = 250
path_to_blacklist = "/path/to/cistopic/dependencies/hg38-blacklist.v2.bed"

consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict,
    peak_half_width,
    chromsizes=chromsizes,
    path_to_blacklist=path_to_blacklist,
)

consensus_peaks.to_bed(
    path=Path(PEAK_PATH, "consensus_regions.bed"),
    keep=True,
    compression="infer",
    chain=False,
)

metadata_bc = {}
bc_passing_filters = {}

for sample in fragments_path_dict.keys():
    with open(
        f"{PUMATAC_PATH}/cistopic_qc_outputs/{sample}_metadata_bc.pkl",
        "rb",
    ) as f:
        metadata_bc_sample = pickle.load(f)

        metadata_bc_sample["full_barcode"] = metadata_bc_sample.index.values

        metadata_bc_sample.index = metadata_bc_sample.index.str.split(
            "_", expand=True
        ).get_level_values(0)

        metadata_bc[sample] = metadata_bc_sample

    with open(
        f"{PUMATAC_PATH}/selected_barcodes/{sample}_bc_passing_filters_otsu.pkl",
        "rb",
    ) as f:
        bc_passing_filters_sample = pickle.load(f)

        bc_passing_filters_sample = bc_passing_filters_sample.str.split(
            "_", expand=True
        ).get_level_values(0)

        bc_passing_filters[sample] = bc_passing_filters_sample.to_list()


ray.shutdown()

path_to_regions = Path(PEAK_PATH, "consensus_regions.bed")

cistopic_obj_list = [
    create_cistopic_object_from_fragments(
        path_to_fragments=fragments_path_dict[key],
        path_to_regions=str(path_to_regions),
        path_to_blacklist=str(path_to_blacklist),
        metrics=metadata_bc[key],
        valid_bc=bc_passing_filters[key],
        n_cpu=6,
        project=key,
    )
    for key in fragments_path_dict.keys()
]

cistopic_obj = merge(cistopic_obj_list)

scrub = scr.Scrublet(cistopic_obj.fragment_matrix.T, expected_doublet_rate=0.1)

doublet_scores, predicted_doublets = scrub.scrub_doublets()

scrublet = pd.DataFrame(
    [scrub.doublet_scores_obs_, scrub.predicted_doublets_],
    columns=cistopic_obj.cell_names,
    index=["Doublet_scores_fragments", "Predicted_doublets_fragments"],
).T

cistopic_obj.add_cell_data(scrublet, split_pattern="-")
cistopic_obj.cell_data["Predicted_doublets_fragments"] = cistopic_obj.cell_data[
    "Predicted_doublets_fragments"
].astype("str")

ray.shutdown()

models = run_cgs_models(
    cistopic_obj,
    n_topics=[70], # We anecdotally have found that using higher numbers of topics (e.g. 70) works better in our hands, without need for further tuning (which can be computationally expensive)
    n_cpu=6,
    n_iter=100,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    save_path=None,
    _temp_dir=str(TMP_DIR),
)

### Because this is Multiome data, we can filter based on RNA barcodes as well. 

ad_rna.obs["barcode"] = ad_rna.obs_names.str.split(
    "_", expand=True
).get_level_values(0)

ad_rna.obs["full_barcode"] = (
    ad_rna.obs["barcode"].astype("str")
    + "___"
    + ad_rna.obs["sample_id"].astype("str")
)

ad_rna.obs["full_barcode"] = ad_rna.obs["full_barcode"].astype("str")

intersect_cells = set(ad_rna.obs["full_barcode"]).intersection(
    cistopic_obj.cell_names
)

cistopic_obj.subset(cells=intersect_cells)

cistopic_obj.cell_data["Level 2"] = (
    ad_rna.obs[["full_barcode", "Level 2"]]
    .set_index("full_barcode")
    .loc[cistopic_obj.cell_names, "Level 2"]
)

binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target="cell",
    method="ntop",
    ntop=1_000,
    plot=True,
    num_columns=5,
    nbins=100,
)
