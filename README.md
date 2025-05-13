# Hegde_Myelopoiesis_Epigenetics

This repository contains data and code relevant to analysis of human and mouse single-cell datasets presented in Hegde S., et al. Nature, in revision (2025).

### Repository contents

This repository contains two main folders which accompany the uploaded sequencing and processed objects:

1. `Accessing manuscript single-cell annotations` for re-analysis of raw sequencing data using new pipelines. 
2. `Code for figure panels` provides information on how data was visualized. 

### Downloading manuscript data

Raw human and mouse sequencing data and main processed objects are available on NCBI GEO under accession number: [GSE270148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270148).

## 1. Accessing manuscript single-cell annotations 

For convenience, filtered cell barcode to annotation csvs are provided for each single-cell analyses for the following datasets in the `1. Accessing manuscript single-cell annotations` folder:

| Dataset | Filename | 
| --- | --- |
| Mouse scRNA Naive vs. Tumor, lung | Hegde_2025_Mouse_Lung_KP_scRNA_bc_cluster_annot.csv
| Mouse scRNA Naive vs. Tumor, BM | Hegde_2025_Mouse_BM_KP_scRNA_bc_cluster_annot.csv
| Mouse scATAC Naive vs. Tumor, lung | Hegde_2025_Mouse_Lung_ATAC_bc_cluster_annot.csv |
| Mouse scATAC Naive vs. Tumor, BM | Hegde_2025_Mouse_BM_ATAC_bc_cluster_annot.csv |
| Mouse scRNA NFE2L2 KO vs. WT, lung | Hegde_2025_Mouse_Lung_Nfe2l2_scRNA_bc_cluster_annot.csv
| Mouse scRNA NFE2L2 KO vs. WT, BM | Hegde_2025_Mouse_BM_Nfe2l2_scRNA_bc_cluster_annot.csv
| Mouse scATAC NFE2L2 KO vs. WT, BM | Hegde_2025_Mouse_BM_Nfe2l2_scATAC_bc_cluster_annot.csv
| Human scMultiome (ATAC), LUAD lung tissue | Hegde_2025_Human_Multiome_ATAC_bc_cluster_annot.csv |
| Human scMultiome (RNA), LUAD lung tissue | Hegde_2025_Human_Multiome_RNA_bc_cluster_annot.csv |
| Human scMultiome (RNA), CD34+ enriched LUAD PBMCs | Hegde_2025_HumanCD34enrich_Multiome_RNA_bc_cluster_annot.csv |

## 2. Code for figure panels

This repository contains code for reproducing figure panels based on single-cell RNA, ATAC, and bulk CUT&RUN data analyses. To help navigate the wide variety of data and code, please refer to the following tables:

| Folder | Manuscript Panels | Assay | Sample | Condition | Contents | 
| --- | --- | --- | --- | --- | --- |
| `Figure 1/fig_1c_edf_1f` | Fig 1c, EDF 1f | scRNA | mouse BM | tumor and naive | pre-processing, annotation, UCell scoring |
| `Figure 1/fig_1e_edf_1ij` | Fig 1e, EDF 1i-j | scATAC | mouse BM | tumor and naive | ArchR pre-processing, gene scoring, motif enrichment |
| `Figure 1/fig_1fg_edf_1k` | Fig 1f-g, EDF 1k | bulk CUT&RUN | mouse BM | preprocessing, tumor and naive | plotting profiles, heatmaps, GREAT |
| `Figure 1/fig_1i` | Fig 1i | scMultiome, RNA | human CD34-enriched PBMCs | lung cancer and healthy | preprocessing, UCell scoring |
| `Figure 1/fig_1j_edf_1p` | Fig 1j, EDF 1p | scMultiome, ATAC | human CD34-enriched PBMCs | lung cancer and healthy | preprocessing, pycisTopic scoring, GREAT |
| `Figure 2/fig_2bce_edf_4efgh` | Fig 2b-c,e, EDF 4e-h | scMultiome | human tumor and PBMCs | NSCLC patients | pre-processing, cell marker plot, motif enrichment, differential analysis, WCGNA on myeloid cells |
| `Figure 2/fig_2d_edf_4ab` | Fig 2d, EDF 4ab | scRNA | mouse lung | tumor-bearing | pre-processing, cell marker plot, differential analysis, module analysis |
| `Figure 4/fig_4bef_edf_8ijklmnopq` | Fig 4b,e-f, EDF 8i-q | scRNA + scATAC | mouse BM and lung | tumor-bearing, Nfe2l2 KO and wild-type | pre-processing, annotation, differential analysis, scoring |
| `EDF 1/edf_1lmn` | EDF 1l-n | bulk CUT&RUN | mouse BM | tumor and naive | plotting profiles, heatmaps, and chromVAR scores |
| `EDF 1/edf_1o` | EDF 1o | scMultiome | human CD34-enriched PBMCs | lung cancer and healthy | cell marker plot |
| `EDF 3/edf_3i` | EDF 3i | scRNA | public data, human NSCLC | pre/post treatment | pre-processing, annotation, and plotting of NRF2 scores and ChIP-X enrichment analysis |
| `EDF 4/edf_4cd` | EDF 4c-d | scATAC | mouse lung | tumor and naive | pre-processing, annotation, cell marker plot, motif enrichment |
| `EDF 10` | EDF 10a-d | scRNA | public data, human solid tumors | pre/post treatment | pre-processing, annotation, and plotting of ChIP-X enrichment analysis |

## Disclaimer

Reproducing these single-cell analyses with different/updated software versions or random generator seeds may result in slight variations. Please see the Methods section for details, and contact the corresponding authors if you have any questions regarding specific packages or data-preprocessing steps. 