# Hegde_Myelopoiesis_Epigenetics

This repository contains data and code relevant to analysis of human and mouse single-cell datasets presented in Hegde S., et al. Nature, in revision (2025).

# Downloading the data

Raw human and mouse sequencing data and main processed objects are available on NCBI GEO under accession number: [GSE270148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270148).

# Contents

This repository contains code for reproducing figure panels based on bulk + single-cell RNA, ATAC, and CUT&RUN data analyses. To help align the wide variety of data and code with the manuscript, please refer to the following tables:

## Main Figures
| Panel | Assay | Sample | Condition | Contents | 
| --- | --- | --- | --- | --- |
| 1c | scRNA | mouse BM | tumor and naive | pre-processing, annotation, UCell scoring |
| 1e-f | bulk CUT&RUN | mouse BM | tumor and naive | plotting profiles, heatmaps |
| 1h | scMultiome | human CD34-enriched PBMCs | lung cancer and healthy | UCell scoring |
| 1i | scMultiome |  human CD34-enriched PBMCs | lung cancer and healthy | pycisTopic scoring |
| 2b | scMultiome | human tumor and PBMCs | NSCLC patients | pre-processing, cell marker plot, motif enrichment, differential analysis, WCGNA on myeloid cells |
| 4 | scRNA + scATAC | mouse BM and lung | tumor-bearing, Nfe2l2 KO and wild-type | pre-processing, annotation, differential analysis, scoring |

## Extended Data Figures
| Panel | Assay | Sample | Condition | Contents | 
| --- | --- | --- | --- | --- |
| 1f | scRNA | mouse BM | tumor and naive | pre-processing, annotation, UCell scoring |
| 1i-j | scATAC | mouse BM | tumor and naive | ArchR pre-processing, gene scoring, motif enrichment |
| 1k-m | bulk CUT&RUN | human CD34-enriched PBMCs | lung cancer and healthy | plotting profiles, heatmaps, and chromVAR scores  |
| 1i | scMultiome | human CD34-enriched PBMCs | lung cancer and healthy | cell marker plot |
| 3i | scRNA | public data, human NSCLC | pre/post treatment | pre-processing, annotation, and plotting of NRF2 scores and ChIP-X enrichment analysis |
| 4a | scRNA | mouse lung tumor | tumor-bearing | pre-processing, annotation, marker plot |
| 4c-d | scATAC | mouse lung tumor | tumor-bearing | pre-processing, annotation, motif enrichment, market plot | 
| 4e-h | scMultiome | human tumor and PBMCs | NSCLC patients | pre-processing, cell marker plot, motif enrichment, differential analysis, WCGNA on myeloid cells |
| 10a-d | scRNA | public data, human solid tumors | pre/post treatment | pre-processing, annotation, and plotting of ChIP-X enrichment analysis |


# Disclaimer

Due to the stochasticity and differing versions of many single-cell analysis pipelines used, figure panels may not reproduce exactly as depicted in the manuscript. Please contact the study authors if there are any questions regarding specific package installations, input data requirements, data pre-processing, or other concerns that you may have! 

