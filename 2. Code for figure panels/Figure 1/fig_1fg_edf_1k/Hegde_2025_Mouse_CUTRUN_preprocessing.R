### Imports

library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(readr)
library(dplyr)
library(tidyr)
library(readxl)
library(GenomicRanges)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BuenRTools)
library(parallel)
library(ggrepel)
library(qs)
library(stringr)
library(fastcluster)
library(readr)

### Setting filepaths

register(MulticoreParam(8, progressbar = TRUE))

DATA_PATH <- file.path("/path/to/input/data")
OUT_PATH <- file.path("/path/to/output/data")

COUNTS_PATH <- file.path(OUT_PATH, "counts_rds")
PLOTS_PATH <- file.path(OUT_PATH, "plots")
DAR_PATH <- file.path(OUT_PATH, "DAR")
CHROMVAR_PATH <- file.path(OUT_PATH, "chromVAR")

setwd(OUT_PATH)

### Processing CUT&RUN data

OCR_PATH <- file.path("/path/to/dependencies/ImmGenATAC18_AllOCRsInfo.csv") # Downloaded from Yoshida et al., 2019 supplementary data (https://www.cell.com/cell/pdf/S0092-8674(18)31650-7.pdf)

df_ocrs_raw <- read_csv(OCR_PATH)

# Establish a ±500bp window around the summit of each OCR
window <- 500

df_ocrs <- df_ocrs_raw %>%
  filter(Included.in.systematic.analysis == 1) %>%
  select(c("chrom", "Summit", "genes.within.100Kb")) %>%
  mutate(start = Summit - 500, end = Summit + 500)

gr_ocrs <- makeGRangesFromDataFrame(
  df_ocrs, keep.extra.columns = TRUE
)

# Set seqlevels/seqinfo for our GRanges so trim works 
seqlevels(gr_ocrs) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(gr_ocrs) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)
gr_ocrs <- trim(gr_ocrs)

gr_ocrs$ucsc_name <- paste0(
  seqnames(gr_ocrs), ":", start(gr_ocrs), "-", end(gr_ocrs)
)

# Get paths of files from metadata
METADATA_PATH <- file.path("/path/to/metadata/metadata.xlsx")

df_metadata <- read_excel(METADATA_PATH)

cell_subsets <- c("Monocyte", "GMP")
labs <- c("Merad", "Buenrostro")

df_metadata_filt <- df_metadata %>%
  filter(
    `Cell Subset` %in% cell_subsets,
    `Lab of Origin` %in% labs
  )

BAM_FILES <- df_metadata_filt %>%
  pull(`Sample Name`)

BAM_PATHS <- list.files(
  path = c(DATA_PATH),
  recursive = TRUE,
  full.names = TRUE
)

BAM_PATHS <- BAM_PATHS[basename(BAM_PATHS) %in% BAM_FILES]

create_counts <- function(path) {
  sample <- getCounts(
    alignment_files = path,
    peaks = gr_ocrs,
    paired = TRUE
  )

  saveRDS(sample, file.path(COUNTS_PATH, paste0(basename(path), ".rds")))

  print("Saved object")

  return(sample)
}

mclapply(BAM_PATHS, create_counts, mc.cores = 8)

# Read files
targets <- c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3")

names(targets) <- targets

count_paths <- list.files(
  path = COUNTS_PATH,
  pattern = "*rds$",,
  full.names = TRUE
)

merge_per_target <- function(target) {
  print(target)

  count_files <- df_metadata_filt %>%
    filter(Target == target) %>%
    mutate(`counts_path` = paste0(`Sample Name`, ".rds")) %>%
    pull(`counts_path`)

  target_merged_obj <- list()

  for (path in count_paths[basename(count_paths) %in% count_files]) {
    target_merged_obj[[basename(path)]] <- readRDS(path)
  }

  target_merged_obj <- do.call(cbind, target_merged_obj)

  merged_metadata <- colData(target_merged_obj) %>%
    merge(df_metadata, by.x = 0, by.y = "Sample Name", sort = FALSE)

  rownames(merged_metadata) <- merged_metadata$Row.names

  colData(target_merged_obj) <- merged_metadata[colnames(target_merged_obj), ]

  target_merged_obj <- filterPeaks(target_merged_obj, non_overlapping = FALSE)
  target_merged_obj <- addGCBias(
    target_merged_obj,
    genome = BSgenome.Mmusculus.UCSC.mm10
  )
}

merged_obj <- mclapply(targets, merge_per_target, mc.cores = length(targets))

saveRDS(merged_obj, file.path(OUT_PATH, paste0(paste0(targets, collapse = "_"), "_merged_obj.rds")))

### chromVAR
set.seed(123)

motif_matching <- function(target){
  print(target)

  target_obj <- merged_obj[[target]]

  motif_ix <- matchMotifs(
    pwms = mouse_pfms_v4,
    subject = target_obj,
    genome = BSgenome.Mmusculus.UCSC.mm10
  )

  saveRDS(motif_ix, file.path(OUT_PATH, paste0(target, "_motif_ix.rds")))

  print("Saved object")

  return(motif_ix)
}

mclapply(targets, motif_matching, mc.cores = length(targets))

compute_deviations <- function(target) {
  print(target)

  target_obj <- merged_obj[[target]]

  motif_ix <- readRDS(file.path(OUT_PATH, paste0(target, "_motif_ix.rds")))

  bg <- getBackgroundPeaks(target_obj, niterations = 250)

  dev_motif <- computeDeviations(
    object = target_obj,
    annotations = motif_ix,
    background_peaks = bg
  )

  saveRDS(dev_motif, file.path(OUT_PATH, paste0(target, "_dev_motif.rds")))

  print("Saved object")

  return(dev_motif)
}

mclapply(targets, compute_deviations, mc.cores = length(targets))

compute_deviations_all_peaks <- function(target) {
  print(target)

  target_obj <- merged_obj[[target]]

  bg <- getBackgroundPeaks(target_obj, niterations = 250)

  dev_peaks <- computeDeviations(
    object = target_obj,
    background_peaks = bg
  )

  saveRDS(
    dev_peaks, file.path(OUT_PATH, paste0(target, "_dev_all_peaks.rds"))
  )

  print("Saved object")

  return(dev_peaks)
}   

mclapply(targets, compute_deviations_all_peaks, mc.cores = length(targets))


### Output top differential chromVAR TFs for each mark
targets <- c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3")

merged_obj <- readRDS(file.path(OUT_PATH, paste0(paste0(targets, collapse = "_"), "_merged_obj.rds")))

cell_subsets <- c("GMP", "Monocyte")

varp <- function(x) {
  return(mean((x - mean(x)) ^ 2))
}

for (cell_subset in cell_subsets) {
  for (target in targets) {
    print(target)

    target_obj <- merged_obj[[target]]

    dev_motif <- readRDS(file.path(OUT_PATH, paste0(target, "_dev_motif.rds")))

    desired_sample_names <- df_metadata_filt %>%
      filter(
        `Cell Subset` == cell_subset,
        `Lab of Origin` == "Merad",
        Target == target
      ) %>%
      pull(`Sample Name`)

    target_obj <- target_obj[, desired_sample_names]
    dev_motif <- dev_motif[, desired_sample_names]

    z_score_motif <- deviationScores(dev_motif)

    rownames(z_score_motif) <- extractTFNames(rownames(z_score_motif))

    df_annot <- data.frame(
      celltype = colData(target_obj)$`Cell Subset`,
      treatment = colData(target_obj)$`Treatment`,
      row.names = colnames(z_score_motif),
      stringsAsFactors = FALSE
    )

    df_annot$treatment[df_annot$treatment == "Kras mutant and p53 mutant"] <- "KP"
    df_annot$treatment[is.na(df_annot$treatment)] <- "Naive"

    group_1 <- "KP"
    group_2 <- "Naive"

    group_1_ind <- (df_annot$treatment == group_1)
    group_2_ind <- (df_annot$treatment == group_2)

    group_1_scores <- z_score_motif[, group_1_ind]
    group_2_scores <- z_score_motif[, group_2_ind]

    group_1_var <- apply(group_1_scores, 1, varp) 
    group_2_var <- apply(group_2_scores, 1, varp) 

    group_1_mean <- apply(group_1_scores, 1, mean)
    group_2_mean <- apply(group_2_scores, 1, mean)

    mean_diff <- group_1_mean - group_2_mean

    var_weighted_avg <- (
      (group_1_var * dim(group_1_scores)[2]) + 
      (group_2_var * dim(group_2_scores)[2])
    ) / dim(z_score_motif)[2]

    z_scored_mean_diff <- mean_diff / sqrt(var_weighted_avg)

    F <- var_weighted_avg / apply(z_score_motif, 1, varp)

    df_results <- tibble(
      mean_diff = mean_diff,
      z_scored_mean_diff = z_scored_mean_diff,
      group_1_mean = group_1_mean,
      group_2_mean = group_2_mean,
      F = F,
      TF = names(mean_diff)
    )

    write.csv(
      df_results,
      file.path(CHROMVAR_PATH, sprintf(
        "Top_TFs_chromVAR_%s_%s_%s_VS_%s.csv",
        cell_subset,
        target,
        group_1,
        group_2
      ))
    )

### Clustering peak chromVAR scores for promoter specific signal, on GMPs and monocytes

targets <- c("H3K4me3")
names(targets) <- targets
cell_subsets <- c("GMP", "Monocyte")
k_clusters <- 24
n_top_variable <- 20000

get_avg_scores <- function(target, n_top_variable) {
  print(target)

  dev_all_peaks <- readRDS(
    file.path(OUT_PATH, paste0(target, "_dev_all_peaks.rds"))
  )

  z_score_motif <- deviationScores(dev_all_peaks)

  z_score_motif_no_na <- na.omit(z_score_motif)
  var <- apply(z_score_motif_no_na, 1, var)

  metadata <- colData(dev_all_peaks)

  colData(dev_all_peaks)$disease <- "KP"

  colData(dev_all_peaks)[
    is.na(colData(dev_all_peaks)$Treatment), "disease"
  ] <- "Naive"

  mask <- (metadata$`Lab of Origin` == "Merad")

  z_score_motif_no_na <- z_score_motif_no_na[!is.na(var), mask]
  dev_all_peaks_no_na <- dev_all_peaks[rownames(z_score_motif_no_na), mask]

  results <- list()

  for (cell_subset in cell_subsets) {
    print(cell_subset)

    cell_subset_mask <- (
      colData(dev_all_peaks_no_na)$`Cell Subset` == cell_subset
    )

    kp_mask <- cell_subset_mask & (colData(dev_all_peaks_no_na)$disease == "KP")
    naive_mask <- cell_subset_mask & (colData(dev_all_peaks_no_na)$disease == "Naive")

    kp_dev <- dev_all_peaks_no_na[, kp_mask]
    naive_dev <- dev_all_peaks_no_na[, naive_mask]

    group_1_scores <- deviationScores(kp_dev)
    group_2_scores <- deviationScores(naive_dev)
    all_scores <- cbind(group_1_scores, group_2_scores)

    group_1_statistic <- apply(group_1_scores, 1, mean)
    group_2_statistic <- apply(group_2_scores, 1, mean)

    all_statistics <- cbind(group_1_statistic, group_2_statistic)

    colnames(all_statistics) <- c(paste0(cell_subset, "_KP_statistic"), paste0(cell_subset, "_Naive_statistic"))

    results[[cell_subset]] <- all_statistics
  }

  df_results <- data.frame(do.call(cbind, results))

  if (!is.null(n_top_variable)) {
    peak_var <- apply(df_results, 1, sd)

    top_variable_index <- order(peak_var, decreasing = TRUE)[1:n_top_variable]

    df_results <- df_results[top_variable_index, ]
  }

  return(df_results)
}

avg_scores <- mcmapply(
  get_avg_scores,
  target = targets,
  n_top_variable = n_top_variable,
  mc.cores = length(targets),
  SIMPLIFY = FALSE
)

cluster_peaks <- function(df_results) {
  print(names(df_results))

  df_results_scaled <- scale(as.data.frame(df_results))

  set.seed(123)
  kmeans_output <- kmeans(
    df_results_scaled,
    centers = k_clusters,
    nstart = k_clusters * 2,
    iter.max = 30
  )

  return(kmeans_output)
}

cluster_outputs <- mcmapply(
  cluster_peaks,
  avg_scores,
  mc.cores = length(avg_scores),
  SIMPLIFY = FALSE
)

output_chosen_cluster <- function(target) {
  df_results <- avg_scores[[target]]
  kmeans_output <- cluster_outputs[[target]]
  chosen_clusters <- chosen_clusters_input[[target]]

  df_output_list <- list()

  for (chosen_cluster in chosen_clusters) {
    desired_peaks <- names(
      kmeans_output$cluster[kmeans_output$cluster %in% chosen_cluster]
    )

    motif_ix <- readRDS(file.path(OUT_PATH, paste0(target, "_motif_ix.rds")))

    desired_ranges <- rowRanges(motif_ix)[as.integer(desired_peaks)]

    desired_ucsc_names <- paste0(
      seqnames(desired_ranges),
      ":",
      start(desired_ranges),
      "-",
      end(desired_ranges)
    )

    gr_ocrs_subset <- values(gr_ocrs)[
      gr_ocrs$ucsc_name %in% desired_ucsc_names,
    ]

    df_output <- cbind(
      df_results[desired_peaks, ],
      gr_ocrs_subset[c("genes.within.100Kb", "ucsc_name")]
    )

    df_output <- df_output[
      order(df_output$GMP_KP_statistic, decreasing = TRUE),
    ]

    df_output_list[[chosen_cluster]] <- df_output
  }

  return(df_output_list)
}

write_outputs <- function(target) {
  chosen_clusters <- chosen_clusters_input[[target]]

  for (chosen_cluster in chosen_clusters) {
    df_output <- data.frame(cluster_output[[target]][[chosen_cluster]])

    write.table(
      df_output,
      file.path(
        DAR_PATH,
        sprintf(
          "%s_peaks_k_means_%s_clustered_grouped_cluster_%s_top_variable_%s.csv",
          target, k_clusters, chosen_cluster, n_top_variable
        )
      ),
      sep = "\t",
      row.names = FALSE
    )

    genes_output <- unique(unlist(strsplit(df_output$genes.within.100Kb, ",")))

    write.table(
      genes_output,
      file.path(
        DAR_PATH,
        sprintf(
          "%s_peaks_k_means_%s_clustered_grouped_cluster_%s_genes_top_variable_%s.csv",
          target, k_clusters, chosen_cluster, n_top_variable
        )
      ),
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )

    ucsc_regions <- df_output$ucsc_name

    # Extract information for GRanges
    regions <- data.frame(
      chr = sapply(strsplit(ucsc_regions, ":"), "[", 1),
      start = str_extract(ucsc_regions, "(?<=:).*?(?=-)"),
      end = sapply(strsplit(ucsc_regions, "-"), "[", 2)
    )

    # Convert rownames into GRanges
    regions <- makeGRangesFromDataFrame(
      regions,
      ignore.strand = TRUE
    )

    df <- data.frame(
      seqnames = seqnames(regions),
      starts = start(regions) - 1,
      ends = end(regions),
      ucsc_regions = ucsc_regions
    )

    write.table(
      df,
      file = file.path(
        DAR_PATH,
        sprintf("%s_peaks_k_means_%s_clustered_grouped_cluster_%s_regions_bed_top_variable_%s.csv", target, k_clusters, chosen_cluster, n_top_variable)
      ),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
  }
}

chosen_clusters_input <- list (
  "H3K4me3" = as.character(unique(cluster_outputs[["H3K4me3"]]$cluster))
)

cluster_output <- mcmapply(
  output_chosen_cluster,
  targets,
  mc.cores = length(targets),
  SIMPLIFY = FALSE
)

written_outputs <- mcmapply(
  write_outputs,
  targets,
  mc.cores = length(targets),
  SIMPLIFY = FALSE
)


### Clustering peak chromVAR signal for other marks, only on GMPs

n_top_variable <- 20000
targets <- c("H3K4me1", "H3K27ac")
names(targets) <- targets
cell_subsets <- c("GMP")
k_clusters <- 8

avg_scores <- mcmapply(
  get_avg_scores,
  target = targets,
  n_top_variable = list(n_top_variable, n_top_variable, n_top_variable),
  mc.cores = length(targets),
  SIMPLIFY = FALSE
)

### Concatenate marks into one matrix
df_H3K4me1 <- as.data.frame(avg_scores["H3K4me1"])
df_H3K27ac <- as.data.frame(avg_scores["H3K27ac"])

motif_H3K4me1 <- readRDS(file.path(OUT_PATH, paste0("H3K4me1", "_motif_ix.rds")))
motif_H3K27ac <- readRDS(file.path(OUT_PATH, paste0("H3K27ac", "_motif_ix.rds")))

peaks_H3K4me1 <- rowRanges(motif_H3K4me1)[as.integer(rownames(df_H3K4me1))]
peaks_H3K27ac <- rowRanges(motif_H3K27ac)[as.integer(rownames(df_H3K27ac))]

intersect_ranges <- Reduce(function(x, y) subsetByOverlaps(x, y, type = "equal"), list(peaks_H3K27ac, peaks_H3K4me1))

common_H3K4me1 <- peaks_H3K4me1 %in% intersect_ranges
common_H3K27ac <- peaks_H3K27ac %in% intersect_ranges

df_avg_scores_concat <- cbind(
  df_H3K4me1[common_H3K4me1, ], 
  df_H3K27ac[common_H3K27ac, ]
)

peak_var <- apply(df_avg_scores_concat, 1, sd)

top_variable_index <- order(peak_var, decreasing = TRUE)[1:n_top_variable]

df_avg_scores_concat <- df_avg_scores_concat[top_variable_index, ]

### Clustering
df_results_scaled <- scale(df_avg_scores_concat)

set.seed(123)
kmeans_output <- kmeans(
  df_results_scaled,
  centers = k_clusters,
  nstart = k_clusters * 2,
  iter.max = 100
)

chosen_clusters <- as.character(unique(kmeans_output$cluster))

# Get cluster info
df_output_list <- list()

for (chosen_cluster in chosen_clusters) {
  print(chosen_cluster)
  desired_peaks <- names(
    kmeans_output$cluster[kmeans_output$cluster %in% chosen_cluster]
  )

  motif_ix <- readRDS(file.path(OUT_PATH, paste0("H3K4me1", "_motif_ix.rds")))

  desired_ranges <- rowRanges(motif_ix)[as.integer(desired_peaks)]

  desired_ucsc_names <- paste0(
    seqnames(desired_ranges),
    ":",
    start(desired_ranges),
    "-",
    end(desired_ranges)
  )

  gr_ocrs_subset <- values(gr_ocrs)[
    gr_ocrs$ucsc_name %in% desired_ucsc_names,
  ]

  df_output <- cbind(
    df_avg_scores_concat[desired_peaks, ],
    gr_ocrs_subset[c("genes.within.100Kb", "ucsc_name")]
  )

  df_output <- df_output[
    order(df_output$ucsc_name, decreasing = FALSE),
  ]

  df_output_list[[chosen_cluster]] <- df_output
}

# Write outputs
n_top_variable <- "20000"
target <- "H3K4me1_H3K27ac_GMP_only_"

for (chosen_cluster in chosen_clusters) {
  df_output <- data.frame(df_output_list[[chosen_cluster]])

  write.table(
    df_output,
    file.path(
      DAR_PATH,
      sprintf(
        "%s_peaks_k_means_%s_clustered_grouped_cluster_%s_top_variable_%s.csv",
        target, k_clusters, chosen_cluster, n_top_variable
      )
    ),
    sep = "\t",
    row.names = FALSE
  )

  genes_output <- unique(unlist(strsplit(df_output$genes.within.100Kb, ",")))

  write.table(
    genes_output,
    file.path(
      DAR_PATH,
      sprintf(
        "%s_peaks_k_means_%s_clustered_grouped_cluster_%s_genes_top_variable_%s.csv",
        target, k_clusters, chosen_cluster, n_top_variable
      )
    ),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  ucsc_regions <- df_output$ucsc_name

  # Extract information for GRanges
  regions <- data.frame(
    chr = sapply(strsplit(ucsc_regions, ":"), "[", 1),
    start = str_extract(ucsc_regions, "(?<=:).*?(?=-)"),
    end = sapply(strsplit(ucsc_regions, "-"), "[", 2)
  )

  # Convert rownames into GRanges
  regions <- makeGRangesFromDataFrame(
    regions,
    ignore.strand = TRUE
  )

  df <- data.frame(
    seqnames = seqnames(regions),
    starts = start(regions) - 1,
    ends = end(regions),
    ucsc_regions = ucsc_regions
  )

  write.table(
    df,
    file = file.path(
      DAR_PATH,
      sprintf("%s_peaks_k_means_%s_clustered_grouped_cluster_%s_regions_bed_top_variable_%s.csv", target, k_clusters, chosen_cluster, n_top_variable)
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
}