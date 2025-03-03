library(UCell)
library(qs)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
BPPARAM <- BiocParallel::MulticoreParam(workers = 16, progressbar = TRUE)

setwd("/folder/containing/genesets/")

sce_pbmc <- qread(
  "[/path/to/PBMC/object]"
)

sce_tumor <- qread(
  "[/path/to/tumor/object]"
)

genesets_paths <- list.files("genesets", full.names = TRUE)

get_file_name <- function(x) {
  return(tools::file_path_sans_ext(basename(x)))
}

get_gene_set <- function(x) {
  gene_list <- as.list(read.table(x, header = TRUE))
  return(lapply(gene_list, toupper))
}

genesets <- lapply(genesets_paths, get_gene_set)

names(genesets) <- lapply(genesets, names)

rowData(sce_pbmc) <- NULL
rowData(sce_tumor) <- NULL

patients_keep <- c("Lu927", "Lu941", "Lu954", "Lu979")
pbmc_celltype_cols <- c("Level 2", "sample_id")
tumor_celltype_cols <- c("Level 3", "sample_id")

sce_pbmc_mask <- (
  grepl(
    pattern = paste0(patients_keep, collapse = "|"),
    x = sce_pbmc$sample_id
  )
)

sce_tumor_mask <- (
  grepl(
    pattern = paste0(patients_keep, collapse = "|"),
    x = sce_tumor$sample_id
  ) & (
    sce_tumor$disease == "Tumor"
  )
)

common_genes <- intersect(rownames(sce_pbmc), rownames(sce_tumor))

celltype_col <- "celltype"
tissue_col <- "tissue"
sample_col <- "sample_id"
new_metadata_cols <- c(celltype_col, sample_col)

sce_pbmc_subset <- sce_pbmc[common_genes, sce_pbmc_mask]
colData(sce_pbmc_subset) <- colData(sce_pbmc_subset)[pbmc_celltype_cols]
colnames(colData(sce_pbmc_subset)) <- new_metadata_cols
sce_pbmc_subset$tissue <- "PBMC"

sce_tumor_subset <- sce_tumor[common_genes, sce_tumor_mask]
colData(sce_tumor_subset) <- colData(sce_tumor_subset)[tumor_celltype_cols]
colnames(colData(sce_tumor_subset)) <- new_metadata_cols
sce_tumor_subset$tissue <- "Tumor"

sce <- cbind(sce_pbmc_subset, sce_tumor_subset)

names(assays(sce)) <- c("counts")

sce <- ScoreSignatures_UCell(
  sce,
  features = genesets,
  assay = "counts",
  name = NULL,
  BPPARAM = BPPARAM
)

df_scores <- as.data.frame(t(assay(altExp(sce, "UCell"))))

df_scores[, new_metadata_cols] <- as.matrix(
  colData(sce)[rownames(df_scores), new_metadata_cols]
)

df_scores[, tissue_col] <- as.matrix(
  colData(sce)[rownames(df_scores), tissue_col]
)

tib_scores <- tibble(df_scores)

# https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(
      data,
      xminv <- x - violinwidth * (x - xmin),
      xmaxv = x + violinwidth * (xmax - x)
    )
    grp <- data[1, "group"]
    newdata <- plyr::arrange(
      transform(data, x = if (grp %% 2 == 1) xminv else xmaxv),
      if (grp %% 2 == 1) y else -y
    )
    newdata <- rbind(
      newdata[1, ],
      newdata,
      newdata[nrow(newdata), ],
      newdata[1, ]
    )
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(
      newdata[1, "x"]
    )

    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[
        rep(1, nrow(quantiles)),
        setdiff(names(data), c("x", "y")),
        drop = FALSE
      ]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname(
        "geom_split_violin",
        grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob)
      )
    } else {
      ggplot2:::ggname(
        "geom_split_violin", GeomPolygon$draw_panel(newdata, ...)
      )
    }
  }
)

geom_split_violin <- function(
  mapping = NULL,
  data = NULL,
  stat = "ydensity",
  position = "identity",
  ...,
  draw_quantiles = NULL,
  trim = TRUE,
  scale = "area",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm,
      ...
    )
  )
}

plot_celltype <- "CD14+ monocyte"
colors <- c("PBMC" = "black", "Tumor" = "black")
fills <- colors <- c("PBMC" = "black", "Tumor" = "tan2")
order <- c("PBMC", "Tumor")
patient_col <- "patient_id"

# grouping by patient with pairing from HSPC -> CD14+ monocyte
plot_patient <- function(tib) {
  p <- tib %>%
    filter(.data[[celltype_col]] %in% plot_celltype) %>%
    mutate(!!patient_col := str_sub(.data[[sample_col]], end = -2)) %>%
    select(-c(!!celltype_col, !!sample_col)) %>%
    group_by(.data[[patient_col]], .data[[tissue_col]]) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(
      cols = names(genesets),
      names_to = "pathway",
      values_to = "pathway_mean"
    ) %>%
    ggplot(
      aes(
        x = factor(.data[[tissue_col]], levels = order),
        y = pathway_mean,
        fill = .data[[tissue_col]]
      )
    ) +
    geom_boxplot() +
    geom_line(aes(group = .data[[patient_col]]), color = "black") +
    geom_point(size = 2) +
    theme_classic() +
    xlab(NULL) +
    ylab("Normalized averaged patient score") +
    scale_fill_manual(values = colors) +
    facet_wrap(~pathway, scales = "free", ncol = 2) +
    theme(
      strip.text = element_text(size = 4),
      axis.line.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.ticks.x = element_blank()
    )

  return(p)
}

ggsave(
  "[/path/to/your/output].pdf",
  plot_patient(tib_scores),
  width = 5,
  height = 17
)

