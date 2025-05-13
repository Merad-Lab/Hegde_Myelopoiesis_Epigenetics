library(UCell)
library(qs)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
BPPARAM <- BiocParallel::MulticoreParam(workers = 32, progressbar = TRUE)

setwd("/folder/containing/genesets/")

sce <- qread("[/path/to/object]")

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

names(assays(sce)) <- c("counts")

sce <- ScoreSignatures_UCell(
  sce,
  features = genesets,
  assay = "counts",
  name = NULL,
  BPPARAM = BPPARAM
)

df_scores <- as.data.frame(t(assay(altExp(sce, "UCell"))))
df_scores$`Level 2` <- colData(sce)[rownames(df_scores), "Level 2"]
df_scores$disease <- colData(sce)[rownames(df_scores), "disease"]
df_scores$sample_id <- colData(sce)[rownames(df_scores), "sample_id"]
df_scores$prog_enrich <- colData(sce)[rownames(df_scores), "prog_enrich"]

tib_scores <- tibble(df_scores)
tib_scores_subset <- tibble(df_scores_subset)

# https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

celltype_col <- "Level 2"
colors <- c("HD" = "gray", "LUAD" = "red")
order <- c("HSPC", "LMPP", "CD14+ monocyte")

plot_violinplot_hd_vs_luad <- function(tib) {
  p <- tib %>%
    filter(`Level 2` %in% order) %>%
    filter(prog_enrich == "Yes") %>%
    pivot_longer(
      cols = names(genesets),
      names_to = "pathway",
      values_to = "Normalized UCell score"
    ) %>%
    ggplot(
      aes(
        x = factor(.data[[celltype_col]], levels = order),
        y = `Normalized UCell score`,
        fill = disease
      )
    ) +
    geom_split_violin(scale = "width") +
    scale_fill_manual(values = colors) +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    ylab("Normalized UCell score") +
    facet_wrap(~pathway, scales = "free", ncol = 2) +
    theme(strip.text = element_text(size = 7))

  return(p)
}

ggsave(
  "[/path/to/your/output/]",
  plot_violinplot_hd_vs_luad(tib_scores),
  width = 8,
  height = 16
)