library(dplyr)
library(broom)
library(survival)
library(extrafont)
library(ggplot2)
library(ggrepel)
font_import() 
loadfonts()    
theme_set(theme_minimal(base_family = "Arial"))

setwd("../yourpath")

enr_ChIPX_up_results <- read.csv('../yourpath/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.Human.enrichr.reports.txt', sep = "\t")
enr_ChIPX_up_results <- enr_ChIPX_up_results %>%
  mutate(label = sapply(strsplit(as.character(Term), " "), `[`, 1))
print(enr_ChIPX_up_results)

enr_ChIPX_up_results %>%
  subset(Term != "TCF3 ENCODE") -> enr_ChIPX_up_results

ggplot(enr_ChIPX_up_results, aes(y = Odds.Ratio, x = -log10(Adjusted.P.value))) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "lightgrey", linewidth = 0.4) +
  geom_point(aes(color = Adjusted.P.value < 0.05), size = 2) + 
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  geom_text_repel(
    data = subset(enr_ChIPX_up_results, label %in% c("SMAD4", "NFE2L2", "TP63", "PPARG", "AR",'SOX2','GATA1','GATA2','TCF3','UBTF','ESR1','SALL4')), 
    aes(label = label),
    nudge_y = -0.25,
    nudge_x = 0.35,
    color = "black",
    size = 5.3,
    box.padding = 0.35,       
    point.padding = 0.3,      
    segment.color = "black"
  ) +
  geom_point(data = subset(enr_ChIPX_up_results, label %in% c("NFE2L2")), 
             aes(x = -log10(Adjusted.P.value), y = Odds.Ratio), 
             color = "red", size = 3) +
  labs(
    title = "Hu et al., 2023 (PMID: 36869384)\nNR vs R to NAIT in NSCLC\nChEA TF regulators in mo-macs",
    y = "Odds ratio",
    x = expression(~-log[10](p.adj))
  ) +
  theme_classic(base_size = 15) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 17, color = "black",family = "Arial"),
    axis.title.x = element_text(size = 18,family = "Arial"),
    axis.title.y = element_text(size = 18,family = "Arial"),
    axis.line.x.bottom = element_line(linewidth = 0.4),
    axis.line.y.left = element_line(linewidth = 0.4),
    legend.position = "none" 
  )

ggsave(paste0("../yourpath.pdf"), plot = last_plot(), dpi = 300, width = 4.5, height = 5)

