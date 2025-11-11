library(tidyverse)
library(magrittr)
library(ggrepel)
library(patchwork)

load("intDat.RData")
load("features.rds")
load("pooled_testing_Results/a_univariate_diff_analysis_res.RData")

testSet <- c("NCAPG2", "TACC1", "KIFC1", "ORC1", "TP53", "RPS6KB1", 
             "ATM", "SNX9", "RPS6KA3", "BIRC5", "CDK1", "MAPRE2", 
             "RCC2", "AURKB", "POGZ", "UHRF1") 

cor_pool_arr <- c()
ftr2Title <- data.frame(prop = c("area_ch4", 
                                 "totalH2AX", "totalpH3Ser10"),
                        title = c(NA,
                                  "Nuclear gamma-H2AX signal (aribtrary units)",
                                  "Nuclear H3S10P signal (aribtrary units)"))
fdrBreaks <- list(area_ch4 = c(1, 5, 10, 15),
                  totalH2AX = c(1, 5, (1:5)*10),
                  totalpH3Ser10 = c(1, 5, 10, 20, 30))
for (prop in c("area_ch4", 
               "totalH2AX", "totalpH3Ser10")) {
  
  ntcCntPool <- paste0("res_gene_level_", prop) %>%
    get() %$% 
    max(cntB)
  poolRes <- paste0("res_gene_level_", prop) %>%
    get() %>% 
    subset(., gene %in% testSet)
  arrayRes <- paste0("res_gene_level_wrt_untransduced_", prop) %>%
    get() %>% 
    subset(., gene %in% c(testSet, "NTC"))
  
  mergeRes <- merge(bind_rows(poolRes[, c("gene", "rawEst", "FDR", "cntS")],
                              data.frame(gene = "NTC", rawEst = 0, FDR = NA,
                                         cntS = ntcCntPool)),
                    arrayRes[, c("gene", "rawEst", "FDR", "cntS")], 
                    by = "gene", suffixes = c("_pool", "_array"))
  
  assign(paste0("mergeRes_", prop), mergeRes)
  
  toPlot <- mergeRes[order(mergeRes$FDR_pool), ]
  toPlot %<>% subset(., ((FDR_pool < 0.05) & (FDR_array < 0.05)) |
                       (gene == "NTC"))
  toPlot$gene %<>% factor(., levels = rev(.))
  poolDat <- toPlot[, c("gene", "rawEst_pool", "FDR_pool", "cntS_pool")] %>% 
    set_colnames(., c("gene", "estimate", "FDR", "cntS")) %>% 
    mutate(set = "arrayed\nminipools")
  arrDat <- toPlot[, c("gene", "rawEst_array", "FDR_array", "cntS_array")] %>% 
    set_colnames(., c("gene", "estimate", "FDR", "cntS")) %>% 
    mutate(set = "arrayed\nindividuals")
  toPlot <- bind_rows(poolDat, arrDat) %>% 
    mutate(set = factor(set, levels = c("arrayed\nminipools",
                                        "arrayed\nindividuals")))
  if (prop == "area_ch4") toPlot$estimate <- toPlot$estimate*0.36
  
  p1 <- ggplot(toPlot, 
               aes(x = 0, y = gene, fill = estimate)) +
    facet_grid(. ~ set) +
    geom_tile() +
    scale_fill_gradient2(low = "#fc8d59", high = '#5ab4ac') +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()
          ) +
    xlab(NULL) + ylab(NULL) +
    guides(fill = guide_legend(title = "effect size"))
  
  if (prop == "area_ch4") p1 <- p1 + 
    scale_fill_gradient2(low = "#fc8d59", high = '#5ab4ac',
                         limits = c(-25, 100))
  
  p0 <- ggplot(toPlot %>% subset(set == "arrayed\nminipools"),
               aes(x = 0, y = gene, size = -log10(FDR))) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "left",
          axis.line = element_blank()
          ) +
    xlab(NULL) + ylab(NULL) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_size_area(breaks = fdrBreaks[[prop]],
                    limits = c(1, -log10(min(mergeRes$FDR_pool)))) +
    guides(size = guide_legend(title = expression(-log["10"]*"(FDR)"))) 
  
  yRange <- toPlot %>% 
    subset(set == "arrayed\nminipools") %$%
    max(cntS, na.rm = TRUE) %>% 
    divide_by(1000) %>% 
    floor()
  p2 <- ggplot(toPlot %>% subset(set == "arrayed\nminipools"),
               aes(x = gene, y = cntS/1000)) +
    geom_col() +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
    ) + coord_flip() +
    xlab(NULL) + ylab(expression("Coverage (10"^3*")")) +
    scale_y_continuous(
                       expand = c(0, 0),
                       breaks = 0:yRange,
                       labels = c("0", rep("", yRange - 1), 
                                  as.character(yRange))) +
    guides(size = guide_legend(title = expression(-log["10"]*"(FDR)"))) 
  
  p <- (p0 | p2 | p1) + plot_layout(widths = c(0.3, 0.5, 2), guides = "collect") +
    plot_annotation(title = ifelse(prop == "area_ch4",
                                   expression("Nuclear area (in "*mu*"m"^2*")"),
                                   ftr2Title$title[ftr2Title$prop == prop]),
                    theme = theme(plot.title = element_text(hjust = 0.5),
                                  plot.margin = margin(c(0,0,0,0), unit = "cm")))
  
  assign(paste0("p_", prop), p)
  
  pdf(paste0("Compare_pooled_vs_individuals_", 
             prop, ".pdf"), width = 4.6, height = 7)
  print(p)
  dev.off()
  
}
