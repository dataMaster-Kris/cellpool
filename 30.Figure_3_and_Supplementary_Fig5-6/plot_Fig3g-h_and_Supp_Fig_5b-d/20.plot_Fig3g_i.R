setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Documents/Ronin_chapters/projects/mcmanus_lab_chapters/Analysis/KC018_MCF10A-203_20250417/20.phenotype/10_plate_P1_clean_run")
library(tidyverse)
library(magrittr)
library(ggrepel)

load("normalized_phenotypic_data.rds")
load("a_univariate_diff_analysis_res.RData")

phenoDat$totalH2AX <- log10(phenoDat$area_ch3*phenoDat$intensity_mean_ch3)
phenoDat$totalpH3Ser10 <- log10(phenoDat$area_ch2*phenoDat$intensity_mean_ch2)
phenoDat$totalATub <- log10(phenoDat$intensity_mean_ch1)
phenoDat %<>% subset(., !is.na(gene))

#----------------------------------------
#Nuclear area
labGenes <- "MAPRE2, BIRC5, CDCA8, ESPL1, PLK1, AURKB, INCENP, TTK, PCNT, ORC1, NCAPG2, NCAPH2, SMC4, ECT2, CDC20, CDK7, KIF20A, CETN2, CDK1" %>%
  #NCAPD2, NCAPD3, NCAPH2, SMC4, NCAPG, SMC3, SMC1A, STAG1, STAG2" %>% 
  strsplit(., split = ", ") %>% unlist()
res_gene_level_area_ch4$labText <- ""
res_gene_level_area_ch4$labText[res_gene_level_area_ch4$gene %in% labGenes] <- 
  res_gene_level_area_ch4$gene[res_gene_level_area_ch4$gene %in% labGenes]
res_gene_level_area_ch4$labText[res_gene_level_area_ch4$gene == "non-targeting_00005"] <-
  "non-targeting"
res_gene_level_area_ch4$labText[res_gene_level_area_ch4$gene %>% 
                                  startsWith(., "non-targeting")] <-
  "non-targeting"
res_gene_level_area_ch4$color <- "not significant"
res_gene_level_area_ch4$color[res_gene_level_area_ch4$nonTargeting] <- "non-targeting sgRNA"
res_gene_level_area_ch4$color[(res_gene_level_area_ch4$estimate > 3) &
                                (res_gene_level_area_ch4$FDR < 0.05)] <- "larger nuclei"
res_gene_level_area_ch4$color[(res_gene_level_area_ch4$estimate < -3) & 
                                (res_gene_level_area_ch4$FDR < 0.05)] <- "smaller nuclei"
res_gene_level_area_ch4$FDRcapped <- replace(-log10(res_gene_level_area_ch4$FDR), 
                                             -log10(res_gene_level_area_ch4$FDR) > 6, 6)

p <- ggplot(res_gene_level_area_ch4, aes(x = estimate, y = FDRcapped,
                                         color = color)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             linewidth = 0.1) +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed",
             linewidth = 0.1) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("black", "#f8766d", "#f0f0f0", "black")) +
  scale_size_manual(values = c(3, 6)) +
  theme_classic() +
  geom_text_repel(aes(label = labText), size = 25/16, #25/14) +
                  max.overlaps = 200,
                  min.segment.length = 0, seed = 0291875,
                  segment.size = 0.1, show.legend = FALSE) + 
  guides(shape = "none", 
         color = "none"
  ) +
  theme(legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), 
        legend.box.background = element_rect(colour = "black"),
        legend.margin = margin(0, 3, 0, 3),
        legend.text = element_text(margin = margin(l = 0, r = 0)),
        legend.position.inside = c(0.65, 0.1),
        text = element_text(size = 7),
        plot.margin = margin(t = 0, b = 0, l = 0, r = 0)) +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  ylab(expression(-log[10]*(FDR))) +
  xlab("Nuclear area (scaled and centered wrt non-targeting sgRNAs)") 

ggsave(paste0("area_ch4_ttest_", Sys.time(), ".pdf"), width = 3.4, height = 3.25)

#----------------------------------------
#H2AX vs H3Ser10
h2axCutoff <- res_gene_level_totalH2AX %>%
  subset(., grepl("non-", gene)) %$% range(estimate, na.rm = TRUE)
h3S10Cutoff <- res_gene_level_totalpH3Ser10 %>%
  subset(., grepl("non-", gene)) %$% range(estimate, na.rm = TRUE)
goAnn <- read.table("genes_GO_terms.txt", header = TRUE, sep = "\t", quote = "") %>% 
  subset(., namespace_1003 == "biological_process")
toPlot <- phenoDat %>% 
  group_by(gene) %>% 
  summarize(meanH2AX = mean(totalH2AX),
            meanH3Ser10 = mean(totalpH3Ser10), 
            coverage = n()) %>% 
  mutate(nonTargeting = grepl("non-", gene)) %>% 
  subset(.,  coverage > 100) %>% 
  merge(., res_gene_level_totalH2AX[, c("gene", "estimate", "FDRcapped")] %>% 
          set_colnames(., c("gene", "estimate_H2AX", "FDRcapped_H2AX")), by = "gene", all.x = TRUE) %>% 
  merge(., res_gene_level_totalpH3Ser10[, c("gene", "estimate", "FDRcapped")] %>% 
          set_colnames(., c("gene", "estimate_pH3Ser10", "FDRcapped_H3Ser10")), 
        by = "gene", all.x = TRUE) %>% 
  mutate(labText = ifelse(((FDRcapped_H2AX > -log10(0.05)) &
                             ((FDRcapped_H3Ser10 > -log10(0.05)) & 
                                (!between(estimate_pH3Ser10, h3S10Cutoff[1], h3S10Cutoff[2])) &
                                (coverage > 1000))) |
                            ((estimate_pH3Ser10 > h3S10Cutoff[2]) &
                               (FDRcapped_H3Ser10 > -log10(0.05)) &
                               (coverage > 1000)) |
                            ((estimate_H2AX < -7) & #h2axCutoff[1]) &
                               (FDRcapped_H2AX > -log10(0.05))) |
                            (gene %in% c("BRCA1", "NCAPH", 
                                         "PDS5A", "BUB1", 
                                         "BUB1B", "MAD2L1")), gene, ""))

toPlot$pointCol <- ifelse(toPlot$labText == "",
                          "grey", "black")
toPlot$pointCol[toPlot$nonTargeting] <- "red"
toPlot$labText[toPlot$nonTargeting] <- "non-targeting"

toPlot <- bind_rows(subset(toPlot, !grepl("non-", gene)), subset(toPlot, grepl("non-", gene)))

toPlot <- bind_rows(subset(toPlot, pointCol == "grey"),
                    subset(toPlot, pointCol != "grey"))

toRmv <- "ATM, AURKC, BRCA2, CCNA2, CDKN1B, MYC" %>%
  strsplit(., split = ", ") %>% unlist()
toPlot$labText[toPlot$labText %in% toRmv] <- ""

toPlot$pointCol <- ifelse(toPlot$labText == "",
                          "grey", "black")
toPlot$pointCol[toPlot$nonTargeting] <- "red"

p <- ggplot(toPlot, aes(x = estimate_H2AX, y = estimate_pH3Ser10, 
                        # shape = nonTargeting, 
                        color = pointCol)) +
  geom_point(size = 0.5) +
  theme_classic() + 
  scale_color_manual(values = c("black", "#f0f0f0", "#f8766d")) +
  scale_shape_manual(values = c(16, 4)) + guides(color = "none", shape = "none") +
  geom_text_repel(aes(label = labText), min.segment.length = 0,
                  size = 25/16, max.overlaps = 25, seed = 987752,
                  segment.size = 0.1, box.padding = 0.2, 
                  vjust = 0.25) +
  scale_x_continuous(breaks = c(-20, -15, -10, -5, 0, 5, 10, 15)) +
  theme(axis.title = element_text(size = 5), text = element_text(size = 5),
        plot.margin = margin(t = 0, b = 0, l = 0, r = 0)) +
  xlab(expression(gamma*"H2AX (scaled and centered wrt non-targeting sgRNAs)")) + 
  ylab("H3S10P (scaled and centered wrt non-targeting sgRNAs)") +
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 7))

ggsave(paste0("H2AX_vs_H3Ser10_", Sys.time(), ".pdf"), p, width = 3.4, height = 3.25)


