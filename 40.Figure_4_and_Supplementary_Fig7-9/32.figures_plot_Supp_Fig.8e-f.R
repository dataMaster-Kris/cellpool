setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Documents/Ronin_chapters/projects/mcmanus_lab_chapters/Analysis/KC018_MCF10A-203_20250417/20.phenotype/dbl_pert_validations_integrated/b_arrayed_minipool_dbls_summary")
library(magrittr)
library(tidyverse)
library(ggrepel)

load("comboAnalysisWithNullDists.RData")
normDat <- readRDS("../../20_30_40_50_plates_P2-P5_take2/normalized_phenotypic_data_take3.rds")
repPlts <- data.frame(plate = c("P2", "P3", "P4", "P5"),
                      rep = c("rep1", "rep1", "rep2", "rep2"))
normDat %<>% merge(., repPlts)

plToGene2 <- c("ATM_P1P2", "KIF15_P1P2", "KIF22_P1", 
               "KIF22_P2", "KIF4A_P1P2", "NTC_1")

lowCntNTC <- 
  c("non-targeting_00002_1_1", "non-targeting_00004_2_1",
    "non-targeting_00006_3_1", "non-targeting_00008_4_1")
res_gene_level_totalH2AX %<>%
  subset(., !(gene1 %in% lowCntNTC))

#------------------------------------------
#-------------------------------------
#Explore agreement of H2AFX between guides in the same context
toPlot <- list(sg1 = c(), sg2 = c(), estimate1 = c(), 
               estimate2 = c(),
               count1 = c(), count2 = c(), arrGene = c())

for (gn2 in plToGene2) {
  thisRes <- res_gene_level_totalH2AX %>% 
    subset(., gene2 == gn2)
  
  for (gn1 in unique(normDat$gene)) {
    thisSgs <- subset(thisRes, startsWith(gene1, gn1))
    if (length(thisSgs) < 2) next
    
    toPlot$sg1 %<>% c(., thisSgs$gene1[1])
    toPlot$sg2 %<>% c(., thisSgs$gene1[2])
    toPlot$estimate1 %<>% c(., thisSgs$estimate[1])
    toPlot$estimate2 %<>% c(., thisSgs$estimate[2])
    toPlot$count1 %<>% c(., thisSgs$count1[1])
    toPlot$count2 %<>% c(., thisSgs$count1[2])
    toPlot$arrGene %<>% c(., gn2)
    toPlot$poolGene %<>% c(., gn1)
  }
}
toPlot %<>% as.data.frame()
toPlot$labText <- 
  replace(toPlot$poolGene, toPlot$poolGene != "H2AFX", "")
# toPlot$labText[toPlot$labText == "H2AFX"] <- "H2AX" 

toPlot$arrGene %<>% replace(., equals(., "NTC_1"), 
                            "non-targeting")
set.seed(24253098)
pdf("Guide_level_estimate_comparison_H2AFX.pdf",
    width = 2.61, height = 1.52)
p <- ggplot(toPlot[sample(nrow(toPlot), nrow(toPlot)), ], 
            aes(x = estimate1, y = estimate2, color = arrGene)) +
  geom_point(size = 0.3#, shape = 1, stroke = 0.2
  ) +
  theme_classic() +
  geom_text_repel(aes(label = labText), size = 25/16, max.overlaps = 50,
                  show.legend = FALSE, segment.size = 0.1) +
  # xlab("Guide 1") +
  # ylab("Guide 2") +
  theme(
    # legend.title = element_blank(),
    # legend.direction = "vertical",
    # legend.position = "right",
    text = element_text(size = 7), 
    # legend.key.size = unit(0.5, "lines"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    # legend.text = element_text(margin = margin(l = 2, r = 0,
    #                                            angle = 60)),
    # legend.spacing = unit(0, "pt"),
    axis.text = element_text(color = "black")) +
  scale_color_brewer(palette = "Dark2") +
  xlab(NULL) + ylab(NULL) +
  coord_fixed(ratio = 1) +
  guides(color = "none") 

p
dev.off()

pdf("Guide_level_estimate_comparison_H2AFX_with_legend.pdf",
    width = 5, height = 5)
p <- ggplot(toPlot[sample(nrow(toPlot), nrow(toPlot)), ], 
            aes(x = estimate1, y = estimate2, color = arrGene)) +
  geom_point(size = 0.3#, shape = 1, stroke = 0.2
  ) +
  theme_classic() +
  geom_text_repel(aes(label = labText), size = 25/16, max.overlaps = 50,
                  show.legend = FALSE, segment.size = 0.1) +
  # xlab("Guide 1") +
  # ylab("Guide 2") +
  theme(
    legend.title = element_blank(),
    legend.direction = "horizontal",
    legend.position = "top",
    text = element_text(size = 7), 
    legend.key.size = unit(0.5, "lines"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(margin = margin(l = 2, r = 0)),
    legend.spacing = unit(0, "pt"),
    axis.text = element_text(color = "black")) +
  scale_color_brewer(palette = "Dark2") +
  xlab(NULL) + ylab(NULL)

p
dev.off()
#----------------------------------
#----------------------------------
fdrCap <- 6
res_gene_level_totalH2AX %<>%
  mutate(FDR_capped = -log10(FDR) %>% 
           replace(., is_greater_than(., fdrCap), fdrCap),
         ntcSg = grepl("non", gene1) & (gene2 == "NTC_1"))
res_gene_level_totalH2AX <- 
  subset(res_gene_level_totalH2AX, !ntcSg) %>%
  sample_n(., nrow(.)) %>%
  bind_rows(., subset(res_gene_level_totalH2AX, ntcSg))

# res_gene_level_totalH2AX$labText <- ""
# res_gene_level_totalH2AX$labText[
#   (res_gene_level_totalH2AX$gene1 %>% startsWith(., "STAG2_")) &
#     (res_gene_level_totalH2AX$gene2 %>% startsWith(., "ATM_")) 
# ] <- "STAG2"

res_gene_level_totalH2AX$labText <- ""
res_gene_level_totalH2AX$labText[
  startsWith(res_gene_level_totalH2AX$gene1, "non-") &
    (res_gene_level_totalH2AX$gene2 == "NTC_1")
] <- "non-targeting"
p <- ggplot(res_gene_level_totalH2AX, 
            aes(x = estimate, y = FDR_capped, color = gene2))  +
  geom_vline(xintercept = 0, linetype = "dotted") +
  # geom_hline(yintercept = 0.69897, linetype = "dotted") +
  geom_point(alpha = 1, size = 0.6) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.direction = "horizontal",
    legend.position = "top",
    text = element_text(size = 7), 
    legend.key.size = unit(0.5, "lines"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(margin = margin(l = 2, r = 0)),
    legend.spacing = unit(0, "pt"),
    axis.text = element_text(color = "black")) +
  geom_text_repel(mapping = aes(label = labText),
                  max.overlaps = 1000, size = 25/12, 
                  show.legend = FALSE, segment.size = 0.1) +
  guides(color = "none") +
  xlab(NULL) + ylab(NULL) #+
  # geom_point(data = subset(res_gene_level_totalH2AX, ntcSg) %>% 
  #              select(estimate, FDR_capped), # %>% 
  #            # bind_rows(., data.frame(estimate = -25,
  #            #                         FDR_capped = 0.6)),
  #            aes(x = estimate, y = FDR_capped), shape = 10,
  #            color = "black") 
# geom_text_repel(size = 5, max.overlaps = 1000)
p
ggsave("Volcano_plot_all.pdf", height = 4, width = 4)

