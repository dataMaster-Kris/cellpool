library(magrittr)
library(tidyverse)
library(ggrepel)
library(patchwork)

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
fdrCap <- 2.15
res_gene_level_totalH2AX %<>%
  mutate(FDR_capped = -log10(FDR) %>% 
           replace(., is_greater_than(., fdrCap), fdrCap),
         ntcSg = grepl("non", gene1) & (gene2 == "NTC_1"))
res_gene_level_totalH2AX <- 
  subset(res_gene_level_totalH2AX, !ntcSg) %>%
  sample_n(., nrow(.)) %>%
  bind_rows(., subset(res_gene_level_totalH2AX, ntcSg))

res_gene_level_totalH2AX$labText <- ""
res_gene_level_totalH2AX$labText[
  startsWith(res_gene_level_totalH2AX$gene1, "non-") &
    (res_gene_level_totalH2AX$gene2 == "NTC_1")
] <- "NTC"

#-----------------------------------
#Get tables
#-----------------------------------
filteredCstH2AX %<>% 
  bind_rows(., kif22_filteredCstH2AX, 
            subset(cst_gene_level_totalH2AX, 
                   startsWith(geneA, "KIF11_") &
                     (geneB == "KIF15_P1P2")) %>%
              mutate(geneAx = strsplit(geneA, split = "_") %>%
                       plyr::laply(., magrittr::extract, i = 1)))
filteredCstH2AX %<>%
  mutate(group = geneB, 
         comboId = paste0(geneA, "__", geneB)) %>%
  bind_rows(., mutate(., geneB = "NTC_1"))

toPlot2 <- res_gene_level_totalH2AX %>% 
  subset(., paste0(gene1, "__", gene2) %>% 
           is_in(., filteredCstH2AX %$% 
                   c(comboId, 
                     strsplit(comboId, split = "__") %>% 
                       plyr::laply(., magrittr::extract, i = 1) %>%
                       paste0(., "__NTC_1"),
                     "KIF11_3_1__KIF15_P1P2",
                     "KIF11_3_2__KIF15_P1P2"
                   ))
         
  ) %>% 
  subset(., is.na(FDR) %>% not()) %>%
  mutate(comboId = paste0(gene1, "__", gene2))


#Fix this to include plate replication based hits
cst_gene_level_totalH2AX %<>%
  mutate(comboId = paste0(geneA, "__", geneB),
         group = geneB)
toPlot2 <- merge(cst_gene_level_totalH2AX, toPlot2,
                 by = "comboId", all.y = TRUE)

toPlot2$labText <- toPlot2$gene1 %>%
  strsplit(., split = "_") %>%
  plyr::laply(., magrittr::extract, i = 1)

toPlot2$guideNum <- toPlot2$gene1 %>%
  strsplit(., split = "_") %>%
  plyr::laply(., magrittr::extract, i = 3)

toPlot2$labText[toPlot2$geneB == "NTC_1"] <- ""
toPlot2$guideNum[toPlot2$geneB == "NTC_1"] <- ""

toPlot2$group %<>%
  as.character() %>%
  strsplit(., split = "_") %>% 
  plyr::laply(., magrittr::extract, i = 1)

#-----------------------------------
#ATM and KIF4A interactors
#-----------------------------------
toPlot4 <- subset(toPlot2, !(group %in% c("KIF15", "KIF22", 
                                          "NTC")))
kif4a_atm_ntc_dat <- 
  subset(res_gene_level_totalH2AX, 
         (gene1 == "non-targeting") & 
           (gene2 %in% c("KIF4A_P1P2", "ATM_P1P2")))
toPlot4 <- 
  bind_rows(toPlot4, 
            (toPlot4[1:4, ]) %>% 
              mutate(group = rep(kif4a_atm_ntc_dat$gene2 %>%
                                   strsplit(., split = "_") %>%
                                   plyr::laply(., magrittr::extract, i = 1),
                                 2),
                     labText = rep(c("", "non-targeting"), each = 2),
                     estimate = c(0, 0, kif4a_atm_ntc_dat$estimate),
                     FDR_capped = c(0, 0, kif4a_atm_ntc_dat$FDR_capped),
                     gene2 = c("NTC_1", "NTC_1", kif4a_atm_ntc_dat$gene2),
                     geneB = c("NTC_1", "NTC_1", 
                               kif4a_atm_ntc_dat$gene2),
                     estimateAB = c(0, 0,
                                    kif4a_atm_ntc_dat$estimate),
                     estimateA = c(0, 0, 0, 0),
                     FDRAB = c(0, 0, kif4a_atm_ntc_dat$FDR),
                     FDRA = c(1, 1, 1, 1)))
toPlot5 <- subset(toPlot4, geneB != "NTC_1")

set.seed(82223342)
p1 <- ggplot(data = toPlot5 %>% subset(., group == "ATM") %>%
              mutate(xPos = estimateA, yPos = -log10(FDRA))) +
  facet_grid(. ~ group, scales = "free") +
  geom_point(data = res_gene_level_totalH2AX %>%
               subset(., gene2 %in% c("ATM_P1P2", "NTC_1")) %>%
               mutate(xPos = estimate, yPos = FDR_capped),
             mapping = aes(x = xPos, y = yPos),
             size = 0.1, color = "gray", alpha = 0.4,
             inherit.aes = FALSE) +
  coord_cartesian(xlim = c(-4, 4.3), ylim = c(0, 2.15)) +
  geom_point(
    mapping = aes(x = xPos, y = yPos),
    color = "black",
    size = 0.02, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dotted",
             linewidth = 0.1) +
  geom_hline(yintercept = 0.69897, linetype = "dotted",
             linewidth = 0.1) +
  theme_classic() +
  geom_text_repel(
    mapping = aes(label = labText,
                  x = xPos,
                  y = yPos),
    segment.size = 0.2,
    size = 25/14, hjust = -0.2,
    box.padding = 0.1,
    show.legend = FALSE, 
    seed = 108) +
  theme(text = element_text(size = 7)) +
  geom_segment(data = toPlot5  %>% subset(., group == "ATM"), 
               aes(xend = estimateAB, x = estimateA,
                   yend = -log10(FDRAB) %>%
                     replace(., is_greater_than(., fdrCap), 
                             fdrCap),
                   y = -log10(FDRA) %>%
                     replace(., is_greater_than(., fdrCap), fdrCap)),
               inherit.aes = FALSE,
               linewidth = 0.1,
               arrow = arrow(length = unit(0.06, "inches")),
               color = "black") +
  scale_color_manual(values = c("black", "black", "red")) +
  guides(color = "none")

p2 <- ggplot(data = toPlot5 %>% subset(., group == "KIF4A") %>%
               mutate(xPos = estimateA, yPos = -log10(FDRA))) +
  facet_grid(. ~ group, scales = "free") +
  geom_point(data = res_gene_level_totalH2AX %>%
               subset(., gene2 %in% c("KIF4A_P1P2", "NTC_1")) %>%
               mutate(xPos = estimate, yPos = FDR_capped),
             mapping = aes(x = xPos, y = yPos),
             size = 0.1, color = "gray", alpha = 0.4,
             inherit.aes = FALSE) +
  coord_cartesian(xlim = c(-8, 7.1), ylim = c(0, 2.15)) +
  geom_point(
    mapping = aes(x = xPos, y = yPos),
    color = "black",
    size = 0.02, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dotted",
             linewidth = 0.1) +
  geom_hline(yintercept = 0.69897, linetype = "dotted",
             linewidth = 0.1) +
  theme_classic() +
  geom_text_repel(
    mapping = aes(label = labText,
                  x = xPos,
                  y = yPos),
    segment.size = 0.2,
    size = 25/14, hjust = -0.2, 
    box.padding = 0,
    show.legend = FALSE, 
    seed = 1198) + 
  theme(text = element_text(size = 7),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  geom_segment(data = toPlot5  %>% subset(., group == "KIF4A"), 
               aes(xend = estimateAB, x = estimateA,
                   yend = -log10(FDRAB) %>%
                     replace(., is_greater_than(., fdrCap), 
                             fdrCap),
                   y = -log10(FDRA) %>%
                     replace(., is_greater_than(., fdrCap), fdrCap)),
               inherit.aes = FALSE,
               linewidth = 0.1,
               arrow = arrow(length = unit(0.06, "inches")),
               color = "black") +
  scale_color_manual(values = c("black", "black", "red")) +
  guides(color = "none")

p2

p <- p1 | p2

ggsave("ATM_and_KIF4A_interactions_on_volcano.pdf", width = 3.5, height = 2.6)

#-----------------------------------
#KIF15 interactors
#-----------------------------------
kif15_ntc_dat <- 
  subset(res_gene_level_totalH2AX, 
         (gene1 == "non-targeting") & (gene2 == "KIF15_P1P2"))
toPlot2$labText[toPlot2$guideNum == 2] <- ""
toPlot4 <- subset(toPlot2, (group == "KIF15") &
                    (endsWith(geneA, "_1") 
                     | (strsplit(geneA, split = "_") %>%
                          plyr::laply(., magrittr::extract, i = 1) %>%
                          is_in(., c("KIF11"#, #"BRD7",
                                     #"CASC5", "KIFC1"
                          )))
                    )) %>%
  subset(., geneA != "KIF11_3_1") #%>%
toPlot4$labText[startsWith(toPlot4$gene1, "KIF11") &
                  (toPlot4$geneB != "NTC_1")] <- "KIF11"
toPlot4 <- 
  bind_rows(toPlot4, 
            (toPlot4[1:2, ]) %>% 
              mutate(labText = c("", "non-targeting"),
                     estimate = c(0, kif15_ntc_dat$estimate),
                     FDR_capped = c(0, kif15_ntc_dat$FDR_capped),
                     gene2 = c("NTC_1", "KIF15_P1P2"),
                     geneB = c("NTC_1", "KIF15_P1P2"),
                     estimateAB = c(0, kif15_ntc_dat$estimate),
                     estimateA = c(0, 0),
                     FDRAB = c(0, kif15_ntc_dat$FDR),
                     FDRA = c(0, 1)))
toPlot5 <- subset(toPlot4, geneB != "NTC_1")

toPlot5$FDRAB %<>% replace(., is_less_than(., 0.007079458),
                           0.007079458)

set.seed(82223342)
p3 <- ggplot() +
  facet_grid(. ~ group, scales = "free") +
  geom_point(data = res_gene_level_totalH2AX %>%
               subset(., gene2 %in% c("KIF15_P1P2", "NTC_1")) %>%
               mutate(xPos = estimate, yPos = FDR_capped),
             mapping = aes(x = xPos, y = yPos),
             size = 0.1, color = "gray", alpha = 0.4,
             inherit.aes = FALSE) +
  coord_cartesian(xlim = c(-10.5, 6), ylim = c(0, 2.15)) +
  geom_point(data = toPlot5 %>%
               mutate(xPos = estimateA, yPos = -log10(FDRA)),
             mapping = aes(x = xPos, y = yPos),
             color = "black",
             size = 0.02, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dotted",
             linewidth = 0.1) +
  geom_hline(yintercept = 0.69897, linetype = "dotted",
             linewidth = 0.1) +
  theme_classic() +
  geom_text_repel(data = toPlot5 %>%
                    mutate(xPos = estimateA, yPos = -log10(FDRA)), 
                  mapping = aes(label = labText,
                                x = xPos,
                                y = yPos),
                  segment.size = 0.2,
                  size = 25/14, 
                  hjust = -0.2,
                  box.padding = 0.05,
                  show.legend = FALSE, 
                  seed = 1928, inherit.aes = FALSE) +
  theme(text = element_text(size = 7),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  geom_segment(data = toPlot5, 
               aes(xend = estimateAB, x = estimateA,
                   yend = -log10(FDRAB) %>%
                     replace(., is_greater_than(., fdrCap), 
                             fdrCap),
                   y = -log10(FDRA) %>%
                     replace(., is_greater_than(., fdrCap), 
                             fdrCap)),
               inherit.aes = FALSE,
               linewidth = 0.1,
               arrow = arrow(length = unit(0.06, "inches")),
               color = "black") +
  scale_color_manual(values = c("black", "black", "red")) +
  guides(color = "none")
p3

p <- (p1 | p3 | p2) + plot_layout(widths = c(1, 2, 1))

ggsave("all_interactions_with_gene_labels_on_volcano.pdf", width = 7.2, height = 2.6)
