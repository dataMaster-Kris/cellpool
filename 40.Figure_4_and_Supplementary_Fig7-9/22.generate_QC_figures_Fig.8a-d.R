setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Documents/Ronin_chapters/projects/mcmanus_lab_chapters/Analysis/KC018_MCF10A-203_20250417/20.phenotype/20_30_40_50_plates_P2-P5_take2")
library(tidyverse)
library(magrittr)

normDat <- readRDS("normalized_phenotypic_data_take3.rds")

FOI <- c(area = "area_ch4", hoechst = "totalHoechst",
         gh2ax = "totalH2AX",
         ph3ser10 = "totalpH3Ser10", 
         aTubulin = "totalATub")

normDat <- subset(normDat, !is.na(gene))

repPlts <- data.frame(plate = c("P2", "P3", "P4", "P5"),
                      rep = c("rep1", "rep1", "rep2", "rep2"))

label <- c(area = "Nuclear\narea", hoechst = "Integrated nuclear\nHoechst 33342",
           gh2ax = "Integrated nuclear\ngamma-H2AX",
           ph3ser10 = "Integrated nuclear\nH3S10P", 
           aTubulin = "Mean cytoplasmic\nalpha-tubulin")
for (ftr in names(FOI)) {
  med_ftr <- normDat %>%
    subset(., 
           (gene != "non-targeting_00002")) %>%
    subset(.,
           (gene != "non-targeting_00004")) %>%
    subset(.,
           (gene != "non-targeting_00006")) %>%
    subset(.,
           (gene != "non-targeting_00008")) %>%
    group_by(gene, plate, poolId) %>%
    summarize(med = get(FOI[ftr]) %>% 
                median(., na.rm = TRUE), count = n())
  med_ftr$comboId <- apply(med_ftr[, c("gene", "poolId")],
                           1, paste0, collapse = "_") 
  med_ftr %<>% 
    subset(., comboId != "NA_NA")
  med_ftr %<>% merge(., repPlts)
  
  med_ftr$med[med_ftr$count < 100] <- NA
  toPlot <- subset(med_ftr, count >= 100, 
                   select = c(-count, 
                              -gene,
                              -poolId,
                              -plate)) %>% 
    spread(., rep, med) %>%
    mutate(color = c("test", "NTC")[1 + startsWith(comboId, "non-")])
  
  toPlot$plotTitle <- label[ftr]
  
  p <- ggplot(toPlot, aes(x = rep1, 
                          y = rep2)) +
    geom_point(size = 1) + 
    facet_grid(. ~ plotTitle) +
    theme_classic() +
    coord_fixed() +
    geom_abline(intercept = 0, slope = 1, color = "red",
                linetype = "dotted") +
    theme(text = element_text(size = 7), 
          plot.margin = margin(0, 0, 0, 0, "cm"),
          panel.spacing = unit(0, "cm"))
  
  ggsave(paste0("Supplementary_Figure_7_", ftr, ".pdf"), 
         width = 2, height = 2)  
  
}
