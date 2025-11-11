setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Documents/Ronin_chapters/projects/mcmanus_lab_chapters/Analysis/KC018_MCF10A-203_20250417/20.phenotype/20_30_40_50_plates_P2-P5_take2")
library(tidyverse)
library(magrittr)

phenoDat <- readRDS("integrated_phenoDat.rds")
phenoDat$well <- substr(phenoDat$label, 1, 3)
phenoDat$column <- substr(phenoDat$label, 2, 3)
phenoDat %<>% subset(., !(well %in% c("H02", "H09"))) #wt-only wells

phenoDat$totalHoechst <- log10(phenoDat$area_ch4 * phenoDat$intensity_mean_ch4)
phenoDat$totalH2AX <- log10(phenoDat$area_ch3*phenoDat$intensity_mean_ch3)
phenoDat$totalpH3Ser10 <- log10(phenoDat$area_ch2*phenoDat$intensity_mean_ch2)
phenoDat$totalATub <- log10(phenoDat$intensity_mean_ch1)

FOI <- c(area = "area_ch4", hoechst = "totalHoechst",
         gh2ax = "totalH2AX",
         ph3ser10 = "totalpH3Ser10", 
         aTubulin = "totalATub")

#----------------------------------
#Normalize
normRefs <- phenoDat %>% 
  subset(., is.na(gene)) %>% 
  group_by(column, plate) %>% 
  summarize(med_area = median(area_ch4, na.rm = TRUE),
            mad_area = mad(area_ch4, na.rm = TRUE),
            med_hoechst = median(totalHoechst, na.rm = TRUE),
            mad_hoechst = mad(totalHoechst, na.rm = TRUE),
            med_gh2ax = median(totalH2AX, na.rm = TRUE),
            mad_gh2ax = mad(totalH2AX, na.rm = TRUE),
            med_ph3ser10 = median(totalpH3Ser10, na.rm = TRUE),
            mad_ph3ser10 = mad(totalpH3Ser10, na.rm = TRUE),
            med_aTubulin = median(totalATub, na.rm = TRUE),
            mad_aTubulin = mad(totalATub, na.rm = TRUE))

normDat <- list()
for (col in unique(phenoDat$column)) {
  for (plt in unique(phenoDat$plate)) {
    thisRef <- subset(normRefs, (plate == plt) & (column == col))
    thisTbl <- phenoDat[(phenoDat$plate == plt) & (phenoDat$column == col), ]
    
    for (ftr in names(FOI)) {
      if (ftr == "area") next
      print(c(ftr, col, plt))
      thisTbl[[FOI[ftr]]] <- 
        (thisTbl[[FOI[ftr]]] - 
           as.numeric(thisRef[, paste0("med_", ftr)]))/as.numeric(thisRef[, paste0("mad_", ftr)])
    }
    normDat[[paste0(plt, "_", col)]] <- thisTbl
  }
}
normDat %<>% bind_rows()

checkNorm <- normDat %>% 
  subset(., is.na(gene)) %>% 
  group_by(column, plate) %>% 
  summarize(med_area = median(area_ch4, na.rm = TRUE),
            mad_area = mad(area_ch4, na.rm = TRUE),
            med_hoechst = median(totalHoechst, na.rm = TRUE),
            mad_hoechst = mad(totalHoechst, na.rm = TRUE),
            med_gh2ax = median(totalH2AX, na.rm = TRUE),
            mad_gh2ax = mad(totalH2AX, na.rm = TRUE),
            med_ph3ser10 = median(totalpH3Ser10, na.rm = TRUE),
            mad_ph3ser10 = mad(totalpH3Ser10, na.rm = TRUE),
            med_aTubulin = median(totalATub, na.rm = TRUE),
            mad_aTubulin = mad(totalATub, na.rm = TRUE))

repPlts <- data.frame(plate = c("P2", "P3", "P4", "P5"),
                      rep = c("rep1", "rep1", "rep2", "rep2"))

#----------------------------
#Plot coverage
#----------------------------
comboFreq <- phenoDat %>% 
  subset(., bcd != "None") %>% #remove untransduced cells
  dplyr::select(sgId, poolId) %>% 
  mutate(comboId = paste0(sgId, poolId)) %$% 
  table(comboId) %>% 
  data.frame()
p <- ggplot(comboFreq, aes(log10(Freq))) +
  geom_histogram(bins = 30) +
  theme_classic() +
  xlab(NULL) + ylab(NULL) +
  theme(axis.text = element_text(size = 5),
        plot.margin = margin(0,0,0,0)) +
  coord_cartesian(xlim = c(1, 4)) +
  scale_x_continuous(breaks = c(1:4)) +
  geom_vline(xintercept = median(log10(comboFreq$Freq)),
             linetype = "dashed", color = "red")

ggsave("figures/Counts_per_guideCombo_histogram.pdf", p, width = 1.34, height = 1.21)

#----------------------------------------

#Outlier to remove (see the rep1 to rep2 comparison for aTubulin)
#subset(toPlot, rep2 < -0.8) => non-targeting_00002_21
#Very low counts for "non-targeting_00002", "non-targeting_00004",
#... "non-targeting_00006", "non-targeting_00008" across minipools.

normDat <- subset(normDat, !is.na(gene))

repPlts <- data.frame(plate = c("P2", "P3", "P4", "P5"),
                      rep = c("rep1", "rep1", "rep2", "rep2"))
saveRDS(normDat[, colnames(normDat) %>% 
                  subset(., grepl("_ch", .) %>% not()) %>% 
                  c("area_ch4")],
        "normalized_phenotypic_data_take3.rds")


