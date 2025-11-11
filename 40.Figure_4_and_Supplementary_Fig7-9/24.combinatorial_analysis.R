setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Documents/Ronin_chapters/projects/mcmanus_lab_chapters/Analysis/KC018_MCF10A-203_20250417/20.phenotype/20_30_40_50_plates_P2-P5_take2")
library(magrittr)
library(tidyverse)
library(ggrepel)

normDat <- readRDS("normalized_phenotypic_data_take3.rds")
repPlts <- data.frame(plate = c("P2", "P3", "P4", "P5"),
                      rep = c("rep1", "rep1", "rep2", "rep2"))
normDat %<>% merge(., repPlts)
normDat$sgId[(normDat$poolId < 60) &
               startsWith(normDat$sgId, "non-")] <- "non-targeting" 

plToGene2 <- c("ATM_P1P2", "KIF15_P1P2", "KIF22_P1", 
               "KIF22_P2", "KIF4A_P1P2", "NTC_1")

prop <- "totalH2AX"
lowCntNTC <- 
  c("non-targeting_00002_1", "non-targeting_00004_2",
    "non-targeting_00006_3", "non-targeting_00008_4")
res <- list(gene1 = c(), gene2 = c(), statistic = c(),
            parameter = c(), p.value = c(), estimate = c(),
            cntS = c(), cntB = c())
for (sg in unique(normDat$sgId)) {
  for (pl in 1:6) {
    if (is.na(sg) | (sg %in% lowCntNTC)) next
    print(sg)
    s1 <- subset(normDat, (sgId == sg) & (round(poolId/10) == pl),
                 select = prop, drop = TRUE) %>%
      subset(., is.na(.) %>% not())
    s2 <- subset(normDat, startsWith(gene, "non-") &
                   (sgId != sg) & (poolId %in% 61:64), 
                 select = prop, drop = TRUE) %>%
      subset(., is.na(.) %>% not())
    
    if (length(s1) < 100) {
      thisRes <- list(statistic = NA, parameter = NA, p.value = NA, estimate = c(NA, NA))
    } else  thisRes <- t.test(s1, s2)
    
    res$gene1 %<>% c(., sg)
    res$gene2 %<>% c(., plToGene2[pl])
    res$statistic %<>% c(., thisRes$statistic)
    res$parameter %<>% c(., thisRes$parameter)
    res$p.value %<>% c(., thisRes$p.value)
    res$estimate %<>% c(., thisRes$estimate[1] - thisRes$estimate[2])
    res$cntS %<>% c(., length(s1))
    res$cntB %<>% c(., length(s2))
  }
}
res %<>% as.data.frame()
res$FDR <- p.adjust(res$p.value, "BH")
res$FWER <- p.adjust(res$p.value, "bonferroni")

ntcDat <- subset(normDat, startsWith(gene, "non-") & 
                   (poolId %in% 61:64), 
                 select = prop, drop = TRUE) %>% 
  subset(., is.na(.) %>% not())
semNTC <- sd(ntcDat)/sqrt(length(ntcDat))
res$rawEst <- res$estimate
res$estimate <- res$estimate/semNTC
assign(paste0("res_gene_level_", prop), res)

#----------------------------------------------
cst <- list(geneA = c(), geneB = c(),
            statisticA = c(), statisticB = c(), statisticAB = c(),
            parameterA = c(), parameterB = c(), parameterAB = c(),
            cntSA = c(), cntSB = c(), cntSAB = c(),
            cntBA = c(), cntBB = c(), cntBAB = c(),
            p.valueA = c(), p.valueB = c(), p.valueAB = c(),
            estimateA = c(), estimateB = c(), estimateAB = c(),
            rawEstA = c(), rawEstB = c(), rawEstAB = c(),
            FDRA = c(), FDRB = c(), FDRAB = c(),
            FWERA = c(), FWERB = c(), FWERAB = c())
for (gn1 in unique(normDat$sgId)) {
  for (gn2 in plToGene2) {
    print(c(gn1, gn2))
    cst$geneA %<>% c(., gn1)
    cst$geneB %<>% c(., gn2)
    thisA <- paste0("res_gene_level_", prop) %>%
      get() %>%
      subset(., (gene1 == gn1) & (gene2 == "NTC_1"))
    
    #Fix this by doing a collective nonTC vs gene2 in above
    thisB <- paste0("res_gene_level_", prop) %>%
      get() %>%
      subset(., (gene1 %>% startsWith(., "non-")) & (gene2 == gn2))
    
    thisAB <- paste0("res_gene_level_", prop) %>%
      get() %>%
      subset(., (gene1 == gn1) & (gene2 == gn2))
    
    for (cx in c("statistic", "parameter", "p.value", "estimate",
                 "FDR", "FWER", "cntS", "cntB", "rawEst")){
      cst[[paste0(cx, "A")]] %<>% 
        c(., thisA[1, cx, drop = TRUE])
      cst[[paste0(cx, "B")]] %<>% 
        c(., thisB[1, cx, drop = TRUE])
      cst[[paste0(cx, "AB")]] %<>% 
        c(., thisAB[1, cx, drop = TRUE])
    }
  }
}
assign(paste0("cst_gene_level_", prop), as.data.frame(cst))
cst_gene_level_totalH2AX$comboScore <- cst_gene_level_totalH2AX %$%
  (estimateAB - estimateA - estimateB)

cst_gene_level_totalH2AX$rawEstComboScore <- cst_gene_level_totalH2AX %$%
  (rawEstAB - rawEstA - rawEstB)

cst_gene_level_totalH2AX$comboScore[cst_gene_level_totalH2AX$geneB == "NTC_1"] <- NA
#-----------------------------------
filteredCstH2AX <- 
  subset(cst_gene_level_totalH2AX, (estimateA > 0) & (estimateB > 0)  &
           (estimateAB < 0) & (FDRAB < 0.2)) %>%
  mutate(geneAx = strsplit(geneA, "_") %>% 
           plyr::laply(magrittr::extract, i = 1))

toKeep <- filteredCstH2AX[, c("geneAx", "geneB")] %>% 
  group_by(geneAx, geneB) %>% 
  summarize(count = n()) %>%
  subset(count == 2) %>%
  ungroup() %>%
  dplyr::select(geneAx, geneB) %>%
  apply(., 1, paste0, collapse = "_")

filteredCstH2AX %<>% 
  subset(., apply(filteredCstH2AX[, c("geneAx", "geneB")], 
                  1, paste0, collapse = "_") %in% toKeep)

#No KIF22 interactors found with above criterion
save(list = ls() %>% 
       subset(., endsWith(., "_gene_level_totalH2AX")) %>%
       c(., "filteredCstH2AX"),
     file = "combo_analysis_res.Rdata")

