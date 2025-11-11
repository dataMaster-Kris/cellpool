library(tidyverse)
library(magrittr)
library(ggrepel)

load("intDat.RData")
load("features.rds")

phenoDat <- merge(features, intDat[, c("well", "objId", "treatment", "GFP_state")],
                  by.x = "label", by.y = "objId")
phenoDat$col <- substr(phenoDat$well, 2, 3)
phenoDat$row <- substr(phenoDat$well, 1, 1) 

phenoDat$totalHoechst <- log10(phenoDat$area_ch4 * phenoDat$intensity_mean_ch4)
phenoDat$totalH2AX <- log10(phenoDat$area_ch3*phenoDat$intensity_mean_ch3)
phenoDat$totalpH3Ser10 <- log10(phenoDat$area_ch2*phenoDat$intensity_mean_ch2)
phenoDat$totalATub <- log10(phenoDat$intensity_mean_ch1)

#-----------------------------------------
phenoDat$treatment %<>% replace(., startsWith(., "NTC"), "NTC")
testSet <- c("NCAPG2", "TACC1", "KIFC1", "ORC1", "TP53", "RPS6KB1", 
             "ATM", "SNX9", "RPS6KA3", "BIRC5", "CDK1", "MAPRE2", 
             "RCC2", "AURKB", "POGZ", "UHRF1") 

for (prop in c("area_ch4", 
               "totalH2AX", "totalpH3Ser10")) {
  res <- list(gene = c(), statistic = c(), parameter = c(), 
              p.value = c(), estimate = c(), cntS = c(), cntB = c())
  for (sg in c(testSet, "NTC")) { 
    if (is.na(sg)) next
    print(sg)
    s1 <- subset(phenoDat, startsWith(treatment, sg) & (GFP_state == 1), 
                 select = c(prop, "well"), drop = TRUE) %>% 
      subset(., get(prop) %>% is.na() %>% not())
    
    thisWells <- subset(phenoDat, startsWith(treatment, sg) & 
                          (GFP_state == 1)) %$%
      unique(well)
    
    if (prop == "area_ch") {
      s2 <- subset(phenoDat, (startsWith(treatment, "NTC")) &
                     (GFP_state == 1),
                   select = prop, drop = TRUE)
    } else {
      s2 <- subset(phenoDat, (GFP_state == 0) &
                     (well %in% thisWells),
                   select = c(prop, "well"), drop = TRUE) %>%
        subset(., get(prop) %>% is.na() %>% not())
      
      for (nxWl in thisWells) {
        wlCnt <- min(subset(s1, well == nxWl) %>% nrow(),
                     subset(s2, well == nxWl) %>% nrow())
        
        if ((subset(s1, well == nxWl) %>% nrow()) > wlCnt) {
          s1 <- bind_rows(subset(s1, well != nxWl),
                          subset(s1, well == nxWl) %>% 
                            sample_n(., size = wlCnt))
        }
        
        if ((subset(s2, well == nxWl) %>% nrow()) > wlCnt) {
          s2 <- bind_rows(subset(s2, well != nxWl),
                          subset(s2, well == nxWl) %>% 
                            sample_n(., size = wlCnt))
        }
      }
      s2 <- s2[, prop, drop = TRUE]
    }
    
    s1 <- s1[, prop, drop = TRUE]
    
    if (length(s1) < 100) {
      thisRes <- list(statistic = NA, parameter = NA, p.value = NA, estimate = c(NA, NA))
    } else  thisRes <- t.test(s1, s2)
    
    res$gene %<>% c(., sg)
    res$statistic %<>% c(., thisRes$statistic)
    res$parameter %<>% c(., thisRes$parameter)
    res$p.value %<>% c(., thisRes$p.value)
    res$estimate %<>% c(., (thisRes$estimate[1] - thisRes$estimate[2]))
    res$cntS %<>% c(., length(s1)) 
    res$cntB %<>% c(., length(s2)) 
  }
  res %<>% as.data.frame()
  res$FDR <- p.adjust(res$p.value, "BH")
  res$FWER <- p.adjust(res$p.value, "bonferroni")
  
  res$labText <- res$gene
  
  res$FDRcapped <- -log10(res$FDR)
  res$FDRcapped[res$FDRcapped > 50] <- 50
  
  ntcDat <- subset(phenoDat, (startsWith(treatment, "NTC")) &
                     (GFP_state == 1),
                   select = prop, drop = TRUE) %>% subset(., is.na(.) %>% not())
  semNTC <- sd(ntcDat, na.rm = TRUE)/sqrt(
    length(ntcDat %>% subset(., is.na(.) %>% not())))
  res$rawEst <- res$estimate 
  res$estimate <- (res$estimate - res$estimate[res$gene == "NTC"])/semNTC
  assign(paste0("res_gene_level_wrt_untransduced_", prop), res)
  
}



