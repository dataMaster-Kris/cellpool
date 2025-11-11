setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Documents/Ronin_chapters/projects/mcmanus_lab_chapters/Analysis/KC018_MCF10A-203_20250417/20.phenotype/10_plate_P1_clean_run")
library(tidyverse)
library(magrittr)
library(openxlsx)

load("a_univariate_diff_analysis_res.RData")
wb <- createWorkbook("de")
shtNm <- c(totalH2AX = "gamma-H2AX", 
           totalpH3Ser10 = "H3S10P", 
           totalATub = "alpha-tubulin",
           area_ch4 = "nuclear area")
for (prop in c("totalH2AX", "totalpH3Ser10", "totalATub", "area_ch4")) {
  addWorksheet(wb, 
             sheetName=shtNm[prop])
  writeData(wb, sheet = shtNm[prop], paste0("res_gene_level_", prop) %>%
              get() %>% 
              select(-FWER, -nonTargeting, -labText, -FDRcapped) %>%
              arrange(FDR))
}
saveWorkbook(wb, "Supplementary_Data_5.xlsx", overwrite = TRUE)

