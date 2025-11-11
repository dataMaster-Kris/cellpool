setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Documents/Ronin_chapters/projects/mcmanus_lab_chapters/Analysis/KC018_MCF10A-203_20250417/20.phenotype/20_30_40_50_plates_P2-P5")
library(tidyverse)
library(magrittr)

for (ix in 2:5) {
  load(paste0("phenotypic_data_P", ix, ".rds"))
  phenoDat$plate <- paste0("P", ix)
  assign(paste0("phenoDat_P", ix), phenoDat)
  rm(phenoDat)
}

phenoDat <- bind_rows(phenoDat_P2, phenoDat_P3,
                      phenoDat_P4, phenoDat_P5)
saveRDS(phenoDat, "integrated_phenoDat.rds")


