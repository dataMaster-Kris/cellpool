library(tidyverse)
library(magrittr)
library(rjson)

params <- fromJSON(file = "../00.config/config.json")
wells <- file.path(params$analysisDir, params$listOfWells) %>%
        readLines()

load("barcodes_to_cells.rds")
load("features.rds")
phenoDat <- merge(features, bcdCalls[, c("bcd", "objId", "sgId", "gene", "poolId")],
                by.x = "label", by.y = "objId")

save(phenoDat, file = "phenotypic_data.rds")
save(phenoDat, file = "normalized_phenotypic_data.rds")

