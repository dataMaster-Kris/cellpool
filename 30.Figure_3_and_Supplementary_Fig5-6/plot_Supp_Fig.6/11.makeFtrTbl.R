library(tidyverse)
library(magrittr)
library(rjson)

params <- fromJSON(file = "../00.config/config.json")
wells <- file.path(params$analysisDir, params$listOfWells) %>%
        readLines()

allFiles <- list.files("../40.segmentation/iter1.step3.features/", full.names = TRUE) %>%
        subset(., grepl("-cyc2", .) & endsWith(., "_phenotypic.txt"))
features <- map(1:4, function(ch) {
       map(subset(allFiles, grepl(paste0("-cyc2.ch", ch), allFiles)), function(x) {
                well <- basename(x) %>% substr(., 1, 3)
                print(well)
                read.table(x, header = TRUE, sep = "\t") %>%
                       mutate(label = paste0(well, "n", label)) %>%
                       set_colnames(., colnames(.) %>%
                               paste0(., "_ch", ch) %>%
                               replace(., startsWith(., "label_ch"), "label"))
        }) %>%
        bind_rows()
}) %>%
       Reduce(merge, .)

trackedLabels <- map(wells, function(well) {
        print(well)
        file.path(params$analysisDir, params$trackObjsDir) %>%
                        list.files() %>%
                        subset(., startsWith(., well) &
                                endsWith(., paste0("_tracking",
                                        params$debarcodingParams$finalTrackingSuffix))) %>%
                        map(.,
                                function(x) file.path(params$analysisDir,
                                                params$trackObjsDir, x) %>%
                                        read.table(., header = TRUE, sep = "\t") %>%
                                        mutate(well = well)) %>%
                        Reduce(merge, .)
}) %>%
        set_names(wells)

trackedLabels %<>%
        bind_rows() %>%
        mutate(cycle1 = paste0(well, "n", cycle1), cycle2 = paste0(well, "n", cycle2))
features %<>% merge(., trackedLabels[, c("cycle1", "cycle2")], by.x = "label", by.y = "cycle2")
features$label <- features$cycle1

save(features, file = "features.rds")


