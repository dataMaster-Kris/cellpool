library(tidyverse)
library(magrittr)
library(rjson)

params <- fromJSON(file = "../00.config/config.json")
wells <- file.path(params$analysisDir, params$listOfWells) %>%
        readLines()
dat <- map(wells, function(well) {
        print(well)
        cycFiles <- file.path(params$analysisDir, params$barcodeAnalysisDir) %>%
                        list.files() %>%
                        subset(., startsWith(., well) & grepl("-cyc", .) & endsWith(., ".csv"))

        uniqCycs <- strsplit(cycFiles, split = "-cyc") %>%
                plyr::laply(., magrittr::extract, i = 2) %>%
                str_remove_all(".csv") %>%
                unique()

        map(uniqCycs, function(cyc) {
                thisFiles <- cycFiles %>%
                        subset(., grepl(paste0("-cyc", cyc), .))
                map(thisFiles, function(f) {
                        file.path(params$analysisDir, params$barcodeAnalysisDir, f) %>%
                                read.table(., header = TRUE, sep = ",") %>%
                                select(label,
                                        paste0("q", params$debarcodingParams$useQntl4Ftrs)) %>%
                                set_colnames(c("label", str_remove(f, ".csv")))
                }) %>%
                        Reduce(merge, .)
        }) %>%
                set_names(uniqCycs)
}) %>%
        set_names(wells)

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
                                        read.table(., header = TRUE, sep = "\t")) %>%
                        Reduce(merge, .)
}) %>%
        set_names(wells)

intDat <- map(wells, function(well) {
        print(well)
        thisDat <- trackedLabels[[well]]
        for (cyc in colnames(trackedLabels[[well]])) {
                thisDat <- merge(thisDat, dat[[well]][[str_remove_all(cyc, "cycle")]],
                                by.x = cyc, by.y = "label")
        }

        cycCh <- colnames(thisDat) %>%
                        subset(., grepl("-cyc", .))
        colNames <- cycCh %>%
                str_remove_all(., paste0(well, "-")) %>%
                map_chr(., function(x) {
                        x <- strsplit(x, "-cyc") %>%
                                unlist() %>%
                                str_remove_all(., "C")
                        paste0("cycle", x[2], ".ch", x[1])
                }) %>%
                c("objId", .)

        thisDat <- thisDat[, c("cycle1", cycCh)] %>%
                set_colnames(., colNames) %>%
                mutate(objId = paste0(well, "n", objId))
        thisDat
}) %>%
        bind_rows()

intDat <- subset(intDat, apply(intDat[, -1], 1, function(x) !any(is.na(x))))

outFile <- gzfile("integrated_intensity_table.txt.gz", "w")
write.table(format(intDat, digits = ceiling(log10(2^params$bitDepth)) + 3, nsmall = 3),
        outFile,
        sep = "\t", row.names = FALSE, quote = FALSE)
close(outFile)

ftrs <- map(wells, function(w) {
        file.path(params$analysisDir, params$segmentationDir,
                paste0("iter", params$segmentationParams$segIter,
                        params$segmentationParams$featureDirSuffix),
                paste0(w, "-cyc1.",
                        str_replace_all(params$segmentationParams$cellposeParams$useStain,
                                " ", "_"), ".txt")) %>%
                read.table(., header = TRUE, sep = "\t") %>%
                mutate(label = paste0(w, "n", label))
}) %>%
        bind_rows()

ftrFile <- gzfile("feature_profiles.txt.gz", "w")
write.table(ftrs, ftrFile,
        sep = "\t", row.names = FALSE, quote = FALSE)
close(ftrFile)

tcMap <- file.path(params$analysisDir, params$tcMapFile) %>%
        read.table(., header = TRUE, sep = "\t")
intDat$well <- substr(intDat$objId, 1, 3)
intDat <- merge(intDat, tcMap[, c("well", "poolId", "treatment")])

gfp_pos_thresh <- 11
gfp_neg_thresh <- 8

pdf(paste0("GFP_intensity_histogram_",
        str_replace_all(Sys.time(), " ", "_"), ".pdf"),
        height = 5, width = 7)
p <- ggplot(intDat, aes(log2(cycle1.ch2))) +
	geom_histogram(bins = 100) +
	coord_cartesian(xlim = c(0, params$bitDepth)) +
	xlab("log2(intensity)") +
	ylab("Count") + 
	geom_vline(xintercept = gfp_pos_thresh, color = "red", linetype = "dotted") +
	geom_vline(xintercept = gfp_neg_thresh, color = "red", linetype = "dotted") +
	ggtitle("GFP")

print(p)
dev.off()

intDat$GFP_state <- NA
intDat$GFP_state <- replace(intDat$GFP_state, intDat$`cycle1.ch2` > 2^gfp_pos_thresh, 1)
intDat$GFP_state <- replace(intDat$GFP_state, intDat$`cycle1.ch2` < 2^gfp_neg_thresh, 0)

ls() %>% setdiff(., c("outFile", "p", "ftrFile")) %>% 
        save(file = "intDat.RData", list = .)





