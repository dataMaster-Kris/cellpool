library(tidyverse)
library(magrittr)
library(rjson)
library(patchwork)

args <- commandArgs(trailingOnly=TRUE)
cyc2 <- "ashlar"
params <- fromJSON(file = args[1])

cyc1 <- "cellpool"
iter <- params$segmentationParams$segIter
nghbrsToConsider <- params$trackingParams$nNgbrsToSavePerImage
trackDir <- file.path(params$analysisDir, params$ashlarAnalysisDir, 
	"45.track_cellpool_to_ashlar")
wells <- file.path(params$analysisDir, params$listOfWells) %>%
	readLines()

trackDat <- map(wells, function(w) {
        file.path(trackDir, paste0(w, ".", cyc1, "_", cyc2, 
			"_tracking", params$trackingParams$stg1OutSuffix)) %>%
                read.table(., sep = "\t", header = TRUE) %>%
                mutate(well = w)
}) %>%
        bind_rows()

cycFtrs <- map(c("cellpool", "ashlar"), function(tool) {
	map(wells, function(w) {
		if (tool == "cellpool") {
			fDir <- params$analysisDir
			cyc <- paste0("-cyc1.", 
				str_replace_all(params$segmentationParams$cellposeParams$useStain, 
					" ", "_"),
				".txt")
			sep <- "\t"
		}
		if (tool == "ashlar") {
			fDir <- file.path(params$analysisDir, params$ashlarAnalysisDir)
			cyc <- paste0("_", params$mosaickingParams$ashlarParamSuffix, ".ch0.csv")
			sep <- ","
		}
                file.path(fDir, params$segmentationDir, 
			paste0("iter", iter, params$segmentationParams$featureDirSuffix), 
			paste0(w, cyc)) %>%
                read.table(., sep = sep, header = TRUE) %>%
                mutate(label = paste0(w, "n", label))
        }) %>% bind_rows()
}) %>% set_names(c("cellpool", "ashlar"))

cellLocs <- map(cycFtrs, function(x) {
        x[, c("label", "centroid.0", "centroid.1")] %>%
                mutate(well = substr(label, 1, 3),
                        objId = strsplit(label, "n") %>%
                                plyr::laply(., magrittr::extract, i = 2) %>%
                                as.numeric())
})

#---------------------
intDat <- trackDat 
intDat[, cyc1] %<>% paste0(intDat$well, "n", .)
intDat[, cyc2] %<>% paste0(intDat$well, "n", .)
intDat %<>%
        merge(., cellLocs[[cyc1]][, c("label", "centroid.0", "centroid.1")],
                by.x = cyc1, by.y = "label") %>%
        merge(., cellLocs[[cyc2]][, c("label", "centroid.0", "centroid.1")],
                by.x = cyc2, by.y = "label", 
		suffixes = paste0(".", c(cyc1, cyc2)))

intDat$dispY <- intDat[, paste0("centroid.0.", cyc2)] - 
			intDat[, paste0("centroid.0.", cyc1)]
intDat$dispX <- intDat[, paste0("centroid.1.", cyc2)] - 
			intDat[, paste0("centroid.1.", cyc1)]
intDat$dispDir <- intDat$dispY/intDat$dispX
intDat$disp <- intDat[, c("dispX", "dispY")] %>% 
	`^`(., 2) %>%
	apply(., 1, function(x) sum(x) %>% sqrt())

#--------------------------
#Find nearest neighbors
cpTrackDir <- file.path(params$analysisDir, params$trackObjsDir)
ngbrs <- map(wells, function(w) {
        file.path(cpTrackDir, paste0(w, "-cyc1", 
		"_nearest_nghbrs_stats.txt")) %>%
                read.table(., sep = "\t", header = TRUE) %>%
                mutate(well = w)
}) %>%
        bind_rows()

#---------------------------
v <- merge(intDat, apply(ngbrs[, grepl("ngbrLabel", colnames(ngbrs))], 2,
                        function(x) paste0(ngbrs$well, "n", x)),
        by.x = cyc1, by.y = "ngbrLabel0", all.y = TRUE)

for (nxbr in 1:(nghbrsToConsider - 1)) {
        print(nxbr)
        thisNxbr <- merge(v[, c(cyc1, paste0("ngbrLabel", nxbr))],
                        v[, c(cyc1, "dispX", "dispY", "dispDir", "disp")],
                        by.x = paste0("ngbrLabel", nxbr), by.y = cyc1) %>%
                select(all_of(cyc1), "dispX", "dispY", "dispDir", "disp") %>%
                set_colnames(c(cyc1, paste0(c("dispXNgbr", "dispYNgbr", 
				"dispDirNgbr", "dispNgbr"), nxbr)))
        v <- merge(v, thisNxbr, all.x = TRUE)
}

v$expDispX <- apply(v[, grepl("dispXNgbr", colnames(v))], 1,
        median, na.rm = TRUE)
v$expDispY <- apply(v[, grepl("dispYNgbr", colnames(v))], 1,
        median, na.rm = TRUE)
v$expDispDir <- apply(v[, grepl("dispDirNgbr", colnames(v))], 1,
	median, na.rm = TRUE)
v$expDisp <- apply(v[, grepl("dispNgbr", colnames(v))], 1,
        median, na.rm = TRUE)

untracked <- setdiff(ngbrs %$% paste0(well, "n", label), 
			paste0(trackDat$well, "n", trackDat[, cyc1]))
whichDispIncnstntWrtLocal <- v[(abs(atan(v$dispDir) - atan(v$expDispY/v$expDispX)) %>%
                log10() %>% is_greater_than(-1)) | 
		(abs(log10(v$disp/v$expDisp)) > 0.05) | 
		(v[, cyc1] %in% untracked), c(cyc1, "expDispX", "expDispY")] 
whichDispIncnstntWrtLocal %<>% 
		subset(., !is.na(expDispX))
map(wells, function(w) {
	this <- whichDispIncnstntWrtLocal %>%
		subset(., get(cyc1) %>% startsWith(., w))
	this[, cyc1] %<>% str_remove_all(., paste0(w, "n"))

	write.table(this, 
			file = paste0(w, ".", cyc1, "_", cyc2, 
						"_potential_errors", 
						params$trackingParams$stg1OutSuffix),
			row.names = FALSE, quote = FALSE, sep = "\t")
}) %>% unlist()

toPlot <- data.frame(dispDir = atan(v$dispDir) - atan(v$expDispY/v$expDispX),
        disp = v$disp/v$expDisp)
p1 <- ggplot(toPlot, aes(log10(abs(dispDir)))) +
        geom_histogram(bins = 100) +
        geom_vline(xintercept = -0.5, linetype = "dotted", color = "red") +
	xlab("Log10 of absolute difference in\ndisplacement direction (radians)")
p2 <- ggplot(toPlot, aes(abs(log10(disp)))) +
        geom_histogram(bins = 100) +
        scale_y_log10() +
	xlab("Absolute value of Log10 of ratio of displacement")

ggsave(paste0(cyc1, "_", cyc2, "_stg1.consistency_of_mate_displacement_wrt_local_average.", 
	str_replace_all(Sys.time(), " ", "_"), ".png"), p1 / p2)

