library(tidyverse)
library(magrittr)
library(rjson)
library(patchwork)

args <- commandArgs(trailingOnly=TRUE)
cyc2 <- "ashlar"
trackStage <- args[1]
params <- fromJSON(file = args[2])

cyc1 <- "cellpool"
iter <- params$segmentationParams$segIter
plotWidth <- 10
plotHeight <- 10
trackSuffix <- paste0("_tracking", params$trackingParams[[paste0("stg", trackStage, "OutSuffix")]])
wells <- file.path(params$analysisDir, params$listOfWells) %>%
        readLines()
dat <- map(wells, function(well) {
	paste0(well, ".", cyc1, "_", cyc2, trackSuffix) %>% 
		file.path(params$analysisDir, params$ashlarAnalysisDir,
			"45.track_cellpool_to_ashlar", .) %>% 
		map(., read.table, header = TRUE, sep = "\t") %>% 
		Reduce(merge, .) %>% 
		apply(., 2, function(x) paste0(well, "n", x)) %>% 
		as.data.frame()
	})

dat %<>% bind_rows()

ftrs <- map(c("cellpool", "ashlar"), function(tool) {
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
                mutate(label = paste0(w, "n", label), tool = tool)
        }) %>% bind_rows()
}) %>% bind_rows()

ornt <- dat
for (tl in unique(ftrs$tool)) {
	ornt <- subset(ftrs, tool == tl, 
			select = c("label", "orientation")) %>%
		set_colnames(., c("label", paste0("orientation_", tl))) %>% 
		merge(ornt, ., by.x = tl, by.y = "label") 
}

toPlot <- ornt %$%
	data.frame(cyc12 = get(paste0("orientation_", cyc1)) - 
			get(paste0("orientation_", cyc2)))

pdf(paste0("featureDiffs_", cyc1, "_", cyc2, str_remove(trackSuffix, ".txt"), ".", 
	   str_replace_all(Sys.time(), " ", "_"), ".pdf"),
	width = plotWidth, height = plotHeight)
p1 <- ggplot(toPlot, aes(cyc12)) +
	geom_histogram(bins = 100) +
	ggtitle("Orientation difference (in radians)") 
print(p1)

ftrsToCompare <- c("area")
for (prop in ftrsToCompare) {

	comp <- dat
	for (tl in unique(ftrs$tool)) {
		comp <- subset(ftrs, tool == tl, 
				select = c("label", prop)) %>%
			set_colnames(., c("label", paste0(prop, "_", tl))) %>% 
			merge(comp, ., by.x = tl, by.y = "label") 
	}

	toPlot <- comp %$%
		data.frame(cyc12 = log10(get(paste0(prop, "_", cyc1))/
					get(paste0(prop, "_", cyc2))),
			row = substr(get(cyc1), 1, 1) %>% 
				map_int(., function(x) equals(LETTERS, x) %>% which()),
			col = substr(get(cyc1), 2, 3) %>% as.numeric())
	p1 <- ggplot(toPlot, aes(cyc12)) +
		geom_histogram(bins = 100) +
		ggtitle(prop) +
		facet_grid(row ~ col) #Turn off this layer to get a single histogram for the plate
	print(p1)
}

euclidDist <- dat
for (tl in unique(ftrs$tool)) {
	euclidDist <- subset(ftrs, tool == tl,
			select = c("label", "centroid.0", "centroid.1")) %>%
		set_colnames(., c("label", paste0(c("centroid.0", "centroid.1"), "_", tl))) %>%
		merge(euclidDist, ., by.x = tl, by.y = "label")
}

euclidDist$cyc12.0 <- (euclidDist[, paste0("centroid.0_", cyc1)] - 
			euclidDist[, paste0("centroid.0_", cyc2)])^2
euclidDist$cyc12.1 <- (euclidDist[, paste0("centroid.1_", cyc1)] - 
			euclidDist[, paste0("centroid.1_", cyc2)])^2

euclidDist %<>%
	mutate(cyc12 = (cyc12.0 + cyc12.1) %>% sqrt())
p <- ggplot(euclidDist, aes(cyc12)) +
	geom_histogram(bins = 100) +
	ggtitle("Euclidean distance between mates (in pixels)")
print(p)
dev.off()

