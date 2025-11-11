library(tidyverse)
library(magrittr)
library(rjson)
library(patchwork)

args <- commandArgs(trailingOnly=TRUE)
cyc2 <- args[1]
trackStage <- args[2]
params <- fromJSON(file = args[3])

cyc1 <- params$registrationParams$registrationPairs[[paste0('cycle', cyc2)]]
iter <- params$segmentationParams$segIter
plotWidth <- 10
plotHeight <- 10
trackSuffix <- paste0("_tracking", params$trackingParams[[paste0("stg", trackStage, "OutSuffix")]])
wells <- file.path(params$analysisDir, params$listOfWells) %>%
        readLines()
dat <- map(wells, function(well) {
	paste0(well, ".", cyc1, "_cycle", cyc2, trackSuffix) %>% 
		file.path(params$analysisDir, params$trackObjsDir, .) %>% 
		map(., read.table, header = TRUE, sep = "\t") %>% 
		Reduce(merge, .) %>% 
		apply(., 2, function(x) paste0(well, "n", x)) %>% 
		as.data.frame()
	})

dat %<>% bind_rows()

ftrDir <- file.path(params$analysisDir, params$segmentationDir, 
	paste0("iter", iter, params$segmentationParams$featureDirSuffix)) 
ftrFiles <- list.files(ftrDir)
ftrs <- map(wells, function(well) {
	cycles <- ftrFiles %>% 
		subset(., startsWith(., well) & 
		       	(grepl(str_replace(cyc1, "cycle", "cyc"), .) | 
				grepl(paste0("cyc", cyc2), .)) & 
		       (grepl(str_replace_all(params$segmentationParams$cellposeParams$useStain, 
					      " ", "_"), .)))
	cycles %>% 
		map(., function(x) file.path(ftrDir, x) %>% 
			read.table(., header = TRUE, sep = "\t") %>%
			mutate(label = paste0(well, "n", label), 
				cycle = str_remove_all(x, ".txt") %>% 
					strsplit(., '-cyc') %>%
					unlist() %>% 
					magrittr::extract(., i = 2) %>%
					strsplit(., '[.]') %>%
					unlist() %>% 
					magrittr::extract(., i = 1))) %>% 
		bind_rows()
	})

ftrs %<>% bind_rows()

ornt <- dat
for (cyc in unique(ftrs$cycle)) {
	ornt <- subset(ftrs, cycle == cyc, 
			select = c("label", "orientation")) %>%
		set_colnames(., c("label", paste0("orientation_", paste0("cycle", cyc)))) %>% 
		merge(ornt, ., by.x = paste0("cycle", cyc), by.y = "label") 
}

toPlot <- ornt %$%
	data.frame(cyc12 = get(paste0("orientation_", cyc1)) - 
			get(paste0("orientation_", paste0("cycle", cyc2))))

pdf(paste0("featureDiffs_", cyc1, "_cycle", cyc2, str_remove(trackSuffix, ".txt"), ".", 
	   str_replace_all(Sys.time(), " ", "_"), ".pdf"),
	width = plotWidth, height = plotHeight)
p1 <- ggplot(toPlot, aes(cyc12)) +
	geom_histogram(bins = 100) +
	ggtitle("Orientation difference (in radians)") 
print(p1)

ftrsToCompare <- c("area", "perimeter", "eccentricity", "feret_diameter_max")
for (prop in ftrsToCompare) {

	comp <- dat
	for (cyc in unique(ftrs$cycle)) {
		comp <- subset(ftrs, cycle == cyc, 
				select = c("label", prop)) %>%
			set_colnames(., c("label", paste0(prop, "_cycle", cyc))) %>% 
			merge(comp, ., by.x = paste0("cycle", cyc), by.y = "label") 
	}

	toPlot <- comp %$%
		data.frame(cyc12 = log10(get(paste0(prop, "_", cyc1))/
					get(paste0(prop, "_cycle", cyc2))),
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
for (cyc in unique(ftrs$cycle)) {
	euclidDist <- subset(ftrs, cycle == cyc,
			select = c("label", "centroid.0", "centroid.1")) %>%
		set_colnames(., c("label", paste0(c("centroid.0", "centroid.1"), "_cycle", cyc))) %>%
		merge(euclidDist, ., by.x = paste0("cycle", cyc), by.y = "label")
}

euclidDist$cyc12.0 <- (euclidDist[, paste0("centroid.0_", cyc1)] - 
			euclidDist[, paste0("centroid.0_", paste0("cycle", cyc2))])^2
euclidDist$cyc12.1 <- (euclidDist[, paste0("centroid.1_", cyc1)] - 
			euclidDist[, paste0("centroid.1_", paste0("cycle", cyc2))])^2

euclidDist %<>%
	mutate(cyc12 = (cyc12.0 + cyc12.1) %>% sqrt())
p <- ggplot(euclidDist, aes(cyc12)) +
	geom_histogram(bins = 100) +
	ggtitle("Euclidean distance between mates (in pixels)")
print(p)
dev.off()

