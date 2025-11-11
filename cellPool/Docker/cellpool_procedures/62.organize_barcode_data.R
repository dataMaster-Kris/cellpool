library(tidyverse)
library(magrittr)
library(rjson)

args <- commandArgs(trailingOnly = TRUE)
params <- fromJSON(file = args[1])
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
#-------------------------------------
#Imaging metadata
#-------------------------------------
chToElem <- file.path(params$analysisDir, params$chToElemId) %>% 
		read.table(., header = TRUE, sep = "\t")
chToStain <- file.path(params$analysisDir, params$chToStainId) %>%
		read.table(., header = TRUE, sep = "\t")
bcdElems <- file.path(params$analysisDir, params$barcodes) %>% 
		read.table(., header = TRUE, sep = "\t") %>%
		select(colnames(.) %>% subset(., grepl("epi", .))) %>% 
		unlist() %>% 
		unique() %>% 
		c(params$debarcodingParams$backboneProt, .)

elemToCh <- map(bcdElems, function(i) apply(chToElem, 1, function(x) i %in% x) %>%
        magrittr::extract(chToElem$cycle, .) %>%
        paste0(., ".ch", apply(chToElem, 2, function(x) i %in% x) %>%
                magrittr::extract(colnames(chToElem), .) %>%
                substr(., nchar(.), nchar(.)))) %>%
        set_names(bcdElems)

chToElem <- list()
for (elem in names(elemToCh)) {
        chToElem[[elemToCh[[elem]]]] <- elem
}

chToElem <- data.frame(channel = names(chToElem), element = unlist(chToElem))
chToElem$fluorophore <- map_chr(chToElem$channel, function(x) {
	chToStain[map_lgl(chToStain$cycle, function(y) startsWith(x, y)), 
		strsplit(x, ".ch") %>% unlist() %>% `[`(., 2) %>% paste0("Channel.", .)]
})
chToElem <- chToElem[order(chToElem$channel), ]
write.table(chToElem, "chToElem.txt", 
	sep = "\t", row.names = FALSE, quote = FALSE)

tcMap <- file.path(params$analysisDir, params$tcMapFile) %>%
        read.table(., header = TRUE, sep = "\t")

#-------------------------------------
#-------------------------------------
#Plot histograms of each cycle intensity
#-------------------------------------
#-------------------------------------
intDat$well <- substr(intDat$objId, 1, 3)
intDat <- merge(intDat, tcMap[, c("well", "poolId")])

toPlot <- (intDat[, chToElem$channel %>% c(., "poolId")]) %>%
        set_colnames(paste0(chToElem$element, "\n", chToElem$fluorophore) %>% c(., "poolId"))

pdf(paste0("PoolId_level_intensity_histograms_",
        str_replace_all(Sys.time(), " ", "_"), ".pdf"),
        height = 5, width = 7)
for (elem in setdiff(colnames(toPlot), "poolId")) {
	p <- ggplot(toPlot, aes(get(elem) %>% log2())) +
                geom_histogram(bins = 100) +
                facet_grid(poolId ~ .) +
		coord_cartesian(xlim = c(0, params$bitDepth)) +
		xlab("log2(intensity)") +
		ylab("Count") +
		ggtitle(elem)

	print(p)
}
dev.off()

toPlot <- (intDat[, chToElem$channel]) %>%
        set_colnames(paste0(chToElem$element, "\n", chToElem$fluorophore)) %>%
        gather(Element, Intensity)

toPlot$Element %<>% factor(.,
        levels = chToElem[order(chToElem$channel), ] %$%
                paste0(element, "\n", fluorophore))

p <- ggplot(toPlot, aes(log2(Intensity))) +
                geom_histogram(bins = 100) +
                facet_wrap(vars(Element), nrow = 2) +
        coord_cartesian(xlim = c(0, params$bitDepth)) +
        xlab("log2(intensity)") +
        ylab("Count")

ggsave(paste0("Population_level_intensity_histograms_",
        str_replace_all(Sys.time(), " ", "_"), ".png"),
        height = 5, width = 7, p)



