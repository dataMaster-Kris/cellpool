library(magrittr)
library(tidyverse)
library(rjson)

args <- commandArgs(trailingOnly = TRUE)
params <- fromJSON(file = args[1])
outFile <- paste0('mosaicQC.', str_replace_all(Sys.time(), " ", "_"), ".pdf")
mosaicDir <- file.path(params$analysisDir, params$mosaicDir)
wells <- readLines(file.path(params$analysisDir, params$listOfWells))
chToStain <- file.path(params$analysisDir, params$chToStainId) %>%
	read.table(., sep = "\t", header = TRUE) %>% 
	gather("ch", "coordSource", -cycle) %>% 
	subset(., !is.na(coordSource)) %>% 
	mutate(ch = str_replace_all(ch, "hannel.", ""))
#--------------------------------------------
#Heatmap of SEMs for all tile coordinates in all cycles
#--------------------------------------------
filesToRead <- list.files(mosaicDir) %>%
	subset(., endsWith(., ".tileCentroidCoords.txt")) %>%
	data.frame(fileName = ., well = substr(., 1, 3))	
dat <- map(wells, function(w) {
	subset(filesToRead, well == w) %$% 
	map(fileName, function(cyc) {
		file.path(mosaicDir, cyc) %>%
			read.table(., sep = "\t", header = TRUE) %>%
			add_column(well = w, 
				cycle = str_remove_all(cyc, ".tileCentroidCoords.txt") %>% 
						strsplit(., split = "-cyc") %>% 
						unlist() %>% 
						magrittr::extract(., i = 2) %>%
						paste0("cycle", .)) %>%
			mutate(gridX = divide_by(medX, abs(medX) %>% 
					subset(., is_greater_than(., 
							params$stitchingParams$maxLtrlShift)) %>% 
							min(.)) %>% 
					round(),
				gridY = divide_by(medY, abs(medY) %>%
					subset(., is_greater_than(., 
							params$stitchingParams$maxLtrlShift)) %>%
							min(.)) %>%
					round())
	}) %>% 
		bind_rows()
}) %>% 
	bind_rows()

dat$row <- dat$well %>% 
	map_int(., function(x) substr(x, 1, 1) %>% 
		equals(LETTERS, .) %>% which())
dat$col <- dat$well %>% substr(., 2, 3) %>% as.numeric()

nGridLenInWell <- params$stitchingParams$nTilesPerAxisInPalette
maxGridPt <- floor(nGridLenInWell/2)
dat$gridX %<>% 
	replace(., is_greater_than(., maxGridPt), maxGridPt)
dat$gridX %<>% 
	replace(., equals(., -maxGridPt), -maxGridPt)
dat$gridY %<>% 
	replace(., equals(., maxGridPt), maxGridPt)
dat$gridY %<>% 
	replace(., equals(., -maxGridPt), -maxGridPt)

toPlot <- dat %>%
	mutate(seX = stdX/sqrt(nPaths), seY = stdY/sqrt(nPaths)) %>% 
	select(cycle, seX, seY, gridX, gridY, row, col) %>% 
	gather(imgAxis, SEM, -cycle, -gridX, -gridY, -row, -col)
toPlot$imgAxis %<>% str_remove_all("se")
nGridLenInWell <- toPlot %$% 
	subset(gridX, !is.nan(gridX)) %>% #Assuming equal radii along the X and Y axes
	unique() %>% 
	length()
pdf(outFile, height = unique(toPlot$cycle) %>% 
		length() %>%
		multiply_by(., 5), width = 15)
ggplot(toPlot, aes(x = col + (1/nGridLenInWell)*gridX, 
		y = -row + (1/nGridLenInWell)*gridY, fill = SEM)) + #row numbers increase downwards
	theme_classic() + 
	geom_tile() + 
	scale_fill_continuous(low = "white", high = "steelblue") +
	facet_grid(cycle ~ imgAxis) +
	theme(axis.line = element_blank(),
		axis.text = element_blank(), 
		axis.ticks = element_blank(),
		axis.title = element_blank())
dev.off()

#--------------------------------------------
#Which channel informed the estimates of tile coordinates
#--------------------------------------------

p <- ggplot(dat, aes(x = medX, y = medY,
                shape = cycle, color = coordSource)) +
        geom_point() +
	scale_shape_manual(values = c(0:25)[1:length(unique(dat$cycle))]) +
        facet_grid(row ~ col) +
        theme_classic() +
	theme(axis.text.x = element_text(angle = 90),
		legend.position = "top",  
		legend.direction = "horizontal") +
	xlab("X coordinate") + ylab("Y coordinate") +
	guides(color = guide_legend(title = "Channel"))

ggsave(paste0("Color_channels_underlying_tileCoords.", 
	str_replace_all(Sys.time(), " ", "_"), ".png"), p, width = 12, height = 12)

#--------------------------------------------
#Distributions of numbers of paths per field coordinate estimate
#--------------------------------------------
maxPaths <- params$mosaickingParams$maxPaths
toPlot <- subset(dat, (nPaths != maxPaths) & (fieldID != params$mosaickingParams$rootForStitching))
p <- ggplot(toPlot, aes(x = nPaths)) +
        geom_histogram(bins = maxPaths + 1) +
        facet_grid(well ~ cycle) +
        theme_classic() +
	xlab("Number of paths")

ggsave(paste0("Distribution_of_nPaths_", str_replace_all(Sys.time(), " ", "_"), ".png"), p,
        width = 4, height = 0.7 * length(unique(toPlot$well)),
        limitsize = FALSE)

#--------------------------------------------
#Displacement patterns per cycle
#--------------------------------------------
ngbrsDir <- params %$%
	file.path(params$analysisDir, params$mosaicDir,
		params$mosaickingParams$mosaicIntermediateFilesDir)
allCycFiles <- ngbrsDir %>%
	list.files() %>%
	subset(., endsWith(., ".coordEsts_est2_filtered.intmdt.txt")) %>%
	data.frame(fileName = ., well = substr(., 1, 3))
datNgbrs <- map(wells, function(w) {
        print(w)
	subset(allCycFiles, well == w) %$%
        map(fileName, function(x) {
                file.path(ngbrsDir, x) %>%
                        read.table(., sep = "\t", header = TRUE) %>%
                        select(FieldID, locEstYWrtNghbrToNorth, locEstXWrtNghbrToNorth,
                                locEstYWrtNghbrToWest, locEstXWrtNghbrToWest) %>%
                        filter_all(any_vars(!is.na(.))) %>%
                        mutate(cycle = str_remove_all(x, ".coordEsts_est2_filtered.intmdt.txt") %>%
                                        strsplit(., split = "-") %>%
                                        unlist() %>% magrittr::extract(., i = 3) %>% 
					str_replace('cyc', 'cycle'),
                                ch = str_remove_all(x, ".txt") %>%
                                        strsplit(., split = "-") %>%
                                        unlist() %>% magrittr::extract(., i = 2),
                                well = w)
        }) %>% bind_rows()
}) %>% bind_rows()

datNgbrs$row <- datNgbrs$well %>% 
	map_int(., function(x) substr(x, 1, 1) %>% 
		equals(LETTERS, .) %>% which())
datNgbrs$col <- datNgbrs$well %>% substr(., 2, 3) %>% as.numeric()

datNgbrs %<>% merge(., chToStain, by = c("cycle", "ch"))

colnames(datNgbrs)[colnames(datNgbrs) == "FieldID"] <- "fieldID"

datNgbrs %<>% 
	merge(dat, ., 
		by = c("cycle", "fieldID", "well", "coordSource", "row", "col"), 
		all.x = TRUE)

p <- ggplot(datNgbrs, aes(x = locEstXWrtNghbrToNorth, y = locEstYWrtNghbrToNorth,
                shape = cycle, color = coordSource)) +
        geom_point() +
        facet_grid(row ~ col) +
        theme(axis.text.x = element_text(angle = 90),
                legend.position = "top",
                legend.direction = "horizontal") +
	guides(color = guide_legend(title = "Channel")) +
	xlab("X coordinate wrt northern neighbor") +
	ylab("Y coordinate wrt northern neighbor")

ggsave(paste0("camera_lateral_and_traversal_displacements_up_down_", 
	str_replace_all(Sys.time(), " ", "_"), ".png"), 
	p, width = 12, height = 12,
        limitsize = FALSE)

p <- ggplot(datNgbrs, aes(x = locEstXWrtNghbrToWest, y = locEstYWrtNghbrToWest,
                shape = cycle, color = coordSource)) +
        geom_point() +
        facet_grid(row ~ col) +
        theme(axis.text.x = element_text(angle = 90),
                legend.position = "top",
                legend.direction = "horizontal") +
        guides(color = guide_legend(title = "Channel")) +
	xlab("X coordinate wrt western neighbor") +
	ylab("Y coordinate wrt western neighbor")

ggsave(paste0("camera_lateral_and_traversal_displacements_left_right_", 
	str_replace_all(Sys.time(), " ", "_"), ".png"), p, width = 12, height = 12, 
	limitsize = FALSE)

