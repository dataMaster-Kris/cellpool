library(magrittr)
library(tidyverse)

mosaicDir <- "../30.mosaics"
#--------------------------------------------
#Heatmap of SEMs for all tile coordinates in all cycles
#--------------------------------------------
filesToRead <- list.files(mosaicDir) %>%
	subset(., endsWith(., ".tileCentroidCoords.txt")) %>%
	data.frame(fileName = ., well = substr(., 1, 3))	
wells <- unique(filesToRead$well)
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
							30)) %>% 
							min(.)) %>% 
					round(),
				gridY = divide_by(medY, abs(medY) %>%
					subset(., is_greater_than(., 
							30)) %>%
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

nGridLenInWell <- 11 
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
	dplyr::select(cycle, seX, seY, gridX, gridY, row, col) %>% 
	gather(imgAxis, SEM, -cycle, -gridX, -gridY, -row, -col)
toPlot$imgAxis %<>% str_remove_all("se")
nGridLenInWell <- toPlot %$% 
	subset(gridX, !is.nan(gridX)) %>% #Assuming equal radii along the X and Y axes
	unique() %>% 
	length()

set.seed(34028766)
smpldWells <- "A07" 
p <- ggplot(toPlot %>%
         mutate(well = 
                  paste0(magrittr::extract(LETTERS, row), 
                         str_pad(col, pad = "0", width = 2)
                    )
                ) %>%
         subset(., well %in% smpldWells), 
       aes(x = col + (1/nGridLenInWell)*gridX, 
		y = #-row + 
		  (1/nGridLenInWell)*gridY, fill = SEM)) + #row numbers increase downwards
	theme_classic() + 
	geom_tile() + 
	scale_fill_continuous(low = "white", high = "steelblue") +
	facet_grid(cycle ~ imgAxis) +
	theme(axis.line = element_blank(),
		axis.text = element_blank(), 
		axis.ticks = element_blank(),
		axis.title = element_blank(),
		legend.position = "top", 
		legend.direction = "horizontal")

ggsave("Supp_fig_2a.pdf", width = 2.7, height = 4.89)
