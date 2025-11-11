library(tidyverse)
library(magrittr)
library(rjson)

args <- commandArgs(trailingOnly = TRUE)
params <- fromJSON(file = args[1])

mosaicDir <- file.path(params$analysisDir, params$mosaicDir) 
wells <- file.path(params$analysisDir, params$listOfWells) %>%
        readLines()

dat <- map(wells, function(x) {
	file.path(mosaicDir, paste0(x, ".rawIntensityQuantiles.txt")) %>%
		read.table(., sep = "\t", header = TRUE)
	}) %>%
	bind_rows() %>%
        group_by(cycle, ch) %>%
        group_map(~ apply(.x, 2, max) %>% 
		data.frame() %>% 
		t() %>% 
		cbind(., .y) %>% 
		data.frame() %>% 
		set_rownames(NULL)) %>% 
	bind_rows()

write.table(dat, 
	"wholePlate.rawIntensityQntls.txt",
	sep = "\t", row.names = FALSE, quote = FALSE)

