library(magrittr)
library(tidyverse)
library(rjson)

args <- commandArgs(trailingOnly=TRUE)
params <- fromJSON(file = args[1])
iter <- as.numeric(params$segmentationParams$segIter)
set.seed(params$unpoolSegQcParams$rndSeed)

ftrDir <- file.path(params$analysisDir, params$segmentationDir, 
	paste0("iter", iter, params$segmentationParams$featureDirSuffix))
wells <- file.path(params$analysisDir, params$listOfWells) %>%
	readLines()

dat <- map(wells, function(well) {
	map(params$cycles, function(cyc) {
		file.path(ftrDir, paste0(well, "-cyc", cyc, ".", 
				str_replace(params$segmentationParams$cellposeParams$useStain, 
				" ", "_"), ".txt")) %>% 
			read.table(., sep = "\t", header = TRUE) %>% 
			mutate(label = paste0(well, "n", label), cycle = cyc)
	}) %>% 
		bind_rows() %>% 
		mutate(well = well)
}) %>% 
	bind_rows()

#--------------------------
#Subsample in different bins
#--------------------------
nBins <- params$unpoolSegQcParams$nFtrBins
nPerClass <- params$unpoolSegQcParams$nSamplesPerBin

ftrs4QC <- params$unpoolSegQcParams$ftrs4Qc

bin_cutoffs <- list()
for (ftr in ftrs4QC) {
	bin_cutoffs[[ftr]] <- log10(dat[, ftr]) %>% 
		quantile(., (0:nBins)/nBins)
	bin_cutoffs[[ftr]] <- data.frame(pctlClass = (0:(nBins - 1) + 0.5)*100/nBins,
			rangeStart = bin_cutoffs[[ftr]][-nBins - 1], 
			rangeEnd = bin_cutoffs[[ftr]][-1]) %>% 
		set_rownames(NULL)

	dat[, paste0(ftr, "PctlClass")] <- NA
	for (i in 1:nBins) {
		dat[(dat %$% get(ftr) %>% 
			between(., 10^bin_cutoffs[[ftr]]$rangeStart[i], 
				10^bin_cutoffs[[ftr]]$rangeEnd[i])) | 
			(dat %$% get(ftr) %>% 
				equals(., round(10^bin_cutoffs[[ftr]]$rangeEnd[i]))) | 
			((i == 1) & (dat[, ftr] == round(10^bin_cutoffs[[ftr]]$rangeStart[i]))),
			paste0(ftr, "PctlClass")] <- 
				paste0(bin_cutoffs[[ftr]]$pctlClass[i], "_", 
					round(10^bin_cutoffs[[ftr]]$rangeStart[i]), "_",
					round(10^bin_cutoffs[[ftr]]$rangeEnd[i]))
	}
}

nominations <- dat %>%
        group_by(areaPctlClass) %>%
        sample_n(nPerClass) %>%
        arrange(label) %>%
	mutate(objId = paste0(label, "_", cycle))

write.table(nominations, 
        file = "nominations.txt",
        sep = "\t",
        row.names = FALSE, quote = FALSE)

map(unique(nominations$well), function(x) {
	map(unique(nominations$cycle), function(cyc) {
        	write.table(subset(nominations, (well == x) & 
					(cycle == cyc), select = "label", drop = FALSE), 
                	file = paste0(x, "-cyc", cyc, ".nominations.txt"),
                	row.names = FALSE, quote = FALSE, col.names = FALSE)
	}) %>% unlist()
}) %>% unlist()



















