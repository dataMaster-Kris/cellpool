library(tidyverse)
library(magrittr)
library(rjson)

args <- commandArgs(trailingOnly = TRUE)
params <- fromJSON(file = args[1])

classLabels <- file.path(params$analysisDir, params$unpoolClustParams$listOfClasses) %>%
	read.table(., header = TRUE, sep = "\t")
tcMap <- file.path(params$analysisDir, params$tcMapFile) %>%
        read.table(., header = TRUE, sep = "\t")
libMap <- file.path(params$analysisDir, params$libraryMapFile) %>%
        read.table(., header = TRUE, sep = "\t")
ftrs <- file.path(params$analysisDir, params$barcodeAnalysisDir, "feature_profiles.txt.gz") %>%
        read.table(., header = TRUE, sep = "\t")
clsTbl <- file.path(params$analysisDir, params$barcodeAnalysisDir, 
		"barcode_classifications.txt.gz") %>%
	read.table(., header = TRUE, sep = "\t") %>% 
	subset(., is.na(barcodeSep) | (barcodeSep >= params$debarcodingParams$barcodeSep)) %>%
	select(objId, bcdElems)
barcodes <- file.path(params$analysisDir, params$barcodes) %>% 
	read.table(., header = TRUE, sep = "\t")
barcodes$bcdElems <- barcodes[, colnames(barcodes) %>% startsWith("epi")] %>% 
	apply(., 1, function(x) paste0(x, collapse = "_"))
nPerFile <- params$unpoolClustParams$nPerFile %>% prod()

dat <- merge(clsTbl, ftrs, by.x = "objId", by.y = "label")
dat$well <- substr(dat$objId, 1, 3)
dat %<>% merge(., tcMap, by = "well", all.x = TRUE)
dat %<>% merge(., barcodes[, c("id", "bcdElems")] %>% 
		set_colnames(., c("bcdId", "bcdElems")), by = "bcdElems", all.x = TRUE)
dat %<>% merge(., libMap, by = c("poolId", "bcdId"), all.x = TRUE)
dat$bcdId[dat$bcdElems == "None"] <- "None"
dat$sgId[dat$bcdElems == "None"] <- "None"
dat$gene[dat$bcdElems == "None"] <- "None"

for (ix in 1:nrow(classLabels)) {
	thisGrp <- classLabels[ix, ]
	thisAll <- dat
	thisGrpId <- thisGrp$groupId
	thisGrp <- thisGrp %>% select(-groupId)
	for (cTag in colnames(thisGrp)) {
		if (is.na(thisGrp[, cTag])) next
		thisAll %<>% subset(., get(cTag) == thisGrp[, cTag])
	}
		
	thisAll <- thisAll[order(thisAll[, params$unpoolClustParams$sortByFeature],
		decreasing = params$unpoolClustParams$orderDecreasing), ]
	for (j in 1:ceiling(nrow(thisAll)/nPerFile)) {
		write.table(thisAll[(nPerFile * j - nPerFile + 1):min(c(nPerFile * j, 
				nrow(thisAll))), 
				c("objId", params$unpoolClustParams$sortByFeature), 
				drop = FALSE], 
			file = paste0(thisGrpId, "_", j, ".toUnpool.txt"), sep = "\t", 
			col.names = FALSE, row.names = FALSE, quote = FALSE)
	}
}

