library(tidyverse)
library(magrittr)
library(rjson)

args <- commandArgs(trailingOnly=TRUE)
params <- fromJSON(file = args[1])

emTblDir <- file.path(params$analysisDir, params$barcodeAnalysisDir)
cuts <- 1:params$debarcodingParams$nCuts
barcodes <- file.path(emTblDir, "barcodes.refmt.txt") %>%
        read.table(., sep = "\t", header = TRUE) %>%
        mutate(bcdIx = 0:(n() - 1))

dat <- file.path(emTblDir, "integrated_intensity_table.txt.gz") %>% 
	read.table(., sep = "\t", header = TRUE)
chToElem <- file.path(emTblDir, "chToElem.txt") %>%
	read.table(., sep = "\t", header = TRUE)
calls <- file.path(emTblDir, "calls_with_pop_level_epitope_gating.txt.gz") %>%
	read.table(., sep = "\t", header = TRUE)

nCuts <- params$debarcodingParams$nCuts
maxLogLDat <- file.path(emTblDir, "init_and_iter_with_max_logL.txt") %>%
	read.table(., sep = "\t", header = TRUE)
scores <- map(1:nCuts, function(x) {
	paste0("cut_", x, "_init_", 
			maxLogLDat$initMax[maxLogLDat$cut == paste0("C", x)], "_iter_", 
			maxLogLDat$iterMax[maxLogLDat$cut == paste0("C", x)], "_r.txt") %>%
		file.path(emTblDir, .) %>%
		read.table(., sep = "\t", header = TRUE)
}) %>% bind_rows()
barcode_cols <- file.path(params$analysisDir, params$barcodeAnalysisDir, "barcodes.refmt.txt") %>% 
	read.table(., sep = "\t", header = TRUE)
barcode_cols$colNumber <- 1:nrow(barcode_cols)
barcode_cols$bcdElems <- apply(barcode_cols[, startsWith(colnames(barcode_cols), "epitope")], 
			       1, paste0, collapse = "_")

#-----------------------------------------
newCalls <- merge(calls[, "objId", drop = FALSE], scores)
newCalls[, c("calls", "barcodeSep")] <- 
	apply(newCalls[, !grepl("objId", colnames(newCalls))], 1, 
		function(x) {
			orderX <- order(as.numeric(x), decreasing = TRUE)
			xMax <- x[orderX[1]]
			x2nd <- x[orderX[2]]
			return(c(orderX[1], xMax - x2nd))
}) %>% t()

newCalls$bcdElems <- barcode_cols$bcdElems[newCalls$calls]
newCalls <- merge(newCalls, calls[, "objId", drop = FALSE], sort = FALSE, all.y = TRUE)

calls <- merge(calls, newCalls[, c("objId", "bcdElems", "barcodeSep")], 
		all.x = TRUE, sort = FALSE)
calls[!is.na(calls$bcdElems), 
      chToElem$channel[chToElem$element != params$debarcodingParams$backboneProt]] <- 
	map_dfr(calls$bcdElems[!is.na(calls$bcdElems)], function(y) {
		data.frame(rbind(map_int(setdiff(chToElem$element, 
						 params$debarcodingParams$backboneProt), 
					 grepl, x = y)))
	})

calls[which(calls[, chToElem$channel[chToElem$element == 
	params$debarcodingParams$backboneProt]] == 0), "bcdElems"] <- "None"
gz1 <- gzfile("barcode_classifications.txt.gz", "w")
write.table(calls, gz1, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz1)

#Show yield for different barcodeSep cutoffs
cutoffs <- seq(75, 99, 1)/100
yields <- map_dbl(cutoffs, function(x) sum(calls$barcodeSep >= x, na.rm = TRUE))
nPosEpis <- apply(calls[, chToElem$channel], 1, sum, na.rm = TRUE) %>% 
	table() %>% 
	data.frame() %>% 
	set_colnames(., c("numberOfEpitopesFoundInCellPlusBackboneProt", "Freq"))
yieldPctWrtFreqOneInfctn <- 
	yields/nPosEpis$Freq[
		nPosEpis$numberOfEpitopesFoundInCellPlusBackboneProt == 
			params$debarcodingParams$nEpiPerCombo + 1]
yields <- data.frame(barcodeSepCutoff = cutoffs, 
			countDebarcoded = yields,
			pctDebarcodedWrtTotalSinglyInfected = yieldPctWrtFreqOneInfctn)

write.table(yields, "Yield_vs_barcodeSepCutoff.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(nPosEpis, "nCells_vs_nEpitopes.txt", sep = "\t", row.names = FALSE, quote = FALSE)



