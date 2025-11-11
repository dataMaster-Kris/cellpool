library(tidyverse)
library(magrittr)
library(rjson)

args <- commandArgs(trailingOnly = TRUE)
params <- fromJSON(file = args[1])

notes <- c()
dat <- read.table("integrated_intensity_table.txt.gz", header = TRUE, sep = "\t")
chToElem <- read.table("chToElem.txt", sep = "\t", header = TRUE)
tcMap <- file.path(params$analysisDir, params$tcMapFile) %>%
        read.table(., header = TRUE, sep = "\t")

calls <- dat[, c("objId", chToElem$channel)]
calls[, names(calls) %>% subset(., startsWith(., "cycle"))] <- NA

bkbnThresh <- readLines("backbone_thresholds.txt") %>% 
	as.numeric()
bkbnCh <- chToElem$channel[chToElem$element == params$debarcodingParams$backboneProt] %>%
	str_replace(., "-", ".")
 
calls[, bkbnCh] <- ifelse(dat[, bkbnCh] > bkbnThresh[2], 1, NA)

def_neg <- dat[, bkbnCh] < bkbnThresh[1]
calls[which(def_neg), chToElem$channel] <- 0

p <- ggplot(dat, aes(get(bkbnCh) %>% log2())) +
	geom_histogram(bins = 100) +
        coord_cartesian(xlim = c(0, params$bitDepth)) +
        xlab(paste0("log2(", params$debarcodingParams$backboneProt, " intensity)")) +
        ylab("Count") +
	geom_vline(xintercept = log2(bkbnThresh), color = "red", linetype = "dotdash")

ggsave(paste0("Population_level_", params$debarcodingParams$backboneProt, "_intensity_histograms_",
        str_replace_all(Sys.time(), " ", "_"), ".png"),
        height = 5, width = 7, p)

notes %<>% c(., "Performed backbone protein thresholding.", 
		paste0("Number barcode positive: ", sum(calls[, bkbnCh], na.rm = TRUE)),
		paste0("Number barcode negative: ", sum(calls[, bkbnCh] == 0, na.rm = TRUE)),
		paste0("Number barcode NA: ", sum(is.na(calls[, bkbnCh]))),
		paste0("Total objects: ", nrow(calls)), "\n\n\n")
bkbnPos <- subset(calls, get(bkbnCh) == 1) %$%
        objId
densityRatioThresh <- params$debarcodingParams$thresholdDensityRatio
minRunLength <- params$debarcodingParams$thresholdRunLengthMin

#----------------------------------------------------
#Thresholding function to compare signal-background mixture with the background distribution
#----------------------------------------------------
findThreshold <- function(mixture, signalEnriched, background, densityRatioThresh, minRunLength,
		scaleDensities = FALSE) {
	mixHist <- hist(c(mixture, background, signalEnriched), nclass = 100, plot = FALSE) %$%
                data.frame(xstart = breaks[-length(breaks)],
                                xend = breaks[-1], density = density)
	
	thisBreaks <- mixHist %$%
                c(xstart, xend) %>%
                unique() %>%
                sort()
	
	sigEnrichedHist <- hist(signalEnriched, breaks = thisBreaks, plot = FALSE) %$%
                data.frame(xstart = breaks[-length(breaks)],
                                xend = breaks[-1], densitySigEnriched = density)	

	bckHist <- hist(background, breaks = thisBreaks, plot = FALSE) %$%
                data.frame(xstart = breaks[-length(breaks)],
                                xend = breaks[-1], densityBckgd = density)
	
	this <- merge(sigEnrichedHist, bckHist,
                        by = c("xstart", "xend"), all = TRUE)
        this$densitySigEnriched[this$densitySigEnriched %>% is.na()] = 0
        this$densityBckgd[this$densityBckgd %>% is.na()] = 0
        this %<>%
                mutate(densityRatio = densitySigEnriched/densityBckgd) %>%
                arrange(xstart)

	if (scaleDensities) this$densityRatio %<>% 
			multiply_by(., max(this$densityBckgd)/max(this$densitySigEnriched))

	thisRle <- this %$%
                is_weakly_greater_than(densityRatio, densityRatioThresh) %>%
                rle()
        firstGoodRun <- which(thisRle$values & (thisRle$lengths >= minRunLength)) %>%
                min()
        threshRow <- sum(thisRle$lengths[1:(firstGoodRun - 1)]) + 1

        this[threshRow, "xstart"]

}

#----------------------------------------------------
#Threshold the signal-background mixture
#----------------------------------------------------
chToElem$channel %<>% str_replace(., "-", ".")
epiThresh <- map(chToElem$channel %>% setdiff(., bkbnCh), function(x) {
	print(x)
	mixture <- log2(dat[, x])

	signalEnriched <- subset(dat, (objId %in% bkbnPos)) %$%
		get(x) %>%
		log2()

	background <- subset(dat, !(objId %in% bkbnPos) & !is.na(get(bkbnCh))) %$%
		get(x) %>%
		log2()

	data.frame(channel = x, threshold = findThreshold(mixture, signalEnriched, background, 
			densityRatioThresh, minRunLength, scaleDensities = FALSE))
}) %>% 
	bind_rows() 

intDat <- dat
intDat$well <- substr(intDat$objId, 1, 3)
intDat <- merge(intDat, tcMap[, c("well", "poolId")])

toPlot <- (intDat[, chToElem$channel %>% c(., "poolId")]) %>%
        set_colnames(paste0(chToElem$element, "\n", chToElem$fluorophore) %>% c(., "poolId"))

pdf(paste0("PoolId_level_intensity_histograms_with_thresholds_",
        str_replace_all(Sys.time(), " ", "_"), ".pdf"),
        height = 5, width = 7)
for (elem in setdiff(colnames(toPlot), "poolId")) {
	thisCh <- chToElem %>%
		subset(., paste0(element, "\n", fluorophore) == elem)
	thisThresh <- subset(epiThresh, channel == thisCh$channel)$threshold
        p <- ggplot(toPlot, aes(get(elem) %>% log2())) +
                geom_histogram(bins = 100) +
                facet_grid(poolId ~ .) +
                coord_cartesian(xlim = c(0, params$bitDepth)) +
                xlab("log2(intensity)") +
                ylab("Count") +
		geom_vline(xintercept = thisThresh, color = "red", linetype = "dotdash") +
                ggtitle(elem)
        print(p)
}
dev.off()

#----------------------------------------------------
epiGates <- list()
sepUp <- file.path(params$analysisDir, params$debarcodingParams$epiThreshUp) %>% 
	read.table(., header = FALSE, sep = "\t")
for (i in chToElem$channel) {
	if (i == bkbnCh) next
	thisSepUp <- sepUp$V2[sepUp$V1 == chToElem$element[chToElem$channel == i]]
	posGate <- epiThresh[(epiThresh$channel == i), "threshold"] + thisSepUp
	found_pos <- subset(dat, get(i) > 2^posGate, 
			select = "objId", drop = TRUE)
	calls[calls$objId %in% found_pos, i] <- 1

	negGate <- subset(dat, calls[, bkbnCh] == 0) %$%
			quantile(get(i) %>% log2(), 
				params$debarcodingParams$bckgdQntlCutoff, na.rm = TRUE)
	if (negGate > posGate - params$debarcodingParams$minSigVsBckgdGateLog2Sep) 
		negGate <- posGate - params$debarcodingParams$minSigVsBckgdGateLog2Sep
	found_neg <- subset(dat, get(i) <= 2^negGate, select = "objId", drop = TRUE)
	calls[calls$objId %in% found_neg, i] <- 0
	notes %<>% c(., paste0("Number ", chToElem$element[chToElem$channel == i], " positive: ",
				length(found_pos)),
			paste0("Number ", chToElem$element[chToElem$channel == i], " negative: ",
				length(found_neg), "\n\n"))
	epiGates[[i]] <- c(neg = as.numeric(negGate), pos = posGate)
}

nBkbnPos <- length(bkbnPos)
bkbnPos <- apply(calls[, chToElem$channel %>% setdiff(., bkbnCh)], 1,
        function(x) sum(x == 1, na.rm = TRUE) >= params$debarcodingParams$minEpiPosForBkbnPos)
calls[bkbnPos, bkbnCh] <- 1
notes %<>% c(., paste0("\n\nNumber of objects deemed positive because >", 
			params$debarcodingParams$minEpiPosForBkbnPos, " epitopes positive: ",
			sum(calls[, bkbnCh], na.rm = TRUE) - nBkbnPos), "\n\n\n")

#-----------------------------------------
#If assumeVeryLowMOI, resolve the ambiguities in cells with > nEpiPerCombo
#-----------------------------------------
if (params$debarcodingParams$assumeVeryLowMOI) {
	nEpiPerCombo <- params$debarcodingParams$nEpiPerCombo
	newCalled <- 0
	calls[, chToElem$channel] <- apply(calls[, chToElem$channel], 1, function(x) {
		out <- x
		if (sum(x[chToElem$channel != bkbnCh], na.rm = TRUE) >= nEpiPerCombo) {
			out[which(out != 1)] <- 0
			out[is.na(out)] <- 0
			newCalled <<- newCalled + 1
		}
		out
	}) %>%
		t()

	notes %<>% c(., paste0("Number of objects with ambiguities resolved under assumption of very low MOI: ", newCalled), "\n\n\n")
}

gz1 <- gzfile("calls_with_pop_level_epitope_gating.txt.gz", "w")
write.table(calls, gz1,
	sep = "\t", row.names = FALSE, quote = FALSE)
close(gz1)
write.table(epiGates %>% data.frame(), "epitopes_population_based_gates.txt",
	sep = "\t", row.names = FALSE, quote = FALSE)

#---------------------------------------------------
#Save notes if verbose is true
#---------------------------------------------------
if (params$debarcodingParams$verbose) {
	fileConn <- file(params$debarcodingParams$saveNotes)
	writeLines(notes, fileConn)
	close(fileConn)
}

#---------------------------------------------------
#Plot epitopes vs backbone protein signal
#---------------------------------------------------
for (ch in chToElem$channel[chToElem$element != params$debarcodingParams$backboneProt]) {
	print(ch)
        toPlot <- data.frame(element = dat[, ch],
                        backbone = dat[, bkbnCh], 
                        calls = as.factor(calls[, ch]))

	thisElement <- chToElem$element[chToElem$channel == ch]        

        p <- ggplot(toPlot, aes(x = log2(backbone), y = log2(element), color = calls)) +
                geom_point(alpha = 0.05) +
                geom_hline(yintercept = epiGates[[ch]],
                        color = "red", linetype = "dotdash") +
                geom_vline(xintercept = log2(bkbnThresh[1]), 
                                color = "red", linetype = "dotdash") +
                geom_vline(xintercept = log2(bkbnThresh[2]), 
                                color = "red", linetype = "dotdash") +
                scale_color_discrete() +
		xlab(paste0("log2(", params$debarcodingParams$backboneProt, ")")) +
		ylab(paste0("log2(", thisElement, ")")) +
		guides(color = guide_legend(override.aes = list(alpha = 1)))

        ggsave(paste0("scatter_backbone_vs_", thisElement, "_", 
		str_replace_all(Sys.time(), " ", "_"), ".png"))
}


