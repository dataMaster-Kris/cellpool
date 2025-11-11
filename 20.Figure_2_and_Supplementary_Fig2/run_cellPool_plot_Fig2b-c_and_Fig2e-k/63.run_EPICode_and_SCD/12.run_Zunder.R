library(tidyverse)
library(magrittr)

bkbnThr <- readLines("../60.barcodes/backbone_thresholds.txt") %>% as.numeric()

dat <- read.table("../60.barcodes/integrated_intensity_table.txt.gz", header = TRUE)
calls <- read.table("../60.barcodes/barcode_classifications.txt.gz", header = TRUE)
tcMap <- read.table("../00.config/tcMap.txt", header = TRUE)

gfpPos <- dat$objId[dat$`cycle1.ch2` > bkbnThr[2]]
outDir <- "tables"

bcds <- read.table("../00.config/barcodes.txt", header = TRUE)
chToElem <- read.table("../00.config/chToElemId.txt", header = TRUE, sep = "\t")
chToStain <- read.table("../00.config/chToStainId.txt", header = TRUE, sep = "\t")

bcds <- bcds[, c("epi1", "epi2", "epi3")]
uniqEpis <- unique(bcds %>% unlist())
epiToCh <- map(uniqEpis, 
	       function(x) {
	       		intmd <- which(chToElem == x, arr.ind = TRUE)
	       		paste0("cycle", intmd[1, 1], ".ch", intmd[1, 2])
	       }) %>% 
	set_names(uniqEpis)
chToEpi <- unlist(epiToCh) %>% names() %>% set_names(epiToCh)

#---------------------------------------
#Call epitrain next
scores <- dat

#Zunder et al. use probs = 1 in the following but it led to only ~10 cells debarcoded.
#Assuming that this is due to outliers, using probs = 0.99 in the following instead.
scores[, setdiff(names(dat), "objId")] <- apply(dat[, setdiff(names(dat), "objId")], 2,
                function(x) x/quantile(x, probs = 0.99, na.rm = TRUE))

chs <- setdiff(names(dat), c("cycle1.ch1", "cycle1.ch2", "cycle2.ch4", "cycle3.ch4", "cycle3.ch3",
			     "cycle4.ch1", "cycle4.ch4", "objId"))
bcdCalls <- apply(scores[, chs], 1, function(x) {
			if (any(is.na(x))) return(c(NA, NA))
			c(chToEpi[unlist(chs[order(x) %>% tail(., n = 3)])], 
			  magrittr::extract(x[order(x)], i = 4:5) %>%
				  diff() %>% set_names(NULL))
		})

bcdCalls %<>% t() %>% as.data.frame() %>% 
	set_colnames(c(paste0("epi", 1:3), "bcdSep")) %>% 
	mutate(bcdSep = as.numeric(bcdSep))
bcdCalls$bcdEpis <- apply(bcdCalls[, paste0("epi", 1:3)], 1, 
			  function(x) sort(x) %>% paste0(., collapse = "_"))

bcdCalls$objId <- dat$objId

tcMap <- read.table("../00.config/tcMap.txt", header = TRUE)
barcodes <- read.table("../00.config/barcodes.txt", header = TRUE)
wells <- read.table("../00.config/wells.txt")
tcMap %<>% subset(., well %in% wells$V1)
bcdId <- data.frame(treatmentId = c("X2", "X3", "X4", "X6", "X7", "X9",
                                    "X10", "X11", "X12", "X13", "X15",
                                    "X16", "X17", "X18", "X19", "X20", "A1", "C1", "E2"),
                    id = 1:19)
tcMap$treatmentId <- strsplit(tcMap$treatment, split = "_") %>%
        plyr::laply(., magrittr::extract, i = 2)

tcMap <- merge(tcMap, bcdId, all.x = TRUE)

barcodes$bcdElems <- apply(barcodes[, paste0("epi", 1:3)], 1, paste0, collapse = "_")
tcMap <- merge(tcMap, barcodes[, c("id", "bcdElems")], all.x = TRUE)

toMerge <- tcMap[, c("well", "bcdElems")] %>%
        set_colnames(., c("well", "bcdGrndTruth"))

toMerge$bcdEpiTruths <- toMerge$bcdGrndTruth %>% 
	map_chr(., function(x) strsplit(x, split = "_") %>% unlist() %>% sort() %>% paste0(., collapse = "_"))

outDat <- bcdCalls
outDat$well <- substr(outDat$objId, 1, 3)
outDat %<>% merge(., toMerge, all.x = TRUE)

write.table(outDat, "Zunder_et_al_2025_05_05.txt", row.names = FALSE, quote = FALSE, sep = "\t")

barcodeSepCutoffs <- 1 - (0:100)*0.01
yield <- c()
error_rate <- c()
for (ct in barcodeSepCutoffs) {
        bcdCalls <- subset(outDat, bcdSep >= ct)
        yield %<>% c(., sum(nrow(bcdCalls)/nrow(outDat)))
        error_rate %<>% c(., sum(bcdCalls$bcdEpis != bcdCalls$bcdEpiTruths, na.rm = TRUE)/nrow(bcdCalls))
}

toPlot <- data.frame(cutoff = barcodeSepCutoffs, yield = yield, fdr = error_rate)

p <- ggplot(toPlot, aes(x = fdr, y = yield)) + geom_point()
ggsave("pr_take1_zunder.png")


p <- ggplot(toPlotDf, aes(x = 1-fdr, y = yield, shape = meth)) + geom_point() + 
	theme_classic()
ggsave("pr_take1_kudo_cellpool_zunder.png")






