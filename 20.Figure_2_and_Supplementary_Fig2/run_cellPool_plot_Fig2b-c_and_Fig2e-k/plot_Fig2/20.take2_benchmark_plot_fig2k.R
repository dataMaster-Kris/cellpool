library(tidyverse)
library(magrittr)

bcdScr <- read.table("../60.barcodes/barcode_classifications.txt.gz", header = TRUE, sep = "\t")
kudo <- read.table("Kudo_et_al_2025_05_05.txt", header = TRUE, sep = "\t")
zunder <- read.table("Zunder_et_al_2025_05_05.txt", header = TRUE, sep = "\t")

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

bcdScr$well <- substr(bcdScr$objId, 1, 3)
bcdScr <- merge(bcdScr, toMerge, all.x = TRUE)
bcdScr$bcdSep <- bcdScr$barcodeSep
kudo$bcdElems <- kudo$bcd
zunder$bcdElems <- zunder$bcdEpis
zunder$bcdGrndTruth <- zunder$bcdEpiTruths

barcodeSepCutoffs <- (1 - (0:100)*0.01) %>% c(0.9999, 0.999, 0.99 + 0.001*(1:8))
yield <- c()
error_rate <- c()
meth <- c()
ctO <- c()
for (ct in barcodeSepCutoffs) {
	for (nxDat in c("bcdScr", "kudo", "zunder")) {
		meth %<>% c(., nxDat)
		nxDat <- get(nxDat)
		bcdCalls <- subset(nxDat, bcdSep >= ct)
		yield %<>% c(., sum(nrow(bcdCalls)/nrow(nxDat)))
		error_rate %<>% 
			c(., sum(bcdCalls$bcdElems != bcdCalls$bcdGrndTruth, 
				 na.rm = TRUE)/nrow(bcdCalls))
		ctO %<>% c(., ct)
	}
}

toPlot <- data.frame(cutoff = ctO, yield = yield, fdr = error_rate, meth = meth)
toPlot %<>% merge(., data.frame(meth = c("bcdScr", "kudo", "zunder"), 
				methId = c("cellPool", "EPICode", "SCD")))
write.table(toPlot, file = "dat_cellpool_vs_epicode_vs_scd_benchmark.txt",
	    row.names = FALSE, quote = FALSE, sep = "\t")
p <- ggplot(toPlot, aes(x = fdr, y = yield, color = methId)) + 
	#geom_point() +
	theme_classic() +
	coord_cartesian(xlim = c(0, 0.1)) +
	geom_line() +
	scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
	theme(legend.direction = "horizontal", legend.position = "top",
	      legend.title = element_blank())
ggsave(paste0("pr_", Sys.time(), ".pdf"), width = 3, height = 3)






