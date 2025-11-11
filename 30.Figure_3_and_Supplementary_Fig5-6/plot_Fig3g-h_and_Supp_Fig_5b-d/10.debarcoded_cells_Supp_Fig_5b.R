library(tidyverse)
library(magrittr)
library(rjson)

params <- fromJSON(file = "../00.config/config.json")
wells <- file.path(params$analysisDir, params$listOfWells) %>%
        readLines()
bkbnThr <- readLines("../60.barcodes/backbone_thresholds.txt") %>% as.numeric()
barcodeSepCutoff <- 0.95
barcodeSepCutoffCtrl <- 0.99
intDat <- read.table("../60.barcodes/integrated_intensity_table.txt.gz", header = TRUE, sep = "\t")
intDat$well <- substr(intDat$objId, 1, 3)

bcdCalls <- read.table("../60.barcodes/barcode_classifications.txt.gz",
                                header = TRUE, sep = "\t")
nGFPpos <- sum(bcdCalls$cycle1.ch2 == 1, na.rm = TRUE)
nMOIgt1 <- sum(apply(bcdCalls[, colnames(bcdCalls) %>% 
		startsWith("cycle")], 1, sum, na.rm = TRUE) > 4)
nNone <- sum(bcdCalls$bcdElems == "None", na.rm = TRUE)
highGFPThr <- subset(intDat, cycle1.ch2 > bkbnThr[2]) %$% quantile(cycle1.ch2, 0.5)
highGFP <- intDat$objId[intDat$cycle1.ch2 > highGFPThr]

barcodes <- read.table("../00.config/barcodes.txt", header = TRUE, sep = "\t")
libMap <- read.table("../00.config/libraryMap.txt", header = TRUE, sep = "\t")
tcMap <- read.table("../00.config/tcMap.txt", header = TRUE, sep = "\t")

barcodes$bcd <- apply(barcodes[, paste0("epi", 1:3)], 1, paste0, collapse = "_")
libMap <- merge(libMap, barcodes[, c("id", "bcd")], by.x = "bcdId", by.y = "id")
libMap <- merge(libMap, tcMap, by = "poolId")
ntcBarcodes <- subset(libMap, startsWith(sgId, "non-")) %$% unique(bcd)
bcdCalls$ntcBcd <- bcdCalls$bcdElems %in% ntcBarcodes
bcdCalls %<>% subset(., (bcdElems == "None") |
                ((barcodeSep >= barcodeSepCutoff) & !(ntcBcd)) |
                ((barcodeSep >= barcodeSepCutoffCtrl) & (ntcBcd) & (objId %in% highGFP)))

bcdCalls$well <- substr(bcdCalls$objId, 1, 3)
bcdCalls %<>% subset(., !is.na(bcdElems), select = c("objId", "bcdElems", "well")) %>%
        set_colnames(., c("objId", "bcd", "well"))
bcdCalls <- merge(bcdCalls, libMap[, c("bcd", "well", "sgId", "gene", "poolId")], all.x = TRUE)

bcdFreq <- table(bcdCalls$sgId) %>% data.frame()
p <- ggplot(bcdFreq, aes(log10(Freq))) +
        geom_histogram(bins = 30) +
        theme_classic() +
        xlab(NULL) + ylab(NULL) +
        theme(axis.text = element_text(size = 5),
                plot.margin = margin(0,0,0,0)) +
	coord_cartesian(xlim = c(1, 4)) +
	scale_x_continuous(breaks = c(1:4)) +
        geom_vline(xintercept = median(log10(bcdFreq$Freq)),
                        linetype = "dashed", color = "red")

ggsave("Counts_per_guide_histogram.pdf", p, width = 1.2, height = 1.85)

save(bcdCalls, file = "barcodes_to_cells.rds")

