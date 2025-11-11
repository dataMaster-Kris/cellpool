library(tidyverse)
library(magrittr)

cp <- read.table("../60.barcodes/integrated_intensity_table.txt.gz",
                 header = TRUE, sep = "\t")
wells <- substr(cp$objId, 1, 3) %>% unique()

cpFtrProfiles <- read.table("../60.barcodes/integrated_intensity_table.txt.gz",
                            header = TRUE, sep = "\t")
cpBcdCalls <- read.table("../60.barcodes/barcode_classifications.txt.gz",
                         header = TRUE, sep = "\t")

elemGrps <- data.frame(cp = c("cycle1.ch1", "cycle2.ch1", "cycle2.ch2", "cycle2.ch3",
                              "cycle3.ch1", "cycle3.ch2", "cycle4.ch2", "cycle4.ch3"),
                       elem = c("EGFP", "OLLAS", "V5", "HA", "ALFA", "FLAG", "MYC", "VSVG"))
elemGrps$cycle <- substr(elemGrps$cp, 1, 6)

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

cp$well <- substr(cp$objId, 1, 3)
elm1 <- "VSVG"
elm2 <- "FLAG"
elm3 <- c("ALFA", "OLLAS")
thisWell <- subset(tcMap, grepl(elm1, bcdElems)) %>%
	subset(., grepl(elm2, bcdElems)) %>%
	subset(., grepl(elm3[1], bcdElems) | grepl(elm3[2], bcdElems))
thisDat <- subset(cp, well %in% thisWell$well) %>% 
	merge(., tcMap[, c("well", "bcdElems")])

thisElms <- subset(elemGrps, elem %in% c(elm1, elm2))
thisDat <- thisDat[, c(thisElms$cp, "bcdElems")] %>%
	set_colnames(c(thisElms$elem, "bcd"))
write.table(thisDat, file = "FLAG_vs_VSVG_in_two_barcodes.txt", 
	    sep = "\t", row.names = FALSE, quote = FALSE)

p <- ggplot(thisDat, aes(x = get(elm1) %>% log2(),
			 y = get(elm2) %>% log2(), color = bcd)) +
	geom_point(alpha = 0.1) +
	coord_cartesian(xlim = c(7, 16), ylim = c(7, 16)) +
	xlab(elm1) + ylab(elm2) +
	theme_classic()
ggsave(paste0("./", 
	      elm1, "_vs_", elm2, ".png"))






