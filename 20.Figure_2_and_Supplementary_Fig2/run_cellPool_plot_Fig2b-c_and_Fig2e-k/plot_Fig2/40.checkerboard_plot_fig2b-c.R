library(tidyverse)
library(magrittr)
library(patchwork)
library(ggh4x)

bcdCalls <- read.table("../60.barcodes/barcode_classifications.txt.gz", header = TRUE)
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
bcdCalls$well <- substr(bcdCalls$objId, 1, 3)

bcdCalls <- merge(bcdCalls, tcMap, by = "well", suffixes = c(".predicted", ".truth"))

#For every well, summarize how many cells with barcode profiles are in the well, how many with GFP+, how many with GFP-, for each epitope, do the same. Also, add truth status of the well.
bcdScr <- subset(bcdCalls, (barcodeSep > 0.99))
toPlot <- data.frame(epitope = c(), well = c(), pos = c(), truth = c(), total = c())
for (ep in elemGrps$elem) {
	for (wl in tcMap$well) {
		print(c(ep, wl))
		thisDat <- subset(bcdScr, grepl(ep, bcdElems.predicted) & 
					(well == wl))
		thisInfo <- data.frame(epitope = ep, 
				       well = wl,
				       pos = nrow(thisDat), 
				       truth = subset(bcdCalls, grepl(ep, bcdElems.truth) & 
						      		(well == wl) & 
								(cycle1.ch1 == 1)) %>% 
				       				nrow(),
				       total = subset(bcdCalls, well == wl) %>% nrow())
		if (ep == "EGFP") {
			thisInfo <- data.frame(epitope = ep,
                                       well = wl,
				       pos = subset(bcdScr, (cycle1.ch1 == 1) &  (well == wl)) %>%
					       nrow(),
				       truth = subset(bcdCalls, well == wl) %>% nrow(),
				       total = subset(bcdCalls, well == wl) %>% nrow()) 
		}
		toPlot %<>% bind_rows(., thisInfo)
	}
}

toPlot$fracPosInPosWells <- toPlot$pos/toPlot$truth
toPlot$fracPosInNegWells <- toPlot$pos/toPlot$total
toPlot$fracPos <- toPlot$fracPosInPosWells
toPlot$fracPos[toPlot$truth == 0] <- toPlot$fracPosInNegWells[toPlot$truth == 0]
toPlot$truthState <- toPlot$truth > 0
toPlot$truthState[(toPlot$well %in% subset(tcMap, treatment == "no_virus")$well) &
		  (toPlot$epitope == "EGFP")] <- FALSE
toPlot %<>% merge(., elemGrps, by.x = "epitope", by.y = "elem")
toPlot$row <- toPlot %$% substr(well, 1, 1)
toPlot$col <- toPlot %$% substr(well, 2, 3) %>% as.numeric()

toPlot$row %<>% factor(., levels = rev(LETTERS[1:8]))
toPlot %<>% subset(., total > 10)
write.table(toPlot, "dat_for_checkerboard_plot_05122025.txt", row.names = FALSE, 
	    quote = FALSE, sep = "\t")

p1 <- ggplot(toPlot, aes(x = col, y = row, fill = fracPos)) +
	geom_tile() + 
	facet_nested_wrap(vars(cycle, epitope), nrow = 2) +
	theme_classic() +
	theme(legend.position = "top", legend.direction = "horizontal",
	      legend.title = element_blank(), legend.frame = element_rect(color = "black"), 
	      axis.text.y = element_blank(), axis.ticks.y = element_blank()) + #,
	      #text = element_text(family = "Arial", size = 5, color = "black")) + 
	ylab(NULL) + xlab(NULL) +
	scale_fill_gradient(low = "white", #"#f45f5a", 
			    high = "#18b2b7", breaks = c(0, 0.5, 1), limits = c(0, 1))

p2 <- ggplot(toPlot, aes(x = col, y = row, fill = truthState)) +
        geom_tile() +
        facet_nested_wrap(vars(cycle, epitope), nrow = 2) +
        theme_classic() + 
        theme(legend.position = "top", legend.direction = "horizontal",
              legend.title = element_blank(), legend.key = element_rect(color = "black")) + #,
	      #text = element_text(family = "Arial", size = 5, color = "black")) + 
	ylab(NULL) + xlab(NULL) +
	scale_fill_manual(values = c("white", #"#f45f5a", 
				     "#18b2b7"))

p <- p2 | p1
ggsave(paste0("checkerboard_plot_", Sys.time(), ".pdf"), width = 5.5*2, height = 4.4*2)









