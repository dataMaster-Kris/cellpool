library(tidyverse)
library(magrittr)
library(ggh4x)

cp <- read.table("../../../60.barcodes/integrated_intensity_table.txt.gz",
		 header = TRUE, sep = "\t")
wells <- substr(cp$objId, 1, 3) %>% unique()
cpFtrProfiles <- read.table("../../../60.barcodes/integrated_intensity_table.txt.gz", 
			    header = TRUE, sep = "\t")
cpBcdCalls <- read.table("../../../60.barcodes/barcode_classifications.txt.gz",
			 header = TRUE, sep = "\t")

allFtrFiles <- list.files("../40.segmentation/iter1.step3.features/") %>%
	subset(., grepl("_maxShift_15_sigma_3.", .))
ashlr_ftrs <- map(wells, function(x) {
	 thisFiles <- subset(allFtrFiles, startsWith(allFtrFiles, x)) %>%
		 paste0("../40.segmentation/iter1.step3.features/", .)
	 thisDat <- map(thisFiles, 
			function(y) read.table(y, header = TRUE, sep = ",") %>%
				mutate(label = basename(y) %>% 
				       substr(., 1, 3) %>% 
				       paste0(., "n", label)) %>% 
				dplyr::select(label, area, orientation, mean, median) %>% 
				set_colnames(., c("label", 
						  paste0(c("area", "orientation", 
							   "mean", "median"), "_ch", 
							 basename(y) %>% 
								 str_remove_all(., ".csv") %>% 
								 strsplit(., split = ".ch") %>% 
								 unlist() %>% 
								 magrittr::extract(., i = 2))))) %>% 
		Reduce(merge, .)
	 }) %>% 
	bind_rows()

allTrackingFiles <- list.files("../45.track_cellpool_to_ashlar/") %>% 
	subset(., grepl(".cellpool_ashlar_tracking.stage2_euclid_dist.txt", .))
cp2AshlrTracking <-  map(wells, function(x) {
                 thisFile <- subset(allTrackingFiles, startsWith(allTrackingFiles, x)) %>%
                         paste0("../45.track_cellpool_to_ashlar/", .)
                 thisDat <- read.table(thisFile, header = TRUE, sep = "\t") %>%
                                        mutate(well = basename(thisFile) %>% substr(., 1, 3))
                 }) %>%
	bind_rows()
cp2AshlrTracking %<>% mutate(cellpool = paste0(well, "n", cellpool),
               ashlar = paste0(well, "n", ashlar))
v <- merge(cp, cp2AshlrTracking, by.x = "objId", by.y = "cellpool")  %>%
 	merge(ashlr_ftrs, ., by.x = "label", by.y = "ashlar")
elemGrps <- data.frame(cp = c("cycle1.ch1", "cycle2.ch1", "cycle2.ch2", "cycle2.ch3",
			      "cycle3.ch1", "cycle3.ch2", "cycle4.ch2", "cycle4.ch3"),
		       ashlr = c("_ch1", "_ch3", "_ch4", "_ch5", "_ch7", "_ch8", "_ch12", "_ch13"),
		       elem = c("EGFP", "OLLAS", "V5", "HA", "ALFA", "FLAG", "MYC", "VSVG"))
elemGrps$cycle <- substr(elemGrps$cp, 1, 6)
elems <- bind_rows(elemGrps[, c("cp", "elem")] %>% set_colnames(., c("meth", "elem")),
		   elemGrps[, c("ashlr", "elem")] %>% mutate(ashlr = paste0("mean", ashlr)) %>% 
			  set_colnames(., c("meth", "elem")))

tcMap <- read.table("../../../00.config/tcMap.txt", header = TRUE)
barcodes <- read.table("../../../00.config/barcodes.txt", header = TRUE)
wells <- read.table("../../../00.config/wells.txt")
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


wSet <- unique(v$well)
toPlot <- list()
gfpPos <- cpBcdCalls$objId[cpBcdCalls$cycle1.ch1 == 1]
for (nxCh in 1:nrow(elemGrps)) {
	wSet <- subset(tcMap, grepl(elemGrps$elem[nxCh], bcdElems) | is.na(bcdElems))$well
	if (elemGrps$elem[nxCh] == "EGFP") wSet <- unique(v$well)
	nx_cp <- log2(v[(v$well %in% wSet), elemGrps$cp[nxCh]])
	nx_ashlr <- log2(v[(v$well %in% wSet), paste0("mean", elemGrps$ashlr[nxCh])]*65535)
	toPlot[[nxCh]] <- data.frame(elem = elemGrps$elem[nxCh],
				     cellPool = nx_cp,
				     ASHLAR = nx_ashlr,
				     category = c("background", "signal")[
					1 + (v$objId[(v$well %in% wSet)] %in% gfpPos)])
}
toPlot %<>% bind_rows()
toPlot %<>% merge(., elemGrps[, c("elem", "cycle")]) %>% 
	gather(., "tool", "intensity", -elem, -category, -cycle)

df1 <- subset(toPlot, category == "signal")
df2 <- data.frame()
for (nxElm in unique(toPlot$elem)) {
	nxDat <- subset(toPlot, (elem == nxElm) & (category == "background")) %>% 
		group_by(tool) %>% 
		sample_n(., size = floor(0.3*sum(df1$elem == nxElm)), replace = TRUE)
	df2 %<>% bind_rows(., nxDat)
}
toPlot1 <- bind_rows(df1, df2)
write.table(toPlot1, file = "dat_for_s_vs_b_ashlar_vs_cellpool.txt",
	    sep = "\t", row.names = FALSE, quote = FALSE)
p <- ggplot(toPlot1, aes(color = tool, x = intensity, y = after_stat(ndensity))) +
	geom_density() + 
	theme_classic() +
	scale_color_manual(values = c("grey", "black")) +
	facet_nested_wrap(vars(cycle, elem), nrow = 2) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) + xlab(NULL) + ylab(NULL) +
	theme(legend.position = "top", legend.direction = "horizontal",
	      text = element_text(family = "Arial", size = 5),
	      legend.key.size = unit(0.2, "lines"),
	      panel.background = element_rect(fill = "transparent", colour = NA),
	strip.text = element_text(margin = margin(t = 2, b = 2, l = 0, r = 0)),
	panel.grid = element_blank(),
	panel.border = element_blank(),
	legend.title = element_blank(), 
	legend.text = element_text(margin = margin(l = 2, r = 0)),
	legend.spacing = unit(0, "pt"),
	axis.text = element_text(color = "black")) 

ggsave(paste0("cpVsAshlr_", Sys.time(), ".png"), width = 3.4, height = 2.3)


