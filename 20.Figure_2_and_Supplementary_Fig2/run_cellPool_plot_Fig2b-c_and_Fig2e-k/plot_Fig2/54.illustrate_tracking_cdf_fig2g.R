library(tidyverse)
library(magrittr)

well <- "A08"
cyc1 <- read.table("../40.segmentation/iter1.step3.features/A08-cyc1.HOECHST_33342.txt",
		   header = TRUE, sep = "\t")
cyc2 <- read.table("../40.segmentation/iter1.step3.features/A08-cyc2.HOECHST_33342.txt",
                   header = TRUE, sep = "\t")
trackDat <- read.table("../50.trackObjs/A08.cycle1_cycle2_tracking.stage3_local_refs.txt", 
		       header = TRUE, sep = "\t") 
ngbrDat <- read.table("../50.trackObjs/A08-cyc1_nearest_nghbrs_stats.txt", 
		      header = TRUE, sep = "\t")
mateDat <- merge(cyc1, trackDat, by.x = "label", by.y = "cycle1") %>%
	merge(., cyc2, by.x = "cycle2", by.y = "label")
rndmDat <- bind_cols(ngbrDat[, "label", drop = FALSE], 
		     data.frame(rndmMate = apply(ngbrDat[, paste0("ngbrLabel", 1:5)], 1, 
			   function(x) sample(x, 1)))) %>%
	merge(., cyc1, by = "label") %>%
	merge(., cyc1, by.x = "rndmMate", by.y = "label")

toPlot <- bind_rows(data.frame(type = "tracking", 
			       areaDiff = log10(mateDat$area.x/mateDat$area.y),
			       orntDiff = (mateDat$orientation.x - mateDat$orientation.y)),
		    data.frame(type = "random",
			       areaDiff = log10(rndmDat$area.x/rndmDat$area.y),
			       orntDiff = (rndmDat$orientation.x - rndmDat$orientation.y)))
toPlot <- gather(toPlot, "prop", "val", -type)

p <- ggplot(toPlot, aes(x = val, color = type)) +
	ggh4x::facet_grid2(col = vars(prop), scales = "free_x", independent = "x") +
	stat_ecdf(geom = "step", pad = FALSE) + theme_classic() +
       scale_color_manual(values = c("grey", "black"))+
	theme(legend.position = "top", legend.direction = "horizontal",
	      legend.title = element_blank(),
	      text = element_text(size = 7),
	      legend.text = element_text(margin = margin(l = 0, r = 0)),
        legend.spacing = unit(0, "pt"),
        panel.grid = element_blank(),
        #legend.key.size = unit(0.2, "cm"),
        panel.border = element_blank(),
        plot.margin = margin(t = 0, b = 0, l = 0, r = 0)) + xlab(NULL) + ylab(NULL) +
	scale_y_continuous(breaks = c(0, 0.5, 1))

ggsave("ecdf_features_with_legend.pdf", width = 1.94, height = 1.46)
ggsave("ecdf_features.pdf", width = 1.94, height = 1.46, p + guides(color = "none"))


p <- ggplot(toPlot, aes(x = val, color = type)) +
	ggh4x::facet_grid2(col = vars(prop), scales = "free_x", independent = "x") +
	geom_density(aes(y = after_stat(scaled))) +
	theme_classic() +
       scale_color_manual(values = c("grey", "black"))+
	theme(legend.position = "top", legend.direction = "horizontal",
	      legend.title = element_blank(),
	      text = element_text(size = 7),
	      legend.text = element_text(margin = margin(l = 0, r = 0)),
        legend.spacing = unit(0, "pt"),
        panel.grid = element_blank(),
        #legend.key.size = unit(0.2, "cm"),
        panel.border = element_blank(),
        plot.margin = margin(t = 0, b = 0, l = 0, r = 0)) + xlab(NULL) + ylab(NULL) +
	scale_y_continuous(breaks = c(0, 0.5, 1))

ggsave("pdf_features_with_legend.pdf", width = 1.94, height = 1.46)
ggsave("pdf_features.pdf", width = 1.94, height = 1.46, p + guides(color = "none"))








