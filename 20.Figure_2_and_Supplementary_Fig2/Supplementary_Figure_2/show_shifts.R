setwd("~/Documents/Analysis/10.U2OS-tinypool/Supplementary_Figure_2")
library(magrittr)
library(tidyverse)

mosaicDir <- "../30.mosaics"
filesToRead <- list.files(mosaicDir) %>%
  subset(., endsWith(., ".tileCentroidCoords.txt")) %>%
  data.frame(fileName = ., well = substr(., 1, 3))	
wells <- c("F11", "G08", "A09")

toPlot <- data.frame()
for (nxWl in wells) {
  for (cycI in 1:3) {
    cycJ <- cycI + 1
    ftrI <- file.path("../40.segmentation/iter1.step3.features", 
                      paste0(nxWl, "-cyc", cycI, ".HOECHST_33342.txt")) %>%
      read.table(., header = TRUE, sep = "\t") %>%
      dplyr::select(label, centroid.0, centroid.1)
    
    ftrJ <- file.path("../40.segmentation/iter1.step3.features", 
                      paste0(nxWl, "-cyc", cycJ, ".HOECHST_33342.txt")) %>%
      read.table(., header = TRUE, sep = "\t") %>%
      dplyr::select(label, centroid.0, centroid.1)
    
    trkDat <- file.path("../50.trackObjs", 
                        paste0(nxWl, ".cycle", cycI, "_cycle", cycJ,
                               "_tracking.stage3_local_refs.txt")) %>%
      read.table(., header = TRUE, sep = "\t")
    
    newPlot <- merge(trkDat, ftrI, by.x = paste0("cycle", cycI),
                    by.y = "label") %>%
      merge(., ftrJ, by.x = paste0("cycle", cycJ), by.y = "label",
            suffixes = c(paste0("_cycle", cycI), paste0("_cycle", cycJ))) %>%
      set_colnames(., colnames(.) %>%
                     str_replace_all(., paste0("cycle", cycJ), "cycleJ") %>%
                     str_replace_all(., paste0("cycle", cycI), "cycleI")) %>%
      mutate(shift = sqrt((centroid.1_cycleI - centroid.1_cycleJ)^2 +
                            (centroid.0_cycleI - centroid.0_cycleJ)^2),
             transition = paste0("cycle", cycJ, "/\ncycle", cycI),
             well = nxWl) 
    
    scaling <- newPlot %$%
      max(centroid.1_cycleJ, centroid.0_cycleJ,
          centroid.1_cycleI, centroid.0_cycleI)
   newPlot %<>%
      mutate(centroid.1_cycleJ = centroid.1_cycleJ/scaling,
             centroid.0_cycleJ = centroid.0_cycleJ/scaling,
             centroid.1_cycleI = centroid.1_cycleI/scaling,
             centroid.0_cycleI = centroid.0_cycleI/scaling)
   
   toPlot %<>% bind_rows(newPlot, .)
  }
}

pdf("show_cell_dist2.pdf", height = 5, width = 3.7)
ggplot(toPlot %>%
         group_by(well, transition) %>%
         sample_n(., 100), 
       aes(x = centroid.1_cycleJ, y = centroid.0_cycleJ,
           color = shift)) +
  facet_grid(transition ~ well) +
  geom_segment(aes(x = centroid.1_cycleI, xend = centroid.1_cycleJ,
                   y = centroid.0_cycleI, yend = centroid.0_cycleJ),
               arrow = arrow(length = unit(0.15, "cm"))) +
  theme_classic() +
  xlab("position X of nucleus") +
  ylab("position Y of nucleus") +
  theme(legend.position = "top", 
        legend.direction = "horizontal") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1))
dev.off()
