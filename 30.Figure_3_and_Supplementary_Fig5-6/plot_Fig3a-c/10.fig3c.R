setwd("~/Documents/Analysis/30.MCF10A-203_M71-M74_Fig3b-c/figures_for_paper/3b_illustrate_imaging_panel")

library(tidyverse)
library(magrittr)
library(knitr)

well <- "E10"
pool <- 73
objId <- 12200

ftrFiles = list.files("../../40.segmentation/iter1.step3.features/") %>%
  subset(., startsWith(., well)) %>%
  subset(., endsWith(., ".HOECHST_33342.txt"))
ftr = list()
for (fl in ftrFiles) {
  ftr[[substr(fl, 1, 8)]] <-
    read.table(file.path("../../40.segmentation/iter1.step3.features/",
                          fl), sep = '\t', header = TRUE)
}

r <- as.integer(ftr[[paste0(well, '-cyc1')]][
  ftr[[paste0(well, '-cyc1')]]$label == objId, 'centroid.0'])
c <- as.integer(ftr[[paste0(well, '-cyc1')]][
  ftr[[paste0(well, '-cyc1')]]$label == objId, 'centroid.1'])
saturated_pixels <- 0.3 #Percent of saturated pixels in the image
height <- 940
width <- 600

nucInROI <- ftr[[paste0(well, '-cyc1')]] %>%
  subset(., centroid.0 > r - as.integer(height/2)) %>%
  subset(., centroid.0 < r + as.integer(height/2)) %>%
  subset(., centroid.1 > c - as.integer(width/2)) %>%
  subset(., centroid.1 < c + as.integer(width/2))


intDat <- read.table('../../60.barcodes/integrated_intensity_table.txt.gz', 
                     sep = '\t', header = TRUE)
intROI <- intDat %>%
  subset(., objId %in% paste0(well, "n", nucInROI$label))

ggplot(nucInROI, aes(x = centroid.1, y = centroid.0)) +
  geom_point() +
  scale_y_reverse() +
  theme_classic() +
  coord_fixed()

bcdIds <- read.table("../../00.config/barcodes.txt",
                     header = TRUE, sep = "\t") %>%
  mutate(bcdElems = paste(epi1, epi2, epi3, sep = "_")) %>%
  select(id, bcdElems)

libMap <- read.table("../../00.config/libraryMap.txt",
                     header = TRUE, sep = "\t") %>%
  subset(., poolId == pool) %>%
  merge(., bcdIds, by.x = "bcdId", by.y = "id")

bcds <- read.table('../../60.barcodes/barcode_classifications.txt.gz', 
                   sep = '\t', header = TRUE)
bcds %<>%
  subset(objId %in% paste0(well, "n", nucInROI$label)) %>%
  merge(., libMap, by = "bcdElems", all.x = TRUE) %>%
  select(objId, bcdElems, barcodeSep, bcdId, poolId, sgId, gene) %>%
  merge(., nucInROI %>% 
          select(label, centroid.0, centroid.1) %>%
          mutate(label = paste0(well, "n", label)),
        by.x = "objId", by.y = "label", all.x = TRUE) #%>%
  # mutate(gene = replace(gene, is.na(gene), "X")) %>%
  # mutate(barcodeSep = replace(barcodeSep, gene == "X", 0.2))

p <- ggplot(bcds, aes(x = centroid.1, y = centroid.0, label = gene,
                 alpha = barcodeSep)) +
  geom_text(size = 1, show.legend = FALSE) +
  scale_y_reverse() +
  theme_void() +
  coord_fixed(expand = FALSE) + 
  scale_alpha(range = c(0, 1)) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_rect(color = "black", fill = NA, 
        #                             linewidth = 0.2),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave("fig3c.pdf", p, height = 1.57, width = 1) 

ggsave("fig3c.pdf", p, device = cairo_pdf, height = 1.57, width = 1)  
knitr::plot_crop("fig3c.pdf")

