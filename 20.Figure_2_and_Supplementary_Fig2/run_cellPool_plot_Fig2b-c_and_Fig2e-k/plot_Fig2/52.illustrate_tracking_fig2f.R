library(tidyverse)
library(magrittr)

well <- "A08" 
xstart = 4000
xend = 4400
ystart = 5000
yend = 5400

c1 <- read.table("../40.segmentation/iter1.step3.features/A08-cyc1.HOECHST_33342.txt", 
		 header = TRUE, sep = "\t")

c2 <- read.table("../40.segmentation/iter1.step3.features/A08-cyc2.HOECHST_33342.txt", 
                 header = TRUE, sep = "\t") 

t1 <- subset(c1, (ystart < centroid.0) & (yend > centroid.0) & 
	     (xstart < centroid.1) & (xend > centroid.1))
t2 <- subset(c2, (ystart < centroid.0) & (yend > centroid.0) & 
	     (xstart < centroid.1) & (xend > centroid.1))

p <- ggplot(t1, aes(x = centroid.1, y = centroid.0)) +
	geom_point() + scale_y_reverse() + 
	theme_classic() +
	geom_text(aes(label = label))

ggsave("setA_cyc1.pdf", p)

p <- ggplot(t2, aes(x = centroid.1, y = centroid.0)) +
	geom_point() + scale_y_reverse() + 
	theme_classic() +
	geom_text(aes(label = label))

ggsave("setA_cyc2.pdf", p)

ilstrtDat <- data.frame(dpLab = c(1:5, "x"), label = c(8991, 9011, 8931, 8934, 8968, 8536))
xLab <- 8536

toPlot <- bind_rows(subset(t1, label %in% ilstrtDat$label, select = c("area", "label")),
	subset(t2, label == xLab, select = c("area", "label"))) %>% 
	merge(., ilstrtDat)
toPlot$dpLab %<>% factor()
p <- ggplot(toPlot, aes(x = dpLab, size = area)) +
	geom_point(shape = 1, y = 0.5) + scale_y_continuous(expand = c(0, 0)) +
	theme_classic() +
	guides(size = "none") + theme(plot.margin = margin(0, 0, 0, 0)) +
	scale_size(range = c(1, 1.75))
ggsave("area_ilstrt_tracking_with_axis.pdf", p)
ggsave("area_ilstrt_tracking.pdf", p + theme_void(), width = 1.2, height = 0.2)

toPlot <- bind_rows(subset(t1, label %in% ilstrtDat$label, select = c("orientation", "label", 
								      "centroid.0", "centroid.1")),
        subset(t2, label == xLab, 
	       select = c("orientation", "label", "centroid.0", "centroid.1"))) %>%
        merge(., ilstrtDat)
toPlot$dpLab %<>% factor()
toPlot$centroid.1 <- as.numeric(toPlot$dpLab)
toPlot$centroid.0 <- 0.5
fctr <- 1
toPlot$x1 <- toPlot$centroid.1 + cos(toPlot$orientation) * 0.5 *fctr
toPlot$y1 <- toPlot$centroid.0 - sin(toPlot$orientation) * 0.5 *fctr
toPlot$x2 <- toPlot$centroid.1 - sin(toPlot$orientation) * 0.5 *fctr
toPlot$y2 <- toPlot$centroid.0 - cos(toPlot$orientation) * 0.5 *fctr
p <- ggplot(toPlot, aes(x = centroid.1, y = centroid.0)) +
	geom_segment(aes(x = centroid.1 + sin(orientation) * 0.5 *fctr, 
			 y = centroid.0 + cos(orientation) * 0.5 *fctr, 
			 xend = x2, yend = y2)) + coord_fixed(ratio = 1) +
	scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) 

ggsave("ornt_ilstrt_tracking_with_axis.pdf", p)
ggsave("ornt_ilstrt_tracking.pdf", p + theme_void(), width = 1.2, height = 0.2)


trin1 <- 9021
trin2 <- 8589
xTrin1 <- subset(t1, label == trin1)$`centroid.1` - subset(t2, label == trin2)$`centroid.1`
yTrin1 <- subset(t1, label == trin1)$`centroid.0` - subset(t2, label == trin2)$`centroid.0`
toPlot <- bind_rows(subset(t1, label %in% ilstrtDat$label, select = c("label",
                                                                      "centroid.0", "centroid.1")),
        subset(t2, label == xLab,
               select = c("label", "centroid.0", "centroid.1"))) %>%
        merge(., ilstrtDat)
toPlot$dpLab %<>% factor()
toPlot$xi <- as.numeric(toPlot$dpLab)
toPlot$yi <- 0.5
toPlot$orientation <- atan((toPlot$centroid.0 - subset(t2, label == xLab)$centroid.0)/
			   (toPlot$centroid.1 - subset(t2, label == xLab)$centroid.1))
toPlot$orientation[toPlot$xi == 6] <- atan(yTrin1/xTrin1)

fctr <- 1
toPlot$len <- (toPlot$centroid.0 - subset(t2, label == xLab)$centroid.0)^2 +
	(toPlot$centroid.1 - subset(t2, label == xLab)$centroid.1)^2
toPlot$len[toPlot$xi == 6] <- yTrin1^2 + xTrin1^2
toPlot$len <- toPlot$len/max(toPlot$len)
toPlot$x2 <- toPlot$xi + cos(toPlot$orientation) * 0.5 *toPlot$len
toPlot$y2 <- toPlot$yi + sin(toPlot$orientation) * 0.5 *toPlot$len
toPlot$xend = toPlot$xi - cos(toPlot$orientation) * 0.5 *toPlot$len
toPlot$yend = toPlot$yi - sin(toPlot$orientation) * 0.5 *toPlot$len
toPlot$swap <- toPlot$xend
toPlot$xend[as.character(toPlot$dpLab) %in% c("4", "5", "x")] <- 
	toPlot$x2[as.character(toPlot$dpLab) %in% c("4", "5", "x")]
toPlot$x2[as.character(toPlot$dpLab) %in% c("4", "5", "x")] <-
	toPlot$swap[as.character(toPlot$dpLab) %in% c("4", "5", "x")]
toPlot$swap <- toPlot$yend
toPlot$yend[as.character(toPlot$dpLab) %in% c("4", "5", "x")] <-
        toPlot$y2[as.character(toPlot$dpLab) %in% c("4", "5", "x")]
toPlot$y2[as.character(toPlot$dpLab) %in% c("4", "5", "x")] <-
        toPlot$swap[as.character(toPlot$dpLab) %in% c("4", "5", "x")]
p <- ggplot(toPlot) +
        geom_segment(aes(xend = xend,
                         yend = yend,
                         x = x2, y = y2),
		     arrow = arrow(length = unit(0.3, "npc"))) + coord_fixed(ratio = 1) +
	scale_y_reverse(expand = c(0.05, 0.05)) + scale_x_continuous(expand = c(0.05, 0.05))

ggsave("shift_ilstrt_tracking_with_axis.pdf", p)
ggsave("shift_ilstrt_tracking.pdf", p + theme_void(), width = 1.2, height = 0.2)


trin1 <- 9021
trin2 <- 8589
xTrin1 <- subset(t1, label == trin1)$`centroid.1` - subset(t2, label == trin2)$`centroid.1`
yTrin1 <- subset(t1, label == trin1)$`centroid.0` - subset(t2, label == trin2)$`centroid.0`
toPlot <- bind_rows(subset(t1, label %in% ilstrtDat$label, select = c("label",
                                                                      "centroid.0", "centroid.1")),
        subset(t2, label == xLab,
               select = c("label", "centroid.0", "centroid.1"))) %>%
        merge(., ilstrtDat)
toPlot$dpLab %<>% factor()
toPlot$xi <- as.numeric(toPlot$dpLab)
toPlot$yi <- 0.5
toPlot$orientation <- atan((toPlot$centroid.0 - subset(t1, label == trin1)$centroid.0)/
			   (toPlot$centroid.1 - subset(t1, label == trin1)$centroid.1))
toPlot$orientation[toPlot$xi == 6] <- atan((subset(t2, label == xLab)$centroid.0 -
					    subset(t2, label == trin2)$centroid.0)/
	(subset(t2, label == xLab)$centroid.1 - subset(t2, label == trin2)$centroid.1))

fctr <- 1
toPlot$len <- (toPlot$centroid.0 - subset(t1, label == trin1)$centroid.0)^2 +
	(toPlot$centroid.1 - subset(t1, label == trin1)$centroid.1)^2
toPlot$len[toPlot$xi == 6] <- (subset(t2, label == xLab)$centroid.0 -
                                            subset(t2, label == trin2)$centroid.0)^2 +
	(subset(t2, label == xLab)$centroid.1 - subset(t2, label == trin2)$centroid.1)^2
toPlot$len <- toPlot$len/max(toPlot$len)
toPlot$x2 <- toPlot$xi + cos(toPlot$orientation) * 0.5 * toPlot$len
toPlot$y2 <- toPlot$yi + sin(toPlot$orientation) * 0.5 * toPlot$len
p <- ggplot(toPlot) +
        geom_segment(aes(xend = xi - cos(orientation) * 0.5 * len,
                         yend = yi - sin(orientation) * 0.5 * len,
                         x = x2, y = y2),
		     arrow = arrow(length = unit(0.3, "npc"))) + coord_fixed(ratio = 1) +
        scale_y_reverse(expand = c(0.05, 0.05)) + scale_x_continuous(expand = c(0.05, 0.05))

ggsave("pos_ilstrt_tracking_with_axis.pdf", p)
ggsave("pos_ilstrt_tracking.pdf", p + theme_void(), width = 1.2, height = 0.2)

