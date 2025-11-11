library(tidyverse)
library(magrittr)
library(rjson)

params <- fromJSON(file = "../00.config/config.json") #args[1])

dat <- read.table("../60.barcodes/integrated_intensity_table.txt.gz", header = TRUE, sep = "\t")
chToElem <- read.table("../60.barcodes/chToElem.txt", sep = "\t", header = TRUE)
tcMap <- file.path(params$analysisDir, params$tcMapFile) %>%
  read.table(., header = TRUE, sep = "\t")

calls <- dat[, c("objId", chToElem$channel)]
calls[, names(calls) %>% subset(., startsWith(., "cycle"))] <- NA

bkbnThresh <- readLines("../60.barcodes/backbone_thresholds.txt") %>%
  as.numeric()
bkbnCh <- chToElem$channel[chToElem$element == params$debarcodingParams$backboneProt] %>%
  str_replace(., "-", ".")
calls[, bkbnCh] <- ifelse(dat[, bkbnCh] > bkbnThresh[2], 1, NA)

def_neg <- dat[, bkbnCh] < bkbnThresh[1]
calls[which(def_neg), chToElem$channel] <- 0

bkbnPos <- subset(calls, get(bkbnCh) == 1) %$%
  objId
densityRatioThresh <- params$debarcodingParams$thresholdDensityRatio
minRunLength <- params$debarcodingParams$thresholdRunLengthMin

chToElem$channel %<>% str_replace(., "-", ".")
x <- "cycle2.ch1"
scaleDensities <- FALSE
print(x)
mixture <- log2(dat[, x])

signalEnriched <- subset(dat, (objId %in% bkbnPos)) %$%
  get(x) %>%
  log2()

background <- subset(dat, !(objId %in% bkbnPos) & !is.na(get(bkbnCh))) %$%
  get(x) %>%
  log2()

mixHist <- hist(c(mixture, background, signalEnriched), nclass = 100, plot = FALSE) %$%
  data.frame(xstart = breaks[-length(breaks)],
             xend = breaks[-1], density = density)

thisBreaks <- mixHist %$%
  c(xstart, xend) %>%
  unique() %>%
  sort()

sigEnrichedHist <- hist(signalEnriched, breaks = thisBreaks, plot = FALSE) %$%
  data.frame(xstart = breaks[-length(breaks)],
             xend = breaks[-1], densitySigEnriched = density)

bckHist <- hist(background, breaks = thisBreaks, plot = FALSE) %$%
  data.frame(xstart = breaks[-length(breaks)],
             xend = breaks[-1], densityBckgd = density)

this <- merge(sigEnrichedHist, bckHist,
              by = c("xstart", "xend"), all = TRUE)
this$densitySigEnriched[this$densitySigEnriched %>% is.na()] = 0
this$densityBckgd[this$densityBckgd %>% is.na()] = 0
this %<>%
  mutate(densityRatio = densitySigEnriched/densityBckgd) %>%
  arrange(xstart)

if (scaleDensities) this$densityRatio %<>%
  multiply_by(., max(this$densityBckgd)/max(this$densitySigEnriched))

thisRle <- this %$%
  is_weakly_greater_than(densityRatio, densityRatioThresh) %>%
  rle()
firstGoodRun <- which(thisRle$values & (thisRle$lengths >= minRunLength)) %>%
  min()
threshRow <- sum(thisRle$lengths[1:(firstGoodRun - 1)]) + 1

this[threshRow, "xstart"]

toPlot <- gather(this, "type", "val", -xstart, -xend)
toPlot$type[toPlot$type == "densityBckgd"] <- "background\ndensity"
toPlot$type[toPlot$type == "densitySigEnriched"] <- "mixture\ndensity"
toPlot$type[toPlot$type == "densityRatio"] <- "mixture /\nbackground"
toPlot$type %<>% factor(levels = c("background\ndensity","mixture\ndensity", "mixture /\nbackground"))
p <- ggplot(toPlot, aes(x = (xstart + xend)/2, y = val)) + #, color = type)) +
  geom_line() +
  theme_classic() + ggh4x::facet_grid2(row = vars(type),
                                       scales = "free_y",  independent = "y") +
  theme(legend.position = "top", legend.direction = "horizontal",
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0)) +
  geom_vline(xintercept = this[threshRow, "xstart"], color = "red", linetype = "dashed") +
  xlab(NULL) + ylab(NULL)

ggsave(paste0("illustrate_thresholding_", x, "_",
              str_replace_all(Sys.time(), " ", "_"), ".pdf"), height = 2, width = 1.44)


