library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly=TRUE)
minL <- as.numeric(args[1])
minLDelta <- as.numeric(args[2])
bkbnThr <- readLines("../60.barcodes/backbone_thresholds.txt") %>% as.numeric()

dat <- read.table("../60.barcodes/integrated_intensity_table.txt.gz", header = TRUE)
calls <- read.table("../60.barcodes/barcode_classifications.txt.gz", header = TRUE)
tcMap <- read.table("../00.config/tcMap.txt", header = TRUE)

gfpPos <- dat$objId[dat$`cycle1.ch2` > bkbnThr[2]]
outDir <- "tables"

#steps: 1. get the background distribution from pool = 0 wells
#steps: 2. convert the raw intensity values to percentile score based on negative wells
#steps: 3. for every barcode calulate the score as average of its epitopes
#steps: 4. for every cell, save the most "likely" barcode based on step 3 and barcodeSep.

#---------------------------------------
#Call epitrain next
neg_ctrl_wells <- subset(tcMap, treatment == "no_virus")$well
chs <- setdiff(names(dat), "objId")
bck_vals <- map(chs,
        function(x) subset(dat[, x], substr(dat$objId, 1, 3) %in% neg_ctrl_wells) %>%
                subset(., is.na(.) %>% not())) %>%
        set_names(chs)

scores <- dat[, c("objId", names(bck_vals))]
scores[, "cycle1.ch1"] <- as.numeric(scores$objId %in% gfpPos)
for (epi in chs) {
        if (epi %in% c("cycle1.ch1", "cycle1.ch2", "cycle2.ch4", "cycle3.ch4", "cycle4.ch4")) next
        print(epi)
        pctlThisBckgd <- quantile(bck_vals[[epi]], seq(0, 1, 0.005), na.rm = TRUE)
        toScore <- which(!is.na(scores[, epi]))
        scores[toScore, epi] <- map_dbl(dat[toScore, epi],
                function(x) {
                        lessPctls <- subset(pctlThisBckgd, pctlThisBckgd < x)
                        ifelse(!length(lessPctls), 0, lessPctls %>% names() %>%
                                str_remove_all("%") %>%
                                as.numeric() %>%
                                max())/100

                })
}

bcds <- read.table("../00.config/barcodes.txt", header = TRUE)
chToElem <- read.table("../00.config/chToElemId.txt", header = TRUE, sep = "\t")
chToStain <- read.table("../00.config/chToStainId.txt", header = TRUE, sep = "\t")

bcds <- bcds[, c("epi1", "epi2", "epi3")]
uniqEpis <- unique(bcds %>% unlist())
epiToCh <- map(uniqEpis, 
	       function(x) {
	       		intmd <- which(chToElem == x, arr.ind = TRUE)
	       		paste0("cycle", intmd[1, 1], ".ch", intmd[1, 2])
	       }) %>% 
	set_names(uniqEpis)

scoresCols <- setdiff(names(scores), 
		      c("cycle1.ch2", "cycle2.ch4", "cycle3.ch4", "cycle4.ch4", "objId"))
allTrainScores <- apply(scores[, scoresCols], 1,
        function(x) {
                if (is.na(x[scoresCols == "cycle1.ch1"])) return(rep(NA, 19))
                if (x[scoresCols == "cycle1.ch1"] == 0) return(rep(NA, 19))
                apply(bcds, 1, function(epis) {
                        mean((x[scoresCols %in% unlist(epiToCh[unlist(epis)])]) %>% as.numeric())
                })
        })
allTrainScores %<>% t()
allTrainScores %<>%
        set_colnames(., apply(bcds, 1, paste0, collapse = "_"))
callsBasedOnMaxScore <- apply(allTrainScores, 1,
        function(x) c(colnames(allTrainScores)[which.max(x)], x[which.max(x)])) %>% 
	t() %>%
	set_colnames(., c("bcd", "bcdSep")) %>% 
	as.data.frame() %>%
	mutate(bcdSep = as.numeric(bcdSep))

outDat <- data.frame(objId = scores$objId, callsBasedOnMaxScore)


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

outDat$well <- substr(outDat$objId, 1, 3)
outDat %<>% merge(., toMerge, all.x = TRUE)

write.table(outDat, "Kudo_et_al_2025_05_05.txt", row.names = FALSE, quote = FALSE, sep = "\t")

barcodeSepCutoffs <- 1 - (0:30)*0.01
yield <- c()
error_rate <- c()
for (ct in barcodeSepCutoffs) {
        bcdCalls <- subset(outDat, bcdSep >= ct)
        yield %<>% c(., sum(nrow(bcdCalls)/nrow(outDat)))
        error_rate %<>% c(., sum(bcdCalls$bcd != bcdCalls$bcdGrndTruth, na.rm = TRUE)/nrow(bcdCalls))
}

toPlot <- data.frame(cutoff = barcodeSepCutoffs, yield = yield, fdr = error_rate)

p <- ggplot(toPlot, aes(x = fdr, y = yield)) + geom_point()
ggsave("pr_take1_kudo.png")


p <- ggplot(toPlotDf, aes(x = fdr, y = yield, shape = meth)) + geom_point() + 
	theme_classic()
ggsave("pr_take1_kudo_cellpool.png")






