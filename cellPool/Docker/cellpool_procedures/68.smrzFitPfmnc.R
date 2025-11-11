library(tidyverse)
library(magrittr)
library(purrr)
library(rjson)

args <- commandArgs(trailingOnly = TRUE)
params <- fromJSON(file = args[1])

emTblDir <- file.path(params$analysisDir, params$barcodeAnalysisDir)
cuts <- 1:params$debarcodingParams$nCuts
barcodes <- file.path(emTblDir, "barcodes.refmt.txt") %>% 
	read.table(., sep = "\t", header = TRUE) %>%
	mutate(bcdIx = 0:(n() - 1))
logLFiles <- list.files(emTblDir) %>% 
	subset(., grepl("logL", .))

logLDat <- expand.grid(cut = cuts,
        init = strsplit(logLFiles, "_init_") %>%
                        plyr::laply(., magrittr::extract, i = 2) %>%
                        strsplit(., split = "_") %>%
                        plyr::laply(., magrittr::extract, i = 1) %>%
                        str_remove_all(".txt") %>%
                        unique() %>%
                        subset(., !is.na(.)))
logLDat <- apply(logLDat, 1, function(x) {
	thisDat <- data.frame(cut = x[1], init = x[2])
	thisFiles <- subset(logLFiles, grepl(paste0("cut_", trimws(as.character(x[1])), 
			"_.*init_", x[2], "_"), logLFiles)) %>% 
		set_names(., strsplit(., split = "iter_") %>% 
			plyr::laply(., magrittr::extract, i = 2) %>%
			strsplit(., split = "_") %>%
			plyr::laply(., magrittr::extract, i = 1) %>%
			str_remove_all(., ".txt"))
	thisLogLs <- map(file.path(emTblDir, thisFiles), read.table, sep = "\t", header = TRUE) %>%
		bind_rows() %>% 
		bind_cols(thisDat, data.frame(iter = names(thisFiles)), .)
}) %>% 
	bind_rows() %>%
	set_rownames(NULL)

logLDat$totalLogL <- logLDat[, paste0("X", barcodes$bcdIx)] %>% 
	apply(., 1, sum)
logLDat$iter %<>% as.numeric()
logLDat$init %<>% as.numeric() 
logLDat$cut %<>% as.character() %>% trimws() %>% paste0("C", .)

toPlot <- logLDat %>%
	gather(., "variable", "value", -cut, -init, -iter) %>% 
	subset(., variable == "totalLogL")
toPlot %<>%
	group_by(cut) %>% 
	mutate(value = value %>% 
		subtract(., subset(., is.infinite(.) %>% 
					not()) %>% 
					min(., na.rm = TRUE)) %>%
					divide_by(., max(.)))
toPlot$cut %<>% 
	factor(., levels = 
		paste0("C", as.character(1:params$debarcodingParams$nCuts) %>% trimws()))	
p <- ggplot(toPlot, aes(x = iter + 1, y = value)) +
	geom_line() +
	geom_point() +
	facet_grid(cut ~ init) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) +
	scale_x_continuous(breaks = c(1, max(toPlot$iter) + 1)) +
	ylab("Likelihood gain") +
	xlab("Iteration of Expectation-Maximization")

ggsave(paste0("EM_totalLogL_", str_replace_all(Sys.time(), " ", "_"), ".png"), 
	p, width = 10, height = 8)

pdf(paste0("EM_barcodeWiseLogL_", str_replace_all(Sys.time(), " ", "_"), ".pdf"), 
	width = 10, height = 8)
for (idx in barcodes$bcdIx) {
	print(idx)
	barcode_elems <- subset(barcodes, bcdIx == idx,
		select = colnames(barcodes) %>% grepl("epitope", .)) %>%
		paste0(., collapse = "_")
	toPlot <- logLDat %>% 
		gather(., "variable", "value", -cut, -init, -iter) %>% 
		subset(., variable == paste0("X", idx))
	toPlot %<>%
		group_by(cut) %>% 
		mutate(value = value %>% 
			subtract(., subset(., is.infinite(.) %>% 
						not()) %>% 
						min(., na.rm = TRUE)) %>%
						divide_by(., max(.)))
		
	p <- ggplot(toPlot, aes(x = iter + 1, y = value)) +
		geom_line() +
		geom_point() +
		facet_grid(cut ~ init) +
		scale_y_continuous(breaks = c(0, 0.5, 1)) +
		scale_x_continuous(breaks = c(1, max(toPlot$iter) + 1)) +
		ylab("Likelihood gain") +
		xlab("Iteration of Expectation-Maximization") +
		ggtitle(barcode_elems)
	print(p)	
}
dev.off()

pdf(paste0("EM_cutWiseLogL_", str_replace_all(Sys.time(), " ", "_"), ".pdf"), 
	width = 10, height = nrow(barcodes))
for (qn in unique(logLDat$cut)) {
	print(qn)
	barcode_elems <- apply(barcodes[, colnames(barcodes) %>% grepl("epitope", .)], 1,
				paste0, collapse = "_")
	toPlot <- logLDat
	colnames(toPlot)[colnames(toPlot) %in% paste0("X", barcodes$bcdIx)] <- barcode_elems
	toPlot <- toPlot %>% 
		gather(., "variable", "value", -cut, -init, -iter) %>% 
		subset(., variable %in% barcode_elems) %>%
		subset(., cut == qn)
	toPlot %<>%
		group_by(variable) %>% 
		mutate(value = value %>% 
			subtract(., subset(., is.infinite(.) %>% 
						not()) %>% 
						min(., na.rm = TRUE)) %>%
						divide_by(., max(.)))
		
	p <- ggplot(toPlot, aes(x = iter + 1, y = value)) +
		geom_line() +
		geom_point() +
		facet_grid(variable ~ init) +
		scale_y_continuous(breaks = c(0, 0.5, 1)) +
		scale_x_continuous(breaks = c(1, max(toPlot$iter) + 1)) +
		ylab("Likelihood gain") +
		xlab("Iteration of Expectation-Maximization") +
		ggtitle(qn) +
		theme(strip.text.y.left = element_text(angle = 0))
	print(p)	
}
dev.off()

logLDat %>% 
	group_by(cut) %>% 
	summarize(iterMax = iter[which.max(totalLogL)], initMax = init[which.max(totalLogL)]) %>%
	write.table(., file = "init_and_iter_with_max_logL.txt", 
		sep = "\t", row.names = FALSE, quote = FALSE)

#--------------------------------------
#Plot mixing proportions
#--------------------------------------
piFiles <- list.files(emTblDir) %>%
        subset(., grepl("_pi.txt", .))
piDat <- expand.grid(cut = cuts,
        init = strsplit(piFiles, "_init_") %>%
                        plyr::laply(., magrittr::extract, i = 2) %>%
                        strsplit(., split = "_") %>%
                        plyr::laply(., magrittr::extract, i = 1) %>%
                        str_remove_all(".txt") %>%
                        unique() %>%
                        subset(., !is.na(.)))
piDat <- apply(piDat, 1, function(x) {
        thisDat <- data.frame(cut = x[1], init = x[2])
        thisFiles <- subset(piFiles, grepl(paste0("cut_", trimws(as.character(x[1])),
                        "_.*init_", x[2], "_iter"), piFiles)) %>%
                set_names(., strsplit(., split = "iter_") %>%
                        plyr::laply(., magrittr::extract, i = 2) %>%
                        strsplit(., split = "_") %>%
                        plyr::laply(., magrittr::extract, i = 1) %>%
                        str_remove_all(., ".txt"))
        thisLogLs <- map(file.path(emTblDir, thisFiles), read.table, sep = "\t", header = TRUE) %>%
                bind_rows() %>%
                bind_cols(thisDat, data.frame(iter = names(thisFiles)), .)
        thisLogLs
}) %>%
        bind_rows() %>%
        set_rownames(NULL)

toPlot <- piDat
barcode_elems <- apply(barcodes[, colnames(barcodes) %>% grepl("epitope", .)], 1,
                                paste0, collapse = "_")
colnames(toPlot)[colnames(toPlot) %in% paste0("X", barcodes$bcdIx)] <- barcode_elems
toPlot <- gather(toPlot, "variable", "value", -cut, -init, -iter)
toPlot$iter %<>% as.numeric()
toPlot$init %<>% as.numeric()
toPlot$cut %<>% as.character() %>% trimws() %>% paste0("C", .)

pdf(paste0("EM_mixing_proportions_", str_replace_all(Sys.time(), " ", "_"), ".pdf"), 
	width = 2*sqrt(nrow(barcodes)), height = 2*sqrt(nrow(barcodes)))
for (thisCut in cuts) {
        p <- ggplot(toPlot %>% subset(., cut == paste0("C", thisCut)),
                        aes(x = iter + 1, y = value, group = init)) +
                geom_line() +
                facet_wrap(vars(variable)) +
                scale_y_log10() +
                scale_x_continuous(breaks = c(1, max(toPlot$iter) + 1)) +
                ylab("Mixing proportion") +
                xlab("Iteration of Expectation-Maximization") +
		ggtitle(paste0("Cut: ", thisCut))
        print(p)
}
dev.off()


