library(tidyverse)
library(magrittr)
library(rjson)
library(ggrepel)

load("normalized_phenotypic_data.rds")

phenoDat$totalHoechst <- log10(phenoDat$area_ch4 * phenoDat$intensity_mean_ch4)
phenoDat$totalH2AX <- log10(phenoDat$area_ch3*phenoDat$intensity_mean_ch3)
phenoDat$totalpH3Ser10 <- log10(phenoDat$area_ch2*phenoDat$intensity_mean_ch2)
phenoDat$totalATub <- log10(phenoDat$intensity_mean_ch1)
lowCntNTC <- phenoDat %>% 
	subset(., grepl("non-", gene)) %>% 
	group_by(gene) %>% 
	summarize(count = n()) %>% 
	subset(., count < 100) %$%
	gene

for (prop in c("area_ch4", "totalH2AX", "totalpH3Ser10", "totalATub")) {
	res <- list(gene = c(), statistic = c(), parameter = c(), 
	            p.value = c(), estimate = c(), coverage_target = c(),
	            coverage_nontargeting = c())
	for (sg in unique(phenoDat$gene)) {
		if (is.na(sg)) next
		print(sg)
		thisPoolId <- subset(phenoDat, gene == sg) %$% unique(poolId)
		s1 <- subset(phenoDat, gene == sg, select = prop, drop = TRUE)
		s2 <- subset(phenoDat, (startsWith(gene, "non-")) &
					(poolId %in% c(71:74)) & 
					(gene != sg) & 
					!(gene %in% lowCntNTC ),
				select = prop, drop = TRUE)

		if (length(s1) < 100) {
			thisRes <- list(statistic = NA, parameter = NA, p.value = NA, estimate = c(NA, NA))
		} else  thisRes <- t.test(s1, s2)

		res$gene %<>% c(., sg)
		res$statistic %<>% c(., thisRes$statistic)
		res$parameter %<>% c(., thisRes$parameter)
		res$p.value %<>% c(., thisRes$p.value)
		res$estimate %<>% c(., thisRes$estimate[1] - thisRes$estimate[2])
		res$coverage_target %<>% c(., length(s1))
		res$coverage_nontargeting %<>% c(., length(s2))
	}
	res %<>% as.data.frame()
	res$FDR <- p.adjust(res$p.value, "BH")
	res$FWER <- p.adjust(res$p.value, "bonferroni")
	res$nonTargeting <- res$gene %>% startsWith(., "non")

	res$labText <- res$gene

	res$FDRcapped <- -log10(res$FDR)
	res$FDRcapped[res$FDRcapped > 3.5] <- 3.5
	res <- bind_rows(subset(res, !grepl("non-", gene)), subset(res, grepl("non-", gene)))

	ntcDat <- subset(phenoDat, (is.na(gene) | startsWith(gene, "non-")) & 
				(poolId %in% c(71:74)), 
				select = prop, drop = TRUE) %>% subset(., is.na(.) %>% not())
	semNTC <- sd(ntcDat)/sqrt(length(ntcDat))
	res$estimate <- res$estimate/semNTC
	res$labText[(abs(res$estimate) < 3) | (res$FDR > 0.05)] <- ""
	assign(paste0("res_gene_level_", prop), res)

	p <- ggplot(res, aes(x = estimate, y = FDRcapped,
			shape = nonTargeting, color = nonTargeting)) +
		geom_point() +
		scale_shape_manual(values = c(16, 4)) +
		scale_color_manual(values = c("black", "red")) +
		scale_size_manual(values = c(3, 6)) +
		theme_classic() +
		geom_text_repel(aes(label = labText), size = 25/16, max.overlaps = 20) + #25/14) +
		guides(shape = "none", color = "none") +
		theme(legend.title = element_blank(),
			axis.title = element_text(size = 5),
			text = element_text(size = 5)) +
		ylab(expression(-log[10]*(FDR))) +
		xlab(prop)

	ggsave(paste0(prop, "_ttest_", Sys.time(), ".png"), width = 30, height = 30)
	ggsave(paste0(prop, "_ttest_", Sys.time(), ".png"), width = 7, height = 7)

}

ls() %>% subset(., grepl("gene_level", .)) %>% 
	save(file = "a_univariate_diff_analysis_res.RData", list = .)







