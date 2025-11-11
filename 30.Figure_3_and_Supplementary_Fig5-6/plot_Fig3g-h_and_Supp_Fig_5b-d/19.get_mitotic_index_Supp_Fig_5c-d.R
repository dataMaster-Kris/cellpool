library(tidyverse)
library(magrittr)

prop <- "totalpH3Ser10"
cutoff <- 6.25 
p <- ggplot(phenoDat, aes(x = get(prop))) +
  geom_histogram(bins = 500) +
  scale_y_log10() +
  geom_vline(xintercept = cutoff, linetype = "dotted", 
             color = "red") +
  theme_classic() + 
  xlab("Integrated nuclear H3S10P") +
  ylab("Frequency") +
  theme(text = element_text(size = 5))

ggsave("Thresholding_mitotic_cells.pdf", width = 2.15, height = 2.15)

mitIndex <- phenoDat %>% 
  mutate(sgId = replace(sgId, startsWith(sgId, "non-"), "non-targeting")) %>%
  group_by(sgId) %>% 
  summarize(mitotic = sum(get(prop) > cutoff, na.rm = TRUE),
            sgCov = n()) %>%
  mutate(mitIndex = mitotic/sgCov)

nullMI <- c(subset(mitIndex, startsWith(sgId, "non")) %$% 
              (sum(mitotic/sum(sgCov))), 
            mitIndex$mitIndex[is.na(mitIndex$sgId)]) %>% 
  mean()
nullMitotic <- subset(mitIndex, startsWith(sgId, "non") | is.na(sgId)) %$% 
  sum(mitotic)
nullAll <- subset(mitIndex, startsWith(sgId, "non") | is.na(sgId)) %$% 
  sum(sgCov)

toPlot <- list(gene = c(), sg1 = c(), sg2 = c(), 
               MI1 = c(), MI2 = c(), coverage = c())
for (gn in unique(phenoDat$gene) %>% subset(., !is.na(.)) %>%
     subset(., !startsWith(., "non"))) {
  toPlot$gene %<>% c(., gn)
  thisSgs <- phenoDat$sgId %>% unique %>% 
    subset(., startsWith(., paste0(gn, "_")))
  toPlot$sg1 %<>% c(., subset(thisSgs, endsWith(thisSgs, "1")))
  toPlot$sg2 %<>% c(., subset(thisSgs, endsWith(thisSgs, "2")))
  toPlot$MI1 %<>% 
    c(., mitIndex$mitIndex[which(mitIndex$sgId == tail(toPlot$sg1, 1))])
  toPlot$MI2 %<>% 
    c(., mitIndex$mitIndex[which(mitIndex$sgId == tail(toPlot$sg2, 1))])
  toPlot$coverage %<>% c(., subset(phenoDat, gene == gn) %$%
                           table(sgId) %>% mean())
  
  ok <- plyr::laply(toPlot, length) %>% 
    unique() %>% length() %>% equals(., 1)
  
  if (!ok) break
}

for (gn in c("non-targeting")) {
  toPlot$gene %<>% c(., gn)
  toPlot$sg1 %<>% c(., NA)
  toPlot$sg2 %<>% c(., NA)
  toPlot$MI1 %<>% 
    c(., mitIndex$mitIndex[which(mitIndex$sgId %>% startsWith(., gn))])
  toPlot$MI2 %<>% 
    c(., NA)
  toPlot$coverage %<>% 
    c(., mitIndex$sgCov[which(mitIndex$sgId == "non-targeting")])
  
  ok <- plyr::laply(toPlot, length) %>% 
    unique() %>% length() %>% equals(., 1)
  
  if (!ok) break
}
toPlot %<>% bind_cols()

p <- ggplot(toPlot %>% subset(., is_greater_than(coverage, 1500)), 
            aes(x = MI1 *100,
                y = MI2 * 100,
            )) +
  geom_point(size = 0.5) + 
  theme_classic() + 
  theme(legend.position = "top", legend.direction = "horizontal",
        text = element_text(size = 5)) +
  geom_abline(slope = 1, intercept = 0, color = "red", 
              linetype = "dotted") +
  xlab("sgRNA #1") +
  ylab("sgRNA #2")

p

ggsave("Mitotic_index_in_mean_cov_gt_1500.pdf", width = 2.15, height = 2.15,
       p)

v <- mitIndex %>% 
  subset(., !is.na(sgId) & (sgId != "non-targeting")) #%>% 
v %<>% subset(sgCov > 1500)
v$p.value <- 
  apply(v[, c("mitotic", "sgCov")], 1, 
        function(x) binom.test(x[1], x[2], p = nullMI)$p.value)
v$FDR <- p.adjust(v$p.value, method = "BH")
v1 <- subset(v, FDR < 0.1)
toKeep <- table(v1 %>% subset(., sgCov > 1500) %$% #gene) %>% 
                  sgId %>% strsplit(., split = "_") %>%
                  plyr::laply(., magrittr::extract, i = 1)) %>%
  subset(., equals(., 2)) %>%
  names()

v2 <- subset(v1, 
             sgId %>% strsplit(., split = "_") %>%
               plyr::laply(., magrittr::extract, i = 1) %>%
               is_in(., toKeep))

subset(v2, mitIndex > nullMI) %$% #gene
  plyr::laply(sgId, function(x) strsplit(x, split = "_")) %>%
  plyr::laply(., function(x) magrittr::extract(x, i = 1))

#Reject SMC3 because of opposite direction of changes between guides
subset(v2, mitIndex < nullMI) %$% #gene %>%
  plyr::laply(sgId, function(x) strsplit(x, split = "_")) %>%
  plyr::laply(., function(x) magrittr::extract(x, i = 1)) %>%
  setdiff(., "SMC3") %>%
  unique() %>% paste0(., collapse = ",")

write.table(v, "Differential_analysis_mitotic_index.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)



