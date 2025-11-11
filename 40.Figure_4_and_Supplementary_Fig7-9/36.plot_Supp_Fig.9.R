library(tidyverse)
library(magrittr)

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Documents/Ronin_chapters/projects/mcmanus_lab_chapters/Analysis/KC018_MCF10A-203_20250417/20.phenotype/dbl_pert_validations_integrated/b_arrayed_minipool_dbls_summary")
library(magrittr)
library(tidyverse)
library(ggrepel)

load("comboAnalysisWithNullDists.RData")

toPlot <- 
  data.frame(comboScore = c(cst_gene_level_totalH2AX$comboScore,
                            nullResAll$comboScore),
             set = c(rep("signal", 
                         nrow(cst_gene_level_totalH2AX)),
                     rep("null", nrow(nullResAll))))

toPlot$set %<>% 
  replace(., equals(., "null"), "null set") %>%
  replace(., equals(., "signal"), "tested set")

normFactor <- 
  subset(res_gene_level_totalH2AX, 
         (gene1 %>% startsWith(., "H2AFX_2_1")) & 
           (gene2 == "NTC_1")) %$% abs(estimate)
p <- ggplot(toPlot, aes(x = comboScore/normFactor, color = set)) +
  geom_density(adjust = 1.5) +
  theme_classic() +
  scale_color_manual(values = c("gray", "black")) +
  geom_vline(xintercept = median(nullResAll$comboScore/normFactor),
             color = "gray", linetype = "dotted") +
  # geom_vline(xintercept = median(cst_gene_level_totalH2AX$comboScore),
  #            color = "black", linetype = "dotted") +
  geom_point(data = filteredCstH2AX, y = -0.1, 
             mapping = aes(x = comboScore/normFactor), 
             inherit.aes = FALSE,
             size = 0.5, shape = 4) +
  coord_cartesian(ylim = c(-0.1, 2.5)) +#0.15)) +
  # guides(color = "none") + 
  theme(text = element_text(size = 7),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = "top",
        legend.key.size = unit(0.2, "cm"),
        legend.margin = margin()) +
  annotate(geom = "point", size = 0.5, shape = 4,
           x = -18/normFactor, y = 2.5) + #0.15) +
  annotate(geom = "text", size = 2, label = "hits",
           x = -16/normFactor, y = 2.5) + #0.15) +
  xlab("interaction score") +
  ylab("density")
p

ggsave("a_interaction_score_distributions_arrayed_minipools.pdf", 
       width = 3, height = 2)

# toPlot <- 
#   data.frame(rawEstComboScore = 
#                c(cst_gene_level_totalH2AX$rawEstComboScore,
#                  nullResAll$rawEstComboScore),
#              set = c(rep("signal", 
#                          nrow(cst_gene_level_totalH2AX)),
#                      rep("null", nrow(nullResAll))))
# 
# ggplot(toPlot, aes(x = rawEstComboScore, color = set)) +
#   geom_density(adjust = 1.5) +
#   theme_classic() +
#   geom_vline(xintercept = median(nullResAll$rawEstComboScore),
#              color = "#f8766d", linetype = "dotted") +
#   geom_vline(xintercept = 
#                median(cst_gene_level_totalH2AX$rawEstComboScore),
#              color = "#00bfc4", linetype = "dotted") +
#   geom_point(data = filteredCstH2AX, y = -0.07, 
#              mapping = aes(x = rawEstComboScore), 
#              inherit.aes = FALSE,
#              size = 1) +
#   coord_cartesian(ylim = c(-0.07, 4))

#-------------------
#Load validation data
# load("../comboScore_arr_individuals.RData")
# comboScores %<>% subset(., is.na(rawEstCombo) %>% not())
# toPlot <- data.frame(
#   rawEstComboScore = 
#     c(cst_gene_level_totalH2AX$rawEstComboScore,
#       nullResAll$rawEstComboScore,
#       comboScores$rawEstCombo),
#   set = c(rep("arrMiniS", 
#               nrow(cst_gene_level_totalH2AX)),
#           rep("arrMiniB", nrow(nullResAll)),
#           rep("arrIndS", nrow(comboScores))))
# 
# ggplot(toPlot, aes(x = rawEstComboScore, y = after_stat(scaled), 
#                    color = set)) +
#   geom_density(adjust = 1.5) +
#   theme_classic() +
#   geom_vline(xintercept = median(nullResAll$rawEstComboScore),
#              color = "#f8766d", linetype = "dotted") +
#   geom_vline(xintercept = 
#                median(cst_gene_level_totalH2AX$rawEstComboScore),
#              color = "#00bfc4", linetype = "dotted") +
#   geom_point(data = filteredCstH2AX, y = -0.07, 
#              mapping = aes(x = rawEstComboScore), 
#              inherit.aes = FALSE,
#              size = 1) +
#   coord_cartesian(ylim = c(-0.07, 1)) 

comboScores %<>% 
  subset(., geneB %>% is_in(c("BRD7", "CASC5", "KIFC1")) %>% not())
toPlot <- data.frame(
  comboScore = 
    c(cst_gene_level_totalH2AX$comboScore,
      nullResAll$comboScore,
      comboScores$estimateCombo - mean(arrIndNullRes$comboScore),
      arrIndNullRes$comboScore - mean(arrIndNullRes$comboScore)),
  set = c(rep("arrMiniS", 
              nrow(cst_gene_level_totalH2AX)),
          rep("arrMiniB", nrow(nullResAll)),
          rep("arrIndS", nrow(comboScores)),
          rep("arrIndB", nrow(arrIndNullRes))))
toPlot$expSetting <- substr(toPlot$set, 1, 
                            nchar(toPlot$set) - 1)
toPlot$type <- substr(toPlot$set, 
                      nchar(toPlot$set),
                      nchar(toPlot$set)) %>%
  replace(., equals(., "B"), "null set") %>%
  replace(., equals(., "S"), "interaction hits") %>%
  factor(., levels = c("null set", "interaction hits"))
# 
# ggplot(toPlot, aes(x = comboScore, y = after_stat(scaled), 
#                    color = set)) +
#   geom_density(adjust = 1.5) +
#   theme_classic() +
#   geom_vline(xintercept = median(nullResAll$comboScore),
#              color = "#f8766d", linetype = "dotted") +
#   geom_vline(xintercept = 
#                median(cst_gene_level_totalH2AX$comboScore),
#              color = "#00bfc4", linetype = "dotted") +
#   geom_point(data = filteredCstH2AX, y = -0.07, 
#              mapping = aes(x = comboScore), 
#              inherit.aes = FALSE,
#              size = 1) +
#   coord_cartesian(ylim = c(-0.07, 1))

toPlot <- subset(toPlot, expSetting == "arrInd")

normFactor2 <- res_gene_level_wrt_untransduced_totalH2AX %>% 
  subset(gene %>% startsWith(., "NA__H2AFX")) %$% 
  mean(estimate) %>% abs()
p <- ggplot(toPlot, aes(y = comboScore/normFactor2, 
                   x = type)) +
  geom_boxplot() +
  theme_classic() +
  ylab("interaction score") +
  geom_hline(yintercept = 0, linetype = "dotted", 
             linewidth = 0.5) +
  xlab(NULL) + theme(text = element_text(size = 7),
                     plot.margin = margin(0, 0, 0, 0, "cm"),)
  # facet_grid(expSetting ~ ., scales = "free_y")

ggsave("b_interaction_scores_arrayed_individuals.pdf", 
       width = 2, height = 2)

#------------------------------
toPlot <- subset(res_gene_level_wrt_untransduced_totalH2AX,
                 startsWith(gene, "NA__H2AFX_set") #| 
                   # startsWith(gene, "NTC__NTC_set")
                 )
toPlot$batch <- "batch #1"
toPlot$batch[toPlot$gene %>% endsWith(., "setB")] <- "batch #2"

p <- ggplot(toPlot, aes(x = batch, y = estimate)) +
  geom_col(width = 0.5) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab("H2AFX targeting control") +
  theme(text = element_text(size = 7),
        plot.margin = margin(0, 0, 0, 0, "cm"))

ggsave("c_h2ax_controls_barplot.pdf", 
       width = 2, height = 2)

#-----------------------------


#-----------------------------
# rm(list= ls())
# load("../phenoDat.rds")
# phenoDat$treatment %<>% 
#   str_replace_all(., "__NTC_2", "__NTC") %>% 
#   str_replace_all(., "__NTC_3", "__NTC") %>% 
#   str_replace_all(., "__NTC_5", "__NTC")
# phenoDat$treatment[phenoDat$treatment %>% endsWith(., "__NTC")] <-
#   paste0(phenoDat$treatment[phenoDat$treatment %>% endsWith(., "__NTC")],
#          "_", phenoDat$set[phenoDat$treatment %>% endsWith(., "__NTC")])
# 
# phenoDat$treatment[phenoDat$treatment %>% endsWith(., "__H2AFX")] <-
#   paste0(phenoDat$treatment[phenoDat$treatment %>% endsWith(., "__H2AFX")],
#          "_", phenoDat$set[phenoDat$treatment %>% endsWith(., "__H2AFX")])
# 
# phenoDat$treatment[phenoDat$treatment %>% endsWith(., "__MKI67_1")] %<>%
#   substr(., 1, nchar(.) - 2)
# phenoDat$treatment[phenoDat$treatment %>% endsWith(., "__MKI67_2")] %<>%
#   substr(., 1, nchar(.) - 2)
# # phenoDat$treatment %<>% replace(., startsWith(., "NTC__NTC"), "NTC__NTC")
# 
# #Removing sample that was accidentally lost except for few cells
# phenoDat %<>% subset(., !(grepl("ZAK", treatment) & 
#                             (plate == "Pb")))
# testSet <- unique(phenoDat$treatment)
# 
# toPlot <- subset(phenoDat, treatment %in%
#                    c("NA__H2AFX_setA", "NA__H2AFX_setB",
#                      "NTC__NTC_setA", "NTC__NTC_setB"))
# toPlot$treatment %<>%
#   str_remove_all(., "_setA") %>%
#   str_remove_all(., "_setB")
# 
# toPlot$sample <- "knockdown"
# for (nxSet in c("setA", "setB")) {
#   medFact <- subset(toPlot, (set == nxSet) &
#                       startsWith(treatment, "NTC__NTC")) %$%
#     median(totalH2AX, na.rm = TRUE)
#   
#   madFact <- subset(toPlot, 
#                     (set == nxSet) &
#                       startsWith(treatment, "NTC__NTC")) %$%
#     mad(totalH2AX, na.rm = TRUE)
#   
#   toPlot$totalH2AX[(toPlot$set == nxSet)] <- 
#     (toPlot$totalH2AX[(toPlot$set == nxSet)] - medFact)/madFact
# }
# 
# toPlot$sample <- "knockdown"
# toPlot$sample[toPlot$treatment == "NTC__NTC"] <- "non-\ntargeting"
# toPlot$sample %<>% factor(., levels = c("non-\ntargeting", "knockdown"))
# 
# toPlot$set %<>% replace(., equals(., "setA"), "batch #1")
# toPlot$set %<>% replace(., equals(., "setB"), "batch #2")
# 
# p <- ggplot(toPlot, aes(x = sample, y = totalH2AX)) +
#   geom_boxplot(outlier.shape = NA) + theme_classic() +
#   coord_cartesian(ylim = c(-3.5, 3.5)) +
#   facet_grid(. ~ set) +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   xlab(NULL) +
#   theme(text = element_text(size = 7),
#         plot.margin = margin(0, 0, 0, 0, "cm"))
# ggsave("c_h2ax_control.pdf", width = 2, height = 2)
# 


