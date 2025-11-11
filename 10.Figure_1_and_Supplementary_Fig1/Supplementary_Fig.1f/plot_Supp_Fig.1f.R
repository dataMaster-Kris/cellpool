library(patchwork)
library(tidyverse)
library(magrittr)

dC9 <- read.csv("control_vectors_KD_efficiency_06-Mar-2025.wsp FlowJo table.csv") %>%
  subset(., startsWith(X, "pV2TUd"), 
         select = c("X", "Cells.Single.Cells...FSC.Single.Cells...SSC.GFP..CD81....Freq..of.Parent....")) %>%
  set_colnames(., c("smp", "freqNoCD81")) %>%
  mutate(cas = "dCas9-KRAB", smp = str_remove_all(smp, "_anti-CD81-stained.fcs") %>% 
           str_remove_all(., "pV2TUdCtrl"))
wtC9 <- read.csv("control_vectors_KO_efficiency_06-Mar-2025.wsp FlowJo table.csv") %>%
  subset(., startsWith(X, "pV2TUd"), 
         select = c("X", "Cells.Single.Cells...FSC.Single.Cells...SSC.GFP..CD81....Freq..of.Parent....")) %>%
  set_colnames(., c("smp", "freqNoCD81")) %>%
  mutate(cas = "Cas9", 
         smp = str_remove_all(smp, "_anti-CD81-stained.fcs") %>% 
           str_remove_all(., "pV2TUdCtrl"))

promoters <- data.frame(vec = c(2, 4, 5, 7, 9), 
                        desc = c("hU6", "mU6", "TATAlox-mU6",
                                 "TATAlox-hU6-1",
                                 "TATAlox-hU6-2"))

dC9$vec <- dC9$smp %>% strsplit(., split = "_") %>% 
  plyr::laply(., magrittr::extract, i = 1) %>% 
  as.numeric()
dC9$insert <- dC9$smp %>% strsplit(., split = "_") %>% 
  plyr::laply(., magrittr::extract, i = 2)

dC9 <- merge(dC9, promoters)
dC9$smpId <- apply(dC9[, c("desc", "insert")], 1, 
                   function(x) paste0(x[1], "_", x[2]))

dC9$smpId %<>% 
  factor(., levels = c(
    "hU6_CD81i", "mU6_CD81i", 
    "TATAlox-hU6-1_CD81i", "TATAlox-mU6_CD81i", 
    "TATAlox-hU6-2_CD81i", 
    "hU6_NTC", "mU6_NTC"))
dC9$freqNoCD81 %<>% divide_by(., 100)

wtC9$vec <- wtC9$smp %>% strsplit(., split = "_") %>% 
  plyr::laply(., magrittr::extract, i = 1) %>% 
  as.numeric()
wtC9$insert <- wtC9$smp %>% strsplit(., split = "_") %>% 
  plyr::laply(., magrittr::extract, i = 2)

wtC9 <- merge(wtC9, promoters)
wtC9$smpId <- apply(wtC9[, c("desc", "insert")], 1, 
                    function(x) paste0(x[1], "_", x[2]))

wtC9$smpId %<>% 
  factor(., levels = c(
    "hU6_CD81ko", "mU6_CD81ko", 
    "TATAlox-hU6-1_CD81ko", "TATAlox-mU6_CD81ko", 
    "TATAlox-hU6-2_CD81ko", 
    "hU6_NTC", "mU6_NTC"))
wtC9$freqNoCD81 %<>% divide_by(., 100)

p1 <- ggplot(dC9, 
             aes(x = smpId, y = freqNoCD81)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) + 
  xlab(NULL) + ylab("Fraction of cells with CD81\nknockdown by dCas9-KRAB")

p2 <- ggplot(wtC9, 
             aes(x = smpId, y = freqNoCD81)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) + 
  xlab(NULL) + ylab("Fraction of cells with CD81\nknockout by Cas9")

p <- (p2 + 
        theme(axis.text = element_text(size = 5), 
              axis.title = element_text(size = 5)) +
        scale_y_continuous(breaks = c(0, 0.5, 1))) | (p1 + 
  theme(axis.text = element_text(size = 5), 
        axis.title = element_text(size = 5)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)))
ggsave("wtCas9_and_dCas9_CD81_KD_efficiency.pdf", p, width = 3.26, height = 2)
