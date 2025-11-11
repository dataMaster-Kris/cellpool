library(tidyverse)
library(magrittr)
library(biomaRt)
library(wordcloud2)
library(ggrepel)
library(ggwordcloud)
library(patchwork)

genes1 <- "MAD2L2, PTTG1, SMC3, MAPRE1, STAG1, CHAF1B, SNX33, MITD1, DMTF1, CCNA2, CHEK2, PAK4, LATS1, CIT, KIF4A, CDCA8, ATM, MAPRE2, NCAPH2, UHRF1, CDKN3, MAP3K8, ZWINT, CDK20, RPS6KB1, RACGAP1, CNTRL, SMC4, CDKN1B, NCAPG2, KIF2B, SMC1A, CCNB1, HAUS8, CENPE, ORC1, NSUN2, SON, MASTL, MAD1L1, MAPK13, KIF18B, PLK3, EIF2AK4, AJUBA, MAPK6, CDK2AP1, TACC1, STAG2, CENPF, ANLN, CDK6, PCNT, FOXM1, CDKN2A, E2F1, BRD7, ANAPC5, MKI67, NDC80, WAPAL, CDC7, KIF2A, ZFYVE26, BUB1, RPS6KA3, PRC1, CENPC, SETDB2, BUB1B, NCAPD3, CETN2, MAP9, BIRC5, CEP55, CDK1, NUF2, NEK11, TERF2, ZAK, H2AFX, KIFC1, LATS2, SEPT6, AURKB, KIF20A, NUSAP1, CUL3, SUV39H1, PIK3C3, NCAPH, CDC20, DSN1, SEPT1, KIF23, ECT2, KIF13A, CDK7, BRCA1, PDS5A, PRKCE, PDS5B, CSNK2A1, MAD2L1, TP53, CHEK1, SNX9, SIRT2, ESCO2, STK10, CKAP2, PLK1, KIF11, NSL1, CDK5, INCENP, CSNK2A2, RIF1, RAD21, TAF1L, NCAPD2, CASC5, MAPK7, MYC, RBL1, KIF22, AURKC, NPM1, CDKN2C, TERF1, PDCD6IP, CYLD, BRCA2, CENPA, RB1, KIF2C, CETN1, KIF20B, BUB3, BORA, NEK1, PKMYT1, TTK, CDKN1C, NASP, USP8, CAMK1, SMARCB1, RCC2, CETN3, SMC2, CDK14, PIM2, GSG2, NEK9, BRSK1, ESPL1, CDKN2B, PRKCD, POGZ, TEX14, NEK3, MAPK4, NCAPG" %>%
  strsplit(., split = ", ") %>% 
  unlist()
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes2 <- getBM(
  attributes = c("name_1006", "go_id", "hgnc_symbol", "namespace_1003"
                 ),
  filters = "hgnc_symbol",
  values = genes1,
  mart = ensembl
) %>%
  subset(., name_1006 != "")

goTerm2Dom <- distinct(genes2[, c("name_1006", "namespace_1003")])
dfWord <- table(genes2$name_1006) %>% 
  as.data.frame() %>%
  subset(., Var1 != "") %>% 
  subset(., Freq > 19) %>%
  merge(., goTerm2Dom, by.x = "Var1", by.y = "name_1006")
dfWord$Var1 %<>% as.character()
dfWord <- subset(dfWord, 
                 !((Var1 == c("identical protein binding")) |
                     (startsWith(Var1, "histone") &
                        endsWith(Var1, "kinase activity"))))

dfWord$Var1 %<>% factor(levels = dfWord$Var1[order(dfWord$Freq)])
keepTerms <- c("protein binding", "ATP binding", "protein kinase activity",
               "centrosome", "kinetochore", "spindle", "microtubule binding",
               "cell division", #"protein phosphorylation", 
               "chromatin remodeling", "nucleolus", "apoptotic process",
               "DNA repair")
for (ix in unique(dfWord$namespace_1003)) {
  thisDf <- dfWord %>% subset(., (namespace_1003 == ix) &
                                as.character(Var1) %in% keepTerms)
  thisDf <- thisDf[order(thisDf$Freq), ]
  p <- ggplot(thisDf,
              aes(x = Freq, y = Var1)) +
    geom_col() +
    theme_classic() +
    ylab(NULL) +
    theme(axis.text = element_text(size = 5),
          plot.margin = margin(0,0,0,0))
  assign(paste0("p_", ix), p)
}

thisDf <- dfWord %>% subset(., 
                              as.character(Var1) %in% keepTerms)
thisDf <- thisDf[order(thisDf$Freq), ]
thisDf[thisDf$namespace_1003 == "biological_process", "namespace_1003"] <- 
  "biological\nprocess"
thisDf[thisDf$namespace_1003 == "cellular_component", "namespace_1003"] <- 
  "cellular\ncomponent"
thisDf[thisDf$namespace_1003 == "molecular_function", "namespace_1003"] <- 
  "molecular\nfunction"
p <- ggplot(thisDf,
            aes(x = Freq, y = Var1)) +
  geom_col() +
  theme_classic() +
  ylab(NULL) + xlab(NULL) +
  theme(axis.text = element_text(size = 5),
        strip.text = element_text(size = 5)) +
  facet_grid(namespace_1003 ~ ., scales = "free")

ggsave("fig3a.pdf", width = 1.7, height = 2)
