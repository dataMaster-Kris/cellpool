library(tidyverse)
library(magrittr)
library(patchwork)

dat <- list.files() %>% 
  subset(., endsWith(., "csv")) %>% 
  map(., function(x) {
    read.csv(x) %>% 
      mutate(filename = x)
  }) %>% 
  bind_rows()

smpIds <- data.frame(
  filename = c("export_no-virus-anti-CD81-stained_Single Cells - SSC.csv", 
               "export_no-virus-unstained_Single Cells - SSC.csv", 
               "export_pV1TUd1_CD81ko_Cre_anti-CD81-stained_Single Cells - SSC.csv", 
               "export_pV2TUd7_CD81ko_Cre_anti-CD81-stained_Single Cells - SSC.csv", 
               "export_pV2TUd8_CD81ko_Cre_anti-CD81-stained_Single Cells - SSC.csv", 
               "export_pV2TUdCtrl2_NTC_Cre_anti-CD81-stained_Single Cells - SSC.csv", 
               "export_pVcTUd4_CD81ko_Cre_anti-CD81-stained_Single Cells - SSC.csv"),
  id = c("no-virus-anti-CD81-stained",
         "no-virus-unstained",
         "pV1TUd1_CD81ko",
         "pV2TUd7_CD81ko",
         "pV2TUd8_CD81ko",
         "pV2TUdCtrl2_NTC",
         "pVcTUd4_CD81ko")
)

dat <- merge(dat, smpIds)
sigma_gfp <- sd(subset(dat, id == "pV2TUdCtrl2_NTC")$BL1.A)
q_gfp <- quantile(subset(dat, id == "no-virus-anti-CD81-stained")$BL1.A, 0.98)
q_cd81 <- quantile(subset(dat, id == "no-virus-anti-CD81-stained")$RL1.A, 0.01)
sigma_cd81 <- sd(subset(dat, id == "no-virus-anti-CD81-stained")$RL1.A)

dat$gfp <- (dat$BL1.A - q_gfp)/sigma_gfp
dat$cd81 <- (dat$RL1.A - q_cd81)/sigma_cd81

#GFP+_fraction is the fraction of CD81- cells in GFP+ cells.
pct_cd81KD_by_FLOWJO_gating <- 
  read.csv("FLOWJO_bivariate_CD81_vs_GFP_based_gating.csv") %>%
  set_colnames(., c("filename", "GFP+_fraction", "NA")) %>%
  subset(., !(filename %in% c("Mean", "SD")), 
         select = c("filename", "GFP+_fraction")) 

smpIds$filename %<>% str_remove_all(., "export_") %>% 
  str_remove_all(., "_Single Cells - SSC.csv")
pct_cd81KD_by_FLOWJO_gating$filename %<>% str_remove_all(., ".fcs")
pct_cd81KD_by_FLOWJO_gating %<>% merge(., smpIds)
pct_cd81KD_by_FLOWJO_gating$`GFP+_fraction` %<>% 
  sprintf("%.1f", .) %>%
  paste0(., "%")

GFP_ePDF <- list()
for (id in unique(dat$id)) {
  thisHist <- hist(dat$gfp[dat$id == id], plot = FALSE, nclass = 30)
  GFP_ePDF[[id]] <- data.frame(id = id, 
                               density = thisHist$density,
                               mid = thisHist$mid)
  GFP_ePDF[[id]]$density <- GFP_ePDF[[id]]$density/max(GFP_ePDF[[id]]$density)
}
GFP_ePDF %<>% bind_rows()

GFP_ePDF$id %<>%
  factor(., 
         levels = c("no-virus-unstained",
                    "no-virus-anti-CD81-stained",
                    "pV2TUdCtrl2_NTC",
                    "pVcTUd4_CD81ko",
                    "pV1TUd1_CD81ko",
                    "pV2TUd7_CD81ko",
                    "pV2TUd8_CD81ko"))

CD81_ePDF <- list()
for (id in unique(dat$id)) {
  thisHist <- hist(dat$cd81[(dat$id == id) & 
                              (dat$gfp > 0)], plot = FALSE, nclass = 50)
  CD81_ePDF[[id]] <- data.frame(id = id, 
                                density = thisHist$density,
                                mid = thisHist$mid)
  CD81_ePDF[[id]]$density <- CD81_ePDF[[id]]$density/max(CD81_ePDF[[id]]$density)
}
CD81_ePDF %<>% bind_rows()

CD81_ePDF$id %<>%
  factor(., 
         levels = c("no-virus-unstained",
                    "no-virus-anti-CD81-stained",
                    "pV2TUdCtrl2_NTC",
                    "pVcTUd4_CD81ko",
                    "pV1TUd1_CD81ko",
                    "pV2TUd7_CD81ko",
                    "pV2TUd8_CD81ko"))

p1 <- ggplot(GFP_ePDF %>% 
         subset(., is_in(id, c("pVcTUd4_CD81ko",
                               "pV1TUd1_CD81ko",
                               "pV2TUd7_CD81ko",
                               "pV2TUd8_CD81ko"))), aes(x = mid, y = density)) +
  annotate(x = GFP_ePDF[GFP_ePDF$id == "no-virus-unstained", 
                        "mid"], 
           y = GFP_ePDF[GFP_ePDF$id == "no-virus-unstained", 
                        "density"],
           geom = "line", color = "grey") +
  geom_smooth(se = FALSE, color = "black", span = 0.3, linewidth = 0.5) +
  facet_grid(id ~ .) +
  ylab(NULL) + xlab("GFP") + 
  coord_cartesian(xlim = c(-4, 4)) +
  theme_classic() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 5),
    strip.background = element_blank(),
    strip.text = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(-4, 0, 4)) +
  geom_hline(yintercept = -0.1)

p2 <- 
  ggplot(CD81_ePDF %>% 
               subset(., is_in(id, c("pVcTUd4_CD81ko",
                                     "pV1TUd1_CD81ko",
                                     "pV2TUd7_CD81ko",
                                     "pV2TUd8_CD81ko"))), 
             aes(x = mid, y = density)) +
  annotate(x = CD81_ePDF[CD81_ePDF$id == "pV2TUdCtrl2_NTC", 
                         "mid"], 
           y = CD81_ePDF[CD81_ePDF$id == "pV2TUdCtrl2_NTC", 
                         "density"],
           # linetype = "dotdash", 
           geom = "line", color = "grey") +
  geom_smooth(se = FALSE, color = "black", span = 0.1, linewidth = 0.5) +
  facet_grid(id ~ .) +
  ylab(NULL) + xlab("CD81") + 
  theme_classic() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 5),
    strip.background = element_blank(),
    strip.text = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(-4, 0, 4)) +
  geom_hline(yintercept = -0.1)

p1NoText <- (p1 + 
               theme(axis.text.x = element_blank(), 
                     axis.title.x = element_blank()))
p2NoText <- (p2 + 
               theme(axis.text.x = element_blank(), 
                     axis.title.x = element_blank()))

ggsave(paste0("p1_enAsCas12a_CD81ko_comparisons_", Sys.time(),".pdf"),
       plot = p1NoText,
       width = 0.6*2.5, height = 1.7*2.5)
ggsave(paste0("p2_enAsCas12a_CD81ko_comparisons_", Sys.time(),".pdf"),
       plot = p2NoText,
       width = 0.6*2.5, height = 1.7*2.5)
