library(tidyverse)
library(magrittr)
library(patchwork)

dat <- list.files() %>% 
  subset(., startsWith(., "FLOWJO") %>% not()) %>%
  subset(., endsWith(., "csv")) %>% 
  map(., function(x) {
    read.csv(x) %>% 
      mutate(filename = x)
  }) %>% 
  bind_rows()

smpIds <- data.frame(
  filename = c("export_no-virus-anti-CD81-stained_Single Cells - FSC.csv", 
               "export_pV2TUd12_CD81i_and_Cre_anti-CD81-stained_Single Cells - FSC.csv", 
               "export_pV2TUdCtrl2_NTC_anti-CD81-stained_Single Cells - FSC.csv"),
  id = c("no-virus-anti-CD81-stained",
         "pV2TUd12_CD81i",
         "pV2TUdCtrl2_NTC")
)

dat <- merge(dat, smpIds)
sigma_gfp <- sd(subset(dat, id == "pV2TUdCtrl2_NTC")$BL1.A)
q_gfp <- quantile(subset(dat, id == "no-virus-anti-CD81-stained")$BL1.A, 0.98)
mu_cd81 <- mean(subset(dat, id == "no-virus-anti-CD81-stained")$RL1.A)
q_cd81 <- quantile(subset(dat, id == "no-virus-anti-CD81-stained")$RL1.A, 0.01)
sigma_cd81 <- sd(subset(dat, id == "no-virus-anti-CD81-stained")$RL1.A)

dat$gfp <- (dat$BL1.A - q_gfp)/sigma_gfp
dat$cd81 <- (dat$RL1.A - q_cd81)/sigma_cd81

#GFP+_fraction is the fraction of CD81- cells in GFP+ cells.
#Overall_fraction is the fraction of CD81- cells in all single cells.
pct_cd81KD_by_FLOWJO_gating <- 
  read.csv("FLOWJO_bivariate_CD81_vs_GFP_based_gating_dCas9-KRAB.csv") %>%
  set_colnames(., c("filename", "Overall_fraction", "GFP+_fraction", "NA")) %>%
  subset(., !(filename %in% c("Mean", "SD", "no-virus-unstained.fcs",
                              "pV2TUd12_CD81i_no_Cre_anti-CD81-stained.fcs")), 
         select = c("filename", "Overall_fraction", "GFP+_fraction")) 

smpIds$filename %<>% str_remove_all(., "export_") %>% 
  str_remove_all(., "_Single Cells - FSC.csv")
pct_cd81KD_by_FLOWJO_gating$filename %<>% str_remove_all(., ".fcs")
pct_cd81KD_by_FLOWJO_gating %<>% merge(., smpIds)
#Using Overall_fraction for no-virus-anti-CD81-stained because very low count of GFP+ cells in this sample.
pct_cd81KD_by_FLOWJO_gating$`GFP+_fraction`[
  pct_cd81KD_by_FLOWJO_gating$id == "no-virus-anti-CD81-stained"] <-
  pct_cd81KD_by_FLOWJO_gating$Overall_fraction[
    pct_cd81KD_by_FLOWJO_gating$id == "no-virus-anti-CD81-stained"
  ]
pct_cd81KD_by_FLOWJO_gating$`GFP+_fraction` %<>% 
  sprintf("%.1f", .) %>%
  paste0(., "%")

GFP_ePDF <- list()
for (id in unique(dat$id)) {
  thisHist <- hist(dat$gfp[(dat$id == id) & 
                             ((dat$gfp > 0.05) | 
                                (dat$id != "pV2TUd12_CD81i"))], plot = FALSE, nclass = 100)
  GFP_ePDF[[id]] <- data.frame(id = id, 
                               density = thisHist$density,
                               mid = thisHist$mid)
  GFP_ePDF[[id]]$density <- GFP_ePDF[[id]]$density/max(GFP_ePDF[[id]]$density)
}
GFP_ePDF %<>% bind_rows()

GFP_ePDF$id %<>%
  factor(., 
         levels = c("no-virus-anti-CD81-stained",
                    "pV2TUdCtrl2_NTC",
                    "pV2TUd12_CD81i"))

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
         levels = c("no-virus-anti-CD81-stained",
                    "pV2TUdCtrl2_NTC",
                    "pV2TUd12_CD81i"))

p1 <- ggplot(GFP_ePDF %>% 
               subset(., id == "pV2TUd12_CD81i"),
             aes(x = mid, y = density)) +
  geom_smooth(se = FALSE, color = "black", span = 0.1, linewidth = 0.3) +
  facet_grid(id ~ .) +
  annotate(x = GFP_ePDF[(GFP_ePDF$id == "no-virus-anti-CD81-stained") &
                          (GFP_ePDF$density != 0), 
                        "mid"], 
           y = GFP_ePDF[(GFP_ePDF$id == "no-virus-anti-CD81-stained") &
                          (GFP_ePDF$density != 0), 
                        "density"],
           geom = "line", color = "grey", linewidth = 0.3) +
  ylab(NULL) + xlab("GFP") + 
  coord_cartesian(xlim = c(-0.5, 5.5)) +
  theme_classic() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    axis.line = element_line(linewidth = 0.3),
    axis.title.x = element_text(size = 7),
    axis.text = element_text(size = 5),
    strip.background = element_blank(),
    strip.text = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) 

p2 <- 
  ggplot(CD81_ePDF %>% 
           subset(., id == "pV2TUd12_CD81i"), 
         aes(x = mid, y = density)) +
  geom_line(linewidth = 0.3) +
  facet_grid(id ~ .) +
  annotate(x = CD81_ePDF[CD81_ePDF$id == "pV2TUdCtrl2_NTC", 
                         "mid"], 
           y = CD81_ePDF[CD81_ePDF$id == "pV2TUdCtrl2_NTC", 
                         "density"],
           geom = "line", color = "grey", linewidth = 0.3) +
  ylab(NULL) + xlab("CD81") + 
  coord_cartesian(xlim = c(-4, 8)) +
  theme_classic() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    axis.line = element_line(linewidth = 0.3),
    axis.title.x = element_text(size = 7),
    axis.text = element_text(size = 5),
    strip.background = element_blank(),
    strip.text = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(-2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 0.5, 1))

p1NoText <- (p1 + 
               theme(axis.text.x = element_blank(), 
                     axis.title.x = element_blank()))
p2NoText <- (p2 + 
               theme(axis.text.x = element_blank(), 
                     axis.title.x = element_blank()))

ggsave(paste0("p1_dCas9_CD81i_comparisons_", Sys.time(),".pdf"),
       plot = p1NoText, width = 2.18, height = 1.14)
ggsave(paste0("p2_dCas9_CD81i_comparisons_", Sys.time(),".pdf"),
       plot = p2NoText, width = 2.18, height = 1.14)

ggsave(paste0("p1_dCas9_CD81i_comparisons_", Sys.time(),".pdf"),
       plot = p1, width = 2.18, height = 1.14)
ggsave(paste0("p2_dCas9_CD81i_comparisons_", Sys.time(),".pdf"),
       plot = p2, width = 2.18, height = 1.14)
