#Code to generate supplemental figures

#load libraries
library(tidyverse)
library(ctmm)
library(lubridate)
library(pals)
library(grid)
library(gridExtra)
library(cowplot)

bearpairs <- read_csv("bearpairsto500m.csv") #when given permission, upload full bearpairs file with meta info
bears_meta <- read_csv("/meta/file/here") #need author permission to post meta file
uds_btwninds <- read_csv("uds_btwninds.csv")
uds_withininds <- read_csv("uds_withininds.csv")

#Fig. S1 (overlap between homeranges of individuals in pair, per pair, colored by subset)
btwnind_vals <- bearpairs %>%
  filter(est < 100) %>%
  select(c("bear1_id", "bear2_id", "time", "sex", "season", "season_breeding")) %>% 
  distinct(bear1_id, bear2_id, time, .keep_all=T) %>%
  left_join(uds_btwninds, by = c("bear1_id"="bear1_id", "bear2_id"="bear2_id", "time"="time")) %>%
  mutate(pair = paste(bear1_id, " ", bear2_id, "\n", year(time), sep="")) %>%
  pivot_longer(c("uds_before_low", "uds_after_low", "uds_before_mean", "uds_after_mean", "uds_before_high", "uds_after_high"), 
               names_pattern = "^((?:[^_]*_){1}[^_]*)_(.*)$", 
               names_to = c("encounter", "CI")) %>%
  pivot_wider(names_from = CI, values_from = value) %>%
  mutate(encounter=sub("^[^_]*_", "", encounter),
         encounter = factor(encounter, levels = c("before","after"))) %>%
  unite(c("sex", "season_breeding"), col="order", sep = " ", remove=FALSE) %>%
  mutate(order = fct_recode(order, M_M = "M_M hunting", M_M = "M_M breeding", 
                            F_F = "F_F hunting",  F_F = "F_F breeding", F_F = "F_F other",
                            `F_M Spring/Summer` = "F_M breeding",  
                            `F_M Spring/Summer` = "F_M other", 
                            `F_M Late Fall` = "F_M hunting"),
         order = fct_relevel(order, c("M_M", "F_F", "F_M Spring/Summer", "F_M Late Fall"))) %>%
  mutate(lab_color = ifelse(order=="M_M", "#ff7678", 
                            ifelse(order=="F_F", "#d9d2e9", 
                                   ifelse(order=="F_M Spring/Summer", "#cbfcbd", "#fedcba")))) %>%
  arrange(match(order, levels(order))) %>%
  mutate(pair = factor(pair, levels = c(unique(pair))))

ind_plots <- list()
p <- NULL

for (i in seq(1, length(btwnind_vals$bear1_id), 2)) {
  
  #red = M_M, purple = F_F, green = F_M, yellow = F_M hunting
  order <- btwnind_vals %>%
    slice(i) %>%
    select(order)
  
  lab_color <- ifelse(order=="M_M", "#ff7678", 
                      ifelse(order=="F_F", "#d9d2e9", 
                             ifelse(order=="F_M Spring/Summer", "#cbfcbd", "#fedcba")))
  
  p <- btwnind_vals %>%
    slice(i:(i+1)) %>%
    ggplot(mapping=aes(x=encounter, y = mean, color=encounter)) +
    geom_point() +
    geom_errorbar(mapping=aes(ymin=low,ymax=high)) +
    scale_color_manual(values=c("purple", "orange"), 
                       labels = c("Before Encounter", "After Encounter")) +
    ylim(0,1)+
    labs(x=NULL, y=NULL, color = NULL) +
    facet_wrap(~pair) +
    theme_bw() + #get rid of background
    theme(panel.grid.major = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth=0.5),
          text=element_text(size=20), strip.text = element_text(size=11), 
          axis.text.y = element_text(size=10),
          strip.background = element_rect(fill=lab_color),
          legend.position = "none",
          plot.background = element_rect(fill='transparent', color=NA))
  
  ind_plots <- c(ind_plots, list(p))
  
}

legend <- btwnind_vals %>%
  ggplot(aes(pair, mean, color = encounter, fill = order)) +
  geom_col(size=10) +
  labs(color=NULL, fill=NULL) +
  scale_fill_manual(values=c("#ff7678", "#d9d2e9", "#cbfcbd", "#fedcba"),
                    labels=c("MM", "FF", "FM Spring/Summer", "FM Late Fall")) +
  scale_color_manual(values=c("purple", "orange"),
                     labels = c("Before Encounter", "After Encounter")) +
  guides(fill=guide_legend(nrow = 1, byrow = T, order =2,
                           override.aes=list(size = 15, shape = 15, color = NULL)),
         color=guide_legend(nrow = 1, byrow = T, order=1,
                            override.aes=list(size = 15, shape = 16, fill=c("purple", "orange")))) +
  theme(legend.background=element_blank(), legend.position = "right", legend.justification="left")

#plot and save
png(filename = "bears_udsbtwn.png", res=95, width=1840, height = 840, units="px")

grid.draw(grid.arrange(arrangeGrob(grobs = c(ind_plots, 
                                             list(as_grob(get_legend(legend)))), #convert legend to grob
                                   ncol=10, left=textGrob("Range Distribution Overlap (Before vs. After)", rot=90, 
                                                          gp = gpar(col = "black", fontsize = 20)))))

dev.off()

#################
#Fig. S2 (overlap within individual homerange before/after encounter, per pair, colored by subset and sex)
bears_ind_sex <- bears_meta %>%
  select(c("Capture ID", "Bear #", "Sex")) %>%
  rename("ind_sex"=Sex)

ind_vals <- bearpairs %>%
  filter(est < 100) %>%
  select(c("bear1_id", "bear2_id", "time", "sex", "season", "season_breeding")) %>%
  distinct(bear1_id, bear2_id, time, .keep_all=T) %>%
  left_join(uds_withininds, by = c("bear1_id"="bear1_id", "bear2_id"="bear2_id", "time"="time")) %>%
  mutate(pair = paste(bear1_id, bear2_id, year(time), sep=" ")) %>%
  pivot_longer(c("uds_ind1_low", "uds_ind2_low", "uds_ind1_mean", "uds_ind2_mean", "uds_ind1_high", "uds_ind2_high"), 
               names_pattern = "^((?:[^_]*_){1}[^_]*)_(.*)$", 
               names_to = c("who", "CI")) %>%
  pivot_wider(names_from = CI, values_from = value) %>%
  mutate(who=sub("^[^_]*_", "", who),
         who = factor(who, levels = c("ind1","ind2"))) %>%
  unite(c("sex", "season_breeding"), col="order", sep = " ", remove=FALSE) %>%
  mutate(order = fct_recode(order, M_M = "M_M hunting", M_M = "M_M breeding", F_F = "F_F hunting", 
                            F_F = "F_F breeding", F_F = "F_F other",
                            `F_M Spring/Summer` = "F_M breeding",  
                            `F_M Spring/Summer` = "F_M other", 
                            `F_M Late Fall` = "F_M hunting"),
         order = fct_relevel(order, c("M_M", "F_F", "F_M Spring/Summer", "F_M Late Fall"))) %>%
  arrange(match(order, levels(order)))

ind_plots <- list()
p <- NULL
for (i in seq(1, length(ind_vals$bear1_id), 2)) {
  
  #red = M_M, purple = F_F, green = F_M, yellow = F_M hunting
  order <- ind_vals %>%
    slice(i) %>%
    select(order)
  
  lab_color <- ifelse(order=="M_M", "#ff7678", 
                      ifelse(order=="F_F", "#d9d2e9", 
                             ifelse(order=="F_M Spring/Summer", "#cbfcbd", "#fedcba")))
  
  
  p <- ind_vals %>%
    slice(i:(i+1)) %>%
    mutate(who = ifelse(#who=="ind1", pair1, pair2),
      who=="ind1", `bear1_id`, `bear2_id`),
      header = year(time)) %>%
    left_join(bears_ind_sex, by = c("who"="Bear #")) %>%
    ggplot(mapping=aes(x=who, y = mean, color=ind_sex)) +
    geom_point() +
    geom_errorbar(mapping=aes(ymin=low, ymax=high)) +
    ylim(0,1)+
    scale_color_manual(values=c("M" = "red", "F" = "blue")) +
    labs(x=NULL, y=NULL) +
    facet_grid(. ~ header) +
    theme_bw() + #get rid of background
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth=0.5),
          strip.background = element_rect(fill=lab_color), 
          text=element_text(size=17), strip.text = element_text(size=15), 
          axis.text.x = element_text(size = 10), axis.text.y= element_text(size=10),
          legend.position = "none")
  
  ind_plots <- c(ind_plots, list(p))
  
}

legend <- ind_vals %>%
  ggplot(aes(pair, mean, color = who, fill = order)) + #color is literally just mapped to a var with two states
  geom_col(size=10) +
  labs(color=NULL, fill=NULL) +
  scale_fill_manual(values=c("#ff7678", "#d9d2e9", "#cbfcbd", "#fedcba"),
                    labels=c("MM", "FF", "FM Spring/Summer", "FM Late Fall")) +
  scale_color_manual(values=c("red", "blue"),
                     labels = c("Male", "Female"),
                     drop=TRUE) +
  guides(fill=guide_legend(nrow = 1, byrow = T, order=2,
                           override.aes=list(size = 15, shape = 15, color = NULL)),
         color=guide_legend(nrow = 1, byrow = T, order=1,
                            override.aes=list(size = 15, shape = 16, fill=c("red", "blue")))) +
  theme(legend.background=element_blank(), legend.position = "right", legend.justification="left")

#plot/save
png(filename = "bears_udswin.png", res=95, width=1840, height = 840, units="px")

grid.draw(grid.arrange(arrangeGrob(grobs = c(ind_plots, 
                                             list(cowplot::as_grob(get_legend(legend)))), #convert legend to grob
                                   ncol=10, left=textGrob("Range Distribution Overlap Within Individuals", rot=90, 
                                                          gp = gpar(col = "black", fontsize = 20)))))

dev.off()

#Fig. S3 & S4 - need data owner permission to post location of carcass pits

