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

#Fig. S3 & S4 
#includes both code to calculate carcass pit distances and to make figures
#posted code but need data owner permission to post actual location of carcass pits

#get location of pits and date ranges active
carcasspit_locs <- read_csv("carcasspit/file/here") %>%
  mutate(`date active`=replace_na(strptime(`date active`, format="%d-%b-%y"), 
                                  as.Date("2000-01-01"))) %>%
  mutate(`date inactive`=replace_na(strptime(`date inactive`, format="%d-%b-%y"),
                                    as.Date("2023-04-01"))) #time ranges confirmed by data owner 

#get lat long coordinates for encounter locations (to compare to carcass pit locations)
encounterlocations <- data.frame("bear1"=bearpairs$pair1,
                                 "bear2"=bearpairs$pair2,
                                 "time"=bearpairs$time, "dist"=bearpairs$est,
                                 "season_breeding"=bearpairs$season_breeding,
                                 "sex"=bearpairs$sex,
                                 "midpt_lat"=NA, "midpt_long"=NA,
                                 "bear1_lat"=NA, "bear1_long"=NA,
                                 "bear2_lat"=NA, "bear2_long"=NA,
                                 "dist_to_pit"=NA, "nearest_pit"=NA) %>%
  mutate(threshold = plyr::round_any(dist, 100, f=ceiling)) %>%
  mutate(threshold = as.factor(ifelse(dist <= 50, 50, threshold))) %>%
  relocate(threshold, .after=dist)

#calculate distances between encounter and carcass pit
for (i in seq_along(encounterlocations$bear1)) {
  
  dat <- bears[c(which(names==encounterlocations$bear1[[i]]),
                 which(names==encounterlocations$bear2[[i]]))]
  
  #used already saved fits; if re-running fits then change this code
  fits <- list()
  for(i in seq_along(dat)) {
    guess <- ctmm.guess(dat[[i]],interactive=FALSE)
    fits[[i]] <- ctmm.select(dat[[i]],GUESS) #all default ctmm model fit settings, this might take a while to run
  }
  
  #predict tracks and get lat/long of continuous tracks for individuals 
  #(this is just a way to double check the encounter location)
  bear1_sim <- predict(dat[[1]],fits[[1]],t=ymd_hms(encounterlocations$time[[i]]), complete=T)
  bear2_sim <- predict(dat[[2]],fits[[2]],t=ymd_hms(encounterlocations$time[[i]]), complete=T)
  
  encounterlocations$bear1_lat[[i]] <- bear1_sim$latitude
  encounterlocations$bear1_long[[i]] <- bear1_sim$longitude
  
  encounterlocations$bear2_lat[[i]] <- bear2_sim$latitude
  encounterlocations$bear2_long[[i]] <- bear2_sim$longitude
  
  projection(dat[[1]]) <- median(dat)
  projection(dat[[2]]) <- median(dat)
  
  #get "midpoint" of bears at encounter time, i.e locaiton of encounter 
  #(they should be essentially on top of each other)
  encounter <- midpoint(dat, fits, t=ymd_hms(encounterlocations$time[[i]]), complete=T)
  
  encounterlocations$midpt_lat[i] <- encounter$latitude
  encounterlocations$midpt_long[i] <- encounter$longitude
  
  dists <- c()
  pitnames <- c()
  
  for (j in seq_along(carcasspit_locs$pit_name)) {
    
    #check to see if pit was active when encounter happened
    if(as.Date(encounterlocations$time[[i]]) %within% 
       interval(start=carcasspit_locs$`date active`[[j]],
                end=carcasspit_locs$`date inactive`[[j]])) {
      
      #get carcass pit as a telemetry object 
      carcasspit_tel <- as.telemetry(data.frame(location.lat=carcasspit_locs$lat[j],
                                                location.long=carcasspit_locs$long[j],
                                                timestamp=encounterlocations$time[i],
                                                individual.local.identifier="pit",
                                                individual.taxon.canonical.name=NA))
      
      #align projections in ctmm
      projection(carcasspit_tel) = projection(encounter)
      
      #get distance
      dists <- c(dists,
                 ctmm::distance(list(encounter, carcasspit_tel),method="Euclidean",sqrt=T)$CI[6])
      
      pitnames <- c(pitnames, carcasspit_locs$pit_name[j])
      
    } 
    
  }
  
  names(dists) <- pitnames
  
  encounterlocations$dist_to_pit[i] <- min(dists)
  encounterlocations$nearest_pit[i] <- names(which.min(dists))
  
  #with regular encounter tel object
  encounter_tel <- as.telemetry(data.frame(location.lat=encounter$latitude,
                                           location.long=encounter$longitude,
                                           timestamp=encounterlocations$time[i],
                                           individual.local.identifier="encounter",
                                           individual.taxon.canonical.name=NA))
  
  print(paste("done", encounterlocations$bear1[i], encounterlocations$bear2[i]))
  
}

#calculate distance from pits at all gps fixes that are NOT encounters (takes a very long time to run)

#set up df
nonencounterlocations <- data.frame("bear"=NA, "time"=NA, "dist"=NA, "season_breeding"=NA, 
                                    "sex"=NA, "dist_from_pit_1"=NA) %>%
  drop_na()

#for each encounter
for (i in seq_along(bearpairs$pair1)) {
  
  bear1 <- bears[[which(names==bearpairs$pair1[[i]])]]
  bear2 <- bears[[which(names==bearpairs$pair2[[i]])]]
  
  #collapse data to day to make it more computationally feasible and 
  #remove day of encounter from data (figured whole day was ok)
  bear1day <- data.frame(bear1) %>%
    mutate(timestamp = as_date(timestamp)) %>%
    distinct(timestamp, .keep_all=T) %>%
    filter(timestamp != as_date(bearpairs$time[[i]])) %>%
    as.telemetry()
  
  bear2day <- data.frame(bear2) %>%
    mutate(timestamp = as_date(timestamp)) %>%
    distinct(timestamp, .keep_all=T) %>%
    filter(timestamp != as_date(bearpairs$time[[i]])) %>%
    as.telemetry()
  
  allpits_dist1 <- data.frame()
  allpits_dist2 <- data.frame()
  
  for (j in seq_along(carcasspit_locs$pit_name)) {
    
    #if the carcass pit is active during the specified date calculate the distance
    if(as.Date(encounterlocations$time[[i]]) %within% 
       interval(start=carcasspit_locs$`date active`[[j]],
                end=carcasspit_locs$`date inactive`[[j]])) {
      
      #distance calculated here with geosphere::distHaversine,
      #calculations work the same with telemetry obj but it's more efficient with geosphere bc there are so many more 
      dists_1 <- unlist(lapply(seq(nrow(bear1)),function(x){
        pt_bear <- unlist(bear1[x, c("longitude", "latitude")])
        pt_pit <- unlist(carcasspit_locs[j, c("long", "lat")])
        geosphere::distHaversine(pt_bear, pt_pit)}))
    
      dists_2 <- unlist(lapply(seq(nrow(bear2)),function(x){
        pt_bear <- unlist(bear2[x, c("longitude", "latitude")])
        pt_pit <- unlist(carcasspit_locs[j, c("long", "lat")])
        geosphere::distHaversine(pt_bear, pt_pit)}))
      
    }
    
    #else skip it
    else {
      
      dists_1 <- rep(NA, length(bear1$timestamp))
      dists_2 <- rep(NA, length(bear2$timestamp))
    }
    
    allpits_dist1 <- rbind(allpits_dist1, dists_1)
    allpits_dist2 <- rbind(allpits_dist2, dists_2)
    
  }
  
  #returns distance from each pit at the specified timestamp
  allpits_dist1 <- allpits_dist1 %>% 
    `rownames<-`(carcasspit_locs$pit_name) %>% 
    t() %>% 
    data.frame() %>% 
    `rownames<-`(NULL)
  
  allpits_dist2 <- allpits_dist2 %>% 
    `rownames<-`(carcasspit_locs$pit_name) %>% 
    t() %>% 
    data.frame() %>% 
    `rownames<-`(NULL)
  
  #get min dist to ANY pit at each timestamp
  allpits_mindist1 <- apply(allpits_dist1, 1, min, na.rm=T)
  allpits_mindist2 <- apply(allpits_dist2, 1, min, na.rm=T)
  
  #bind fix for each bear
  nonencounterlocations <- rbind(nonencounterlocations,
                                 cbind("bear"=bearpairs$pair1[[i]], 
                                       "time"=bearpairs$time[[i]], 
                                       "dist"=bearpairs$est[[i]], 
                                       "season_breeding"=bearpairs$season_breeding[[i]], 
                                       "sex"=bearpairs$sex[[i]], 
                                       "dist_from_pit"=allpits_mindist1)) %>%
    rbind(nonencounterlocations,
          cbind("bear"=bearpairs$pair2[[i]], 
                "time"=bearpairs$time[[i]], 
                "dist"=bearpairs$est[[i]], 
                "season_breeding"=bearpairs$season_breeding[[i]], 
                "sex"=bearpairs$sex[[i]], 
                "dist_from_pit"=allpits_mindist2))
  
  print(paste("done", bearpairs$pair1[[i]], bearpairs$pair2[[i]]))
  
}

#Make figures!

#pts removed for being >20000m away from any carcass pit
n_removed <- (nonencounterlocations %>%
                filter(dist_from_pit > 50000) %>%
                count() %>%
                summarize(n=n/nrow(nonencounterlocations)))$n

#histograms/density plots of dis from carcass pit
n_removed_fall <- (nonencounterlocations %>%
                     filter(season_breeding=="hunting") %>%
                     filter(dist_from_pit > 50000) %>%
                     count() %>%
                     summarize(n=n/nrow(filter(nonencounterlocations, season_breeding=="hunting"))))$n

n_removed_spring <- (nonencounterlocations %>%
                       filter(season_breeding=="breeding") %>%
                       filter(dist_from_pit > 50000) %>%
                       count() %>%
                       summarize(n=n/nrow(filter(nonencounterlocations, season_breeding=="breeding"))))$n

n_removed_other <- (nonencounterlocations %>%
                      filter(season_breeding=="other") %>%
                      filter(dist_from_pit > 50000) %>%
                      count() %>%
                      summarize(n=n/nrow(filter(nonencounterlocations, season_breeding=="other"))))$n

#Figure 3. Distance from carcass pits based on sex
nonencounterlocations %>%
  filter(dist_from_pit < 50000) %>%
  mutate(sex=factor(sex, levels=c("F_M", "F_F", "M_M"))) %>%
  ggplot(mapping=aes(x=dist_from_pit, y=after_stat(density), color=sex)) +
  geom_density() +
  geom_histogram(data=mutate(encounterlocations, sex=factor(sex, levels=c("F_M", "F_F", "M_M"))),
                 mapping=aes(x=dist_to_pit, fill=sex),
                 alpha=0.6, size=0.1) +
  scale_color_manual(labels = c("FM", "FF", "MM"),
                     values=c("purple", "blue", "red")) +
  scale_fill_manual(labels = c("FM", "FF", "MM"),
                    values=c("purple", "blue", "red")) +
  facet_wrap(~sex, 
             labeller = labeller(sex = c(`F_M`="FM",  `F_F`="FF", `M_M`="MM"))) +
  labs(x="Distance from nearest carcass pit (m)", y="Density of fixes", 
       color="Non-encounter fixes", fill="Encounter fixes") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  guides(color = guide_legend(override.aes = list(fill = "white")))

#Figure 4: Distance from carcass pits based on season
nonencounterlocations %>%
  filter(dist_from_pit < 50000) %>%
  mutate(season_breeding=factor(season_breeding, 
                                levels=c("hunting", "breeding", "other"))) %>%
  ggplot(mapping=aes(x=dist_from_pit, y=after_stat(density), color=season_breeding)) +
  geom_density() +
  geom_histogram(data=mutate(encounterlocations,
                             season_breeding=factor(season_breeding, 
                                                    levels=c("hunting", "breeding", "other"))),
                 mapping=aes(x=dist_to_pit, fill=season_breeding),
                 alpha=0.6, size=0.1) +
  scale_color_manual(labels = c("Late Fall","Late Spring/\nSummer", "Other"),
                     values=c("#00BA38", "#F8766D", "#619CFF")) +
  scale_fill_manual(labels = c("Late Fall", "Late Spring/Summer", "Other"),
                    values=c("#00BA38", "#F8766D", "#619CFF")) +
  facet_wrap(~season_breeding, 
             labeller = labeller(season_breeding = c(hunting="Late Fall", 
                                                     breeding="Late Spring/\nSummer", 
                                                     other="Other"))) +
  labs(x="Distance from nearest carcass pit (m)", y="Density of fixes", 
       color="Non-encounter fixes", fill="Encounter fixes") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  guides(color = guide_legend(order=1,
                              override.aes = list(fill = "white"), size=1),
         fill = guide_legend(order = 2))


