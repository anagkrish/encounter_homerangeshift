#Code to generate supplemental figures

#load libraries
library(tidyverse)
library(ctmm)
library(lubridate)
library(pals)
library(grid)
library(gridExtra)
library(cowplot)

#NOTE:: some of these figures require information on cubs, which has not been explicitly allowed to share by author -- will post once we have full permission!

bearpairs <- read_csv("bearpairsto500m.csv") %>% #when given permission, upload full bearpairs file with all meta info
  filter(pair1_ess_before > 2, pair2_ess_before > 2, pair1_ess_after > 2, pair2_ess_after > 2) #drop all encounters with subsets with ess <2
  
bears_meta <- read_csv("trunc_meta.csv") #truncated file with some meta info on bears, need author permission to post meta file

#uds overlap between individuals before vs after encounter
uds_btwninds <- read_csv("uds_btwninds.csv") %>%
  filter(pair1_ess_before > 2, pair2_ess_before > 2, pair1_ess_after > 2, pair2_ess_after > 2) #drop all encounters with subsets with ess <2

#uds overlap within individual before vs after encounter
uds_withininds <- read_csv("uds_withininds.csv") %>%
  filter(ind1_ess_before > 2, ind2_ess_before > 2, ind1_ess_after > 2, ind2_ess_after > 2) #drop all encounters with subsets with ess <2

##########################
#Fig. S1 (overlap between homeranges of individuals in pair, per pair, colored by subset)
btwnind_vals <- uds_vals %>%
  filter(pair1_ess_before > 2, pair2_ess_before > 2, pair1_ess_after > 2, pair2_ess_after > 2) %>%
  filter(est < 100) %>%
  dplyr::select(-c("low","est","high")) %>%
  mutate(pair = paste(`bear#1`, " ", `bear#2`, "\n", year(time), sep="")) %>%
  pivot_longer(c("uds_before_low", "uds_after_low", "uds_before_mean", "uds_after_mean", "uds_before_high", "uds_after_high"), 
               names_pattern = "^((?:[^_]*_){1}[^_]*)_(.*)$", 
               names_to = c("encounter", "CI")) %>%
  pivot_wider(names_from = CI, values_from = value) %>%
  mutate(encounter=sub("^[^_]*_", "", encounter),
         encounter = factor(encounter, levels = c("before","after"))) %>%
  unite(c("cubs", "season_breeding"), col="order", sep = " ", remove=FALSE) %>%
  mutate(order = fct_recode(order, `Cubs Involved\nLate Fall` = "yes hunting",
                            `Late Fall` = "no hunting",
                            `All Other` = "yes other", `All Other` = "no other", 
                            `All Other` = "yes breeding", `All Other` = "no breeding"),
         order = fct_relevel(order, c("Cubs Involved\nLate Fall", "Late Fall", "All Other"))) %>%
  mutate(lab_color = ifelse(order=="Cubs Involved\nLate Fall", "#ff7678", 
                            ifelse(order=="Late Fall", "#d9d2e9", "#cbfcbd"))) %>%
  arrange(match(order, levels(order))) %>%
  mutate(pair = factor(pair, levels = c(unique(pair))))

ind_plots <- list()
p <- NULL

for (i in seq(1, length(btwnind_vals$pair1), 2)) {
  
  #red = M_M, purple = F_F, green = F_M, yellow = F_M hunting
  order <- btwnind_vals %>%
    slice(i) %>%
    dplyr::select(order)
  
  lab_color <- ifelse(order=="Cubs Involved\nLate Fall", "#ff7678", 
                      ifelse(order=="Late Fall", "#d9d2e9", "#cbfcbd"))
  
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
  scale_fill_manual(values=c("#ff7678", "#d9d2e9", "#cbfcbd"),
                    labels=c("Cubs Involved\nLate Fall", "Late Fall", "All Other")) +
  scale_color_manual(values=c("purple", "orange"),
                     labels = c("Before Encounter", "After Encounter")) +
  guides(fill=guide_legend(nrow = 1, byrow = T, order =2,
                           override.aes=list(size = 15, shape = 15, color = NULL)),
         color=guide_legend(nrow = 1, byrow = T, order=1,
                            override.aes=list(size = 15, shape = 16, fill=c("purple", "orange")))) +
  theme(legend.background=element_blank(), legend.position = "right", legend.justification="left")

#plot and save
grid.draw(grid.arrange(arrangeGrob(grobs = c(ind_plots, 
                                             #list(grid.rect(gp=gpar(col="white"))), #blank to space plot/ legend
                                             list(as_grob(get_legend(legend)))), #convert legend to grob
                                   ncol=6, left=textGrob("Range Distribution Overlap (Before vs. After)", rot=90, 
                                                         gp = gpar(col = "black", fontsize = 20)))))

#################
#Fig. S2 (overlap within individual homerange before/after encounter, per pair, colored by subset and sex)
bears_ind_sex <- bears_meta %>%
  select(c("Capture ID", "Bear #", "Sex")) %>%
  rename("ind_sex"=Sex)

ind_vals <- uds_ind_vals %>%
  filter(ind1_ess_before > 2, ind2_ess_before > 2, ind1_ess_after > 2, ind2_ess_after > 2) %>%
  dplyr::select(-c("low","est","high")) %>%
  mutate(pair = paste(`bear_id1`, " ", `bear_id2`, "\n", year(time), sep="")) %>%
  pivot_longer(c("uds_ind1_low", "uds_ind2_low", "uds_ind1_mean", 
                 "uds_ind2_mean", "uds_ind1_high", "uds_ind2_high"), 
               names_pattern = "^((?:[^_]*_){1}[^_]*)_(.*)$", 
               names_to = c("who", "CI")) %>%
  pivot_wider(names_from = CI, values_from = value) %>%
  mutate(who=sub("^[^_]*_", "", who),
         who = factor(who, levels = c("ind1","ind2"))) %>%
  unite(c("cubs", "season_breeding"), col="order", sep = " ", remove=FALSE) %>%
  mutate(order = fct_recode(order, `Cubs Involved\nLate Fall` = "yes hunting",
                            `Late Fall` = "no hunting",
                            `All Other` = "yes other", `All Other` = "no other", 
                            `All Other` = "yes breeding", `All Other` = "no breeding"),
         order = fct_relevel(order, c("Cubs Involved\nLate Fall", "Late Fall", "All Other"))) %>%
  arrange(match(order, levels(order)))

ind_plots <- list()
p <- NULL
for (i in seq(1, length(ind_vals$time), 2)) {
  
  #red = M_M, purple = F_F, green = F_M, yellow = F_M hunting
  order <- ind_vals %>%
    slice(i) %>%
    dplyr::select(order)
  
  lab_color <- ifelse(order=="Cubs Involved\nLate Fall", "#ff7678", 
                      ifelse(order=="Late Fall", "#d9d2e9", "#cbfcbd"))
  
  
  p <- ind_vals %>%
    slice(i:(i+1)) %>%
    mutate(who = ifelse(#who=="ind1", ind1, ind2),
      who=="ind1", `bear#1`, `bear#2`),
      header = year(time)) %>%
    #left_join(bears_ind_sex, by = c("who"="Capture ID")) %>%
    left_join(bears_ind_sex, by = c("who"="Bear #")) %>%
    #mutate(ind_sex=factor(ind_sex, labels=c("F","M")))
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
  scale_fill_manual(values=c("#ff7678", "#d9d2e9", "#cbfcbd"),
                    labels=c("Cubs Involved\nLate Fall", "Late Fall", "All Other")) +
  scale_color_manual(values=c("red", "blue"),
                     labels = c("Male", "Female"),
                     drop=TRUE) +
  guides(fill=guide_legend(nrow = 1, byrow = T, order=2,
                           override.aes=list(size = 15, shape = 15, color = NULL)),
         color=guide_legend(nrow = 1, byrow = T, order=1,
                            override.aes=list(size = 15, shape = 16, fill=c("red", "blue")))) +
  theme(legend.background=element_blank(), legend.position = "right", legend.justification="left")

#plot/save
grid.draw(grid.arrange(arrangeGrob(grobs = c(ind_plots, 
                                             #list(grid.rect(gp=gpar(col="white"))), #blank to space plot/ legend
                                             list(as_grob(get_legend(legend)))), #convert legend to grob
                                   ncol=6, left=textGrob("Range Distribution Overlap Within Individuals", rot=90, 
                                                         gp = gpar(col = "black", fontsize = 20)))))


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
encounterlocations <- data.frame("bear1"=bearpairs_meta$pair1,
                                 "bear2"=bearpairs_meta$pair2,
                                 "time"=bearpairs_meta$time, "dist"=bearpairs_meta$est,
                                 "season_breeding"=bearpairs_meta$season_breeding,
                                 "sex"=bearpairs_meta$sex,
                                 "cubs"=bearpairs_meta$cubs,
                                 "pair1_ess_before"=bearpairs_meta$pair1_ess_before,
                                 "pair2_ess_after"=bearpairs_meta$pair2_ess_before,
                                 "pair1_ess_after"=bearpairs_meta$pair1_ess_after,
                                 "pair2_ess_before"=bearpairs_meta$pair2_ess_after,
                                 "midpt_lat"=NA, "midpt_long"=NA,
                                 "bear1_lat"=NA, "bear1_long"=NA,
                                 "bear2_lat"=NA, "bear2_long"=NA,
                                 "dist_to_pit"=NA, "nearest_pit"=NA) %>%
  filter(pair1_ess_before > 2, pair2_ess_before > 2, pair1_ess_after > 2, pair2_ess_after > 2) %>%
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
  
  #remove day of encounter from data (figured whole day was ok) and filter to just that encounter year/season
  bear1day <- data.frame(bear1) %>%
    mutate(individual.local.identifier=bearpairs$pair1[[i]]) %>%
    mutate(timestamp = as_date(timestamp)) %>%
    distinct(timestamp, .keep_all=T) %>%
    filter(timestamp != as_date(bearpairs$time[[i]])) %>% #remove day of encounter
    filter(year(timestamp) == year(bearpairs$time[[i]])) %>%  #filter to just that year
    #add in season_breeding and filter to that
    mutate("season_breeding" = ifelse(parse_number(str_replace(str_sub(timestamp,6,10),"-",""))>600 &
                                        parse_number(str_replace(str_sub(timestamp,6,10),"-","")) < 750,
                                      "breeding",
                                      ifelse(parse_number(str_replace(str_sub(timestamp,6,10),"-",""))>900 &
                                               parse_number(str_replace(str_sub(timestamp,6,10),"-","")) < 1130,
                                             "hunting",
                                             "other"))) %>%
    filter(season_breeding==bearpairs$season_breeding[[i]]) %>%
    as.telemetry()
  
  bear2day <- data.frame(bear2) %>%
    mutate(individual.local.identifier=bearpairs$pair2[[i]]) %>%
    mutate(timestamp = as_date(timestamp)) %>%
    distinct(timestamp, .keep_all=T) %>%
    filter(timestamp != as_date(bearpairs$time[[i]])) %>%#remove day of encounter
    filter(year(timestamp) == year(bearpairs$time[[i]])) %>%  #filter to just that year
    #add in season_breeding and filter to that
    mutate("season_breeding" = ifelse(parse_number(str_replace(str_sub(timestamp,6,10),"-",""))>600 &
                                        parse_number(str_replace(str_sub(timestamp,6,10),"-","")) < 750,
                                      "breeding",
                                      ifelse(parse_number(str_replace(str_sub(timestamp,6,10),"-",""))>900 &
                                               parse_number(str_replace(str_sub(timestamp,6,10),"-","")) < 1130,
                                             "hunting",
                                             "other"))) %>%
    filter(season_breeding==bearpairs$season_breeding[[i]]) %>%
    as.telemetry()
  
  allpits_dist1 <- data.frame()
  allpits_dist2 <- data.frame()
  
  for (j in seq_along(carcasspit_locs$pit_name)) {
    
    # print(as.Date(encounterlocations$time[[i]]))
    # print(interval(start=carcasspit_locs$`date active`[[j]],
    #             end=carcasspit_locs$`date inactive`[[j]]))
    
    if(as.Date(bearpairs$time[[i]]) %within% 
       interval(start=carcasspit_locs$`date active`[[j]],
                end=carcasspit_locs$`date inactive`[[j]])) {
      
      #dist calculations turn out to be essentially the same with the distances() telemetry obj,
      #but it's more efficient with geosphere bc there are so many more nonencounter fixes  so just go with that here
      
      dists_1 <- unlist(lapply(seq(nrow(bear1day)), function(x) {
        
        # pt_bear <- as.telemetry(data.frame(location.lat=bear1day$latitude[x],
        #                  location.long=bear1day$longitude[x],
        #                  timestamp=encounterlocations$time[i],
        #                  individual.local.identifier="bear",
        #                  individual.taxon.canonical.name=NA))
        # 
        # carcasspit_tel <- as.telemetry(data.frame(location.lat=carcasspit_locs$lat[j],
        #                  location.long=carcasspit_locs$long[j],
        #                  timestamp=encounterlocations$time[i],
        #                  individual.local.identifier="pit",
        #                  individual.taxon.canonical.name=NA))
        # 
        # projection(pt_bear) = projection(carcasspit_tel)
        
        #return(ctmm::distance(list(pt_bear, carcasspit_tel),method="Euclidean",sqrt=T)$CI[6])
        
        return(geosphere::distHaversine(unlist(bear1day[x, c("longitude", "latitude")]),
                                        unlist(carcasspit_locs[j, c("long", "lat")])))
        
      }))
      
      dists_2 <- unlist(lapply(seq(nrow(bear2day)), function(x) {
        
        # pt_bear <- as.telemetry(data.frame(location.lat=bear2day$latitude[x],
        #                  location.long=bear2day$longitude[x],
        #                  timestamp=encounterlocations$time[i],
        #                  individual.local.identifier="bear",
        #                  individual.taxon.canonical.name=NA))
        # 
        # carcasspit_tel <- as.telemetry(data.frame(location.lat=carcasspit_locs$lat[j],
        #                  location.long=carcasspit_locs$long[j],
        #                  timestamp=encounterlocations$time[i],
        #                  individual.local.identifier="pit",
        #                  individual.taxon.canonical.name=NA))
        # 
        # projection(pt_bear) = projection(carcasspit_tel)
        
        #return(ctmm::distance(list(pt_bear, carcasspit_tel),method="Euclidean",sqrt=T)$CI[6])
        
        return(geosphere::distHaversine(unlist(bear2day[x, c("longitude", "latitude")]),
                                        unlist(carcasspit_locs[j, c("long", "lat")])))
        
      }))
      
    }
    
    else {
      
      #print("nah")
      dists_1 <- rep(NA, length(bear1day$timestamp))
      dists_2 <- rep(NA, length(bear2day$timestamp))
    }
    
    allpits_dist1 <- rbind(allpits_dist1, c(carcasspit_locs$pit_name[[j]], dists_1))
    allpits_dist2 <- rbind(allpits_dist2, c(carcasspit_locs$pit_name[[j]], dists_2))
    
  }
  
  allpits_dist1 <- allpits_dist1 %>% 
    column_to_rownames("X.Olsen.") %>% 
    t() %>% 
    data.frame() %>% 
    `rownames<-`(NULL)
  
  allpits_dist2 <- allpits_dist2 %>% 
    column_to_rownames("X.Olsen.") %>% 
    t() %>% 
    data.frame() %>% 
    `rownames<-`(NULL)
  
  #get min dist to ANY pit at each timestamp
  allpits_mindist1 <- apply(allpits_dist1, 1, min, na.rm=T)
  allpits_mindist2 <- apply(allpits_dist2, 1, min, na.rm=T)
  
  nonencounterlocations <- rbind(nonencounterlocations,
                                 cbind("bear"=bearpairs$pair1[[i]], 
                                       "time"=bearpairs$time[[i]], 
                                       "dist"=bearpairs$est[[i]], 
                                       "season_breeding"=bearpairs$season_breeding[[i]], 
                                       "sex"=bearpairs$sex[[i]], "cubs"=bearpairs_meta$cubs[[i]], 
                                       "dist_from_pit"=allpits_mindist1)) %>%
    rbind(cbind("bear"=bearpairs$pair2[[i]], 
                "time"=bearpairs$time[[i]], 
                "dist"=bearpairs$est[[i]], 
                "season_breeding"=bearpairs$season_breeding[[i]], 
                "sex"=bearpairs$sex[[i]], "cubs"=bearpairs$cubs[[i]],
                "dist_from_pit"=allpits_mindist2))
  
}

#Make figures!

#pts removed for being >20000m away from any carcass pit
#pts removed for being >20000m awat
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

n_removed_cubs_fall <- (nonencounterlocations %>%
                          filter(cubs=="yes", season_breeding=="hunting") %>%
                          filter(dist_from_pit > 50000) %>%
                          count() %>%
                          summarize(n=n/nrow(filter(nonencounterlocations, cubs=="yes", season_breeding=="hunting"))))$n

#Figure 3. Distance from carcass pits based on sex and cubs
nonencounterlocations %>%
  filter(dist_from_pit < 50000) %>%
  unite(c("cubs", "season_breeding"), col="order", sep = " ", remove=FALSE) %>%
  mutate(order = fct_recode(order, `Cubs Involved\nLate Fall` = "yes hunting",
                            `Late Fall` = "no hunting",
                            `All` = "yes other", `All` = "no other", 
                            `All` = "yes breeding", `All` = "no breeding"),
         order = fct_relevel(order, c("All", "Late Fall", "Cubs Involved\nLate Fall"))) %>%
  #scale density to make it visible with counts: https://stackoverflow.com/questions/72071765/how-to-add-a-second-variable-to-histogram-ggplot-and-plot-on-top-current-histogr
  ggplot(mapping=aes(x=dist_from_pit, y=after_stat(density)*40000)) +
  geom_density(aes(color=order), show.legend=F) +
  stat_density(aes(color=order), geom="line") +
  geom_histogram(encounterlocations %>%
                   filter(dist_to_pit < 50000) %>%
                   unite(c("cubs", "season_breeding"), col="order", sep = " ", remove=FALSE) %>%
                   mutate(order = fct_recode(order, `Cubs Involved\nLate Fall` = "yes hunting",
                                             `Late Fall` = "no hunting",
                                             `All` = "yes other", `All` = "no other", 
                                             `All` = "yes breeding", `All` = "no breeding"),
                          order = fct_relevel(order, c("All", "Late Fall", 
                                                       "Cubs Involved\nLate Fall"))),
                 mapping=aes(x=dist_to_pit, y=after_stat(count), fill=order),
                 alpha=0.6, size=0.1) +
  scale_color_manual(values=c("lightgreen", "#C688E3", "#ff7678")) +
  scale_fill_manual(values=c("lightgreen", "#C688E3", "#ff7678")) +
  facet_wrap(~order, scales="free_y") +
  #labeller = labeller(order = c(hunting="Late Fall FM"))) +
  labs(x="Distance from nearest carcass pit (m)", y="Count of encounter fixes", 
       fill="Encounter fixes", color="Non-encounter fixes") +
  scale_y_continuous(sec.axis = sec_axis(~.x/40000, name = "Density of non-encounter fixes")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90),
        legend.position="bottom",
        legend.direction = "horizontal", 
        legend.box = "horizontal") + 
  guides(fill = guide_legend(order=1, title.position="top",
                             override.aes = list(alpha=1)),
         color = guide_legend(order=2, title.position="top",
                              override.aes = list(alpha=1, fill=NA, lwd=1)))



