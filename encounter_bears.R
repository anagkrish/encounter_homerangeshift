#This file contains the code to fit home range models and calculate overlap for all 29 bear encounters 
#and also to run the meta population level analysis showing changes in overlap following encounters.

#lodd libraries
library(tidyverse)
library(ctmm)
library(lubridate)
library(doParallel)
library(parallel)
library(pals)
library(ggpattern)

#function to calculate overlap between individuals before/after encounter
get_overlap = function(data, name1, name2, date_encounter, timestart, timestop) {
  
  names <- c()
  for (i in data) {
    names <- c(names,i[1]@info$identity)
  } #get list of individual names to subset data with 
  
  id1 <- dat[[which(names==name1)]]
  id1 <- id1[id1$timestamp > timestart,]
  id1 <- id1[id1$timestamp < timestop,] #subset individual 1 between specified dates
  
  id2 <- dat[[which(names==name2)]]
  id2 <- id2[id2$timestamp > timestart,]
  id2 <- id2[id2$timestamp < timestop,] #subset individual 2 between specified dates
  
  #split data at time and get before/after
  before <- list(id1[id1$timestamp < date_encounter,],
                 id2[id2$timestamp < date_encounter,]) 
  
  after <- list(id1[id1$timestamp > date_encounter,],
                id2[id2$timestamp > date_encounter,])
  
  #correct all projections
  ctmm:::projection(before[1]) <- median(before)  
  ctmm:::projection(before[2]) <- median(before)  
  ctmm:::projection(after[1]) <- median(before)  
  ctmm:::projection(after[2]) <- median(before)
  
  #fit home range models for before and after encounter for each individuals 
  akde_all <- function(d){
    guess <- map(d, ctmm.guess, interactive = FALSE)
    fit <- map2(d, guess, ctmm.select, cores = 2, trace = 2)
    return(fit)
  }
  
  seq_tele <- list(before = before, after = after)
  
  #fit ctmm models to each home range
  model_list <- foreach(i = seq_along(seq_tele),
                        .packages = c("ctmm","lubridate","tidyverse","purrr","stringr")) %dopar% {
                          akde_all(seq_tele[[i]])}
  
  before_model_fit <- model_list[[1]]
  after_model_fit <- model_list[[2]]
  
  #calculate UDS overlap of home ranges before and after encounter
  UDS_before <- akde(seq_tele[[1]], before_model_fit)
  UDS_after <- akde(seq_tele[[2]], after_model_fit)
  
  #calculate overlap of home ranges before and after encounter
  before_overlap <- overlap(UDS_before)
  after_overlap <- overlap(UDS_after)
  
  return(list(before_overlap, after_overlap))
  
}

#-----------------------------------------------------------------------------------------
#part 1: calculate overlaps for all 29 encounter pairs
#info on all 29 pairs, including time of encounter and time range, included in file bearpairs_final.csv

#initialize folder to store overlap objects in
folder <- "/your/folder/name/here"

#load file with selected bear pairs (n=29)
bearpairs <- read_csv("bearpairs_final.csv")

#load bears data
#data should be in format of a large list where each element is a telemetry object for each individual 
load("/your/grizzly/data/here")

#run overlap function on all individuals per pair 
foreach (i = seq_along(bearpairs$time),
         .packages = c("tidyverse","ctmm","lubridate","foreach","doSNOW","parallel"),
         .combine = c) %dopar% {
           
           print(paste(bearpairs$pair1[i], bearpairs$pair2[i],
                       ymd_hms(bearpairs$trackstart[i]), ymd_hms(bearpairs$trackstop[i]),
                       ymd_hms(bearpairs$time[i])))
           
           overlap <- get_overlap(bears, bearpairs$pair1[i], bearpairs$pair2[i], ymd_hms(bearpairs$time[i]),
                                  ymd_hms(bearpairs$trackstart[i]), ymd_hms(bearpairs$trackstop[i]))
           #creates a list with two items where the first item is overlap of homeranges before the encounter
           #and the second item is overlap of homeranges after the encounter
           
           save(overlap, file=(paste(folder, "/", bearpairs$pair1[i], bearpairs$pair2[i], 
                                     bearpairs$time[i], "overlaps.rda"))) #save to specified folder
           
           }

#-----------------------------------------------------------------------------------------
#part 2: load all overlaps from folder and get UDS overlaps

before <- list()
after <- list() #initialize empty lists

uds_vals <- data.frame(pair=NA, 
                       uds_before_low=NA, uds_before_mean=NA, uds_before_high=NA,
                       uds_after_low=NA, uds_after_mean=NA, uds_after_high=NA)

for (i in seq_along(list.files(folder))) {

  #load objects and individually append before/after overlap distributions from each overlap object
  load(paste(folder, list.files(folder)[[i]], sep ="/"))
  before <- c(before, list(overlap[[1]])) 
  after <- c(after, list(overlap[[2]]))
  
  uds_vals <- add_row(uds_vals, 
                      pair = str_replace(list.files(folder)[[i]][1], "2.*",""), 
                      uds_before_low=overlap[[1]][["CI"]][[2]], 
                      uds_before_mean=overlap[[1]][["CI"]][[6]], 
                      uds_before_high=overlap[[1]][["CI"]][[10]],
                      uds_after_low=overlap[[2]][["CI"]][[2]], 
                      uds_after_mean=overlap[[2]][["CI"]][[6]], 
                      uds_after_high=overlap[[2]][["CI"]][[10]])
  
}

#rearrange/neaten UDS values data frame
uds_vals <- uds_vals %>%
  drop_na() %>%
  mutate(pair = gsub('([[:upper:]])', ' \\1', pair)) %>%
  separate(pair, into=c("drop","pair1", "pair2"), sep = " ") %>%
  mutate(pair2 = replace(pair2, pair2 == "E", "EVGF84")) %>%
  dplyr::select(-c("drop"))

#use meta function to compare Bhattacharryya distances of before/after overlap of all grizzlies in population
meta <- meta(list(before = before, after = after),
             level=0.95, col=c("red", "purple"))

meta

#-----------------------------------------------------------------------------------------
#part 3: code to generate Fig. 2

#load bhattacharyya distances for various subsets:
#All encounters,encounters between same sex individuals, encounters between diff sex individuals,
#encounters between diff sex individuals during late fall
bhattacharyya <- read_csv("bhattacharyya.csv")

bhattacharyya %>%
  mutate(split = factor(split, levels = c("before", "after")), #make variable as factor so bars are in the correct order
         subset = factor(subset, levels = c("all", "same", "diff", "diffhunting"))) %>%
  ggplot(mapping = aes(x=subset, y=est, group=split, pattern=split)) +
  geom_col_pattern(position=position_dodge(0.85), width=0.85, pattern_color = "#7A7A7A",
                   fill = "white", color = "black", pattern_spacing = 0.015,
                   pattern_frequency = 5, pattern_angle = 45) +
  geom_errorbar(aes(ymin=low, ymax=high), width = .25, position = position_dodge(0.85)) +
  scale_pattern_manual(values=c('none', 'stripe'), 
                       labels=c("Before Encounter","After Encounter")) +
  labs(x=NULL, y="Bhattacharyya Distance", pattern=NULL) +
  scale_x_discrete(labels=c('ALL', 'SAME SEX', 'DIFF SEX', 'DIFF SEX \nLATE FALL'),
                   expand=c(-1,0)) +
  scale_y_continuous(expand = c(0, 0.001), limits=c(0,5.5)) + #removes space between bars and axis line
  geom_text(aes(label=c("n=29",NA,"n=13",NA,"n=16",NA,"n=8",NA)), 
            y = rep(c(4), times = 8)) +
  theme_bw() + #get rid of background
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth=0.5),
        text=element_text(size=20), legend.position =c(.15,.90))


