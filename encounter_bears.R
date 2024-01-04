#This file contains the code to fit home range models and calculate overlap for all bear encounters 
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
#returns a list of overlap before and overlap after
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
  
  #bootstrap implememted for low ESS (2 < ess < 5) encounter subsets
  bootstrap <- function(dat, mod) {
    
    #bootstrap when ess < 5 and > 2
    purrr::map2(dat, mod, function(dat, mod) {
      if(summary(mod)$DOF[["area"]]<5&summary(mod)$DOF[["area"]]>2) {
        return(ctmm.boot(dat, mod, error=0.1, clamp=0.001, iterate=T, trace=2))
      }
      else {
        return(mod)
      }
    })
  }
  
  #run bootstrap in parallel
  model_boot <- foreach(i = seq_along(model_list),
                        .packages = c("ctmm","lubridate","tidyverse","purrr","stringr",
                                      "snow","doSNOW","foreach","parallel","doParallel")) %dopar% {
                                        bootstrap(seq_tele[[i]], model_list[[i]])}
  
  before_model_boot <- model_boot[[1]]
  after_model_boot <- model_boot[[2]]
  
  #calculate UDS overlap of home ranges before and after encounter
  UDS_before <- akde(seq_tele[[1]], before_model_fit)
  UDS_after <- akde(seq_tele[[2]], after_model_fit)
  
  #bootstrapped fits UDS
  UDS_before_boot <- akde(seq_tele[[1]],before_model_boot)
  UDS_after_boot <- akde(seq_tele[[2]],after_model_boot)
  
  #calculate overlap of home ranges before and after encounter
  before_overlap <- overlap(UDS_before)
  after_overlap <- overlap(UDS_after)
  
  #bootstrapped overlap UDS
  boot_before_overlap <- overlap(UDS_before_boot)
  boot_after_overlap <- overlap(UDS_after_boot)
  
  #return all info
  return(list(before, after, before_model_fit, after_model_fit, before_model_boot, after_model_boot,
              UDS_before, UDS_after, UDS_before_boot, UDS_after_boot, before_overlap, after_overlap,
              boot_before_overlap, boot_after_overlap))
  
}

#function to calculate overlap within individual home range before/after encounter 
#returns a list of overlap for individual 1 and overlap for individual 2
get_overlap_individual = function(dat, name1, name2, date_encounter, timestart, timestop) {
  
  names <- c()
  for (i in dat) {
    names <- c(names,i[1]@info$identity)
  }
  
  id1 <- dat[[which(names==name1)]]
  id1 <- id1[id1$timestamp > timestart,]
  id1 <- id1[id1$timestamp < timestop,] #subset individual 1 between specified dates
  
  id2 <- dat[[which(names==name2)]]
  id2 <- id2[id2$timestamp > timestart,]
  id2 <- id2[id2$timestamp < timestop,] #subset individual 2 between specified dates
  
  ind1 <- list(id1[id1$timestamp < date_encounter,],
               id1[id1$timestamp > date_encounter,])
  
  ind2 <- list(id2[id2$timestamp < date_encounter,],
               id2[id2$timestamp > date_encounter,])
  
  #standardize projections
  ctmm:::projection(ind1[1]) <- median(ind1)   #correct the projection --Qianru
  ctmm:::projection(ind1[2]) <- median(ind1)   #correct the projection --Qianru
  ctmm:::projection(ind2[1]) <- median(ind1)   #correct the projection --Qianru
  ctmm:::projection(ind2[2]) <- median(ind1)   #correct the projection --Qianru
  
  akde_all <- function(d){
    guess <- map(d, ctmm.guess, interactive = FALSE)
    fit <- map2(d, guess, ctmm.select, cores = 4, trace = 2) #use all cores
    #akde_final <- map2(data, fit, akde)
    return(fit)
    #, akde = akde_final))
  }
  
  seq_tele <- list(ind1 = ind1, ind2 = ind2)
  model_list <- foreach(i = seq_along(seq_tele), 
                        .packages = c("ctmm","lubridate","tidyverse","purrr","stringr")) %dopar% {
                          akde_all(seq_tele[[i]])}
  
  ind1_fit <- model_list[[1]]
  ind2_fit <- model_list[[2]]
  
  #bootstrap function for encounter subsets with low ESS
  bootstrap <- function(dat, mod) {
    
    #bootstrap when ess < 5 and > 2
    purrr::map2(dat, mod, function(dat, mod) {
      
      if(summary(mod)$DOF[["area"]]<5&summary(mod)$DOF[["area"]]>2) {
        return(ctmm.boot(dat, mod, error=0.1, clamp=0.001, iterate=T, trace=2))
      }
      else {
        return(mod)
      }
    })
    
  }
  
  model_boot <- foreach(i = seq_along(model_list),
                        .packages = c("ctmm","lubridate","tidyverse","purrr","stringr",
                                      "snow","doSNOW","foreach","parallel","doParallel")) %dopar% {
                                        bootstrap(seq_tele[[i]], model_list[[i]])}
  
  
  #calculate UDS overlap of home ranges before and after encounter
  UDS_ind1 <- akde(seq_tele[[1]], ind1_fit)
  UDS_ind2 <- akde(seq_tele[[2]], ind2_fit)
  
  UDS_boot_ind1 <- akde(seq_tele[[1]], ind1_boot)
  UDS_boot_ind2 <- akde(seq_tele[[2]], ind2_boot)
  
  #calculate overlap of home ranges before and after encounter
  ind1_overlap <- overlap(UDS_ind1)
  ind2_overlap <- overlap(UDS_ind2)
  
  ind1_boot_overlap <- overlap(UDS_boot_ind1)
  ind2_boot_overlap <- overlap(UDS_boot_ind2)
  
  return(list(ind1, ind2, ind1_fit, ind2_fit, ind1_boot, ind2_boot,
                UDS_ind1, UDS_ind2, UDS_boot_ind1, UDS_boot_ind2, ind1_overlap, ind2_overlap,
                ind1_boot_overlap, ind2_boot_overlap))
  
}

#-----------------------------------------------------------------------------------------
#part 1: get overlaps from original movement data

#load bear movement data
bears <- "/movebank/file/here"

#filter for individuals with overlapping time
done <- list()
pair1 <- c()
pair2 <- c()
overlap <- c() 

for (i in bears) {
  done <- c(done,i@info$identity) #keep track of individuals you've already checked
  
  for (j in bears) {
    
    if (j@info$identity != i@info$identity && j@info$identity %notin% done) {
        
      #check for overlap between time periods
      if(int_overlaps(as.interval(i$timestamp[1],tail(i$timestamp,1)), 
                      as.interval(j$timestamp[1],tail(j$timestamp,1))) == TRUE) {
        
        pair1 <- c(pair1, i@info$identity)
        pair2 <- c(pair2, j@info$identity)
        
        #calculate overlap time in weeks
        overlap <- c(overlap, as.duration(intersect(as.interval(i$timestamp[1],
                                                                tail(i$timestamp,1)), 
                                                    as.interval(j$timestamp[1],
                                                                tail(j$timestamp,1))))/604800)
        
      }
    }
  }
}

candidates <- data.frame(pair1,pair2,overlap)

#check distance between all candidates

threshold = 100 #encounter threshold in meters

dists <- data.frame() #initialize data frame

for (i in seq_along(candidates$pair1)) { 
  
  dat <- bears[c(which(names==candidates$pair1[i]),
                 which(names==candidates$pair2[i]))] 
  
  #fit ctmm models to data
  fits <- list()
  for(i in seq_along(dat)) {
    guess <- ctmm.guess(dat[[i]],interactive=FALSE)
    fits[[i]] <- ctmm.select(dat[[i]],guess, trace=2,cores=2)
  }
  
  dat[[2]]@info$projection <- dat[[1]]@info$projection #correct projection
  
  #calculates point distances with uncertainty between two movement tracks
  diffs <- ctmm::distances(dat, fits) %>%
    filter(est<=threshold) %>%
    mutate(pair1=dat[[1]]@info$identity, pair2=dat[[2]]@info$identity) %>%
    dplyr::select(c("pair1","pair2","timestamp","low","est","high"))
  
  #contains ALL distances for each pair that are estimated below threshold; 
  #may be in same time period across multiple years depending on if there were multiple encounters/individual 
  dists <- rbind(dists, diffs)
  
}

#summarize individual pairs (by year for pairs with multiple encounters) 
#and lowest estimated encounter distance for each pair

summary <- dists %>%
  group_by(pair1, pair2, year(timestamp)) %>%
  summarize(est=min(est))

#-----------------------------------------------------------------------------------------
#part 2: calculate overlaps across selected individuals

#initialize folder to store overlap objects in
folder <- "/your/folder/name/here"

#load file with candidate pairs
#for threshold=100, n=32
bearpairs <- read_csv("bearpairsto500m.csv")

#load bears data
#data should be in format of a large list where each element is a telemetry object for each individual 
load("/your/grizzly/data/here")

#run overlap function on all individuals per pair 
#written to run in parallel (but can be modified)
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
#part 3: load all overlaps from folder and get UDS overlaps

before <- list()
after <- list() #initialize empty lists

#(optional) can extract and store UDS values as well
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
  mutate(pair = gsub('([[:upper:]])', ' \\1', pair)) %>% #space between words (ideally each name but formatting is a bit iffy, hence the other steps)
  separate(pair, into=c("pair1", "pair2"), sep = "^\\s*\\S+\\K\\s+") %>% #split only at first instance of space
  rowwise() %>%
  mutate(pair1 = gsub(" ", "", pair1),
         pair2 = gsub(" ", "", pair2)) %>%
  ungroup() %>%
  mutate(pair1 = ifelse(pair1=="E", "EVGF84", pair1),
         pair2 = ifelse(pair2=="VGF84Stu", "Stu", pair2)) #correct one instancce where split messes up


#use meta function to compare Bhattacharryya distances of before/after overlap of all grizzlies in population
meta <- meta(list(before = before, after = after),
             level=0.95, col=c("red", "purple"))

meta

#calculate relative percent overlap values (translate BD into physical overlap change)
#value represents what percentage the after overlap is of the the before overlap

#get percent reduction in overlap
BM <- meta(before,level=NA)[1,2:3] # mean and variance
AM <- meta(after,level=NA)[1,2:3] # mean and variance

# difference point estimate
D <- AM[1] - BM[1]
# difference variance
VAR <- AM[2] + BM[2]

# distance difference CI (can be negative)
CI <- ctmm:::norm.ci(D,VAR=VAR,level=0.90)
# overlap ratio CI
CI <- rev(exp(-CI))
names(CI) <- ctmm:::NAMES.CI

# On average, AFTER overlap is only
sigfig(100*CI)

#-----------------------------------------------------------------------------------------
#part 4: code to generate Fig. 2

#load relative overlap for various subsets of encounters at 100m:
#all individuals, same sex individuals, diff sex individuals, individuals during late fall, individuals during spring, diff sex individuals during late fall
rel_overlap <- read_csv("bd_reloverlap_500m.csv") %>%
  dplyr::select(subset, threshold, rel_overlap_low, rel_overlap_est, rel_overlap_high)

rel_overlap %>%
  filter(threshold==100,
         subset %in% c("all", "late fall", "diff sex late fall", "cubs involved late fall")) %>%
  mutate(subset=factor(subset, levels = c("all", "late fall", 
                                          "diff sex late fall",
                                          "cubs involved late fall"))) %>%
  drop_na(subset) %>%
  ggplot(mapping=aes(x=subset, y=rel_overlap_est)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin=rel_overlap_low, ymax=rel_overlap_high), width=0.2) +
  geom_hline(mapping=aes(yintercept=1), linetype="dashed") +
  scale_x_discrete(breaks = c("all", "late fall",  "diff sex late fall", "cubs involved late fall"),
                   labels=c('ALL', "LATE FALL",
                            "DIFF SEX \nLATE FALL", "CUBS INVOLVED \nLATE FALL"),
                   expand=c(0.1,0.1)) +
  labs(x=NULL, y = "Change in Relative Overlap Between Individuals \n(After / Before Encounter)") +
  scale_y_continuous(breaks=c(0.25, 0.50, 0.75, 1.00, 1.25), limits=c(0.2,1.4)) +
  geom_text(aes(label=c("n=32\n", "n=9\n**", "n=12\n",  "n=19\n"), 
                #all, cubs involved late fall, diff sex late fall, late fall
                y = rep(c(1.32), times = 4)), size=6) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth=0.5),
        text=element_text(size=20), legend.direction="vertical")

#-----------------------------------------------------------------------------------------
#part 5: code to generate Fig 3
#showing difference in relative overlap after/before for encounters at different thresholds, subsetted by all vs different sex

all <- rel_overlap %>%
  filter(subset %in% c("all", "late fall")) %>%
  mutate(subset=factor(subset, levels = c("all", "late fall"))) %>%
  ggplot(mapping=aes(x=threshold, y=rel_overlap_est, group=subset, shape=subset)) +
  geom_point(size = 3, position=position_dodge(width=50)) +
  geom_errorbar(aes(ymin=rel_overlap_low, ymax=rel_overlap_high), position=position_dodge(width=50)) +
  geom_hline(mapping=aes(yintercept=1), linetype="dashed") +
  scale_shape_manual(values = c(1, 4),
                     labels=c("All Encounters", "Late Fall \nEncounters")) +
  labs(x=NULL, y = NULL, shape=NULL) +
  scale_x_continuous(breaks=c(50, 100, 200, 300, 400, 500),
                     labels=c("50", "100", "200", "300", "400", "500"),
                     limits=c(0,550)) +
  scale_y_continuous(breaks=c(0.25, 0.50, 0.75, 1.00, 1.25), limits=c(0.2,1.3)) +
  geom_text(aes(label=c("","","","","","", #all 50 100 200 300 400 500 #pval
                        "*","*","*","","",""), #all hunting 50 100 200 300 400 500 #pval
                y = rep(c(1.2), times = 12)), size=8, position=position_dodge(width=50)) +
  geom_text(aes(label= c("(n)", "", "", "", "", "", "", "", "", "", "", ""), #all hunting 50 100 200 300 400 500 pval
                y = rep(c(0.25), times = 12)), size=5,
            nudge_x= -39) +
  geom_text(aes(label= c("20", "32", "51", "66", "74", "77", #all 50 100 200 300 400 500 pval
                         "11", "19", "29", "33", "40", "41"), #all hunting 50 100 200 300 400 500 pval
                y = rep(c(0.25), times = 12)), size=5,
            position=position_dodge(width=50)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth=0.5),
        text=element_text(size=20), legend.direction="vertical")

diff <- rel_overlap %>%
  filter(subset %in% c("diff sex late fall", "cubs involved late fall")) %>%
  filter(!c(threshold == 500 & subset == "diff sex late fall")) %>%
  mutate(subset=factor(subset, levels = c("diff sex late fall", "cubs involved late fall"))) %>%
  ggplot(mapping=aes(x=threshold, y=rel_overlap_est, group=subset, shape=subset)) +
  geom_point(size = 3, position=position_dodge(width=50)) +
  geom_errorbar(aes(ymin=rel_overlap_low, ymax=rel_overlap_high), position=position_dodge(width=50)) +
  geom_hline(mapping=aes(yintercept=1), linetype="dashed") +
  scale_shape_manual(values = c(19, 8),
                     labels=c("Different Sex \nLate Fall Encounters", "Cubs Involved \nLate Fall Encounters")) +
  labs(x="Encounter Threshold (m)", y = NULL, shape=NULL) +
  scale_x_continuous(breaks=c(50, 100, 200, 300, 400, 500),
                     labels=c("50", "100", "200", "300", "400", "500"),
                     limits=c(0, 550)) +
  scale_y_continuous(breaks=c(0.25, 0.50, 0.75, 1.00, 1.25, 1.5, 1.75), limits=c(0.2,1.77)) +
  geom_text(aes(label= c("**", "**", "", "", "", "", #cubs hunting50 100 200 300 400 500  
                         "", "", "**", "", ""), #diff hunting 50 100 200 300 400 500
                y = rep(c(1.51), times = 11)), size=8, position=position_dodge(width=50)) +
  geom_text(aes(label= c("(n)", "", "", "", "", "", "", "", "", "", ""), #all hunting 50 100 200 300 400 500 pval
                y = rep(c(0.25), times = 11)), size=5,
            nudge_x= -39) +
  geom_text(aes(label= c("8", "9", "15", "18", "21", "22",#cubs hunting50 100 200 300 400 500 pval
                         "7", "12 ", "19", "21", "22"), #diff hunting 50 100 200 300 400 500 pval
                y = rep(c(0.25), times = 11)), size=5,
            position=position_dodge(width=50)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth=0.5),
        text=element_text(size=20), legend.direction="vertical")

yaxis <- cowplot::ggdraw() + 
  cowplot::draw_label("Change in Relative Overlap Between Individuals \n(After / Before Encounter)",
                      size=25,angle=90, x=0.5)

cowplot::plot_grid(yaxis, 
                   cowplot::plot_grid(all, diff, ncol=1,
                                      align = "v",
                                      labels = c("a)","b)")),
                   ncol=2, nrow=1, rel_widths=c(0.07,1))

