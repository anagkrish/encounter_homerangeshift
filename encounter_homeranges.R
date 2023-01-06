#This file contains the code to calculate the home range shifts before and after the encounter
#It can also be used to calculate pseudoencounters by changing the date of the encounter

#Example of implementing functions is at the bottom of the file
#Can be modified to run multiple individuals at a time in parallel (code not shown because it's super customizable)

#load libraries
library(tidyverse)
library(lubridate)
library(ctmm)
library(doSNOW)
library(foreach)
library(gridExtra)
library(purrr)

#-----------------------------------------------------------------------------------------
# functions to fit models and obtain home ranges before and after encounter

#first version by Maria Nikolaitchik 05/2022
#second version by Qianru Liao --08/2022
#third version by me 12/22

#-- function: compare overlap
#Takes in two lists. 
#Before_model takes a list containing two fitted models of the individuals before encounter.
#After_model is the list of fitted models after encounter
#Returns a data.frame that lets you compare values easier

compare_overlap <- function(before_model, after_model){
  
  before_overlap <- overlap(before_model)
  after_overlap <- overlap(after_model)
  low_comp <- c(before_overlap$CI[2], after_overlap$CI[2]) 
  est_comp <- c(before_overlap$CI[6], after_overlap$CI[6]) 
  high_comp <- c(before_overlap$CI[10], after_overlap$CI[10])
  return(data.frame( rows = c("before","after"), low = low_comp, est = est_comp, high = high_comp, 
                     row.names = "rows"))
}

#-- function: run_akde_analysis
#The following is a script that runs all the AKDE analysis and saves all the data into a folder it creates

#This splits the data EXACTLY before and after the encounter time, it does not factor in "transition times"

#You need to input the following:
#dat: movement data from movebank
#name1: the name of individual one (make sure it's spelled exactly)
#name2: the name of individual two (make sure it's spelled exactly)
#study: the study the pair is from
#date_encounter: the date of encounter as a lubridate datetime (use ymd_hms("2016-04-23 14:30:00") to create them)

###NOTE: this code works if dat is imported directly from movebank as a dataframe of multiple individuals. 
#To run instead on a list of telemetry objects, make the following modifications:
  #uncomment lines 67 to 70
  #replace lines 79-83 with the following: 
    #id1 <- dat[[which(names==name1)]]
    #id2 <- dat[[which(names==name2)]]
  #replace lines 85-89 with the following:
    # before <- list(id1[id1$timestamp < date_encounter,],
    #                id2[id2$timestamp < date_encounter,]) 
    # after <- list(id1[id1$timestamp > date_encounter,],
    #               id2[id2$timestamp > date_encounter,]) ######

run_akde_analysis = function(dat, name1, name2, study, date_encounter) {
  #folder_name <- str_c(name1, "_", name2, "_Homeranges",date_encounter, sep = "")
  folder_name <- str_c(gsub("/","",name1), "_", name2, "encounters", 
                       as_date(date_encounter), sep = "")
  
  # names <- c()
  # for (i in dat) {
  #   names <- c(names,i[1]@info$identity)
  # }
  
  #Change the folder_dir string to the folder where you want all your results to be stored
  #Do not add any slashes to the end of the string
  #You ONLY need to change this variable
  folder_dir <- "/Users/anankekrishnan/Documents/faganlab2022"
  folder_path <- str_c(folder_dir, folder_name, sep = "/")
  dir.create(folder_path)
  
  id1 <- dat %>% filter(study.id == study,
                        individual.local.identifier == name1)
  
  id2 <- dat %>% filter(study.id == study,
                        individual.local.identifier == name2)
  
  before <- list(filter(id1, timestamp < date_encounter),
                 filter(id2, timestamp < date_encounter)) %>% purrr::map(., as.telemetry)

  after <- list(filter(id1, timestamp > date_encounter),
                filter(id2, timestamp > date_encounter)) %>% map(., as.telemetry)
  
  save(before, file = str_c(folder_path, "/", folder_name, "_before_data.rda"))
  save(after, file = str_c(folder_path, "/", folder_name, "_after_data.rda"))
  
  #Sets the projections to be the same as before. Change both to after is you want them standardized to after
  ctmm:::projection(before[1]) <- median(before)   #correct the projection --Qianru
  ctmm:::projection(before[2]) <- median(before)   #correct the projection --Qianru
  ctmm:::projection(after[1]) <- median(before)   #correct the projection --Qianru
  ctmm:::projection(after[2]) <- median(before)   #correct the projection --Qianru
  
  #Print the number of entries in each split dataset, gives you a sense of scale
  print(map_dbl(before, nrow))
  print(map_dbl(after, nrow))
  
  akde_all <- function(d){
    guess <- map(d, ctmm.guess, interactive = FALSE)
    fit <- map2(d, guess, ctmm.select, cores = 4, trace = 2) #use all cores
    return(fit)
  }
  
  seq_tele <- list(before = before, after = after)
  model_list <- foreach(i = seq_along(seq_tele), 
                        .packages = c("ctmm","lubridate","tidyverse","purrr","stringr")) %dopar% {
                          akde_all(seq_tele[[i]])}
  
  before_model_fit <- model_list[[1]]
  after_model_fit <- model_list[[2]]
  
  save(before_model_fit, file = str_c(folder_path, "/", folder_name, "_before_model.rda", sep = ""))
  save(after_model_fit, file = str_c(folder_path, "/", folder_name, "_after_model.rda", sep = ""))
  
  #fit overlap analysis table
  overlap_analysis_fit <- compare_overlap(before_model_fit, after_model_fit)
  png(filename = str_c(folder_path, "/", folder_name,"_overlap_table_fit.png", sep = ""),
      height = 50*nrow(overlap_analysis_fit), width = 200*ncol(overlap_analysis_fit))
  grid.table(overlap_analysis_fit)
  dev.off()
  
  #Save Overlap Table as dataframe
  save(overlap_analysis_fit, file = str_c(folder_path, "/", folder_name,"_overlap_table_fit.rda", sep = ""))
  
  #UDS overlap analysis table ---------------added by Qianru
  UDS_before <- akde(seq_tele[[1]],before_model_fit)
  UDS_after <- akde(seq_tele[[2]],after_model_fit)
  overlap_analysis_akde <- compare_overlap(UDS_before, UDS_after)
  png(filename = str_c(folder_path, "/", folder_name,"_overlap_table_UDS.png", sep = ""), 
      height = 50*nrow(overlap_analysis_akde), width = 200*ncol(overlap_analysis_akde))
  grid.table(overlap_analysis_akde)
  dev.off()
  
  #Save UDS Overlap Table as dataframe ---------------added by Qianru
  save(overlap_analysis_akde, file = str_c(folder_path, "/", folder_name,"_overlap_table_UDS.rda", sep = ""))
  
  pair_names <- c(name1, name2)
  
  #Make the effective home range table for reporting  ---------------edited by Qianru
  eff_home_range <- data.frame(x = c(UDS_before[[1]]$DOF.area[1], 
                                     UDS_after[[1]]$DOF.area[1]),
                               y = c(UDS_before[[2]]$DOF.area[1], 
                                     UDS_after[[2]]$DOF.area[1]))
  colnames(eff_home_range) <- pair_names
  rownames(eff_home_range) <- c("Before","After")
  
  png(filename = str_c(folder_path, "/", folder_name,"_effective_homerange_table.png", sep = ""), 
      height = 50*nrow(eff_home_range), width = 200*ncol(eff_home_range))
  grid.table(eff_home_range)
  dev.off()
  
  save(eff_home_range, file = str_c(folder_path, "/", folder_name,"_effective_homerange_table.rda", sep = ""))
  
  #Prepare to plot home ranges before and after
  extents <- rbind(ctmm::extent(before_model_fit[[1]]), ctmm::extent(before_model_fit[[2]]),
                   ctmm::extent(UDS_before[[1]]), ctmm::extent(UDS_before[[2]]),
                   ctmm::extent(after_model_fit[[1]]), ctmm::extent(after_model_fit[[2]]),
                   ctmm::extent(UDS_after[[1]]), ctmm::extent(UDS_after[[2]]))
  
  
  x_range <- c((min(extents$x) - 500 %#% "meters"),(max(extents$x) + 500 %#% "meters"))
  y_range <- c((min(extents$y) - 500 %#% "meters"),(max(extents$y) + 500 %#% "meters"))
  
  sorted_names <- sort(c(name1, name2), decreasing = F)
  
  #Before encounter home range estimates (Fig. 1B)
  png(filename = str_c(folder_path, "/", folder_name,"_before.png", sep = ""),
      width = 13, height = 6, units = "in", res = 75)
  plot(before, col = color(before, by = "individual"), UD = UDS_before,
       ylim = y_range, xlim = x_range, level.UD = 0.95,
       main = str_c(name1, " and ", name2, " before encounter: ",date_encounter))
  legend("bottomleft", legend = sorted_names, col = c("red","blue"), pch = 1, cex =2.5) 
  
  dev.off()
  
  #After encounter home range estimates (Fig. 1C)
  png(filename = str_c(folder_path, "/", folder_name,"_after.png", sep = ""), 
      width = 13, height = 6, units = "in", res = 75)
  plot(after, col = color(after, by = "individual"), UD = UDS_after, 
       ylim = y_range, xlim = x_range, 
       main = str_c(name1, " and ", name2, " after encounter: ",date_encounter))
  legend("bottomleft", legend = sorted_names, col = c("red","blue"), pch = 1, cex = 2.5) 
  dev.off()
  
  #return output to calculate AIC values for home range shift with meta function
  UDS_before <- akde(seq_tele[[1]],before_model_fit)
  UDS_after <- akde(seq_tele[[2]],after_model_fit)
  
  #calculate overlap of home ranges before and after encounter
  before_overlap <- overlap(UDS_before)
  after_overlap <- overlap(UDS_after)
  
  UDS_for_meta <- list(before_overlap, after_overlap)
  save(UDS_for_meta, file=  file = str_c(folder_path, "/", folder_name,"_UDS_for_meta.rda", sep = ""))
  
}

#-----------------------------------------------------------------------------------------
#Example of running the functions above with coyotes (PEC068 and PEC088) in Fig. 1

#import movement data
load("~/allmovement.rda") #change to name of your data file
ind1 = "PEC068"
ind2 = "PEC088"
study = "Wheeldon"
encounter_date = ymd_hms("2012-05-26 03:10:37") #change all info to that of corresponding encounter of interest

run_akde_analysis(data, ind1, ind2, study, encounter_date)

#-----------------------------------------------------------------------------------------
#Get physical distances over time (data to generate Fig. 1A)

#get data
ind1_dat <- data %>%
  filter(individual.local.identifier==ind1) %>%
  as.telemetry()

ind2_dat <- data %>%
  filter(individual.local.identifier==ind2) %>%
  as.telemetry()

projection(ind2_dat) <- projection(ind1_dat) #align projections

#fit WHOLE timeseries to a ctmm model
ind1_fit <- ctmm.select(ind1_dat, ctmm.guess(ind1_dat, interactive=FALSE), trace=2, cores=2)
ind2_fit <- ctmm.select(ind2_dat, ctmm.guess(ind2_dat, interactive=FALSE), trace=2, cores=2)

dists <- distances(list(ind1_dat, ind2_dat),list(ind1_fit, ind2_fit))

#-----------------------------------------------------------------------------------------
#run meta analysis to get AIC values associated with home range shift

#get meta object from folder output (loads object UDS_for_meta)
load("/uds/for/meta/object")

#must coerce overlap objects into a list format in order to get correct comparison
meta <- meta(list(before = list(UDS_for_meta[[1]]), after = list(UDS_for_meta[[2]])),
             level=0.95, col=c("red", "purple"))

meta

#this code was also used to generate the pseudoencounter plot (Fig 1D) by simply changing the encounter_date to the pseudoencounter dates
