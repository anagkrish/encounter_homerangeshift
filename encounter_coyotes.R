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
  folder_dir <- "folder/name/here"
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
load("movement/data/here") #change to name of your data file
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

dists <- ctmm::distances(list(ind1_dat, ind2_dat), list(ind1_fit, ind2_fit))

#-----------------------------------------------------------------------------------------
#run meta analysis to get AIC values associated with home range shift

#get meta object from folder output (loads object UDS_for_meta)
load("/uds/for/meta/object")

#must coerce overlap objects into a list format in order to get correct comparison
meta <- meta(list(before = list(UDS_for_meta[[1]]), after = list(UDS_for_meta[[2]])),
             level=0.95, col=c("red", "purple"))

meta

#this code was also used to generate the pseudoencounter plot (Fig 1D) by simply changing the encounter_date to the pseudoencounter dates

#code to generate components of fig. 1

#A-----------
"distance/data/here" %>% 
  mutate(day=format(as.POSIXct(cut(timestamp, breaks='days')), format="%Y-%m-%d %H:%M:%S")) %>%
  group_by(day) %>%
  summarize(low=min(low, na.rm=T),
            est=min(est, na.rm=T),
            high=min(high, na.rm=T)) %>%
  mutate(days_to_encounter = difftime(day, 
                                      as.POSIXct(encounter_date, format="%Y-%m-%d %H:%M:%S"), units="days"),
         day=ymd_hms(day)) %>%
  filter(days_to_encounter > -181,
         days_to_encounter < 63) %>%
  ggplot(mapping=aes(x=day,y=est)) +
  geom_line() +
  geom_ribbon(mapping=aes(ymin=low, ymax=high),fill="grey",alpha=0.5) +
  geom_vline(mapping=aes(xintercept= encounter_date), linetype="dotdash") +
  labs(x=NULL, y = "Distance (m)") +
  scale_x_datetime(date_breaks = "1 month",
                   date_labels = "%Y\n%b")+ #manually added so they'll be aligned
  scale_y_log10() +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth=0.5),
        text=element_text(size=18))

#B/C-----------
plot("telemetry/objects/here", col = c("red","blue"),
     UD = "UD/objects/here", col.DF = "#969ca1",
     ylim = y_range, xlim = x_range, level.UD = 0.95, 
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

#D-----------
UDS_all <- #get UD objects from run_akde_analysis function

UDS_all_aligned <- ctmm:::same.grids(UDS_all)

MEAN1 <- mean(list(UDS_all_aligned[[1]], UDS_all_aligned[[2]]),sample=FALSE) # pec068 and pec088 before
MEAN2 <- mean(list(UDS_all_aligned[[3]], UDS_all_aligned[[4]]),sample=FALSE) # pec068 and pec088 after
HDIST <- sqrt(MEAN2$PDF) - sqrt(MEAN1$PDF)

lev=0.95
UD_all <- rbind(SpatialPolygonsDataFrame.UD(UDS_all[[1]], level.UD=lev),
                SpatialPolygonsDataFrame.UD(UDS_all[[2]], level.UD=lev),
                SpatialPolygonsDataFrame.UD(UDS_all[[3]], level.UD=lev),
                SpatialPolygonsDataFrame.UD(UDS_all[[4]], level.UD=lev))
UD_all <- UD_all[c(2,5,8,11),] #subset just 95% est, not low/high

diff_ras <- raster(HDIST) #hellinger distance to show change in space use
projection(diff_ras) <- ctmm::projection(UD_all)
diff_ras <- raster::setExtent(diff_ras, UD_all)
diff_ras <- raster::mask(x = diff_ras, mask = UD_all)
diff_ras <- crop(x = diff_ras, y = extent(UD_all))

breaks = round(seq(min(HDIST), max(HDIST), length.out=21), 5)
em = #merge(extent(diff_ras)+5,extent(diff_ras)+0.2)
  extent(c(-10.2, 7, -7.5, 7.5))

plot(em, lwd=0, xlab="x (kilometers)", ylab = "y (kilometers)", legend=FALSE,
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(diff_ras, col=pals::coolwarm(21), legend=F,
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, add=T)
plot(UD_all, col=NA, border=c("black", "black", "darkgrey", "darkgrey"), lwd=2, add=T) #after
plot(diff_ras, breaks=breaks, 
     col=pals::coolwarm(21),
     legend.only=TRUE, horizontal = T,
     legend.width=0.75, legend.shrink=0.4,
     smallplot=c(0.15,.9, .95,.97),
     #c(min % from left, max % from left, min % from bottom, max % from bottom)
     axis.args=list(at=breaks,
                    labels=format(breaks, scientific=F),
                    cex.axis=0.6),
     legend.args = list(text="Difference of Mean Coyote UDs Before/After Encounter", side = 1, line=-1.9),
     cex.lab=3,
     add=T)

#E-----------
#calculate pseudoencounters and encounters with above described method
"between/individual/pseudoencounter/data/here" %>% 
  mutate(date = as.POSIXct(date, format="%Y-%m-%d %H:%M:%S"),
         days_to_encounter = difftime(date, as.POSIXct(encounter_date, format="%Y-%m-%d %H:%M:%S"), units="days"),
         encounter=ifelse(round(days_to_encounter)==0, "yes","no"),
         shapecolor=factor(paste(`before/after`,encounter), 
                           levels=c("before no", "after no", "before yes", "after yes"))) %>%
  ggplot(mapping = aes(x = days_to_encounter, y = est, color = shapecolor, shape = shapecolor)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  geom_vline(xintercept = 0, linetype="dotdash") + 
  labs(x = "Days to Encounter", y = "Range Distribution Overlap \nBetween Individuals") +
  scale_color_manual(name=" ",
                     labels=c('Before Null Encounter', 'After Null Encounter', 
                              'Before Encounter', 'After Encounter'),
                     values = c("purple","orange", "purple", "orange")) +
  scale_shape_manual(name=" ",
                     labels=c('Before Null Encounter', 'After Null Encounter', 
                              'Before Encounter', 'After Encounter'),
                     values= c(1, 1, 16, 16)) +
  scale_x_continuous(limits=c(-181, 63),
                     breaks = c(-240, -210, -180, -150, -120, -90, -60, -30, 0, 30, 60)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth=0.5),
        text=element_text(size=18),
        legend.text = element_text(size = 18),
        egend.box.margin = margin(0, 0, 0, 0),
        legend.direction = "horizontal",
        legend.position = "top",
        legend.justification = c(0.31,0.2))

#F-----------
#calculate pseudoencounters and encounters with above described method
"within/individual/pseudoencounter/data/here" %>% 
  mutate(dt = as.POSIXct(dt, format="%Y-%m-%d"),
         days_to_encounter = round(difftime(dt, as.POSIXct(encounter_date, format="%Y-%m-%d %H:%M:%S"), units="days"))) %>%
  mutate_at(c("overlap_low", "overlap_est", "overlap_high"), parse_number) %>%
  #encounter= factor(encounter)) %>%
  ggplot(mapping = aes(x = days_to_encounter, y = overlap_est, color=ind)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = overlap_low, ymax = overlap_high)) +
  geom_vline(xintercept = 0, linetype="dotdash") + 
  labs(x = "Days to Encounter", 
       y = "Range Distribution Overlap \nWithin Individuals") +
  scale_y_continuous(limits=c(0.65,1)) +
  scale_x_continuous(limits=c(-181, 63),
                     breaks = c(-240, -210, -180, -150, -120, -90, -60, -30, 0, 30, 60)) +
  scale_color_manual(values=c("red","blue")) +
  theme_bw() + 
  #facet_wrap(~ind, ncol=1) +
  labs(color=NULL) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth=0.5),
        text=element_text(size=18), legend.position = c(0.89, 0.12),
        legend.key = element_blank(), legend.background=element_blank())


#G-----------
#calculate bls from encounter_bls.R file

pec068_bls %>%
  filter(days_to_encounter > -181,
         days_to_encounter < 63) %>%
  ggplot(mapping=aes(x=days_to_encounter, y=b_l_s_mid)) +
  geom_line(aes(color="PEC068")) +
  geom_line(mapping=aes(x=days_to_encounter, y=b_l_s_mid, color="PEC088"), data=pec088_bls) +
  geom_ribbon(aes(ymin = b_l_s_high, ymax = b_l_s_low), alpha = 0.2) +
  geom_ribbon(aes(ymin = b_l_s_high, ymax = b_l_s_low), alpha = 0.2, data=pec088_bls) +
  geom_vline(mapping=aes(xintercept=0), linetype="dotdash") +
  scale_x_continuous(expand = c(0, 0), limits = c(-184, 65)) +
  labs(x="Days to Encounter", y = "Ballistic Length Scale (m)") +
  scale_color_manual(name=NULL,
                     breaks=c('PEC068', 'PEC088'),
                     values=c('PEC068'='red', 'PEC088'='blue')) +
  scale_x_continuous(limits=c(-181, 63),
                     breaks = c(-240, -210, -180, -150, -120, -90, -60, -30, 0, 30, 60)) +
  theme_bw() + #get rid of background
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth=0.5),
        text=element_text(size=18), legend.position = c(0.89, 0.12), 
        legend.key = element_blank(), legend.background=element_blank())
