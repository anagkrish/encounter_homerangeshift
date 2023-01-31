# Code and Data Associated with <<Paper name here>>

This repo contains the code and additional data required to reproduce the results of <<paper name here>>. More detail on each file is provided below:

## R Files
- encounter_homeranges.R: Code to calculate home range shifts before and after the encounter for both individuals involved in encounter. Used to generate data for Fig. 1A-D.
- encounter_bls.R: Code to calculate ballistic length scale over sliding windows for a single individual. Used to generate data for Fig. 1E.
- encounter_bears.R: Code to calculate home range shifts for bears and run population level comparisons on individuals' home ranges before and after encounters. Used to generate data for Fig. 2.

## Data Files
- bearpairs100m.csv: selected pairs that came within 100m of each other (n=35). Columns are as follows:
  - pair1: individual 1
  - pair2: individual 2
  - time: year-month-day hour-minute-second of encounter
  - dist: minimum distance between individuals at encounter (in meters)
  - trackstart: beginning of overall track for both individuals
  - trackstop: ending of overall track for both individuals.
- bearpairs500m.csv: same as bearpairs100m but for all pairs that came within up to 500m of each other.
- bhattacharyya.csv: bhattacharyya distances between home ranges before and after encounters within 100m,calculated from population-level analysis of homeranges. Columns are as follows:
  - split: before vs after encounter
  - subset: describes selected group of individuals (all encounters, encounters btwn same sex individuals, encounters btwn diff sex individuals, encounters in late fall, encounters in spring, encounters btwn diff sex individuals in late fall)
  - low: low est for bhattacharyya distance
  - est: mean est for bhattacharyya distance
  - high: high est for bhattacharyya distance
- percentoverlap100m.csv: population-level overlap between home ranges after encounters within 100m as a percent of the overlap of home ranges before encounter. 
  - subset: describes selected group of individuals (all encounters, encounters btwn same sex individuals, encounters btwn diff sex individuals, encounters in late fall, encounters in spring, encounters btwn diff sex individuals in late fall)
  - low: low est for percent overlap (after/before)
  - est: mean est for percent overlap (after/before)
  - high: high est for percent overlap (after/before)
- percentoverlap500m.csv: same as percentoverlap100m.csv but for encounters within 500m.
- udsvalsbearsbtwn100.csv: UDS percent overlap for each pair (n=35) of encounters with 100m before and after the encounter. Shows change between overlap of both individuals' home ranges before and after encounter. 
- udsvalsbearswin100.csv: UDS percent overlap for each individual in pair (n=35) of encounters with 100m before and after the encounter. Shows change within overlap of individual home ranges for each individual in pair before and after encounter. 
      
*All other movement data is available on Movebank (insert dataset IDs here)
      
      
      
      