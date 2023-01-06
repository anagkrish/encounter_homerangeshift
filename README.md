# Code and Data Associated with <<Paper name here>>

This repo contains the code and additional data required to reproduce the results of <<paper name here>>. More detail on each file is provided below:

## R Files
- encounter_homeranges.R: Code to calculate home range shifts before and after the encounter for both individuals involved in encounter. Used to generate data for Fig. 1A-D.
- encounter_bls.R: Code to calculate ballistic length scale over sliding windows for a single individual. Used to generate data for Fig. 1E.
- encounter_bears.R: Code to calculate home range shifts for bears and run population level comparisons on individuals' home ranges before and after encounters. Used to generate data for Fig. 2.

## Data Files
- bearpairs_final.csv: selected pairs that came within 100m of each other (n=29). Columns are as follows:
  - pair1: individual 1
  - pair2: individual 2
  - time: year-month-day hour-minute-second of encounter
  - dist: minimum distance between individuals at encounter (in meters)
  - trackstart: beginning of overall track for both individuals
  - trackstop: ending of overall track for both individuals.
- bhattacharyya.csv: bhattacharyya distances between home ranges before and after encounters,calculated from population-level analysis of homeranges. Columns are as follows:
  - split: before vs after encounter
  - subset: describes selected group of individuals (all encounters, encounters btwn same sex individuals, encounters btwn diff sex individuals, encounters btwn diff sex individuals in late fall)
  - low: low est for bhattacharyya distance
  - est: mean est for bhattacharyya distance
  - high: high est for bhattacharyya distance
      
*All other movement data is available on Movebank (insert dataset IDs here)
      
      
      
      