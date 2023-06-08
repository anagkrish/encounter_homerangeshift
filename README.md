# Intraspecific encounters can induce home-range shifts: Associated R Scripts

Zenodo archive v1.0.0 (for initial review/preprint)

[![DOI](https://zenodo.org/badge/586006515.svg)](https://zenodo.org/badge/latestdoi/586006515)

<b>Authors:</b> William F. Fagan<sup>1,*</sup>, Ananke Krishnan<sup>1§</sup>, Qianru Liao<sup>1§</sup>, Christen H. Fleming<sup>1,2,3</sup>, Daisy Liao<sup>1</sup>, Clayton Lamb<sup>4</sup>, Brent Patterson<sup>5</sup>, Tyler Wheeldon<sup>5</sup>, Ricardo Martinez-Garcia<sup>3,6</sup>, Jorge F. S. Menezes<sup>3</sup>, Michael J. Noonan<sup>7</sup>, Eliezer Gurarie<sup>8</sup>, and Justin M. Calabrese<sup>1,3,9</sup>

<sup>1</sup>Department of Biology, University of Maryland, College Park, MD, United States

<sup>2</sup>Smithsonian Conservation Biology Institute, National Zoological Park, Front Royal, VA, United States 

<sup>3</sup>Center for Advanced Systems Understanding (CASUS), Helmholtz-Zentrum Dresden-Rosendorf (HZDR), Görlitz, Germany

<sup>4</sup>Department of Biology, University of British Columbia, Kelowna, BC, Canada

<sup>5</sup>Ontario Ministry of Natural Resources and Forestry, Trent University, Peterborough, ON, Canada

<sup>6</sup>ICTP - South American Institute for Fundamental Research & Instituto de Física Teórica, Universidade Estadual Paulista - UNESP, São Paulo, SP, Brazil.

<sup>7</sup>Department of Biology, The University of British Columbia Okanagan, Kelowna, BC, Canada

<sup>8</sup>Department of Environmental Biology, SUNY Environmental Science and Forestry, Syracuse, NY, United States

<sup>9</sup>Department of Ecological Modelling, Helmholtz Centre for Environmental Research-UFZ, Leipzig, Germany


<sup>*</sup>Correspondence Author

<sup>§</sup>Joint Second Authors

This repo contains the code and additional data required to reproduce the results in this manuscript. More detail on each file is provided below:

## R Files
- encounter_homeranges.R: Code to calculate home range shifts before and after the encounter for both individuals involved in encounter. Used to generate data for Fig. 1A-D.
- encounter_bls.R: Code to calculate ballistic length scale over sliding windows for a single individual. Used to generate data for Fig. 1E.
- encounter_bears.R: Code to calculate home range shifts for bears and run population level comparisons on individuals' home ranges before and after encounters. Used to generate data for Fig. 2.
- supplementalfigs.R: Code to generate supplemental figures from data files (code to generate data files shown in encounter_bears.R). Fig. S3/S4 waiting on data owner permission.

## Data Files
- bearpairsto500m.csv: all pairs that came within up to 500m of each other (n=103). Columns are as follows:
  - bear1_id: individual 1
  - bear2_id: individual 2
  - time: year-month-day hour-minute-second of encounter
  - dist: minimum distance between individuals at encounter (in meters)
  - trackstart: beginning of overall track for both individuals
  - trackstop: ending of overall track for both individuals.
- bhattacharyya.csv: bhattacharyya distances between home ranges before and after encounters within 500m,calculated from population-level analysis of homeranges. Columns are as follows:
  - split: before vs after encounter
  - subset: describes selected group of individuals (all encounters, encounters btwn same sex individuals, encounters btwn diff sex individuals, encounters in late fall, encounters in spring, encounters btwn diff sex individuals in late fall)
  - low: low est for bhattacharyya distance
  - est: mean est for bhattacharyya distance
  - high: high est for bhattacharyya distance
- percentoverlap500m.csv: population-level overlap between home ranges after encounters up to 500m as a percent of the overlap of home ranges before encounter. 
  - subset: describes selected group of individuals (all encounters, encounters btwn same sex individuals, encounters btwn diff sex individuals, encounters in late fall, encounters in spring, encounters btwn diff sex individuals in late fall)
  - low: low est for percent overlap (after/before)
  - est: mean est for percent overlap (after/before)
  - high: high est for percent overlap (after/before)
- uds_btwninds.csv: UDS percent overlap for each individual in pair before and after all encounters (n=103). Shows change between overlap of both individuals' home ranges before and after encounter. 
- uds_withininds.csv: UDS percent overlap for each individual in pair before and after encounters up to 100m (n=44) . Shows change within overlap of individual home ranges for each individual in pair before and after encounter.
      
*Movement data is available on Movebank.org data as datasets 1614661371 (coyotes) and 1044288582 (grizzly bears).

      