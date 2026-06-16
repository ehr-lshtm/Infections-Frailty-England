
#'load libraries
library(tidyverse);library(janitor);library(magrittr);
library(purrr);library(pryr);library(glue);library(rio);
library(arrow);library(duckdb);library(collapse)

#' STEP 1 : load exposed cohort data; obtained from Type 1 HES-APC data request; Type 1 request is for defining a cohort in HES
#' STEP 2 : load denominator data and create eligible cohort based on pre-specified eligibility in study protocol
#' STEP 3 : join cohort_eligible to cohort_exposed to import exposure variable 
#' STEP 4 : create matchable cohort; allow exposed patients to be available as controls prior to their date of first exposure
#' STEP 5 : create cohort_matched

# NOTE; Data sets are very large; so to use big data analytic tools like arrow and or duckdb
#' STEP 1 : load raw datasets; 1. exposed cohort data

#' cohort_exposed data
#' loaded from paths file

#' STEP 2 : load denominator data and create eligible cohort based on pre-specified eligibility in study protocol
#' loaded from paths file

cohort_eligible <- denominator |> # denominator dataset
  # keep observations eligible for HES linkage
  filter(patid %in% eligible_linkage_hes_apc_e) |>
  collect() |>
  transmute(patid, yob, gender, pracid,
            regstartdate=as.Date(regstartdate,format = "%d/%m/%Y"),
            regenddate=as.Date(regenddate,format = "%d/%m/%Y"),
            lcd=as.Date(lcd,format = "%d/%m/%Y"),
            cprd_ddate= as.Date(cprd_ddate,format = "%d/%m/%Y")) |>
  # Adult date: 65+ years old at study start in 2006-04-01
  mutate(birthday_65th=as.Date(paste0(yob+65, "/06/01"),format = "%Y/%m/%d")) |> 
  # Date of Eligibility for matching: latest of: : 
  # study start (2006-04-01); Was initially 2004 in protocol, we agreed to start from 2006 when ethnicity data is more complete
  # 65th birthday; or 
  # one year after practice registration 
  mutate(startdate=pmax(regstartdate+365,as.Date("2006-04-01","%Y-%m-%d"),birthday_65th)) |> 
  # Enddate: earliest of: study_end (2019-12-31); death; registration end; or no further data collected from practice.
  # lcd = Date of the most recent CPRD data collection for the practice
  mutate(enddate=pmin(regenddate, cprd_ddate, lcd, as.Date("2019-12-31","%Y-%m-%d"), na.rm = TRUE)) |>
  #Keep only people with follow-up: are eligible 
  filter(startdate < enddate) |> 
  select(patid,startdate, enddate,birthday_65th, yob, gender, pracid) |>
  mutate(across(matches("date"),~as.Date(as.character(.x)))) |>
  to_duckdb(con=my_con_integer64)

#' STEP 3 : join cohort_eligible to cohort_exposed to import exposure information
cohort_eligible <- cohort_eligible |>
  left_join(cohort_exposed,join_by(patid)) |>
  # filter indexdates that fall within study period
  filter((indexdate>=startdate & indexdate<=enddate) | is.na(indexdate)) 

#' STEP 4 : create match-able cohort; allow exposed patients 
# to be available as controls prior to their date of first exposure

# variables required in memory before running (correctly named):
# patid:          CPRD patient id [nb if using Aurum data, patid MUST be stored as a "double" precision variable]
# indexdate:      date of "exposure" for exposed patients (missing for potential controls)
# gender:         gender, numerically coded (e.g. 1=male, 2=female)
# startdate:      date of start of CPRD follow-up
# enddate:        date of end of follow-up as a potential control, generally = end of CPRD follow-up, but see "important note" below
# exposed:        indicator: 1 for exposed patients, 0 for potential controls
# yob:            year of birth
# IMPORTANT NOTE: in most cases it is desirable to allow exposed patients 
# to be available as controls prior to their date of first exposure. Such 
# patients should be included in the dataset twice (i.e. two separate rows):
# once as exposed (exposed = 1, startdate = start of CPRD follow-up, 
# indexdate = date of first exposure, enddate = end of CPRD follow-up),
# and once as a potential control (exposed = 0, startdate = start of CPRD 
# follow-up, indexdate = missing, enddate = date of first exposure-1)

create_cohort_matchable <- function(data) {
  
  cohort_matchable <- data |> 
    mutate(exposed=if_else(is.na(indexdate), 0, 1))
  
# Get exposed patients to be available as controls prior to their date of first exposure
  pre_exposure_ppl <- cohort_matchable |> 
    filter(exposed==1) |> 
    mutate(enddate=indexdate-1, indexdate=NA_Date_, exposed=0,ICD=NA_character_,spno=NA) |> 
    filter(enddate > startdate)
  
  bind_rows(cohort_matchable, pre_exposure_ppl) 
}

# create_cohort_matchable
cohort_matchable <- create_cohort_matchable(data=cohort_eligible)|>
  to_duckdb(con=my_con_integer64)

#' STEP 5 : create cohort_matched

#' Create matched cohort using sequential trials matching
#' @description each daily trial includes all n eligible people who 
#' become exposed on that day (exposed=1) and
#' a sample of n eligible controls (exposed=0) who:
#' - had not been exposed on or before that day (still at risk of becoming exposed);
#' - still at risk of an outcome (not left the study); 
#' - had not already been selected as a control in a previous trial
#' @param data A data frame with one or two rows per participant with:
#'  1. $patid: The patient ID
#'  2. $startdate: the date people become eligible 
#'  3. $indexdate: the date eligible people got exposed 
#'  4. $enddate: the date people leave the study 
#'  @param pracid_list A vector of practice IDs to iterate over 
#'  @param dayspriorreg days prior registration required for controls

create_cohort_matched <- function(dayspriorreg=0,data) {
  #Map across every practice group
  library(future)
  plan(multisession, workers = 8)
  
  pracid_list <- data |> as_tibble() |> 
    fmutate(pracid=as.character(pracid)) |> pull(pracid) |> unique()
  
  purrr::map_df(pracid_list, .progress = TRUE, \(x) {
    
    data_pracid <- data |> fsubset(pracid %in% x) 
    
    # Arrange by indexdate, and then randomly
    cohort_matchable <- data_pracid |>
      ftransform(sortunique=runif(nrow(data_pracid))) |> 
      roworder(exposed, indexdate, sortunique) |> 
      fselect(-sortunique)
    
    exposed <- cohort_matchable |> fsubset(exposed==1) |> ftransform(setid=patid)
    unexposed <- cohort_matchable |> fsubset(exposed==0)
    
    matched <- exposed[0,] #Make empty dataframe with same columns to be filled
    
    #Loop through all people that ever get exposed (each one gets matched to people who are unexposed at the same time)
    for (i in 1:nrow(exposed)) {
      
      exposed_pat <- exposed[i,]
      matchday <- exposed_pat$indexdate
      
      #Drop people that can't be matched anymore 
      #(either because they have already been matched or they have passed the study end date)
      unexposed <- unexposed |> 
        fsubset(enddate > matchday & !(patid %in% matched$patid))
      
      #Perform matching
      new <- unexposed |> 
        fsubset(gender==exposed_pat$gender) |> 
        fsubset((startdate+dayspriorreg) <= matchday) |> 
        ftransform(age_difference=abs(yob-exposed_pat$yob)) |> 
        fsubset(age_difference<=2) |> 
        roworder(age_difference) |> # the closest matches are given priority
        slice(1:5) |> 
        ftransform(setid=exposed_pat$patid,
                   indexdate=matchday) |>  #Set the indexdate for everyone to the day the exposed individual got exposed
        fselect(-age_difference)
      
      if (nrow(new)>0) matched <- bind_rows(matched, exposed_pat, new)
      
    }
    
    matched |> 
      roworder(setid, -exposed) }) # put exposed first ; thus exposed in ascending order
}

# run matched_cohort function to create matched cohort

length(pracid_list) 

# create matched cohort
cohort_matched <- create_cohort_matched(data = cohort_matchable) |>
  to_duckdb(con=my_con_integer64)


                    # FURTHER DATA REQUESTS

# TYPE-2 HES APC DATA REQUEST - 1: Type 2 data request is for requesting data for a defined cohort
# In our case we need this to redefine the cohort; determin diagnostic position of infection during hospitalisation
# 1. for hospitalisation data including diagnostic position of infection codes: exposed individuals only
 
# 2. Linked IMD data requestL for indices of multiple deprivation (individual and practice level): all individuals

# What is required: patid, pracid; NOTE: unique patids

# unique patids and pracids
cohort_matched |>
  select(patid,pracid) |>
  distinct(patid,.keep_all = T) |>
  write.table("./data/type2_linkage_request_patid.txt",quote = F,row.names = F)

# unique patids
cohort_matched_patids <- cohort_matched |>
  select(patid,pracid) |>
  distinct(patid,.keep_all = T)

length(cohort_matched_patids |> pull(patid)) 

write.table(cohort_matched_patids,"./data/cohort_matched_patids.txt",quote = F,col.names = F,row.names = F)
write.table(cohort_matched_patids,"./data/imd_request.txt",quote = F,col.names = F,row.names = F)

# unique patids: exposed ohort
cohort_matched_exposed_patids <- cohort_matched |>
  filter(exposed==1) |>
  pull(patid) |> unique()

length(cohort_matched_exposed_patids) 

write.table(cohort_matched_exposed_patids,"./data/exposed_patids.txt",quote = F,col.names = F,row.names = F)
write.table(cohort_matched_exposed_patids,"./data/hes_request.txt",quote = F,col.names = F,row.names = F)

# CREATE SECOND MATCHED COHORT
# RE-MATCHING after obtaining diagnostic position variable to fully define exposure
# Update the exposed population: Filter exposed cohort to include only individuals with an infection in the first diagnostic position 
# (select the first if more than one)
# then perform re-matching: REPEAT THE WHOLE PROCESS:
# NOTE: Only new thing is the new exposed cohort exposed_cohort_2

cohort_exposed_2 <- hes_primary_diag_hosp |>
  # keep observations eligible for HES linkage
  filter(patid %in% eligible_linkage_hes_apc_e) |> 
  rename(indexdate=admidate,ICD=ICD_PRIMARY) |> 
  # select the first diagnosis of all first diagnoses within hospitalizations
  group_by(patid) |>
  arrange(indexdate,.by_group = T) |>
  slice(1) |>
  ungroup()|>
  to_duckdb(con=my_con_integer64)

# cohort_eligible_2
cohort_eligible_2 <- collect(denominator) |> # denominator dataset
  # keep observations eligible for HES linkage
  filter(patid %in% eligible_linkage_hes_apc_e) |>
  transmute(patid, yob, gender, pracid,
            regstartdate=as.Date(regstartdate,format = "%d/%m/%Y"),
            regenddate=as.Date(regenddate,format = "%d/%m/%Y"),
            lcd=as.Date(lcd,format = "%d/%m/%Y"),
            cprd_ddate= as.Date(cprd_ddate,format = "%d/%m/%Y")) |>
  # Adult date: 65+ years old at study start in 2006-04-01
  mutate(birthday_65th=as.Date(paste0(yob+65, "/06/01"),format = "%Y/%m/%d")) |> 
  # Date of Eligibility for matching: latest of: : 
  # study start (2006-04-01); Was initially 2004 in protocol, we agreed to start from 2006 when ethnicity data is more complete
  # 65th birthday; or 
  # one year after practice registration 
  mutate(startdate=pmax(regstartdate+365,as.Date("2006-04-01","%Y-%m-%d"),birthday_65th)) |> 
  # Enddate: earliest of: study_end (2019-12-31); death; registration end; or no further data collected from practice.
  # lcd = Date of the most recent CPRD data collection for the practice
  mutate(enddate=pmin(regenddate, cprd_ddate, lcd, as.Date("2019-12-31","%Y-%m-%d"), na.rm = TRUE)) |>
  #Keep only people with follow-up: are eligible 
  filter(startdate < enddate) |> 
  select(patid,pracid,startdate, enddate,birthday_65th, yob, gender, cprd_ddate) |>
  mutate(across(matches("date"),~as.Date(as.character(.x))))  |>
  # join with new cohort_exposed_2
  full_join(collect(cohort_exposed_2),join_by(patid)) |>
  # filter indexdates that fall within study period
  filter(between(indexdate,startdate, enddate) | is.na(indexdate)) |>
  to_duckdb(con=my_con_integer64)

# Create cohort_matchable_2
# allow exposed patients to be available as controls prior to their date of first exposure

cohort_matchable_2 <- create_cohort_matchable(collect(cohort_eligible_2))|>
  to_duckdb(con=my_con_integer64)

# Create cohort_matched_2
cohort_matched_2 <- create_cohort_matched(data = collect(cohort_matchable_2)) |>
  to_duckdb(con=my_con_integer64)

# PRIMARY CARE DATA REQUEST
# all other variables : for all individuals
# what is required: unique patids

cohort_matched_2_distint_patids <- cohort_matched_2 |>
  select(patid) |>
  distinct(patid,.keep_all = T) |>
  collect()

length(cohort_matched_2_distint_patids |> pull(patid))        

cohort_matched_2_all_patids <- cohort_matched_2 |>
  select(patid)

length(cohort_matched_2_all_patids |> pull(patid)) 

