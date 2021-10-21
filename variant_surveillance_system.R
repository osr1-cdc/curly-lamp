options(java.parameters = "-Xmx8000m")
# Background

#Emerging variants of SARS-CoV-2 with increasing share among cases are of potential public health concern. 
# A survey design to track variants, based on the sequencing data residing at CDC and weighted to adjust 
# for multiplicity of sources as well as variations in sequencing rates, is discussed here.

# Data and methods

#An event recorded in the surveillance system is the collection of a sample for (PCR) testing (labelled by 
# date and place of collection) that results in a positive test and is selected for sequencing and 
# submission to the system. Prevalence of variants is estimated from sequence data from multiple sources, 
# weighted to adjust for variations in sequencing and testing rates over time, source and jurisdiction, 
# using a design-based approach. In this analysis:
# * The event time is date of collection of sample for testing
# * The event location is currently set to the state reporting the test result
# * The protocol for selection of positive test samples for sequencing is assumed to be random, as the 
#   commercial laboratories are unable to report their selection protocol, and the NS3 instructions neglect 
#   to specify any
# * The strata are U.S. states and territories.
# * The clusters are the data streams: CDC NS3 surveillance, and (currently) three contractors -- Laboratory 
#   Corporation of America, Quest Diagnostics, and Helix. 
# * Inverse-probability-of-selection weights are generated for estimates of prevalence among those sequenced 
#   (unweighted), among test positives, as well as among all infections.

## Current data streams and workflow

#The current workflow is a provisional implementation of the methods. It begins with tapping the current 
# data streams:
# * Survey sample data: Genomic data (lineage and S gene mutation list), by collection date, reporting 
#   state, data source (lab, NS3, etc.) and any known sampling bias (targeted oversampling of "S gene 
#   target failure" specimens, etc.), from a server internal to CDC. Other variables available that may 
#   be used for stratified analysis or more careful weighting, but not used currently, are demographic 
#   information (limited completeness) on the individual tested, and ZIP code of specimen collection. 
# * Supplementary sample data: Pangolin lineage data table, from the CDC internal server, for current 
#   lineage information
# * Weighting data: Two tables, updated on or close to the date of analysis, are downloaded from HHS Protect: 
#   count of PCR tests results (positive, negative, etc.), aggregated by collection date and reporting state, 
#   and count of PCR tests results (positive, negative, etc.) from selected labs (currently lab names that 
#   contain substrings "labcorp", "laboratory corporation", "helix", "illumina", "quest"), aggregated by 
#   collection date, reporting state and lab name. 
# * Supplementary weighting data: state population estimates (ACS 2018 5 year set)

# ver: 1Sep202

#--------------- Packages Filepaths and Settings ---------------
library(implyr)
library(odbc)
library(survey)
library(dplyr)
library(RJDBC)
library(optparse)

# optparse option list
option_list <- list(
  make_option(c("-u", "--user"), type = "character", default= NA,
              help="User Name",
              metavar = "character"),
  make_option(c("-p", "--password"), type = "character", default= NA,
              help="Password",
              metavar = "character")
  
)

# parseing options list
opts <- parse_args(OptionParser(option_list = option_list))

# get base directory and driver location
initial.options = commandArgs(trailingOnly = FALSE)
file.arg.name = "--file="
script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename = dirname(script.name)
# correction for if running from base repo interactively
if(length(script.basename) == 0) {
  script.basename = "."
}
# create data dir
dir.create(paste0(script.basename,"/data"))

jdbc_driver = paste0(script.basename, "/jdbc/ClouderaImpalaJDBC-2.6.20.1024/ClouderaImpalaJDBC41-2.6.20.1024/")

options(survey.adjust.domain.lonely=T,survey.lonely.psu="average")

#load("variant_survey_dat2021-08-06.RData")

make.svy.dat = !(paste0("svydat_", Sys.Date(), ".RData") %in% list.files("surveillance"))

#define a not in logical statement
`%notin%`=Negate(`%in%`)

#variable identifying which Hadoop table to read in. Options are:
# ~ "deduplication_cdcncbigisaid_auto": this data table is updated regularly but may be subject to cleaning issues
# ~  "analytics_metadata_frozen": this data table is updated less frequently, but is the cleanest version of the sequence data
# NOTE - if running the official Friday analysis use the "analytics_metadata_frozen" table
seq_table = "sc2_archive.analytics_metadata_frozen"
node="11" #you can use any node from 08 - 13. Sometimes nodes fail so if you get an error you can try another node
#-------------- Import Data ---------------

## Get the sequence data

impala_classpath <- list.files(path = jdbc_driver,
                               pattern = "\\.jar$", full.names = TRUE)
.jinit(classpath = impala_classpath)

drv <- JDBC(
  driverClass = "com.cloudera.impala.jdbc41.Driver",
  classPath = impala_classpath,
  identifier.quote = "`"
)

impala <- dbConnect(drv, paste0("jdbc:impala://cdp-", node, ".biotech.cdc.gov:21050/sc2_air;AuthMech=3;useSasl=1;SSL=1;AllowSelfSignedCerts=1;CAIssuedCertNamesMismatch=1"),
                     opts$user, opts$password)

# Get names if needed: tbls = dbListTables(impala)
all.vars = dbGetQuery(impala, 'DESCRIBE sc2_archive.analytics_metadata_frozen')[,1]
if(seq_table=="sc2_archive.analytics_metadata_frozen"){
get.vars = grep("(^csid)|(^primary)|(^covv)|(^ct)|(zip$)|(targeted)|(vendor)", all.vars, value=TRUE)
#drop the contractor_vendor_id
get.vars = c(get.vars, "age") } else {
get.vars= c( "primary_virus_name","primary_nt_id","primary_country" ,              
             "primary_state","primary_state_abv","primary_host" ,                 
             "primary_specimen","primary_collection_date_dt","csid" ,                         
             #"contractor_vendor_id",
             "contractor_vendor_name","contractor_targeted_sequencing",
             "covv_accession_id","covv_orig_lab","covv_add_host_info" ,           
             "primary_collection_date")#,"patient_age_cl","covv_patient_age")
}

# Subset to US sequences for faster download
query = paste('SELECT', paste0('A.', get.vars, collapse=', '),
              paste0(' FROM sc2_archive.analytics_metadata_frozen as A
INNER JOIN
(SELECT max(date_frozen) as max_frozen
  FROM sc2_archive.analytics_metadata_frozen) as M
ON A.date_frozen = M.max_frozen'),'WHERE A.primary_country in ("United States", "USA")') # if unavailable, use test_deduplication_cdcncbigisaid_auto for testing
dat = dbGetQuery(impala, query) #saves sequencing data as data frame in R- includes GISAID data

# Current Pangolin lineages:
#pangolin = dbReadTable(impala, "sc2_src.pangolin") # Lineage master list
pangolin = dbGetQuery(impala, 'SELECT distinct * FROM sc2_src.pangolin') # 
# S gene mutation lists, and source (for NS3 and labs)
baseline = dbGetQuery(impala, 'SELECT nt_id, source, primary_virus_name, s_mut, collection_date FROM sc2_dev.baselineseq') # state, zip available, also in dedup; S1 slower, built at query
#baseline = dbGetQuery(impala, "SELECT nt_id, source, primary_virus_name, collection_date FROM sc2_dev.baselineseq") # state, zip available, also in dedup; S1 slower, built at query
tests = dbGetQuery(impala, 'SELECT to_date(H.collection_date) as collection_date,
  H.reporting_state,
  H.INDETERMINATE,
  H.INVALID,
  H.NEGATIVE,
  H.POSITIVE
  FROM sc2_archive.hhs_protect_testing_frozen as H
  INNER JOIN
  (SELECT max(date_frozen) as max_date_frozen
    FROM sc2_archive.hhs_protect_testing_frozen
  ) as F
  ON H.date_frozen = F.max_date_frozen
  WHERE H.collection_date is NOT NULL')

colnames(tests) = c("collection_date", "reporting_state", "INDETERMINATE", "INVALID" , "NEGATIVE" , "POSITIVE")

dbDisconnect(impala)

pops = read.delim(paste0(script.basename,"/resources/ACStable_B01001_40_2018_5.txt"))[, c("STUSAB", "Total.")] #population data may need to be updated

save(dat, pangolin, baseline, tests, pops, file=paste0(script.basename, "/data/", "variant_survey_dat",Sys.Date(),".RData"))#"surveillance/variant_survey_dat.RData")

#get capture of just analytics_metadata_frozen
#write.csv(dat, paste0("Analytics_data_",Sys.Date(),"_MS.csv"),row.names=FALSE)


#--------------- Data Cleaning ---------------

#The genomic data are subset to U.S. specimens among human hosts. Some light cleaning of the sequence data includes dropping records where the state abbreviation is incorrect, or the collection data is prior to October 2019. An HHS Region variable is merged in.

## Some general parameters:
week0day1 = as.Date("2020-01-05") # First Sunday in 2020
# HHS regions
HHS_reg = list(HHS1=c("CT", "ME", "MA", "NH", "RI", "VT"),
               HHS2=c("NJ", "NY", "PR", "VI"),
               HHS3=c("DE", "DC", "MD", "PA", "VA", "WV"),
               HHS4=c("AL", "FL", "GA", "KY", "MS", "NC", "SC", "TN"),
               HHS5=c("IL", "IN", "MI", "MN", "OH", "WI"),
               HHS6=c("AR", "LA", "NM", "OK", "TX"),
               HHS7=c("IA", "KS", "MO", "NE"),
               HHS8=c("CO", "MT", "ND", "SD", "UT", "WY"),
               HHS9=c("AZ", "CA", "HI", "NV", "AS", "MP", "GU", "MH"),
               HHS10=c("AK", "ID", "OR", "WA")
)
# Territory 2020-07-01 populations from https://en.wikipedia.org/wiki/List_of_states_and_territories_of_the_United_States_by_population (2021-04-19)
# Marshall Islands 2018 estimate https://en.wikipedia.org/wiki/Marshall_Islands (2021-04-19)
pops = rbind(pops, data.frame(STUSAB=c("AS", "GU", "MP", "VI", "MH"), `Total.`=c(49437, 168485, 51433, 106235, 58413)))
hhs = data.frame(STUSAB=toupper(pops$STUSAB))
hhs$HHS = sapply(hhs$STUSAB, grep, HHS_reg)

## Subset by time and place: USA, states with two letter abbreviations, human host, reasonable collection date
# U.S., legitimately coded states, human host; dat may already be subset to US, but retained for backward compatibility
us.dat = subset(dat, primary_country %in% c("United States", "USA") & primary_host=="Human") 
us.dat = subset(us.dat, nchar(primary_state_abv)==2) 
## Get unique records in Hadoop table (2021-04-26)
us.dat = distinct(us.dat)
pangolin = distinct(pangolin)
baseline = distinct(baseline)
# Disambiguate and remove unreasonable dates
us.dat$collection_date = as.Date(us.dat$primary_collection_date)
us.dat$collection_date = with(us.dat, as.Date(ifelse(is.na(collection_date), as.Date(primary_collection_date_dt), collection_date), origin="1970-01-01"))
if ("covv_subm_date" %in% names(us.dat)) us.dat$covv_subm_date = as.Date(us.dat$covv_subm_date)
if ("contractor_receive_date_to_cdc" %in% names(us.dat)) us.dat$contractor_receive_date_to_cdc = as.Date(us.dat$contractor_receive_date_to_cdc)
us.dat$yr_wk = as.character(us.dat$collection_date - as.numeric(strftime(us.dat$collection_date, format="%w"))) 
us.dat$DAY = as.numeric(us.dat$collection_date - week0day1) 
us.dat = subset(us.dat, collection_date >= as.Date("2019-10-01"))
# drops any smaples collected prior to 10/1/2020
us.dat = merge(us.dat, pangolin[, c("nt_id", "lineage")], by.x="primary_nt_id", by.y="nt_id", all.x=TRUE)
if ("covv_lineage" %in% names(us.dat)) us.dat$lineage = with(us.dat, ifelse(is.na(lineage), covv_lineage, lineage))#this line of code only works if using the dedupe table
# NS3 identifier in baseline$source:
us.dat = merge(us.dat, baseline[, c("nt_id", "source", "primary_virus_name", "s_mut")],
              by.x=c("primary_nt_id", "primary_virus_name"), by.y=c("nt_id", "primary_virus_name"), all.x=TRUE)
#us.dat = merge(us.dat, baseline[, c("nt_id", "source", "primary_virus_name")],
#                             by.x=c("primary_nt_id", "primary_virus_name"), by.y=c("nt_id", "primary_virus_name"), all.x=TRUE)
## Add in HHS regions (2021-03-09)
us.dat = merge(us.dat, hhs, by.x="primary_state_abv", by.y="STUSAB", all.x=TRUE)

#--------------- Aggregate Test Data ---------------
#Test data are aggregated by week for weighting.


## Testing data by state and date (HHS Protect) for weighting
# Done: tests = read.csv("surveillance/tests_by_state_collection_date.csv")
names(tests) = gsub("^[^a-zA-Z]+", "", names(tests)) # Odd first name
names(tests) = gsub("^[[:print:]]+state$", "STUSAB", names(tests)) # Clean an odd name
#tests$collection_date = as.Date(tests$collection_date)
tests$collection_date = as.Date(tests$collection_date,"%Y-%m-%d")
tests = subset(tests, collection_date >= as.Date("2019-10-01") & collection_date <= Sys.Date())
tests$yr_wk = as.character(tests$collection_date - as.numeric(strftime(tests$collection_date, format="%w"))) 
tests$TOTAL = rowSums(tests[, c("INDETERMINATE", "INVALID", "NEGATIVE", "POSITIVE")], na.rm=TRUE)
tests$POSITIVE = ifelse(is.na(tests$POSITIVE), 0, tests$POSITIVE)
# Aggregate by week for weighting
tests_wk = expand.grid(STUSAB=unique(tests$STUSAB), yr_wk=unique(tests$yr_wk), stringsAsFactors=FALSE)
for (cc in c("POSITIVE", "TOTAL")) {
  tests_wk = merge(tests_wk, data.frame(xtabs(tests[, cc] ~ STUSAB + yr_wk, tests)), all.x=TRUE)
  names(tests_wk) = gsub("[Ff]req", cc, names(tests_wk)) 
}


#---------------Weighting---------------

#Weights are estimated by treating each source (contracting lab, or CDC surveillance) as a cluster 
#   nested with strata (states), calculated for each week (i.e., to ensure representation weekly). 
#Two sets of weights are estimated: the first, $w_p$, for representation among (PCR) test positive 
#   indivduals; and, second, $w_i$, for representation among all prevalent infections. Assumptions 
#   (modifiable as data become available) in estimating weights are:
  
#  * Each positive (PCR) test is in the sampling frame of one of the source streams, so that for 
#    each state, week and source
#$$w_p = \frac{\mbox{number of positive PCR test results}}{\mbox{number of sequences submitted}}$$
#  * Oversampling of SGTF samples by one source results in a reduction in weights of SGTF sequences 
#   from that source by a factor that is estimated using a logistic regression model relating the 
#   odds of finding an "SGTF variant" by source, state and week.
#   There is no reliable and precise method for this yet, so these weights are subject to 
#   considerable uncertainties. Here, a [strategy](https://www.medrxiv.org/content/10.1101/2020.10.07.20208504v2.full-text) 
#   based on test positivity is used (other [strategies](https://covid19-projections.com/estimating-true-infections-revisited/) 
#   are available, too):
#  $$\frac{\mbox{number of prevalent infections}}{\mbox{number of test positives}} = \sqrt{\frac{\mbox{population of jurisdiction}}{\mbox{number of tests}}}$$
#  If each source (lab) stream is assumed to sample from a base population with the same prevalence 
#   of infection as the jurisdiction (state), it can be shown that the weight specific to each source, 
#   based on the test positivity of each source, is:
#  $$w_i = \frac{\mbox{number of prevalent infections in source}}{\mbox{number of source positives}}
#= \frac{\mbox{number of source positives}}{\mbox{number of source tests}}
#{\frac{\sqrt{\mbox{population of jurisdiction}\times\mbox{number of tests}}}{\mbox{number of positives}}}$$
  
#  The "infection" weight for each sequence is $w_p\times w_i$, and depends on stratum (state), time (week of collection) and cluster (source).

## Test tally denominator streams
test_tallies_wk = merge(merge(tests_wk, hhs, all.x=TRUE), data.frame(STUSAB=toupper(pops$STUSAB), state_population=pops$Total.), all.x=TRUE)
#estimates of infections based on: https://www.medrxiv.org/content/10.1101/2020.10.07.20208504v2.full-text
#sqrt(state_population/TOTAL) this estimates the geometric mean of the total population versus the state population
#to try to get around the bias of assuming either total or state population
test_tallies_wk = within(test_tallies_wk, INFECTIONS <- ifelse(TOTAL>0, POSITIVE * sqrt(state_population/TOTAL), 0))
test_tallies_wk$week = as.numeric(as.Date(test_tallies_wk$yr_wk) - week0day1)%/%7
incidence_by_region = merge(
  aggregate(INFECTIONS ~ yr_wk + HHS, data=test_tallies_wk, sum),
  aggregate(state_population ~ yr_wk + HHS, data=test_tallies_wk, sum)
)
incidence_by_region$HHS_INCIDENCE = incidence_by_region$INFECTIONS/incidence_by_region$state_population

## SGTF over-sampling weights -> only correcting for the Helix upsampling (other labs were not specifically SGTF up-sampling)
sgtf.1 = table(us.dat$lineage, paste(us.dat$contractor_vendor_name, us.dat$contractor_targeted_sequencing)) #lineages by contractor name
sgtf.vars = rownames(sgtf.1)[ sgtf.1[, grep("dropout$", colnames(sgtf.1))] > sgtf.1[, grep("Illumina $", colnames(sgtf.1))] ]
# Smoothed weights: use set of states and weeks where targeted samples were sequenced
sgtf.sub = unique(subset(us.dat, contractor_targeted_sequencing %in% "Screened for S dropout")[, c("primary_state_abv", "yr_wk")])
# [2021-08-31]: if () {} else {} wrapper added to handle data without SGTF flag:
if (nrow(sgtf.sub)>0) {
#linear model that fits sgtf lineage weights
sgtf.glm = glm(
  I(lineage %in% sgtf.vars)~ I(contractor_vendor_name %in% "Helix/Illumina") + primary_state_abv + yr_wk,
  family="binomial",
  subset(us.dat, yr_wk %in% sgtf.sub$yr_wk & primary_state_abv %in% sgtf.sub$primary_state_abv))
sgtf.sub$sgtf_weights = predict(sgtf.glm, cbind(contractor_vendor_name=c("Helix/Illumina"), sgtf.sub), type="response")/
  predict(sgtf.glm, cbind(contractor_vendor_name=c("Other"), sgtf.sub), type="response")
} else {sgtf.sub$sgtf_weights = numeric(0)}
# End modified code (2021-08-31)
test_tallies_wk = merge(test_tallies_wk, sgtf.sub, by.x=c("STUSAB", "yr_wk"), by.y=c("primary_state_abv", "yr_wk"), all.x=TRUE)
test_tallies_wk$sgtf_weights = ifelse(is.na(test_tallies_wk$sgtf_weights), 1, test_tallies_wk$sgtf_weights) # Inconsequential, but avoids disruptions due to NA
##!!!!!!!!!!!! 8/6/2021: bandaid fix for no SGTF flags
#test_tallies_wk$sgtf_weights =1

## Creating a trimmed down survey dataset
# Removing submission/receive dates for now, for consistency with frozen dataset
svy.dat = data.frame(STUSAB=us.dat$primary_state_abv, HHS=us.dat$HHS, yr_wk=us.dat$yr_wk, DAY=us.dat$DAY,
                     # SUBM_DT=us.dat$covv_subm_date, CDC_DT = us.dat$contractor_receive_date_to_cdc,
                     #AGE=us.dat$age_group,
                     LAB=toupper(us.dat$source),
                     SGTF_UPSAMPLING=(us.dat$contractor_targeted_sequencing %in% "Screened for S dropout"), SOURCE=us.dat$source,
                     VARIANT=us.dat$lineage, S_MUT=us.dat$s_mut)
svy.dat$LAB[is.na(svy.dat$LAB)] = "OTHER"
svy.dat$week = as.numeric(as.Date(svy.dat$yr_wk) - week0day1)/7
svy.dat = merge(svy.dat, data.frame(STUSAB=toupper(pops$STUSAB), state_population=pops$Total.), all.x=TRUE)#merge in state population

#---------------Cleaning state tagged sequences----------------------
unique(svy.dat$LAB)
svy.dat$count=1
aggregate(count~LAB,data=svy.dat,FUN=sum)

#clean LAB names
svy.dat$LAB2=svy.dat$LAB

svy.dat[svy.dat$LAB %in% unique(c(sort(unique(svy.dat$LAB)[grep(ignore.case=T,"Maryland",unique(svy.dat$LAB))]),unique(svy.dat$LAB)[grep("MD",unique(svy.dat$LAB))])),"LAB2"] <- "MD-DPH"
svy.dat[svy.dat$LAB %in% unique(c(sort(unique(svy.dat$LAB)[grep(ignore.case=T,"New Jersey",unique(svy.dat$LAB))]),unique(svy.dat$LAB)[grep("NJ",unique(svy.dat$LAB))])),"LAB2"] <- "NJ-DPH"
svy.dat[svy.dat$LAB %in% unique(c(sort(unique(svy.dat$LAB)[grep(ignore.case=T,"Texas",unique(svy.dat$LAB))]),unique(svy.dat$LAB)[grep("TX",unique(svy.dat$LAB))])),"LAB2"] <- "TX-DPH"
svy.dat[svy.dat$LAB %in% sort(unique(svy.dat$LAB)[grep(ignore.case=T,"Centers for Disease Control and Prevention",unique(svy.dat$LAB))]),"LAB2"] <- "CDC"
svy.dat[svy.dat$LAB %in% sort(unique(svy.dat$LAB)[grep(ignore.case=T,"LSU",unique(svy.dat$LAB))]),"LAB2"] <- "LSU LAB"
svy.dat[svy.dat$LAB %in% sort(unique(svy.dat$LAB)[grep(ignore.case=T,"Orange",unique(svy.dat$LAB))]),"LAB2"] <- "Orange County PHL"
svy.dat[svy.dat$LAB %in% sort(unique(svy.dat$LAB)[grep(ignore.case=T,"LAURING LAB",unique(svy.dat$LAB))]),"LAB2"] <- "LAURING LAB, UNIVERSITY OF MICHIGAN"
svy.dat[svy.dat$LAB %in% sort(unique(svy.dat$LAB)[grep(ignore.case=T,"HELIX",unique(svy.dat$LAB))]),"LAB2"] <- "HELIX"
svy.dat[svy.dat$LAB %in% sort(unique(svy.dat$LAB)[grep(ignore.case=T,"SOUTH DAKOTA",unique(svy.dat$LAB))]),"LAB2"] <- "SD-DPH"
svy.dat[svy.dat$LAB %in% sort(unique(svy.dat$LAB)[grep(ignore.case=T,"MOUNES",unique(svy.dat$LAB))]),"LAB2"] <- "OMEGA DIAGNOSTICS AT MOUNES"

aggregate(count~LAB2,data=svy.dat,FUN=sum)

#write.csv(aggregate(count~STUSAB,data=svy.dat,FUN=sum),"counts_by_state_without_vendor_id.csv",row.names=F)
#check which state-level sources have less than 100 samples
check_count <- aggregate(count~LAB2,data=svy.dat,FUN=sum)
low_lab <- check_count$LAB2[check_count$count<100]

#Drop sequences from sources with less than 100 samples or were listed as CDC sequenced
`%notin%` <- Negate(`%in%`)
dim(svy.dat)
svy.dat <- svy.dat[svy.dat$LAB2 %notin% c("CDC",low_lab),]

#aggregate state-level sources that have small numbers of samples and can be aggregated (i.e., county-level sources aggregated to state)
#svy.dat[svy.dat$LAB %in% sort(unique(svy.dat$LAB)[grep("Illinois",unique(svy.dat$LAB))]),"LAB2"] <- "IL-DPH"


#--------------- Create survey weights ---------------

#Three separate survey designs are used in this analysis:
  
# * Unweighted, to estimate variant prevalence among sequenced samples
# * Weighted for estimation among test positives
# * Weighted for estimation among infections, but using strategies of unproven reliability in this context


## Weights
# Infections to test-positive ratio
# Geometric mean strategy: https://www.medrxiv.org/content/10.1101/2020.10.07.20208504v2.full-text
#   infections/test positives = sqrt(population/number tested)
svy.dat = merge(svy.dat, cbind(test_tallies_wk[, c("STUSAB", "yr_wk")],
                               I_over_POSITIVE_unadj=sqrt(test_tallies_wk$state_population/test_tallies_wk$TOTAL)), all.x=TRUE)#using same bias correcting method as above get ratio of infections to positives

# State-week totals of total and positive test counts, and population, added in to enable alternate weighting (2021-03-18)
# sgtf_weights merged in to enable alternate weighting (2021-03-19)
svy.dat = merge(svy.dat, test_tallies_wk[, c("STUSAB", "yr_wk", "POSITIVE", "TOTAL", "INFECTIONS", "sgtf_weights")], all.x=TRUE)
svy.dat$sgtf_weights[is.na(svy.dat$sgtf_weights) | !svy.dat$SGTF_UPSAMPLING] = 1
svy.dat = merge(svy.dat, incidence_by_region[, c("HHS", "yr_wk", "HHS_INCIDENCE")], all.x=TRUE)#for states missing testing data (OH), using the region level incidence

#
save(svy.dat, file=paste0(script.basename, "/data/", "svydat_", Sys.Date(), ".RData"))



#--------------- Limitations, etc.--------------- 

#Merging multiple disparate data streams into a single unified survey design is achieved 
# here by assuming each stream is a cluster in a PPS design. The analysis is limited by 
# the limitations of this assumption, and by incomplete information about the sampling 
# design within each stream. Owing to incomplete and imprecise test data, sources with 
# substantial contribution to the sequence database may be underweighted.  

# * Weights are derived using count of test results (including negatives) by state and 
#   collection date, currently from HHS Protect. There are a few states missing, and it 
#   has proved to be difficult to assign tests unambiguoisly to their source.
# * The protocols used by different sources to select for sequences are unknown; an 
#   assumption of random selection is made here (with modifications for oversampling of 
#   SGTF samples, where noted), leading to the possibility of biased estimates.
# * Current [NS3 guidelines for sampling](https://www.aphl.org/programs/preparedness/Crisis-Management/Documents/FEB%202021%20Revised%20NS3%20Submission%20Guidance_02052021%20FNL.pdf) are inadequate:

#  > Ideally, specimens should represent geographic, demographic (e.g., age), and clinical 
#   (e.g., disease severity or outcome) diversity from across the jurisdiction. This can be
#   achieved through random selection of specimens collected within the last 7 days.


# Changelog

# * 2021-03-07: updated lineage with pangolin as master-list
# * 2021-03-07: updated test stats by collection date (some states missing!)
# * 2021-03-07: updated designation of year-week to the date of first Sunday of week
# * 2021-03-08: calibrated infections-over-positive weights to test-positivity by cluster
# * 2021-03-08: switched to a PPS design to more accurately reflect subsampling within clusters
# * 2021-03-09: added HHS region 
# * 2021-04-XX: added populations for territories; removed lab based weights; infection rates 
#   imputed using HHS Regional rates for states with missing testing data; removed 
#   submission/receive dates from svy.dat; sgtf_weights are now weights even for SGTF_UPSAMPLING = FALSE

# Known issues/next steps

# * Smooth out positve-over-sequenced weights, if needed
# * Find test counts by collection date for states consistently missing at HHS Protect 
# * Adjust multinomial logistic growth rate models for R_t
# * Locate all SGTF sequences using in silico PCR [row 5 (TCAACTCAGGACTTGTTCTTACCT) https://www.protocols.io/view/multiplexed-rt-qpcr-to-screen-for-sars-cov-2-b-1-1-br9vm966/materials]? (Probably not needed any longer)
