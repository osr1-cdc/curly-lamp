# Overview -------------------------------------------------------------
# Setting for each run:
# the command-line argument "run_number" contributes to the R variable "tag",
# with options:
# Run1, Run2, Run3
# - Run1 calculates the variant share/proportion and confidence interval using
#   survey design ("myciprop" function) for fortnights (HHS regions & nationally)
#   and weeks  (HHS regions & nationally). It produces 2 files as output:
#   - "results/variant_share_weighted_",       ci.type,"CI_",svy.type,"_",data_date,tag,".csv"
#   - "results/variant_share_weekly_weighted_",ci.type,"CI_",svy.type,"_",data_date,tag,".csv"
# - Run2 does everything that "Run1" does and also fits the "nowcast" model,
#   which estimated smoothed trends in variant share: National and by HHS Region
#   NOTE: Run2 also produces output for Run1! (to get complete Run1 output, you
#         must run Run1 AND Run2.)
#   It creates 38 output files:
#   - "/results/wow_growth_variant_share",data_date,tag,".csv"
#   - "/results/updated_nowcast_fortnightly_",data_date,"_state_tag_included_Run1.csv"
#   - "/results/updated_nowcast_fortnightly_",data_date,"_state_tag_included_Run2.csv"
#   - "/results/updated_nowcast_weekly_",data_date,"_state_tag_included_Run1.csv"
#   - "/results/updated_nowcast_weekly_",data_date,"_state_tag_included_Run2.csv"
#   - and 33 images
#      - "/results/wtd_shares_",data_date,"_","barplot_US",tag,".jpg"
#      - "/results/wtd_shares_",data_date,"_","projection_US",tag,".jpg"
#      - "/results/wtd_shares_",data_date,"_","growthrate_US",tag,".png"
#      Create the same 3 figures for each HHS region
#      - "/results/wtd_shares_",data_date,"_","barplot_HHS",hhs,tag,".jpg"
#      - "/results/wtd_shares_",data_date,"_","projection_HHS",hhs,tag,".jpg"
#      - "/results/wtd_shares_",data_date,"_","growthrate_HHS", hhs,tag,".jpg"
# - Run3 createws state-level estimates based on rolling 4 wk bins. It saves a
#   single file:
#   - "/results/state_weighted_roll4wk_",ci.type,"CI_svyNEW_",data_date,tag,".csv"


# Update the following each run ----------------------------------------
custom_lineages = TRUE

# specify survey design type (NO NEED TO CHANGE)
svy.type <- "svyNEW"
# Survey design options:
# svyREG: stratified by HHS;
#         Primary sampling unit = state
#         Secondary sampling unit = source
# [anything else]: stratified by state & week;
#                  Sampling unit = source

# specify confidence interval type for weighted proportion estimates
# (NO NEED TO CHANGE)
ci.type <- "KG"
# ci.type options include:
#   xlogit: uses function "svyciprop" with method='xlogit'
#   [anything else]: uses function "svycipropkg"

# Update the following each run ----------------------------------------
# set end date for national and regional survey estimates
time_end <- as.Date("2021-11-13")
# this is generally the end of the previous week:   Sys.Date() - as.numeric(format(Sys.Date(), '%w')) - 1
# alternatively, just set this manually:            as.Date("2021-11-13")

# set end dates for state-level estimates (Run3)
# (end of week = Saturday)
state_time_end = c( # as.Date("2021-09-18"),
  # as.Date("2021-09-25"),
  # as.Date("2021-10-02"),
  as.Date("2021-10-09"),
  as.Date("2021-10-16"),
  as.Date("2021-10-23"),
  as.Date("2021-10-30"),
  as.Date("2021-11-06")
)
# this is generally the end of 5 of the 6 most recent weeks (doesn't include most recent week):    (Sys.Date() - as.numeric(format(Sys.Date(), '%w')) - 1) - (7*5:1)

# set date for data creation
# (generally set to current date to allow more portability)
# data_date <- Sys.Date()
data_date <- as.Date('2021-11-18')

# Use data from the frozen data created on this date
date_frozen <- if(data_date == Sys.Date()){
  "to_date(date_frozen)"
} else {
  paste0('"', data_date, '"')
}
# Options:
# - default = newest data available: "to_date(date_frozen)"
# - alternative = set a date:        '"2021-11-04"'

#variable for whether or not to include the state tagged sequences
state_source <- "state_tag_included" #argument indicating whether to include state tagged data
# options:
# state_tag_included: include samples tagged by states as surveillance quality
# [anything else]: excludes samples tagged by states as surveillance quality (i.e. only NS3 & CDC sampling)

# a tag for the filename to indicate which run from Lab TF request the results are for
tag <- paste0("_",state_source,"_Run", opts$run_number)
# potentially adjust if there are custom lineages to include
custom_tag = ifelse(custom_lineages == TRUE, "_custom", "")
# options:
# paste0("_",state_source,"_Run1"): calc proportions USING SURVEY DESIGN for reduced set of VOCs
# paste0("_",state_source,"_Run2"): calc proportions for extended set of VOCs (many Delta subclades)
#                                   this is the only run that includes the multinomial "nowcast" model
# Note: Nowcast runs best with a large set of variants (if you group everything into delta, then the model breaks)
# paste0("_",state_source,"_Run3"): calc proportions for another reduced set of VOCs
#                                   Specific to state-level runs & VOC's generally won't change

# arguments to indicate whether sublineages should be aggregated to parent lineage
P.1_agg     = TRUE
B.1.351_agg = TRUE
AY_agg      = TRUE
Q.1_3_agg   = TRUE
B.1.621_agg = TRUE
B429_7_agg  = TRUE

# List of variants to track (not just VOC or VOI):
# VOCs
# The current branch has the following custom defined lineages:
# - AY.4.2+ - AY.4.2 with Y145H and A222V
# - AY.35+ - AY.35 with E484Q
#
# All other lineages (including AY.4.2 and AY.35) are from default pangolin calls.

# Set custom lineages
custom_lineage_names = c("AY.35+",
                         "AY.4.2+")
# Set voc's for Run1
voc1 = c("AY.1",
         "AY.2",
         "B.1.617.2")
voc1_custom = c(voc1,
                custom_lineage_names)
# Set voc's for Run2
voc2 = c(
  "AY.100",
  "AY.103",
  "AY.117",
  "AY.118",
  "AY.119",
  "AY.14",
  "AY.20",
  "AY.25",
  "AY.26",
  "AY.3",
  "AY.3.1",
  "AY.39",
  "AY.44",
  "AY.47",
  "AY.75",
  "B.1.617.2"
)
voc2_custom = c(voc2,
                custom_lineage_names)
# Set voc's for Run3
voc3 = c("B.1.1.7",# with Q.1 to 8*
         "B.1.351", #and B.1.351.*
         "P.1", #and P.*
         "B.1.617.2", #and AY.3-AY.25*
         "AY.1",
         "AY.2",
         "B.1.427",#/B.1.429* # B.1.427 & B.1.429 are aggregated on line 137 of "weekly_variant_report_nowcast.R"
         "B.1.525",
         "B.1.526",
         "B.1.617.1",
         "B.1.617.3",
         "P.2",
         "B.1.621")# and B.1.621.1*
voc3_custom = c(voc3,
                custom_lineage_names)

# Choose which list of vocs to use based on the run number
if( grepl("Run1",tag) ){
  if(custom_lineages == FALSE) {
    voc = voc1
  } else {
    voc = voc1_custom
  }
}
if( grepl("Run2",tag) ) {
  if(custom_lineages == FALSE) {
    voc = voc2
  } else {
    voc = voc2_custom
  }
}
if( grepl("Run3",tag) ) {
  if(custom_lineages == FALSE) {
      voc = voc3
  } else {
      voc = voc3_custom
  }
}


# Argument determining whether figures should be output as jpgs
fig_gen_run = TRUE

# define "not in" function
`%notin%` <- Negate(`%in%`)

# Some parameters defining what is modeled and displayed ---------------
# Top n variants by variant share that must be included in output
n_top = 10
# Only a subset of variants that are/were common *in recent weeks* are included
# in the multinomial model. (Note that the multinomial model itself is run on
# more weeks of data than are used for defining common variants.) Focal
# variants are determined using sample weights and "n_recent_weeks" of data.

# Window for estimates (focus to top variants in this model)
n_recent_weeks = 7

# Multinomial model includes current week + model_weeks weeks of previous data
model_weeks = 20

# Criterion for inclusion in model (i.e to be included in model, weighted share
# must be at least 0.01 in the n_recent_weeks)
# share_cutoff = 0.01

# start of the first week
week0day1 = get0("week0day1",
                 ifnotfound = as.Date("2020-01-05"))

# current week
current_week = as.numeric(as.Date(data_date) - week0day1) %/% 7

# Specify the variants to include in plots. Either "top7" or "voc"
display_option = c("top7", "voc")[1]

# number of weeks (up to current_week) to include in plots
display_lookback = 8
