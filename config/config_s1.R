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

# custom_lineages = FALSE
# set date for data creation
# (generally set to current date to allow more portability)
# data_date <- Sys.Date()
data_date <- as.Date('2023-06-06')
# This needs to be a date on which data were frozen in the CDP database, which is often Thursdays.

# results folder name inherits from data_date for auto completion, however the set name needs to be edited to 
# the specified run set before each set is run
results_tag <- "s1"
results_folder <- paste0("results_", data_date, '_', results_tag)

# do not need to change these on a regular basis -------------------------------
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


# set end date for national and regional survey estimates
# this is generally the end of the previous week.
time_end <- data_date - as.numeric(format(data_date, '%w')) - 1
# time_end <- as.Date('2021-12-11') # VERY little data for this past week. Not worth including.
# otherwise, set manually:
# time_end <- as.Date("2021-10-30")

# set end dates for state-level estimates (Run3)
# (end of week = Saturday)
# this is generally the end of 5 of the 6 most recent weeks (doesn't include most recent week):
state_time_end = (data_date - as.numeric(format(data_date, '%w')) - 1) - (7*5:1)
# otherwise set manually:
# state_time_end=c(as.Date("2021-09-25"),as.Date("2021-10-02"),as.Date("2021-10-09"),as.Date("2021-10-16"),as.Date("2021-10-23"))


# Use data from the frozen data created on data_date
if(data_date == Sys.Date()){
  date_frozen <- "to_date(date_frozen)" # "date_frozen" is a column in pangolin table
  # flag for whether or not current data is being used
  current_data = TRUE
} else {
  date_frozen <- paste0('"', data_date, '"')
  current_data = FALSE
}

# Options:
# - default = newest data available: "to_date(date_frozen)"
# - alternative = set a date:        '"2021-11-04"'

#variable for whether or not to include the state tagged sequences
state_source <- "state_tag_included" #argument indicating whether to include state tagged data
# options:
# state_tag_included: include samples tagged by states as surveillance quality
# [anything else]: excludes samples tagged by states as surveillance quality (i.e. only NS3 & CDC sampling)

# arguments for s1 proportion
include_other = TRUE

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
# start of the first week
week0day1 = get0("week0day1",
                 ifnotfound = as.Date("2020-01-12"))
# Multinomial model includes current week + model_weeks weeks of previous data
model_weeks = 8 # early on the model ended up including 1 more week than was set here. Now it includes the number set here.
# FOR S1 species proprotion ONLY ---- modified 2022-11-15
# model weeks only include the 8 weeks included in geni analysis
# modified 2023-01-25 to chang max to -2nd week, corresponding to change Nick made using geni report from -2nd week
model_weeks = 8
model_week_max = as.numeric(as.Date(time_end-14) - week0day1) %/% 7

# Criterion for inclusion in model (i.e to be included in model, weighted share
# must be at least 0.01 in the n_recent_weeks)
# share_cutoff = 0.01



# current week
current_week = as.numeric(as.Date(data_date) - week0day1) %/% 7

# Specify the variants to include in plots. Either "top7" or "voc"
display_option = c("top7", "voc")[2]

# number of weeks (up to current_week) to include in plots
display_lookback = 8

# define a start time to filter out old data
# (this is intended to speed up processing of the overall dataset)
# number of weeks to produce "weighted"/"thencast" estimates
weighted_weeks <- 12
# start-time for the weighted estimates
# (this speeds up calculations by only calculating weighted variant proportions for the most recent "weighted_weeks")
# time_start_weights <- time_end - 6 - weighted_weeks*7
time_start_weights <- as.Date('2021-05-09') # keep using week of (2021-05-02 to 2021-05-08) for consistency
# start time will be the earlier of: 1) time_start_weights; 2) "model_weeks" before time_end
time_start <- min(time_start_weights, time_end - model_weeks*7 + 1) # +1 to start on Monday

# Option to calculate the number of confirmed cases attributable to each variant
# This is done by simplying multiplying the proportions estimated in this script
# by the number of confirmed positive cases from ICATT testing data.
calc_confirmed_infections <- TRUE

# Option to just fit the nowcast model and avoid the slower parts of the script
# (this is only valid if the run number == 2)
# this can be removed eventually. It's here to make it easier to make updates to the Nowcast model.
nowcast_only = FALSE

# This is an option that probably won't be used, but I don't want to delete it yet
# so it's hiding here just in case I want to use it again.
# grouped weights = 3 most recent weeks (up to "time_end") are grouped together for weighting; 2 weeks prior are grouped; all weeks before that are single weeks. The purpose was to avoid extreme weights, but we went with weight trimming over this option.
use_group_weights <- FALSE

# rescale the weights that are used in the multinomial Nowcast model.
# this can help avoid numerical overflow when trying to calculate prediction intervals.
rescale_model_weights <- TRUE
# how to rescale model weights
rescale_model_weights_by <- 'mean'
# options: "max", "mean", [number]

# optionally remove UTAH PHL sequences (b/c they were causing issues with Region 8 estimates in January, 2022)
remove_utahphl <- FALSE

# optionally remove BROAD sequences (b/c they were having trouble with dropout on the Omicron spike protein, resulting in an inability to distinguish between BA.1 and BA.1+R346K in Jan/Feb 2022)
remove_broad <- FALSE

# optionally remove Quest sequences (b/c there seems to be some XBB sequences from Quest that were not received, so need to make sure the overall proportion is not skewed. Dec. 12, 2022)
remove_Quest <- TRUE
remove_Quest_cutoff <- "2022-10-08"
remove_Quest_cutoff_end <- data_date
received_Quest_cutoff <- "2023-01-17"

# new option for testing data exclustion
exclude_testing_data_portion <- TRUE
exclusion_states <- c("TX")
testing_exclusion_cutoff <- "2022-06-26"
testing_exclusion_cutoff_end <- "2022-07-02"
