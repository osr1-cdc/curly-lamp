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
data_date <- Sys.Date()
# data_date <- as.Date('2022-05-26')
# This needs to be a date on which data were frozen in the CDP database, which is often Thursdays.

## List of variants to track (not just VOC or VOI, but we name them voc in these scripts):
# These variables (custom_lineage_names, voc*) are *only* used in the weekly_variant_report_nowcast.R script. They are not used in the variant_surveillance_system.R script.

# The current branch has the following custom defined lineages:
# - BA.1+ = BA.1 with R346K
# FORMERLY "custom" lineages included:
# - AY.4.2+ - AY.4.2 with Y145H and A222V
# - AY.35+ - AY.35 with E484Q
#
# All other lineages (including AY.4.2 and AY.35) are from default pangolin calls.

# Set custom lineages
custom_lineage_names <- c("BA.1+")
# NOTE! If you change the custom lineages, you much also change the "custom"
#       pangolin sql query (lines 305-320) in variant_surveillance_system.R to match!

# Set voc's for Run1
voc1 = c(# "AY.1", "AY.2",
         # "BA.1",
         "BA.1.1",
         "BA.2",
         "BA.2.12.1",
         'BA.4',
         'BA.5',
         "B.1.617.2", # Delta
         "B.1.1.529") # Omicron
# define an alternate set of vocs
# (the reason for including two sets instead of just redefining the first set is
#  to make it easier to run both sets simultaneously just changing an option
#  passed to weekly_variant_report_nowcast.R)
voc1_reduced = c(
  'B.1.1.529', # Omicron
  'B.1.617.2'  # Delta
)

# Set voc's for Run2
# THESE ARE NOW DOWNLOADED IN "variant_surveillance_system.R".
# ONLY SET VALUES HERE TO OVERRIDE THE SQL QUERY IN "variant_surveillance_system.R".
# Leave set to "NA" to use the variants automatically identified (pulled to
# "voc2_df" in variant_surveillance_system.R)
voc2_manual = c(NA)

# optionally specify additional variants that will be added on to the lineages
# from the SQL query in "variant_surveillance_system.R"
# (this will not have any effect if "voc2_manual" is used)
voc2_additional = c(#"AY.1", "AY.2",
                    "BA.1",
                    "BA.1.1",
                    "BA.2",
                    # "BA.2.12",
                    "BA.2.12.1",
                    # "BA.3",
                    'BA.4',
                    'BA.5',
                    "B.1.617.2", # Delta 
                    "B.1.1.529" # Omicron
                    )
# voc2_custom = c(voc2,
#                 custom_lineage_names)

# define an alternate set of vocs
voc2_reduced = voc1_reduced

# Set voc's for Run3
voc3 = c("B.1.1.7",   # Alpha  # and Q.1 to 8*
         "B.1.351",   # Beta   # and B.1.351.*
         "P.1",       # Gamma  # and P.*
         "B.1.617.2", # Delta # and AY.3-AY.25*
         "AY.1",
         "AY.2",
         "B.1.427",   # Epsilon  # B.1.427 & B.1.429 are aggregated on line 137 of "weekly_variant_report_nowcast.R"
         "B.1.525",   # Eta
         "B.1.526",   # Iota
         "B.1.617.1", # Kappa
         "B.1.617.3", # (unnamed)
         "P.2",       # Zeta
         "B.1.621",   # Mu
         "B.1.1.529", # Omicron # and BA.*
         "BA.2",
         "BA.2.12.1",
         "BA.4",
         "BA.5")
# define an alternate set of vocs
voc3_reduced = voc1_reduced



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

# arguments to indicate whether sublineages should be aggregated to parent lineage
P.1_agg     = TRUE
B.1.351_agg = TRUE
AY_agg      = TRUE
Q.1_3_agg   = TRUE
B.1.621_agg = TRUE
B429_7_agg  = TRUE
B.1.1.529_agg = TRUE  # aggregate omicrons


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
model_weeks = 21 # early on the model ended up including 1 more week than was set here. Now it includes the number set here.

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

# define a start time to filter out old data
# (this is intended to speed up processing of the overall dataset)
# number of weeks to produce "weighted"/"thencast" estimates
weighted_weeks <- 12
# start-time for the weighted estimates
# (this speeds up calculations by only calculating weighted variant proportions for the most recent "weighted_weeks")
# time_start_weights <- time_end - 6 - weighted_weeks*7
time_start_weights <- as.Date('2021-05-02') # keep using week of (2021-05-02 to 2021-05-08) for consistency
# start time will be the earlier of: 1) time_start_weights; 2) "model_weeks" before time_end
time_start <- min(time_start_weights, time_end - model_weeks*7 + 1) # +1 to start on Monday
