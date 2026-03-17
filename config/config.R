# Overview for configuration settings
#
# Update the following each run ----------------------------------------
# set date for data creation
data_date <- as.Date('2026-01-08')
# This needs to be a date on which data were frozen in the CDP database
date_frozen_toread <- data_date
# If the data was already pulled and you want to just use that data instead of re-pulling it, set here.
# This is useful if you aggregate some lab names at the end of this code and then want to re-run the
# script after changing which labs get aggregated.
use_previously_imported_data <- FALSE

# results folder name inherits from data_date for auto completion, however the set name needs to be edited to
# the specified run set before each set is run
results_tag <- "CDT"
results_folder <- paste0("results_", data_date, '_', results_tag)

# If pre_aggregation is TRUE, force aggregate sublineages to voc1 list, no need to generate run1 postaggregated nowcast results in run2.
pre_aggregation <- FALSE
## List of variants to track (not just VOC or VOI, but we name them voc in these scripts):

voc1 = c(
         "BA.1.1",
         "BA.2",
         "BA.2.12.1",
         "BA.2.75",
         'BA.2.86',
         "BA.2.75.2",
         "CH.1.1",
         "BF.7",
         "BF.11",
         'BA.4',
         'BA.4.6',
         'BA.5',
         'BA.5.2.6',
         "BQ.1",
         "BQ.1.1",
         "BN.1",
         "XBB",
         'XBB.1.5',
         'XBB.1.5.1',
         'XBB.1.9.1',
         'XBB.1.5.10',
         'FD.1.1',
         'FD.2',
         'EU.1.1',
         'XBB.1.5.59',
         'XBB.1.5.68',
         'XBB.1.5.70',
         'XBB.1.5.72',
         'JD.1.1',
         'GK.1.1',
         'GK.2',
         'HV.1',
         'XBB.1.9.2',
         'JG.3',
         'HK.3',
         'EG.5',
         'EG.5.1.8',
         'EG.6.1',
         'XBB.1.16',
         'XBB.1.16.1',
         'XBB.1.16.6',
         'JF.1',
         'XBB.1.16.11',
         'XBB.1.16.15',
         'XBB.1.16.17',
         'HF.1',
         'FE.1.1',
         'FL.1.5.1',
         'XBB.2.3',
         'XBB.2.3.8',
         'GE.1',
         'XBB.1.42.2',
         'JN.1',
         'JN.1.7',
         'JN.1.8.1',
         'JN.1.18',
         'JN.1.13',
         'JN.1.13.1',
         'JN.1.11.1',
         'JN.1.16',
         'KP.3',
         'JN.1.16.1',
         'JN.1.32',
         'KS.1',
         'KW.1.1',
         'KV.2',
         'KQ.1',
         'KP.2',
         'KP.2.3',
         'KP.3.1.1',
         'LF.3.1',
         'KP.1.1',
         'KP.1.2',
         'XDP',
         'LB.1',
         'XDV.1',
         'B.1.617.2', # Delta
         'B.1.1.529',
         'JN.1.4.3',
         'KP.4.1',
         'KP.1.1.3',
         'KP.2.15',
         'XEC',
         'LP.1',
	       'MC.1',
	       'LB.1.3.1',
         'LF.7',
         'MC.10.1',
         'LP.8.1',
         'MC.19',
         'XEK',
         'XEQ',
         'MC.28.1',
         'XEC.4',
         'XFC',
         'LF.7.2.1',
         'PA.1',
         'LF.7.7.2',
         'LF.7.7.1',
         'JN.1.18.6',
         'XFG',
         'NB.1.8.1',
         'LF.7.9',
         'NW.1',
         'XFG.1',
         'XFG.14.1',
         'XFZ',
         'XFV',
         'XFG.6',
         'XFY') # Omicron
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
voc2_additional = c(
                    'BA.2.75.2',
                    'BQ.1',
                    'BQ.1.1',
                    'BA.1.1',
                    'BA.2',
                    'BA.2.12.1',
                    'BA.2.75',
                    'BA.2.86',
                    'CH.1.1',
                    'BF.7',
                    'BF.11',
                    'BA.4',
                    'BA.4.6',
                    'BA.5',
                    'BA.5.2.6',
                    'BN.1',
                    'XBB',
                    'XBB.1.5',
                    'XBB.1.16.1',
                    'XBB.1.5.1',
                    'XBB.1.9.1',
                    'XBB.1.5.10',
                    'FD.1.1',
                    'FD.2',
                    'EU.1.1',
                    'XBB.1.5.59',
                    'XBB.1.5.68',
                    'XBB.1.5.70',
                    'XBB.1.5.72',
                    'JD.1.1',
                    'GK.1.1',
                    'GK.2',
                    'HV.1',
                    'XBB.1.9.2',
                    'JG.3',
                    'HK.3',
                    'EG.5',
                    'EG.5.1.8',
                    'EG.6.1',
                    'XBB.1.16',
                    'XBB.1.16.6',
                    'JF.1',
                    'XBB.1.16.11',
                    'XBB.1.16.15',
                    'XBB.1.16.17',
                    'HF.1',
                    'FE.1.1',
                    'FL.1.5.1', 
                    'XBB.2.3',
                    'XBB.2.3.8',
                    'GE.1',
                    'XBB.1.42.2',
                    'B.1.617.2', # Delta
                    'B.1.1.529',
                    'KQ.1',
                    'KP.1.1',
                    'KP.1.2',
                    'KP.2',
                    'KP.2.3',
                    'KP.3.1.1',
                    'LF.3.1',
                    'KS.1',
                    'JN.1',
                    'JN.1.7',
                    'JN.1.8.1',
                    'JN.1.16',
                    'JN.1.18',
                    'KP.3',
                    'JN.1.16.1',
                    'JN.1.32',
                    'KW.1.1',
                    'KV.2',
                    'XDP',
                    'LB.1',
                    'XDV.1',
                    'JN.1.13.1',
                    'JN.1.11.1',
                    'JN.1.13',
                    'JN.1.4.3',
                    'KP.4.1',
                    'LP.1',
                    'KP.2.15',
                    'XEC',
                    'KP.1.1.3',
                    'LF.7',
                    'MC.10.1',
                    'LP.8.1',
                    'MC.19',
	                    'XEK',
	                    'XEQ',
                    'MC.28.1',
                    'XEC.4',
                    'XFC',
                    'LF.7.2.1',
                    'PA.1',
                    'LF.7.7.2',
                    'LF.7.7.1',
                    'JN.1.18.6',
                    'XFG',
                    'NB.1.8.1',
                    'LF.7.9',
                    'NW.1',
                    'XFG.1',
                    'XFG.14.1',
                    'XFZ',
                    'XFV',
                    'XFG.6',
                    'XFY'
                    )

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

date_frozen <- paste0('"', data_date, '"')
if(data_date == Sys.Date()){
  # flag for whether or not current data is being used
  current_data = TRUE
} else {
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
XBB_agg = TRUE # aggregate XBBs
# preaggregation for CDT run
force_preaggregate_XBB = FALSE
force_aggregate_XBB_except <- c("XBB.1", "XBB.1.5")
force_preaggregate_BN.1 = FALSE

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
model_weeks = 21
model_time_end = time_end
# Criterion for inclusion in model (weighted share must be at least 0.01 in the n_recent_weeks)

# start of the first week
week0day1 = get0("week0day1",
                 ifnotfound = as.Date("2020-01-12"))
                 # for 2023-5-9 to be implemented biweek estimates, change to 2020-01-12

# current week
current_week = as.numeric(as.Date(data_date) - week0day1) %/% 7

# Specify the variants to include in plots. Either "top7" or "voc"
display_option = c("top7", "voc")[2]

# number of weeks (up to current_week) to include in plots
display_lookback = 8

# define a start time to filter out old data
# (this is intended to speed up processing of the overall dataset)
# start-time for the weighted estimates
time_start_weights <- as.Date('2021-05-09')
# for 2023-5-9 to be implemented biweek estimates, change to 2021-05-09
# start time will be the earlier of: 1) time_start_weights; 2) "model_weeks" before time_end
time_start <- min(time_start_weights, time_end - model_weeks*7 + 1) # +1 to start on Monday

# Option to calculate the number of confirmed cases attributable to each variant
# This is done by simplying multiplying the proportions estimated in this script
# by the number of confirmed positive cases from ICATT testing data.
calc_confirmed_infections <- FALSE

# Option to just fit the nowcast model and avoid the slower parts of the script
# (this is only valid if the run number == 2)
# this can be removed eventually. It's here to make it easier to make updates to the Nowcast model.
nowcast_only = FALSE

# This is an option that probably won't be used, but I don't want to delete it yet
# so it's hiding here just in case I want to use it again.
# grouped weights = 3 most recent weeks (up to "time_end") are grouped together for weighting; 2 weeks prior are grouped; all weeks before that are single weeks. The purpose was to avoid extreme weights, but we went with weight trimming over this option.
use_group_weights <- FALSE

# force_aggregate_R346T will force all custom R346T lineages aggregate to a total "R346T" lineage, and include R346T in voc
# This is to accompany the custom pull to have main lineages with R346T separately, and would like to have a combined R346T lineage analysis
# With this option turned on, all Omicron with R346T will be grouped together before modeling analysis.
force_aggregate_R346T <- FALSE

# force_aggregate_xxx will REMOVE variants from the voc list, which will result
# in their subsequent aggregation into a parent lineage *IF* the parent lineage
# is included in VOC.
# This option will likely go away with time, but for now it's here to override
# the automated voc2 selection in some cases. This will prevent splitting omicron
# into more sublineages than we want.
# if 'B.1.1.529' is listed in the vocs, then force-aggregate omicron sublineages
# even if sublineages are also included in the vocs. To avoid aggregating a
# specific sublineage, include the sublineage in both "voc" and
# "force_aggregate_omicron_except".
# THIS WILL LIKELY NEED TO BE REPLACED IN THE FUTURE, BUT IT'S HERE TO AVOID
# SPLITTING OUT BA.1, WHICH IS OFTEN AUTOMATICALLY INCLUDED IN VOC2 B/C IT'S > 1% NATIONALLY.
force_aggregate_omicron <- FALSE
# list omicron sublineages that will not be aggregated (if they are also in voc) (these are the only Omicron sublineages that will be permitted)
force_aggregate_omicron_except <- c(
         "BA.1.1",
         "BA.2",
         "BA.2.12.1",
         "BA.2.75",
         "BA.2.75.2",
         "BF.7",
         'BA.4',
         'BA.4.6',
         'BA.5',
         'BA.5.2.6',
         "BQ.1",
         "BQ.1.1")


# force-aggregate Delta sublineages
# this option will force any/all Delta sublineages that show up in the vocs to be
# aggregated into "B.1.617.2". This was introduced on 2022-02-01 to help stabalize
# Nowcast estimates of the BA.2 sublineage.
# Another option would be to use the "voc_manual" setting in config/config.R to prevent
# AY sublineages from being split out from Delta. The difference between using "voc_manual"
# and using "force_aggregate_delta" is that using "voc_manual" requires one to look up
# all other lineages > 1% for inclusion.
# another option to try to control this problem is to set "n_top" to a low number.
force_aggregate_delta <- FALSE

# force-aggregate "B" into "other"
# variant "B" most likely indicates trouble sequencing, rather than an actual variant, so don't split it out.
# this also prevents "B" from being included in the Nowcast model.
force_aggregate_B <- TRUE
# same for "B.1"
force_aggregate_B.1 <- TRUE

# rescale the weights that are used in the multinomial Nowcast model.
# this can help avoid numerical overflow when trying to calculate prediction intervals.
rescale_model_weights <- TRUE
# how to rescale model weights.
# code will try all options (in order) until a value works for both the national and regional nowcast model ("works" means that the multinomial model fit successfully estimates the Hessian, allowing for estimates of uncertainty)
rescale_model_weights_by <- c('mean', "max", 'sum', 100, 'min', 2^seq(1:14))
# options: "max", "mean", [number]
# note that R vectors can only have 1 data type, so combining strings with numbers results in all values being treated as string values. That's ok. The code converts number strings to numeric.

# optionally remove UTAH PHL sequences (b/c they were causing issues with Region 8 estimates in January, 2022)
remove_utahphl <- FALSE

# optionally remove BROAD sequences (b/c they were having trouble with dropout on the Omicron spike protein, resulting in an inability to distinguish between BA.1 and BA.1+R346K in Jan/Feb 2022)
# Added 2023-02-07, removing broad sequences submited on specific dates since Danny Park from Broad notified us these sequences need to be retracted
remove_broad <- TRUE
received_broad_dates <- c('2023-01-12', '2023-02-02')

# optionally remove Quest sequences (b/c there seems to be some XBB sequences from Quest that were not received, so need to make sure the overall proportion is not skewed. Dec. 12, 2022)
remove_Quest <- TRUE
remove_Quest_cutoff <- "2022-10-08"
remove_Quest_cutoff_end <- data_date
# added 2023-01-31, quest data received after 2023-01-17 is supposed to be error-free
# updated 2023-02-07, changed the cutoff date to 2023-01-23
received_Quest_cutoff <- "2023-01-23"

# new option for testing data exclustion
exclude_testing_data_portion <- TRUE
exclusion_states <- c("TX")
testing_exclusion_cutoff <- "2022-06-26"
testing_exclusion_cutoff_end <- "2022-07-02"

# Location of JDBC Driver 

jdbc_driver <- '/scicomp/groups-pure/OID/NCIRD/CORVD/CRVLB/MIST/apps/jdbc/ClouderaImpalaJDBC-2.6.20.1024/ClouderaImpalaJDBC41-2.6.20.1024/'
