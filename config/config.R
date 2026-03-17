# Run-specific settings ---------------------------------------------------------
data_date <- as.Date('2026-01-08')
date_frozen_toread <- data_date
use_previously_imported_data <- FALSE

results_tag <- "CDT"
results_folder <- paste0("results_", data_date, '_', results_tag)

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
         
voc1_reduced = c(
  'B.1.1.529', # Omicron
  'B.1.617.2'  # Delta
)

# Run 2 variants. Leave `voc2_manual` as `NA` to use the downloaded list.
voc2_manual = c(NA)

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

# Stable settings ---------------------------------------------------------------

svy.type <- "svyNEW"

ci.type <- "KG"

# set end date for national and regional survey estimates
# this is generally the end of the previous week.
time_end <- data_date - as.numeric(format(data_date, '%w')) - 1

date_frozen <- paste0('"', data_date, '"')
if(data_date == Sys.Date()){
  current_data = TRUE
} else {
  current_data = FALSE
}

state_source <- "state_tag_included" #argument indicating whether to include state tagged data

# Parent-lineage aggregation switches ------------------------------------------
P.1_agg     = TRUE
B.1.351_agg = TRUE
AY_agg      = TRUE
Q.1_3_agg   = TRUE
B.1.621_agg = TRUE
B429_7_agg  = TRUE
B.1.1.529_agg = TRUE  # aggregate omicrons
XBB_agg = TRUE # aggregate XBBs
force_preaggregate_XBB = FALSE
force_aggregate_XBB_except <- c("XBB.1", "XBB.1.5")
force_preaggregate_BN.1 = FALSE

# Argument determining whether figures should be output as jpgs
fig_gen_run = TRUE

# define "not in" function
`%notin%` <- Negate(`%in%`)

# Model and display parameters --------------------------------------------------
n_top = 10
n_recent_weeks = 7
model_weeks = 21
model_time_end = time_end
# Criterion for inclusion in model (weighted share must be at least 0.01 in the n_recent_weeks)

# start of the first week
week0day1 = get0("week0day1",
                 ifnotfound = as.Date("2020-01-12"))

current_week = as.numeric(as.Date(data_date) - week0day1) %/% 7

# Specify the variants to include in plots. Either "top7" or "voc"
display_option = c("top7", "voc")[2]

# number of weeks (up to current_week) to include in plots
display_lookback = 8

# Data window ------------------------------------------------------------------
time_start_weights <- as.Date('2021-05-09')
time_start <- min(time_start_weights, time_end - model_weeks*7 + 1) # +1 to start on Monday

calc_confirmed_infections <- FALSE

# Option to just fit the nowcast model and avoid the slower parts of the script (this is only valid if the run number == 2)
nowcast_only = FALSE

# This is an option that probably won't be used
use_group_weights <- FALSE

# force_aggregate_R346T will force all custom R346T lineages aggregate to a total "R346T" lineage, and include R346T in voc
force_aggregate_R346T <- FALSE

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
# this option will force any/all Delta sublineages that show up in the vocs to be aggregated into "B.1.617.2". 
force_aggregate_delta <- FALSE

# force-aggregate "B" into "other"
# variant "B" most likely indicates trouble sequencing, rather than an actual variant, so don't split it out.
force_aggregate_B <- TRUE
# same for "B.1"
force_aggregate_B.1 <- TRUE

# rescale the weights that are used in the multinomial Nowcast model - this can help avoid numerical overflow when trying to calculate prediction intervals.
rescale_model_weights <- TRUE

# how to rescale model weights.
rescale_model_weights_by <- c('mean', "max", 'sum', 100, 'min', 2^seq(1:14))

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
received_Quest_cutoff <- "2023-01-23"

# new option for testing data exclustion
exclude_testing_data_portion <- TRUE
exclusion_states <- c("TX")
testing_exclusion_cutoff <- "2022-06-26"
testing_exclusion_cutoff_end <- "2022-07-02"

# Location of JDBC Driver 
jdbc_driver <- '/scicomp/groups/covlab/apps/jdbc/ClouderaImpalaJDBC-2.6.20.1024/ClouderaImpalaJDBC41-2.6.20.1024/'
