# Updates ----------------------------------------------------------------------
#
#   Variant share estimates and nowcasts
#   created by Prabasaj Paul
#
# Updates:
#   ~ weight trimming to avoid overly-influential sequences when samples are small
#   ~ Error handling added to multinomial model fitting to handle non-invertible
#     Hessian
#   ~ Multinomial model run on aggregated data to speed up model fitting
#   ~ Survey design updated to have PSU=Source and strata=week and state
#   ~ Less reliable CIs are now flagged based on NCHS data presentation
#     standards for proportions
#   ~ Code now generates state-level weighted estimates for rolling 4 week
#     time bins
#   ~ Multinomial model code updated to formally account for survey design
#     in variance estimation via svyrecvar
#   ~ Function to calculate binomial CI for nowcast proportion

# Setup ------------------------------------------------------------------------
library(optparse)   # to get arguments from the command line
library(survey)     # package with survey desgin functions
library(nnet)       # package with multinomial regression for nowcast
library(data.table) # package for speeding up calculation of simple adjusted weights

# set global options:
#  - tell the survey package how to handle surveys where only a single primary
#    sampling units has observations from a particular domain or subpopulation.
#    See more here:
#    https://r-survey.r-forge.r-project.org/survey/exmample-lonely.html
#  - do not convert strings to factors (as was default in R < 4.0)
options(survey.adjust.domain.lonely = T,
        survey.lonely.psu = "average",
        stringsAsFactors = FALSE)

## optparse option list --------------------------------------------------------
{
  # (get the run number from the command line)
  option_list <- list(

    # Run number
    optparse::make_option(
      opt_str = c("-r", "--run_number"),
      type    = "character",
      default = "1",
      help    = "Run number",
      metavar = "character"),
    # options:
    # Run1: calc proportions USING SURVEY DESIGN for VOCs specified in voc1
    # Run2: calc proportions for VOCs with >= 1% unweighted share (many Delta subclades)
    #       this is the only run that includes the multinomial "nowcast" model
    #     Note: Nowcast runs best with a large set of variants (or very few variants with Omicron)
    #           (if you group everything into Delta, then the model breaks)
    # Run3: calc state-level proportions in 4-week bins for VOCs specified in voc3
    #       runs & VOC's generally won't change

    # whether or not to use Custom Lineages
    optparse::make_option(
      opt_str = c("-c", "--custom_lineages"),
      type    = "character",
      default = "F",
      help    = "Whether or not to use custom lineages (character value of T or F)",
      metavar = "character"
    ),

    # whether or not to use reduced vocs
    optparse::make_option(
      opt_str = c("-v", "--reduced_vocs"),
      type    = "character",
      default = "F",
      help    = "Whether or not to use reduced set of vocs (character value of T or F)",
      metavar = "character"
    ),

    # whether (and how) to trim weights.
    # options include:
    #   - "quantile_99" = 99th percentile. This must contain the letters "quant" followed by a number.
    #   - "F", "FALSE", "N", or "NO" = no weight trimming
    #   - "IQR" = median + 6 * IQR
    optparse::make_option(
      opt_str = c("-t", "--trim_weights"),
      type    = "character",
      default = "quantile_99",
      help    = "Maximum weight for any individual sequence",
      metavar = "character"
    ),

    # optionally save datasets that are produced by this script
    # setting to true will save:
    #  - src.dat (processed data used for weighted variant proportions)
    #  - data used for the nowcast models
    optparse::make_option(
      opt_str = c("-s", "--save_datasets_to_file"),
      type    = "character",
      default = "F",
      help    = "Maximum weight for any individual sequence",
      metavar = "character"
    )
  )

  # parsing options list
  opts <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

  # use the specified flags to set several variables
  # convert custom_lineages flag to a logical value
  if(toupper(opts$custom_lineages) %in% c('T', 'TRUE', 'Y', 'YES')){
    custom_lineages = TRUE
    custom_tag = "_custom"
  } else {
    if(toupper(opts$custom_lineages) %in% c('F', 'FALSE', 'N', 'NO')){
      custom_lineages = FALSE
      custom_tag = ""
    } else {
      errorCondition(message = paste0('custom_lineages must be "T" or "F". Argument provide: ', opts$custom_lineages))
    }
  }
  # convert reduced_vocs flag to a logical value
  if(toupper(opts$reduced_vocs) %in% c('T', 'TRUE', 'Y', 'YES')){
    reduced_vocs = TRUE
    reduced_voc_tag = "_reduced_vocs"
  } else {
    if(toupper(opts$reduced_vocs) %in% c('F', 'FALSE', 'N', 'NO')){
      reduced_vocs = FALSE
      reduced_voc_tag = ""
    } else {
      errorCondition(message = paste0('reduced_vocs must be "T" or "F". Argument provide: ', opts$reduced_vocs))
    }
  }

  # convert trim_weights flag to a logical value
  if( toupper(opts$trim_weights) %in% c('F', 'FALSE', 'N', 'NO')){
    trim_weights = FALSE
  } else {
    trim_weights = TRUE
  }

  # convert save_datasets_to_file flag to a logical value
  if( toupper(opts$save_datasets_to_file) %in% c('F', 'FALSE', 'N', 'NO')){
    save_datasets_to_file = FALSE
  } else {
    save_datasets_to_file = TRUE
  }

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

  # another option that will go away with time, but for now it's here to avoid
  # splitting omicron into multiple sublineages.
  # if 'B.1.1.529' is listed in the vocs, then force-aggregate omicron sublineages
  # even if sublineages are also included in the vocs. To avoid aggregating a
  # specific sublineage, include the sublineage in both "voc" and 
  # "force_aggregate_omicron_except". 
  # THIS WILL LIKELY NEED TO BE REPLACED IN THE FUTURE, BUT IT'S HERE TO AVOID
  # SPLITTING OUT BA.1, WHICH IS OFTEN AUTOMATICALLY INCLUDED IN VOC2 B/C IT'S > 1% NATIONALLY.
  force_aggregate_omicron <- TRUE
  # list omicron sublineages that will not be aggregated (if they are also in voc) (these are the only Omicron sublineages that will be permitted)
  force_aggregate_omicron_except <- c('BA.1', 'BA.2', 'BA.3') # , 'BA.1.1')

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
}


# get base directory
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name   <- "--file="
script.name     <- sub(pattern = file.arg.name,
                       replacement = "",
                       x = initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
# set script.basename when running from Rstudio interactively
if(length(script.basename) == 0) {
  script.basename = "."
}

# create results dir
dir.create(paste0(script.basename,"/results"), showWarnings = F)


## Data prep -------------------------------------------------------------------
#capture system time
tstart = proc.time()

# source options from config.R file
source(paste0(script.basename, "/config/config.R"))

#Source the functions used in this script.
# source(paste0(script.basename, "/svycipropkg.R"))
source(paste0(script.basename, "/weekly_variant_report_functions.R"))

# Load output from variant_surveillance_system.r
# (filtered genomic surveillance data)
load(paste0(script.basename, "/data/svydat_", data_date, custom_tag, ".RData"))


# create a tag for the filenames to differentiate results from different runs
tag <- paste0("_",state_source,"_Run", opts$run_number, reduced_voc_tag, custom_tag)

### choose vocs ----------------------------------------------------------------
# - run number
# - custom_lineage
# - reduced_vocs
if( grepl("Run1", tag) ){
  if(reduced_vocs){
    voc = voc1_reduced
  } else {
    voc = voc1
  }
}
if( grepl("Run2", tag) ) {
  # get the list of vocs from either the manually specified values in "config/config.R"
  # or from the automatically calculated values in "variant_surveillance_system.R"
  if(reduced_vocs){
    voc = voc2_reduced
    # also need the voc1 vocs for the aggregation matrix from the nowcast model
    voc1 = voc1_reduced
  } else {
    if( is.na(voc2_manual) ) {
    # read in the list of voc created in "variant_surveillance_system.R"
    # and add on any variants specified in voc2_additional
    voc = unique( # make sure there are no duplicates
      c(
        readRDS(file = paste0(script.basename,
                              "/data/voc2_auto_", data_date, custom_tag, ".RDS")),
        # add in the "additional" vocs that need to be included
        voc2_additional
      ))
    # copy the 'voc2_auto' file over to the results folder
    file.copy(
      from = paste0(script.basename,
                    "/data/voc2_auto_", data_date, custom_tag, ".RDS"),
      to = paste0(script.basename,
                  "/results/voc2_auto_", data_date, custom_tag, ".RDS")
    )

    } else {
      voc = voc2_manual
    }
  }
}
if ( grepl("Run3", tag) ) {
  if (reduced_vocs){
    voc = voc3_reduced
  } else {
    voc = voc3
  }
}

# optionally add on the custom lineages
if (custom_lineages == TRUE) {
  voc = unique(c(voc, custom_lineage_names))
}

# force-aggregate omicron
# (even if BA.1 is listed in voc2, this will force all BA sublineages of Omicron
#  to have the same VARIANT name (B.1.1.529)) (unless explicitly excluded in "force_aggregate_omicron_except")
if ( force_aggregate_omicron & ('B.1.1.529' %in% voc) ){

  # omicron sub-lineages to exclude/aggregate
  omicron_sublineages <- voc[ grepl('(BA\\.[0-9])', voc, ignore.case = TRUE) ]

  # but don't force-aggregate any sublineages in "force_aggregate_omicron_except"
  # if they are also in the vocs.
  omicron_sublineages <- omicron_sublineages[omicron_sublineages %notin% force_aggregate_omicron_except]

  # remove omicron sublineages (leaving "B.1.1.529")
  voc <- voc[ voc %notin% omicron_sublineages ]
}

if (force_aggregate_delta){
  # Delta sublineages to exclude
  delta_sublineages <- voc[ grepl('^AY', voc, ignore.case = TRUE)]

  # remove Delta sublineages
  voc <- voc[ voc %notin% delta_sublineages ]
}

# force-aggregate B into "Other"
if (force_aggregate_B & ('B' %in% voc)) {
  voc <- voc[ voc %notin% 'B' ]
}

############# REMOVE specified labs
if (remove_utahphl){
  # grep(pattern = 'utah', x = unique(svy.dat$SOURCE), ignore.case = T, value = T)
  svy.dat <- subset(x = svy.dat,
                    subset = SOURCE != 'UTAH PUBLIC HEALTH LABORATORY')
}
if (remove_broad){
  svy.dat <- subset(x = svy.dat,
                    subset = SOURCE %notin% c('BROAD INSTITUTE', 'INFECTIOUS DISEASE PROGRAM, BROAD INSTITUTE OF HARVARD AND MIT'))
}

### subset data ----------------------------------------------------------------

# Convert any factor to string
# (redundant with "stringsAsFactors = FALSE")
column_classes = sapply(svy.dat, class) # get the class of each column
factor_columns = names(column_classes[column_classes == "factor"]) # names of the columns that are factors
# convert columns that are factors to strings
for (vv in factor_columns) svy.dat[, vv] = as.character(svy.dat[, vv])

# subset the survey data based on whether or not state-tagged data should be included
if(state_source == "state_tag_included"){

  # only include samples where both the lab and the variant are defined
  src.dat = subset(x = svy.dat,
                   SOURCE != "OTHER" &
                     !is.na(VARIANT) &
                     VARIANT != "None")
} else {
  # only include samples from these labs (i.e. no state-tagged data)
  src.dat = subset(x = svy.dat,
                   SOURCE %in% c("UW VIROLOGY LAB",
                                 "FULGENT GENETICS", # NOTE: there are ~700 FULGENT sequences that were NOT part of CDC contractor sequencing...
                                 "HELIX",
                                 "HELIX/ILLUMINA",
                                 "LABORATORY CORPORATION OF AMERICA",
                                 "AEGIS SCIENCES CORPORATION",
                                 "QUEST DIAGNOSTICS INCORPORATED",
                                 "BROAD INSTITUTE",
                                 "INFINITY BIOLOGIX",
                                 "NS3",
                                 "MAKO MEDICAL"))
  # exclude samples where the variant has not been identified
  src.dat = subset(x = src.dat,
                   !is.na(VARIANT) &
                     VARIANT != "None")
  # this *might* be adequate for identifying the FULGENT sequences that were NOT
  # part of NS3 or CDC sequencing
  #    & !is.na(src.dat$covv_accession_id)
}


### SGTF weights ---------------------------------------------------------------

# set s-gene upsampling weights to 1 until I move the SGTF weights here...
src.dat$sgtf_weights = 1
# set this so that default SAW_ALT will be 1
# src.dat$HHS_INCIDENCE = 1/src.dat$state_population

# SGTF WEIGHT CALCULATIONS SHOULD BE MOVED FROM VARIANT_SURVEILLANCE_SYSTEM.R TO HERE!



# (re)calculate the survey weights
# Note: these weights are already calculated in "variant_surveillance_system.R",
#       but they're recalculated here after subsampling data based on lab.
# For info on how these weights are calculated, see notes in "variant_surveillance_system.R"
# and: https://www.cdc.gov/mmwr/volumes/70/wr/mm7023a3.htm
# (using data.table saves about 99% of the time compared to for-loop. The for-loop
#  is still here b/c it's easier to understand step-by-step calculations.)

# # (Weighted) count of sequences in each state & week
# seq.tbl = with(src.dat, xtabs((1/sgtf_weights) ~ STUSAB + yr_wk))

# for (rr in 1:nrow(src.dat)) {
#   # (estimated) number of infections represented by each positive test result
#   w_i = sqrt(src.dat$state_population[rr]/src.dat$TOTAL[rr])
#
#   # (estimated) total number of infections
#   n_infections = w_i * src.dat$POSITIVE[rr]
#
#   # number of sequences (samples) in a given week & state
#   n_sequences = seq.tbl[src.dat$STUSAB[rr], src.dat$yr_wk[rr]]
#
#   # note: n_infections / n_sequences = w_p = # of infections represented by each sequence
#
#   # survey/sample weight
#   src.dat$SIMPLE_ADJ_WT[rr] = n_infections / n_sequences / src.dat$sgtf_weights[rr]
#
#   # Impute by HHS region for states with missing testing data
#   if (is.na(src.dat$SIMPLE_ADJ_WT[rr])) {
#     # (estimated) total number of infections
#     # (rate assumed to be same as HHS region as a whole)
#     n_infections = src.dat$state_population[rr] * src.dat$HHS_INCIDENCE[rr]
#
#     # number of sequences (samples) in a given week & state
#     n_sequences = seq.tbl[src.dat$STUSAB[rr], src.dat$yr_wk[rr]]
#
#     # survey/sample weight
#     src.dat$SIMPLE_ADJ_WT[rr] = n_infections / n_sequences / src.dat$sgtf_weights[rr]
#   }
# }

### survey weights -------------------------------------------------------------
# Replacing for loop with faster data.table code
# see: https://git.biotech.cdc.gov/sars2seq/sc2_proportion_modeling/-/issues/6#note_86543
src.dat = data.table::data.table(src.dat)

# calculate simple adjusted weights using formula from https://www.cdc.gov/mmwr/volumes/70/wr/mm7023a3.htm
src.dat[, "SAW" := sqrt(state_population/TOTAL) * POSITIVE/sum(1/sgtf_weights)/sgtf_weights,
        .(STUSAB, yr_wk)]
# calculate simple adjusted weights for states that lack testing data
# (use average values from HHS region)
src.dat[, "SAW_ALT" := state_population*HHS_INCIDENCE / sum(1/sgtf_weights)/sgtf_weights,
        .(STUSAB, yr_wk)]
# use "SAW_ALT" when SAW is NA
src.dat[, "SIMPLE_ADJ_WT" := ifelse(is.na(SAW), SAW_ALT, SAW)]
# remove "SAW" and "SAW_ALT" columns
src.dat[, c("SAW", "SAW_ALT") := .(NULL, NULL)]

# optionally use group weights (where some recent weeks are grouped together)
if(use_group_weights){
  # calculate simple adjusted weights using formula from https://www.cdc.gov/mmwr/volumes/70/wr/mm7023a3.htm
  src.dat[, "SAW" := sqrt(group_population/TOTAL_gp) * POSITIVE_gp/sum(1/sgtf_weights_gp)/sgtf_weights_gp,
          .(STUSAB, group)]
  # calculate simple adjusted weights for states that lack testing data
  # (use average values from HHS region)
  src.dat[, "SAW_ALT" := group_population*HHS_INCIDENCE_gp / sum(1/sgtf_weights_gp)/sgtf_weights_gp,
          .(STUSAB, group)]
  # use "SAW_ALT" when SAW is NA
  src.dat[, "SIMPLE_ADJ_WT" := ifelse(is.na(SAW), SAW_ALT, SAW)]
  # remove "SAW" and "SAW_ALT" columns
  src.dat[, c("SAW", "SAW_ALT") := .(NULL, NULL)]
}

# check for weights of 0 (possible if there was very little testing in a state and no tests were positive)
# If using weight trimming, these are replaced with the minimum weight > 0.
zero_weights <- src.dat[SIMPLE_ADJ_WT == 0,]
if(nrow(zero_weights)>0 & !trim_weights) {
  warning('Some sequences have weights of 0: ')
}
if(nrow(zero_weights)>0 & !trim_weights){
  print(zero_weights)
}

# check for weights of NA weights
na_weights <- src.dat[is.na(SIMPLE_ADJ_WT),]
if(nrow(na_weights)>0) {
  warning('Some sequences have weights of NA: ')
}
if(nrow(na_weights)>0){
  print(na_weights)
}

# check for Infinite weights
inf_weights <- src.dat[is.infinite(SIMPLE_ADJ_WT),]
if(nrow(inf_weights)>0) {
  warning('Some sequences have infinite weights: ')
}
if(nrow(inf_weights)>0){
  print(inf_weights)
}

# Remove NA and INF weights
src.dat = subset(x = src.dat,
                 !is.na(SIMPLE_ADJ_WT) &
                   SIMPLE_ADJ_WT < Inf)



### aggregate sublineages ------------------------------------------------------
# make sure "VARIANT" is a character (rather than factor)
# (redundant with "stringsAsFactors = FALSE")
src.dat$VARIANT = as.character(src.dat$VARIANT)

#Identify all the clades/lineages to aggregate in the surveillance dataset
# all the AY variants
AY = sort(unique(src.dat$VARIANT)[grep("AY",unique(src.dat$VARIANT))])
# just the AY variants that are to be aggregated (i.e. not listed in "voc")
AY = AY[which(AY %notin% voc)]
P1=sort(unique(src.dat$VARIANT)[grep("P\\.1.",unique(src.dat$VARIANT))])
P1=P1[which(P1 %notin% voc)] #vector of the P1s to aggregate
Q=sort(unique(src.dat$VARIANT)[grep("Q\\.",unique(src.dat$VARIANT))])
Q=Q[which(Q %notin% voc)] #vector of the Qs to aggregate
B351=sort(unique(src.dat$VARIANT)[grep("B\\.1\\.351\\.",unique(src.dat$VARIANT))])
B351=B351[which(B351 %notin% voc)] #vector of the B351s to aggregate
B621=sort(unique(src.dat$VARIANT)[grep("B\\.1\\.621\\.",unique(src.dat$VARIANT))])
B621=B621[which(B621 %notin% voc)] #vector of the B621s to aggregate
B429=sort(unique(src.dat$VARIANT)[grep("B\\.1\\.429",unique(src.dat$VARIANT))])
B429=B429[which(B429 %notin% voc)] #vector of the B429s to aggregate
# omicrons [including lots of sublineages]
# aggregate the BA sublineages [NOTE! this *WILL* aggregate all BA.1.1.x sub-sublineages into BA.1.1 *even* if BA.1.1.x is listed in voc.]
if('BA.1' %in% voc) B529.BA1 <- sort(grep("(BA\\.1)(?!(\\.1$)|(\\.1\\.))",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA1 <- NULL
if('BA.1.1' %in% voc) B529.BA1.1 <- sort(grep("(BA\\.1\\.1)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T))     else B529.BA1.1 <- NULL
if('BA.2' %in% voc) B529.BA2 <- sort(grep("(BA\\.2)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T))             else B529.BA2 <- NULL
if('BA.3' %in% voc) B529.BA3 <- sort(grep("(BA\\.3)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T))             else B529.BA3 <- NULL
if('BA.4' %in% voc) B529.BA4 <- sort(grep("(BA\\.4)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T))             else B529.BA4 <- NULL
if('BA.5' %in% voc) B529.BA5 <- sort(grep("(BA\\.5)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T))             else B529.BA5 <- NULL
# safety check: make sure that no variants are in the multiple sublineage groups 
if(any(duplicated(c(B529.BA1, B529.BA1.1, B529.BA2, B529.BA3, B529.BA4, B529.BA5)))) errorCondition(message = paste0(c(B529.BA1, B529.BA1.1, B529.BA2, B529.BA3, B529.BA4, B529.BA5)[duplicate(c(B529.BA1, B529.BA1.1, B529.BA2, B529.BA3, B529.BA4, B529.BA5))], ' appear in multiple BA sublineage groups. Check B529.BA1, B529.BA1.1, B529.BA2, B529.BA3, B529.BA4, B529.BA5.'))
B529=sort(grep("(B\\.1\\.1\\.529)|(BA\\.[0-9])",unique(src.dat$VARIANT), value = T))
B529=B529[ B529 %notin% c(voc, B529.BA1, B529.BA1.1, B529.BA2, B529.BA3, B529.BA4, B529.BA5) ] #vector of the B529s to aggregate


# Aggregate sublineages to the parent lineage
if(P.1_agg==TRUE)     {src.dat[src.dat$VARIANT %in% P1,  "VARIANT"]    <- "P.1"}
if(B.1.351_agg==TRUE) {src.dat[src.dat$VARIANT %in% B351,"VARIANT"]    <- "B.1.351"}
if(B.1.621_agg==TRUE) {src.dat[src.dat$VARIANT %in% B621,"VARIANT"]    <- 'B.1.621' }
if(Q.1_3_agg==TRUE)   {src.dat[src.dat$VARIANT %in% Q,   "VARIANT"]    <- "B.1.1.7"}
if(AY_agg==TRUE)      {src.dat[src.dat$VARIANT %in% AY,  "VARIANT"]    <- "B.1.617.2"}
if(B429_7_agg==TRUE)  {src.dat[src.dat$VARIANT %in% B429,"VARIANT"]    <- "B.1.427"}
if(B.1.1.529_agg==TRUE)  {
   src.dat[src.dat$VARIANT %in% B529,"VARIANT"] <- "B.1.1.529"
   src.dat[src.dat$VARIANT %in% B529.BA1,"VARIANT"] <- "BA.1"
   src.dat[src.dat$VARIANT %in% B529.BA1.1,"VARIANT"] <- "BA.1.1"
   src.dat[src.dat$VARIANT %in% B529.BA2,"VARIANT"] <- "BA.2"
   src.dat[src.dat$VARIANT %in% B529.BA3,"VARIANT"] <- "BA.3"
   src.dat[src.dat$VARIANT %in% B529.BA4,"VARIANT"] <- "BA.4"
   src.dat[src.dat$VARIANT %in% B529.BA5,"VARIANT"] <- "BA.5"
}

# create another column for the varients of interest
# this is only used to get (unweighted) counts of the sequences by lineage (used in all runs)
src.dat$VARIANT2 = as.character(src.dat$VARIANT)
# group all non-"variants of interest" together
src.dat[src.dat$VARIANT %notin% voc, "VARIANT2"] <- "Other"



### survey design --------------------------------------------------------------
# (used in "Not Run3" and "Run3" below. The trimmed weights that are calculated
#  from svyDES are (optionally) used in Run2.)
if(svy.type=="svyREG"){
  svyDES = survey::svydesign(ids     = ~STUSAB+SOURCE,
                             strata  = ~HHS,
                             weights = ~SIMPLE_ADJ_WT,
                             nest    = TRUE, # TRUE = disable checking for
                             #                 duplicate cluster ID's across strata
                             data    = src.dat)
  # fpc not specified = sampling with replacement, which is equivalent to:
  # ids = ~ STUSAB; it's most appropriate when sampling relatively small
  # proportion of the population
} else {
  svyDES = survey::svydesign(ids     = ~SOURCE,
                             strata  = ~STUSAB + yr_wk,
                             weights = ~SIMPLE_ADJ_WT,
                             nest    = TRUE,
                             data    = src.dat)
}

# Add a column for weights that can be either SIMPLE_ADJ_WT or trimmed weights
src.dat$weights = src.dat$SIMPLE_ADJ_WT

# trim weights using survey package:
if(trim_weights){

  # get the max value to trim weights down to
  if( grepl(pattern = 'quant', x = opts$trim_weights, ignore.case = TRUE) ){
    quantile <- as.numeric(sub(pattern = '([^0-9]*)([0-9]+)',
                               replacement = '\\2',
                               x = opts$trim_weights))
    max_weight <- quantile(weights(svyDES),
                           probs = ifelse(quantile>1, quantile/100,quantile))
  } else if( grepl(pattern = 'IQR', x = opts$trim_weights, ignore.case = TRUE)){
    max_weight <- median(weights(svyDES)) + 6 * IQR(weights(svyDES))
  }

  # get the minimum weight to replace weights of 0
  wts <- weights(svyDES)
  min_wt <- min(wts[wts > 0])

  # trim weights using the min & max weights calculated above
  svyDES <- survey::trimWeights(design = svyDES,
                                upper = max_weight,
                                lower = min_wt)

  # add weights back into the dataframe
  myweights <- weights(svyDES)
  src.dat$wts_trimmed = myweights

  # Set "weights" to be the trimmed weights
  src.dat$weights = myweights
}



### model weeks ----------------------------------------------------------------
# add in another column for model_week to use in the model
# the reason for this is b/c large values of "week" were causing non-invertible hessians in the nowcast model.
# could scale "week" using mean and SD. This is just another option.
# NOTE! When switching to this method, I also switched to fitting the model to
#       20 weeks instead of 21 weeks.
# maximum week (not maximum "model_week") included in the model
model_week_max = as.numeric(as.Date(time_end) - week0day1) %/% 7
# first week (not first "model_week") included in the model
model_week_min = model_week_max - model_weeks
# a midpoint week that will be used to center week values
# this is also the maximum value of "model_week"
model_week_mid = round(model_weeks/2) # center it around 0 instead of using 1:model_weeks
# create a dataframe of old and new values to help visualize things
model_week_df = data.frame(week = 1:model_weeks + model_week_min,
                           model_week = 1:model_weeks - model_week_mid,
                           week_start = (1:model_weeks + model_week_min)*7 + as.Date(week0day1),
                           week_mid   = (1:model_weeks + model_week_min)*7 + as.Date(week0day1) + 3,
                           week_end   = (1:model_weeks + model_week_min)*7 + as.Date(week0day1) + 6)
# add the new model_week info to the dataframe
src.dat$model_week = src.dat$week - model_week_min - model_week_mid

# a function to convert dates to model_week
# (used in making predictions with Nowcast model)
date_to_model_week = function(date){
  # whole week
  # week = as.numeric(as.Date(date) - week0day1) %/% 7

  # fractional week (centered on Wednesday)
  week = as.numeric(as.Date(date) - (week0day1+3)) / 7
  return(week - model_week_min - model_week_mid)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # Background
#
# This brief summarizes findings from the variant surveillance system for variants
#  of concern and variants of interest.
#
# The sequences are those submitted to NS3 and by contract with lab vendors. Both
# raw numbers and weighted estimates are displayed. Weights allow estimates to be
# representative of all infected persons, aggregated by week of specimen collection
# and by state. Currently, the probability of selection for each sample is assumed
# to be independent of data source, but that may be modified as more data on
# sampling protocol become available. An adjustment is applied to account for
# oversampling of SGTF specimens in certain data streams.


# The general process for the runs below:
# 1) subset the data & survey design
# 2) generate estimates for individual regions/states
#    - create a dataframe with all combinations of time period, region, and
#      variant for which you want estimates
#    - calculate proportions & CI using survey design
# 3) generate estimates for entire US
#    - create a dataframe with all combinations of time period & variant for
#      which you want estimates
#    - calculate proportions & CI using survey design
# 4) combine regional & national estimates
# 5) add raw counts of sequences
# 6) add NCHS flags to identify estimates that may not be reliable
# 7) save results


# Runs 1 & 2 (i.e. not 3) ------------------------------------------------------
# calculate the variant share/proportion and confidence interval using survey
# design ("myciprop" function) for:
# - fortnights (HHS regions & nationally)
# - weeks  (HHS regions & nationally)
# creates 2 output files:
# - "results/variant_share_weighted_",       ci.type,"CI_",svy.type,"_",data_date,tag,".csv"
# - "results/variant_share_weekly_weighted_",ci.type,"CI_",svy.type,"_",data_date,tag,".csv"
if ( !grepl("Run3", tag) ){ # fortnight and weekly estimates

  # subset data to only include
  # - data since May 8, 2021
  # - older than "time_end"
  dat2 <- subset(x = src.dat,
                 as.Date(FORTNIGHT_END) >= as.Date("2021-05-08") &
                   as.Date(FORTNIGHT_END) <= time_end)

  # get the relevant fortnights from the data
  # (should this be sorted?)
  ftnts = (unique(dat2$FORTNIGHT_END))

  # skip running this if only running/testing the nowcast model.
  # (This section takes 90+% of the time required to run this script.)
  if(!nowcast_only){
    # create a dataframe with all unique combinations of variants, fortnights, and regions
    # (and then reverse column order for convenience)
    all.ftnt = expand.grid(Variant = voc,
                           Fortnight_ending = ftnts,
                           USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

    # get the proportion estimates & CI
    ests = apply(X = all.ftnt,
                 MARGIN = 1,
                 FUN = function(rr) myciprop(voc = rr[3],
                                             geoid = rr[1],
                                             svy = subset(x = svyDES,
                                                          FORTNIGHT_END == rr[2]),
                                             str = FALSE))

    # add in the estimates to the dataframe
    all.ftnt = cbind(all.ftnt,
                     Share    = ests[1,],
                     Share_lo = ests[2,],
                     Share_hi = ests[3,],
                     DF       = ests[4,],
                     eff.size = ests[5,],
                     cv.mean  = ests[6,],
                     deff     = ests[7,])

    ## make predictions for the "other" variants (following the same steps)
    others = expand.grid(Variant          = "Other",
                         Fortnight_ending = ftnts,
                         USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

    # get the proportion estimates & CI (for "other" variants)
    ests.others = apply(X = others,
                        MARGIN = 1,
                        FUN = function(rr) myciprop(voc = voc,
                                                    geoid = rr[1],
                                                    svy = subset(svyDES,
                                                                 FORTNIGHT_END == rr[2]),
                                                    str = FALSE))

    # add in the estimates to the dataframe (for "other" variants)
    others = cbind(others,
                   Share    = 1-ests.others[1,],
                   Share_lo = 1-ests.others[3,],
                   Share_hi = 1-ests.others[2,],
                   DF       = ests.others[4,],
                   eff.size = ests.others[5,],
                   cv.mean  = ests.others[6,],
                   deff     = ests.others[7,])

    # combine the estimates for the vocs with the estimates for "other" variants
    all.ftnt = rbind(all.ftnt,
                     others)

    # create a table of counts by variant, time period, and HHS region
    raw_counts_REG <- aggregate(count ~ VARIANT2 + FORTNIGHT_END + HHS,
                                data = dat2,
                                FUN  = sum,
                                drop = FALSE) # drop = FALSE: keep all combinations, even if no observations

    # convert HHS region to character
    raw_counts_REG$HHS <- as.character(raw_counts_REG$HHS)

    # create a table of counts by variant and time period for the whole US
    raw_counts_US <- aggregate(count ~ VARIANT2 + FORTNIGHT_END,
                               data = dat2,
                               FUN  = sum,
                               drop = FALSE)

    # add a column for HHS region
    raw_counts_US <- cbind(raw_counts_US[,1:2],
                           HHS = "USA",
                           count = raw_counts_US[,3])

    # combine dataframe of counts by region with dataframe of counts for US
    raw_counts <- rbind.data.frame(raw_counts_US,
                                   raw_counts_REG)

    # merge weighted proportions estimates with sequence counts
    all.ftnt2 <- merge(x = all.ftnt,
                       y = raw_counts,
                       by.x = c("USA_or_HHSRegion",
                                "Fortnight_ending",
                                "Variant"),
                       by.y = c("HHS",
                                "FORTNIGHT_END",
                                "VARIANT2"),
                       all = T)

    # replace NA counts with 0
    all.ftnt2[is.na(all.ftnt2$count)==T, "count"] <- 0

    #calculate denominator counts by region & time period
    dss <- aggregate(count ~ USA_or_HHSRegion + Fortnight_ending,
                     data = all.ftnt2,
                     FUN  = sum)

    # change the names of the "count" column
    names(dss)[grep("count",names(dss))] <- "denom_count"

    # add the denominator counts into the dataframe of results
    all.ftnt2 <- merge(x = all.ftnt2,
                       y = dss)

    #set the Share 0 and CI limits to NA when the count for a lineage is 0
    all.ftnt2$Share = ifelse(test = all.ftnt2$Share != 0 & all.ftnt2$count == 0,
                             yes = 0,
                             no = all.ftnt2$Share)
    all.ftnt2$Share_lo = ifelse(test = is.na(all.ftnt2$Share_lo) == F & all.ftnt2$count == 0,
                                yes = NA,
                                no = all.ftnt2$Share_lo)
    all.ftnt2$Share_hi = ifelse(test = is.na(all.ftnt2$Share_hi) == F & all.ftnt2$count==0,
                                yes = NA,
                                no = all.ftnt2$Share_hi)

    # calculate absolute CI width
    all.ftnt2$CI_width = all.ftnt2$Share_hi - all.ftnt2$Share_lo

    ## generate NCHS flags
    # flag estimates with Degrees of Freedom < 8
    all.ftnt2$flag_df = as.numeric(all.ftnt2$DF < 8)
    # flag estimates with effective size of < 30 (or NA)
    all.ftnt2$flag_eff.size = ifelse(test = all.ftnt2$eff.size < 30 |
                                       is.na(all.ftnt2$eff.size) == T,
                                     yes = 1,
                                     no = 0)
    # flag estimates with "denominator count" of < 30 (or NA) (denominator = count of sequences of all variants in a given region and time period)
    all.ftnt2$flag_dss = ifelse(test = all.ftnt2$denom_count < 30 |
                                  is.na(all.ftnt2$denom_count) == T,
                                yes = 1,
                                no = 0)
    # flag estimates with wide (absolute) confidence intervals
    all.ftnt2$flag_abs.ciw = ifelse(test = all.ftnt2$CI_width > 0.30 |
                                      is.na(all.ftnt2$CI_width) == T,
                                    yes = 1,
                                    no = 0)
    # flag estimates with wide (relative) confidence intervals
    all.ftnt2$flag_rel.ciw = ifelse(test = ((all.ftnt2$CI_width/all.ftnt2$Share)*100) > 130 |
                                      is.na((all.ftnt2$CI_width/all.ftnt2$Share)*100) == T,
                                    yes = 1,
                                    no = 0)

    # Single identifier for observations that have *any* NCHS flag
    all.ftnt2$nchs_flag = ifelse(test = all.ftnt2$flag_df == 1 |
                                   all.ftnt2$flag_eff.size == 1 |
                                   all.ftnt2$denom_count == 1 |
                                   all.ftnt2$flag_abs.ciw == 1 |
                                   all.ftnt2$flag_rel.ciw == 1,
                                 yes = 1,
                                 no = 0)
    # Single identifier for observations that have any NCHS flag *other* than the
    # degrees of freedom flag.
    all.ftnt2$nchs_flag_wodf = ifelse(test = all.ftnt2$flag_eff.size == 1 |
                                        all.ftnt2$denom_count == 1 |
                                        all.ftnt2$flag_abs.ciw == 1 |
                                        all.ftnt2$flag_rel.ciw == 1,
                                      yes = 1,
                                      no = 0)

    # select columns for the final results
    all.ftnt2 = all.ftnt2[, c("USA_or_HHSRegion",
                              "Fortnight_ending",
                              "Variant",
                              "Share",
                              "Share_lo",
                              "Share_hi",
                              "count",
                              "denom_count",
                              "DF",
                              "eff.size",
                              "CI_width",
                              "nchs_flag",
                              "nchs_flag_wodf")]
    
    # optionally calculate the number of infections attributable to each variant
    if (calc_confirmed_infections){
      test_filepath <- paste0(script.basename, 
                         "/data/backup_", 
                         data_date, "/", 
                         data_date, "_tests_aggregated", 
                         custom_tag, ".RDS")
      
      if (file.exists(test_filepath)){
        test_list <- readRDS(file = test_filepath)
      
        # get the fortnightly test tallies & aggregate them by fn across USA
        tests_fn_us <- test_list$tests_fortnight[,
                                                  .('total_test_positives' = sum(POSITIVE, na.rm = T)), 
                                                  by = 'fortnight_end'][,'HHS' := 'USA']
        # aggregate fortnightly tests by HHS region
        tests_fn_hhs <- test_list$tests_fortnight[,
                                                  .('total_test_positives' = sum(POSITIVE, na.rm = T)), 
                                                  by = c('fortnight_end', 'HHS')]
        
        # merge the positive test results in with the variant proportion estimates
        all.ftnt2 <- merge(
          x = all.ftnt2, 
          y = rbind(tests_fn_us,
                    tests_fn_hhs), 
          by.x = c("USA_or_HHSRegion",
                   "Fortnight_ending"), 
          by.y = c('HHS',
                   'fortnight_end'),
          all.x = TRUE)
        
        # calculate case totals for each variant
        all.ftnt2$cases    <- all.ftnt2$total_test_positives * all.ftnt2$Share
        all.ftnt2$cases_lo <- all.ftnt2$total_test_positives * all.ftnt2$Share_lo
        all.ftnt2$cases_hi <- all.ftnt2$total_test_positives * all.ftnt2$Share_hi
      } else {
        print(paste0('File ', 
                     test_filepath, 
                     ' not found. Not calculating number of infections attributable to each variant for fortnights.'))
      }
      
    }

    # re-order results by HHS region [so that "other" variants are not listed seperately]
    all.ftnt2 <- all.ftnt2[order(all.ftnt2$USA_or_HHSRegion),]

    # write results to file
    write.csv(x = all.ftnt2,
              file = paste0(script.basename,
                            "/results/variant_share_weighted_",
                            ci.type,
                            "CI_",
                            svy.type,
                            "_",
                            data_date,
                            tag,
                            ".csv"),
              row.names = FALSE)


    ### Weekly estimates
    # (repeat of above code, but with weeks instead of fortnights)

    # subset data to only include data since 2 May, 2021
    # and (end of week) older than "time_end"
    dat2 <- subset(src.dat,
                   as.Date(yr_wk) >= as.Date("2021-05-02") &
                     (as.Date(yr_wk) + 6) <= time_end)

    # add a column for the date of the final day of each week
    dat2$WEEK_END = as.Date(dat2$yr_wk) + 6

    # get the unique weeks that are in the appropriate time frame
    wks = sort(unique(dat2$yr_wk))

    # create a dataframe with variant, time period, and region
    # (and reorder columns for convenience)
    all.wkly = expand.grid(Variant = voc,
                           Week_of = wks,
                           USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

    # get predictions & CI for each variant in each week and region
    # (using survey design but NOT multinomial nowcast model)
    ests = apply(X = all.wkly,
                 MARGIN = 1,
                 FUN = function(rr) myciprop(voc   = rr[3],
                                             geoid = rr[1],
                                             svy   = subset(svyDES, yr_wk == rr[2]),
                                             str   = FALSE))

    # add the predictions into the dataframe
    all.wkly = cbind(all.wkly,
                     Share    = ests[1,],
                     Share_lo = ests[2,],
                     Share_hi = ests[3,],
                     DF       = ests[4,],
                     eff.size = ests[5,],
                     cv.mean  = ests[6,])

    # get predictions for the non-focal ("other") variants grouped together
    # (and reorder columns for convenience)
    others = expand.grid(Variant = "Other",
                         Week_of = wks,
                         USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

    # get predictions & CI for "other" (grouped) variant in each week and region
    # (using survey design but NOT multinomial nowcast model)
    ests.others = apply(X = others,
                        MARGIN = 1,
                        FUN = function(rr) myciprop(voc = voc,
                                                    geoid = rr[1],
                                                    svy = subset(svyDES, yr_wk == rr[2]),
                                                    str =  FALSE))

    # add the predictions into the dataframe
    others = cbind(others,
                   Share    = 1-ests.others[1,],
                   Share_lo = 1-ests.others[3,],
                   Share_hi = 1-ests.others[2,],
                   DF       = ests.others[4,],
                   eff.size = ests.others[5,],
                   cv.mean  = ests.others[6,])

    # combine estimates for individual variants with the estimates for the "other"
    # (grouped) variants
    all.wkly = rbind(all.wkly,
                     others)

    # Add a column for the last day of each week
    all.wkly$WEEK_END = as.Date(all.wkly$Week_of) + 6


    ## generate sequence counts by lineage, location and date
    # raw counts of each variant in each week and each region
    raw_counts_REG <- aggregate(count ~ VARIANT2 + WEEK_END + HHS,
                                data = dat2,
                                FUN  = sum,
                                drop = FALSE)

    # convert region names to character
    raw_counts_REG$HHS <- as.character(raw_counts_REG$HHS)

    # raw counts for each variant in each week nationally
    raw_counts_US <- aggregate(count ~ VARIANT2 + WEEK_END,
                               data = dat2,
                               FUN  = sum,
                               drop = FALSE)

    # add in a column for HHS region = "USA"
    raw_counts_US <- cbind(raw_counts_US[,1:2],
                           HHS   = "USA",
                           count = raw_counts_US$count)

    # combine the regional counts with the national counts
    raw_counts <- rbind.data.frame(raw_counts_US,
                                   raw_counts_REG)

    # merge sequence counts with weighted proportions estimates
    all.wkly2 <- merge(x = all.wkly,
                       y = raw_counts,
                       by.x = c("USA_or_HHSRegion",
                                "WEEK_END",
                                "Variant"),
                       by.y = c("HHS",
                                "WEEK_END",
                                "VARIANT2"),
                       all = T)

    # Set counts to 0 if current count is NA and variant != "Other"
    all.wkly2[is.na(all.wkly2$count)==T ,"count"] <- 0

    # calculate denominator counts
    # (i.e. raw counts per region & week)
    dss <- aggregate(count ~ USA_or_HHSRegion + WEEK_END,
                     data = all.wkly2,
                     FUN  = sum)

    # change name of denominator counts
    names(dss)[grep("count",names(dss))] <- "denom_count"

    # add denominator counts into the dataframe of raw counts
    all.wkly2 <- merge(x = all.wkly2,
                       y = dss)

    #calculate absolute CI width
    all.wkly2$CI_width = all.wkly2$Share_hi - all.wkly2$Share_lo

    #set the Share 0 and CI limits to NA when the count for a lineage is 0
    all.wkly2$Share=ifelse(all.wkly2$Share!=0 & all.wkly2$count==0,0,all.wkly2$Share)
    all.wkly2$Share_lo=ifelse(is.na(all.wkly2$Share_lo)==F & all.wkly2$count==0,NA,all.wkly2$Share_lo)
    all.wkly2$Share_hi=ifelse(is.na(all.wkly2$Share_hi)==F & all.wkly2$count==0,NA,all.wkly2$Share_hi)

    ## generate NCHS flags
    # flag estimates with Degrees of Freedom < 8
    all.wkly2$flag_df = ifelse(test = all.wkly2$DF < 8,
                               yes = 1,
                               no = 0)
    # flag estimates with effective size of < 30 (or NA)
    all.wkly2$flag_eff.size = ifelse(test = all.wkly2$eff.size < 30 |
                                       is.na(all.wkly2$eff.size) == T,
                                     yes = 1,
                                     no = 0)
    # flag estimates with "denominator count" of < 30 (or NA)
    # (denominator = count of sequences of all variants in a given region and time period)
    all.wkly2$flag_dss = ifelse(test = all.wkly2$denom_count < 30|
                                  is.na(all.wkly2$denom_count) == T,
                                yes = 1,
                                no = 0)
    # flag estimates with wide (absolute) confidence intervals
    all.wkly2$flag_abs.ciw = ifelse(test = all.wkly2$CI_width > 0.30 |
                                      is.na(all.wkly2$CI_width) == T,
                                    yes = 1,
                                    no = 0)
    # flag estimates with wide (relative) confidence intervals
    all.wkly2$flag_rel.ciw = ifelse(test = ((all.wkly2$CI_width/all.wkly2$Share)*100) > 130 |
                                      is.na((all.wkly2$CI_width/all.wkly2$Share)*100) == T,
                                    yes = 1,
                                    no = 0)

    # Single identifier for observations that have *any* NCHS flag
    all.wkly2$nchs_flag = ifelse(test = all.wkly2$flag_df == 1 |
                                   all.wkly2$flag_eff.size == 1 |
                                   all.wkly2$denom_count == 1 |
                                   all.wkly2$flag_abs.ciw == 1 |
                                   all.wkly2$flag_rel.ciw == 1,
                                 yes = 1,
                                 no = 0)
    # Single identifier for observations that have any NCHS flag *other* than the
    # degrees of freedom flag.
    all.wkly2$nchs_flag_wodf = ifelse(all.wkly2$flag_eff.size == 1 |
                                        all.wkly2$denom_count == 1 |
                                        all.wkly2$flag_abs.ciw == 1 |
                                        all.wkly2$flag_rel.ciw == 1,
                                      yes = 1,
                                      no = 0)

    # flag estimates with raw count of < 20
    all.wkly2$count_LT20 = ifelse(test = all.wkly2$count < 20,
                                  yes = 1,
                                  no = 0)
    # flag estimates with raw count of < 10
    all.wkly2$count_LT10 = ifelse(test = all.wkly2$count < 10,
                                  yes = 1,
                                  no = 0)

    # select the columns to save
    all.wkly2 = all.wkly2[,c("USA_or_HHSRegion",
                             "WEEK_END",
                             "Variant",
                             "Share",
                             "Share_lo",
                             "Share_hi",
                             "count",
                             "denom_count",
                             "DF",
                             "eff.size",
                             "CI_width",
                             "nchs_flag",
                             "nchs_flag_wodf",
                             "count_LT20",
                             "count_LT10")]

    # optionally calculate the number of infections attributable to each variant
    if (calc_confirmed_infections){
      test_filepath <- paste0(script.basename, 
                              "/data/backup_", 
                              data_date, "/", 
                              data_date, "_tests_aggregated", 
                              custom_tag, ".RDS")
      
      if (file.exists(test_filepath)){
        test_list <- readRDS(file = test_filepath)
        
        # get the fortnightly test tallies & aggregate them by fn across USA
        tests_wk_us <- test_list$tests_weekly[,
                                              .('total_test_positives' = sum(POSITIVE, na.rm = T)), 
                                              by = 'yr_wk'][,'HHS' := 'USA']
        # aggregate fortnightly tests by HHS region
        tests_wk_hhs <- test_list$tests_weekly[,
                                               .('total_test_positives' = sum(POSITIVE, na.rm = T)), 
                                               by = c('yr_wk', 'HHS')]
        
        # merge the positive test results in with the variant proportion estimates
        all.wkly2 <- merge(
          x = all.wkly2, 
          y = rbind(tests_wk_us,
                    tests_wk_hhs)[,'WEEK_END' := as.Date(yr_wk) + 6][, 'yr_wk' := NULL], 
          by.x = c("USA_or_HHSRegion",
                   "WEEK_END"), 
          by.y = c('HHS',
                   'WEEK_END'),
          all.x = TRUE)
        
        # calculate case totals for each variant
        all.wkly2$cases    <- all.wkly2$total_test_positives * all.wkly2$Share
        all.wkly2$cases_lo <- all.wkly2$total_test_positives * all.wkly2$Share_lo
        all.wkly2$cases_hi <- all.wkly2$total_test_positives * all.wkly2$Share_hi
      } else {
        print(paste0('File ', 
                     test_filepath, 
                     ' not found. Not calculating number of infections attributable to each variant for fortnights.'))
      }
      
    }
    
    # sort by HHS region
    all.wkly2 <- all.wkly2[order(all.wkly2$USA_or_HHSRegion),]

    # save the results to file
    write.csv(x = all.wkly2,
              file = paste0(script.basename,
                            "/results/variant_share_weekly_weighted_",
                            ci.type,
                            "CI_",
                            svy.type,
                            "_",
                            data_date,
                            tag,
                            ".csv"),
              row.names=FALSE)
  }
} # end run (not 3)

# Run2 (nowcast) ---------------------------------------------------------------
# Model-based smoothed trends in variant share: National and by HHS Region
# creates 5 + 33 output files:
# - "/results/wow_growth_variant_share",data_date,tag,".csv"
# - "/results/updated_nowcast_fortnightly_",data_date,"_state_tag_included_Run1.csv"
# - "/results/updated_nowcast_fortnightly_",data_date,"_state_tag_included_Run2.csv"
# - "/results/updated_nowcast_weekly_",data_date,"_state_tag_included_Run1.csv"
# - "/results/updated_nowcast_weekly_",data_date,"_state_tag_included_Run2.csv"
# - and 33 images
#   - "/results/wtd_shares_",data_date,"_","barplot_US",tag,".jpg"
#   - "/results/wtd_shares_",data_date,"_","projection_US",tag,".jpg"
#   - "/results/wtd_shares_",data_date,"_","growthrate_US",tag,".png"
#     Create the same 3 figures for each HHS region
#   - "/results/wtd_shares_",data_date,"_","barplot_HHS",hhs,tag,".jpg"
#   - "/results/wtd_shares_",data_date,"_","projection_HHS",hhs,tag,".jpg"
#   - "/results/wtd_shares_",data_date,"_","growthrate_HHS", hhs,tag,".jpg"
if ( grepl("Run2",tag) ){

  ## Data prep for nowcast model ----

  ### get vocs to include in the multinomial model ----
  # when using "reduced_vocs", just include the reduced_vocs, NOT the n_top
  if(reduced_vocs){
    # could exclude particular variants
    # all_tops <- all_tops[all_tops %notin% exclude_variants]
    # easier to just set variants based on the reduced vocs
    model_vars <- voc
  } else {

    # Based on sequences where the specimen was collected in the current or previous
    # `n_recent_weeks` weeks, the count of sequences, and variant share (weighted
    #  percent) are:
    #
    # Select display and model variants - this basically selects which variants are
    # the top variants to include in the model
    us_var = sort(
      prop.table(x = xtabs(formula = weights ~ VARIANT,
                           data = src.dat,
                           subset = src.dat$week >= current_week - n_recent_weeks)),
      decreasing = TRUE)

    # names of all the variants
    us_rank = names(us_var)

    # # counts of variants for each hhs region
    # hhs_var = prop.table(
    #   x = xtabs(weights ~ HHS + VARIANT,
    #             subset(src.dat,
    #                    week >= current_week - n_recent_weeks)),
    #   margin = 1)

    # # get the names of the variants for each hhs region (in order of abundance)
    # hhs_rank = apply(hhs_var, 1, function(rr) names(sort(rr, decreasing=TRUE)))

    # names of variants that are either in the n_top or are in "voc"
    # Ordered by national rank
    # these will be included in results
    model_vars = us_rank[us_rank %in% c(us_rank[1:n_top], voc)]

    # optionally make sure B is not in the model_vars
    if (force_aggregate_B) model_vars <- model_vars[model_vars %notin% 'B']
  }

  ### add variant ranks to src.dat ----
  # makes sure each seq is assigned a rank/number based on the weighted proportion
  # in the last few weeks (1 = most common)
  # these are treated as categories and used as the response variable in the
  # multinomial nowcast model
  src.dat$K_US <- match(x = src.dat$VARIANT,
                        table = model_vars)

  # for variants that are not in "model_vars", assign last number,
  # which will correspond to "Other"
  src.dat$K_US[is.na(src.dat$K_US)] = length(model_vars) + 1

  ### subset src.dat (only include the weeks that are in the multinomial model) ----
  src.moddat = subset(src.dat,
                      model_week %in% ((1:model_weeks)-model_week_mid))
  # this filter now incorporates time_end; the previous filter (week >= max(week) - model_weeks) did not.

  # scale the data weights to help avoid numerical overflow in the multinomial model
  if(rescale_model_weights){
    if (tolower(rescale_model_weights_by) == 'max') {
      src.moddat$wts <- src.moddat$weights / max(src.moddat$weights)
    }
    if (tolower(rescale_model_weights_by) == 'mean') {
      src.moddat$wts <- src.moddat$weights / mean(src.moddat$weights)
    }
    if (is.numeric(rescale_model_weights_by)){
      src.moddat$wts <- src.moddat$weights / rescale_model_weights_by
    }
  } else {
    src.moddat$wts <- src.moddat$weights
  }

  ### survey design ----
  # weight using "weights", which is either SIMPLE_ADJ_WT or wts_trimmed
  mysvy = survey::svydesign(ids     = ~SOURCE,
                            strata  = ~STUSAB + yr_wk,
                            weights = ~wts,
                            nest = TRUE, # TRUE = disable checking for duplicate cluster ID's across strata
                            data = src.moddat) # not specifying fpc signifies sampling with replacement

  # optionally save src.moddat that's been prepped for analysis
  if(save_datasets_to_file){
    saveRDS(object = src.moddat,
            file = paste0(script.basename,
                          '/results/src.moddat_', # save to results instead of 'data' folder
                          data_date,
                          tag,
                          '.RDS'))
  }



  ## Fit Nowcast model ----

  # produces results for Run1 & Run2
  # (summarizes the output of the nowcast model to only include variants in voc1 for the Run1 output)

  ### National model ----
  svymlm_us = svymultinom(mod.dat = src.moddat,
                          mysvy = mysvy,
                          fmla = formula("as.numeric(as.factor(K_US))  ~ model_week"),
                          model_vars = model_vars)

  ### Regional model ----
  svymlm_hhs = svymultinom(mod.dat = src.moddat,
                           mysvy = mysvy,
                           fmla = formula("as.numeric(as.factor(K_US)) ~ model_week + as.factor(HHS)"),
                           model_vars = model_vars)

  ## Plot results ----

  ### plot prep -----

  # choose which variants to display in the results
  if (display_option=="top7") {
    # most abundant 7 variants
    display_vars = head(model_vars, 7)
  } else {
    # all the vocs
    display_vars = voc
  }

  # get the indices of the variants that will be displayed in the results
  display_indices = which(model_vars %in% display_vars)
  # get the names of the variants to display (in correct order)
  display_vars = model_vars[display_indices]

  # Normalize survey weights
  # (used for making bar plots)
  src.dat$NORM_WTS = with(src.dat,
                          weights/sum(weights,
                                      na.rm = TRUE))

  # normalize so the total weights = # of sequences (instead of 1)
  src.dat$NORM_WTS = src.dat$NORM_WTS * sum(!is.na(src.dat$NORM_WTS))

  # Set up colors
  col.dk = hcl.colors(length(display_vars),
                      palette = "TealRose",
                      alpha = 0.8)
  names(col.dk) = display_vars

  # base filename for figures
  stub = paste0(script.basename,
                "/results/wtd_shares_",
                data_date,
                "_")


  # Weighted variant shares of the top variants in the past `display_lookback`
  # weeks (number of sequences collected weekly above each bar), and model-based
  # smoothed estimates, nationwide:

  ### barplot of weighted share (national) ----

  # tabulate/count the normalized survey weights by week, and variant
  # (only used for barplot below)
  bp_us = xtabs(formula = NORM_WTS ~ model_week + VARIANT,
                data = src.dat,
                subset = src.dat$model_week %in% ((model_week_mid - display_lookback + 1):model_week_mid))

  # Normalize values to proportions
  bp_us = prop.table(bp_us, 1)

  # subset the weights to only include the variants to be displayed
  bp_us = bp_us[, display_vars]

  # give the table prettier column names
  # rownames(bp_us) = week_label(as.numeric(rownames(bp_us)) - current_week)
  rownames(bp_us) = format(
    x = ((as.numeric(rownames(bp_us)) + model_week_mid) + model_week_min) * 7 + as.Date(week0day1),
    format = '%m-%d'
  )

  # create plot
  if (fig_gen_run) jpeg(filename  = paste0(stub, "barplot_US", tag, ".jpg"),
                        width     = 1500,
                        height    = 1500,
                        pointsize = 40)
  # create a barplot of "observed" values (i.e. weighted counts)
  bp = barplot(height = 100 * t(bp_us),
               xlab = "Week beginning",
               ylab = "Weighted variant share (%)",
               main = "Nationwide",
               border = NA,
               ylim = 110 * 0:1,
               col = col.dk,
               names.arg = rownames(bp_us),
               legend.text = display_vars,
               args.legend = list(x = "topleft",
                                  bty = "n",
                                  border = NA))

  # add text to the barplot
  text(x = bp,
       y = 3 + colSums(100 * t(tail(bp_us, 12))),
       labels = with(subset(src.dat,
                            week < current_week &
                              week >= current_week - display_lookback),
                     table(week)),
       cex = 0.7)
  if (fig_gen_run) dev.off()


  ### barplot of model-predicted data (national) ----
  # create a dataframe for predictions for each week
  pred_us.df = expand.grid(model_week = seq(from = -display_lookback,
                                            to = 2,
                                            by = (1/7)) + (model_weeks - model_week_mid))

  #add a column for (predicted) each variant proportion for each timepoint
  pred_us.df = cbind(pred_us.df,
                     predict(object = svymlm_us$mlm,
                             newdata = pred_us.df,
                             type = "probs"))

  # add in dates
  # pred_us.df$date       = as.Date((pred_us.df$model_week + model_week_mid + model_week_min) * 7 + week0day1)
  # pred_us.df$week_start = pred_us.df$date - as.numeric(format(pred_us.df$date, format = '%w'))
  # pred_us.df$week_end   = pred_us.df$week_start + 6
  # pred_us.df$week_mid   = pred_us.df$week_start + 3

  # Barplot of Nowcast predicted values
  if (fig_gen_run) jpeg(filename  = paste0(stub, "projection_US", tag, ".jpg"),
                        width     = 1500,
                        height    = 1500,
                        pointsize = 40)

  bp = barplot(height = 100 * t(pred_us.df[, 1 + display_indices]),
               xlab = "Week beginning",
               ylab = "Weighted variant share (%)",
               main = "Nationwide",
               space = 0,
               border = NA,
               ylim = 110 * 0:1,
               col = col.dk,
               names.arg = ifelse(test = pred_us.df$model_week %% 1 == 0,
                                  yes = format(((pred_us.df$model_week + model_week_mid) + model_week_min) * 7 + as.Date(week0day1), format = '%m-%d'),
                                  no = NA),
               legend.text = display_vars,
               args.legend = list(x = "topleft",
                                  bty = "n",
                                  border = NA))

  # predicted percent contributions of some variants for the current week
  # pc = unlist(100 * subset(pred_us.df,
  #                          week==current_week)[, 1+display_indices])
  pc = unlist(100 * subset(pred_us.df,
                           model_week == model_week_mid)[, 1+display_indices])

  # define a y-value for the text to be added to the plot
  y = cumsum(pc) - pc/2

  # add text to the plot
  text(x = 1.02 * tail(bp, 1),
       y = y,
       labels = round(pc, 1),
       cex = 0.7,
       xpd = TRUE,
       adj = c(0, 0.5))

  # define x-values for grey boxes signifying that recent data is likely incomplete
  # x = bp[which(pred_us.df$week %in% (current_week + c(-2, 0, 2)))]
  x = bp[which(pred_us.df$model_week %in% (model_week_mid + c(-2, 0, 2)))]

  # add grey rectangles to plot
  rect(xleft   = x[1],
       ybottom = 0,
       xright  = x[2],
       ytop    = 100,
       border  = NA,
       col = "#00000020")
  rect(xleft   = x[2],
       ybottom = 0,
       xright  = x[3],
       ytop    = 100,
       border  = NA,
       col = "#00000040")
  if (fig_gen_run) dev.off()

  ### WoW growth rate vs. transmission ----
  #  - the vertical axis depicts the variant share growth rate
  #    (derivative of log of variant proportion with respect to time).
  #  - The axis on the right shows an estimate of variant transmissibility with
  #    respect to the overall mean transmissibility.

  # Get the SE of the national estimate
  us.summary = se.multinom(mlm   = svymlm_us$mlm,
                           newdata_1row = data.frame(
                             model_week = model_week_mid,
                             HHS = "USA"))

  # calculate the SE of the estimated growth rate
  se.gr = with(data = us.summary,
               expr = 100 * exp((sqrt(se.b_i^2 * (1 - 2 * p_i) + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2)))) - 100) # coeff of se.b_i^2 corrected [2022-01-13]

  # p_i    = predicted probability (proportional representation/frequency)
  # se.p_i = SE of predicted frequency
  # b_i    = coefficient value (for intercept & week)
  # se.b_i = SE of coefficient value

  # calculate the estimated growth rate
  gr = with(data = us.summary,
            expr = 100 * exp(b_i - sum(p_i * b_i)) - 100)

  # calculate doubling times
  se.gr_link = with(data = us.summary,
                    expr = sqrt(se.b_i^2 * (1 - 2 * p_i) + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2)))
  gr_link = with(data = us.summary,
                 expr = (b_i - sum(p_i * b_i)))
  gr_lo_link = gr_link - 1.96 * se.gr_link
  gr_hi_link = gr_link + 1.96 * se.gr_link

  gr_lo = 100 * exp(gr_lo_link) - 100
  gr_hi = 100 * exp(gr_hi_link) - 100

  # calculate doubling time
  doubling_time = log(2)/gr_link * 7
  doubling_time_lo = log(2)/gr_lo_link * 7
  doubling_time_hi = log(2)/gr_hi_link * 7


  if (fig_gen_run) png(filename = paste0(stub, "growthrate_US", tag, ".png"),
                       width = 8,
                       height = 8,
                       units = "in",
                       pointsize = 16,
                       res = 1000)

  # make enough room on the right side of the plot for a secondary axis.
  orpar <- par()
  par(mar = c(5.1, 4.1, 4.1, 4.1))

  # plot "nowcast" groth rates by variant
  plot(x = 100 * us.summary$p_i,
       y = gr,
       log = "x",
       type = "n",
       ylim = range(gr + 1.96 * se.gr, gr - 1.96 * se.gr),
       xaxt = "n",
       xlim = c(0.0005,110),
       xlab = "Nowcast Estimated Proportion (%)",
       ylab = "Week over week growth rate (%)",
       main = "Nationwide")

  # Add an explanation for the lack of confidence intervals in the event that
  # the Hessian was non-invertible
  if(is.null(svymlm_us$SE)){
    mtext(text = "*No SE estimates b/c of non-invertible Hessian in multinomial model fit.",
          side = 3,
          line = 0,
          cex = 0.75,
          font = 4,
          col = 'red' )
  }

  # add an x axis
  axis(side = 1,
       at     = c(0.001, 0.01, 0.1, 1, 10, 100),
       labels = c(0.001, 0.01, 0.1, 1, 10, 100))

  # add a horizontal line at 0
  abline(h = 0,
         col = "grey65")

  # add lines for each variant
  for (vv in seq(model_vars)) {

    # vertical lines for uncertainty in growth rate
    lines(x = 100 * rep(us.summary$p_i[vv], 2),
          y = gr[vv] + 1.96 * c(1,-1) * se.gr[vv],
          col = "blue",
          lwd = 2)

    # horizontal lines for uncertainty in variant proportion
    lines(x = pmax(0.0001, 100 * us.summary$p_i[vv] + 196 * c(1,-1) * us.summary$se.p_i[vv]),
          y = rep(gr[vv], 2),
          col = "blue",
          lwd = 2)
  }

  # add the name of each variant
  text(x = 100 * us.summary$p_i,
       y = gr,
       labels = c(model_vars,""),
       cex = 0.85,
       col = "grey25",
       adj = 1.15)

  # add text for doubling time
  plot_growth_rates <- axTicks(2)
  plot_doubling_times <- (log(2) / log((100 + plot_growth_rates)/100)) * 7

  # Add second axis
  axis(side = 4,
       at = plot_growth_rates,
       labels = round(plot_doubling_times,1))
  # Add second axis label
  mtext("Doubling time (days)",
        side = 4,
        line = 3)

  if (fig_gen_run)  dev.off()

  # create a dataframe of variant shares & growth rates
  gr_tab = data.frame(variant          = c(model_vars, "OTHER"),
                      variant_share    = (100 * us.summary$p_i),
                      growth_rate      = gr,
                      growth_rate_lo   = gr_lo,
                      growth_rate_hi   = gr_hi,
                      doubling_time    = doubling_time,
                      doubling_time_lo = doubling_time_lo,
                      doubling_time_hi = doubling_time_hi,
                      model_week       = model_week_mid)

  # merge in date
  gr_tab <- merge(gr_tab,
                  model_week_df,
                  by = 'model_week')

  # save growth rates to file
  write.csv(x = gr_tab,
            file = paste0(script.basename,
                          "/results/wow_growth_variant_share",
                          data_date,
                          tag,
                          ".csv"),
            row.names = FALSE)


  ### regional barplots (weighted) ----

  # Same as above for each HHS region:
  # Weighted variant shares of the top variants in the past `display_lookback`
  # weeks (number of sequences collected weekly above each bar), and model-based
  # smoothed estimates, for each HHS region:


  # tabulate/count the normalized survey weights by region, week, and variant
  # (just used for barplots)
  bp_hhs = xtabs(formula = NORM_WTS ~ HHS + model_week + VARIANT,
                 data = src.dat,
                 subset = src.dat$model_week %in% ((model_week_mid - display_lookback + 1):model_week_mid))

  # subset the weights to only include the variants to be displayed
  bp_hhs = prop.table(bp_hhs, 1:2)[,, display_vars]

  # give the table prettier column names (week start date instead of model_week)
  # dimnames(bp_hhs)[[2]] = week_label(as.numeric(dimnames(bp_hhs)[[2]]) - current_week)
  dimnames(bp_hhs)[[2]] = format(
    x = ((as.numeric(dimnames(bp_hhs)[[2]]) + model_week_mid) + model_week_min) * 7 + as.Date(week0day1),
    format = '%m-%d'
  )


  # add region to the growth rate table (and then add on the HHS regional growth rates)
  gr_tab$region = 'USA'

  for (hhs in sort(unique(src.dat$HHS))) {

    # Barplot of regional data
    if (fig_gen_run) jpeg(filename = paste0(stub, "barplot_HHS", hhs,tag,".jpg"),
                          width = 1500,
                          height = 1500,
                          pointsize = 40)

    # create a barplot of "observed values (weighted proportions)
    bp = barplot(height = 100 * t(bp_hhs[hhs,,]),
                 xlab = "Week beginning",
                 ylab = "Weighted variant share (%)",
                 main = paste("HHS Region", hhs),
                 border = NA,
                 ylim = 110 * 0:1,
                 col = col.dk,
                 # names.arg = rownames(bp_hhs[hhs,,]) - current_week,
                 # names.arg = format((as.numeric(as.Date(rownames(bp_hhs[hhs,,]))) + model_week_mid + model_week_min) * 7 + week0day1, format = '%m-%d'),
                 names.arg = rownames(bp_hhs[hhs,,]),
                 legend.text = display_vars,
                 args.legend = list(x = "topleft",
                                    bty = "n",
                                    border = NA))

    # add text to the barplot
    text(x = bp,
         y = 3 + colSums(100 * t(tail(bp_hhs[hhs,,], 12))),
         labels = with(subset(src.dat,
                              week < current_week &
                                week >= current_week - display_lookback),
                       table(week)),
         cex = 0.7)
    if (fig_gen_run) dev.off()

    ### regional barplot (Nowcast) ----

    # create a dataframe for predictions for each week & HHS region
    pred_hhs.df = expand.grid(model_week = seq(from = -display_lookback,
                                               to = 2,
                                               by = (1/7)) + (model_weeks - model_week_mid), #  + current_week,
                              HHS = sort(unique(src.moddat$HHS)))

    # predict a sequence's clade (ID'd by rank abundance) given only time & hhs region
    # (i.e. predicted proportion for each clade in each hhs region & time period)
    pred_hhs.df = cbind(pred_hhs.df,
                        predict(object  = svymlm_hhs$mlm,
                                newdata = pred_hhs.df,
                                type    = "probs"))

    # add in dates
    # pred_hhs.df$date = as.Date((pred_hhs.df$model_week + model_week_mid + model_week_min) * 7 + week0day1)
    # pred_hhs.df$week_start = pred_hhs.df$date - as.numeric(format(pred_hhs.df$date, format = '%w'))
    # pred_hhs.df$week_end   = pred_hhs.df$week_start + 6
    # pred_hhs.df$week_mid   = pred_hhs.df$week_start + 3


    if (fig_gen_run) jpeg(filename = paste0(stub, "projection_HHS", hhs,tag, ".jpg"),
                          width = 1500,
                          height = 1500,
                          pointsize = 40)

    # subset predicted data to only include the region of interest
    pred.df = subset(pred_hhs.df,
                     HHS==hhs)[, -2]

    # Barplot of Nowcast predicted values
    bp = barplot(height = 100 * t(pred.df[, 1 + display_indices]),
                 xlab = "Week beginning",
                 ylab = "Weighted variant share (%)",
                 main = paste("HHS Region", hhs),
                 space = 0,
                 border = NA,
                 ylim = 110 * 0:1,
                 col = col.dk,
                 # names.arg = ifelse(pred.df$week %% 1 == 0,
                 #                    week_label(pred.df$week - current_week),
                 #                    NA),
                 names.arg = ifelse(pred.df$model_week %% 1 == 0,
                                    format(((pred.df$model_week + model_week_mid) + model_week_min) * 7 + as.Date(week0day1), format = '%m-%d'),
                                    NA),
                 legend.text = display_vars,
                 args.legend = list(x = "topleft",
                                    bty = "n",
                                    border = NA))

    # predicted percent contributions of some variants for the current week
    # pc = unlist(100 * subset(pred.df, week==current_week)[, 1+display_indices])
    pc = unlist(100 * subset(pred.df, model_week==model_week_mid)[, 1+display_indices])

    # define a y-value for the text on the plot
    y = cumsum(pc) - pc/2

    # add text to the plot
    text(x = 1.02 * tail(bp, 1),
         y = y,
         labels = round(pc, 1),
         cex = 0.7,
         xpd = TRUE,
         adj = c(0, 0.5))

    # define x-values for grey boxes signifying that recent data is likely incomplete
    # x = bp[which(pred.df$week %in% (current_week + c(-2, 0, 2)))]
    x = bp[which(pred.df$model_week %in% (model_week_mid + c(-2, 0, 2)))]

    # add boxes/rectangles to the plot
    rect(xleft = x[1],
         ybottom = 0,
         xright = x[2],
         ytop = 100,
         border = NA,
         col = "#00000020")
    rect(xleft = x[2],
         ybottom = 0,
         xright = x[3],
         ytop = 100,
         border = NA,
         col = "#00000040")
    if (fig_gen_run) dev.off()


    ### regional WoW growth rate -----

    # get the SE for the given model, week, and HHS region
    hhs.summary = se.multinom(mlm = svymlm_hhs$mlm,
                              newdata_1row = data.frame(
                                model_week = model_week_mid,
                                HHS = hhs
                              ),
                              composite_variant = NULL)

    # calculate the SE of the growth rate
    se.gr = with(data = hhs.summary,
                 expr = 100 * exp(sqrt(se.b_i^2 * (1 - 2 * p_i) + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2))) - 100)

    # calculate the growth rate
    gr = with(hhs.summary,
              100 * exp(b_i - sum(p_i * b_i)) - 100)

    # add in doubling times
    se.gr_link = with(data = hhs.summary,
                      expr = sqrt(se.b_i^2 * (1 - 2 * p_i) + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2)))
    gr_link = with(data = hhs.summary,
                   expr = (b_i - sum(p_i * b_i)))
    gr_lo_link = gr_link - 1.96 * se.gr_link
    gr_hi_link = gr_link + 1.96 * se.gr_link

    gr_lo = 100 * exp(gr_lo_link) - 100
    gr_hi = 100 * exp(gr_hi_link) - 100

    # calculate doubling time
    doubling_time = log(2)/gr_link * 7
    doubling_time_lo = log(2)/gr_lo_link * 7
    doubling_time_hi = log(2)/gr_hi_link * 7

    # create a dataframe of variant shares & growth rates
    gr_tab_hhs = data.frame(variant          = c(model_vars, "OTHER"),
                            variant_share    = (100 * hhs.summary$p_i),
                            growth_rate      = gr,
                            growth_rate_lo   = gr_lo,
                            growth_rate_hi   = gr_hi,
                            doubling_time    = doubling_time,
                            doubling_time_lo = doubling_time_lo,
                            doubling_time_hi = doubling_time_hi,
                            region           = hhs,
                            model_week       = model_week_mid
    )

    # merge in date
    gr_tab_hhs <- merge(gr_tab_hhs,
                        model_week_df,
                        by = 'model_week')

    # combine national and regional growth rates
    gr_tab = rbind(
      gr_tab,
      gr_tab_hhs
    )

    # create plot of growth rates by region
    if (fig_gen_run) jpeg(filename = paste0(stub, "growthrate_HHS", hhs, tag,".jpg"),
                          width = 1500,
                          height = 1500,
                          pointsize = 40)

    # make enough room on the right side of the plot for a secondary axis.
    orpar <- par()
    par(mar = c(5.1, 4.1, 4.1, 4.1))

    # plot "nowcast" growth rates by variant
    plot(x = 100 * hhs.summary$p_i,
         y = gr,
         log = "x",
         type = "n",
         ylim = c(min(range(gr + 1.96 * se.gr, gr - 1.96 * se.gr)[1], -75),
                  max(range(gr + 1.96 * se.gr, gr - 1.96 * se.gr)[2], 100)),
         xaxt = "n",
         xlim = c(0.0005,110),
         xlab = "Weighted share (%)",
         ylab = "Week over week growth rate (%)",
         main = paste("HHS Region", hhs))

    # Add an explanation for the lack of confidence intervals in the event that
    # the Hessian was non-invertible
    if(is.null(svymlm_hhs$SE)){
      mtext(text = "*No SE estimates b/c of non-invertible Hessian in multinomial model fit.",
            side = 3,
            line = 0,
            cex = 0.75,
            font = 4,
            col = 'red' )
    }

    # add an x-axis
    axis(side = 1,
         at     = c(0.001, 0.01, 0.1, 1, 10, 100),
         labels = c(0.001, 0.01, 0.1, 1, 10, 100))

    # add a horizontal line at 0
    abline(h = 0,
           col = "grey65")

    # add lines for each variant
    for (vv in seq(model_vars)) {
      # vertical lines for uncertainty in growth rate
      lines(x = 100 * rep(hhs.summary$p_i[vv], 2),
            y = gr[vv] + 1.96 * c(1,-1) * se.gr[vv],
            col = "blue",
            lwd = 2)

      # horizontal lines for uncertainty in variant proportions
      lines(x = pmax(0.0001, 100 * hhs.summary$p_i[vv] + 1.96 * c(1,-1) * hhs.summary$se.p_i[vv]),
            y = rep(gr[vv], 2),
            col = "blue",
            lwd = 2)
    }


    # identify variants with share > 1% or growth rate > 0%
    labels = unique(c(which(100 * hhs.summary$p_i>1),
                      which(gr > 1)))

    # label variants
    text(x = 100 * hhs.summary$p_i,
         y = gr,
         labels = c(model_vars,""),
         cex = 0.85,
         col = "grey25",
         adj = 1.15)

    # add text for doubling time
    plot_growth_rates <- axTicks(2)
    plot_doubling_times <- (log(2) / log((100 + plot_growth_rates)/100)) * 7

    # Add second axis
    axis(side = 4,
         at = plot_growth_rates,
         labels = round(plot_doubling_times,1))
    # Add second axis label
    mtext("Doubling time (days)",
          side = 4,
          line = 3)
    if (fig_gen_run) dev.off()
  } # end loop over HHS regions

  # save growth rates to file
  write.csv(x = gr_tab,
            file = paste0(script.basename,
                          "/results/wow_growth_variant_share",
                          data_date,
                          tag,
                          "_HHS.csv"),
            row.names = FALSE)


  ## Model-smoothed estimates ("updated nowcast") ----
  # (for runs 1 & 2)

  ## Model-smoothed estimates (needs Hessian for regional multinomial model)
  # (uses confidence intervals from "svymultinom" > "se.multinom" > "svyCI")

  # This section of code is only run for Run 2, but it also produces results for
  # Run 1. Therefore, it needs to aggregate some of the variants in the Nowcast
  # model to match those that are reported for Run 1.

  ### aggregation matrix ----
  # Check to see which lineages are in model_vars
  # this returns all variants with "AY" in the name
  AY_vars = model_vars[grep("AY", model_vars)]
  # this returns all variants with BA. in the name (Omicron sublineages)
  BA_vars = model_vars[grep("BA\\.", model_vars)]

  # get the names of the lineages included in Run1
  if(custom_lineages){
    run1_lineages = c(voc1, custom_lineage_names)
  } else {
    run1_lineages = voc1
  }

  # this returns all "AY" variants to be aggregated for Run 1 results
  # (i.e. not listed in run1_lineages)
  AY_agg = AY_vars[AY_vars %notin% run1_lineages]
  # this returns all "BA" variants to be aggregated for Run 1 results
  BA_agg = BA_vars[BA_vars %notin% run1_lineages]

  # all other variants to be aggregated (used for Run 1 & Run 2 results)
  Other_agg = model_vars[model_vars %notin% voc]

  # generate a matrix that indicates which lineages to aggregate for the nowcast
  # Columns are the lineages in the nowcast model, so all the defined lineages
  #  plus the "other" lineage
  # Rows are the aggregated lineages desired
  agg_var_mat <- matrix(data = 0,
                        nrow = 3,
                        ncol = (length(model_vars)+1))
  colnames(agg_var_mat) <- c(model_vars,"Other")

  # Fill in matrix values: if lineage is to be aggregated to parent lineage in
  # given row, then value = 1, else value = 0
  agg_var_mat[1,] <- ifelse(colnames(agg_var_mat) %in% c("B.1.617.2", AY_agg),1,0)
  agg_var_mat[2,] <- ifelse(colnames(agg_var_mat) %in% c("B.1.1.529", BA_agg),1,0)
  agg_var_mat[3,] <- ifelse(colnames(agg_var_mat) %in% c(Other_agg, "Other"),1,0)
  row.names(agg_var_mat) <-c("Delta Aggregated", "Omicron Aggregated", "Other Aggregated")

  ### Fortnightly estimates -----
  #define fortnights and regions to get nowcasts for
  proj_ftnts = as.Date(tail(ftnts, 2))
  # add in 2 fortnights into the future
  proj_ftnts = sort(unique(c(proj_ftnts,
                             proj_ftnts + 14,
                             proj_ftnts + 28)))

  # create a dataframe with all regions for predictions
  dfs = expand.grid(USA_or_HHSRegion = c("USA", as.character(1:10)))

  # create an empty object to hold predicted values for each region
  proj.res = c()

  # cycle over regions
  for (rgn in dfs$USA_or_HHSRegion){
    # cycle over time period
    for (ftn in proj_ftnts) {

      # get the model fit object & geiod
      if (rgn=="USA") {
        mlm   = svymlm_us$mlm
        geoid = rgn
      } else {
        mlm   = svymlm_hhs$mlm
        geoid = as.numeric(rgn)
      }

      # get the week for the timepoint
      wk = as.numeric(as.Date(ftn, origin="1970-01-01") - week0day1) %/% 7
      # convert week to model_week
      wk = wk - model_week_min - model_week_mid

      # get the estimates (and SE) for the given place & time
      ests = se.multinom(mlm = mlm,
                         newdata_1row = data.frame(
                           model_week = wk,
                           HHS = geoid
                         ),
                         composite_variant = agg_var_mat)

      # calculate the SE of the growth rate
      se.gr = with(data = ests,
                   expr = 100 * exp(sqrt(se.b_i^2 * (1 - 2 * p_i) + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2))) - 100)

      # calculate the growth rate
      gr = with(ests,
                100 * exp(b_i - sum(p_i * b_i)) - 100)

      # add in doubling times
      se.gr_link = with(data = ests,
                        expr = sqrt(se.b_i^2 * (1 - 2 * p_i) + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2)))
      gr_link = with(data = ests,
                     expr = (b_i - sum(p_i * b_i)))
      gr_lo_link = gr_link - 1.96 * se.gr_link
      gr_hi_link = gr_link + 1.96 * se.gr_link

      gr_lo = 100 * exp(gr_lo_link) - 100
      gr_hi = 100 * exp(gr_hi_link) - 100

      # calculate doubling time
      doubling_time    = log(2)/gr_link * 7
      doubling_time_lo = log(2)/gr_lo_link * 7
      doubling_time_hi = log(2)/gr_hi_link * 7

      # format estimates into dataframe with relevant info
      ests = data.frame(
        USA_or_HHSRegion = rgn,
        Fortnight_ending = as.Date(ftn, origin="1970-01-01"),
        Variant = c(model_vars,
                    "Other",
                    row.names(ests$composite_variant$matrix)),
        Share = c(ests$p_i,
                  ests$composite_variant$p_i),
        se.Share = c(ests$se.p_i,
                     ests$composite_variant$se.p_i),
        growth_rate    = c(gr,    rep(NA, length(ests$composite_variant$p_i))),
        growth_rate_lo = c(gr_lo, rep(NA, length(ests$composite_variant$p_i))),
        growth_rate_hi = c(gr_hi, rep(NA, length(ests$composite_variant$p_i))),
        doubling_time    = c(doubling_time,    rep(NA, length(ests$composite_variant$p_i))),
        doubling_time_lo = c(doubling_time_lo, rep(NA, length(ests$composite_variant$p_i))),
        doubling_time_hi = c(doubling_time_hi, rep(NA, length(ests$composite_variant$p_i)))
      )

      # Get binomial CI from p_i and se.p_i
      binom.ci = apply(X = ests,
                       MARGIN = 1,
                       FUN = function(rr) svyCI(p = as.numeric(rr[4]),
                                                s = as.numeric(rr[5])))

      # add the CI into the estimates dataframe
      ests$Share_lo = binom.ci[1,]
      ests$Share_hi = binom.ci[2,]

      # add the estimates for this specific place & time to the results
      proj.res = rbind(proj.res,
                       ests)
    } # end loop over fortnights
  } # end loop over regions

  # select just the columns we want
  proj.res = proj.res[,c('USA_or_HHSRegion',
                         'Fortnight_ending',
                         'Variant',
                         'Share',
                         'Share_lo',
                         'Share_hi',
                         'se.Share',
                         'growth_rate',
                         'growth_rate_lo',
                         'growth_rate_hi',
                         'doubling_time',
                         'doubling_time_lo',
                         'doubling_time_hi')]

  # optionally calculate the number of infections attributable to each variant
  if (calc_confirmed_infections){
    test_filepath <- paste0(script.basename, 
                            "/data/backup_", 
                            data_date, "/", 
                            data_date, "_tests_aggregated", 
                            custom_tag, ".RDS")
    
    if (file.exists(test_filepath)){
      test_list <- readRDS(file = test_filepath)
      
      # get the fortnightly test tallies & aggregate them by fn across USA
      tests_fn_us <- test_list$tests_fortnight[,
                                               .('total_test_positives' = sum(POSITIVE, na.rm = T)), 
                                               by = 'fortnight_end'][,'HHS' := 'USA']
      # aggregate fortnightly tests by HHS region
      tests_fn_hhs <- test_list$tests_fortnight[,
                                                .('total_test_positives' = sum(POSITIVE, na.rm = T)), 
                                                by = c('fortnight_end', 'HHS')]
      
      # merge the positive test results in with the variant proportion estimates
      proj.res <- merge(
        x = proj.res, 
        y = rbind(tests_fn_us,
                  tests_fn_hhs), 
        by.x = c("USA_or_HHSRegion",
                 "Fortnight_ending"), 
        by.y = c('HHS',
                 'fortnight_end'),
        all.x = TRUE)
      
      # calculate case totals for each variant
      proj.res$cases    <- proj.res$total_test_positives * proj.res$Share
      proj.res$cases_lo <- proj.res$total_test_positives * proj.res$Share_lo
      proj.res$cases_hi <- proj.res$total_test_positives * proj.res$Share_hi
    } else {
      print(paste0('File ', 
                   test_filepath, 
                   ' not found. Not calculating number of infections attributable to each variant for fortnights.'))
    }
  }
  
  
  # Format output for the run 1 lineage list
  # only include Variants that are NOT in the list provided (to avoid duplicates)
  run_1 = proj.res[proj.res$Variant %notin% c(AY_agg,
                                              "B.1.617.2",
                                              BA_agg,
                                              'B.1.1.529',
                                              Other_agg,
                                              "Other"),]

  # change the name of "Other Aggregated" to "Other" to match other output files
  run_1[run_1$Variant == "Other Aggregated","Variant"] <- "Other"

  # save the results to file
  write.csv(x = run_1,
            file = paste0(script.basename,
                          "/results/updated_nowcast_fortnightly_",
                          data_date,
                          sub(pattern = '2', replacement = '1', x = tag),
                          ".csv"),
            row.names = FALSE)

  # Format output for the run2 lineage list
  # get all the lineages names that are *NOT* "Other Aggregated"
  drop_lin <- row.names(agg_var_mat)[row.names(agg_var_mat) %notin% "Other Aggregated"]

  # Only include variants that are NOT in the list provided
  run_2 = proj.res[proj.res$Variant %notin% c(drop_lin,
                                              Other_agg,
                                              "Other"),]

  # change the name of "Other Aggregated" to "Other" to match other output files
  run_2[run_2$Variant=="Other Aggregated","Variant"] <- "Other"

  # save the results to file
  write.csv(x = run_2,
            file = paste0(script.basename,
                          "/results/updated_nowcast_fortnightly_",
                          data_date,
                          tag,
                          ".csv"),
            row.names = FALSE)


  ### Weekly estimates ----
  # same as above, but for weekly estimates instead of fortnightly

  # define and use a function to get the end-of-week dates
  cast_wks = (function(dd) as.Date(seq(from = dd[1],
                                       to = dd[2],
                                       by = 7),
                                   origin = "1970-01-01"))(range(as.Date(proj_ftnts)) + c(-7, 0))

  # optionally make predictions on a daily basis instead of weekly basis

  # cast weeks are based on fortnights, which are end-of-week dates
  cast_wks <- seq(from = min(cast_wks) - 6,
                  to   = max(cast_wks),
                  by   = 1)

  # create an empty object to hold predicted values for each region
  proj.res = c()

  # cycle over regions
  for (rgn in dfs$USA_or_HHSRegion){
    # cycle over weeks
    for (cwk in cast_wks) {

      # get the model fit & geoid
      if (rgn=="USA") {
        mlm = svymlm_us$mlm
        geoid = rgn
      } else {
        mlm = svymlm_hhs$mlm
        geoid = as.numeric(rgn)
      }

      # get the week for the given timepoint
      wk_date = as.Date(cwk, origin="1970-01-01")
      # convert date to model_week
      wk = date_to_model_week(wk_date)

      # Sunday of week
      week_start  = wk_date - as.numeric(format(wk_date, format = '%w'))
      # Saturday of week
      week_ending = week_start + 6


      # get the estimates (and SE) for the given place & time
      ests = se.multinom(mlm = mlm,
                         newdata_1row = data.frame(
                           model_week = wk,
                           HHS = geoid
                         ),
                         composite_variant = agg_var_mat)

      # calculate the SE of the growth rate
      se.gr = with(data = ests,
                   expr = 100 * exp(sqrt(se.b_i^2 * (1 - 2 * p_i) + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2))) - 100)

      # calculate the growth rate
      gr = with(ests,
                100 * exp(b_i - sum(p_i * b_i)) - 100)

      # add in doubling times
      se.gr_link = with(data = ests,
                        expr = sqrt(se.b_i^2 * (1 - 2 * p_i) + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2)))
      gr_link = with(data = ests,
                     expr = (b_i - sum(p_i * b_i)))
      gr_lo_link = gr_link - 1.96 * se.gr_link
      gr_hi_link = gr_link + 1.96 * se.gr_link

      gr_lo = 100 * exp(gr_lo_link) - 100
      gr_hi = 100 * exp(gr_hi_link) - 100

      # calculate doubling time
      doubling_time    = log(2)/gr_link * 7
      doubling_time_lo = log(2)/gr_lo_link * 7
      doubling_time_hi = log(2)/gr_hi_link * 7

      # format estimates into dataframe with relevant info
      ests = data.frame(
        USA_or_HHSRegion = rgn,
        Week_ending = week_ending, # this no longer identifies a single estimate.
        Variant = c(model_vars,
                    "Other",
                    row.names(ests$composite_variant$matrix)),
        Share = c(ests$p_i,
                  ests$composite_variant$p_i),
        se.Share = c(ests$se.p_i,
                     ests$composite_variant$se.p_i),
        growth_rate    = c(gr,    rep(NA, length(ests$composite_variant$p_i))),
        growth_rate_lo = c(gr_lo, rep(NA, length(ests$composite_variant$p_i))),
        growth_rate_hi = c(gr_hi, rep(NA, length(ests$composite_variant$p_i))),
        doubling_time    = c(doubling_time,    rep(NA, length(ests$composite_variant$p_i))),
        doubling_time_lo = c(doubling_time_lo, rep(NA, length(ests$composite_variant$p_i))),
        doubling_time_hi = c(doubling_time_hi, rep(NA, length(ests$composite_variant$p_i))),
        # add in dates
        date        = wk_date,
        week_start  = week_start,
        model_week  = wk
      )

      # Get binomial CI from p_i and se.p_i
      binom.ci = apply(X = ests,
                       MARGIN = 1 ,
                       FUN = function(rr) svyCI(p = as.numeric(rr[4]),
                                                s = as.numeric(rr[5])))

      # add the CI into the estimates dataframe
      ests$Share_lo = binom.ci[1,]
      ests$Share_hi = binom.ci[2,]

      # add the estimates for this specific place & time to the results
      proj.res = rbind(proj.res,
                       ests)
    } # end loop over weeks
  } # end loop over regions

  # select just the columns we want
  proj.res = proj.res[,c('USA_or_HHSRegion',
                         "Week_ending",
                         'Variant',
                         'Share',
                         'Share_lo',
                         'Share_hi',
                         'se.Share',
                         'growth_rate',
                         'growth_rate_lo',
                         'growth_rate_hi',
                         'doubling_time',
                         'doubling_time_lo',
                         'doubling_time_hi',
                         'date',
                         'week_start',
                         'model_week')]

  # optionally calculate the number of infections attributable to each variant
  if (calc_confirmed_infections){
    test_filepath <- paste0(script.basename, 
                            "/data/backup_", 
                            data_date, "/", 
                            data_date, "_tests_aggregated", 
                            custom_tag, ".RDS")
    
    if (file.exists(test_filepath)){
      test_list <- readRDS(file = test_filepath)
      
      # get the daily test tallies & aggregate them by fn across USA
      tests_dy_us <- test_list$tests_daily[,
                                            .('total_test_positives_daily' = sum(POSITIVE_daily, na.rm = T)), 
                                            by = 'date'][,'HHS' := 'USA']
      # aggregate daily tests by HHS region
      tests_dy_hhs <- test_list$tests_daily[,
                                             .('total_test_positives_daily' = sum(POSITIVE_daily, na.rm = T)), 
                                             by = c('date', 'HHS')]
      
      # get the weekly test tallies & aggregate them by fn across USA
      tests_wk_us <- test_list$tests_weekly[,
                                            .('total_test_positives_weekly' = sum(POSITIVE, na.rm = T)), 
                                            by = 'yr_wk'][,'HHS' := 'USA']
      # aggregate weekly tests by HHS region
      tests_wk_hhs <- test_list$tests_weekly[,
                                             .('total_test_positives_weekly' = sum(POSITIVE, na.rm = T)), 
                                             by = c('yr_wk', 'HHS')]
      
      
      # merge the daily positive test results in with the variant proportion estimates
      proj.res <- merge(
        x = proj.res, 
        y = rbind(tests_dy_us,
                  tests_dy_hhs), 
        by.x = c("USA_or_HHSRegion",
                 "date"), 
        by.y = c('HHS',
                 'date'),
        all.x = TRUE)
      
      
      # merge the positive test results in with the variant proportion estimates
      proj.res <- merge(
        x = proj.res, 
        y = rbind(tests_wk_us,
                  tests_wk_hhs)[,'Week_ending' := as.Date(yr_wk) + 6][, 'yr_wk' := NULL], 
        by.x = c("USA_or_HHSRegion",
                 "Week_ending"), 
        by.y = c('HHS',
                 'Week_ending'),
        all.x = TRUE)
      
      # calculate case totals for each variant
      proj.res$cases_daily    <- proj.res$total_test_positives_daily * proj.res$Share
      proj.res$cases_lo_daily <- proj.res$total_test_positives_daily * proj.res$Share_lo
      proj.res$cases_hi_daily <- proj.res$total_test_positives_daily * proj.res$Share_hi
      
      proj.res$cases_weekly    <- proj.res$total_test_positives_weekly * proj.res$Share
      proj.res$cases_lo_weekly <- proj.res$total_test_positives_weekly * proj.res$Share_lo
      proj.res$cases_hi_weekly <- proj.res$total_test_positives_weekly * proj.res$Share_hi
    } else {
      print(paste0('File ', 
                   test_filepath, 
                   ' not found. Not calculating number of infections attributable to each variant for weeks.'))
    }
    
  }
  
  # Format output for the run 1 lineage list
  # only include variants that are NOT in the list provided
  run_1 = proj.res[proj.res$Variant %notin% c(AY_agg,
                                              "B.1.617.2",
                                              BA_agg,
                                              'B.1.1.529',
                                              Other_agg,
                                              "Other"),]

  # change the name of "Other Aggregated" to "Other" to match other output files
  run_1[run_1$Variant == "Other Aggregated", "Variant"] <- "Other"

  # save the weekly results to file
  write.csv(x = run_1[ run_1$model_week %% 1 == 0 , !names(run_1) %in% c('total_test_positives_daily', 'cases_daily', 'cases_lo_daily', 'cases_hi_daily')],
            file = paste0(script.basename,
                          "/results/updated_nowcast_weekly_",
                          data_date,
                          sub(pattern = '2', replacement = '1', x = tag),
                          ".csv"),
            row.names = FALSE)
  # save the daily results to file
  write.csv(x = run_1[,!names(run_1) %in% c('total_test_positives_weekly', 'cases_weekly', 'cases_lo_weekly', 'cases_hi_weekly')],
            file = paste0(script.basename,
                          "/results/updated_nowcast_weekly_",
                          data_date,
                          sub(pattern = '2', replacement = '1', x = tag),
                          "_daily.csv"),
            row.names = FALSE)

  #Format output for the run2 lineage list
  # get all the lineage names that are *NOT* "Other Aggregated"
  drop_lin <- row.names(agg_var_mat)[row.names(agg_var_mat) %notin% "Other Aggregated"]

  # Only include variants that are NOT in the list provided
  run_2 = proj.res[proj.res$Variant %notin% c(drop_lin,
                                              Other_agg,
                                              "Other"),]

  # change the name of "Other Aggregated" to "Other" to match other output files
  run_2[run_2$Variant=="Other Aggregated", "Variant"] <- "Other"

  # save the weekly results to file
  write.csv(x = run_2[run_2$model_week %% 1 == 0 ,!names(run_2) %in% c('total_test_positives_daily', 'cases_daily', 'cases_lo_daily', 'cases_hi_daily')],
            file = paste0(script.basename,
                          "/results/updated_nowcast_weekly_",
                          data_date,
                          tag,
                          ".csv"),
            row.names = FALSE)
  # save the daily results to file
  write.csv(x = run_2[,!names(run_2) %in% c('total_test_positives_weekly', 'cases_weekly', 'cases_lo_weekly', 'cases_hi_weekly')],
            file = paste0(script.basename,
                          "/results/updated_nowcast_weekly_",
                          data_date,
                          tag,
                          "_daily.csv"),
            row.names = FALSE)
} # end Run2



# Run3 - State-level estimates; Rolling 4 wk bins ------------------------------
# saves a single file:
# - "/results/state_weighted_roll4wk_",ci.type,"CI_svyNEW_",data_date,tag,".csv"
if ( grepl("Run3", tag) ){

  # get the week number that corresponds to each date defined in state_time_end
  data_weeks <- as.numeric(state_time_end - week0day1) %/% 7

  # all combinations of states, 4-wk periods, & variants
  # (reverse column order for convenience)
  all.state = expand.grid(Variant = voc,
                          Roll_4wk_end = data_weeks,
                          State = sort(unique(src.dat$STUSAB)))[, 3:1]

  # calculate estimated proportions (and CI) using survey design
  ests = apply(X = all.state,
               MARGIN = 1,
               FUN = function(rr) myciprop(voc = rr[3],
                                           geoid = rr[1],
                                           svy = subset(svyDES,
                                                        week >= (as.numeric(rr[2])- 3) &
                                                          week < (as.numeric(rr[2])+1)),
                                           str = FALSE))

  # add together the combinations of states, time period, & variants with estimates
  all.state = cbind(all.state,
                    Share    = ests[1,],
                    Share_lo = ests[2,],
                    Share_hi = ests[3,],
                    DF       = ests[4,],
                    eff.size = ests[5,],
                    cv.mean  = ests[6,])

  # all combinations of states, 4-wk periods for "Other" variants
  # (column order reversed for convenience)
  others = expand.grid(Variant = "Other",
                       Roll_4wk_end = data_weeks,
                       State = sort(unique(src.dat$STUSAB)))[, 3:1]

  # calculate the estimated proportions (and CI) using survey design
  ests.others = apply(X = others,
                      MARGIN = 1,
                      FUN = function(rr) myciprop(voc = voc,
                                                  geoid = rr[1],
                                                  svy = subset(svyDES,
                                                               week >= (as.numeric(rr[2])- 3) &
                                                                 week < (as.numeric(rr[2])+1)),
                                                  str = FALSE))

  # add together the combinations of states, time period, & variants with estimates
  others = cbind(others,
                 Share    = 1-ests.others[1,],
                 Share_lo = 1-ests.others[3,],
                 Share_hi = 1-ests.others[2,],
                 DF       = ests.others[4,],
                 eff.size = ests.others[5,],
                 cv.mean  = ests.others[6,])

  # combine estimates for individual variants with estimates for "Other" variants
  all.state = rbind(all.state,
                    others)

  # empty object to hold sequence counts
  all.state.out <-c()

  # generate sequence counts by lineage, location, & date
  for(i in 1:length(data_weeks)){

    # subset the data to the relevant time period
    dat2 <- src.dat[src.dat$week >= (as.numeric(data_weeks[i]) - 3) &
                      src.dat$week < (as.numeric(data_weeks[i]) + 1),]

    # raw sample counts by variant & state
    raw_counts_state <- aggregate(count ~ VARIANT2 + STUSAB,
                                  data = dat2,
                                  FUN  = sum,
                                  drop = FALSE)

    # set NA counts to 0
    raw_counts_state$count <- ifelse(test = is.na(raw_counts_state$count) == T,
                                     yes = 0,
                                     no = raw_counts_state$count)


    # merge sequence counts with weighted proportions estimates
    all.state2 <- merge(x = all.state[all.state$Roll_4wk_end == data_weeks[i],],
                        y = raw_counts_state,
                        by.x = c("State",
                                 "Variant"),
                        by.y = c("STUSAB",
                                 "VARIANT2"),
                        all = T)

    # set NA counts to 0 again
    all.state2[is.na(all.state2$count) == T, "count"] <- 0

    # calculate denominator counts
    dss <- aggregate(count ~ State,
                     data = all.state2,
                     FUN = sum)

    # change the column name for the denominator counts
    names(dss)[grep("count",names(dss))] <- "denom_count"

    # add the denominator counts into the results dataframe
    all.state2 <- merge(x = all.state2,
                        y = dss,
                        all = T)

    # add a column for the rolling 4-week period
    all.state2$Roll_Fourweek_ending <- unique(as.Date(src.dat$yr_wk[src.dat$week==data_weeks[i]]) + 6)

    # add the results for the given time period onto the dataframe of results
    all.state.out <- rbind.data.frame(all.state.out,
                                      all.state2)
  } # end loop over weeks

  # set the Share 0 and CI limits to NA when the count for a lineage is 0
  all.state.out$Share = ifelse(test = all.state.out$Share != 0 &
                                 all.state.out$count == 0,
                               yes = 0,
                               no = all.state.out$Share)
  all.state.out$Share_lo = ifelse(test = is.na(all.state.out$Share_lo) == F &
                                    all.state.out$count == 0,
                                  yes = NA,
                                  no = all.state.out$Share_lo)
  all.state.out$Share_hi = ifelse(test = is.na(all.state.out$Share_hi) == F &
                                    all.state.out$count == 0,
                                  yes = NA,
                                  no = all.state.out$Share_hi)

  #calculate absolute CI width
  all.state.out$CI_width = all.state.out$Share_hi - all.state.out$Share_lo

  ## generate NCHS flags
  # flag estimates with Degrees of Freedom < 8
  all.state.out$flag_df = ifelse(all.state.out$DF<8, 1, 0)
  # flag estimates with effective size of < 30 (or NA)
  all.state.out$flag_eff.size = ifelse(test = all.state.out$eff.size < 30 |
                                         is.na(all.state.out$eff.size)==T,
                                       yes = 1,
                                       no = 0)
  # flag estimates with "denominator count" of < 30 (or NA)
  # (denominator = count of sequences of all variants in a given region and time period)
  all.state.out$flag_dss = ifelse(test = all.state.out$denom_count < 30|
                                    is.na(all.state.out$denom_count) == T,
                                  yes = 1,
                                  no = 0)
  # flag estimates with wide (absolute) confidence intervals
  all.state.out$flag_abs.ciw = ifelse(test = all.state.out$CI_width > 0.30 |
                                        is.na(all.state.out$CI_width) == T,
                                      yes = 1,
                                      no = 0)
  # flag estimates with wide (relative) confidence intervals
  all.state.out$flag_rel.ciw = ifelse(test = ((all.state.out$CI_width/all.state.out$Share)*100) > 130 |
                                        is.na((all.state.out$CI_width/all.state.out$Share)*100) == T,
                                      yes = 1,
                                      no = 0)

  # Single identifier for observations that have *any* NCHS flag
  all.state.out$nchs_flag = ifelse(test = all.state.out$flag_df == 1 |
                                     all.state.out$flag_eff.size == 1 |
                                     all.state.out$denom_count == 1 |
                                     all.state.out$flag_abs.ciw == 1 |
                                     all.state.out$flag_rel.ciw==1,
                                   yes = 1,
                                   no = 0)

  # Single identifier for observations that have any NCHS flag *other* than the degrees of freedom flag.
  all.state.out$nchs_flag_wodf = ifelse(test = all.state.out$flag_eff.size == 1 |
                                          all.state.out$denom_count == 1 |
                                          all.state.out$flag_abs.ciw == 1 |
                                          all.state.out$flag_rel.ciw == 1,
                                        yes = 1,
                                        no = 0)

  # select columns for the final results
  all.state.out = all.state.out[,c("State",
                                   "Roll_Fourweek_ending",
                                   "Variant",
                                   "Share",
                                   "Share_lo",
                                   "Share_hi",
                                   "count",
                                   "denom_count",
                                   "DF",
                                   "eff.size",
                                   "CI_width",
                                   "nchs_flag",
                                   "nchs_flag_wodf")]

  
  # optionally calculate the number of infections attributable to each variant
  if (calc_confirmed_infections){
    test_filepath <- paste0(script.basename, 
                            "/data/backup_", 
                            data_date, "/", 
                            data_date, "_tests_aggregated", 
                            custom_tag, ".RDS")
    
    if (file.exists(test_filepath)){
      test_list <- readRDS(file = test_filepath)
      
      # data_weeks <- as.numeric(state_time_end - week0day1) %/% 7
      
      aso_list <- lapply(X = state_time_end, FUN = function(ste){
        
        # get the tests for the given 4-week period
        tests_state <- test_list$tests_4weeks[ names(test_list$tests_4weeks) == ste ][[1]]
        
        data.table::setnames(x = tests_state,
                             old = c('STUSAB', 'POSITIVE'),
                             new = c('State', 'total_test_positives'))
        
        # merge the positive test results in with the variant proportion estimates
        merge(
          x = all.state.out[ all.state.out$Roll_Fourweek_ending == ste, ], 
          y = tests_state[, 'TOTAL' := NULL], 
          by = c("State"), 
          all.x = TRUE)
      })
      
      # combine the elements in the list
      all.state.out <- do.call(rbind, aso_list)
      
      # calculate case totals for each variant
      all.state.out$cases    <- all.state.out$total_test_positives * all.state.out$Share
      all.state.out$cases_lo <- all.state.out$total_test_positives * all.state.out$Share_lo
      all.state.out$cases_hi <- all.state.out$total_test_positives * all.state.out$Share_hi
    } else {
      print(paste0('File ', 
                   test_filepath, 
                   ' not found. Not calculating number of infections attributable to each variant for fortnights.'))
    }
    
  }
  
  
  # write results to file
  write.csv(x = all.state.out,
            file = paste0(script.basename,
                          "/results/state_weighted_roll4wk_",
                          ci.type,
                          "CI_svyNEW_",
                          data_date,
                          tag,
                          ".csv"),
            row.names = FALSE)
} # end Run3


# # save src.dat that's been prepped for analysis
if(save_datasets_to_file){
  saveRDS(object = src.dat,
          file = paste0(script.basename,
                        '/results/src.dat_', # save to results instead of 'data' folder
                        data_date,
                        tag,
                        '.RDS'))
}

# Cleanup ----------------------------------------------------------------------
# capture system time
tend = proc.time()

# calculate amount of time code took to run
tdiff = tend[3] - tstart[3]

# print out the runtime
sprintf(paste0("code for ",tag, " run took %f minutes"), round(tdiff/60, 1))