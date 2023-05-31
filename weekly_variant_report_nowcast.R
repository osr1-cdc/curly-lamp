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

# optionally calculate 99% confidence intervals (in addition to 95% intervals)
calc_99_CI_weighted <- FALSE # these might add a lot of time
calc_99_CI_nowcast  <- TRUE  # these should not add a lot of time

## optparse option list --------------------------------------------------------
{
  # -r run_number
  # -c custom_lineages
  # -n nextclade_pango
  # -v reduced_vocs
  # -t trim_weights
  # -s save_datasets_to_file
  # -p parallel_cores
  # -w weighted_methods
  # -b weight_type
  # -d voc2_extra_preaggregation
  # -e voc_aggregation_method

  # (get the run number from the command line)
  option_list <- list(

    # Run number
    optparse::make_option(
      opt_str = c("-r", "--run_number"),
      type    = "character",
      default = "2",
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

    # Get lineage definition from sc2_src.nextclade instead of the default sc2_src.pangolin
    optparse::make_option(
      opt_str = c("-n", "--nextclade_pango"),
      type    = "character",
      default = "F",
      help    = "Whether or not to swith to get lineage defintion from sc2_src.nextclade instead of sc2_src.pangolin (character value of T or F)",
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
      default = "T",
      help    = "Maximum weight for any individual sequence",
      metavar = "character"
    ),

    # optionally run the weighted estimates in parallel
    # setting to a number will use that many cores (make sure the qsub script reserves adequate cores)
    # options: must be a number or FALSE
    optparse::make_option(
      opt_str = c("-p", "--parallel_cores"),
      type    = "character",
      default = "12",
      help    = "Number of cores for parallel estimation of weighted proportions",
      metavar = "character"
    ),

    # choose how to calculated proportions (and CI)
    # options: weighted   =   weighted estimate with survey-design-based CI
    #          unweighted = unweighted estimate with survey-design-based CI
    #          both       = both weighted & unweighted (saved to different files)
    # (all options will also calculate unweighted proportions with binomial CI)
    optparse::make_option(
      opt_str = c("-w", "--weighted_methods"),
      type    = "character",
      default = "weighted",
      help    = '"weighted", "unweighted" (with survey design CI), or "both"',
      metavar = "character"
    ),

    # type of weights to use
    # - original        (weights used from 2021 to May 2023)
    # - updated         (uses regional positivity rate from NREVSS testing data (intead of state positivity rate from CLERS, as for "original"))
    # - population      (uses population-based weights (no testing data))
    optparse::make_option(
      opt_str = c("-b", "--weight_type"),
      type    = "character",
      default = "updated",
      help    = 'options include "original", "updated", and "population" ',
      metavar = "character"
    ),

    # whether or not to aggregate subvariants when looking for voc2 that are
    # over 0.5% in lag week -3
    # e.g. when TRUE, if variant xyz.123 meets the 0.5% criterion for inclusion in voc2,
    #      then all subvariants of xyz.123 (that are not already included in voc2) will be
    #      aggregated into xyz.123. This scenario is unlikely to occur unless a variant and
    #      its subvariants are split out in at the same time.
    optparse::make_option(
      opt_str = c("-d", "--voc2_extra_preaggregation"),
      type    = "character",
      default = "FALSE",
      help    = 'T/F, whether or not to pre-aggregate subvariants (to myriad parent variants) when looking for extra vocs to include in Nowcast modeling',
      metavar = "character"
    ),
    # voc aggregation method
    # original = the lengthy code based on abbreviated pango lineages
    # updated = shorter/faster code based on lineage_expanded
    optparse::make_option(
      opt_str = c("-e", "--voc_aggregation_method"),
      type    = "character",
      default = "updated", # "original"; "updated"/"lineage_expanded"
      help    = "Whether to use the original voc aggregation method (based on variant name) or the updated method (based on lineage_expanded)",
      metavar = "character"
    )
  )

  # parsing options list
  opts <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

  # use the specified flags to set several variables
  # convert custom_lineages flag to a logical value
  if(toupper(opts$custom_lineages) %in% c('T', 'TRUE', 'Y', 'YES')){
    if(toupper(opts$nextclade_pango) %in% c('F', 'FALSE', 'N', 'NO')){
      custom_lineages = TRUE
      custom_tag = "_custom"
    } else {
      if(toupper(opts$nextclade_pango) %in% c('T', 'TRUE', 'Y', 'YES')){
        custom_lineages = TRUE
        custom_tag = "_custom_nextcladepango"
      } else {
        errorCondition(message = paste0('nextclade_pango must be "T" or "F". Argument provide: ', opts$nextclade_pango))
      }
    }
  } else {
    if(toupper(opts$custom_lineages) %in% c('F', 'FALSE', 'N', 'NO')){
      if(toupper(opts$nextclade_pango) %in% c('F', 'FALSE', 'N', 'NO')){
        custom_lineages = FALSE
        custom_tag = ""
      } else {
        if(toupper(opts$nextclade_pango) %in% c('T', 'TRUE', 'Y', 'YES')){
          custom_lineages = FALSE
          custom_tag = "_nextcladepango"
        }
      }
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
      stop(message = paste0('reduced_vocs must be "T" or "F". Argument provide: ', opts$reduced_vocs))
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

  # whether or not to use parallel processing on weighted estimates
  if( toupper(opts$parallel_cores) %in% c('F', 'FALSE', 'N', 'NO', '1', '0')){
    use_parallel = FALSE
  } else {
    use_parallel = TRUE
    ncores = as.numeric(opts$parallel_cores)
  }

  # whether or not to aggregate subvariants to various parent levels when looking for extra voc to include in Nowcast modeling
  if( toupper(opts$voc2_extra_preaggregation) %in% c('F', 'FALSE', 'N', 'NO', '0')){
    voc2_extra_preaggregation = FALSE
  } else {
    voc2_extra_preaggregation = TRUE
  }


  # type of estimates to produce
  if( toupper(opts$weighted_methods) %in% c('BOTH')){
    weighted_methods <- c('weighted', 'unweighted')
  } else if( toupper(opts$weighted_methods) %in% c('WEIGHTED')){
    weighted_methods <- c('weighted')
  } else if( toupper(opts$weighted_methods) %in% c('UNWEIGHTED')){
    weighted_methods <- c('unweighted')
  } else {
    # if opts$weighted_methods is not a valid value, throw an error.
    stop(paste('"weighted_method" must be one of "weighted", "unweighted", "both".', opts$weighted_methods, 'is not an option.'))
  }
} # end optparse options list


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
if(date_frozen_toread != data_date){
  load(paste0(script.basename, "/data/", "svydat_", data_date, custom_tag, "_", date_frozen_toread, "_frozendata",".RData"))
} else {
  load(paste0(script.basename, "/data/svydat_", data_date, custom_tag, ".RData"))
  # load(paste0('/scicomp/groups-pure/Projects/SARS2Seq/repos/sc2_proportion_modeling', "/data/svydat_", data_date, custom_tag, ".RData"))
}

# create results dir and output folder definition
dir.create(paste0(script.basename,"/results"), showWarnings = F)
dir.create(paste0(script.basename,"/results/",results_folder), showWarnings = F)
output_folder <- paste0("/results/", results_folder)

# # filter out data that's older than we actually use
# # NOTE: this will prevent calculation of # of old sequences removed b/c of invalid lab name, invalid variant name, invalid weight
# svy.dat <- subset(svy.dat,
#                   yr_wk >= time_start)

# create a tag for the filenames to differentiate results from different runs
tag <- paste0("_",state_source,"_Run", opts$run_number, reduced_voc_tag, custom_tag)


# define some functions for automated lineage aggregation
# these will be needed for both "voc2_extra_preaggregation" and "voc_aggregation_method" == 'updated',
{
  # define a function to find a variant's nearest parent from a list of
  # options. It defines "parent" as the longest string that is a subset
  # of the child variant name. (So it only works on extended names.)
  # If no match is found, it returns "Other"
  # takes 3 arguments:
  # arg x  = (string) variant name to look up
  # arg y  = (string) vector of voc (extended)
  # arg no_match = (string) what to return if there is no match
  # returns string with the parent's name
  np = function(x, y, no_match = 'Other') {

    # only look for parents in y, so each y MUST be shorter than x
    y <- y[ nchar(y) <= nchar(x) ]

    # add a "." to all y that are shorter than x (to make sure that the parent variant IS a parent)
    # e.g. without the extra ".", the parent of BA.5.22 might be identified as BA.5.2 (shorter & a perfect subset)
    y. <- y
    y.[ nchar(y) < nchar(x) ] <- paste0(y[ nchar(y) < nchar(x) ], '.')

    # split out each character in x & y
    sx = strsplit(x, "", fixed = TRUE)
    sy = strsplit(y., "", fixed = TRUE)

    # calculate an array of match proportions (based on: https://stackoverflow.com/a/33120009/3174566)
    match_array <- array(
      # cycle over each string in X and Y
      data = mapply(function(X, Y) {
        # get the length of the shorter (parent) string
        #slen = seq_len(min(length(X), length(Y)))
        ylen <- seq_len(length(Y))

        # test if the first slen characters are the same in both strings
        # if they're not the same, calculate the proportion of characters (starting from the beginning) that are the same (i.e. how far you get through the string before the first mismatch)
        # wh <- (X[slen] == Y[slen])
        # if(all(wh)) return(1) else (which.min(wh) - 1) / length(slen)
        wh <- (X[ylen] == Y[ylen])
        if(all(wh)) return(1) else (which.min(wh) - 1) / length(ylen)
      },
      # values of X to pass to mapply function (just repeat x for each item in y)
      rep(sx, each = length(sy)),
      # values of Y to pass to mapply function
      sy),
      # dimensions for the array
      dim = c(length(x), length(y)),
      # names of the dimensions for the array
      dimnames = list(x, y)
    )
    # match_array has x strings on the x-axis and y strings on the
    # y-axis. The numbers in the array are the proportion of the 2 strings
    # that match (starting from beginning; so the proportion is how far through
    # the string you get before the first mismatch). A value of 1
    # means that the shorter string is a subset of the longer string.

    # id which matches are complete matches
    match100 <- colnames(match_array)[ match_array[1,] == 1 ]

    # if there are complete matches, get the parent
    if( length(match100) > 0 ) {

      # the longest complete match is the closest parent
      return(match100[ nchar(match100) == max(nchar(match100)) ])

      # # return the longest complete match (that is still shorter than x) as the nearest parent
      # match100_2 <- match100[ (nchar(match100) < nchar(x)) ]
      #
      # # if there are complete matches where
      # if( length(match100_2) > 0 ) return(match100_2[ nchar(match100_2) == max(nchar(match100_2)) ])
      # else return(no_match)
    } else {
      # if there are no complete matches, return "Other"
      return(no_match)
    }
  }
  # sapply(sort(unname(unlist(extra_voc_to_consider))), function(x) nearest_parent(x, vocxl))

  # convert function "np" to work on a vector of x (reuses the same y vector for each element in x)
  nearest_parent <- function(x, y, no_match = 'Other') {
    # this just uses "sapply" to look for the nearest parent of each x in y
    sapply( x, function(x_i) np(x_i,  y ) )
  }

  # # This is a shorter & easier-to-understand version of the "np" function
  # # (but it runs about 8 times slower than "np")
  # np2 <- function(x, y, no_match = 'Other'){
  #
  #    # y is a parent of x if x & y start with identical characters, then the next character in x is either "\\." or end of string "$"
  #    # all the parent variants in y
  #    all_pv <- y[unlist(lapply(y, function(p){
  #       grepl(pattern = paste0(gsub(pattern = '\\.', replacement = '\\\\.', paste0('^', p)), '((\\.)|($))'),
  #            x = x)
  #    }))]
  #
  #    # if there are complete matches, get the parent
  #    if( length(all_pv) > 0 ) {
  #
  #       # only the max-length parent variant
  #       all_pv[ nchar(all_pv) == max(nchar(all_pv))]
  #
  #    } else {
  #       # if there are no complete matches, return "Other"
  #       return(no_match)
  #    }
  # }
}




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
  } else if(pre_aggregation) {
    # do all aggregation before analysis ~line 500, only voc1 lineages (already aggregated)
    # will be included weighted and nowcast estimates. This might not work when voc1 has a very small list.
    voc = voc1
  } else {
    if( is.na(voc2_manual) ) {
      # read in the list of voc created in "variant_surveillance_system.R"
      # and add on any variants specified in voc2_additional
      voc = unique( # make sure there are no duplicates
        c(
          readRDS(file = paste0(script.basename, "/data/voc2_auto_", data_date, custom_tag, ".RDS")),
          # add in the "additional" vocs that need to be included
          voc2_additional
        ))
      # copy the 'voc2_auto' file over to the results folder
      file.copy(
        from = paste0(script.basename,
                      "/data/voc2_auto_", data_date, custom_tag, ".RDS"),
        to = paste0(script.basename,
                    output_folder, "/voc2_auto_", data_date, custom_tag, ".RDS")
      )

      # these will be needed for both "voc2_extra_preaggregation" and "voc_aggregation_method" == 'updated',
      if (voc2_extra_preaggregation | (tolower(opts$voc_aggregation_method) == 'updated')) {
        # create a table of alias & expanded lineages
        lut <- svy.dat[, .(expanded_lineage = unique(expanded_lineage)), by = 'VARIANT']

        # create a table of all samples with short & long names
        dat.voc <- svy.dat[,.(week, VARIANT, expanded_lineage)]

        # all the expanded lineages
        xl <- unique(svy.dat$expanded_lineage)
        # voc list (using expanded names)
        vocxl <- lut[ VARIANT %in% voc ][['expanded_lineage']]
        # sort by the number of "." in the name
        # (this is so that subvariants are analyzed before their parents)
        vocxl <- vocxl[order(nchar(as.character(vocxl)) - nchar( gsub("\\.", "", vocxl)), decreasing = T)]
      }

      # voc2_extra_preaggregation looks for sublineages of rare variants.
      # currently lineages can be added to voc2 if they're > 0.5% in a particular
      # timeframe, but that's only 0.5% without ANY aggregation.
      # Instead of looking through every possible aggregation among all variants in
      # the dataset, this takes the shortcut of only looking at possible aggregations
      # among sublineages of the variants that are in "voc".
      if(voc2_extra_preaggregation){

        # subset of svy.dat including only the weeks used for 1p
        # dat.voc.1p  <- svy.dat[ week > (current_week - 14) & week < (current_week - 2) , .(yr_wk, FORTNIGHT_END, VARIANT, expanded_lineage)]
        dat.voc.1p  <- svy.dat[ week > (current_week - 16) & week <= (current_week - 2) , .(yr_wk, FORTNIGHT_END, VARIANT, expanded_lineage)]
        # subset of svy.dat including only the weeks used for .5p
        # dat.voc.05p <- svy.dat[ week == (current_week - 2) , , .(yr_wk, FORTNIGHT_END, VARIANT, expanded_lineage)]
        dat.voc.05p <- svy.dat[ FORTNIGHT_END == (svy.dat[week == (current_week - 4),unique(FORTNIGHT_END)]) , .(yr_wk, FORTNIGHT_END, VARIANT, expanded_lineage)]
        # These filters need to match this query:   https://cdp-01.biotech.cdc.gov:8889/hue/editor?editor=72645

        # vector to hold the extra vocs that I may want to add to voc2
        # to speed things up, I'll only look for extra sublineages of the vocs, NOT extra sublineages within "OTHER"
        # (any individual sublineage that meets the criteria will already be included)
        extra_voc_to_consider <- vector(mode = 'list', length = length(vocxl))
        # each item in the list corresponds to one of the parent variants in voc
        # each item contains 0 or more subvariants of that voc that will be included
        # (split out into a list like this to make double-checking easier)
        names(extra_voc_to_consider) <- vocxl

        # find sublineages to aggregate
        # and check their (unweighted) proportion
        for(v in vocxl){
          # all subvariants of "v" in the data
          # identified b/c they start with identical characters, followed by "."
          sv <- grep(pattern = paste0(gsub(pattern = '\\.', replacement = '\\\\.', paste0('^',v)), '\\.'),
                     x = xl,
                     value = T)
          # identify subvariants of v that are already included in voc
          svv_in_voc <- base::intersect(vocxl, sv)

          # remove lineages from sve that are already included in voc
          sve <- setdiff(sv, vocxl)

          # sub-lineages of "v" should not include any lineages that are sub-lineages
          # of a child lineage that's included in voc
          for( cl in svv_in_voc ){
            # identify sub-lineage of (sublineage of v that's in voc)
            sv_sv <- grep(pattern = paste0(gsub(pattern = '\\.', replacement = '\\\\.', paste0('^',cl)), '\\.'),
                       x = xl,
                       value = T)
            # exclude those sublineages from sv
            sve <- setdiff(sve,
                           sv_sv)
          }
          # sve

          # only proceed if there are sub-lineages
          if(length(sve)>0){
            # count the "." in each name
            svedf <- data.frame(
              lin = sve,
              np = nchar(as.character(sve)) - nchar( gsub("\\.", "", sve)) # count the "." by removing all "." and seeing how much the total length changes
            )

            # order the subvariants by number of periods and then lineage
            svedf <- svedf[order(-svedf$np, svedf$lin, decreasing = F),]
            row.names(svedf) <- 1:nrow(svedf)

            # add columns for whether or not the subvariant meets the qualification
            # for inclusion in Nowcast model based on unweighted proportion
            svedf$over0.5 <- 0
            svedf$over1   <- 0

            # cycle through the lineages, aggregate any sub-lineages into their
            #  nearest parent, and see if they're over 0.5 or 1 in the relevant weeks
            for (r in 1:nrow(svedf)){

              # all subvariants for row r
              # identified b/c they start with identical characters, followed by "." and something else
              svr <- grep(pattern = paste0(gsub(pattern = '\\.', replacement = '\\\\.', paste0('^', svedf$lin[r])), '\\.'),
                          x = xl, # all extended lineages in the data
                          value = T)

              # exclude subvariants that are already included in some other subvariant of v
              # (svr above includes all subvariants of variant v, but we want to exclude
              #   subvariants (along with their children) of v that are themselves voc)
              #   e.g. when looking for subvariants of BA.2, if BA.2.12.1 is also in voc,
              #        then we don't want to include subvariants of BA.2.12.1 in the list
              #        of subvariants of BA.2
              {
                # exclude subvariants in svr that are subvariants of (a voc that is a subvariant of v)
                # e.g. when looking for subvariants of BA.2, if BA.2.12.1 is also in voc, then we
                #      don't want to include subvariants of BA.2.12.1 in the list of subvariants of BA.2
                # for XBB, XBB.1.18.1 should be excluded b/c of XBB.1, why is it not filtered out??
                svr <- setdiff(svr,
                               # get all the subvariants of the vocs that are subvariants of svr (they will be excluded)
                               unlist(lapply(base::intersect(vocxl, svr), # for each variant in vocxl or in svr
                                             # look for any child variants (starts the same followed by ".")
                                             function(v) grep(pattern = paste0(gsub(pattern = '\\.', replacement = '\\\\.', paste0('^', v)), '\\.'),
                                                              x = xl, value = T)))
                )

                # exclude subvariants in svr that have already been identified for inclusion in voc2
                # example: if BA.2.12.1 meets criteria for inclusion in voc2, then do NOT include BA.2.12.1
                #          in calculations of BA.2.12 frequencies
                svr <- setdiff(svr,
                               # get all the rows of svedf that are subvariants of svr & get all their subvariants (they will be excluded)
                               unlist(lapply(intersect(svedf[ (svedf$over0.5 == 1) | (svedf$over1 == 1), 'lin'], svr),
                                             function(v) grep(pattern = paste0(gsub(pattern = '\\.', replacement = '\\\\.', paste0('^', v)), '\\.'),
                                                              x = xl, value = T)))
                )

                # exclude subvariants in svr that are subvariants of (anything in extra_voc_to_consider that is a subvariant of v)
                svr <- setdiff(svr,
                               # get all the rows of svedf that are subvariants of svr & get all their subvariants (they will be excluded)
                               unlist(lapply(intersect(unlist(extra_voc_to_consider), svr),
                                             function(v) grep(pattern = paste0(gsub(pattern = '\\.', replacement = '\\\\.', paste0('^', v)), '\\.'),
                                                              x = xl, value = T)))
                )
              }

              # determine whether or not the variant meets the criteria for inclusion in voc2
              {
                # subset the data to get the number of sequences in each time frame that's in the specified lineage
                # whether or not the lineage is in "svr"
                dat.voc.1p[, 'in_SVR' := expanded_lineage %in% svr]
                # the number of observations of svr in each fortnight
                counts_by_var_and_fn <- dat.voc.1p[, .(fortnight_count = sum(in_SVR)), by = 'FORTNIGHT_END']
                # total counts in each time period
                counts_by_fn <- dat.voc.1p[, .(count = .N), by = 'FORTNIGHT_END']
                counts_by_var_and_fn <- merge(counts_by_var_and_fn, counts_by_fn, by = 'FORTNIGHT_END')
                counts_by_var_and_fn[, 'proportion' := fortnight_count / count]

                # identify variants that meet inclusion critera
                if( any(counts_by_var_and_fn$proportion > 0.01) )                              svedf[r, 'over1']   <- 1
                if( nrow(dat.voc.05p[expanded_lineage %in% svr ]) / nrow(dat.voc.05p)> 0.005)  svedf[r, "over0.5"] <- 1

                # clean up
                rm('counts_by_var_and_fn', 'counts_by_fn')
                dat.voc.1p[, 'in_SVR' := NULL]
              }
            } # end for loop over rows in svedf
            # extra VOCs to consider aggregating into their nearest parent
            # svedf[ (svedf$over0.5 == 1) | (svedf$over1 == 1), ]
            # columns are
            # - lin = lineage
            # - np = number of periods in the expanded lineage
            # - over0.5 = the lineage (plus any child variants) is over 0.5% in the 0.5% time period
            # - over1 = the lineage (plus any child variants) is over 1% in the 1% time period

            # add the variants to the list of vocs to consider adding to voc2
            # extra_voc_to_consider <- c(extra_voc_to_consider, svedf[ (svedf$over0.5 == 1) | (svedf$over1 == 1), 'lin'])
            extra_voc_to_consider[[ which(v == vocxl) ]] <- svedf[ (svedf$over0.5 == 1) | (svedf$over1 == 1), 'lin']
          } # end if statement (only proceed if there are subvariants)
        } # end loop over vocs to look for all subvariants
        # get all levels of all expanded_lineages
        # extra_voc_to_consider[sapply(extra_voc_to_consider, function(x) !is.null(x) & length(x) > 0)]

        # add the extra_voc_to_consider to voc (after converting back to abbreviated variant names)
        extra_voc_to_consider_shortname <- setdiff(
          unname(setNames(lut$VARIANT, lut$expanded_lineage)[unname(unlist(extra_voc_to_consider))]),
          voc
        )
        voc <- c(voc, extra_voc_to_consider_shortname)

        # save the extra voc's to file
        write.csv(x = data.frame('Added because of voc_extra_preaggregation' = extra_voc_to_consider_shortname ),
                  file = paste0(script.basename, output_folder, '/voc_extra_preaggregation_variants_', data_date, tag, '.csv'),
                  row.names = F)
      } # end if (voc2_extra_preaggregation)

    } else {
      voc = voc2_manual
    }
  } # end (NOT reduced_vocs) and (NOT pre_aggregation)
} # end Run2 voc definitions
if ( grepl("Run3", tag) ) {
  if (reduced_vocs){
    voc = voc3_reduced
  } else {
    voc = voc3
  }
} # end Run3 voc definitions

# optionally add on the custom lineages
if (custom_lineages == TRUE) {
  voc = unique(c(voc, custom_lineage_names))
}
# force-aggregate R346T lineages
if (force_aggregate_R346T==TRUE)  {
  R346T_sublineages <- voc[grepl('R346T', voc, ignore.case = TRUE)]
  voc <- voc[ voc %notin% R346T_sublineages ]
  voc <- c(voc, "R346T")
}

# force-aggregate omicron
# (even if BA.1 is listed in voc2, this will force all BA sublineages of Omicron
#  to have the same VARIANT name (B.1.1.529)) (unless explicitly excluded in "force_aggregate_omicron_except")
if ( force_aggregate_omicron & ('B.1.1.529' %in% voc) ){

  # omicron sub-lineages to exclude/aggregate
  omicron_sublineages <- voc[ grepl('(^B[AC-HJ-NP-WYZ]\\.[0-9])|(^C[A-G]\\.)', voc, ignore.case = TRUE) ]

  # but don't force-aggregate any sublineages in "force_aggregate_omicron_except"
  # if they are also in the vocs.
  omicron_sublineages <- omicron_sublineages[omicron_sublineages %notin% force_aggregate_omicron_except]

  # remove omicron sublineages (leaving "B.1.1.529")
  voc <- voc[ voc %notin% omicron_sublineages ]
}

# force aggregate XBB to other before analysis
# if (XBB_agg_to_other==TRUE) {
#   voc <- setdiff(voc, c('XBB'))
# }

# preaggregate certain lineages for CDT run
if(force_preaggregate_BN.1) {
  BN.1sub_in_voc <- sort(grep("(^BN\\.1\\.)", voc, perl = T, value = T))
  voc <- setdiff(voc, BN.1sub_in_voc)
}
if(force_preaggregate_XBB) {
  XBBsub_in_voc <- sort(grep("(^XBB\\.)", voc, perl = T, value = T))
  #don't force-aggregate any sublineagese in "force_aggregate_XBB_except"
  XBBsub_in_voc <- XBBsub_in_voc[XBBsub_in_voc %notin% force_aggregate_XBB_except]
  voc <- setdiff(voc, XBBsub_in_voc)
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
# force-aggregate B.1 into "Other"
if (force_aggregate_B.1 & ('B.1' %in% voc)) {
  voc <- voc[ voc %notin% 'B.1' ]
}


# Define voc_lut (VOC look-up table) after voc is finalized
if (voc2_extra_preaggregation | (tolower(opts$voc_aggregation_method) == 'updated')){

  # if these are not already defined, define them
  # (these will already be defined for only run2)
  if ( !exists('lut') ) {
    # create a table of alias & expanded lineages
    lut <- svy.dat[, .(expanded_lineage = unique(expanded_lineage)), by = 'VARIANT']

    # create a table of all samples with short & long names
    dat.voc <- svy.dat[,.(week, VARIANT, expanded_lineage)]

    # all the expanded lineages
    xl <- unique(svy.dat$expanded_lineage)
    # voc list (using expanded names)
    vocxl <- lut[ VARIANT %in% voc ][['expanded_lineage']]
    # sort by the number of "." in the name
    # (this is so that subvariants are analyzed before their parents)
    vocxl <- vocxl[order(nchar(as.character(vocxl)) - nchar( gsub("\\.", "", vocxl)), decreasing = T)]
  }

  # lineage_expanded for each voc
  voc_expanded <- lut[ VARIANT %in% voc ][['expanded_lineage']]

  # all variants in the data
  unique_vars <- na.omit(unique(svy.dat$VARIANT))
  # lineage_expanded for each of unique_vars
  unique_lineage_expanded <- unname(setNames(lut$expanded_lineage, lut$VARIANT)[ unique_vars ])

  # create a lookup table for all the variants in the data
  # (originally "lut" included all Pango lineages in the database)
  # variant = abreviated variant name
  # lineage_expanded
  # parent_lineage_expanded = lineage_expanded that this variant will be aggregated into
  # parent_variant = variant that this variant will be aggregated into
  voc_lut <- data.frame(
    'variant'                  = unique_vars,
    'lineage_expanded'         = unique_lineage_expanded,
    'parent_lineage_expanded'  = nearest_parent(unique_lineage_expanded, voc_expanded)
  )

  # add in the parent variant
  voc_lut$parent_variant <- setNames(voc_lut$variant, voc_lut$lineage_expanded)[voc_lut$parent_lineage_expanded]
  # update the row names
  row.names(voc_lut) <- 1:nrow(voc_lut)

  # save the voc aggregation look-up table to file
  write.csv(x = voc_lut,
            file = paste0(script.basename, output_folder, '/voc_aggregation_table_', data_date, tag, '.csv'),
            row.names = F)
}




############# REMOVE specified labs
if (remove_utahphl){
  # grep(pattern = 'utah', x = unique(svy.dat$SOURCE), ignore.case = T, value = T)
  svy.dat <- subset(x = svy.dat,
                    subset = SOURCE != 'UTAH PUBLIC HEALTH LABORATORY')
}
if (remove_broad){
  svy.dat <- subset(x = svy.dat,
                    subset = (SOURCE %notin% c('BROAD INSTITUTE', 'INFECTIOUS DISEASE PROGRAM, BROAD INSTITUTE OF HARVARD AND MIT') |
                                as.Date(received_date) %notin% as.Date(received_broad_dates)))
}
if (remove_Quest){
  svy.dat <- subset(x = svy.dat,
                    subset = (SOURCE %notin% c('QUEST DIAGNOSTICS INCORPORATED','Quest Diagnostics Incorporated', 'Infectious Diseases,  Quest Diagnostics', 'Quest Diagnostics') |
                                as.Date(yr_wk) < as.Date(remove_Quest_cutoff)-7 |
                                as.Date(yr_wk) > as.Date(remove_Quest_cutoff_end)+1 |
                                as.Date(received_date) >= as.Date(received_Quest_cutoff)))
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

  # index of sequences to be excluded b/c of invalid lab name
  # NOTE! This does not include labs that are explicity excluded above by "remove_utahphl" and "remove_broad"
  invalid_labname <- svy.dat$SOURCE == 'OTHER'
  # invalid_labname <- is.na(svy.dat$SOURCE)
  # index of sequences to be excluded b/c of invalid variant name
  invalid_variant <- is.na(svy.dat$VARIANT) | svy.dat$VARIANT == "None" | svy.dat$VARIANT == "Unassigned"

  # count sequences that are excluded by week
  if(file.exists(paste0(script.basename, "/data/backup_",data_date, custom_tag, "/dropped_sequence_counts_", data_date, custom_tag, "_v1.csv"))){
    # read in the counts of dropped sequences
    dropped_sequences <- as.data.table(read.csv(file = paste0(script.basename, "/data/backup_",data_date, custom_tag, "/dropped_sequence_counts_", data_date, custom_tag, "_v1.csv")))[,'week' := as.Date(week)]

    # invalid lab names
    iln_by_wk <- svy.dat[ invalid_labname, .(count = .N), by = yr_wk]
    iln_by_wk[,'yr_wk' := as.Date(yr_wk)]
    # merge in the counts of invalid lab names
    if(nrow(iln_by_wk) > 0)
      dropped_sequences <- rbind(
        dropped_sequences[!(week %in% iln_by_wk$yr_wk & reason == 'n_dropped_invalid_lab_name')],
        iln_by_wk[, .('week' = yr_wk, 'reason' = 'n_dropped_invalid_lab_name', 'count' = count)]
      )
    # fill in 0's for rows that didn't have any sequences dropped
    dropped_sequences[is.na(count) & reason == 'n_dropped_invalid_lab_name', 'count' := 0]



    # invalid variant names
    iv_by_wk <- svy.dat[ invalid_variant, .(count = .N), by = yr_wk]
    iv_by_wk[,'yr_wk' := as.Date(yr_wk)]
    # merge in the counts of invalid variant name
    if(nrow(iv_by_wk) > 0)
      dropped_sequences <- rbind(
        dropped_sequences[!(week %in% iv_by_wk$yr_wk & reason == 'n_dropped_invalid_variant_name')],
        iv_by_wk[, .('week' = yr_wk, 'reason' = 'n_dropped_invalid_variant_name', 'count' = count)]
      )
    # fill in 0's for rows that didn't have any sequences dropped
    dropped_sequences[is.na(count) & reason == 'n_dropped_invalid_variant_name', 'count' := 0]
  }
  # only include samples where both the lab and the variant are defined
  src.dat = subset(x = svy.dat,
                   !invalid_labname & # SOURCE != "OTHER"
                     !invalid_variant & # !is.na(VARIANT) & VARIANT != "None"
                     yr_wk >= time_start) # filter out old sequences to speed everything up
} else {
  # only include samples from these labs (i.e. no state-tagged data)
  # it would be better to update this to use svy.dat$source_type == 'Contractor' if we want to exclude tagged data again. (because the contractors have changed over time and hopefully source_type accounts for that?)
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
                     VARIANT != "None" &
                     VARIANT != "Unassigned")
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
# Note that SGTF weights are calculated for every contractor & week, but they should ONLY be applied to sequences that were identified as being oversampled! (e.g. join on c('STUSAB', 'yr_wk', 'contractor_targeted_sequencing') )



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
# calculate the "original" survey weights
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
src.dat[, "SIMPLE_ADJ_WT0" := ifelse(is.na(SAW), SAW_ALT, SAW)]
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
  src.dat[, "SIMPLE_ADJ_WT0" := ifelse(is.na(SAW), SAW_ALT, SAW)]
  # remove "SAW" and "SAW_ALT" columns
  src.dat[, c("SAW", "SAW_ALT") := .(NULL, NULL)]
}

# population-based weighting
# (number of people represented by each sequence)
src.dat[, "population_weight" := state_population / sum(count), by = c('STUSAB', 'yr_wk')]
# src.dat[, "population_weight.hhs"   := population_reporting.HHS / sum(count), by = c('HHS',    'yr_wk')] # this uses CLERS testing data to get the regional population (total population of states within each region that reported testing data to CLERS)

# UPDATED weights using NREVSS testing data
# Based on Prabasaj's "weighty_matters.rmd".
# proxy_infections    = estimated number of infections in a given state-week
# state_population    = state population
# POSITIVE.HHS.nrevss = number of positive tests reported to NREVSS in a given HHS region in a given week
# TOTAL.HHS.nrevss    = number of total    tests reported to NREVSS in a given HHS region in a given week
# POSITIVE.HHS        = number of positive tests reported to CLERS in a given HHS region in a given week
# population_reporting.HHS.nrevss = total population of the states in a given HHS region that reported tests to NREVSS in a given week
src.dat[
  ,
  "proxy_infections" := state_population * sqrt(POSITIVE.HHS.nrevss / TOTAL.HHS.nrevss) *
    sqrt( POSITIVE.HHS / population_reporting.HHS )]
# this formula is equivalent to calculating the regional infections as: sqrt(POSITIVE.HHS.nrevss / TOTAL.HHS.nrevss * population_reporting.HHS * POSITIVE.HHS) and then splitting it up among states in the region based on population (i.e. multiplying by): (state_population / population_reporting.HHS)

# calculate the updated weight (number of infections represented by each sequence) based on "proxy_infections"
src.dat[, 'updated_weight' := proxy_infections / sum(count), by = c('STUSAB', 'yr_wk')] # this still works for states that are missing POSITIVE in a given week
# choose which weights to use
src.dat$SIMPLE_ADJ_WT <- switch(
  EXPR = tolower(opts$weight_type),
  'original'       = {src.dat$SIMPLE_ADJ_WT0},
  'population'     = {src.dat$population_weight},
  'updated'        = {src.dat$updated_weight}
)


# check for weights of 0 (possible if there was very little testing in a state and no tests were positive)
# If using weight trimming, these are replaced with the minimum weight > 0.
zero_weights <- src.dat[SIMPLE_ADJ_WT == 0,]
if(nrow(zero_weights)>0 & !trim_weights) {
  warning('Some sequences have weights of 0: ')
  print(zero_weights)
}

# check for weights of NA weights
na_weights <- src.dat[is.na(SIMPLE_ADJ_WT),]
if(nrow(na_weights)>0) {
  warning('Some sequences have weights of NA: ')
  print(na_weights)
}

# check for Infinite weights
inf_weights <- src.dat[is.infinite(SIMPLE_ADJ_WT),]
if(nrow(inf_weights)>0) {
  warning('Some sequences have infinite weights: ')
  print(inf_weights)
}



# Remove NA and INF weights
{
  # index of sequences to be excluded b/c of invalid weights
  invalid_weight <- is.na(src.dat$SIMPLE_ADJ_WT) | is.infinite(src.dat$SIMPLE_ADJ_WT)

  # count sequences excluded b/c of invalid weights
  if(exists('dropped_sequences')){
    iw_by_wk <- src.dat[ invalid_weight, .(count = .N), by = yr_wk]
    # merge in the counts of invalid weights
    if(nrow(iw_by_wk) > 0){
      iw_by_wk[,'yr_wk' := as.Date(yr_wk)]

      dropped_sequences <- rbind(
        dropped_sequences[!(week %in% iw_by_wk$yr_wk & reason == 'n_dropped_invalid_weight')],
        iw_by_wk[, .('week' = yr_wk, 'reason' = 'n_dropped_invalid_weight', 'count' = count)]
      )
    }

    # fill in 0's for rows that didn't have any sequences dropped
    dropped_sequences[week %in% as.Date(unique(src.dat$yr_wk)) & is.na(count) & reason == 'n_dropped_invalid_weight', 'count' := 0]

    # save the dropped_sequence data.table again
    write.csv(x = dropped_sequences[order(week, count, decreasing = T)],
              file = paste0(script.basename, "/data/backup_",data_date, custom_tag, "/dropped_sequence_counts_", data_date, custom_tag, ".csv"),
              row.names = F)
  }
}
# remove sequences excluded b/c of invalid weights
src.dat = subset(x = src.dat,
                 !invalid_weight)

### aggregate sublineages ------------------------------------------------------
# make sure "VARIANT" is a character (rather than factor)
# (redundant with "stringsAsFactors = FALSE")
src.dat$VARIANT <- src.dat$lineage <- as.character(src.dat$VARIANT)

if( tolower(opts$voc_aggregation_method) %in% c('updated', 'lineage_expanded') ){

  # use the voc_lut to get the aggregation variant ("parent_variant") for each variant
  src.dat$VARIANT <- setNames(voc_lut$parent_variant, voc_lut$variant)[ src.dat$VARIANT ]

} else {
  # use the original variant aggregation code
  # IMPORTANT! LAST UPDATED on 2023-05-31 for pangolin-data v1.19, pangolin v4.2, usher v0.6.2
  # This version would not be suitable for usage unless manually updated

  #Identify all the clades/lineages to aggregate in the surveillance dataset
  # all the AY variants
  AY = sort(unique(src.dat$VARIANT)[grep("^AY",unique(src.dat$VARIANT), perl = T)])
  # just the AY variants that are to be aggregated (i.e. not listed in "voc")
  AY = AY[which(AY %notin% voc)]
  P1=sort(unique(src.dat$VARIANT)[grep("^P\\.1",unique(src.dat$VARIANT), perl = T)])
  P1=P1[which(P1 %notin% voc)] #vector of the P1s to aggregate
  Q=sort(unique(src.dat$VARIANT)[grep("^Q\\.",unique(src.dat$VARIANT), perl = T)])
  Q=Q[which(Q %notin% voc)] #vector of the Qs to aggregate
  B351=sort(unique(src.dat$VARIANT)[grep("^B\\.1\\.351",unique(src.dat$VARIANT), perl = T)])
  B351=B351[which(B351 %notin% voc)] #vector of the B351s to aggregate
  B621=sort(unique(src.dat$VARIANT)[grep("^B\\.1\\.621|^BB\\.",unique(src.dat$VARIANT), perl = T)])
  B621=B621[which(B621 %notin% voc)] #vector of the B621s to aggregate
  B429=sort(unique(src.dat$VARIANT)[grep("^B\\.1\\.429",unique(src.dat$VARIANT), perl = T)])

  # XBB recombinant lineages
  if('XBB.1.5.1' %in% voc) XBB.1.5.1 <- sort(grep("^XBB\\.1\\.5\\.1(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.1 <- NULL
  if('XBB.1.5.2' %in% voc) XBB.1.5.2 <- sort(grep("^XBB\\.1\\.5\\.2(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.2 <- NULL
  if('XBB.1.5.4' %in% voc) XBB.1.5.4 <- sort(grep("^XBB\\.1\\.5\\.4(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.4 <- NULL
  if('XBB.1.5.5' %in% voc) XBB.1.5.5 <- sort(grep("^XBB\\.1\\.5\\.5(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.5 <- NULL
  if('XBB.1.5.10' %in% voc) XBB.1.5.10 <- sort(grep("^XBB\\.1\\.5\\.10(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.10 <- NULL
  if('XBB.1.5.11' %in% voc) XBB.1.5.11 <- sort(grep("^XBB\\.1\\.5\\.11(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.11 <- NULL
  if('XBB.1.5.13' %in% voc) XBB.1.5.13 <- sort(grep("^XBB\\.1\\.5\\.13(?![0-9])|(^EK\\.)", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.13 <- NULL
  if('FD.2' %in% voc) FD.2 <- sort(grep("^FD\\.2(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else FD.2 <- NULL
  if('XBB.1.5.15' %in% voc) {
    XBB.1.5.15 <- sort(grep("^XBB\\.1\\.5\\.15(?![0-9])|(^FD\\.)", unique(src.dat$VARIANT), perl = T, value = T))
    XBB.1.5.15 <- setdiff(XBB.1.5.15, c(FD.2))
  } else XBB.1.5.15 <- NULL
  if('XBB.1.5.16' %in% voc) XBB.1.5.16 <- sort(grep("^XBB\\.1\\.5\\.16(?![0-9])|(^FG\\.)", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.16 <- NULL
  if('XBB.1.5.17' %in% voc) XBB.1.5.17 <- sort(grep("^XBB\\.1\\.5\\.17(?![0-9])|(^FH\\.)", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.17 <- NULL
  if('XBB.1.5.19' %in% voc) XBB.1.5.19 <- sort(grep("^XBB\\.1\\.5\\.19(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.19 <- NULL
  if('XBB.1.5.20' %in% voc) XBB.1.5.20 <- sort(grep("^XBB\\.1\\.5\\.20(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.20 <- NULL
  if('XBB.1.5.21' %in% voc) XBB.1.5.21 <- sort(grep("^XBB\\.1\\.5\\.21(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.21 <- NULL
  if('XBB.1.5.30' %in% voc) XBB.1.5.30 <- sort(grep("^XBB\\.1\\.5\\.30(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.30 <- NULL
  if('XBB.1.5.31' %in% voc) XBB.1.5.31 <- sort(grep("^XBB\\.1\\.5\\.31(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.31 <- NULL
  if('XBB.1.5.32' %in% voc) XBB.1.5.32 <- sort(grep("^XBB\\.1\\.5\\.32(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.32 <- NULL
  if('XBB.1.5.33' %in% voc) XBB.1.5.33 <- sort(grep("^XBB\\.1\\.5\\.33(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.33 <- NULL
  if('XBB.1.5.35' %in% voc) XBB.1.5.35 <- sort(grep("^XBB\\.1\\.5\\.35(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.5.35 <- NULL
  if('XBB.1.5' %in% voc){
    XBB.1.5 <- sort(grep("^XBB\\.1\\.5(?![0-9])|^E[KLMU]|^F[DGH]\\.", unique(src.dat$VARIANT), perl = T, value = T))
    XBB.1.5 <- setdiff(XBB.1.5, c(XBB.1.5.1, XBB.1.5.2, XBB.1.5.4, XBB.1.5.5, XBB.1.5.10, XBB.1.5.11, XBB.1.5.13, XBB.1.5.15, FD.2, XBB.1.5.16, XBB.1.5.17,
                                  XBB.1.5.19, XBB.1.5.20, XBB.1.5.21, XBB.1.5.30, XBB.1.5.31, XBB.1.5.32, XBB.1.5.33, XBB.1.5.35))
  } else XBB.1.5 <- NULL
  if('XBB.1.9.1' %in% voc) XBB.1.9.1 <- sort(grep("^XBB\\.1\\.9\\.1(?![0-9])|(^FL\\.)", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.9.1 <- NULL
  if('EG.1' %in% voc) EG.1 <- sort(grep("^EG\\.1(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else EG.1 <- NULL
  if('XBB.1.9.2' %in% voc) {
    XBB.1.9.2 <- sort(grep("^XBB\\.1\\.9\\.2(?![0-9])|(^EG\\.)", unique(src.dat$VARIANT), perl = T, value = T))
    XBB.1.9.2 <- setdiff(XBB.1.9.2, EG.1)
  } else XBB.1.9.2 <- NULL
  if('XBB.1.15' %in% voc) XBB.1.15 <- sort(grep("^XBB\\.1\\.15(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.15 <- NULL
  if('XBB.1.16.1' %in% voc) XBB.1.16.1 <- sort(grep("^XBB\\.1\\.16\\.1(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.1.16.1 <- NULL
  if('XBB.1.16' %in% voc) {
    XBB.1.16 <- sort(grep("^XBB\\.1\\.16(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T))
    XBB.1.16 <- setdiff(XBB.1.16, XBB.1.16.1)
  } else XBB.1.16 <- NULL
  if('FE.1' %in% voc) FE.1 <- sort(grep("^FE\\.1(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else FE.1 <- NULL
  if('XBB.1' %in% voc){
    XBB.1 <- sort(grep("^XBB\\.1(?![0-9])|^E[GKLMU]|^F[DEGHL]\\.", unique(src.dat$VARIANT), perl = T, value = T))
    XBB.1 <- setdiff(XBB.1, c(XBB.1.5.1, XBB.1.5.2, XBB.1.5.4, XBB.1.5.5, XBB.1.5.10, XBB.1.5.11, XBB.1.5.13, XBB.1.5.15, FD.2, XBB.1.5.16, XBB.1.5.17,
                              XBB.1.5.19, XBB.1.5.20, XBB.1.5.21, XBB.1.5.30, XBB.1.5.31, XBB.1.5.32, XBB.1.5.33, XBB.1.5.35, XBB.1.9.1, XBB.1.9.2, EG.1,
                              XBB.1.15, XBB.1.16, XBB.1.16.1, FE.1))
  } else XBB.1 <- NULL
  if('XBB.2.3' %in% voc) XBB.2.3 <- sort(grep("^XBB\\.2\\.3(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else XBB.2.3 <- NULL
  if('XBB.2' %in% voc) {
    XBB.2 <- sort(grep("^XBB\\.2(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T))
    XBB.2 <- setdiff(XBB.2, XBB.2.3)
  } else XBB.2 <- NULL
  if('XBB' %in% voc) {
    XBB <- sort(grep("(^XBB\\.)|^E[GKLMU]|^F[DEGHL]\\.", unique(src.dat$VARIANT), perl = T, value = T))
    XBB <- setdiff(XBB, c(XBB.1, XBB.1.5.1, XBB.1.5.2, XBB.1.5.4, XBB.1.5.5, XBB.1.5.10, XBB.1.5.11, XBB.1.5.13, XBB.1.5.15, FD.2, XBB.1.5.16, XBB.1.5.17,
                          XBB.1.5.19, XBB.1.5.20, XBB.1.5.21, XBB.1.5.30, XBB.1.5.31, XBB.1.5.32, XBB.1.5.33, XBB.1.5.35, XBB.1.9.1, XBB.1.9.2, EG.1, XBB.1.15, XBB.1.16, XBB.1.16.1,FE.1,
                          XBB.2, XBB.2.3))
  } else XBB <- NULL

  B429=B429[which(B429 %notin% voc)] #vector of the B429s to aggregate
  # omicrons [including lots of sublineages]
  # aggregate the BA sublineages [NOTE! this *WILL* aggregate all BA.1.1.x sub-sublineages into BA.1.1 *even* if BA.1.1.x is listed in voc.]
  if('BA.1.1' %in% voc) B529.BA1.1 <- sort(grep("^BC\\.|(^BA\\.1\\.1)(?![0-9])",  unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA1.1 <- NULL
  if('BA.1.15' %in% voc) B529.BA1.15 <- sort(grep("(^BA\\.1\\.15)(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA1.15 <- NULL
  if('BA.1' %in% voc){
    B529.BA1 <- sort(grep("^BC\\.|^BD\\.|(^BA\\.1)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T))
    # BA.1 subvariants do not include subvariants already included in B529.BA1.1, B529.BA1.15
    B529.BA1 <- setdiff(B529.BA1, c(B529.BA1.1, B529.BA1.15))
  } else B529.BA1 <- NULL
  if('BA.2.9' %in% voc) B529.BA2.9 <- sort(grep("(^BA\\.2\\.9)(?![0-9])",  unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA2.9 <- NULL
  if('BA.2.10.1' %in% voc) B529.BA2.10.1 <- sort(grep("(^BA\\.2\\.10\\.1)(?![0-9])|^BJ\\.", unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA2.10.1 <- NULL
  if('BA.2.10' %in% voc) {
    B529.BA2.10 <- sort(grep("(^BA\\.2\\.10)(?![0-9])|^BJ\\.", unique(src.dat$VARIANT), perl = T, value = T))
    B529.BA2.10<- setdiff(B529.BA2.10, B529.BA2.10.1)
  } else B529.BA2.10 <- NULL
  if('BA.2.18' %in% voc) B529.BA2.18 <- sort(grep("(^BA\\.2\\.18)(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA2.18 <- NULL
  if('BA.2.12.1' %in% voc) B529.BA2.12.1 <- sort(grep("^BG\\.|(^BA\\.2\\.12\\.1)(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA2.12.1 <- NULL
  if('BA.2.3.20' %in% voc) B529.BA2.3.20 <- sort(grep("(^BA\\.2\\.3\\.20)(?![0-9])|(^CM\\.)", unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA2.3.20 <- NULL
  if('BA.2.3' %in% voc) {
    B529.BA2.3 <- sort(grep("(^BA\\.2\\.3)(?![0-9])|(^DD\\.)|(^CM\\.)",  unique(src.dat$VARIANT), perl = T, value = T))
    B529.BA2.3 <- setdiff(B529.BA2.3, B529.BA2.3.20)
  } else B529.BA2.3 <- NULL
  if('BA.2.75.2' %in% voc) B529.BA2.75.2 <- sort(grep("(^BA\\.2\\.75\\.2)(?![0-9])|(^CA\\.)", unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA2.75.2 <- NULL
  if('CH.1.1.1' %in% voc) B529.CH.1.1.1 <- sort(grep("(^CH\\.1\\.1\\.1)(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else B529.CH.1.1.1<- NULL
  if('CH.1.1' %in% voc) {
    B529.CH.1.1 <- sort(grep("(^CH\\.1\\.1)(?![0-9])|(^DV\\.)|(^F[JK]\\.)", unique(src.dat$VARIANT), perl = T, value = T))
    B529.CH.1.1 <- setdiff(B529.CH.1.1, B529.CH.1.1.1)
  } else B529.CH.1.1 <- NULL
  if('BN.1.3' %in% voc) B529.BN.1.3 <- sort(grep("(^BN\\.1\\.3)(?![0-9])|(^DS\\.)|(^EJ\\.)", unique(src.dat$VARIANT), perl = T, value = T)) else B529.BN.1.3 <- NULL
  if('BN.1.5' %in% voc) B529.BN.1.5 <- sort(grep("(^BN\\.1\\.5)(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T)) else B529.BN.1.5 <- NULL
  if('BN.1' %in% voc){
    B529.BN.1 <- sort(grep("(^BN\\.1)(?!([0-9]))|(^DS\\.)|(^EJ\\.)", unique(src.dat$VARIANT), perl = T, value = T))
    B529.BN.1 <- setdiff(B529.BN.1, c(B529.BN.1.3, B529.BN.1.5))
  } else B529.BN.1 <- NULL
  if('BA.2.75' %in% voc) {
    B529.BA2.75 <- sort(grep("(^BA\\.2\\.75)(?![0-9])|^B[LMNRY]\\.|^C[ABHJV]\\.|(^D[SV]\\.)|(^E[JP]\\.)|(^F[JK]\\.)", unique(src.dat$VARIANT), perl = T, value = T))
    B529.BA2.75 <- setdiff(B529.BA2.75, c(B529.BA2.75.2, B529.BN.1, B529.BN.1.3, B529.BN.1.5, B529.CH.1.1.1, B529.CH.1.1))
  }else B529.BA2.75 <- NULL
  if('BA.2.12' %in% voc) {
    B529.BA2.12 <- sort(grep("^BG\\.|(^BA\\.2\\.12)(?![0-9])", unique(src.dat$VARIANT), perl = T, value = T))
    B529.BA2.12 <- setdiff(B529.BA2.12, c(B529.BA2.12.1))
  } else B529.BA2.12 <- NULL
  if('BA.2' %in% voc){
    B529.BA2 <- sort(grep("^B[GHJLMNPRSY]\\.|(^BA\\.2)(?![0-9])|^C[ABHJVM]\\.|(^D[DSV]\\.)|(^E[JP]\\.)|(^F[JK]\\.)", unique(src.dat$VARIANT), perl = T, value = T))
    # BA.2 subvariants do not include subvariants already included in B529.BA2.3, B529.BA2.9, B529.BA2.10, B529.BA2.12.1, B529.BA2.12
    B529.BA2 <- setdiff(B529.BA2, c(B529.BA2.3, B529.BA2.3.20, B529.BA2.9, B529.BA2.10, B529.BA2.10.1, B529.BA2.12.1, B529.BA2.12, B529.BA2.18, B529.BA2.75, B529.BA2.75.2, B529.BN.1, B529.BN.1.3, B529.BN.1.5, B529.CH.1.1.1, B529.CH.1.1))
  } else B529.BA2 <- NULL
  if('BA.3' %in% voc) B529.BA3 <- sort(grep("(^BA\\.3)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA3 <- NULL
  if('BA.4.4' %in% voc) B529.BA4.4 <- sort(grep("(^BA\\.4\\.4)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA4.4 <- NULL
  if('BA.4.6' %in% voc) B529.BA4.6 <- sort(grep("(^BA\\.4\\.6)(?![0-9])|(^DC\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA4.6 <- NULL
  if('BA.4.1' %in% voc) B529.BA4.1 <- sort(grep("(^BA\\.4\\.1)(?![0-9])|(^CS\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA4.1 <- NULL
  if('BA.4' %in% voc) {
    B529.BA4 <- sort(grep("(^BA\\.4)(?![0-9])|(^CS\\.)|(^DC\\.)",unique(src.dat$VARIANT), perl = T, value = T))
    # BA.4 subvariants do not include subvariants already included in B529.BA4.1
    B529.BA4 <- setdiff(B529.BA4, c(B529.BA4.1, B529.BA4.4, B529.BA4.6))
  } else B529.BA4 <- NULL
  # BA.5 sublineages includes BE.x [NOTE! Change this if any BE sublineages are added to the VOCs]
  if('BQ.1.1.1' %in% voc) B529.BQ.1.1.1 <- sort(grep("(^BQ\\.1\\.1\\.1)(?![0-9])|(^CZ\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.1.1 <- NULL
  if('BQ.1.1.3' %in% voc) B529.BQ.1.1.3 <- sort(grep("(^BQ\\.1\\.1\\.3)(?![0-9])|(^DR\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.1.3 <- NULL
  if('BQ.1.1.4' %in% voc) B529.BQ.1.1.4 <- sort(grep("(^BQ\\.1\\.1\\.4)(?![0-9])|(^EE\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.1.4 <- NULL
  if('BQ.1.1.5' %in% voc) B529.BQ.1.1.5 <- sort(grep("(^BQ\\.1\\.1\\.5)(?![0-9])|(^DN\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.1.5 <- NULL
  if('BQ.1.1.7' %in% voc) B529.BQ.1.1.7 <- sort(grep("(^BQ\\.1\\.1\\.7)(?![0-9])|(^DK\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.1.7 <- NULL
  if('BQ.1.1.10' %in% voc) B529.BQ.1.1.10 <- sort(grep("(^BQ\\.1\\.1\\.10)(?![0-9])|(^FA\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.1.10 <- NULL
  if('BQ.1.1.13' %in% voc) B529.BQ.1.1.13 <- sort(grep("(^BQ\\.1\\.1\\.13)(?![0-9])|(^E[FY]\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.1.13 <- NULL
  if('BQ.1.1.18' %in% voc) B529.BQ.1.1.18 <- sort(grep("(^BQ\\.1\\.1\\.18)(?![0-9])|(^ED\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.1.18 <- NULL
  if('BQ.1.1.32' %in% voc) B529.BQ.1.1.32 <- sort(grep("(^BQ\\.1\\.1\\.32)(?![0-9])|(^DT\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.1.32 <- NULL
  if('BQ.1.1.41' %in% voc) B529.BQ.1.1.41 <- sort(grep("(^BQ\\.1\\.1\\.41)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.1.41 <- NULL
  if('BQ.1.1' %in% voc) {
    B529.BQ.1.1 <- sort(grep("(^BQ\\.1\\.1)(?![0-9])|(^C[WZ]\\.)|(^D[KMNPRTU]\\.)|(^E[ADEFHNRSTVWYZ]|^F[ACMN]\\.)",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BQ.1.1 <- setdiff(B529.BQ.1.1, c(B529.BQ.1.1.1, B529.BQ.1.1.3, B529.BQ.1.1.4, B529.BQ.1.1.5, B529.BQ.1.1.7, B529.BQ.1.1.10, B529.BQ.1.1.13, B529.BQ.1.1.18, B529.BQ.1.1.32, B529.BQ.1.1.41))
  } else B529.BQ.1.1 <- NULL
  if('BQ.1.10' %in% voc) B529.BQ.1.10 <- sort(grep("(^BQ\\.1\\.10)(?![0-9])|(^EC\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.10 <- NULL
  if('BQ.1.11' %in% voc) B529.BQ.1.11 <- sort(grep("(^BQ\\.1\\.11)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.11 <- NULL
  if('BQ.1.12' %in% voc) B529.BQ.1.12 <- sort(grep("(^BQ\\.1\\.12)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.12 <- NULL
  if('BQ.1.13' %in% voc) B529.BQ.1.13 <- sort(grep("(^BQ\\.1\\.13)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.13 <- NULL
  if('BQ.1.14' %in% voc) B529.BQ.1.14 <- sort(grep("(^BQ\\.1\\.14)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.14 <- NULL
  if('BQ.1.19' %in% voc) B529.BQ.1.19 <- sort(grep("(^BQ\\.1\\.19)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.19 <- NULL
  if('BQ.1.2' %in% voc) B529.BQ.1.2 <- sort(grep("(^BQ\\.1\\.2)(?![0-9])|(^FB\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.2 <- NULL
  if('BQ.1.3' %in% voc) B529.BQ.1.3 <- sort(grep("(^BQ\\.1\\.3)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.3 <- NULL
  if('BQ.1.5' %in% voc) B529.BQ.1.5 <- sort(grep("(^BQ\\.1\\.5)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.5 <- NULL
  if('BQ.1.22' %in% voc) B529.BQ.1.22 <- sort(grep("(^BQ\\.1\\.22)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.22 <- NULL
  if('BQ.1.23' %in% voc) B529.BQ.1.23 <- sort(grep("(^BQ\\.1\\.23)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.23 <- NULL
  if('BQ.1.25.1' %in% voc) B529.BQ.1.25.1 <- sort(grep("(^BQ\\.1\\.25\\.1)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.25.1 <- NULL
  if('BQ.1.25' %in% voc) {
    B529.BQ.1.25 <- sort(grep("(^BQ\\.1\\.25)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BQ.1.25 <- setdiff(B529.BQ.1.25, B529.BQ.1.25.1)
  } else B529.BQ.1.25 <- NULL
  if('BQ.1.28' %in% voc) B529.BQ.1.28 <- sort(grep("(^BQ\\.1\\.28)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BQ.1.28 <- NULL
  if('BQ.1' %in% voc) {
    B529.BQ.1 <- sort(grep("(^BQ\\.1)(?![0-9])|(^C[WZ]\\.)|(^D[KMNPRTU]\\.)|(^E[ACDEFHNRSTVWYZ]\\.)|(^F[ABCFMN]\\.)",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BQ.1 <- setdiff(B529.BQ.1, c(B529.BQ.1.1, B529.BQ.1.1.1, B529.BQ.1.1.3, B529.BQ.1.1.4, B529.BQ.1.1.5, B529.BQ.1.1.7, B529.BQ.1.1.10, B529.BQ.1.1.13, B529.BQ.1.1.18, B529.BQ.1.1.32, B529.BQ.1.1.41,
                                      B529.BQ.1.2, B529.BQ.1.3, B529.BQ.1.5, B529.BQ.1.10, B529.BQ.1.11, B529.BQ.1.12, B529.BQ.1.13, B529.BQ.1.14, B529.BQ.1.19, B529.BQ.1.22, B529.BQ.1.23, B529.BQ.1.25, B529.BQ.1.25.1, B529.BQ.1.28))
  }else B529.BQ.1 <- NULL
  if('BE.1.1' %in% voc) {
    B529.BE.1.1 <- sort(grep("(^BE\\.1\\.1)(?![0-9])|(^BQ\\.)|(^C[CWZ]\\.)|(^D[KMNPRTU]\\.)|(^E[ACDEFHNRSTVWYZ]\\.)|(^F[ABCFMN]\\.)",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BE.1.1 <- setdiff(B529.BE.1.1, c(B529.BQ.1, B529.BQ.1.1, B529.BQ.1.1.1, B529.BQ.1.1.3, B529.BQ.1.1.4, B529.BQ.1.1.5, B529.BQ.1.1.7, B529.BQ.1.1.10, B529.BQ.1.1.13, B529.BQ.1.1.18, B529.BQ.1.1.32, B529.BQ.1.1.41,
                                          B529.BQ.1.2, B529.BQ.1.3, B529.BQ.1.5, B529.BQ.1.10, B529.BQ.1.11, B529.BQ.1.12, B529.BQ.1.13, B529.BQ.1.14, B529.BQ.1.19, B529.BQ.1.22, B529.BQ.1.23, B529.BQ.1.25, B529.BQ.1.25.1, B529.BQ.1.28))
  }else B529.BE.1.1 <- NULL
  if('BE.1' %in% voc) {
    B529.BE.1 <- sort(grep("(^BE\\.1)(?![0-9])|(^BQ\\.)|(^C[CWZ]\\.)|(^D[KMNPRTUW]\\.)|(^E[ACDEFHNRSTVWYZ]\\.)|(^F[ABCFMN]\\.)",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BE.1 <- setdiff(B529.BE.1, c(B529.BE.1.1, B529.BQ.1, B529.BQ.1.1, B529.BQ.1.1.1, B529.BQ.1.1.3, B529.BQ.1.1.4, B529.BQ.1.1.5, B529.BQ.1.1.7, B529.BQ.1.1.10, B529.BQ.1.1.13, B529.BQ.1.1.18, B529.BQ.1.1.32, B529.BQ.1.1.41,
                                      B529.BQ.1.2, B529.BQ.1.3, B529.BQ.1.5, B529.BQ.1.10, B529.BQ.1.11, B529.BQ.1.12, B529.BQ.1.13, B529.BQ.1.14, B529.BQ.1.19, B529.BQ.1.22, B529.BQ.1.23, B529.BQ.1.25, B529.BQ.1.25.1, B529.BQ.1.28))
  } else B529.BE.1 <- NULL
  if('BE.3' %in% voc) B529.BE.3 <- sort(grep("(^BE\\.3)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BE.3 <- NULL
  if('CQ.2' %in% voc) B529.CQ.2 <- sort(grep("(^CQ\\.2)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.CQ.2 <- NULL
  if('BA.5.3.1' %in% voc) {
    B529.BA5.3.1 <- sort(grep("(^BA\\.5\\.3\\.1)(?![0-9])|(^B[EQ]\\.)|(^C[CWZQ]\\.)|(^D[KMNPRTUW]\\.)|(^E[ACDEFHNRSTVWYZ]\\.)|(^F[ABCFMN]\\.)",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BA5.3.1 <- setdiff(B529.BA5.3.1, c(B529.BE.1, B529.BE.1.1, B529.BE.3, B529.BQ.1, B529.BQ.1.1, B529.BQ.1.1.1, B529.BQ.1.1.3, B529.BQ.1.1.4, B529.BQ.1.1.5, B529.BQ.1.1.7, B529.BQ.1.1.10, B529.BQ.1.1.13, B529.BQ.1.1.18, B529.BQ.1.1.32, B529.BQ.1.1.41,
                                            B529.BQ.1.2, B529.BQ.1.3, B529.BQ.1.5, B529.BQ.1.10, B529.BQ.1.11, B529.BQ.1.12, B529.BQ.1.13, B529.BQ.1.14, B529.BQ.1.19, B529.BQ.1.22, B529.BQ.1.23, B529.BQ.1.25, B529.BQ.1.25.1, B529.BQ.1.28, B529.CQ.2))
  } else B529.BA5.3.1 <- NULL
  if('BF.5' %in% voc) B529.BF.5 <- sort(grep("(^BF\\.5)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BF.5 <- NULL
  if('BF.7.4.1' %in% voc) B529.BF.7.4.1 <- sort(grep("(^BF\\.7\\.4\\.1)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BF.7.4.1 <- NULL
  if('BF.7.4' %in% voc) {
    B529.BF.7.4 <- sort(grep("(^BF\\.7\\.4)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BF.7.4 <- setdiff(B529.BF.7.4, B529.BF.7.4.1)
  } else B529.BF.7.4 <- NULL
  if('BF.7' %in% voc) {
    B529.BF.7 <- sort(grep("(^BF\\.7)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BF.7 <- setdiff(B529.BF.7, c(B529.BF.7.4.1, B529.BF.7.4))
  } else B529.BF.7 <- NULL
  if('BF.8' %in% voc) B529.BF.8<- sort(grep("(^BF\\.8)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BF.8 <- NULL
  if('BF.10' %in% voc) B529.BF.10 <- sort(grep("(^BF\\.10)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BF.10 <- NULL
  if('BF.11' %in% voc) B529.BF.11 <- sort(grep("(^BF\\.11)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BF.11 <- NULL
  if('BF.13' %in% voc) B529.BF.13 <- sort(grep("(^BF\\.13)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BF.13 <- NULL
  if('BF.21' %in% voc) B529.BF.21 <- sort(grep("(^BF\\.21)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BF.21 <- NULL
  if('BF.26' %in% voc) B529.BF.26 <- sort(grep("(^BF\\.26)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BF.26 <- NULL
  if('BF.27' %in% voc) B529.BF.27 <- sort(grep("(^BF\\.27)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BF.27 <- NULL
  if('BA.5.2.1' %in% voc) {
    B529.BA5.2.1 <- sort(grep("^BF\\.|(^BA\\.5\\.2\\.1)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BA5.2.1 <- setdiff(B529.BA5.2.1, c(B529.BF.5, B529.BF.7, B529.BF.7.4.1, B529.BF.7.4, B529.BF.8, B529.BF.10, B529.BF.11, B529.BF.13, B529.BF.21, B529.BF.26, B529.BF.27))
  } else B529.BA5.2.1 <- NULL
  if('BA.5.2.6' %in% voc) B529.BA5.2.6 <- sort(grep("(^BA\\.5\\.2\\.6)(?![0-9])|(^CP\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.2.6 <- NULL
  if('BA.5.2.9' %in% voc) B529.BA5.2.9 <- sort(grep("(^BA\\.5\\.2\\.9)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.2.9 <- NULL
  if('BA.5.2.20' %in% voc) B529.BA5.2.20 <- sort(grep("(^BA\\.5\\.2\\.20)(?![0-9])|(^BV\\.)",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.2.20 <- NULL
  if('BA.5.2.21' %in% voc) B529.BA5.2.21 <- sort(grep("(^BA\\.5\\.2\\.21)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.2.21 <- NULL
  if('BA.5.2.23' %in% voc) B529.BA5.2.23 <- sort(grep("(^BA\\.5\\.2\\.23)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.2.23 <- NULL
  if('BA.5.2.31' %in% voc) B529.BA5.2.31 <- sort(grep("(^BA\\.5\\.2\\.31)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.2.31 <- NULL
  if('BA.5.2.34' %in% voc) B529.BA5.2.34 <- sort(grep("(^BA\\.5\\.2\\.34)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.2.34 <- NULL
  if('CK.1' %in% voc) B529.CK.1 <- sort(grep("(^CK\\.1)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.CK.1<- NULL
  if('CR.1.1' %in% voc) B529.CR.1.1 <- sort(grep("(^CR\\.1\\.1)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.CR.1.1<- NULL
  if('BA.5.2' %in% voc) {
    B529.BA5.2 <- sort(grep("(^BA\\.5\\.2)(?![0-9])|^B[FUVZ]\\.|^C[DEFGKNPRTY]\\.|^D[ABGQYZ]\\.",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BA5.2 <- setdiff(B529.BA5.2, c(B529.BA5.2.1, B529.BA5.2.6, B529.BA5.2.9, B529.BA5.2.20, B529.BA5.2.21, B529.BA5.2.23, B529.CK.1, B529.CR.1.1, B529.BA5.2.31, B529.BA5.2.34,
                                        B529.BF.5, B529.BF.7, B529.BF.7.4.1, B529.BF.7.4, B529.BF.8, B529.BF.10, B529.BF.11, B529.BF.13, B529.BF.21, B529.BF.26, B529.BF.27))
  } else B529.BA5.2 <- NULL
  if('BA.5.6' %in% voc) B529.BA5.6 <- sort(grep("(^BA\\.5\\.6)(?![0-9])|^BW\\.",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.6 <- NULL
  if('BA.5.5.1' %in% voc) B529.BA5.5.1 <- sort(grep("(^BA\\.5\\.5\\.1)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.5.1 <- NULL
  if('BA.5.5' %in% voc) {
    B529.BA5.5 <- sort(grep("(^BA\\.5\\.5)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BA5.5 <- setdiff(B529.BA5.5, B529.BA5.5.1)
  }else B529.BA5.5 <- NULL
  if('BA.5.1.1' %in% voc) B529.BA5.1.1 <- sort(grep("(^BA\\.5\\.1\\.1)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.1.1 <- NULL
  if('BA.5.1.10' %in% voc) B529.BA5.1.10 <- sort(grep("(^BA\\.5\\.1\\.10)(?![0-9])|^BK\\.",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.1.10 <- NULL
  if('BA.5.1.18' %in% voc) B529.BA5.1.18 <- sort(grep("(^BA\\.5\\.1\\.18)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.1.18 <- NULL
  if('BA.5.1.22' %in% voc) B529.BA5.1.22 <- sort(grep("(^BA\\.5\\.1\\.22)(?![0-9])|^DH\\.",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.1.22 <- NULL
  if('BA.5.1.23' %in% voc) B529.BA5.1.23 <- sort(grep("(^BA\\.5\\.1\\.23)(?![0-9])|^DE\\.",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.1.23 <- NULL
  if('BA.5.1.27' %in% voc) B529.BA5.1.27 <- sort(grep("(^BA\\.5\\.1\\.27)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.1.27 <- NULL
  if('BA.5.1.2' %in% voc) B529.BA5.1.2 <- sort(grep("(^BA\\.5\\.1\\.2)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.1.2 <- NULL
  if('BA.5.1.5' %in% voc) B529.BA5.1.5 <- sort(grep("(^BA\\.5\\.1\\.5)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.BA5.1.5 <- NULL
  if('BA.5.1' %in% voc) {
    B529.BA5.1 <- sort(grep("(^BA\\.5\\.1)(?![0-9])|^BA\\.5\\.1\\.|^B[KT]\\.|^C[LU]\\.|^D[HEJL]\\.|^E[BQ]\\.",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BA5.1 <- setdiff(B529.BA5.1, c(B529.BA5.1.1, B529.BA5.1.10, B529.BA5.1.2, B529.BA5.1.5, B529.BA5.1.18, B529.BA5.1.22, B529.BA5.1.23, B529.BA5.1.27))
  } else B529.BA5.1 <- NULL
  if('DF.1' %in% voc) B529.DF.1 <- sort(grep("(^DF\\.1)(?![0-9])",unique(src.dat$VARIANT), perl = T, value = T)) else B529.DF.1 <- NULL
  if('BA.5' %in% voc){
    B529.BA5 <- sort(grep("(^B[KFEQTUVWZ]\\.)|(^C[CDEFGLKNPQRTUWYZ]\\.)|(^D[ABEFGHJKLMNPQRTUWYZ])\\.|((^BA\\.5)(?![0-9]))|(^E[ABCDEFHNQRSTVWYZ]\\.)|(^F[ABCFMN]\\.)",unique(src.dat$VARIANT), perl = T, value = T))
    B529.BA5 <- setdiff(B529.BA5, c(B529.BA5.1, B529.BA5.1.1, B529.BA5.1.10, B529.BA5.1.2, B529.BA5.1.5, B529.BA5.1.18, B529.BA5.1.22, B529.BA5.1.23, B529.BA5.1.27,
                                    B529.BA5.2, B529.BA5.2.1, B529.BA5.2.6, B529.BA5.2.9, B529.BA5.2.20, B529.BA5.2.21, B529.BA5.2.23, B529.CK.1, B529.CR.1.1, B529.BA5.2.31, B529.BA5.2.34, B529.BA5.3.1, B529.BA5.5, B529.BA5.5.1, B529.BA5.6,
                                    B529.BE.1, B529.BE.1.1, B529.BE.3, B529.BF.5, B529.BF.7, B529.BF.7.4.1, B529.BF.7.4, B529.BF.8, B529.BF.10, B529.BF.11, B529.BF.13, B529.BF.21, B529.BF.26, B529.BF.27,
                                    B529.BQ.1, B529.BQ.1.1, B529.BQ.1.1.1, B529.BQ.1.1.3, B529.BQ.1.1.4, B529.BQ.1.1.5, B529.BQ.1.1.7, B529.BQ.1.1.10, B529.BQ.1.1.13, B529.BQ.1.1.18, B529.BQ.1.1.32, B529.BQ.1.1.41,
                                    B529.BQ.1.2, B529.BQ.1.3, B529.BQ.1.5, B529.BQ.1.10, B529.BQ.1.11, B529.BQ.1.12, B529.BQ.1.13, B529.BQ.1.14, B529.BQ.1.19, B529.BQ.1.22, B529.BQ.1.23, B529.BQ.1.25, B529.BQ.1.25.1, B529.BQ.1.28,
                                    B529.CQ.2, B529.DF.1))
  } else B529.BA5 <- NULL

  # safety check: make sure that no variants are in the multiple sublineage groups
  B.529.all <- c(B529.BA1, B529.BA1.1, B529.BA1.15, B529.BA2, B529.BA2.3, B529.BA2.3.20, B529.BA2.9,
                 B529.BA2.10, B529.BA2.10.1, B529.BA2.12, B529.BA2.12.1, B529.BA2.18, B529.BA2.75.2, B529.BA2.75,B529.BN.1, B529.BN.1.3, B529.BN.1.5, B529.CH.1.1.1, B529.CH.1.1, B529.BA3, B529.BA4,
                 B529.BA4.1, B529.BA4.4, B529.BA4.6, B529.BA5, B529.BA5.1, B529.BA5.1.1, B529.BA5.1.10, B529.BA5.1.18, B529.BA5.1.22, B529.BA5.1.23, B529.BA5.1.27,
                 B529.BA5.1.2, B529.BA5.1.5, B529.BA5.2, B529.BA5.2.1, B529.BA5.2.6, B529.BA5.2.9, B529.BA5.2.20, B529.BA5.2.21, B529.BA5.2.23, B529.CK.1, B529.CR.1.1, B529.BA5.2.31, B529.BA5.2.34,
                 B529.BA5.3.1, B529.BA5.5, B529.BA5.5.1, B529.BA5.6,
                 B529.BE.1, B529.BE.1.1, B529.BE.3, B529.BF.5, B529.BF.7, B529.BF.7.4.1, B529.BF.7.4, B529.BF.8, B529.BF.10, B529.BF.11, B529.BF.13, B529.BF.21, B529.BF.26, B529.BF.27,
                 B529.BQ.1, B529.BQ.1.1, B529.BQ.1.1.1, B529.BQ.1.1.3, B529.BQ.1.1.4, B529.BQ.1.1.5, B529.BQ.1.1.7, B529.BQ.1.1.10, B529.BQ.1.1.13, B529.BQ.1.1.18, B529.BQ.1.1.32, B529.BQ.1.1.41,
                 B529.BQ.1.2, B529.BQ.1.3, B529.BQ.1.5, B529.BQ.1.10, B529.BQ.1.11, B529.BQ.1.12, B529.BQ.1.13, B529.BQ.1.14, B529.BQ.1.19, B529.BQ.1.22, B529.BQ.1.23, B529.BQ.1.25, B529.BQ.1.25.1, B529.BQ.1.28,
                 B529.CQ.2, B529.DF.1)

  if(any(duplicated(B.529.all))) stop(message = paste0(B.529.all[duplicated(B.529.all)], ' appear in multiple BA sublineage groups. Check B529.BA1, B529.BA1.1, B529.BA.1.15, B529.BA2, B529.BA2.3, B529.BA2.9, B529.BA2.10, B529.BA3, B529.BA4, B529.BA5.'))
  B529=sort(grep("(^B\\.1\\.1\\.529)|(^B[AC-HJ-NP-VYZ]\\.)|(^C[A-HJ-NP-WYZ]\\.)|(^D[A-HJ-NP-WYZ]\\.)|(^E[A-FHJNPQRSTVWYZ]\\.)|(^F[ABCFJKMN]\\.)",unique(src.dat$VARIANT), value = T))
  B529=B529[ B529 %notin% c(voc, B.529.all) ] #vector of the B529s to aggregate


  # Aggregate sublineages to the parent lineage
  if(P.1_agg==TRUE)     {src.dat[src.dat$VARIANT %in% P1,  "VARIANT"]    <- "P.1"}
  if(B.1.351_agg==TRUE) {src.dat[src.dat$VARIANT %in% B351,"VARIANT"]    <- "B.1.351"}
  if(B.1.621_agg==TRUE) {src.dat[src.dat$VARIANT %in% B621,"VARIANT"]    <- 'B.1.621' }
  if(Q.1_3_agg==TRUE)   {src.dat[src.dat$VARIANT %in% Q,   "VARIANT"]    <- "B.1.1.7"}
  if(AY_agg==TRUE)      {src.dat[src.dat$VARIANT %in% AY,  "VARIANT"]    <- "B.1.617.2"}
  if(B429_7_agg==TRUE)  {src.dat[src.dat$VARIANT %in% B429,"VARIANT"]    <- "B.1.427"}
  if(B.1.1.529_agg==TRUE)  {
    src.dat[src.dat$VARIANT %in% B529[B529 %notin% voc],            "VARIANT"] <- "B.1.1.529"
    src.dat[src.dat$VARIANT %in% B529.BA1[   B529.BA1    %notin% voc],"VARIANT"] <- "BA.1"
    src.dat[src.dat$VARIANT %in% B529.BA1.1[ B529.BA1.1  %notin% voc],"VARIANT"] <- "BA.1.1"
    src.dat[src.dat$VARIANT %in% B529.BA1.15[B529.BA1.15 %notin% voc],"VARIANT"] <- "BA.1.15"
    src.dat[src.dat$VARIANT %in% B529.BA2[   B529.BA2    %notin% voc],"VARIANT"] <- "BA.2"
    src.dat[src.dat$VARIANT %in% B529.BA2.3[ B529.BA2.3  %notin% voc],"VARIANT"] <- "BA.2.3"
    src.dat[src.dat$VARIANT %in% B529.BA2.3.20[ B529.BA2.3.20  %notin% voc],"VARIANT"] <- "BA.2.3.20"
    src.dat[src.dat$VARIANT %in% B529.BA2.9[ B529.BA2.9  %notin% voc],"VARIANT"] <- "BA.2.9"
    src.dat[src.dat$VARIANT %in% B529.BA2.10.1[B529.BA2.10.1 %notin% voc],"VARIANT"] <- "BA.2.10.1"
    src.dat[src.dat$VARIANT %in% B529.BA2.10[B529.BA2.10 %notin% voc],"VARIANT"] <- "BA.2.10"
    src.dat[src.dat$VARIANT %in% B529.BA2.12[B529.BA2.12 %notin% voc],"VARIANT"] <- "BA.2.12"
    src.dat[src.dat$VARIANT %in% B529.BA2.75[B529.BA2.75 %notin% voc],"VARIANT"] <- "BA.2.75"
    src.dat[src.dat$VARIANT %in% B529.CH.1.1.1[B529.CH.1.1.1 %notin% voc],"VARIANT"] <- "CH.1.1.1"
    src.dat[src.dat$VARIANT %in% B529.CH.1.1[B529.CH.1.1 %notin% voc],"VARIANT"] <- "CH.1.1"
    src.dat[src.dat$VARIANT %in% B529.BA2.75.2[B529.BA2.75.2 %notin% voc],"VARIANT"] <- "BA.2.75.2"
    src.dat[src.dat$VARIANT %in% B529.BA2.12.1[B529.BA2.12.1 %notin% voc],"VARIANT"] <- "BA.2.12.1"
    src.dat[src.dat$VARIANT %in% B529.BA2.18[B529.BA2.18 %notin% voc],"VARIANT"] <- "BA.2.18"
    src.dat[src.dat$VARIANT %in% B529.BN.1[B529.BN.1 %notin% voc],"VARIANT"] <- "BN.1"
    src.dat[src.dat$VARIANT %in% B529.BN.1.3[B529.BN.1.3 %notin% voc],"VARIANT"] <- "BN.1.3"
    src.dat[src.dat$VARIANT %in% B529.BN.1.5[B529.BN.1.5 %notin% voc],"VARIANT"] <- "BN.1.5"
    src.dat[src.dat$VARIANT %in% B529.BA3[B529.BA3 %notin% voc],"VARIANT"] <- "BA.3"
    src.dat[src.dat$VARIANT %in% B529.BA4[B529.BA4 %notin% voc],"VARIANT"] <- "BA.4"
    src.dat[src.dat$VARIANT %in% B529.BA4.1[B529.BA4.1 %notin% voc],"VARIANT"] <- "BA.4.1"
    src.dat[src.dat$VARIANT %in% B529.BA4.4[B529.BA4.4 %notin% voc],"VARIANT"] <- "BA.4.4"
    src.dat[src.dat$VARIANT %in% B529.BA4.6[B529.BA4.6 %notin% voc],"VARIANT"] <- "BA.4.6"
    src.dat[src.dat$VARIANT %in% B529.BA5[B529.BA5 %notin% voc],"VARIANT"] <- "BA.5"
    src.dat[src.dat$VARIANT %in% B529.BA5.1[B529.BA5.1 %notin% voc],"VARIANT"] <- "BA.5.1"
    src.dat[src.dat$VARIANT %in% B529.BA5.1.1[B529.BA5.1.1 %notin% voc],"VARIANT"] <- "BA.5.1.1"
    src.dat[src.dat$VARIANT %in% B529.BA5.1.2[B529.BA5.1.2 %notin% voc],"VARIANT"] <- "BA.5.1.2"
    src.dat[src.dat$VARIANT %in% B529.BA5.1.5[B529.BA5.1.5 %notin% voc],"VARIANT"] <- "BA.5.1.5"
    src.dat[src.dat$VARIANT %in% B529.BA5.1.10[B529.BA5.1.10 %notin% voc],"VARIANT"] <- "BA.5.1.10"
    src.dat[src.dat$VARIANT %in% B529.BA5.1.18[B529.BA5.1.18 %notin% voc],"VARIANT"] <- "BA.5.1.18"
    src.dat[src.dat$VARIANT %in% B529.BA5.1.22[B529.BA5.1.22 %notin% voc],"VARIANT"] <- "BA.5.1.22"
    src.dat[src.dat$VARIANT %in% B529.BA5.1.23[B529.BA5.1.23 %notin% voc],"VARIANT"] <- "BA.5.1.23"
    src.dat[src.dat$VARIANT %in% B529.BA5.1.27[B529.BA5.1.27 %notin% voc],"VARIANT"] <- "BA.5.1.27"
    src.dat[src.dat$VARIANT %in% B529.BA5.2[B529.BA5.2 %notin% voc],"VARIANT"] <- "BA.5.2"
    src.dat[src.dat$VARIANT %in% B529.BA5.2.1[B529.BA5.2.1 %notin% voc],"VARIANT"] <- "BA.5.2.1"
    src.dat[src.dat$VARIANT %in% B529.BA5.2.6[B529.BA5.2.6 %notin% voc],"VARIANT"] <- "BA.5.2.6"
    src.dat[src.dat$VARIANT %in% B529.BA5.2.9[B529.BA5.2.9 %notin% voc],"VARIANT"] <- "BA.5.2.9"
    src.dat[src.dat$VARIANT %in% B529.BA5.2.20[B529.BA5.2.20 %notin% voc],"VARIANT"] <- "BA.5.2.20"
    src.dat[src.dat$VARIANT %in% B529.CR.1.1[B529.CR.1.1 %notin% voc],"VARIANT"] <- "CR.1.1"
    src.dat[src.dat$VARIANT %in% B529.BA5.2.21[B529.BA5.2.21 %notin% voc],"VARIANT"] <- "BA.5.2.21"
    src.dat[src.dat$VARIANT %in% B529.BA5.2.23[B529.BA5.2.23 %notin% voc],"VARIANT"] <- "BA.5.2.23"
    src.dat[src.dat$VARIANT %in% B529.CK.1[B529.CK.1 %notin% voc],"VARIANT"] <- "CK.1"
    src.dat[src.dat$VARIANT %in% B529.BA5.2.31[B529.BA5.2.31 %notin% voc],"VARIANT"] <- "BA.5.2.31"
    src.dat[src.dat$VARIANT %in% B529.BA5.2.34[B529.BA5.2.34 %notin% voc],"VARIANT"] <- "BA.5.2.34"
    src.dat[src.dat$VARIANT %in% B529.BA5.3.1[B529.BA5.3.1 %notin% voc],"VARIANT"] <- "BA.5.3.1"
    src.dat[src.dat$VARIANT %in% B529.BA5.5[B529.BA5.5 %notin% voc],"VARIANT"] <- "BA.5.5"
    src.dat[src.dat$VARIANT %in% B529.BA5.5.1[B529.BA5.5.1 %notin% voc],"VARIANT"] <- "BA.5.5.1"
    src.dat[src.dat$VARIANT %in% B529.BA5.6[B529.BA5.6 %notin% voc],"VARIANT"] <- "BA.5.6"
    src.dat[src.dat$VARIANT %in% B529.BE.1[B529.BE.1 %notin% voc],"VARIANT"] <- "BE.1"
    src.dat[src.dat$VARIANT %in% B529.BE.1.1[B529.BE.1.1 %notin% voc],"VARIANT"] <- "BE.1.1"
    src.dat[src.dat$VARIANT %in% B529.BE.3[B529.BE.3 %notin% voc],"VARIANT"] <- "BE.3"
    src.dat[src.dat$VARIANT %in% B529.BF.5[B529.BF.5 %notin% voc],"VARIANT"] <- "BF.5"
    src.dat[src.dat$VARIANT %in% B529.BF.7[B529.BF.7 %notin% voc],"VARIANT"] <- "BF.7"
    src.dat[src.dat$VARIANT %in% B529.BF.7.4.1[B529.BF.7.4.1 %notin% voc],"VARIANT"] <- "BF.7.4.1"
    src.dat[src.dat$VARIANT %in% B529.BF.7.4[B529.BF.7.4 %notin% voc],"VARIANT"] <- "BF.7.4"
    src.dat[src.dat$VARIANT %in% B529.BF.8[B529.BF.8 %notin% voc],"VARIANT"] <- "BF.8"
    src.dat[src.dat$VARIANT %in% B529.BF.10[B529.BF.10 %notin% voc],"VARIANT"] <- "BF.10"
    src.dat[src.dat$VARIANT %in% B529.BF.11[B529.BF.11 %notin% voc],"VARIANT"] <- "BF.11"
    src.dat[src.dat$VARIANT %in% B529.BF.13[B529.BF.13 %notin% voc],"VARIANT"] <- "BF.13"
    src.dat[src.dat$VARIANT %in% B529.BF.21[B529.BF.21 %notin% voc],"VARIANT"] <- "BF.21"
    src.dat[src.dat$VARIANT %in% B529.BF.26[B529.BF.26 %notin% voc],"VARIANT"] <- "BF.26"
    src.dat[src.dat$VARIANT %in% B529.BF.27[B529.BF.27 %notin% voc],"VARIANT"] <- "BF.27"
    src.dat[src.dat$VARIANT %in% B529.BQ.1[B529.BQ.1 %notin% voc],"VARIANT"] <- "BQ.1"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.1[B529.BQ.1.1 %notin% voc],"VARIANT"] <- "BQ.1.1"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.1.1[B529.BQ.1.1.1 %notin% voc],"VARIANT"] <- "BQ.1.1.1"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.1.3[B529.BQ.1.1.3 %notin% voc],"VARIANT"] <- "BQ.1.1.3"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.1.4[B529.BQ.1.1.4 %notin% voc],"VARIANT"] <- "BQ.1.1.4"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.1.5[B529.BQ.1.1.5 %notin% voc],"VARIANT"] <- "BQ.1.1.5"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.1.7[B529.BQ.1.1.7 %notin% voc],"VARIANT"] <- "BQ.1.1.7"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.1.10[B529.BQ.1.1.10 %notin% voc],"VARIANT"] <- "BQ.1.1.10"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.1.13[B529.BQ.1.1.13 %notin% voc],"VARIANT"] <- "BQ.1.1.13"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.1.18[B529.BQ.1.1.18 %notin% voc],"VARIANT"] <- "BQ.1.1.18"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.1.32[B529.BQ.1.1.32 %notin% voc],"VARIANT"] <- "BQ.1.1.32"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.1.41[B529.BQ.1.1.41 %notin% voc],"VARIANT"] <- "BQ.1.1.41"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.10[B529.BQ.1.10 %notin% voc],"VARIANT"] <- "BQ.1.10"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.11[B529.BQ.1.11 %notin% voc],"VARIANT"] <- "BQ.1.11"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.12[B529.BQ.1.12 %notin% voc],"VARIANT"] <- "BQ.1.12"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.13[B529.BQ.1.13 %notin% voc],"VARIANT"] <- "BQ.1.13"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.14[B529.BQ.1.14 %notin% voc],"VARIANT"] <- "BQ.1.14"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.19[B529.BQ.1.19 %notin% voc],"VARIANT"] <- "BQ.1.19"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.2[B529.BQ.1.2 %notin% voc],"VARIANT"] <- "BQ.1.2"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.3[B529.BQ.1.3 %notin% voc],"VARIANT"] <- "BQ.1.3"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.5[B529.BQ.1.5 %notin% voc],"VARIANT"] <- "BQ.1.5"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.22[B529.BQ.1.22 %notin% voc],"VARIANT"] <- "BQ.1.22"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.23[B529.BQ.1.23 %notin% voc],"VARIANT"] <- "BQ.1.23"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.25[B529.BQ.1.25 %notin% voc],"VARIANT"] <- "BQ.1.25"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.25.1[B529.BQ.1.25.1 %notin% voc],"VARIANT"] <- "BQ.1.25.1"
    src.dat[src.dat$VARIANT %in% B529.BQ.1.28[B529.BQ.1.28 %notin% voc],"VARIANT"] <- "BQ.1.28"
    src.dat[src.dat$VARIANT %in% B529.CQ.2[B529.CQ.2 %notin% voc],"VARIANT"] <- "CQ.2"
    src.dat[src.dat$VARIANT %in% B529.DF.1[B529.DF.1 %notin% voc],"VARIANT"] <- "DF.1"
  }
  if(XBB_agg) {
    src.dat[src.dat$VARIANT %in% XBB.1.5.1[XBB.1.5.1 %notin% voc],"VARIANT"] <- "XBB.1.5.1"
    src.dat[src.dat$VARIANT %in% XBB.1.5.2[XBB.1.5.2 %notin% voc],"VARIANT"] <- "XBB.1.5.2"
    src.dat[src.dat$VARIANT %in% XBB.1.5.4[XBB.1.5.4 %notin% voc],"VARIANT"] <- "XBB.1.5.4"
    src.dat[src.dat$VARIANT %in% XBB.1.5.5[XBB.1.5.5 %notin% voc],"VARIANT"] <- "XBB.1.5.5"
    src.dat[src.dat$VARIANT %in% XBB.1.5.10[XBB.1.5.10 %notin% voc],"VARIANT"] <- "XBB.1.5.10"
    src.dat[src.dat$VARIANT %in% XBB.1.5.11[XBB.1.5.11 %notin% voc],"VARIANT"] <- "XBB.1.5.11"
    src.dat[src.dat$VARIANT %in% XBB.1.5.13[XBB.1.5.13 %notin% voc],"VARIANT"] <- "XBB.1.5.13"
    src.dat[src.dat$VARIANT %in% XBB.1.5.15[XBB.1.5.15 %notin% voc],"VARIANT"] <- "XBB.1.5.15"
    src.dat[src.dat$VARIANT %in% XBB.1.5.16[XBB.1.5.16 %notin% voc],"VARIANT"] <- "XBB.1.5.16"
    src.dat[src.dat$VARIANT %in% XBB.1.5.17[XBB.1.5.17 %notin% voc],"VARIANT"] <- "XBB.1.5.17"
    src.dat[src.dat$VARIANT %in% XBB.1.5.19[XBB.1.5.19 %notin% voc],"VARIANT"] <- "XBB.1.5.19"
    src.dat[src.dat$VARIANT %in% XBB.1.5.20[XBB.1.5.20 %notin% voc],"VARIANT"] <- "XBB.1.5.20"
    src.dat[src.dat$VARIANT %in% XBB.1.5.21[XBB.1.5.21 %notin% voc],"VARIANT"] <- "XBB.1.5.21"
    src.dat[src.dat$VARIANT %in% XBB.1.5.30[XBB.1.5.30 %notin% voc],"VARIANT"] <- "XBB.1.5.30"
    src.dat[src.dat$VARIANT %in% XBB.1.5.31[XBB.1.5.31 %notin% voc],"VARIANT"] <- "XBB.1.5.31"
    src.dat[src.dat$VARIANT %in% XBB.1.5.32[XBB.1.5.32 %notin% voc],"VARIANT"] <- "XBB.1.5.32"
    src.dat[src.dat$VARIANT %in% XBB.1.5.33[XBB.1.5.33 %notin% voc],"VARIANT"] <- "XBB.1.5.33"
    src.dat[src.dat$VARIANT %in% XBB.1.5.35[XBB.1.5.35 %notin% voc],"VARIANT"] <- "XBB.1.5.35"
    src.dat[src.dat$VARIANT %in% XBB.1.5[XBB.1.5 %notin% voc],"VARIANT"] <- "XBB.1.5"
    src.dat[src.dat$VARIANT %in% XBB.1.9.1[XBB.1.9.1 %notin% voc],"VARIANT"] <- "XBB.1.9.1"
    src.dat[src.dat$VARIANT %in% XBB.1.9.2[XBB.1.9.2 %notin% voc],"VARIANT"] <- "XBB.1.9.2"
    src.dat[src.dat$VARIANT %in% EG.1[EG.1 %notin% voc],"VARIANT"] <- "EG.1"
    src.dat[src.dat$VARIANT %in% XBB.1.15[XBB.1.15 %notin% voc],"VARIANT"] <- "XBB.1.15"
    src.dat[src.dat$VARIANT %in% XBB.1.16[XBB.1.16 %notin% voc],"VARIANT"] <- "XBB.1.16"
    src.dat[src.dat$VARIANT %in% XBB.1.16.1[XBB.1.16.1 %notin% voc],"VARIANT"] <- "XBB.1.16.1"
    src.dat[src.dat$VARIANT %in% XBB.1[XBB.1 %notin% voc],"VARIANT"] <- "XBB.1"
    src.dat[src.dat$VARIANT %in% XBB.2[XBB.2 %notin% voc],"VARIANT"] <- "XBB.2"
    src.dat[src.dat$VARIANT %in% XBB.2.3[XBB.2.3 %notin% voc],"VARIANT"] <- "XBB.2.3"
    src.dat[src.dat$VARIANT %in% XBB[XBB %notin% voc],"VARIANT"] <- "XBB"
  }
} # end old aggregation methods
# src.dat[ , table( VARIANT0 == VARIANT ) ]

# define an "R346T" variant.
# Note: "R346T" is not reflected in "voc_lut", but will be in "lineage_aggregations" table saved below
if(force_aggregate_R346T==TRUE & length(grep("^R346T", unique(src.dat$VARIANT), perl=T))> 0){
  src.dat[src.dat$VARIANT %in% R346T_sublineages, "VARIANT"] <- "R346T"
}

# create another column for the varients of interest
# this is only used to get (unweighted) counts of the sequences by lineage (used in all runs)
src.dat$VARIANT2 = as.character(src.dat$VARIANT)
# group all non-"variants of interest" together
src.dat[src.dat$VARIANT %notin% voc, "VARIANT2"] <- "Other"

# create a table of all the old and new (after aggregating) variant names
# (using "VARIANT", which does NOT have the "OTHER" aggregation)
write.csv(
  x = src.dat[
    , # no filtering
    .(old = lineage, new = VARIANT), # create new columns with clearer names
    by = c('lineage', 'VARIANT') # group by lineage and Variant
  ][
    ,
    .(old, new) # select only these columns
  ][
    order(old), # reorder by the original lineage name
  ],
  file = paste0(script.basename,
                output_folder, "/lineage_aggregations_",
                ci.type,
                "CI_",
                svy.type,
                "_",
                data_date,
                tag,
                ".csv"),
  row.names = F
)

# create a table of all the old and final (after aggregating) variant names
# (using "VARIANT2", which includes the "OTHER" aggregation)
write.csv(
  x = src.dat[
    , # no filtering
    .(original = lineage, aggregated = VARIANT2), # create new columns with clearer names
    by = c('lineage', 'VARIANT2') # group by lineage and Variant
  ][
    ,
    .(original, aggregated) # select only these columns
  ][
    order(original), # reorder by the original lineage name
  ],
  file = paste0(script.basename,
                output_folder, "/lineage_aggregations_summary_",
                ci.type,
                "CI_",
                svy.type,
                "_",
                data_date,
                tag,
                ".csv"),
  row.names = F
)



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

# tally sequences by source
# 1. tagged as surveillance
# 2. CDC contracting labs
# 3. NS3
# alternatively, tally them using SQL on HUE: https://cdp-01.biotech.cdc.gov:8889/hue/editor?editor=37419
# (This isn't calculated on a regular basis; it's just whenever we want to look at it again)
if(FALSE){
  # identify contractor sequences
  src.dat[ SOURCE %in% c("UW VIROLOGY LAB",
                         "FULGENT GENETICS",
                         "HELIX",
                         "HELIX/ILLUMINA",
                         "LABORATORY CORPORATION OF AMERICA",
                         "AEGIS SCIENCES CORPORATION",
                         "QUEST DIAGNOSTICS INCORPORATED",
                         "BROAD INSTITUTE",
                         "INFINITY BIOLOGIX",
                         "MAKO MEDICAL"),
           'source_type' := 'Contractor']
  # identify NS3 sequences
  src.dat[ SOURCE %in% c("NS3"),
           'source_type' := 'NS3']
  # identify tagged sequences
  src.dat[ is.na(source_type), 'source_type' := 'Tagged']

  # tally totals by source_type
  src.dat[ ,
           .('N' = .N,
             'Prop' = .N / nrow(src.dat)),
           by = c('source_type')]

  # compare to results from HUE on 2022-05-31
  data.frame(
    'source_type' = c('Contractor', 'NS3', 'Tagged'),
    'N' = c(1643231, 31781, 529791),
    'Prop' = c(1643231, 31781, 529791)/sum(c(1643231, 31781, 529791)))

  # plot proportion of totals by week & source_type
  merge(
    x = src.dat[ ,
                 .('N' = .N),
                 by = c('source_type', 'week')],
    y = src.dat[ ,
                 .('N_total' = .N),
                 by = c('week')],
    by =  'week'
  ) %>%
    mutate(Prop = N / N_total) %>%
    ggplot(., aes(x = week, y = Prop, color = source_type)) + theme_bw() + geom_line() + scale_y_continuous(labels = scales::percent_format())
}


### model weeks ----------------------------------------------------------------
# add in another column for model_week to use in the model
# the reason for this is b/c large values of "week" were causing non-invertible hessians in the nowcast model.
# could scale "week" using mean and SD. This is just another option.
# NOTE! When switching to this method, I also switched to fitting the model to
#       20 weeks instead of 21 weeks. (then when switching back to 21 weeks, it screwed up this method a little)
# maximum week (not maximum "model_week") included in the model
model_week_max = as.numeric(as.Date(model_time_end) - week0day1) %/% 7
### SHOULD THIS BE BASED ON TIME_END OR ON DATA_DATE???

# first week (not first "model_week") included in the model
model_week_min = model_week_max - model_weeks
# a midpoint week that will be used to center week values
# this is also the maximum value of "model_week"
model_week_mid = round(model_weeks/2) # center it around 0 instead of using 1:model_weeks
# create a dataframe of old and new values to help visualize things
model_week_df = data.frame(week       = 1:model_weeks + model_week_min,
                           model_week = 1:model_weeks - model_week_mid,
                           week_start = (1:model_weeks + model_week_min)*7 + as.Date(week0day1),
                           week_mid   = (1:model_weeks + model_week_min)*7 + as.Date(week0day1) + 3,
                           week_end   = (1:model_weeks + model_week_min)*7 + as.Date(week0day1) + 6)
# add the new model_week info to the dataframe
src.dat$model_week = src.dat$week - model_week_min - model_week_mid

# data week (this might be the same as the last week in "model_week_df", or it might be 1 week after; it depends on "time_end")
data_week_df = data.frame(
  week       = current_week, # "current_week" is defined in config.R: as.numeric(data_date - week0day1) %/% 7
  model_week = current_week - model_week_min - model_week_mid,
  week_start = (data_date - as.numeric(format(data_date, '%w'))),
  week_mid   = (data_date - as.numeric(format(data_date, '%w'))) + 3,
  week_end   = (data_date - as.numeric(format(data_date, '%w'))) + 6
)

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
if ( grepl("Run(1|2)", tag) ){ # fortnight and weekly estimates

  # subset data to only include
  # - data since May 8, 2021
  # - older than "time_end"
  dat2 <- subset(x = src.dat,
                 as.Date(FORTNIGHT_END) >= (time_start_weights - 7*((as.numeric(time_end+1 - time_start_weights)/7) %%2)) & # this is an ugly way to make sure the start date is a multiple of 2 weeks.
                   as.Date(FORTNIGHT_END) <= time_end)

  # get the relevant fortnights from the data
  # this should be sorted
  ftnts = sort(unique(dat2$FORTNIGHT_END))

  # skip running this if only running/testing the nowcast model.
  # (This section takes 90+% of the time required to run this script.)
  if(!nowcast_only){
    # calculate UNweighted estimates with binomial confidence intervals
    if(TRUE){

      # national estimates
      unwt.ftnt_us <- dat2[,
                           .(count = .N,
                             denom_count = sum(dat2$FORTNIGHT_END == FORTNIGHT_END),
                             Unwt_share = .N / sum(dat2$FORTNIGHT_END == FORTNIGHT_END),
                             USA_or_HHSRegion = 'USA'),
                           by = c('FORTNIGHT_END', 'VARIANT2')]

      # regional estimates
      unwt.ftnt_hhs <- dat2[,
                            .(count = .N,
                              denom_count = sum(dat2$HHS == HHS & dat2$FORTNIGHT_END == FORTNIGHT_END),
                              Unwt_share = .N / sum(dat2$HHS == HHS & dat2$FORTNIGHT_END == FORTNIGHT_END),
                              USA_or_HHSRegion = HHS),
                            by = c('HHS', 'FORTNIGHT_END', 'VARIANT2')]

      # join the national and regional estimates
      unwt.ftnt <- rbind(unwt.ftnt_us,
                         unwt.ftnt_hhs[,.SD, .SD = names(unwt.ftnt_us)] )

      # calculate binomial confidence intervals
      unwt.ftnt[,Unwt_share_lo:= NA_real_]
      unwt.ftnt[,Unwt_share_hi:= NA_real_]
      for(i in 1:nrow(unwt.ftnt)) {
        nsuccess <- unwt.ftnt[i][['count']]
        ntry     <- unwt.ftnt[i][['denom_count']]

        if(ntry > 0){
          # perform a proportion test to get the binomial confidence interval
          pt <- prop.test(x = nsuccess,
                          n = ntry,
                          conf.level = 0.95)

          # add the confidence limits back into the data
          unwt.ftnt[i,Unwt_share_lo:= pt$conf.int[1]]
          unwt.ftnt[i,Unwt_share_hi:= pt$conf.int[2]]
        }
      }# end for loop

      # calculate growth rates
      # make sure sort columns are not factors
      unwt.ftnt[, ':='('USA_or_HHSRegion' = factor(as.character(USA_or_HHSRegion), levels = c('USA', 1:10)),
                       'Fortnight_ending' = as.Date(as.character(FORTNIGHT_END)),
                       'Variant' = as.character(VARIANT2))]

      # calculate some estimate of growth rate in a more efficient manner
      # set the key to order the dataframe & speed up some calculations
      data.table::setkey(unwt.ftnt, USA_or_HHSRegion, Variant, Fortnight_ending)

      # calculate the time period between 2 estimates
      unwt.ftnt[,
                time_diff_days := as.numeric(as.Date(Fortnight_ending)) - as.numeric(as.Date(data.table::shift(Fortnight_ending, 1, type = "lag"))),
                by = c("USA_or_HHSRegion", "Variant")
      ]
      # add a column for the starting Variant share
      unwt.ftnt[,start_share := data.table::shift(Unwt_share, 1, type = "lag"),by = c("USA_or_HHSRegion", "Variant")]

      unwt.ftnt[
        ,':='(
          # calculate multinomial growth rate (weekly)
          unwt_growth_rate = (exp((log(Unwt_share) - log(start_share))/(time_diff_days/7))-1)*100,
          # calculate log-odds growth rate
          unwt_growth_rate_logodds = log( (Unwt_share / (1 - Unwt_share)) / (start_share / (1 - start_share))) / (time_diff_days/7))
      ]
      unwt.ftnt[
        ,':='(
          # calculate multinomial doubling time (in days)
          unwt_doubling_time = (log(2) / log((100 + unwt_growth_rate) / 100)) * 7,
          # calculate time for log odds to double (this is equivalent to doubling time calculated from the time coefficient in a logistic model)
          unwt_doubling_time_logodds = log(2) / unwt_growth_rate_logodds * 7
        )
      ]

      # convert to character (for merging)
      unwt.ftnt[, 'Fortnight_ending' := as.character(Fortnight_ending)]
    } # end calculated unweighted proportions

    # calculate weighted and/or unweighted estimates with survey design
    for (meth in weighted_methods){

      # specify the weights in the survey design
      survey_Design_temp <- svyDES
      if (meth == 'unweighted'){
        # set all the weights to be equal
        survey_Design_temp$prob <- rep(1, length(weights(svyDES)))
        # weights(survey_Design_temp, type = 'sampling')
        # weights(survey_Design_temp, type = 'analysis')
      }
      # create a dataframe with all unique combinations of variants, fortnights, and regions
      # (and then reverse column order for convenience)
      all.ftnt = expand.grid(Variant = voc,
                             Fortnight_ending = ftnts,
                             USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

      # get the proportion estimates & CI
      all_ftnt_ests_function <- function(all.ftnt, svyDES, calc_99_CI = calc_99_CI_weighted){
        # calculate estimates with 95% CI
        ests <- apply(X = all.ftnt,
                      MARGIN = 1,
                      FUN = function(rr) myciprop(voc   = rr[3],
                                                  geoid = rr[1],
                                                  svy = subset(x = svyDES,
                                                               FORTNIGHT_END == rr[2]),
                                                  str = FALSE))
        # optionally calculate estimates with 99% CI
        if(calc_99_CI){
          ests_99 <- apply(X = all.ftnt,
                           MARGIN = 1,
                           FUN = function(rr) myciprop(voc = rr[3],
                                                       geoid = rr[1],
                                                       svy = subset(x = svyDES,
                                                                    FORTNIGHT_END == rr[2]),
                                                       str = FALSE,
                                                       level = 0.99))

          # add the 99% CI onto the ests object
          rownames(ests_99) <- sub(pattern = 'ucl', replacement = 'ucl_99', x = sub(pattern = 'lcl', replacement = 'lcl_99', x = rownames(ests_99)))

          # add the 99% CI to the output
          ests <- rbind(ests,
                        ests_99[c('lcl_99', 'ucl_99'),])
        }

        return(ests)
      } # end all_ftnt_ests_function definition

      # do the estimates in parallel
      if(use_parallel){

        # use a tryCatch just in case the parallel operation fails
        ests <- tryCatch(expr = { # "expr" is what we want to run (not in a function form, unlike "error" and "warning")
          # choose the number of cores
          # ncores <- max(1, parallel::detectCores())

          # split the rows of data into "ncores" sets
          cut_list <- split(x = all.ftnt,
                            f = rep(1:ncores,
                                    each = ceiling(nrow(all.ftnt)/ncores),
                                    length.out = nrow(all.ftnt)))

          # make a cluster
          cl <- parallel::makeCluster(ncores)

          # pass everything to each of the cores
          parallel::clusterEvalQ(cl = cl, {
            library(survey)
            library(data.table)
          })
          # export all the other R objects that are used within all_ftnt_ests_function to each cluster node
          parallel::clusterExport(cl = cl,
                                  varlist = c('survey_Design_temp', 'script.basename', 'ci.type', 'voc', 'calc_99_CI_weighted'))
          parallel::clusterEvalQ(cl = cl, {
            source(paste0(script.basename, "/weekly_variant_report_functions.R"))
            options(survey.adjust.domain.lonely = T,
                    survey.lonely.psu = "average",
                    stringsAsFactors = FALSE)
          })

          # perform the calculations on the cluster
          ests_list <- parallel::parLapply(cl = cl,
                                           X = cut_list,
                                           fun = all_ftnt_ests_function,
                                           svyDES = survey_Design_temp)

          # combine the results into a single dataframe
          do.call(what = 'cbind', args = ests_list)
        },
        error = function(cond) { # if "expr" throws an error, do this instead of actually stopping
          message('Parallel execution of biweekly weighted proportions failed. Running in series instead.')
          message("Here's the original error message:")
          message(cond)
          # Choose a return value in case of error
          return(all_ftnt_ests_function(all.ftnt = all.ftnt, svyDES = survey_Design_temp))
        })
      } else ests <- all_ftnt_ests_function(all.ftnt = all.ftnt, svyDES = survey_Design_temp)

      # add in the estimates to the dataframe
      all.ftnt = cbind(all.ftnt,
                       Share    = ests[1,],
                       Share_lo = ests[2,],
                       Share_hi = ests[3,],
                       DF       = ests[4,],
                       eff.size = ests[5,],
                       cv.mean  = ests[6,],
                       deff     = ests[7,])

      # add on 99% CI
      if('lcl_99' %in% row.names(ests)){
        all.ftnt <- cbind(all.ftnt,
                          Share_lo_99 = ests['lcl_99',],
                          Share_hi_99 = ests['ucl_99',])
      }

      ## make predictions for the "other" variants (following the same steps)
      others = expand.grid(Variant          = "Other",
                           Fortnight_ending = ftnts,
                           USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

      # function to get the proportion estimates & CI (for "other" variants)
      others_ftnt_ests_function <- function(others, svyDES, calc_99_CI = calc_99_CI_weighted){
        # calculate estimates with 95% CI
        ests.others = apply(X = others,
                            MARGIN = 1,
                            FUN = function(rr) myciprop(voc = voc, # "Other" isn't in the Nowcast model; if voc is a vector, myciprop will return the aggregated proportion. So feed in all the vocs, then subtract from 1 to get "Other" estimates.
                                                        geoid = rr[1],
                                                        svy = subset(svyDES,
                                                                     FORTNIGHT_END == rr[2]),
                                                        str = FALSE))

        # optionally calculate estimates with 99% CI
        if(calc_99_CI){
          ests.others_99 <- apply(X = others,
                                  MARGIN = 1,
                                  FUN = function(rr) myciprop(voc = voc, # "Other" isn't in the Nowcast model; if voc is a vector, myciprop will return the aggregated proportion. So feed in all the vocs, then subtract from 1 to get "Other" estimates.
                                                              geoid = rr[1],
                                                              svy = subset(svyDES,
                                                                           FORTNIGHT_END == rr[2]),
                                                              str = FALSE,
                                                              level = 0.99))

          # add the 99% CI onto the ests object
          rownames(ests.others_99) <- sub(pattern = 'ucl',
                                          replacement = 'ucl_99',
                                          x = sub(pattern = 'lcl',
                                                  replacement = 'lcl_99',
                                                  x = row.names(ests.others_99)))

          ests.others <- rbind(ests.others,
                               ests.others_99[c('lcl_99', 'ucl_99'),])
        }


        return(ests.others)
      } # end others_ftnt_ests_function definition

      if(use_parallel){
        # use a tryCatch just in case the parallel operation fails
        ests.others <- tryCatch(expr = { # "expr" is what we want to run (not in a function form, unlike "error" and "warning")
          # split the rows of data into "ncores" sets
          cut_list <- split(x = others,
                            f = rep(1:ncores,
                                    each = ceiling(nrow(others)/ncores),
                                    length.out = nrow(others)))

          # perform the calculations on the cluster
          others_list <- parallel::parLapply(cl = cl,
                                             X = cut_list,
                                             fun = others_ftnt_ests_function,
                                             svyDES = survey_Design_temp)

          # combine the results into a single dataframe
          do.call(what = 'cbind', args = others_list)
        },
        error = function(cond) { # if "expr" throws an error, do this instead of actually stopping
          message('Parallel execution of biweekly "other" weighted proportions failed. Running in series instead.')
          message("Here's the original error message:")
          message(cond)
          # Choose a return value in case of error
          return(others_ftnt_ests_function(others = others, svyDES = survey_Design_temp))
        })
      } else ests.others <- others_ftnt_ests_function(others = others, svyDES = survey_Design_temp)

      # add in the estimates to the dataframe (for "other" variants)
      others = cbind(others,
                     Share    = 1-ests.others[1,],
                     Share_lo = 1-ests.others[3,],
                     Share_hi = 1-ests.others[2,],
                     DF       = ests.others[4,],
                     eff.size = ests.others[5,],
                     cv.mean  = ests.others[6,],
                     deff     = ests.others[7,])

      # add on 99% CI
      if('lcl_99' %in% row.names(ests.others)){
        others <- cbind(others,
                        Share_lo_99 = 1-ests.others['ucl_99',],
                        Share_hi_99 = 1-ests.others['lcl_99',])
      }


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

      # estimate week-over-week growth rate
      # NOTE! This is NOT equivalent to the Nowcast WoW growth rates.
      #       Nowcast weekly growth rates are instantaneous rates at the weekly midpoint (Wednesday)
      #       These estimates are averaged over the prior week/fortnight and will therefore always be
      #       slightly higher than the Nowcast growth rates.
      {
        # # tidyverse syntax
        # all.ftnt2 %>%
        #   tibble::as_tibble() %>%
        #    dplyr::arrange(USA_or_HHSRegion, Variant, Fortnight_ending) %>%
        #    dplyr::group_by(USA_or_HHSRegion, Variant) %>%
        #    dplyr::mutate(
        #       # calculate the time difference between the two estimates
        #       time_diff_days = as.numeric(as.Date(Fortnight_ending) - as.Date(dplyr::lag(Fortnight_ending))),
        #       # calculate growth rate that's equivalent to Nowcast (multinomial logistic regression) growth rate
        #       # multinomial model growth rate is the derivative (slope) of the log of the proportion;
        #       # - convert proportion to log proportion
        #       # - get the change in log proportion
        #       # - divide by number of weeks in the time period
        #       # - convert back to proportion scale
        #       # so that it's comparable to our other estimates
        #       growth_rate = (exp((log(Share) - dplyr::lag(log(Share)))/(time_diff_days/7))-1)*100,
        #       # calculate growth rate that's equivalent to logistic regression growth rate
        #       # logistic regression growth rate is the change in log odds per unit time,
        #       # - convert proportion to log odds
        #       # - get the change in log odds
        #       # - divide time period by number of weeks in the time period (to get growth rate per week)
        #       logistic_growth_rate = ((log(Share/(1-Share)) - dplyr::lag(log(Share/(1-Share))))/(time_diff_days/7)),
        #       # multinomial doubling time (in days)
        #       doubling_time = (log(2) / log((100 + growth_rate) / 100)) * 7,
        #       # logistic doubling time (in days)
        #       logistic_doubling_time = log(2) / logistic_growth_rate * 7
        #    ) %>%
        #    dplyr::ungroup()

        # data.table syntax
        all.ftnt2 <- data.table::as.data.table(all.ftnt2)

        # make sure sort columns are not factors
        all.ftnt2[, ':='('USA_or_HHSRegion' = factor(as.character(USA_or_HHSRegion), levels = c('USA', 1:10)),
                         'Fortnight_ending' = as.Date(as.character(Fortnight_ending)),
                         'Variant' = as.character(Variant))]

        # calculate some estimate of growth rate in a more efficient manner
        # set the key to order the dataframe & speed up some calculations
        data.table::setkey(all.ftnt2, USA_or_HHSRegion, Variant, Fortnight_ending)

        # calculate the time period between 2 estimates
        all.ftnt2[
          ,
          time_diff_days := as.numeric(as.Date(Fortnight_ending)) -
            as.numeric(as.Date(data.table::shift(Fortnight_ending, 1, type = "lag"))),
          by = c("USA_or_HHSRegion", "Variant")
        ]
        # add a column for the starting Variant share
        all.ftnt2[
          ,
          start_share := data.table::shift(Share, 1, type = "lag"),
          by = c("USA_or_HHSRegion", "Variant")
        ]
        # calculate multinomial growth rate (weekly)
        # This is an exponential growth rate.
        # https://en.wikipedia.org/wiki/Doubling_time#Cell_culture_doubling_time:~:text=Growth%20rate%3A,or
        # exponential growth rate = log(Share / start_share)/(time_diff_days/7)
        #    continuous time: xt = x0*exp(r*t) (r = proportion of e)
        #                  xt/x0 = exp(r*t)
        #             log(xt/x0) = r*t
        #           log(xt/x0)/t = r
        #
        # But then we take the exponent of it so that the growth rate is on the
        # scale of proportions. This essentially converts it to discrete time.
        # continuous time: xt = x0*exp(r*t) (r = proportion of e)
        #                  xt = x0*exp(r)^t
        # discrete time:   xt = x0*(1+r)^t  (r = proportion of x0)
        # so we can convert the continuous time growth rate, r, to the discrete time growth rate, using the equivalence above: exp(r) = (1+r)
        #
        # discrete time:   xt = x0*(1+r)^t  (r = proportion of x0)
        #               xt/x0 = (1+r)^t
        #       (xt/x0)^(1/t) = (1+r)
        #   (xt/x0)^(1/t) - 1 = r
        #
        # multinomial model growth rate is the derivative (slope) of the log of the proportion;
        # - convert proportion to log proportion
        # - get the change in log proportion
        # - divide by number of weeks in the time period
        # - convert back to proportion scale, (exponent() - 1) * 100
        # so that it's comparable to our other estimates
        # Note! This is the average growth rate over the previous week/fortnight
        #       and will therefore be > the Nowcast growth rate estimate for the same week
        all.ftnt2[
          ,
          growth_rate := (exp((log(Share) - log(start_share))/(time_diff_days/7))-1)*100
          # could calculate using the discrete time formula
          # growth_rate = ((Share/start_share)^(1/(time_diff_days/7))-1)*100
          # continuous time: (then exponentiate, subtract 1, *100)
          # growth_rate = (exp( log(Share / start_share)/(time_diff_days/7))-1) *100 )
        ]
        # all.ftnt2[!is.na(growth_rate),.(Share,start_share,time_diff_days, growth_rate, growth_rate2, growth_rate3)]

        # calculate logistic growth rate (log-odds) (weekly)
        all.ftnt2[
          ,
          growth_rate_logodds := log( (Share / (1 - Share)) / (start_share / (1 - start_share))) / (time_diff_days/7)
        ]
        # all.ftnt2[!is.na(logistic_growth_rate),.(Share,start_share,time_diff_days, growth_rate_logodds)]

        # calculate multinomial doubling time (in days)
        all.ftnt2[
          ,
          # https://en.wikipedia.org/wiki/Doubling_time
          # Countinuous time (https://en.wikipedia.org/wiki/Doubling_time#Cell_culture_doubling_time)
          # doubling time = log(2)/r where r = multiplicative (exponential) growth rate
          # A  = P*e^(r*t)
          # 2P = P*e^(r*t)
          # 2  = e^(r*t)
          # log(2) = r*t
          # log(2)/r = t (r=0 is no growth)
          #
          # or....
          # Discrete time
          # Doubling time = log(2) / log(r+1) where r = (exponential) growth rate
          # y   = y0*(1+r)^t
          # 2y0 = y0*(1+r)^t
          # 2   = (1+r)^t
          # log(2) = log((1+r)^t)
          # log(2) = t*log(1+r)
          # log(2)/log(r+1) = t (r=0 is no growth)
          #
          # log((100+growth_rate)/100) = exponential (multiplicative) growth rate
          # log(2) / log((100+growth_rate)/100) = number of weeks per doubling
          # multiply by days / week to get the number of days per doubling
          doubling_time := (log(2) / log((100 + growth_rate) / 100)) * 7
          # equivalent to
          # doubling_time := log(2)/ (continuous gr)      = log(2) / (log(Share/start_share) / (time_diff_days/7)) * 7
          # doubling_time := log(2)/ log(discrete gr + 1) = log(2) / log((Share/start_share)^(1/time_diff_days) - 1 + 1) [divide by 7 to convert to weeks then multiply by 7 to convert back to days]
        ]
        # all.ftnt2[!is.na(doubling_time1),.(Share,start_share,time_diff_days, doubling_time1, doubling_time3)]

        # calculate time for log odds to double (this is equivalent to doubling time calculated from the time coefficient in a logistic model)
        all.ftnt2[,doubling_time_logodds := log(2) / growth_rate_logodds * 7]

        # # remove infinite growth rates
        # all.ftnt2[is.infinite(growth_rate), growth_rate := NA]


        # remove temporary columns
        all.ftnt2[, time_diff_days := NULL]
        all.ftnt2[, start_share    := NULL]
        # convert time back to character for merging
        all.ftnt2[, Fortnight_ending := as.character(Fortnight_ending)]
      } # end calc w-o-w growth rate

      # select columns for the final results
      all.ftnt2_column_names <- c("USA_or_HHSRegion",
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
                                  "nchs_flag_wodf",
                                  'growth_rate',
                                  'growth_rate_logodds',
                                  'doubling_time',
                                  'doubling_time_logodds')
      # possibly add in 99% CI
      if('Share_lo_99' %in% colnames(all.ftnt2)) all.ftnt2_column_names <- c(all.ftnt2_column_names, 'Share_lo_99', 'Share_hi_99')

      # subset the columns ( the .. tells data.table that "all.ftnt2_column_names" is an object in the global environment, NOT a column in the data.table)
      all.ftnt2 = all.ftnt2[, ..all.ftnt2_column_names]

      # optionally calculate the number of infections attributable to each variant
      if (calc_confirmed_infections){
        test_filepath <- paste0(script.basename,
                                "/data/backup_",
                                data_date, custom_tag, "/",
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

          # fortnightly tests
          tests_fn <- rbind(tests_fn_us,
                            tests_fn_hhs)
          # convert time to character for merging
          tests_fn[, fortnight_end := as.character(fortnight_end)]

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

      } # end calc_confirmed_infections

      # re-order results by HHS region [so that "other" variants are not listed seperately]
      all.ftnt2 <- all.ftnt2[order(all.ftnt2$USA_or_HHSRegion),]

      # make sure both Fortnight_ending columns are the same type
      all.ftnt2[, 'Fortnight_ending' := as.character(Fortnight_ending)]
      unwt.ftnt[, 'Fortnight_ending' := as.character(Fortnight_ending)]

      # merge on the unweighted estimates
      all.ftnt3 <- merge(
        x = all.ftnt2,
        y = unwt.ftnt[,.(USA_or_HHSRegion, Fortnight_ending, Variant,
                         Unwt_share, Unwt_share_lo, Unwt_share_hi,
                         unwt_growth_rate, unwt_growth_rate_logodds,
                         unwt_doubling_time, unwt_doubling_time_logodds)],
        by = c('USA_or_HHSRegion', 'Fortnight_ending', 'Variant'),
        all.x = TRUE
      )


      # output data file depends on NO change in all.ftnt3 column order!!!
      # If new columns are to be added, they need to go to the end.
      # write results to file
      write.csv(x = all.ftnt3,
                file = paste0(script.basename,
                              output_folder, "/variant_share_",
                              meth, # weighted or unweighted
                              "_",
                              ci.type,
                              "CI_",
                              svy.type,
                              "_",
                              data_date,
                              tag,
                              "_",
                              results_tag,
                              ".csv"),
                row.names = FALSE)

      if(meth == 'weighted'){
        # process dataframe and save to a format for direct hadoop upload
        all.ftnt3_hadoop <- data.frame(all.ftnt3[,c("USA_or_HHSRegion",
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
                                                    "nchs_flag_wodf")])
        all.ftnt3_hadoop[is.na(all.ftnt3_hadoop)] = '\\N'
        all.ftnt3_hadoop[,14:15] = '\\N'
        all.ftnt3_hadoop[,16] = 'weighted'
        all.ftnt3_hadoop[,17] = 'biweekly'
        all.ftnt3_hadoop[,18] = data_date
        all.ftnt3_hadoop[,19] = paste0(results_tag, '_Run', opts$run_number)
        all.ftnt3_hadoop[,20] = 1
        all.ftnt3_hadoop[,21:24] = all.ftnt3[,c("cases", "cases_hi", "cases_lo")]
        # all.ftnt3_hadoop[,3] = gsub('Delta Aggregated', 'B.1.617.2', all.ftnt3_hadoop[,3])
        # all.ftnt3_hadoop[,3] = gsub('Omicron Aggregated', 'B.1.1.529', all.ftnt3_hadoop[,3])
        # all.ftnt3_hadoop[,3] = gsub(' Aggregated', '', all.ftnt3_hadoop[,3])
        write.table(x = all.ftnt3_hadoop,
                    file = paste0(script.basename,
                                  output_folder, "/variant_share_weighted_",
                                  ci.type,
                                  "CI_",
                                  svy.type,
                                  "_",
                                  data_date,
                                  tag,
                                  "_",
                                  results_tag,
                                  "_hadoop.csv"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE,
                    sep =  ",")
      }

      ### Weekly estimates
      # (repeat of above code, but with weeks instead of fortnights)

      include_run1_weekly <- TRUE
      if (include_run1_weekly){

        # subset data to only include data since 2 May, 2021
        # and (end of week) older than "time_end"
        dat2 <- subset(src.dat,
                       as.Date(yr_wk) >= time_start_weights &
                         (as.Date(yr_wk) + 6) <= time_end)

        # add a column for the date of the final day of each week
        dat2[, WEEK_END := as.Date(yr_wk) + 6]

        # get the unique weeks that are in the appropriate time frame
        wks = sort(unique(dat2$yr_wk))

        # calculate raw/unweighted estimates (NO survey design; binomial CI)
        {
          # national estimates
          unwt.wkly_us <- dat2[,
                               .(count = .N,
                                 denom_count = sum(dat2$WEEK_END == WEEK_END),
                                 Unwt_share = .N / sum(dat2$WEEK_END == WEEK_END),
                                 USA_or_HHSRegion = 'USA'),
                               by = c('WEEK_END', 'VARIANT2')]

          # regional estimates
          unwt.wkly_hhs <- dat2[,
                                .(count = .N,
                                  denom_count = sum(dat2$HHS == HHS & dat2$WEEK_END == WEEK_END),
                                  Unwt_share = .N / sum(dat2$HHS == HHS & dat2$WEEK_END == WEEK_END),
                                  USA_or_HHSRegion = HHS),
                                by = c('HHS', 'WEEK_END', 'VARIANT2')]

          # join the national and regional estimates
          unwt.wkly <- rbind(unwt.wkly_us,
                             unwt.wkly_hhs[,.SD, .SD = names(unwt.wkly_us)] )

          # calculate binomial confidence intervals
          unwt.wkly[,Unwt_share_lo:= NA_real_]
          unwt.wkly[,Unwt_share_hi:= NA_real_]
          for(i in 1:nrow(unwt.wkly)) {
            nsuccess <- unwt.wkly[i][['count']]
            ntry     <- unwt.wkly[i][['denom_count']]

            if(ntry > 0){
              # perform a proportion test to get the binomial confidence interval
              pt <- prop.test(x = nsuccess,
                              n = ntry,
                              conf.level = 0.95)

              # add the confidence limits back into the data
              unwt.wkly[i,Unwt_share_lo:= pt$conf.int[1]]
              unwt.wkly[i,Unwt_share_hi:= pt$conf.int[2]]
            }
          }# end for loop

          # calculate growth rates
          # make sure sort columns are not factors
          unwt.wkly[, ':='('USA_or_HHSRegion' = factor(as.character(USA_or_HHSRegion), levels = c('USA', 1:10)),
                           'WEEK_END' = as.Date(as.character(WEEK_END)),
                           'Variant' = as.character(VARIANT2))]

          # calculate some estimate of growth rate in a more efficient manner
          # set the key to order the dataframe & speed up some calculations
          data.table::setkey(unwt.wkly, USA_or_HHSRegion, Variant, WEEK_END)

          # calculate the time period between 2 estimates
          unwt.wkly[,
                    time_diff_days := as.numeric(as.Date(WEEK_END)) - as.numeric(as.Date(data.table::shift(WEEK_END, 1, type = "lag"))),
                    by = c("USA_or_HHSRegion", "Variant")
          ]
          # add a column for the starting Variant share
          unwt.wkly[,start_share := data.table::shift(Unwt_share, 1, type = "lag"),by = c("USA_or_HHSRegion", "Variant")]

          unwt.wkly[
            ,':='(
              # calculate multinomial growth rate (weekly)
              unwt_growth_rate = (exp((log(Unwt_share) - log(start_share))/(time_diff_days/7))-1)*100,
              # calculate log-odds growth rate
              unwt_growth_rate_logodds = log( (Unwt_share / (1 - Unwt_share)) / (start_share / (1 - start_share))) / (time_diff_days/7))
          ]
          unwt.wkly[
            ,':='(
              # calculate multinomial doubling time (in days)
              unwt_doubling_time = (log(2) / log((100 + unwt_growth_rate) / 100)) * 7,
              # calculate time for log odds to double (this is equivalent to doubling time calculated from the time coefficient in a logistic model)
              unwt_doubling_time_logodds = log(2) / unwt_growth_rate_logodds * 7
            )
          ]

          # convert time to character for merging
          unwt.wkly[, WEEK_END := as.character(WEEK_END)]
        } # end calculated unweighted proportions






        # create a dataframe with variant, time period, and region
        # (and reorder columns for convenience)
        all.wkly = expand.grid(Variant = voc,
                               Week_of = wks,
                               USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

        # function to get predictions & CI for each variant in each week and region
        # (using survey design but NOT multinomial nowcast model)
        all_wkly_ests_function <- function(all.wkly, svyDES, calc_99_CI = calc_99_CI_weighted){
          # calculate estimates with 95% CI
          ests <- apply(X = all.wkly,
                        MARGIN = 1,
                        FUN = function(rr) myciprop(voc   = rr[3],
                                                    geoid = rr[1],
                                                    svy   = subset(svyDES,
                                                                   yr_wk == rr[2]),
                                                    str   = FALSE))

          # optionally calculate estimates with 99% CI
          if(calc_99_CI){
            ests_99 <- apply(X = all.wkly,
                             MARGIN = 1,
                             FUN = function(rr) myciprop(voc   = rr[3],
                                                         geoid = rr[1],
                                                         svy   = subset(svyDES,
                                                                        yr_wk == rr[2]),
                                                         str   = FALSE,
                                                         level = 0.99))

            # change some of the names to make them distinct from the 95% CI
            rownames(ests_99) <- sub(pattern = 'ucl',
                                     replacement = 'ucl_99',
                                     x = sub(pattern = 'lcl',
                                             replacement = 'lcl_99',
                                             x = rownames(ests_99)))

            # add the 99% CI onto the ests object
            ests <- rbind(ests,
                          ests_99[c('lcl_99', 'ucl_99'),])
          }

          return(ests)
        } # end all_wkly_ests_function definition

        # do the estimates in parallel
        if(use_parallel){
          # use a tryCatch just in case the parallel operation fails
          ests <- tryCatch(expr = { # "expr" is what we want to run (not in a function form, unlike "error" and "warning")
            # split the rows of data into "ncores" sets
            cut_list <- split(x = all.wkly,
                              f = rep(1:ncores,
                                      each = ceiling(nrow(all.wkly)/ncores),
                                      length.out = nrow(all.wkly)))

            # perform the calculations on the cluster
            ests_list <- parallel::parLapply(cl = cl,
                                             X = cut_list,
                                             fun = all_wkly_ests_function,
                                             svyDES = survey_Design_temp)

            # combine the results into a single dataframe
            do.call(what = 'cbind', args = ests_list)
          },
          error = function(cond) { # if "expr" throws an error, do this instead of actually stopping
            message('Parallel execution of weekly weighted proportions failed. Running in series instead.')
            message("Here's the original error message:")
            message(cond)
            # Choose a return value in case of error
            return(all_wkly_ests_function(all.wkly = all.wkly, svyDES = survey_Design_temp))
          })
        } else ests <- all_wkly_ests_function(all.wkly = all.wkly, svyDES = survey_Design_temp)

        # add the predictions into the dataframe
        all.wkly = cbind(all.wkly,
                         Share    = ests[1,],
                         Share_lo = ests[2,],
                         Share_hi = ests[3,],
                         DF       = ests[4,],
                         eff.size = ests[5,],
                         cv.mean  = ests[6,])

        # add on 99% CI
        if('lcl_99' %in% row.names(ests)){
          all.wkly <- cbind(all.wkly,
                            Share_lo_99 = ests['lcl_99',],
                            Share_hi_99 = ests['ucl_99',])
        }

        # get predictions for the non-focal ("other") variants grouped together
        # (and reorder columns for convenience)
        others = expand.grid(Variant = "Other",
                             Week_of = wks,
                             USA_or_HHSRegion = c("USA", 1:10))[, 3:1]

        # function to get the proportion estimates & CI (for "other" variants)
        # get predictions & CI for "other" (grouped) variant in each week and region
        # (using survey design but NOT multinomial nowcast model)
        others_wkly_ests_function <- function(others, svyDES, calc_99_CI = calc_99_CI_weighted){
          # calculate estimates with 95% CI
          ests.others <- apply(X = others,
                               MARGIN = 1,
                               FUN = function(rr) myciprop(voc = voc,
                                                           geoid = rr[1],
                                                           svy = subset(svyDES, yr_wk == rr[2]),
                                                           str =  FALSE))

          # optionally calculate estimates with 99% CI
          if(calc_99_CI){
            ests.others_99 <- apply(X = others,
                                    MARGIN = 1,
                                    FUN = function(rr) myciprop(voc = voc,
                                                                geoid = rr[1],
                                                                svy = subset(svyDES, yr_wk == rr[2]),
                                                                level = .99,
                                                                str =  FALSE))

            # add the 99% CI onto the ests object
            rownames(ests.others_99) <- sub(pattern = 'ucl',
                                            replacement = 'ucl_99',
                                            x = sub(pattern = 'lcl',
                                                    replacement = 'lcl_99',
                                                    x = rownames(ests.others_99)))

            # add the 99% CI onto the results
            ests.others <- rbind(ests.others,
                                 ests.others_99[c('lcl_99', 'ucl_99'),])
          }

          return(ests.others)
        } # end others_wkly_ests_function definition

        if(use_parallel){
          # use a tryCatch just in case the parallel operation fails
          ests.others <- tryCatch(expr = { # "expr" is what we want to run (not in a function form, unlike "error" and "warning")
            # split the rows of data into "ncores" sets
            cut_list <- split(x = others,
                              f = rep(1:ncores,
                                      each = ceiling(nrow(others)/ncores),
                                      length.out = nrow(others)))

            # perform the calculations on the cluster
            others_list <- parallel::parLapply(cl = cl,
                                               X = cut_list,
                                               fun = others_wkly_ests_function,
                                               svyDES = survey_Design_temp)

            # close the cluster
            parallel::stopCluster(cl)

            # combine the results into a single dataframe
            do.call(what = 'cbind', args = others_list)
          },
          error = function(cond) { # if "expr" throws an error, do this instead of actually stopping
            message('Parallel execution of weekly "other" weighted proportions failed. Running in series instead.')
            message("Here's the original error message:")
            message(cond)
            # Choose a return value in case of error
            return(others_wkly_ests_function(others = others, svyDES = survey_Design_temp))
          })
        } else ests.others <- others_wkly_ests_function(others = others, svyDES = survey_Design_temp)

        # add the predictions into the dataframe
        others = cbind(others,
                       Share    = 1-ests.others[1,],
                       Share_lo = 1-ests.others[3,],
                       Share_hi = 1-ests.others[2,],
                       DF       = ests.others[4,],
                       eff.size = ests.others[5,],
                       cv.mean  = ests.others[6,])

        # add on 99% CI
        if('lcl_99' %in% row.names(ests.others)){
          others <- cbind(others,
                          Share_lo_99 = 1-ests.others['ucl_99',],
                          Share_hi_99 = 1-ests.others['lcl_99',])
        }

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

        # estimate week-over-week growth rate
        # NOTE! This is NOT equivalent to the Nowcast WoW growth rates.
        #       Nowcast weekly growth rates are instantaneous rates at the weekly midpoint (Wednesday)
        #       These estimates are averaged over the prior week/fortnight and will therefore always be
        #       slightly higher than the Nowcast growth rates.
        {
          # # tidyverse syntax
          # all.wkly2 %>%
          #   tibble::as_tibble() %>%
          #   dplyr::arrange(USA_or_HHSRegion, Variant, WEEK_END) %>%
          #   dplyr::group_by(USA_or_HHSRegion, Variant) %>%
          #   dplyr::mutate(
          #     # calculate the time difference between the two estimates
          #     time_diff_days = as.numeric(as.Date(WEEK_END) - as.Date(dplyr::lag(WEEK_END))),
          #     # calculate growth rate that's equivalent to Nowcast (multinomial logistic regression) growth rate
          #     # multinomial model growth rate is the derivative (slope) of the log of the proportion;
          #     # - convert proportion to log proportion
          #     # - get the change in log proportion
          #     # - divide by number of weeks in the time period
          #     # - convert back to proportion scale
          #     # so that it's comparable to our other estimates
          #     growth_rate = (exp((log(Share) - dplyr::lag(log(Share)))/(time_diff_days/7))-1)*100,
          #     # calculate growth rate that's equivalent to logistic regression growth rate
          #     # logistic regression growth rate is the change in log odds per unit time,
          #     # - convert proportion to log odds
          #     # - get the change in log odds
          #     # - divide time period by number of weeks in the time period (to get growth rate per week)
          #     logistic_growth_rate = ((log(Share/(1-Share)) - dplyr::lag(log(Share/(1-Share))))/(time_diff_days/7)),
          #     # multinomial doubling time (in days)
          #     doubling_time = (log(2) / log((100 + growth_rate) / 100)) * 7,
          #     # logistic doubling time (in days)
          #     logistic_doubling_time = log(2) / logistic_growth_rate * 7
          #   ) %>%
          #   dplyr::ungroup()

          # data.table syntax
          all.wkly2 <- data.table::as.data.table(all.wkly2)

          # make sure sort columns are not factors
          all.wkly2[, ':='('USA_or_HHSRegion' = factor(as.character(USA_or_HHSRegion), levels = c('USA', 1:10)),
                           'WEEK_END' = as.Date(as.character(WEEK_END)),
                           'Variant' = as.character(Variant))]

          # calculate some estimate of growth rate in a more efficient manner
          # set the key to order the dataframe & speed up some calculations
          data.table::setkey(all.wkly2, USA_or_HHSRegion, Variant, WEEK_END)

          # calculate the time period between 2 estimates
          all.wkly2[
            ,
            time_diff_days := as.numeric(as.Date(WEEK_END)) -
              as.numeric(as.Date(data.table::shift(WEEK_END, 1, type = "lag"))),
            by = c("USA_or_HHSRegion", "Variant")
          ]
          # add a column for the starting Variant share
          all.wkly2[
            ,
            start_share := data.table::shift(Share, 1, type = "lag"),
            by = c("USA_or_HHSRegion", "Variant")
          ]
          # calculate multinomial growth rate (weekly)
          # This is an exponential growth rate.
          # https://en.wikipedia.org/wiki/Doubling_time#Cell_culture_doubling_time:~:text=Growth%20rate%3A,or
          # exponential growth rate = log(Share / start_share)/(time_diff_days/7)
          #    continuous time: xt = x0*exp(r*t) (r = proportion of e)
          #                  xt/x0 = exp(r*t)
          #             log(xt/x0) = r*t
          #           log(xt/x0)/t = r
          #
          # But then we take the exponent of it so that the growth rate is on the
          # scale of proportions. This essentially converts it to discrete time.
          # continuous time: xt = x0*exp(r*t) (r = proportion of e)
          #                  xt = x0*exp(r)^t
          # discrete time:   xt = x0*(1+r)^t  (r = proportion of x0)
          # so we can convert the continuous time growth rate, r, to the discrete time growth rate, using the equivalence above: exp(r) = (1+r)
          #
          # discrete time:   xt = x0*(1+r)^t  (r = proportion of x0)
          #               xt/x0 = (1+r)^t
          #       (xt/x0)^(1/t) = (1+r)
          #   (xt/x0)^(1/t) - 1 = r
          #
          # multinomial model growth rate is the derivative (slope) of the log of the proportion;
          # - convert proportion to log proportion
          # - get the change in log proportion
          # - divide by number of weeks in the time period
          # - convert back to proportion scale, (exponent() - 1) * 100
          # so that it's comparable to our other estimates
          # Note! This is the average growth rate over the previous week/fortnight
          #       and will therefore be > the Nowcast growth rate estimate for the same week
          all.wkly2[
            ,
            growth_rate := (exp((log(Share) - log(start_share))/(time_diff_days/7))-1)*100
            # could calculate using the discrete time formula
            # growth_rate = ((Share/start_share)^(1/(time_diff_days/7))-1)*100
            # continuous time: (then exponentiate, subtract 1, *100)
            # growth_rate = (exp( log(Share / start_share)/(time_diff_days/7))-1) *100 )
          ]

          # calculate logistic growth rate (log-odds) (weekly)
          all.wkly2[
            ,
            # change in log odds: (log(odds_1) - log(odds_0) / # of weeks
            #                   =  log(odds_1/odds_0) / # of weeks
            growth_rate_logodds := log( (Share / (1 - Share)) / (start_share / (1 - start_share))) / (time_diff_days/7)
          ]
          # calculate multinomial doubling time (in days)
          all.wkly2[
            ,
            # https://en.wikipedia.org/wiki/Doubling_time
            # Countinuous time (https://en.wikipedia.org/wiki/Doubling_time#Cell_culture_doubling_time)
            # doubling time = log(2)/r where r = multiplicative (exponential) growth rate
            # A  = P*e^(r*t)
            # 2P = P*e^(r*t)
            # 2  = e^(r*t)
            # log(2) = r*t
            # log(2)/r = t (r=0 is no growth)
            #
            # or....
            # Discrete time
            # Doubling time = log(2) / log(r+1) where r = (exponential) growth rate
            # y   = y0*(1+r)^t
            # 2y0 = y0*(1+r)^t
            # 2   = (1+r)^t
            # log(2) = log((1+r)^t)
            # log(2) = t*log(1+r)
            # log(2)/log(r+1) = t (r=0 is no growth)
            #
            # log((100+growth_rate)/100) = exponential (multiplicative) growth rate
            # log(2) / log((100+growth_rate)/100) = number of weeks per doubling
            # multiply by days / week to get the number of days per doubling
            doubling_time := (log(2) / log((100 + growth_rate) / 100)) * 7
            # equivalent to
            # doubling_time := log(2)/ (continuous gr)      = log(2) / (log(Share/start_share) / (time_diff_days/7)) * 7
            # doubling_time := log(2)/ log(discrete gr + 1) = log(2) / log((Share/start_share)^(1/time_diff_days) - 1 + 1) [divide by 7 to convert to weeks then multiply by 7 to convert back to days]
          ]
          # calculate time for log odds to double (this is equivalent to doubling time calculated from the time coefficient in a logistic model)
          all.wkly2[,doubling_time_logodds := log(2) / growth_rate_logodds * 7]
          # remove infinite growth rates
          # all.wkly2[is.infinite(growth_rate), growth_rate := NA]


          # remove temporary columns
          all.wkly2[, time_diff_days := NULL]
          all.wkly2[, start_share    := NULL]

          # convert time to character for merging
          all.wkly2[, WEEK_END := as.character(WEEK_END)]
        } # end calculate WoW growth rate

        # select the columns to save
        all.wkly2_column_names <- c("USA_or_HHSRegion",
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
                                    "count_LT10",
                                    'growth_rate',
                                    'growth_rate_logodds',
                                    'doubling_time',
                                    'doubling_time_logodds'
        )
        if('Share_lo_99' %in% colnames(all.wkly2)) all.wkly2_column_names <- c(all.wkly2_column_names, 'Share_lo_99', 'Share_hi_99')

        # subset the columns to save ( the .. tells data.table that "all.wkly2_column_names" is an object in the global environment, not a column name in the data.table)
        all.wkly2 = all.wkly2[,..all.wkly2_column_names]

        # optionally calculate the number of infections attributable to each variant
        if (calc_confirmed_infections){
          test_filepath <- paste0(script.basename,
                                  "/data/backup_",
                                  data_date, custom_tag, "/",
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

            tests_wk <- rbind(tests_wk_us,
                              tests_wk_hhs)
            tests_wk[,'WEEK_END' := as.character(as.Date(yr_wk) + 6)]
            tests_wk[, 'yr_wk' := NULL]

            # merge the positive test results in with the variant proportion estimates
            all.wkly2 <- merge(
              x = all.wkly2,
              y = tests_wk,
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
                         ' not found. Not calculating number of infections attributable to each variant for weeks.'))
          }

        } # end calc confirmed infections

        # sort by HHS region
        all.wkly2 <- all.wkly2[order(all.wkly2$USA_or_HHSRegion),]

        # output data file depends on NO change in all.wkly2 column order!!!
        # If new columns are to be added, they need to go to the end.

        # make sure both WEEK_END columns are the same type
        all.wkly2[, 'WEEK_END' := as.character(WEEK_END)]
        unwt.wkly[, 'WEEK_END' := as.character(WEEK_END)]

        # merge on the unweighted estimates
        all.wkly3 <- merge(
          x = all.wkly2,
          y = unwt.wkly[,.(USA_or_HHSRegion, WEEK_END, Variant,
                           Unwt_share, Unwt_share_lo, Unwt_share_hi,
                           unwt_growth_rate, unwt_growth_rate_logodds,
                           unwt_doubling_time, unwt_doubling_time_logodds)],
          by = c('USA_or_HHSRegion', 'WEEK_END', 'Variant'),
          all.x = TRUE
        )



        # save the results to file
        write.csv(x = all.wkly3,
                  file = paste0(script.basename,
                                output_folder, "/variant_share_weekly_",
                                meth,
                                "_",
                                ci.type,
                                "CI_",
                                svy.type,
                                "_",
                                data_date,
                                tag,
                                "_",
                                results_tag,
                                ".csv"),
                  row.names=FALSE)

        if(meth == 'weighted'){
          # process dataframe and save to a format for direct hadoop upload
          all.wkly3_hadoop = data.frame(all.wkly3[,c(1:15)])
          all.wkly3_hadoop[is.na(all.wkly3_hadoop)] = '\\N'
          all.wkly3_hadoop[,16] = 'weighted'
          all.wkly3_hadoop[,17] = 'weekly'
          all.wkly3_hadoop[,18] = data_date
          all.wkly3_hadoop[,19] = paste0(results_tag, '_Run', opts$run_number)
          all.wkly3_hadoop[,20] = 1
          all.wkly3_hadoop[,21:24] = all.wkly3[,c("cases", "cases_hi", "cases_lo")]
          all.wkly3_hadoop[,3] = gsub('Delta Aggregated', 'B.1.617.2', all.wkly3_hadoop[,3])
          all.wkly3_hadoop[,3] = gsub('Omicron Aggregated', 'B.1.1.529', all.wkly3_hadoop[,3])
          all.wkly3_hadoop[,3] = gsub(' Aggregated', '', all.wkly3_hadoop[,3])
          write.table(x = all.wkly3_hadoop,
                      file = paste0(script.basename,
                                    output_folder, "/variant_share_weekly_weighted_",
                                    ci.type,
                                    "CI_",
                                    svy.type,
                                    "_",
                                    data_date,
                                    tag,
                                    "_",
                                    results_tag,
                                    "_hadoop.csv"),
                      quote = FALSE,
                      row.names = FALSE,
                      col.names = FALSE,
                      sep = ",")
        } # end save hadoop data
      } # end run 1 weekly estimates ("include_run1_weekly")
    } # end loop over weighted_method
  } # end !nowcast_only
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

    ## force-remove delta from the list of most abundant variants
    # (so that it is aggregated into "Other" and not automatically added into the
    #  Nowcast model as 1 of the "n_top" most abundant lineages)
    if(FALSE){
      us_var <- us_var[ names(us_var) != 'B.1.617.2' ]
    }

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

    # optionally make sure B and B.1 are not in the model_vars
    if (force_aggregate_B)   model_vars <- model_vars[model_vars %notin% 'B']
    if (force_aggregate_B.1) model_vars <- model_vars[model_vars %notin% 'B.1']
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
                      model_week %in% ((1:model_weeks)-model_week_mid)) # same as model_week_df$model_week
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

  # If a model_variant is only observed in the final model_week, then throw a warning
  # because there's a good chance that the model will predict that the variant will
  # just take over everything
  # variants only observed in the final model week:
  final_model_week_vars <- setdiff(unique(src.moddat$VARIANT), unique(src.moddat[model_week != max(src.moddat$model_week)][['VARIANT']]))
  if( length(final_model_week_vars) > 0 ){
    warning( paste0(
      'Variant(s) ',
      paste(final_model_week_vars, collapse = ', '),
      ' are only observed in the final model week. Results will likely be misleading!'
    ))
  }

  # optionally save src.moddat that's been prepped for analysis
  if(save_datasets_to_file){

    # consider saving a list of everything that is needed to fit the nowcast model: list(src.moddat, mysvy, model_vars)
    saveRDS(object = src.moddat,
            file = paste0(script.basename,
                          output_folder, '/src.moddat_', # save to results instead of 'data' folder
                          data_date,
                          tag,
                          '.RDS'))
  }


  # loop over weighting methods (weighted & unweighted only)
  for (meth in weighted_methods){
    ### survey design ----
    # weight using "weights", which is either SIMPLE_ADJ_WT or wts_trimmed
    mysvy = survey::svydesign(ids     = ~SOURCE,
                              strata  = ~STUSAB + yr_wk,
                              weights = ~wts,
                              nest = TRUE, # TRUE = disable checking for duplicate cluster ID's across strata
                              data = src.moddat) # not specifying fpc signifies sampling with replacement

    # for unweighted estimates, set all weights to 1
    if (meth == 'unweighted'){
      mysvy$prob <- rep(1, length(mysvy$prob))
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
    #write svymlm_hhs to file
    saveRDS(object = svymlm_hhs,
            file = paste0(script.basename,
                          output_folder, '/svymlm_hhs_', # save to results instead of 'data' folder
                          data_date,
                          tag,
                          '.RDS'))

    saveRDS(object = svymlm_us,
            file = paste0(script.basename,
                          output_folder, '/svymlm_us_', # save to results instead of 'data' folder
                          data_date,
                          tag,
                          '.RDS'))

    ## Plot results ----

    ### plot prep -----

    # choose which variants to display in the results
    if (display_option=="top7") {
      # most abundant 7 variants
      display_vars = head(model_vars, 7)
    } else {
      # all the vocs
      display_vars = head(model_vars,15)
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
                  output_folder, "/wtd_shares_",
                  data_date,
                  "_")


    # Weighted variant shares of the top variants in the past `display_lookback`
    # weeks (number of sequences collected weekly above each bar), and model-based
    # smoothed estimates, nationwide:

    # ### barplot of weighted share (national) ----

    # # tabulate/count the normalized survey weights by week, and variant
    # # (only used for barplot below)
    # bp_us = xtabs(formula = NORM_WTS ~ model_week + VARIANT,
    #               data = src.dat,
    #               subset = src.dat$model_week %in% ((data_week_df$model_week - 2 - display_lookback + 1):(data_week_df$model_week - 2))) #  ((model_week_mid - display_lookback + 1):model_week_mid))
    # # DOUBLE-CHECK THAT THIS ALWAYS WORKS!

    # # Normalize values to proportions
    # bp_us = prop.table(bp_us, 1)

    # # subset the weights to only include the variants to be displayed
    # bp_us = bp_us[, display_vars]

    # # give the table prettier column names
    # # rownames(bp_us) = week_label(as.numeric(rownames(bp_us)) - current_week)
    # rownames(bp_us) = format(
    #   x = ((as.numeric(rownames(bp_us)) + data_week_df$model_week - 2) + model_week_min) * 7 + as.Date(week0day1),
    #   # x = ((as.numeric(rownames(bp_us)) + model_week_mid) + model_week_min) * 7 + as.Date(week0day1), # DOUBLE-CHECK THAT THIS IS ALWAYS GIVING THE CORRECT WEEKS
    #   format = '%m-%d'
    # )

    # # create plot
    # if (fig_gen_run) jpeg(filename  = paste0(stub, "barplot_US", tag, ".jpg"),
    #                       width     = 1500,
    #                       height    = 1500,
    #                       pointsize = 40)
    # # create a barplot of "observed" values (i.e. weighted counts)
    # bp = barplot(height = 100 * t(bp_us),
    #              xlab = "Week beginning",
    #              ylab = "Weighted variant share (%)",
    #              main = "Nationwide",
    #              border = NA,
    #              ylim = 110 * 0:1,
    #              col = col.dk,
    #              names.arg = rownames(bp_us),
    #              legend.text = display_vars,
    #              args.legend = list(x = "topleft",
    #                                 bty = "n",
    #                                 border = NA))

    # # add text to the barplot
    # text(x = bp,
    #      y = 3 + colSums(100 * t(tail(bp_us, 12))),
    #      labels = with(subset(src.dat,
    #                           week < current_week - 1 &
    #                             week >= current_week - display_lookback),
    #                    table(week)),
    #      cex = 0.7)
    # if (fig_gen_run) dev.off()


    # ### barplot of model-predicted data (national) ----
    # # note: the x-axis labels are at the week midpoint, but labeled with week starting date

    # # create a dataframe for predictions for each week
    # pred_us.df = expand.grid(model_week = seq(from = -display_lookback,
    #                                           to = 2, # go 2 weeks into the future
    #                                           by = (1/7)) + data_week_df$model_week)

    # #add a column for (predicted) each variant proportion for each timepoint
    # pred_us.df = cbind(pred_us.df,
    #                    predict(object = svymlm_us$mlm,
    #                            newdata = pred_us.df,
    #                            type = "probs"))

    # # add in dates
    # # pred_us.df$date       = as.Date((pred_us.df$model_week + model_week_mid + model_week_min) * 7 + week0day1)
    # # pred_us.df$week_start = pred_us.df$date - as.numeric(format(pred_us.df$date, format = '%w'))
    # # pred_us.df$week_end   = pred_us.df$week_start + 6
    # # pred_us.df$week_mid   = pred_us.df$week_start + 3

    # # Barplot of Nowcast predicted values
    # if (fig_gen_run) jpeg(filename  = paste0(stub, "projection_US", tag, ".jpg"),
    #                       width     = 1500,
    #                       height    = 1500,
    #                       pointsize = 40)

    # bp = barplot(height = 100 * t(pred_us.df[, 1 + display_indices]),
    #              xlab = "Week beginning",
    #              ylab = "Weighted variant share (%)",
    #              main = "Nationwide",
    #              space = 0,
    #              border = NA,
    #              ylim = 110 * 0:1,
    #              col = col.dk,
    #              names.arg = ifelse(test = pred_us.df$model_week %% 1 == 0,
    #                                 yes = format(((pred_us.df$model_week + (data_week_df$model_week - 2)) + model_week_min) * 7 + as.Date(week0day1), format = '%m-%d'),
    #                                 no = NA),
    #              legend.text = display_vars,
    #              args.legend = list(x = "topleft",
    #                                 bty = "n",
    #                                 border = NA))

    # # predicted percent contributions of some variants for the current week
    # # pc = unlist(100 * subset(pred_us.df,
    # #                          week==current_week)[, 1+display_indices])
    # pc = unlist(100 * subset(pred_us.df,
    #                          model_week == data_week_df$model_week)[, 1+display_indices])

    # # define a y-value for the text to be added to the plot
    # y = cumsum(pc) - pc/2

    # # define x-values for grey boxes signifying that recent data is likely incomplete
    # # x = bp[which(pred_us.df$week %in% (current_week + c(-2, 0, 2)))]
    # x <- bp[which(pred_us.df$model_week %in% (data_week_df$model_week + c(-2, 0, 2)))]

    # # add text to the plot
    # text(x = x[2], # 1.02 * tail(bp, 1),
    #      y = y,
    #      labels = round(pc, 1),
    #      cex = 0.7,
    #      xpd = TRUE,
    #      adj = c(0.5, 0.5))

    # # add grey rectangles to plot
    # rect(xleft   = x[1] + (x[2]-x[1])*.25, # shift it over by 1/2 week so that it is centered on individual weeks
    #      ybottom = 0,
    #      xright  = x[2] + (x[3]-x[2])*.25,
    #      ytop    = 100,
    #      border  = NA,
    #      col = "#00000020")
    # # abline(v = x[2],
    # #        lty = 2, color = 'grey20')

    # # add text to the plot
    # text(x = mean(c(x[1], x[2])) + (x[2]-x[1])*.25,
    #      y = 101,
    #      labels = 'Nowcast',
    #      cex = 0.7,
    #      xpd = TRUE,
    #      adj = c(0.5, 0))
    # rect(xleft   = x[2] + (x[3]-x[2])*.25,
    #      ybottom = 0,
    #      xright  = x[3],
    #      ytop    = 100,
    #      border  = NA,
    #      col = "#00000040")
    # # add text to the plot
    # text(x = mean(c(x[2], x[3]))+ (x[2]-x[1])*.125,
    #      y = 101,
    #      labels = 'Future',
    #      cex = 0.7,
    #      xpd = TRUE,
    #      adj = c(0.5, 0))
    # if (fig_gen_run) dev.off()

    ### WoW growth rate vs. transmission ----
    #  - the vertical axis depicts the variant share growth rate
    #    (derivative of log of variant proportion with respect to time).
    #  - The axis on the right shows an estimate of variant transmissibility with
    #    respect to the overall mean transmissibility.

    # Get the SE of the national estimate
    us.summary = se.multinom(mlm   = svymlm_us$mlm,
                             newdata_1row = data.frame(
                               model_week = data_week_df$model_week, # formerly "model_week_mid", but that doesn't always work
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

    # create a dataframe of variant shares & growth rates
    gr_tab = data.frame(variant          = c(model_vars, "OTHER"),
                        variant_share    = (100 * us.summary$p_i),
                        variant_share_lo = 100 * (us.summary$p_i - 1.96 * us.summary$se.p_i),
                        variant_share_hi = 100 * (us.summary$p_i + 1.96 * us.summary$se.p_i),
                        growth_rate      = gr,
                        growth_rate_lo   = gr_lo,
                        growth_rate_hi   = gr_hi,
                        doubling_time    = doubling_time,
                        doubling_time_lo = doubling_time_lo,
                        doubling_time_hi = doubling_time_hi,
                        model_week       = data_week_df$model_week)

    # merge in date
    gr_tab <- merge(gr_tab,
                    data_week_df,
                    by = 'model_week')

    # save growth rates to file
    write.csv(x = gr_tab,
              file = paste0(script.basename,
                            output_folder, "/wow_growth_variant_share_",
                            meth,
                            "_",
                            data_date,
                            tag,
                            ".csv"),
              row.names = FALSE)


    if (fig_gen_run) png(filename = paste0(stub, "growthrate_US_", meth, tag, ".png"),
                         width = 8,
                         height = 8,
                         units = "in",
                         pointsize = 16,
                         res = 1000)

    # make enough room on the right side of the plot for a secondary axis.
    orpar <- par()
    par(mar = c(5.1, 4.1, 4.1, 4.1))

    # filter out extremely rare variants from the plot
    if (force_aggregate_omicron && custom_lineages == FALSE) {
      gtp <- subset(gr_tab,
                    variant %in% c(voc1, "OTHER"))
    } else {
      gtp <- subset(gr_tab,
                    variant_share >= 0.01) # this is already a percent, so this is filtering out variants with less than 1/1000 of a percent (not 1 percent)
    }

    if (custom_lineages == TRUE) {
      gtp$variant <- gsub("R346T_(B[AF])([1-9])([1-9]{2})([1-9]{1})", "\\1.\\2.\\3.\\4+R346T", gtp$variant)
      gtp$variant <- gsub("R346T_(B[AFEQ])([1-9])([1-9]{1,2})", "\\1.\\2.\\3+R346T", gtp$variant)
      gtp$variant <- gsub("R346T_(B[AFQ])([1-9])", "\\1.\\2+R346T", gtp$variant)
      gtp$variant <- gsub("R346T_B11529", "B.1.1.529+R346T", gtp$variant)
    }

    wow_x_scale <- floor(log(min(gtp$variant_share),10))
    wow_x_min <- 5 * (10 ^ (wow_x_scale - 1))

    # plot "nowcast" groth rates by variant
    plot(x = gtp$variant_share,
         y = gtp$growth_rate,
         log = "x",
         type = "n",
         ylim = range(gtp$growth_rate_lo, gtp$growth_rate_hi),
         xaxt = "n",
         xlim = c(wow_x_min,110),
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
         at     = 10 ^ seq(wow_x_scale,2,by=1),
         labels = 10 ^ seq(wow_x_scale,2,by=1))

    # add a horizontal line at 0
    abline(h = 0,
           col = "grey65")

    # add lines for each variant
    for (vv in 1:nrow(gtp)) {

      # vertical lines for uncertainty in growth rate
      lines(x = rep(gtp$variant_share[vv], 2),
            y = gtp[vv,c('growth_rate_lo', 'growth_rate_hi')],
            col = "blue",
            lwd = 2)

      # horizontal lines for uncertainty in variant proportion
      lines(x = pmax(0.0001, gtp[vv, c('variant_share_lo', 'variant_share_hi')]),
            y = rep(gtp$growth_rate[vv], 2),
            col = "blue",
            lwd = 2)
    }

    # add the name of each variant
    text(x = gtp$variant_share,
         y = gtp$growth_rate,
         labels = gtp$variant,
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


    ### regional barplots (weighted) ----

    # Same as above for each HHS region:
    # Weighted variant shares of the top variants in the past `display_lookback`
    # weeks (number of sequences collected weekly above each bar), and model-based
    # smoothed estimates, for each HHS region:


    # tabulate/count the normalized survey weights by region, week, and variant
    # (just used for barplots)
    bp_hhs = xtabs(formula = NORM_WTS ~ HHS + model_week + VARIANT,
                   data = src.dat,
                   subset = src.dat$model_week %in% ((data_week_df$model_week - 2 - display_lookback + 1):(data_week_df$model_week - 2))) # ((model_week_mid - display_lookback + 1):model_week_mid))
    # DOUBLE-CHECK THIS USE OF MODEL_WEEK_MID

    # subset the weights to only include the variants to be displayed
    bp_hhs = prop.table(bp_hhs, 1:2)[,, display_vars]

    # give the table prettier column names (week start date instead of model_week)
    # dimnames(bp_hhs)[[2]] = week_label(as.numeric(dimnames(bp_hhs)[[2]]) - current_week)
    dimnames(bp_hhs)[[2]] = format(
      x = ((as.numeric(dimnames(bp_hhs)[[2]]) + data_week_df$model_week - 2) + model_week_min) * 7 + as.Date(week0day1),
      # x = ((as.numeric(dimnames(bp_hhs)[[2]]) + model_week_mid) + model_week_min) * 7 + as.Date(week0day1),
      format = '%m-%d'
    )


    # add region to the growth rate table (and then add on the HHS regional growth rates)
    gr_tab$region = 'USA'

    for (hhs in sort(unique(src.dat$HHS))) {

      # Barplot of regional data
      if (fig_gen_run) jpeg(filename = paste0(stub, "barplot_HHS", hhs, "_", meth, tag,".jpg"),
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
                                HHS == hhs &
                                  week < current_week - 1 &
                                  week >= current_week - display_lookback),
                         table(week)),
           cex = 0.7)
      if (fig_gen_run) dev.off()

      ### regional barplot (Nowcast) ----

      # create a dataframe for predictions for each week & HHS region
      pred_hhs.df = expand.grid(model_week = seq(from = -display_lookback,
                                                 to = 2, # go 2 weeks into the future
                                                 by = (1/7)) + data_week_df$model_week, # (model_weeks - model_week_mid), #  + current_week,
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


      if (fig_gen_run) jpeg(filename = paste0(stub, "projection_HHS", hhs, "_", meth, tag, ".jpg"),
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
                                      format(((pred.df$model_week + (data_week_df$model_week - 2)) + model_week_min) * 7 + as.Date(week0day1), format = '%m-%d'),
                                      NA),
                   legend.text = display_vars,
                   args.legend = list(x = "topleft",
                                      bty = "n",
                                      border = NA))

      # predicted percent contributions of some variants for the current week
      # pc = unlist(100 * subset(pred.df, week==current_week)[, 1+display_indices])
      pc = unlist(100 * subset(pred.df, model_week == data_week_df$model_week)[, 1+display_indices])

      # define a y-value for the text on the plot
      y = cumsum(pc) - pc/2

      # define x-values for grey boxes signifying that recent data is likely incomplete
      # x = bp[which(pred.df$week %in% (current_week + c(-2, 0, 2)))]
      x = bp[which(pred.df$model_week %in% (data_week_df$model_week + c(-2, 0, 2)))]

      # add text to the plot
      text(x = x[2], # 1.02 * tail(bp, 1),
           y = y,
           labels = round(pc, 1),
           cex = 0.7,
           xpd = TRUE,
           adj = c(0.5, 0.5))

      # add boxes/rectangles to the plot
      rect(xleft   = x[1] + (x[2]-x[1])*.25, # shift it over by 1/2 week so that it is centered on individual weeks
           ybottom = 0,
           xright  = x[2] + (x[3]-x[2])*.25,
           ytop = 100,
           border = NA,
           col = "#00000020")
      # add text to the plot
      text(x = mean(c(x[1], x[2])) + (x[2]-x[1])*.25,
           y = 101,
           labels = 'Nowcast',
           cex = 0.7,
           xpd = TRUE,
           adj = c(0.5, 0))
      rect(xleft   = x[2] + (x[3]-x[2])*.25,
           ybottom = 0,
           xright = x[3],
           ytop = 100,
           border = NA,
           col = "#00000040")
      # add text to the plot
      text(x = mean(c(x[2], x[3]))+ (x[2]-x[1])*.125,
           y = 101,
           labels = 'Future',
           cex = 0.7,
           xpd = TRUE,
           adj = c(0.5, 0))
      if (fig_gen_run) dev.off()


      ### regional WoW growth rate -----

      # get the SE for the given model, week, and HHS region
      hhs.summary = se.multinom(mlm = svymlm_hhs$mlm,
                                newdata_1row = data.frame(
                                  model_week = data_week_df$model_week, # model_week_mid,
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
                              variant_share_lo = 100 * (hhs.summary$p_i - 1.96 * hhs.summary$se.p_i),
                              variant_share_hi = 100 * (hhs.summary$p_i + 1.96 * hhs.summary$se.p_i),
                              growth_rate      = gr,
                              growth_rate_lo   = gr_lo,
                              growth_rate_hi   = gr_hi,
                              doubling_time    = doubling_time,
                              doubling_time_lo = doubling_time_lo,
                              doubling_time_hi = doubling_time_hi,
                              region           = hhs,
                              model_week       = data_week_df$model_week
      )

      # merge in date
      gr_tab_hhs <- merge(gr_tab_hhs,
                          data_week_df,
                          by = 'model_week')

      # combine national and regional growth rates
      gr_tab = rbind(
        gr_tab,
        gr_tab_hhs
      )

      # create plot of growth rates by region
      if (fig_gen_run) jpeg(filename = paste0(stub, "growthrate_HHS", hhs, "_", meth, tag,".jpg"),
                            width = 1500,
                            height = 1500,
                            pointsize = 40)

      # make enough room on the right side of the plot for a secondary axis.
      orpar <- par()
      par(mar = c(5.1, 4.1, 4.1, 4.1))

      # filter out variants with proportions less than 0.01%
      if (force_aggregate_omicron && custom_lineages == FALSE) {
        gtphhs <- subset(gr_tab_hhs,
                         variant %in% c(voc1, "OTHER"))
      } else {
        gtphhs <- subset(gr_tab_hhs,
                         variant_share >= 0.01) # this is already a percent, so this is filtering out variants with less than 1/1000 of a percent (not 1 percent)
      }

      if (custom_lineages == TRUE) {
        gtphhs$variant <- gsub("R346T_(B[AF])([1-9])([1-9]{2})([1-9]{1})", "\\1.\\2.\\3.\\4+R346T", gtphhs$variant)
        gtphhs$variant <- gsub("R346T_(B[AFEQ])([1-9])([1-9]{1,2})", "\\1.\\2.\\3+R346T", gtphhs$variant)
        gtphhs$variant <- gsub("R346T_(B[AFQ])([1-9])", "\\1.\\2+R346T", gtphhs$variant)
        gtphhs$variant <- gsub("R346T_B11529", "B.1.1.529+R346T", gtphhs$variant)
      }


      wow_x_scale <- floor(log(min(gtphhs$variant_share),10))
      wow_x_min <- 5 * (10 ^ (wow_x_scale - 1))


      # plot "nowcast" growth rates by variant
      plot(x = 100 * gtphhs$variant_share,
           y = gtphhs$growth_rate,
           log = "x",
           type = "n",
           ylim = range(gtphhs$growth_rate_lo, gtphhs$growth_rate_hi),
           xaxt = "n",
           xlim = c(wow_x_min,110),
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
           at     = 10 ^ seq(wow_x_scale,2,by=1),
           labels = 10 ^ seq(wow_x_scale,2,by=1))

      # add a horizontal line at 0
      abline(h = 0,
             col = "grey65")

      # add lines for each variant
      for (vv in 1:nrow(gtphhs)) {
        # vertical lines for uncertainty in growth rate
        lines(x = rep(gtphhs$variant_share[vv], 2),
              y = gtphhs[vv, c('growth_rate_lo', 'growth_rate_hi')],
              col = "blue",
              lwd = 2)

        # horizontal lines for uncertainty in variant proportions
        lines(x = gtphhs[vv, c("variant_share_lo", "variant_share_hi")],
              y = rep(gtphhs$growth_rate[vv], 2),
              col = "blue",
              lwd = 2)
      }


      # identify variants with share > 1% or growth rate > 0%
      labels = unique(c(which(100 * gtphhs$variant_share>1),
                        which(gr > 1)))

      # label variants
      text(x = gtphhs$variant_share,
           y = gtphhs$growth_rate,
           labels = gtphhs$variant,
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
                            output_folder, "/wow_growth_variant_share_",
                            meth,
                            "_",
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
    # this is used to aggregate the results of run 2 to the vocs included in run 1.
    # The aggregation matrix has a column for each of the vocs in run 2 and a row for each aggregated voc.
    # 1's in the matrix specify which sub-variants to aggregate into an aggregate lineage.
    if( tolower(opts$voc_aggregation_method) %in% c('updated', 'lineage_expanded') ){

      # create a column for each voc in the model (plus "OTHER")
      agg_var_mat <- matrix(data = 0,
                            nrow = 0,
                            ncol = (length(model_vars)+1))
      colnames(agg_var_mat) <- c(model_vars,"Other")

      # get the lineage_expanded for the model_vars
      model_vars_expanded <- setNames(voc_lut$lineage_expanded, voc_lut$variant)[ model_vars ]

      # potential parent variants are voc1_expanded
      voc1_expanded <- setNames(voc_lut$lineage_expanded, voc_lut$variant)[ voc1 ]

      # get the parent varients for the model_vars_expanded from the voc1_expanded
      model_var_parents_expanded <- nearest_parent( model_vars_expanded,  voc1_expanded )
      # get the short names for the parent variants of model_vars
      model_var_parents <- setNames( voc_lut$variant, voc_lut$lineage_expanded )[ model_var_parents_expanded ]

      # model_var look-up table
      model_var_lut <- data.frame(
        voc2_model_vars = model_vars,
        voc1_model_vars = model_var_parents
      )

      # only need to convert rows where voc2_model_vars != voc1_model_vars
      mvl_sub <- model_var_lut[ model_var_lut$voc2_model_vars != model_var_lut$voc1_model_vars, ]

      # get the names of the parent variants that will have sub-lineages aggregated into them
      unique_mvl_sub <- unique(mvl_sub$voc1_model_vars)

      # for each unique mvl_sub$voc1_model_vars, fill in values in the matrix for the variants that will be aggregated into it.
      for (i in unique_mvl_sub) {

        # get the voc2 that are aggregated into the given voc1
        # model_var_lut[ model_var_lut$voc1_model_vars == i, 'voc2_model_vars']

        # define an extra row to add onto the agg_var_mat
        extra_row <- ifelse( colnames(agg_var_mat) %in%
                               # the variant itself & all the sub-lineages to be aggregated into it
                               c(i, model_var_lut[ model_var_lut$voc1_model_vars == i, 'voc2_model_vars']) ,
                             1,
                             0)

        # add the new row onto the aggregation matrix
        agg_var_mat <- rbind(
          agg_var_mat,
          extra_row
        )
        row.names(agg_var_mat)[nrow(agg_var_mat)] <- paste(i, 'Aggregated')
      } # end loop over unique_mvl_sub

      # all the variants (not in voc1) AND (not aggregated into something else)
      # should be aggregated into "Other Aggregated"
      other_agg <- base::setdiff(colnames(agg_var_mat)[colSums(agg_var_mat) == 0],
                                 voc1)
      # add the new row onto the aggregation matrix
      agg_var_mat <- rbind(
        agg_var_mat,
        ifelse(colnames(agg_var_mat) %in% other_agg, 1, 0)
      )
      # update the row name
      row.names(agg_var_mat)[nrow(agg_var_mat)] <- 'Other Aggregated'

      # double-check that no variant is aggregated into multiple vocs
      if(max(colSums(agg_var_mat)) > 1){
        warning(message = paste(
          'Aggregated results are invalid! These variants are being aggregated multiple times:',
          names(agg_var_mat)[colSums(agg_var_mat) > 1],
          '. Fix the aggregation matrix.'))
      }

      # NOTE! Setting up the agg_var_mat in this way assumes that variants included in
      # voc2 but NOT in voc1 (i.e: setdiff(voc, voc1) ) will be aggregated into something
      # that IS in voc1 (other than "Other Aggregated"). This is necessary to ensure that the
      # "Other Aggregated" row is the SAME for run1 output and run2 output, which
      # is important because they both use that row for output below.
      # this is a check to make sure that the assumption stated above is being met.
      sub_mat <- agg_var_mat[row.names(agg_var_mat) != 'Other Aggregated',setdiff(voc, voc1), drop = FALSE]
      if(!all(colSums(sub_mat) > 0)){
        problem_vocs <- colnames(sub_mat)[colSums(sub_mat) == 0]
        warning(message = paste0('Not all variants in "voc" are aggregated into something in "voc1". ',
                                 paste(problem_vocs, collapse = ', '),
                                 ' is in voc2 but is not aggregated into anything in voc1. Therefore it will be aggregated into "Other Aggregated" for both run_1 and run_2 output. Either change voc1 and/or voc2 or change the way agg_var_mat works.'))
      }


      # save the aggregation matrix to file for double-checking
      write.csv(
        x = replace(agg_var_mat, agg_var_mat == 0, NA), # it's easier to view in Excel with only the 1's and no 0's
        file = paste0(script.basename, output_folder, "/agg_var_mat_",
                      ci.type,
                      "CI_",
                      svy.type,
                      "_", # meth, '_', # weighted & unweighted agg_var_mat will be identical, so no need to save both
                      data_date,
                      tag,
                      ".csv"),
        row.names = T,
        na = ""
      )
    } else {
      # use the former aggregation

      # Check to see which lineages are in model_vars
      # this returns all variants with "AY" in the name
      AY_vars = model_vars[grep("^AY\\.", model_vars, perl=T)]
      # this returns all variants with BA. in the name (Omicron sublineages)
      BA_vars = model_vars[grep("(^B[AC-HJ-NP-VYZ]\\.)|(^C[A-HJ-NP-WYZ]\\.)|(^D[A-HJ-NP-WYZ]\\.)|(^E[A-FHJNPQRSTVWYZ]\\.)|(^F[ABCFJKMN]\\.)", model_vars, perl=T)]
      # this returns all variants with XBB. in the name
      XBB_vars = model_vars[grep("(^XBB\\.)|^E[GKLMU]\\.|^F[DEGHL]\\.", model_vars, perl=T)]

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
      # Note: this assumes that all BA subvariants will be aggregated into parent Omicron. That's no longer the case, so below I'm adding code to deal with BA subvariants that need their own aggregations.
      BA_agg = BA_vars[BA_vars %notin% run1_lineages]
      # some of these should be aggregated into parent omicron; others should be aggregated into some BA sublineage
      # this returns all "XBB" variants to be aggregated for Run 1 results
      XBB_agg = XBB_vars[XBB_vars %notin% run1_lineages]
      # all other variants to be aggregated (used for Run 1 & Run 2 results)
      Other_agg = model_vars[model_vars %notin% c(voc, 'B.1.617.2', 'B.1.1.529')] # also have to exclude any variants that have their own row in agg_var_mat. If delta is not included in VOC, then it gets included in both the delta row and the Other row.

      # generate a matrix that indicates which lineages to aggregate for the nowcast
      # Columns B529.BF.7.4are the lineages in the nowcast model, so all the defined lineages
      #  plus the "other" lineage
      # Rows are the aggregated lineages desired

      # Check whether any of the large parent lineage aggregation has zero sublineage in the model variant list
      agg_var_mat <- matrix(data = 0,
                            nrow = 1,
                            ncol = (length(model_vars)+1))
      colnames(agg_var_mat) <- c(model_vars,"Other")
      agg_var_mat[1,] <- ifelse(colnames(agg_var_mat) %in% c(Other_agg, "Other"),1,0)
      row.names(agg_var_mat) <- c("Other Aggregated")

      if(length(AY_agg) != 0) {
        extra_row <- ifelse(colnames(agg_var_mat) %in% c("B.1.617.2", AY_agg),1,0)
        agg_var_mat <- rbind(agg_var_mat, extra_row)
        row.names(agg_var_mat)[nrow(agg_var_mat)] <- c("Delta Aggregated")
      }
      if(length(BA_agg) != 0) {
        extra_row <- ifelse(colnames(agg_var_mat) %in% c("B.1.1.529", BA_agg),1,0)
        agg_var_mat <- rbind(agg_var_mat, extra_row)
        row.names(agg_var_mat)[nrow(agg_var_mat)] <- c("Omicron Aggregated")
      }
      if(length(XBB_agg) != 0) {
        extra_row <- ifelse(colnames(agg_var_mat) %in% c("XBB", XBB_agg),1,0)
        agg_var_mat <- rbind(agg_var_mat, extra_row)
        row.names(agg_var_mat)[nrow(agg_var_mat)] <- c("XBB Aggregated")
      }

      # agg_var_mat <- matrix(data = 0,
      #                       nrow = 4,
      #                       ncol = (length(model_vars)+1))
      # colnames(agg_var_mat) <- c(model_vars,"Other")

      # # Fill in matrix values: if lineage is to be aggregated to parent lineage in
      # # given row, then value = 1, else value = 0
      # agg_var_mat[1,] <- ifelse(colnames(agg_var_mat) %in% c("B.1.617.2", AY_agg),1,0)
      # agg_var_mat[2,] <- ifelse(colnames(agg_var_mat) %in% c("B.1.1.529", BA_agg),1,0)
      # agg_var_mat[3,] <- ifelse(colnames(agg_var_mat) %in% c(Other_agg, "Other"),1,0)
      # agg_var_mat[4,] <- ifelse(colnames(agg_var_mat) %in% c("XBB", XBB_agg),1,0)
      # row.names(agg_var_mat) <-c(#"Delta Aggregated",
      #                             "Omicron Aggregated", "Other Aggregated", "XBB Aggregated")

      # add rows to agg_var_mat for each XBB subvariant in run1_lingeages
      {
        XBBs_in_r1l <- XBB_vars[XBB_vars %in% run1_lineages]
        for (ll in XBBs_in_r1l){
          if(exists('ll_agg')) rm(ll_agg)
          if(ll == 'XBB.1.5'){
            ll_agg <- grep("^XBB\\.1\\.5(?![0-9])|^FD\\.", XBB_vars, perl = T, value = T)
            ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
            ll_agg <- c(ll_agg, 'XBB.1.5')
          }
          if(ll == 'XBB.1.16'){
            ll_agg <- grep("^XBB\\.1\\.16(?![0-9])", XBB_vars, perl = T, value = T)
            ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
            ll_agg <- c(ll_agg, 'XBB.1.16')
          }
          if(ll == 'XBB.1.9.2'){
            ll_agg <- grep("^XBB\\.1\\.9\\.2(?![0-9])|^EG\\.", XBB_vars, perl = T, value = T)
            ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
            ll_agg <- c(ll_agg, 'XBB.1.9.2')
          }
          if(ll == 'XBB.2.3'){
            ll_agg <- grep("^XBB\\.2\\.3(?![0-9])", XBB_vars, perl = T, value = T)
            ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
            ll_agg <- c(ll_agg, 'XBB.2.3')

          }
          if(exists('ll_agg')){
            # if ll_agg contains subvariants of ll, then add a new row to agg_var_mat
            # if ll_agg only contains ll, don't add a new row
            if(any(ll_agg != ll)){
              # create an extra row for the agg_var_mat
              extra_row <- ifelse(colnames(agg_var_mat) %in% ll_agg,1,0)

              # add the new row onto the aggregation matrix
              agg_var_mat <- rbind(
                agg_var_mat,
                extra_row
              )
              row.names(agg_var_mat)[nrow(agg_var_mat)] <- paste(ll, 'Aggregated')

              # remove aggregation indices from "XBB Aggregated" if they're in the new row
              om_row <- which(row.names(agg_var_mat) == 'XBB Aggregated')
              agg_var_mat[om_row,which(agg_var_mat[ om_row,] == 1 & extra_row == 1)] <- 0
            }
          }
        }
      }


      # add rows to agg_var_mat for each BA subvariant in run1_lineages
      {
        # get the BA subvariants that are in run1_lineages, and therefore might need their own aggregations
        #  (e.g. aggregate BA.2.12 into BA.2 if both are in voc2, but only BA.2 is in voc1)
        BAs_in_r1l <- BA_vars[BA_vars %in% run1_lineages]
        for( ll in BAs_in_r1l){
          if(exists('ll_agg')) rm(ll_agg)
          # first get the BA_vars that should be aggregated into this "ll" instead of parent omicron
          # grep(pattern = ll, x = BA_vars) # this won't work in all cases. I'll have to use a piece-wise/messy "solution"
          if (ll == 'BA.1'){
            # this always excludes BA.1.1
            # if voc1 includes BA.1, but not BA.1.1, then this is going to result in BA.1.1 aggregated into parent omicron, NOT BA.1
            ll_agg <- grep("(^BA\\.1)(?!(\\.1$))",BA_vars, perl = T, value = T)
          }
          if (ll == 'BA.1.1'){
            ll_agg <- grep("(^BA\\.1\\.1)(?![0-9])",BA_vars, perl = T, value = T)
          }
          if (ll == 'BA.2') {
            ll_agg <- grep("(^BA\\.2)(?![0-9])",BA_vars, perl = T, value = T)
            ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
            ll_agg <- c(ll_agg, 'BA.2')
            # # this will always aggregate BA.2.12 into BA.2
            # if( length(grep("(^BA\\.2)(?![0-9])",run1_lineages, perl = T, value = T)) == 1 ){
            #   ll_agg <- grep("(^BA\\.2)(?![0-9])",BA_vars, perl = T, value = T)
            # }
            # # this will keep BA.2.12.1 seperate ONLY if BA.2.12.1 is *ALSO* listed in run1_lineages
            # if( length(grep("(^BA\\.2\\.12\\.1)",run1_lineages, perl = T, value = T)) == 1 ){
            #   ll_agg <- grep("(^BA\\.2)(?!([0-9])|(\\.12\\.1))",BA_vars, perl = T, value = T)
            # }
            # # this will keep BA.2.75 and BA.2.75.* seperate ONLY if BA.2.75 is *ALSO* listed in run1_lineages
            # if( length(grep("(^BA\\.2\\.75)",run1_lineages, perl = T, value = T)) >= 1 ){
            #   ll_agg <- grep("(^BA\\.2)(?!([0-9])|(\\.75))",BA_vars, perl = T, value = T)
            # }
            # # this will keep ba.2.12.1 AND BA.2.75 and BA.2.75.* seperate ONLY if BOTH BA.2.12.1 and BA.2.75 are both *ALSO* listed in run1_lineages
            # if( length(grep("(^BA\\.2\\.75)",run1_lineages, perl = T, value = T)) >= 1 && length(grep("(BA\\.2\\.12\\.1)",run1_lineages, perl = T, value = T)) == 1 ){
            #   if('BN.1' %in% run1_lineages) {
            #     ll_agg <- grep("(^BA\\.2)(?!([0-9])|(\\.75)|(\\.12\\.1))",BA_vars, perl = T, value = T)
            #   } else {
            #     ll_agg <- grep("(^BA\\.2)(?!([0-9])|(\\.75)|(\\.12\\.1))|(^BN\\.1)(?![0-9])",BA_vars, perl = T, value = T)
            #     }
            # }
          }
          if(ll == 'BA.2.12.1') {
            ll_agg <- grep("(^BA\\.2\\.12\\.1)|(^BG\\.)",BA_vars, perl = T, value = T)
          }

          if(ll == 'BA.2.75') {
            ll_agg <- grep("(^BA\\.2\\.75)(?![0-9])|(^BN\\.)|(^CH\\.)",BA_vars, perl = T, value = T)
            if(length(grep("(^CH\\.1\\.1)",run1_lineages, perl = T, value = T)) == 1) {
              ll_agg_not <- grep("(^CH\\.1\\.1)(?![0-9])",ll_agg, perl = T, value =T)
              ll_agg <- setdiff(ll_agg, ll_agg_not)
            }
            ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
            ll_agg <- c(ll_agg, 'BA.2.75')
          }
          if(ll == 'CH.1.1'){
            ll_agg <- grep("(^CH\\.1\\.1)(?![0-9])",BA_vars, perl = T, value =T)
            ll_agg <- c(ll_agg, 'CH.1.1')
          }
          if (ll == 'BA.3') {
            ll_agg <- grep("(^BA\\.3)(?![0-9])",BA_vars, perl = T, value = T)
          }
          if (ll == 'BA.4') {
            ll_agg <- grep("(^BA\\.4)(?![0-9])",BA_vars, perl = T, value = T)
            # this will keep BA.4.6 seperate ONLY if BA.4.6 is *ALSO* listed in run1_lineages
            if( length(grep("(^BA\\.4\\.6)",run1_lineages, perl = T, value = T)) == 1 ){
              ll_agg <- grep("(^BA\\.4)(?!([0-9])|(\\.6))",BA_vars, perl = T, value = T)
            }
          }
          if (ll == 'BA.4.6'){
            ll_agg <- grep("(^BA\\.4\\.6)(?![0-9])",BA_vars, perl = T, value = T)
          }
          if(ll == 'BQ.1.1') {
            ll_agg <- grep("(^BQ\\.1\\.1)(?![0-9])|(^DK\\.)",BA_vars, perl = T, value = T)
          }
          if(ll == 'BQ.1') {
            if(length(grep("^BQ\\.1\\.1(?![0-9])", run1_lineages, perl=T, value=T)) > 0){
              ll_agg <- grep("(^BQ\\.1)(?![0-9]|(\\.1(?![0-9])))",BA_vars, perl = T, value = T)
            } else {
              ll_agg <- grep("(^BQ\\.1)(?![0-9])",BA_vars, perl = T, value = T)
            }
            ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
            ll_agg <- c(ll_agg, 'BQ.1')
          }
          if(ll == 'BF.11') {
            ll_agg <- grep("(^BF\\.11)(?![0-9])",BA_vars, perl = T, value = T)
          }
          if(ll == 'BF.7') {
            ll_agg <- grep("(^BF\\.7)(?![0-9])",BA_vars, perl = T, value = T)
          }
          if(ll == 'BA.5') {
            if( length(grep("(^BF\\.7)(?![0-9])",run1_lineages, perl = T, value = T)) > 0 ){
              ll_agg <- grep("(^CQ\\.)|(^BA\\.5)(?![0-9])|(^BE\\.)|(^BF\\.(?!7(?![0-9])))|(^C[KR]\\.)|(^DF\\.)",BA_vars, perl = T, value = T)
            } else {
              ll_agg <- grep("(^CQ\\.)|(^BA\\.5)(?![0-9])|(^B[EF]\\.)|(^C[KR]\\.)|(^DF\\.)",BA_vars, perl = T, value = T)
            }
            # this will keep the sublineage seperate if it is ALSO listed in run1 lineages
            ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
            ll_agg <- c(ll_agg, 'BA.5')
          }
          if(exists('ll_agg')){
            # if ll_agg contains subvariants of ll, then add a new row to agg_var_mat
            # if ll_agg only contains ll, don't add a new row
            if(any(ll_agg != ll)){
              # create an extra row for the agg_var_mat
              extra_row <- ifelse(colnames(agg_var_mat) %in% ll_agg,1,0)

              # add the new row onto the aggregation matrix
              agg_var_mat <- rbind(
                agg_var_mat,
                extra_row
              )
              row.names(agg_var_mat)[nrow(agg_var_mat)] <- paste(ll, 'Aggregated')

              # remove aggregation indices from "Omicron Aggregated" if they're in the new row
              om_row <- which(row.names(agg_var_mat) == 'Omicron Aggregated')
              agg_var_mat[om_row,which(agg_var_mat[ om_row,] == 1 & extra_row == 1)] <- 0
            }
          }
        }
      }
      # print a warning if any columns have totals > 1
      if(any(colSums(agg_var_mat)>1)) warning(paste0('agg_var_mat not correctly specified. Some variants are aggregated more than once.', agg_var_mat))

      # QA: make sure that any aggregated variant includes itself!
      agg_var_names <- gsub(' Aggregated', '',
                            sub('Delta', 'B.1.617.2',
                                sub('Omicron', 'B.1.1.529',
                                    row.names(agg_var_mat))))
      if(!all(diag(agg_var_mat[, agg_var_names])==1)) {
        # which variant is not included in itself?
        stop(paste0(
          "In 'agg_var_mat', ",
          agg_var_names[diag(agg_var_mat[, agg_var_names])!=1],
          ' is not included in ',
          row.names(agg_var_mat)[diag(agg_var_mat[, agg_var_names])!=1])
        )
      }

      # save agg_var_mat
      # weighting_methods (meth) will produce identical agg_var_mat for weighted & unweighted, so no need to save both. This will overwrite one with the other if running weighting_methods = "both"
      write.csv(
        x = replace(agg_var_mat, agg_var_mat == 0, NA), # it's easier to view in Excel with only the 1's and no 0's
        file = paste0(script.basename,
                      output_folder, "/agg_var_mat_",
                      ci.type,
                      "CI_",
                      svy.type,
                      "_",
                      data_date,
                      tag,
                      ".csv"),
        row.names = T,
        na = ""
      )
    } # end agg_var_mat old method


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
        wk = as.numeric(as.Date(ftn, origin="1970-01-01") - (week0day1+3) - 6.5) / 7 # this is to get the midpoint of each fortnight
        #wk = as.numeric(as.Date(ftn, origin="1970-01-01") - week0day1) %/% 7
        # convert week to model_week
        wk = wk - model_week_min - model_week_mid
        # same as date_to_model_week(ftn - 3)

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

        # empty dataframe to store growth rates (and doubling times) for aggregated variants
        gr_agg <- data.frame( variant = rownames(agg_var_mat),
                              gr    = NA,
                              se.gr = NA,
                              gr_lo = NA,
                              gr_hi = NA,
                              dt    = NA,
                              dt_lo = NA,
                              dt_hi = NA)
        # extract the growth rates for the aggregated variants
        for(r in 1:nrow(agg_var_mat)){
          # if nothing is actually being aggregated, just get the growth rate of the individual component
          if(unname(rowSums(agg_var_mat)[r]) == 1){
            col_ind <- which(agg_var_mat[r,]>0)
            gr_agg[r,'gr']    <- gr[col_ind]
            gr_agg[r,'se.gr'] <- se.gr[col_ind]
            gr_agg[r,'gr_lo'] <- gr_lo[col_ind]
            gr_agg[r,'gr_hi'] <- gr_hi[col_ind]
            gr_agg[r,'dt']    <- doubling_time[col_ind]
            gr_agg[r,'dt_lo'] <- doubling_time_lo[col_ind]
            gr_agg[r,'dt_hi'] <- doubling_time_hi[col_ind]
          } else {
            # if there are component variants, then take the weighted mean to get the aggregated growth rate
            col_ind <- unname(which(agg_var_mat[r,]>0))
            gr_agg[r,'gr']    <- sum(   gr[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
            gr_agg[r,'gr_lo'] <- sum(gr_lo[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
            gr_agg[r,'gr_hi'] <- sum(gr_hi[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
            gr_agg[r,'dt']    <- sum(   doubling_time[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
            gr_agg[r,'dt_lo'] <- sum(doubling_time_lo[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
            gr_agg[r,'dt_hi'] <- sum(doubling_time_hi[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
          }
        }

        # format estimates into dataframe with relevant info
        ests = data.table::data.table(
          USA_or_HHSRegion = rgn,
          Fortnight_ending = as.Date(ftn, origin="1970-01-01"),
          Variant = c(model_vars,
                      "Other",
                      row.names(ests$composite_variant$matrix)),
          Share = c(ests$p_i,
                    ests$composite_variant$p_i),
          se.Share = c(ests$se.p_i,
                       ests$composite_variant$se.p_i),
          growth_rate    = c(gr,    gr_agg$gr),
          growth_rate_lo = c(gr_lo, gr_agg$gr_lo),
          growth_rate_hi = c(gr_hi, gr_agg$gr_hi),
          doubling_time    = c(doubling_time,    gr_agg$dt),
          doubling_time_lo = c(doubling_time_lo, gr_agg$dt_lo),
          doubling_time_hi = c(doubling_time_hi, gr_agg$dt_hi)
        )

        # Get binomial CI from p_i and se.p_i
        binom.ci = apply(X = ests,
                         MARGIN = 1,
                         FUN = function(rr) svyCI(p = as.numeric(rr[4]),
                                                  s = as.numeric(rr[5])))

        # add the CI into the estimates dataframe
        ests$Share_lo = binom.ci[1,]
        ests$Share_hi = binom.ci[2,]

        # optionally add 99% prediction interval
        if(calc_99_CI_nowcast){
          binom.ci.99 <- array(data = numeric(), dim = c(2, nrow(ests)))
          for (rr in seq(nrow(ests))){
            binom.ci.99[,rr] <- svyCI(p = as.numeric(ests[rr,4]),
                                      s = as.numeric(ests[rr,5]),
                                      conf.level = 0.99)
          }

          # add the CI into the estimates dataframe
          ests$Share_lo_99 = binom.ci.99[1,]
          ests$Share_hi_99 = binom.ci.99[2,]
        }

        # add the estimates for this specific place & time to the results
        proj.res = rbind(proj.res,
                         ests)
      } # end loop over fortnights
    } # end loop over regions

    # select just the columns we want
    proj.res_column_names <- c('USA_or_HHSRegion',
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
                               'doubling_time_hi')
    if(calc_99_CI_nowcast) proj.res_column_names <- c(proj.res_column_names, 'Share_lo_99', 'Share_hi_99')
    proj.res = proj.res[,..proj.res_column_names]

    # optionally calculate the number of infections attributable to each variant
    if (calc_confirmed_infections){
      test_filepath <- paste0(script.basename,
                              "/data/backup_",
                              data_date, custom_tag, "/",
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
                    tests_fn_hhs)[,.(fortnight_end = as.Date(fortnight_end),
                                     total_test_positives = total_test_positives,
                                     HHS = HHS)],
          by.x = c("USA_or_HHSRegion",
                   "Fortnight_ending"),
          by.y = c('HHS',
                   'fortnight_end'),
          all.x = TRUE)

        # calculate case totals for each variant
        proj.res[, cases    := total_test_positives * Share]
        proj.res[, cases_lo := total_test_positives * Share_lo]
        proj.res[, cases_hi := total_test_positives * Share_hi]
      } else {
        print(paste0('File ',
                     test_filepath,
                     ' not found. Not calculating number of infections attributable to each variant for fortnights.'))
      }
    }

    if(pre_aggregation==FALSE){
      # Format output for the run 1 lineage list
      # exclude variants that have been aggregated into other groups (i.e. have value > 0 in the matrix)
      # run_1 = proj.res[proj.res$Variant %notin% colnames(agg_var_mat)[colSums(agg_var_mat)>0],]
      agg_lineages <- colnames(agg_var_mat)[colSums(agg_var_mat)>0] # need to make sure that run1 keeps either "Other" or Other aggregated! (not both or neither)
      if("Other" %notin% agg_lineages) agg_lineages <- c(agg_lineages, "Other Aggregated")
      run_1 = proj.res[Variant %notin% agg_lineages]

      # change the name of "Other Aggregated" to "Other" to match other output files
      run_1[run_1$Variant == "Other Aggregated","Variant"] <- "Other"

      # output data file depends on NO change in run_1 column order!!!
      # If new columns are to be added, they need to go to the end.

      # QA: make sure that the shares add up to 1 each time period to make sure that the aggregation is doing what I want it to:
      if (!all(run_1[, .(total_share = sum(Share)), by = c('USA_or_HHSRegion', 'Fortnight_ending')][,unique(round(total_share, 5))] == 1)){
        warning(paste(
          paste0(script.basename,
                 output_folder, "/updated_nowcast_fortnightly_",
                 data_date,
                 sub(pattern = '2', replacement = '1', x = tag),
                 ".csv"),
          'results are invalid! The total proportion does not add up to 100% in each time period!')
        )
      } else {
        # only save the results to file if the proportions add up to 100% each week
        # save the results to file
        write.csv(x = run_1,
                  file = paste0(script.basename,
                                output_folder, "/updated_nowcast_fortnightly_",
                                meth,
                                '_',
                                data_date,
                                sub(pattern = '2', replacement = '1', x = tag),
                                "_",
                                results_tag,
                                ".csv"),
                  row.names = FALSE)
        if(meth == 'weighted'){
          # process dataframe and save to a format for direct hadoop upload
          run_1_hadoop = data.frame(run_1[,1:6])
          run_1_hadoop[is.na(run_1_hadoop)] = '\\N'
          run_1_hadoop[,7:15] = '\\N'
          run_1_hadoop[,16] = 'smoothed'
          run_1_hadoop[,17] = 'biweekly'
          run_1_hadoop[,18] = data_date
          run_1_hadoop[,19] = paste0(results_tag, '_Run1')
          run_1_hadoop[,20] = 1
          run_1_hadoop[,21:24] = run_1[,14:17]
          run_1_hadoop[,3] = gsub('Delta Aggregated', 'B.1.617.2', run_1_hadoop[,3])
          run_1_hadoop[,3] = gsub('Omicron Aggregated', 'B.1.1.529', run_1_hadoop[,3])
          run_1_hadoop[,3] = gsub(' Aggregated', '', run_1_hadoop[,3])
          write.table(x = run_1_hadoop,
                      file = paste0(script.basename,
                                    output_folder, "/updated_nowcast_fortnightly_",
                                    data_date,
                                    sub(pattern = '2', replacement = '1', x = tag),
                                    "_",
                                    results_tag,
                                    "_hadoop.csv"),
                      quote = FALSE,
                      row.names = FALSE,
                      col.names = FALSE,
                      sep = ",")
        } # end save hadoop data
      } # end save data
    } # end save data when pre-aggregation == FALSE

    # Format output for the run2 lineage list
    # exclude the lineages that were aggregated (other than "Other")
    drop_lin <- row.names(agg_var_mat)[row.names(agg_var_mat) %notin% "Other Aggregated"]
    # alternatively, just include only the variants that are in c(voc, "Other Aggregated")

    # Only include variants that are NOT in the list provided
    # (also use "Other Aggregated" instead of "Other"; this requires dropping the variants that were aggregated into "Other Aggregated")
    # output data file depends on NO change in run_2 column order!!!
    # If new columns are to be added, they need to go to the end.
    run_2 = proj.res[Variant %notin% c(drop_lin, "Other", colnames(agg_var_mat['Other Aggregated',,drop=F])[agg_var_mat['Other Aggregated',,drop=F] > 0])]

    # change the name of "Other Aggregated" to "Other" to match other output files
    run_2[run_2$Variant=="Other Aggregated","Variant"] <- "Other"

    # QA: make sure that the shares add up to 1 each week to make sure that the aggregation is doing what I want it to:
    if (!all(run_2[, .(total_share = sum(Share)), by = c('USA_or_HHSRegion', 'Fortnight_ending')][,unique(round(total_share, 5))] == 1)){
      warning(paste(
        paste0(script.basename,
               output_folder, "/updated_nowcast_fortnightly_",
               data_date,
               tag,
               ".csv"),
        'results are invalid! The total proportion does not add up to 100% in all weeks!')
      )
    } else {
      # only save the results to file if the proportions add up to 100% each week
      write.csv(x = run_2,
                file = paste0(script.basename,
                              output_folder, "/updated_nowcast_fortnightly_",
                              meth,
                              "_",
                              data_date,
                              tag,
                              "_",
                              results_tag,
                              ".csv"),
                row.names = FALSE)
      if(meth == "weighted"){
        # process dataframe and save to a format for direct hadoop upload
        run_2_hadoop = data.frame(run_2[,1:6])
        run_2_hadoop[is.na(run_2_hadoop)] = '\\N'
        run_2_hadoop[,7:15] = '\\N'
        run_2_hadoop[,16] = 'smoothed'
        run_2_hadoop[,17] = 'biweekly'
        run_2_hadoop[,18] = data_date
        run_2_hadoop[,19] = paste0(results_tag, '_Run2')
        run_2_hadoop[,20] = 1
        run_2_hadoop[,21:24] = run_2[,14:17]
        run_2_hadoop[,3] = gsub('Delta Aggregated', 'B.1.617.2', run_2_hadoop[,3])
        run_2_hadoop[,3] = gsub('Omicron Aggregated', 'B.1.1.529', run_2_hadoop[,3])
        run_2_hadoop[,3] = gsub(' Aggregated', '', run_2_hadoop[,3])
        write.table(x = run_2_hadoop,
                    file = paste0(script.basename,
                                  output_folder, "/updated_nowcast_fortnightly_",
                                  data_date,
                                  tag,
                                  "_",
                                  results_tag,
                                  "_hadoop.csv"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE,
                    sep = ",")
      } # end save hadoop data
    } # end save run2 fortnightly data

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

        # empty dataframe to hold growth rates of aggregated variants
        gr_agg <- data.frame( variant = rownames(agg_var_mat),
                              gr    = NA,
                              se.gr = NA,
                              gr_lo = NA,
                              gr_hi = NA,
                              dt    = NA,
                              dt_lo = NA,
                              dt_hi = NA)
        # extract the growth rates for the aggregated variants
        for(r in 1:nrow(agg_var_mat)){
          # if nothing is actually being aggregated, just get the growth rate of the individual component
          if(unname(rowSums(agg_var_mat)[r]) == 1){
            col_ind <- which(agg_var_mat[r,]>0)
            gr_agg[r,'gr']    <- gr[col_ind]
            gr_agg[r,'se.gr'] <- se.gr[col_ind]
            gr_agg[r,'gr_lo'] <- gr_lo[col_ind]
            gr_agg[r,'gr_hi'] <- gr_hi[col_ind]
            gr_agg[r,'dt']    <- doubling_time[col_ind]
            gr_agg[r,'dt_lo'] <- doubling_time_lo[col_ind]
            gr_agg[r,'dt_hi'] <- doubling_time_hi[col_ind]
          } else {
            # if there are component variants, then take the weighted mean to get the aggregated growth rate
            col_ind <- unname(which(agg_var_mat[r,]>0))
            gr_agg[r,'gr']    <- sum(   gr[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
            gr_agg[r,'gr_lo'] <- sum(gr_lo[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
            gr_agg[r,'gr_hi'] <- sum(gr_hi[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
            gr_agg[r,'dt']    <- sum(   doubling_time[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
            gr_agg[r,'dt_lo'] <- sum(doubling_time_lo[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
            gr_agg[r,'dt_hi'] <- sum(doubling_time_hi[col_ind] * ests$p_i[col_ind]) / sum(ests$p_i[col_ind])
          }
        }

        # format estimates into dataframe with relevant info
        ests = data.table::data.table(
          USA_or_HHSRegion = rgn,
          Week_ending = week_ending, # this no longer identifies a single estimate.
          Variant = c(model_vars,
                      "Other",
                      row.names(ests$composite_variant$matrix)),
          Share = c(ests$p_i,
                    ests$composite_variant$p_i),
          se.Share = c(ests$se.p_i,
                       ests$composite_variant$se.p_i),
          growth_rate    = c(gr,    gr_agg$gr),
          growth_rate_lo = c(gr_lo, gr_agg$gr_lo),
          growth_rate_hi = c(gr_hi, gr_agg$gr_hi),
          doubling_time    = c(doubling_time,    gr_agg$dt),
          doubling_time_lo = c(doubling_time_lo, gr_agg$dt_lo),
          doubling_time_hi = c(doubling_time_hi, gr_agg$dt_hi),
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

        if(calc_99_CI_nowcast){
          binom.ci.99 <- array(data = numeric(), dim = c(2, nrow(ests)))
          for (rr in seq(nrow(ests))){
            binom.ci.99[,rr] <- svyCI(p = as.numeric(ests[rr,4]),
                                      s = as.numeric(ests[rr,5]),
                                      conf.level = 0.99)
          }

          # add the CI into the estimates dataframe
          ests$Share_lo_99 = binom.ci.99[1,]
          ests$Share_hi_99 = binom.ci.99[2,]
        }


        # add the estimates for this specific place & time to the results
        proj.res = rbind(proj.res,
                         ests)
      } # end loop over weeks
    } # end loop over regions

    # select column names to save
    proj.res_column_names <- c(
      "USA_or_HHSRegion",
      "Week_ending",
      "Variant",
      "Share",
      "Share_lo",
      "Share_hi",
      "se.Share",
      "growth_rate",
      "growth_rate_lo",
      "growth_rate_hi",
      "doubling_time",
      "doubling_time_lo",
      "doubling_time_hi",
      "date",
      "week_start",
      "model_week"
    )
    if(calc_99_CI_nowcast) proj.res_column_names <- c(proj.res_column_names, 'Share_lo_99', 'Share_hi_99')
    # subset to the columns we want
    proj.res = proj.res[,..proj.res_column_names]

    # optionally calculate the number of infections attributable to each variant
    if (calc_confirmed_infections){
      test_filepath <- paste0(script.basename,
                              "/data/backup_",
                              data_date, custom_tag, "/",
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
        proj.res[, cases_daily    := total_test_positives_daily * Share]
        proj.res[, cases_lo_daily := total_test_positives_daily * Share_lo]
        proj.res[, cases_hi_daily := total_test_positives_daily * Share_hi]

        proj.res[, cases_weekly    := total_test_positives_weekly * Share]
        proj.res[, cases_lo_weekly := total_test_positives_weekly * Share_lo]
        proj.res[, cases_hi_weekly := total_test_positives_weekly * Share_hi]
      } else {
        print(paste0('File ',
                     test_filepath,
                     ' not found. Not calculating number of infections attributable to each variant for weeks.'))
      }

    }

    if(pre_aggregation==FALSE){
      # Format output for the run 1 lineage list
      # exclude variants that have been aggregated into other groups (i.e. have value > 0 in the matrix)
      # run_1 = proj.res[proj.res$Variant %notin% colnames(agg_var_mat)[colSums(agg_var_mat)>0],]
      agg_lineages <- colnames(agg_var_mat)[colSums(agg_var_mat)>0] # need to make sure that run1 keeps either "Other" or Other aggregated! (not both or neither)
      if("Other" %notin% agg_lineages) agg_lineages <- c(agg_lineages, "Other Aggregated")
      run_1 = proj.res[Variant %notin% agg_lineages]

      # change the name of "Other Aggregated" to "Other" to match other output files
      run_1[run_1$Variant == "Other Aggregated","Variant"] <- "Other"

      # weekly results = remove daily results & daily columns
      # output data file depends on NO change in run_1_weekly column order!!!
      # If new columns are to be added, they need to go to the end.
      run_1_weekly <- data.table:::subset.data.table(x = run_1,
                                                     subset = model_week %% 1 == 0,
                                                     select = !names(run_1) %in% c('total_test_positives_daily', 'cases_daily', 'cases_lo_daily', 'cases_hi_daily'))
      # daily results = remove weekly columns
      run_1_daily <- data.table:::subset.data.table(x = run_1,
                                                    select = !names(run_1) %in% c('total_test_positives_weekly', 'cases_weekly', 'cases_lo_weekly', 'cases_hi_weekly'))

      # QA: make sure that the shares add up to 1 each week to make sure that the aggregation is doing what I want it to:
      if (!all(run_1_weekly[, .(total_share = sum(Share)), by = c('USA_or_HHSRegion', 'Week_ending')][,unique(round(total_share, 5))] == 1)){
        warning(paste(
          paste0(script.basename,
                 output_folder, "/updated_nowcast_weekly_",
                 data_date,
                 sub(pattern = '2', replacement = '1', x = tag),
                 ".csv"),
          'results are invalid! The total proportion does not add up to 100% in all weeks!')
        )
      } else {
        # only save the results to file if the proportions add up to 100% each week
        # save the results to file
        write.csv(x = run_1_weekly,
                  file = paste0(script.basename,
                                output_folder, "/updated_nowcast_weekly_",
                                meth,
                                "_",
                                data_date,
                                sub(pattern = '2', replacement = '1', x = tag),
                                "_",
                                results_tag,
                                ".csv"),
                  row.names = FALSE)
        # process dataframe and save to a format for direct hadoop upload
        if(meth == 'weighted'){
          run_1_weekly_hadoop = data.frame(run_1_weekly[,c(1:2,4:7)])
          run_1_weekly_hadoop[is.na(run_1_weekly_hadoop)] = '\\N'
          run_1_weekly_hadoop[,7:15] = '\\N'
          run_1_weekly_hadoop[,16] = 'smoothed'
          run_1_weekly_hadoop[,17] = 'weekly'
          run_1_weekly_hadoop[,18] = data_date
          run_1_weekly_hadoop[,19] = paste0(results_tag, '_Run1')
          run_1_weekly_hadoop[,20] = 1
          run_1_weekly_hadoop[,21:24] = run_1_weekly[,17:20]
          run_1_weekly_hadoop[,3] = gsub('Delta Aggregated', 'B.1.617.2', run_1_weekly_hadoop[,3])
          run_1_weekly_hadoop[,3] = gsub('Omicron Aggregated', 'B.1.1.529', run_1_weekly_hadoop[,3])
          run_1_weekly_hadoop[,3] = gsub(' Aggregated', '', run_1_weekly_hadoop[,3])
          write.table(x = run_1_weekly_hadoop,
                      file = paste0(script.basename,
                                    output_folder, "/updated_nowcast_weekly_",
                                    data_date,
                                    sub(pattern = '2', replacement = '1', x = tag),
                                    "_",
                                    results_tag,
                                    "_hadoop.csv"),
                      quote = FALSE,
                      row.names = FALSE,
                      col.names = FALSE,
                      sep = ",")
        } # end save hadoop data
      } # end save run1 data
      # daily results
      # QA: make sure that the shares add up to 1 each week to make sure that the aggregation is doing what I want it to:
      if (!all(run_1_daily[, .(total_share = sum(Share)), by = c('USA_or_HHSRegion', 'date')][,unique(round(total_share, 5))] == 1)){
        warning(paste(
          paste0(script.basename,
                 output_folder, "/updated_nowcast_weekly_",
                 data_date,
                 sub(pattern = '2', replacement = '1', x = tag),
                 "_daily.csv"),
          'results are invalid! The total proportion does not add up to 100% in all weeks!')
        )
      } else {
        # only save the results to file if the proportions add up to 100% each week
        # save the results to file
        write.csv(x = run_1_daily,
                  file = paste0(script.basename,
                                output_folder, "/updated_nowcast_weekly_",
                                meth,
                                "_",
                                data_date,
                                sub(pattern = '2', replacement = '1', x = tag),
                                "_",
                                results_tag,
                                "_daily.csv"),
                  row.names = FALSE)
      }
    }
    # Format output for the run2 lineage list
    # exclude the lineages that were aggregated (other than "Other")
    drop_lin <- row.names(agg_var_mat)[row.names(agg_var_mat) %notin% "Other Aggregated"]
    # alternatively, just include only the variants that are in c(voc, "Other Aggregated")

    # Only include variants that are NOT in the list provided
    # (also use "Other Aggregated" instead of "Other"; this requires dropping the variants that were aggregated into "Other Aggregated")
    run_2 = proj.res[Variant %notin% c(drop_lin, "Other", colnames(agg_var_mat['Other Aggregated',,drop=F])[agg_var_mat['Other Aggregated',,drop=F] > 0])]

    # change the name of "Other Aggregated" to "Other" to match other output files
    run_2[run_2$Variant=="Other Aggregated","Variant"] <- "Other"

    # weekly results = remove daily results & daily columns
    # output data file depends on NO change in run_2_weekly column order!!!
    # If new columns are to be added, they need to go to the end.
    run_2_weekly <- data.table:::subset.data.table(x = run_2,
                                                   subset = model_week %% 1 == 0,
                                                   select = !names(run_2) %in% c('total_test_positives_daily', 'cases_daily', 'cases_lo_daily', 'cases_hi_daily'))
    # daily results = remove weekly columns
    run_2_daily <- data.table:::subset.data.table(x = run_2,
                                                  select = !names(run_2) %in% c('total_test_positives_weekly', 'cases_weekly', 'cases_lo_weekly', 'cases_hi_weekly'))

    # QA: make sure that the shares add up to 1 each week to make sure that the aggregation is doing what I want it to:
    if (!all(run_2_weekly[, .(total_share = sum(Share)), by = c('USA_or_HHSRegion', 'Week_ending')][,unique(round(total_share, 5))] == 1)){
      warning(paste(
        paste0(script.basename,
               output_folder, "/updated_nowcast_weekly_",
               data_date,
               tag,
               ".csv"),
        'results are invalid! The total proportion does not add up to 100% in all weeks!')
      )
    } else {
      # only save the results to file if the proportions add up to 100% each week
      # save the results to file
      write.csv(x = run_2_weekly,
                file = paste0(script.basename,
                              output_folder, "/updated_nowcast_weekly_",
                              meth,
                              "_",
                              data_date,
                              tag,
                              "_",
                              results_tag,
                              ".csv"),
                row.names = FALSE)
      # process dataframe and save to a format for direct hadoop upload
      if(meth == 'weighted'){
        run_2_weekly_hadoop = data.frame(run_2_weekly[,c(1:2,4:7)])
        run_2_weekly_hadoop[is.na(run_2_weekly_hadoop)] = '\\N'
        run_2_weekly_hadoop[,7:15] = '\\N'
        run_2_weekly_hadoop[,16] = 'smoothed'
        run_2_weekly_hadoop[,17] = 'weekly'
        run_2_weekly_hadoop[,18] = data_date
        run_2_weekly_hadoop[,19] = paste0(results_tag, '_Run2')
        run_2_weekly_hadoop[,20] = 1
        run_2_weekly_hadoop[,21:24] = run_2_weekly[,17:20]
        run_2_weekly_hadoop[,3] = gsub('Delta Aggregated', 'B.1.617.2', run_2_weekly_hadoop[,3])
        run_2_weekly_hadoop[,3] = gsub('Omicron Aggregated', 'B.1.1.529', run_2_weekly_hadoop[,3])
        run_2_weekly_hadoop[,3] = gsub(' Aggregated', '', run_2_weekly_hadoop[,3])
        write.table(x = run_2_weekly_hadoop,
                    file = paste0(script.basename,
                                  output_folder, "/updated_nowcast_weekly_",
                                  data_date,
                                  tag,
                                  "_",
                                  results_tag,
                                  "_hadoop.csv"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE,
                    sep = ",")
      } # end save hadoop data
    } # end save run2 weekly data
    # daily results
    # QA: make sure that the shares add up to 1 each week to make sure that the aggregation is doing what I want it to:
    if (!all(run_2_daily[, .(total_share = sum(Share)), by = c('USA_or_HHSRegion', 'date')][,unique(round(total_share, 5))] == 1)){
      warning(paste(
        paste0(script.basename,
               output_folder, "/updated_nowcast_weekly_",
               data_date,
               tag,
               "_daily.csv"),
        'results are invalid! The total proportion does not add up to 100% in all weeks!')
      )
    } else {
      # only save the results to file if the proportions add up to 100% each week
      # save the results to file
      write.csv(x = run_2_daily,
                file = paste0(script.basename,
                              output_folder, "/updated_nowcast_weekly_",
                              meth,
                              '_',
                              data_date,
                              tag,
                              "_",
                              results_tag,
                              "_daily.csv"),
                row.names = FALSE)
    } # end QA check for daily results
  } # end loop over weighted_methods (unweighted, weighted)
} # end run 2 Nowcast modeling


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
                          State = c("AK", "AL", "AR", "AZ", "CA", "CO", "CT",
                                    "DC", "DE", "FL", "GA", "GU", "HI", "IA",
                                    "ID", "IL", "IN", "KS", "KY", "LA", "MA",
                                    "MD", "ME", "MI", "MN", "MO", "MP", "MS",
                                    "MT", "NC", "ND", "NE", "NH", "NJ", "NM",
                                    "NV", "NY", "OH", "OK", "OR", "PA", "PR",
                                    "RI", "SC", "SD", "TN", "TX", "UT", "VA",
                                    "VT", "WA", "WI", "WV", "WY", "AS", "MH",
                                    "VI", "PW", "FM"))[, 3:1]




  # get the proportion estimates & CI
  all_state_ests_function <- function(all.state, svyDES){
    ests <- apply(X = all.state,
                  MARGIN = 1,
                  FUN = function(rr) myciprop(voc = rr[3],
                                              geoid = rr[1],
                                              svy = subset(svyDES,
                                                           week >= (as.numeric(rr[2])- 3) &
                                                             week < (as.numeric(rr[2])+1)),
                                              str = FALSE))
    return(ests)
  }

  # do the estimates in parallel
  if(use_parallel){
    # use a tryCatch just in case the parallel operation fails
    ests <- tryCatch(expr = { # "expr" is what we want to run (not in a function form, unlike "error" and "warning")
      # choose the number of cores
      # ncores <- max(1, parallel::detectCores())

      # split the rows of data into "ncores" sets
      cut_list <- split(x = all.state,
                        f = rep(1:ncores,
                                each = ceiling(nrow(all.state)/ncores),
                                length.out = nrow(all.state)))

      # make a cluster
      cl <- parallel::makeCluster(ncores)

      # pass everything to each of the cores
      parallel::clusterEvalQ(cl = cl, {
        library(survey)
        library(data.table)
      })
      # export all the other R objects that are used within all_ftnt_ests_function to each cluster node
      parallel::clusterExport(cl = cl,
                              varlist = c('svyDES', 'script.basename', 'ci.type', 'voc'))
      parallel::clusterEvalQ(cl = cl, {
        source(paste0(script.basename, "/weekly_variant_report_functions.R"))
        options(survey.adjust.domain.lonely = T,
                survey.lonely.psu = "average",
                stringsAsFactors = FALSE)
      })

      # perform the calculations on the cluster
      ests_list <- parallel::parLapply(cl = cl,
                                       X = cut_list,
                                       fun = all_state_ests_function,
                                       svyDES = svyDES)

      # combine the results into a single dataframe
      do.call(what = 'cbind', args = ests_list)

      # close the cluster
      # parallel::stopCluster(cl)
    },
    error = function(cond) { # if "expr" throws an error, do this instead of actually stopping
      message('Parallel execution of monthly state weighted proportions failed. Running in series instead.')
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(all_state_ests_function(all.state = all.state, svyDES = svyDES))
    })
  } else ests <- all_state_ests_function(all.state = all.state, svyDES = svyDES)

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

  # get the proportion estimates & CI
  other_state_ests_function <- function(others, svyDES){
    ests.others <- apply(X = others,
                         MARGIN = 1,
                         FUN = function(rr) myciprop(voc = voc,
                                                     geoid = rr[1],
                                                     svy = subset(svyDES,
                                                                  week >= (as.numeric(rr[2])- 3) &
                                                                    week < (as.numeric(rr[2])+1)),
                                                     str = FALSE))
    return(ests.others)
  }

  # do the estimates in parallel
  if(use_parallel){
    # use a tryCatch just in case the parallel operation fails
    ests.others <- tryCatch(expr = { # "expr" is what we want to run (not in a function form, unlike "error" and "warning")
      # split the rows of data into "ncores" sets
      cut_list <- split(x = others,
                        f = rep(1:ncores,
                                each = ceiling(nrow(others)/ncores),
                                length.out = nrow(others)))

      # perform the calculations on the cluster
      ests_list <- parallel::parLapply(cl = cl,
                                       X = cut_list,
                                       fun = other_state_ests_function,
                                       svyDES = svyDES)

      # close the cluster
      parallel::stopCluster(cl)

      # combine the results into a single dataframe
      do.call(what = 'cbind', args = ests_list)
    },
    error = function(cond) { # if "expr" throws an error, do this instead of actually stopping
      message('Parallel execution of "other" monthly state weighted proportions failed. Running in series instead.')
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(other_state_ests_function(others = others, svyDES = svyDES))
    })
  } else ests.others <- other_state_ests_function(others = others, svyDES = svyDES)

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
                            data_date, custom_tag, "/",
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
                          output_folder, "/state_weighted_roll4wk_",
                          ci.type,
                          "CI_svyNEW_",
                          data_date,
                          tag,
                          ".csv"),
            row.names = FALSE)

  # Generate file ready for hadoop upload
  all.state.out_hadoop = data.frame(all.state.out[,c("State",
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
                                                     "nchs_flag_wodf")])
  all.state.out_hadoop[,14] = data_date
  all.state.out_hadoop[,15] = ''
  all.state.out_hadoop[,16:19] = all.state.out[,c("total_test_positives",
                                                  "cases",
                                                  "cases_lo",
                                                  "cases_hi")]
  all.state.out_hadoop[is.na(all.state.out_hadoop)] = '\\N'
  write.table(x = all.state.out_hadoop,
              file = paste0(script.basename,
                            output_folder, "/state_weighted_roll4wk_",
                            ci.type,
                            "CI_svyNEW_",
                            data_date,
                            tag,
                            "_hadoop.csv"),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              sep =  ",")

} # end Run3


# # save src.dat that's been prepped for analysis
if(save_datasets_to_file){
  saveRDS(object = src.dat,
          file = paste0(script.basename,
                        output_folder, '/src.dat_', # save to results instead of 'data' folder
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
