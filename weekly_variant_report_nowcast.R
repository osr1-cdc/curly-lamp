#-------------------------------------------------------------------------------
#
#   Variant share estimates and nowcasts
#   created by Prabasaj Paul
#
# Updates:
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
#

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


## optparse option list
# (get the run number from the command line)
option_list <- list(
  # Run number
  optparse::make_option(
    opt_str = c("-r", "--run_number"),
    type    = "character",
    default = "1",
    help    = "Run number",
    metavar = "character"),
  # whether or not to use Custom Lineages
  optparse::make_option(
    opt_str = c("-c", "--custom_lineages"),
    type    = "character",
    default = "F",
    help    = "Whether or not to use custom lineages (character value of T or F)",
    metavar = "character"
  )
)

# parseing options list
opts <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# use the specified flags to set several variables
# convert custom_lineages flag to logical value
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
dir.create(paste0(script.basename,"/results"))


# Data prep --------------------------------------------------------------------
#capture system time
tstart = proc.time()

# source options from config.R file
source(paste0(script.basename, "/config/config.R"))
#Source the code with svycipropkg function created by Crescent Martin (odb4).
#NOTE: the svyciprop function from the survey package uses the calculated
#     effective sample size for the KG CI (which can be larger than the actual sample size)
#     The svycipropkg function caps the effective sample size at the actual sample size
#CAVEATE: the svycipropkg code hasn't gone through comprehensive testing.
source(paste0(script.basename, "/svycipropkg.R"))

# Load output from variant_surveillance_system.r
# (genomic surveillance data plus survey weights)
load(paste0(script.basename, "/data/svydat_", data_date, custom_tag, ".RData"))

# create a tag for the filenames to differentiate results from different runs
tag <- paste0("_",state_source,"_Run", opts$run_number, custom_tag)
# options:
# Run1: calc proportions USING SURVEY DESIGN for reduced set of VOCs
# Run2: calc proportions for extended set of VOCs (many Delta subclades)
#       this is the only run that includes the multinomial "nowcast" model
#     Note: Nowcast runs best with a large set of variants (if you group everything
#           into delta, then the model breaks)
# Run3: calc proportions for another reduced set of VOCs specific to state-level
#       runs & VOC's generally won't change

# Choose which list of vocs to use based on the run number
if( grepl("Run1",tag) ){
  if(custom_lineages == FALSE) {
    voc = voc1
  } else {
    voc = voc1_custom
  }
}
if( grepl("Run2",tag) ) {
  # get the list of vocs from either the manually specified values in "config/config.R"
  # or from the automatically calculated values in "variant_surveillance_system.R"
  if( is.na(voc2_manual) ) {
    # read in the list of voc created in "variant_surveillance_system.R"
    # and add on any variants specified in voc2_additional
    voc2 = unique( # make sure there are no duplicates
      c(
        readRDS(file = paste0(script.basename,
                              "/data/voc2_auto", data_date, custom_tag, ".RDS")),
        # add in the "additional" vocs that need to be included
        voc2_additional
        ))
  } else {
    voc2 = voc2_manual
  }

  # possibly add on custom lineages
  if(custom_lineages == FALSE) {
    voc = voc2
  } else {
    voc = c(voc2, custom_lineage_names) # "custom_lineage_names" defined in "config/config.R"
  }
}
if( grepl("Run3",tag) ) {
  if(custom_lineages == FALSE) {
    voc = voc3
  } else {
    voc = voc3_custom
  }
}

# Convert any factor to string
fac2str = sapply(svy.dat, class) # get the class of each column
fac2str = names(fac2str[fac2str=="factor"]) # names of the columns that are factors
# convert columns that are factors to strings
for (vv in fac2str) svy.dat[, vv] = as.character(svy.dat[, vv])

# Define "source" as "LAB2"; used for survey design weights
svy.dat$SOURCE = svy.dat$LAB2

# create a variable for aggregating counts of different groupings of samples
svy.dat$count=1

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
                                 "FULGENT GENETICS",
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
}

# set s-gene upsampling weights to 1 for samples from weeks where no labs did any
# upsampling or for samples with weight of NA
src.dat$sgtf_weights[is.na(src.dat$sgtf_weights) |
                       !src.dat$SGTF_UPSAMPLING] = 1

# (Weighted) count of sequences in each state & week
seq.tbl = with(src.dat, xtabs((1/sgtf_weights) ~ STUSAB + yr_wk))


# (re)calculate the survey weights
# Note: these weights are already calculated in "variant_surveillance_system.R",
#       but they're recalculated here after subsampling data based on lab.
# For info on how these weights are calculated, see notes in "variant_surveillance_system.R"
# and: https://www.cdc.gov/mmwr/volumes/70/wr/mm7023a3.htm
# (using data.table saves about 99% of the time compared to for-loop)

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


# Make sure all survey weights are valid numbers
# (this shouldn't be necessary, just precautionary)
src.dat = subset(x = src.dat,
                 !is.na(SIMPLE_ADJ_WT) &
                   SIMPLE_ADJ_WT < Inf)

# Normalize the weights for use in the multinomial model
src.dat$NORM_WTS = with(src.dat,
                        SIMPLE_ADJ_WT/sum(SIMPLE_ADJ_WT,
                                          na.rm = TRUE))
# normalize so the total weights = # of sequences (instead of 1)
src.dat$NORM_WTS = src.dat$NORM_WTS * sum(!is.na(src.dat$NORM_WTS))

# calculate final date of 2-week bins for each observation
src.dat$FORTNIGHT_END = as.character(week0day1 + src.dat$week%/%2 * 14 + 13)
# Create variable that indicates date ends for 4-week time bins
src.dat$FOURWEEK_END = as.character(week0day1 + src.dat$week%/%4 * 28 + 27)

# make sure "VARIANT" is a character (rather than factor)
src.dat$VARIANT = as.character(src.dat$VARIANT)

#Identify all the clades/lineages to aggregate in the surveillance dataset
# all the AY variants
AY = sort(unique(src.dat$VARIANT)[grep("AY",unique(src.dat$VARIANT))])
# just the AY variants that are to be aggregated (i.e. not listed in "voc")
AY = AY[which(AY %notin% voc)]
P1=sort(unique(src.dat$VARIANT)[grep("P.1.",unique(src.dat$VARIANT))])
P1=P1[which(P1 %notin% voc)] #vector of the P1s to aggregate
Q=sort(unique(src.dat$VARIANT)[grep("Q.",unique(src.dat$VARIANT))])
Q=Q[which(Q %notin% voc)] #vector of the Qs to aggregate
B351=sort(unique(src.dat$VARIANT)[grep("B.1.351.",unique(src.dat$VARIANT))])
B351=B351[which(B351 %notin% voc)] #vector of the B351s to aggregate
B621=sort(unique(src.dat$VARIANT)[grep("B.1.621.",unique(src.dat$VARIANT))])
B621=B621[which(B621 %notin% voc)] #vector of the B621s to aggregate
B429=sort(unique(src.dat$VARIANT)[grep("B.1.429",unique(src.dat$VARIANT))])
B429=B429[which(B429 %notin% voc)] #vector of the B621s to aggregate

# Aggregate sublineages to the parent lineage
if(P.1_agg==TRUE)     {src.dat[src.dat$VARIANT %in% P1,  "VARIANT"] <- "P.1"}
if(B.1.351_agg==TRUE) {src.dat[src.dat$VARIANT %in% B351,"VARIANT"] <- "B.1.351"}
if(B.1.621_agg==TRUE) {src.dat[src.dat$VARIANT %in% B621,"VARIANT"] <- 'B.1.621' }
if(Q.1_3_agg==TRUE)   {src.dat[src.dat$VARIANT %in% Q,   "VARIANT"] <- "B.1.1.7"}
if(AY_agg==TRUE)      {src.dat[src.dat$VARIANT %in% AY,  "VARIANT"] <- "B.1.617.2"}
if(B429_7_agg==TRUE)  {src.dat[src.dat$VARIANT %in% B429,"VARIANT"] <- "B.1.427"}



# Functions --------------------------------------------------------------------
# Helper wrapper for svyciprop
# this returns confidence intervals based on survey design only
# (no multinomial model for temporal smoothing)
myciprop = function(voc,
                    geoid,
                    svy,
                    str   = TRUE,
                    range = FALSE,
                    mut   = FALSE,
                    level = 0.95,
                    ...) {
  # Arguments:
  #  ~ voc:   character string of the variants of concern to analyze
  #           (if multiple listed, estimated proportion will be for the group of
  #           variants)
  #  ~ geoid: character string/numeric of the geographic resolution to analyze
  #           (options: 2 letter state code, 1:10 for HHS region, or "USA")
  #  ~ svy:   survey design object
  #  ~ str:   logical argument indicating if output should be formatted in single
  #           character object or as vector/dataframe
  #  ~ range: logical argument indicating whether to present CI range vs CI bounds
  #  ~ mut:   logical argument indicating whether to analyze mutation profiles or
  #           lineages
  #  ~ level: numeric specifying confidence level (1-alpha)
  #  ~ ...:   when method = "asin" or "mean", these arguments are passed on to
  #           svyciprop (which passes them to "confint.svystat")

  # Output: Confidence Intervals in 1 of several formats:
  #  1) if str == FALSE: named numeric vector with 7 values:
  #                      estimate: estimated proportion
  #                      lcl: lower confidence limit
  #                      ucl: upper confidence limit
  #                      DF: degrees of freedom
  #                      n.eff.capped: number of effective samples (capped to actual number sampled)
  #                      CVmean: coefficient of variation of the mean (estimated proportion)
  #                      deffmean: design effect of the mean (estimated proportion)
  #  2) if str == TRUE:
  #         if range == TRUE:
  #               string range: "ll.l-hh.h"
  #         if range == FALSE:
  #               string: "pp.p (ll.l-hh.h) DF=__, Eff.sampsize=__, CVprop=__, DEff=__"

  # Create a binary indicator variable to ID samples to analyze
  # (based on either S_MUT or VARIANT)
  if (mut) {
    VOC = grepl(
      pattern = paste0("^(?=.*\\",
                       paste(voc, collapse="\\b)(?=.*\\"),"\\b)"),
      x = svy$variables$S_MUT,
      perl = T
    )
  } else {
    VOC = (svy$variables$VARIANT %in% voc)
  }

  # extract the relevant survey design
  srv_design = subset(update(svy, VOC = VOC),
                      (STUSAB %in% geoid) |
                        (HHS %in% geoid) |
                        (geoid == "USA"))

  # calculate the confidence interval
  # (of all specified "voc" grouped together and all non-"voc" grouped together)
  if (ci.type == "xlogit"){
    # (note: this will throw a lot of warnings b/c of strata that have single
    # PSU's, so I'll suppress warnings for just this function)
    res = suppressWarnings(
      survey::svyciprop(~VOC,
                        design = srv_design,
                        method = "xlogit",
                        ...)
    )
  } else {
    # use a modified version of svyciprop to limit effective sample size so that
    # it's never > observed sampled size.
    res = suppressWarnings(
      svycipropkg(~VOC,
                  design = srv_design,
                  ...)
    )
  }

  # svyciprop & svycipropkg return a single number (the proportion estimate)
  # with the confidence interval in the attributes.
  # Add the confidence interval directly to the results
  res = c('estimate' = unname(res[1]),
          'lcl' = survey:::confint.svyciprop(res)[1],
          'ucl' = survey:::confint.svyciprop(res)[2])

  # add on the degrees of freedom
  res = c(res,
          'DF' = survey::degf(design = srv_design))


  #add effective sample size to res output
  m <- eval(bquote( # eval and bquote are superfluous b/c there's no variable substitution going on...
    # suppress warnings caused by single PSU in stratum
    suppressWarnings(
      #extracts the mean and SE for each estimate
      survey::svymean(x = ~as.numeric(VOC),
                      design = srv_design)
    )))

  # coefficient of variation
  CVmean <- suppressWarnings(
    # suppress warnings caused by single PSU in stratum
    # then calculate CV of the mean
    survey::cv(object = survey::svymean(x = ~as.numeric(VOC),
                                        design = srv_design))
  )

  # calculate the design effect
  # "The design effect compares the variance of a mean or total to the variance
  #  from a study of the same size using simple random sampling without replacement"
  deffmean <- suppressWarnings(
    # suppress warnings caused by single PSU in stratum
    survey::deff(object = survey::svymean(x = ~as.numeric(VOC),
                                          design = srv_design,
                                          deff = T))
  )

  # n_effective (this comes from the "survey::svyciprop" function when using the
  #              beta method for calculating the CI)
  #              coef(m) = mean     vcov(m) = variance
  n.eff <- coef(m) * (1 - coef(m))/vcov(m)

  # define alpha based on CI level
  alpha <- 1-level

  # this is also from the "survey::svyciprop" function using method == 'beta'
  n.eff.capped <- min(
    n.eff * (qt(p = alpha/2,
                df = nrow(srv_design) - 1) /
               qt(p = alpha/2,
                  df = survey::degf(design = srv_design)))^2,
    nrow(srv_design)
  )

  # add more info into the results
  res = c(res,
          'n.eff.capped' = n.eff.capped,
          'CVmean'= CVmean,
          'deffmean' = unname(deffmean))

  # optionally convert results to a string
  if (str) {
    res = c(round( 100 * res[1:3], 1),
            res[4],
            res[5],
            res[6],
            res[7])

    # optionally format results as a range
    if (range) {
      res = paste0(res[2], "-", res[3])
    } else {
      res = paste0(res[1],
                   " (",res[2],"-",res[3],
                   ") DF=",res[4],
                   ", Eff.sampsize=",res[5],
                   ", CVprop=",res[6],
                   ",DEff=",res[7])
    }
  }

  # return the results
  res
}

# function to calculate & format a date based on # of weeks from data_date
week_label = function(week_from_current,
                      datadate = data_date) {
  # returns the date of the first day of a future week in "%m-%d" format

  # the current day of the week
  day_of_week = as.numeric(strftime(as.Date(datadate), format="%w"))

  # start of current week
  start_of_week = as.Date(datadate) - day_of_week

  # formatted date of a future week
  strftime(x = start_of_week + 7 * week_from_current,
           format = "%m-%d")
}


#Functions for multinomial nowcast model
# this function fits a multinomial model to the data to get temporally smoothed
# proportion estimates while also accounting survey design.
# (only run by "Run2" tag)
# Note: "svymultinom" is used in conjunction with "se.multinom" and "svyCI"
#       functions.
#       - "svymultinom" fits the model and adjusts the variance-covariance
#          matrix (of the model coefficients) for the sampling design.
#       - "se.multinom" uses the model coefficients from "svymultinom" to calculate
#          the linear predictor values (and covariances), and then converts the
#          linear predictor values into proportions and calculates the covariances
#          of those proportions.
#       - "svyCI" uses the proportions and SE from "se.multinom" to calculate
#          confidence intervals on estimated proportions (using score intervals
#          for a Binomial distribution rather than the (less reliable) normal
#          approximation)
svymultinom = function(src.dat,
                       mysvy,
                       fmla = formula("as.numeric(as.factor(K_US)) ~ week + as.factor(HHS)")) {
  # Arguments:
  #  ~  src.dat: source data frame
  #  ~  mysvy:   survey design object
  #  ~  fmla:    multinomial model formula

  # Output: list of 6 objects:
  #  mlm:       nnet multinomial model object (edited to include "svyvcov", which
  #             is the variance-covariance matrix after accounting for survey
  #             design (i.e. "sandwich" below))
  #  estimates: coefficient estimates from the mlm model (named numeric vector)
  #  scores:    gradient of the loglikelihood * survey weights * model matrix
  #             (numeric matrix with row for each observation and column for each coefficient)
  #  invinf:    variance-covariance matrix for coefficients of the mlm model (i.e.
  #             not accounting for survey design)
  #             (numeric matrix with row and column for each coefficient)
  #  sandwich:  variance-covariance matrix for coefficients after accounting for
  #             survey design
  #             (numeric matrix with row and column for each coefficient)
  #  SE:        SE of coefficients after accounting for survey design
  #             (named numeric vector)


  # aggregate data before fitting the multinomial model to vastly improve run time
  # see: https://git.biotech.cdc.gov/sars2seq/sc2_proportion_modeling/-/issues/6#note_86544
  fmla.vars = all.vars(fmla)
  mlm.dat = data.table::data.table(cbind(data.frame(src.dat)[, fmla.vars],
                                         weight = weights(mysvy)))[
                                           ,
                                           .(weight = sum(weight)), # aggregate "weight" column
                                           by = fmla.vars] # by the formula

  # Fit multinomial logistic regression
  # (without survey design, but with survey weights)
  multinom_geoid = nnet::multinom(formula = fmla,
                                  data    = mlm.dat,
                                  weights = weight,
                                  Hess    = TRUE,
                                  maxit   = 1000,
                                  trace   = FALSE)

  ## Format results to fit into svymle-like workflow
  # get the number of variants being modeled
  # (should be length of model_vars plus 1 for others)
  num_var = length(unique(with(src.dat,
                               eval(terms(fmla)[[2]]))))

  # generate the model formula without the response term
  fmla.no.response = formula(delete.response(terms(fmla)))

  # repeat the model formula w/o response term for each variant listed in model_vars
  formulas = rep(list(fmla.no.response),
                 num_var - 1)

  # Add response back to first formula to ensure inclusion as response later on
  formulas[[1]] = fmla

  # sets the names of the formulas to correspond to beta coefficients for each variant
  names(formulas) = paste0("b", 1:length(formulas) + 1)

  # creates a list that contains the mlm object and the coefficient estimates
  # (this will form the output of this function after other things are added on)
  rval = list(mlm = multinom_geoid,
              estimates = coefficients(multinom_geoid))

  # transforms the estimates object to be a list where each element is a vector
  #  of the coefficients for a given variant
  # (hhs model: Intercept, week, HHS 2:10;
  #   us model: Intercept, week)
  rval$estimates = as.list(data.frame(t(rval$estimates)))
  ## End multinomial regression


  # if the Hessian is not invertible, avoid errors by not running the rest of the function
  # (note: the se.multinom function will still run on the output of this function,
  #  but will assume the SE is 0)
  # if(det(multinom_geoid$Hessian) == -Inf) {
  if(det(multinom_geoid$Hessian) < -9e100) {
    # might want to replace -Inf with a very large negative number: e.g. -9e100

    # add empty items to the results to match structure of a successful run
    # (setting to NULL inside of list() *does* create the named item)
    rval <- append(rval,
                   list(
                     'scores' = NULL,
                     'invinf' = NULL,
                     'sandwich' = NULL,
                     'SE' = NULL
                   ))

    # remove the Hessian
    # (setting to NULL outside of list() *removes* item)
    rval$mlm$Hessian = NULL

    # make the estimates a vector instead of a list
    rval$estimates = unlist(rval$estimates)

    # add prettier names to the estimates (to match names that are produced below)
    names(rval$estimates) <- paste0('b',
                                    sub(pattern = ':',
                                        replacement = '\\.',
                                        x = colnames(multinom_geoid$Hessian)))

    # add a warning
    warning_message = 'Hessian is non-invertible. Nowcast estimates will not have confidence intervals. Check for geographic regions with very few samples of a variant.'
    warning(warning_message)
    # also print out an alert (not really necessary)
    print(warning_message)

  } else {

    ## Define likelihood

    # Loglikelihood:
    lmultinom = function(v, ...) {
      # vectorized loglikelihood function
      # Arguments: v  (positive integer vector) is position of response variable
      #                in ordered list of possible values (1=reference)
      #            ...  vectors of  linear predictors

      # create a matrix of linear predictors (with 0's for reference class)
      b = cbind(0, ...)
      b = matrix(b,
                 nrow = length(v))

      # get the predicted probability of the realized (actual) outcomes and
      # subtract the total probability for all classes???
      sapply(X = seq_along(v),
             FUN = function(rr) b[rr, v[rr]])  - log(rowSums(exp(b)))
    }

    # Gradient of loglikelihood:
    gmultinom = function(v, ...) {
      # vectorized partial derivatives of lmultinom(v, ...) with respect to
      #  linear predictors at ...
      # vectorized loglikelihood function
      # Arguments: v  (positive integer vector) is position of response variable
      #                in ordered list of possible values (1=reference)
      #            ...  vectors of linear predictors

      # create a matrix of linear predictors (with 0's for reference class)
      #  column for each observation
      #  row for each possible outcome class
      b = cbind(0, ...)
      b = matrix(b,
                 nrow=length(v))

      # convert to probabilities by normalizing linear predicted values by the
      # total value across all possible outcome classes
      p = exp(b) / rowSums(exp(b))

      # create a matrix of the same dimensions
      # to hold the observed class of each observation
      delta_ij = p * 0

      # set some values of the matrix to 1
      # (for each observation (i.e. row) assign a value of 1 to the column
      #  representing its clade/rank order)
      delta_ij[1:length(v) + (v-1) * length(v)] = 1

      # Column n is partial derivative with respect to ...[n]:
      # (realized outcome minus the predicted probability of belonging to that class)
      (delta_ij - p)[, -1]
    }
    ## End likelihood definitions

    ## Build dataframe to pass to svyrecvar
    # Adapted from code in svymle
    # modify formula names
    nms = c("", names(formulas))

    # logical that identifies which formula has a response variable
    has.response = sapply(formulas, length) == 3

    # creates variable equal to regression formula that has response variable
    ff = formulas[[which(has.response)]]

    #sets predictor terms to 1 so: as.numeric(as.factor(K_US)) ~ 1
    ff[[3]] = 1

    #I think this just gets you what the rankings are for the variants from src.dat
    # (i.e. the response data for the regression model)
    y = eval.parent(model.frame(formula   = ff,
                                data      = src.dat,
                                na.action = na.pass))

    # sets the first formula to a formula without the response term
    formulas[[which(has.response)]] = formula(delete.response(terms(formulas[[which(has.response)]])))

    # creates empty list the same length of the formulas object
    # mf = vector("list", length(formulas))
    # object "mf" gets redefined below. This line doesn't seem to do anything.

    # vector with the names of the regression predictors
    vnms = unique(do.call(c, lapply(formulas, all.vars)))

    # converts the vnms to a formula object w/o the response term
    uformula = make.formula(vnms)

    # gets the values for the model predictors from src.dat
    mf = model.frame(formula   = uformula,
                     data      = src.dat,
                     na.action = na.pass)

    # adds a column for the variant rankings to mf with the column header as
    #  as.numeric(as.factor(K_US))
    mf = cbind(`(Response)` = y,
               mf)
    # this isn't changing/setting the column header to "(Response)"...

    # I think this just ensures that there aren't any duplicated columns
    #  (the drop = FALSE ensures that mf remains a dataframe even if only 1 column
    #   is returned)
    mf = mf[, !duplicated(colnames(mf)),
            drop = FALSE]

    # get the svy weights from the survey design object
    weights = weights(mysvy)

    # defines Y (response variable) as variant group (rank abundance)
    Y = mf[, 1]

    # for each formula listed in "formulas", generate the underlying data for
    # those parameters
    # (model.frame object is required for the model.matrix function)
    mmFrame = lapply(X    = formulas,
                     FUN  = model.frame,
                     data = mf) # argument passed to "model.frame"

    # creates model design objects for each rank class
    mm = mapply(FUN      = model.matrix,
                object   = formulas,
                data     = mmFrame,
                SIMPLIFY = FALSE)

    # get the cumulative number of columns across all dataframes listed in mm
    np = c(0,
           cumsum(sapply(X   = mm,
                         FUN = NCOL)))

    # get the column names from all the dataframes in mm
    parnms = lapply(mm, colnames)

    # format column names in each list element to correspond to the coefficient
    for (i in 1:length(parnms)){
      parnms[[i]] = paste(nms[i + 1],
                          parnms[[i]],
                          sep = ".")
    }

    # convert the parameter names from a list to a vector
    parnms = unlist(parnms)

    # get the estimated coefficients from multinom as a vector
    theta = unlist(rval$estimates)

    # creates empty list the same length as the number of variants plus "other"
    args = vector("list", length(nms))

    # sets first list element to the variant rank data (i.e. response variable)
    args[[1]] = Y

    # assigns the names of the args list to nms
    names(args) = nms

    # get the linear predictor for each response variable class (rank abundance)
    # by multiplying the covariate values in "mm" by the coefficient values, "theta"
    for (i in 2:length(nms)){
      args[[i]] = drop(mm[[i - 1]] %*% theta[(np[i - 1] + 1):np[i]])
    }

    # this takes the args list and uses the gmultinom function above to estimate
    # the gradient of the log likelihood
    deta = matrix(data = do.call(what = "gmultinom",
                                 args = args),
                  ncol = length(args) - 1)

    # creates an empty list element called scores in the rval list
    rval <- append(rval,
                   list('scores' = NULL))

    # creates numerical vector the same length as the names of formulas
    reorder = na.omit(match(x = names(formulas),
                            table = nms[-1]))

    #for each variant this multiplies the gradient of the loglikelihood by the survey weights
    for (i in reorder){
      rval$scores = cbind(rval$scores,
                          deta[, i] * weights * mm[[i]])
    }

    # not 100% sure what this step does. It looks like it's solving the inverse
    # (negative?) multinomial regression object and adds as a list element to rval
    rval$invinf = solve(-multinom_geoid$Hessian)
    # when provided a single matrix, the "solve" function will invert it. So here
    # it's inverting the negative of the Hessian.
    # inverse of the negative Hessian is the estimate of the variance-covariance
    # matrix of model parameters (which is used to estimate parameter SE)

    # assign names to the rval$invinf
    dimnames(rval$invinf) = list(parnms, parnms)

    # Use matrix multiplication to multiply the "scores" (i.e., the gradient
    # loglikelihood * survey weights * model matrix) by the variance-covariance
    # matrix of parameter estimates
    db = rval$scores %*% rval$invinf

    # Everything up until now has been formatting steps to get the estimates/data
    # formatted in way that's compatible with the svyrecvar function

    # Computes the variance of a total under multistage sampling, using a recursive
    # descent algorithm. The object that's returned ("sandwich") is the covariance matrix
    rval$sandwich = survey::svyrecvar(x          = db,
                                      clusters   = mysvy$cluster,
                                      stratas    = mysvy$strata,
                                      fpcs       = mysvy$fpc,
                                      postStrata = mysvy$postStrata)

    # estimate the SE as the square root of the diagonal elements of the
    # variance-covariance matrix (diagonal = variance)
    rval$SE = sqrt(diag(rval$sandwich))

    # convert the "estimates" object from a list to a vector
    rval$estimates = unlist(rval$estimates)

    # assign the names of "estimates" in rval to be the same as the names of "SE"
    names(rval$estimates) = names(rval$SE)

    # adds the variance-covariance matrix to the mlm object (mlm being the results
    # object you get from solving the multinomial model)
    rval$mlm$svyvcov = rval$sandwich
  }

  # return the rval object
  rval
}


# function to extract the predictions and SE of predictions from the multinomial
# nowcast model for a particular week and location
# (and optionally grouped clades/variants)
# Note: This function estimates proportions (predicted values) and their standard
#       errors using the results of "svymultinom"; the "svymultinom" function
#       fits the model and estimates the SE of model coefficients.
se.multinom = function(mlm,
                       week,
                       geoid = "USA",
                       composite_variant = NULL) {
  # Arguments:
  #  ~  mlm:   model output, with Hessian;
  #  ~  week:  focal week to calculate the SE for
  #  ~  geoid: geographic region; options include "USA" or HHS Region (1:10)
  #  ~  composite_variant: (if not NA) is a matrix with one column per model
  #                        lineage and one row per composite variant matrix
  #                        element of 1 marks each component lineage (column)
  #                        for each composite variant (row).
  #                        Example: matrix(c(1, 0, 0, 1, 0, 0), nrow = 1)
  #                        designates a variant comprising the first and fourth
  #                        lineages in the model estimates for composites will
  #                        be appended at end of p_i and se.p_i

  # If mlm is without Hessian, all se's are set to zero
  # mlm not geographically stratified for geoid="USA": ~ (week - current_week)
  # mlm, for HHS Regions: ~ week + HHS (1 as reference level)

  # Output: list of 5 items:
  #  ~  p_i:     predicted values for each clade/variant (numeric vector)
  #  ~  se.p_i:  SE of predicted values (numeric vector)
  #  ~  b_i:     coefficient (beta) values (numeric vector)
  #  ~  se.b_i:  SE of coefficients (numeric vector)
  #  ~  composite_variant: list with 3 element:
  #                        matrix = composite_variant (input matrix)
  #                        p_i    = predicted proportions for composite variants
  #                        se.p_i = SE of the predicted proportions

  # get the model coefficients
  cf = coefficients(mlm)

  # get the variance-covariance matrix
  if ("svyvcov" %in% names(mlm)) {
    vc = mlm$svyvcov
  } else {
    if ("Hessian" %in% names(mlm)) {
      vc = solve(mlm$Hessian)
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      # Why use solve(hessian) instead of solve(-hessian)?
      #   svymle(), through a call to nlm(), minimizes a function. The functions
      #   lmultinom() and gmultinom() are, therefore, negative log-likelihood
      #   and it's gradient, respectively. The Hessian is, therefore, negative
      #   of what you'd expect.
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    } else {
      # If there's no variance-covariance matrix from the svymultinom function
      # and there's no Hessian, just set the variance-covariance to 0
      vc = rep(0, length(cf)^2)
    }
  }

  # Rearrange terms to ease pulling out relevant submatrix
  # convert symmetrical matrix with rows & columns for each coefficient:
  # 4-d array with dim 1 for each covariates
  #                dim 2 for each outcome class
  #                dim 3 for each covariate
  #                dim 4 for each outcome class
  dim(vc) = rep(x = rev(dim(cf)),
                times = 2)

  # Boolean: is the geoid the reference level
  ref_geoid = ((length(mlm$xlevels) == 0) || (geoid == mlm$xlevels[[1]][1]))

  # get the number/index of the geoid/area being modeled
  if (ref_geoid) { # if the geoid is the reference level
    n_geoid = NULL
  } else { # if it's not the reference level
    n_geoid = which(mlm$xlevels[[1]] == geoid)
  }

  # Set value for "hhs1tail"
  # index of column identifying column for the hhs region in question
  if (ref_geoid) { # if the geoid is the reference level
    hhs1tail = NULL
  } else { # if it's not the reference level
    hhs1tail = n_geoid + 1
  }

  # save values to vector
  # only want the coefficients for the intercept, the week in question, and the
  # hhs region in question
  indices = c(1, 2, hhs1tail)

  # reset value for "hhs1tail"
  # indicator (0/1) covariate value (1 for the hhs region in question)
  if (ref_geoid) {
    hhs1tail = NULL
  } else {
    hhs1tail = 1
  }

  # save values to another vector
  # (aren't these really covariate values, not coefficients?)
  coeffs = c(1, week, hhs1tail)

  # b is vector of coefficients of time (week)
  # subset the variance-covariance matrix to the relevant parameters
  sub.vc = vc[indices,,indices,]

  # get the linear predictor values
  # (what's the 0 on the beginning for?)
  y_i = c(0, cf[, indices] %*% coeffs)

  # get the variance-covariance matrix of the linear predictor
  # (after adjusting for survey design)
  vc.y_i = outer(1:dim(sub.vc)[2],
                 1:dim(sub.vc)[4],
                 Vectorize(function(i, j) c(coeffs %*% sub.vc[, i, , j] %*% coeffs)))
  # (the "coeffs" here really covariates)

  # add in row and column for the intercept???
  vc.y_i = rbind(0, cbind(0, vc.y_i))

  # calculate the predicted proportions
  p_i = exp(y_i)/sum(exp(y_i))

  # Taylor series based variance:
  # dp_i/dy_j = delta_ij p_i - p_i * p_j
  dp_dy = diag(p_i) - outer(p_i, p_i, `*`)

  # calculate variance/covariance of the predicted proportions
  p.vcov = dp_dy %*% vc.y_i %*% dp_dy
  # SE = sqrt of variance
  se.p_i = as.vector(sqrt(diag(p.vcov)))

  # if there are composite variants, calculate calculate their proportions & SE
  if (!is.null(composite_variant)) {
    composite_variant = list(
      matrix = composite_variant,
      p_i    = as.vector(composite_variant %*% p_i),
      se.p_i = as.vector(sqrt(diag(composite_variant %*% p.vcov %*% t(composite_variant))))
    )
  }

  # output
  list(p_i = p_i, # predicted values
       se.p_i = se.p_i, # SE of predicted values
       b_i = c(0, cf[, 2]), # coefficient values for ?? (intercept & week?)
       se.b_i = c(0, sqrt(diag(vc[2,,2,]))), # SE of the coefficients
       # (but unadjusted for survey design. Why not take these from svymultinom ouput?)
       composite_variant = composite_variant) # predicted probabilities (and SE) for the composite variants
}

# Function to get binomial confidence interval based on the point estimated
# proportion and associated SE
# (based on output from svymultinom & se.multinom functions)
svyCI = function(p, s) {
  # p = point estimate
  # s = estimated standard error

  # calculate the sample size
  n = p*(1-p)/s^2

  # calculate the confidence interval using a proportion test
  out = prop.test(x = n * p,   # a vector of counts of successes
                  n = n        # a vector of counts of trials
  )$conf.int
  # could calculate a normal approximation to the confidence interval as:
  # logit( p +/- 1.96 * SE )
  # prop.test says: "The confidence interval is computed by inverting the score test."
  # Wikipedia (https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval)
  # says: "The Wilson score interval is an improvement over the normal approximation
  # interval in multiple respects... Unlike the symmetric normal approximation
  # interval (above), the Wilson score interval is asymmetric. It does not suffer
  # from problems of overshoot and zero-width intervals that afflict the normal
  # interval, and it may be safely employed with small samples and skewed
  # observations.[3] The observed coverage probability is consistently closer to
  # the nominal value.

  # Is this incorporating the Korn-Graubard method?
  #   Korn and Graubard (1998) have suggested a method for producing confidence intervals for
  #   proportions estimated from a sample based on a complex sample design where the proportions
  #   are either very small or very large, or the sample size is small. Their method uses the exact
  #   binomial confidence intervals but with the sample size modified by dividing by the estimated
  #   design effect for the proportion in question.

  # return the lower and upper confidence interval limits
  c(out[1],out[2])
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
#
## Results ---------------------------------------------------------------------

# Current top variants by share and variants of concern or interest
# (run for all values of "tag")
{
  # Based on sequences where the specimen was collected in the current or previous
  # `n_recent_weeks` weeks, the count of sequences, and variant share (weighted
  #  percent) are:
  #
  # Select display and model variants - this basically selects which variants are
  # the top variants to include in the model
  us_var = sort(
    prop.table(x = xtabs(SIMPLE_ADJ_WT ~ VARIANT,
                         subset(src.dat,
                                week >= current_week - n_recent_weeks &
                                  VARIANT != "None"))),
    decreasing = TRUE)

  # counts of all variants in the US
  us_seq = table(subset(src.dat,
                        week >= current_week - n_recent_weeks)$VARIANT)

  # names of all the variants
  us_rank = names(us_var)

  # counts of variants for each hhs region
  hhs_var = prop.table(
    x = xtabs(SIMPLE_ADJ_WT ~ HHS + VARIANT,
              subset(src.dat,
                     week >= current_week - n_recent_weeks)),
    # should this also have VARIANT != "None"?
    margin = 1)

  # get the names of the variants for each hhs region (in order of abundance)
  hhs_rank = apply(hhs_var, 1, function(rr) names(sort(rr, decreasing=TRUE)))

  # names of variants that are either in the n_top or are in "voc"
  # Ordered by national rank
  # these will be included in results
  all_tops = us_rank[us_rank %in% c(us_rank[1:n_top], voc)]

  # variants to include in the nowcast model
  model_vars = us_var[all_tops]

  # choose which variants to display in the results
  if (display_option=="top7") {
    display_vars = head(names(model_vars), 7)
  } else {
    display_vars = voc
  }

  # get the indices of the variants that will be displayed in the results
  display_indices = which(names(model_vars) %in% display_vars)
  # get the names of the variants to display (in correct order)
  display_vars = names(model_vars)[display_indices]

  #add variant ranks to src.dat
  # makes sure each seq is assigned a rank/number based on the weighted proportion
  # in the last few weeks (1 = most common)
  src.dat$K_US = sapply(X = 1:nrow(src.dat),
                        FUN = function(nn) which( c(names(model_vars), src.dat$VARIANT[nn]) ==
                                                    src.dat$VARIANT[nn])[1])

  # add the current week to the source data
  src.dat$current_week = current_week


  # reported variants will include
  # - voc for short list
  # - all_tops for long
  # (will differ in what falls under "other")
  # (used in "Not Run3" and "Run3" groups)
  reported_variants = voc

  # specify the survey design
  # (used in "Not Run3" and "Run2". It's redefined for "Run3")
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


  # create another column for the varients of interest
  src.dat$VARIANT2 = as.character(src.dat$VARIANT)

  # group all non-"variants of interest" together
  src.dat[src.dat$VARIANT %notin% voc, "VARIANT2"] <- "Other"
}

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

  # create a dataframe with all unique combinations of variants, fortnights, and regions
  # (and then reverse column order for convenience)
  all.ftnt = expand.grid(Variant = reported_variants,
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
                      FUN = function(rr) myciprop(voc = reported_variants,
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
  all.wkly = expand.grid(Variant = reported_variants,
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
                      FUN = function(rr) myciprop(voc = reported_variants,
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
  # produces results for Run1 & Run2


  #### Run 2.1 Fit the model


  #create a subset of src.dat that only contains the weeks that will be included
  # in multinomial model
  src.moddat = subset(src.dat,
                      week >= max(week) - model_weeks) # Modeling window

  #create survey design object based on subsetted src.dat
  mysvy = svydesign(ids     = ~SOURCE,
                    strata  = ~STUSAB + yr_wk,
                    weights = ~SIMPLE_ADJ_WT,
                    nest = TRUE, # TRUE = disable checking for duplicate cluster ID's across strata
                    data = src.moddat)
  # not specifying fpc signifies sampling with replacement

  #runs the svymultinom function and saves the output
  svymlm_hhs = svymultinom(src.dat = src.moddat,
                           mysvy = mysvy,
                           fmla = formula("as.numeric(as.factor(K_US)) ~ week + as.factor(HHS)"))

  # create a dataframe for predictions for each week & HHS region
  pred_hhs.df = expand.grid(week = seq(from = -display_lookback,
                                       to = 2,
                                       by = 0.02) + current_week,
                            HHS = sort(unique(src.moddat$HHS)))

  # predict a sequence's clade (ID'd by rank abundance) given only time & hhs region
  # (i.e. predict the proportion for each clade in each hhs region & time period)
  pred_hhs.df = cbind(pred_hhs.df,
                      predict(object  = svymlm_hhs$mlm,
                              newdata = pred_hhs.df,
                              type    = "probs"))

  # tabulate/count the normalized survey weights by region, week, and variant
  bp_hhs = xtabs(formula = NORM_WTS ~ HHS + week + VARIANT,
                 data = subset(src.moddat,
                               week < current_week &
                                 week >= current_week - display_lookback))

  # subset the weights to only include the variants to be displayed
  bp_hhs = prop.table(bp_hhs, 1:2)[,, display_vars]

  # give the table prettier column names
  dimnames(bp_hhs)[[2]] = week_label(as.numeric(dimnames(bp_hhs)[[2]]) - current_week)

  # Run another model for the entire US
  svymlm_us = svymultinom(src.dat = src.moddat,
                          mysvy = mysvy,
                          fmla = formula("as.numeric(as.factor(K_US))  ~ week"))

  # create a dataframe for predictions for each week
  pred_us.df = expand.grid(week = seq(from = -display_lookback,
                                      to = 2,
                                      by = 0.02) + current_week)

  #add a column for (predicted) each variant proportion for each timepoint
  pred_us.df = cbind(pred_us.df,
                     predict(object = svymlm_us$mlm,
                             newdata = pred_us.df,
                             type = "probs"))

  # tabulate/count the normalized survey weights by week, and variant
  bp_us = xtabs(formula = NORM_WTS ~ week + VARIANT,
                data = subset(src.moddat,
                              week < current_week &
                                week >= (current_week - display_lookback)))

  # Normalize values to proportions
  bp_us = prop.table(bp_us, 1)

  # subset the weights to only include the variants to be displayed
  bp_us = bp_us[, display_vars]

  # give the table prettier column names
  rownames(bp_us) = week_label(as.numeric(rownames(bp_us)) - current_week)


  #### Run 2.2 Plot results


  # Weighted variant shares of the top variants in the past `display_lookback`
  # weeks (number of sequences collected weekly above each bar), and model-based
  # smoothed estimates, nationwide:

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

  ### create figures and (optionally) save to file
  # Barplot of US national data
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

  # Barplot of Nowcast predicted values
  if (fig_gen_run) jpeg(filename  = paste0(stub, "projection_US", tag, ".jpg"),
                        width     = 1500,
                        height    = 1500,
                        pointsize = 40)

  # create a barplot with model-predicted data
  bp = barplot(height = 100 * t(pred_us.df[, 1 + display_indices]),
               xlab = "Week beginning",
               ylab = "Weighted variant share (%)",
               main = "Nationwide",
               space = 0,
               border = NA,
               ylim = 110 * 0:1,
               col = col.dk,
               names.arg = ifelse(pred_us.df$week %% 1 == 0,
                                  week_label(pred_us.df$week - current_week),
                                  NA),
               legend.text = display_vars,
               args.legend = list(x = "topleft",
                                  bty = "n",
                                  border = NA))

  # predicted percent contributions of some variants for the current week
  pc = unlist(100 * subset(pred_us.df,
                           week==current_week)[, 1+display_indices])

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
  x = bp[which(pred_us.df$week %in% (current_week + c(-2, 0, 2)))]

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


  # Get the SE of the national estimate
  us.summary = se.multinom(mlm   = svymlm_us$mlm,
                           week  = current_week,
                           geoid = "USA")

  # calculate the SE of the estimated growth rate
  se.gr = with(data = us.summary,
               expr = 100 * exp(sqrt(se.b_i^2 + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2))) - 100)
  # p_i    = predicted probability (proportional representation/frequency)
  # se.p_i = SE of predicted frequency
  # b_i    = coefficient value (for intercept & week)
  # se.b_i = SE of coefficient value

  # calculate the estimated growth rate
  gr = with(data = us.summary,
            expr = 100 * exp(b_i - sum(p_i * b_i)) - 100)


  # Plot growth rate vs. transmission
  #  - the vertical axis depicts the variant share growth rate
  #    (derivative of log of variant proportion with respect to time).
  #  - The axis on the right shows an estimate of variant transmissibility with
  #    respect to the overall mean transmissibility.
  if (fig_gen_run) png(filename = paste0(stub, "growthrate_US", tag, ".png"),
                       width = 8,
                       height = 8,
                       units = "in",
                       pointsize = 16,
                       res = 1000)

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
    lines(x = pmax(100 * us.summary$p_i[vv] + 196 * c(1,-1) * us.summary$se.p_i[vv]),
          y = rep(gr[vv], 2),
          col = "blue",
          lwd = 2)
  }

  # add the name of each variant
  text(x = 100 * us.summary$p_i,
       y = gr,
       labels = c(names(model_vars),""),
       cex = 0.85,
       col = "grey25",
       adj = 1.15)
  if (fig_gen_run)  dev.off()

  # create a dataframe of variant shares & growth rates
  gr_tab = cbind(variant = c(names(model_vars), "OTHER"),
                 variant_share = (100 * us.summary$p_i),
                 growth_rate = gr)

  # save growth rates to file
  write.csv(x = gr_tab,
            file = paste0(script.basename,
                          "/results/wow_growth_variant_share",
                          data_date,
                          tag,
                          ".csv"),
            row.names = FALSE)


  # Same as above for each HHS region:
  # Weighted variant shares of the top variants in the past `display_lookback`
  # weeks (number of sequences collected weekly above each bar), and model-based
  # smoothed estimates, for each HHS region:

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
                 names.arg = rownames(bp_hhs[hhs,,] - current_week),
                 legend.text = display_vars,
                 args.legend = list(x = "topleft",
                                    bty = "n",
                                    border = NA))

    # add text to the barplot
    text(x = bp,
         y = 3 + colSums(100 * t(tail(bp_hhs[hhs,,], 12))),
         labels = with(subset(src.dat,
                              HHS == hhs &
                                week < current_week &
                                week >= current_week - display_lookback),
                       table(week)),
         cex = 0.7)
    if (fig_gen_run) dev.off()

    # create another barplot for Nowcast predictions
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
                 names.arg = ifelse(pred.df$week %% 1 == 0,
                                    week_label(pred.df$week - current_week),
                                    NA),
                 legend.text = display_vars,
                 args.legend = list(x = "topleft",
                                    bty = "n",
                                    border = NA))

    # predicted percent contributions of some variants for the current week
    pc = unlist(100 * subset(pred.df, week==current_week)[, 1+display_indices])

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
    x = bp[which(pred.df$week %in% (current_week + c(-2, 0, 2)))]

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


    # create plot of growth rates by region
    if (fig_gen_run) jpeg(filename = paste0(stub, "growthrate_HHS", hhs, tag,".jpg"),
                          width = 1500,
                          height = 1500,
                          pointsize = 40)

    # get the SE for the given model, week, and HHS region
    hhs.summary = se.multinom(mlm = svymlm_hhs$mlm,
                              week = current_week,
                              geoid = hhs)

    # calculate the SE of the growth rate
    se.gr = with(hhs.summary,
                 100 * exp(sqrt(se.b_i^2 + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2))) - 100)

    # calculate the growth rate
    gr = with(hhs.summary,
              100 * exp(b_i - sum(p_i * b_i)) - 100)

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
      lines(x = pmax(100 * hhs.summary$p_i[vv] + 1.96 * c(1,-1) * hhs.summary$se.p_i[vv]),
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
         labels = c(names(model_vars),""),
         cex = 0.85,
         col = "grey25",
         adj = 1.15)

    if (fig_gen_run) dev.off()
  }


  #### Run 2.3 Model-smoothed estimates ("updated nowcast")
       # (for runs 1 & 2)


  ## Model-smoothed estimates (needs Hessian for regional multinomial model)
  # (uses confidence intervals from "svymultinom" > "se.multinom" > "svyCI")

  # This section of code is only run for Run 2, but it also produces results for
  # Run 1. Therefore, it needs to aggregate some of the variants in the Nowcast
  # model to match those that are reported for Run 1.

  # Check to see which lineages are in model_vars
  # this returns all variants with "AY" in the name
  AY_agg = names(model_vars)[grep("AY", names(model_vars))]

  # get the names of the lineages to avoid aggregating for Run1
  if(custom_lineages){
    run1_lineages = voc1_custom
  } else {
    run1_lineages = voc1
  }

  # this returns all "AY" variants to be aggregated for Run 1 results
  # (i.e. not listed in run1_lineages)
  AY_agg = AY_agg[AY_agg %notin% run1_lineages]

  # all other variants to be aggregated (used for Run 1 & Run 2 results)
  Other_agg = names(model_vars)[names(model_vars) %notin% voc]

  # generate a matrix that indicates which lineages to aggregate for the nowcast
  # Columns are the lineages in the nowcast model, so all the defined lineages
  #  plus the "other" lineage
  # Rows are the aggregated lineages wanted
  agg_var_mat <- matrix(data = 0,
                        nrow = 2,
                        ncol = (length(model_vars)+1))
  colnames(agg_var_mat) <- c(names(model_vars),"Other")

  # Fill in matrix values: if lineage is to be aggregated to parent lineage in
  # given row, then value = 1, else value = 0
  agg_var_mat[1,] <- ifelse(colnames(agg_var_mat) %in% c("B.1.617.2", AY_agg),1,0)
  agg_var_mat[2,] <- ifelse(colnames(agg_var_mat) %in% c(Other_agg, "Other"),1,0)
  row.names(agg_var_mat) <-c("Delta Aggregated", "Other Aggregated")

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
        mlm = svymlm_us$mlm
        geoid = rgn
      } else {
        mlm = svymlm_hhs$mlm
        geoid = as.numeric(rgn)
      }

      # get the week for the timepoint
      wk = as.numeric(as.Date(ftn, origin="1970-01-01") - week0day1) %/% 7

      # get the estimates (and SE) for the given place & time
      ests = se.multinom(mlm = mlm,
                         week = wk,
                         geoid = geoid,
                         composite_variant = agg_var_mat)

      # format estimates into dataframe with relevant info
      ests = data.frame(
        USA_or_HHSRegion = rgn,
        Fortnight_ending = as.Date(ftn, origin="1970-01-01"),
        Variant = c(names(model_vars),
                    "Other",
                    row.names(ests$composite_variant$matrix)),
        Share = c(ests$p_i,
                  ests$composite_variant$p_i),
        se.Share = c(ests$se.p_i,
                     ests$composite_variant$se.p_i)
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
    }
  }

  # select just the columns we want
  proj.res = proj.res[, c("USA_or_HHSRegion",
                          "Fortnight_ending",
                          "Variant",
                          "Share",
                          "Share_lo",
                          "Share_hi")]

  # Format output for the run 1 lineage list
  # only include Variants that are NOT in the list provided
  run_1 = proj.res[proj.res$Variant %notin% c(AY_agg,
                                              "B.1.617.2",
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


  ## Weekly estimates
  # same as above, but for weekly estimates instead of fortnightly

  # define and use a function to get the end-of-week dates
  cast_wks = (function(dd) as.Date(seq(from = dd[1],
                                       to = dd[2],
                                       by = 7),
                                   origin = "1970-01-01"))(range(as.Date(proj_ftnts)) + c(-7, 0))

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
      wk = as.numeric(as.Date(cwk, origin="1970-01-01") - week0day1) %/% 7

      # get the estimates (and SE) for the given place & time
      ests = se.multinom(mlm = mlm,
                         week = wk,
                         geoid = geoid,
                         composite_variant = agg_var_mat)

      # format estimates into dataframe with relevant info
      ests = data.frame(
        USA_or_HHSRegion = rgn,
        Week_ending = as.Date(cwk,
                              origin="1970-01-01"),
        Variant = c(names(model_vars),
                    "Other",
                    row.names(ests$composite_variant$matrix)),
        Share = c(ests$p_i,
                  ests$composite_variant$p_i),
        se.Share = c(ests$se.p_i,
                     ests$composite_variant$se.p_i)
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
  proj.res = proj.res[, c("USA_or_HHSRegion",
                          "Week_ending",
                          "Variant",
                          "Share",
                          "Share_lo",
                          "Share_hi")]

  # Format output for the run 1 lineage list
  # only include variants that are NOT in the list provided
  run_1 = proj.res[proj.res$Variant %notin% c(AY_agg,
                                              "B.1.617.2",
                                              Other_agg,
                                              "Other"),]

  # change the name of "Other Aggregated" to "Other" to match other output files
  run_1[run_1$Variant == "Other Aggregated", "Variant"] <- "Other"

  # save the results to file
  write.csv(x = run_1,
            file = paste0(script.basename,
                          "/results/updated_nowcast_weekly_",
                          data_date,
                          sub(pattern = '2', replacement = '1', x = tag),
                          ".csv"),
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

  # save the results to file
  write.csv(x = run_2,
            file = paste0(script.basename,
                          "/results/updated_nowcast_weekly_",
                          data_date,
                          tag,
                          ".csv"),
            row.names = FALSE)
} # end Run2

# Run3 - State-level estimates; Rolling 4 wk bins ------------------------------
# saves a single file:
# - "/results/state_weighted_roll4wk_",ci.type,"CI_svyNEW_",data_date,tag,".csv"
if ( grepl("Run3", tag) ){

  # empty value to hold week numbers
  data_week = c()

  # get the week number that corresponds to each date defined in state_time_end
  for(i in 1:length(state_time_end)){
    data_week[i] <-unique(src.dat$week[src.dat$yr_wk == state_time_end[i]-6])
    # have to subtract 6 days because the yr_wk variable defines week starting
    # on Sunday whereas the state_time_end is defined as week ending Saturday
  }

  # redefine the survey design to the new design for state-level estimates
  svyDES = svydesign(ids     = ~ SOURCE,
                     strata  = ~ STUSAB + yr_wk,
                     weights = ~ SIMPLE_ADJ_WT,
                     nest = TRUE,
                     data = src.dat)

  # all combinations of states, 4-wk periods, & variants
  # (reverse column order for convenience)
  all.state = expand.grid(Variant = reported_variants,
                          Roll_4wk_end = data_week,
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
                       Roll_4wk_end = data_week,
                       State = sort(unique(src.dat$STUSAB)))[, 3:1]

  # calculate the estimated proportions (and CI) using survey design
  ests.others = apply(X = others,
                      MARGIN = 1,
                      FUN = function(rr) myciprop(voc = reported_variants,
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
  for(i in 1:length(data_week)){

    # subset the data to the relevant time period
    dat2 <- src.dat[src.dat$week >= (as.numeric(data_week[i]) - 3) &
                      src.dat$week < (as.numeric(data_week[i]) + 1),]

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
    all.state2 <- merge(x = all.state[all.state$Roll_4wk_end == data_week[i],],
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
    all.state2$Roll_Fourweek_ending <- unique(as.Date(src.dat$yr_wk[src.dat$week==data_week[i]]) + 6)

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


# capture system time
tend = proc.time()

# calculate amount of time code took to run
tdiff = tend[3] - tstart[3]

# print out the runtime
sprintf(paste0("code for ",tag, " run took %f minutes"), round(tdiff/60, 1))