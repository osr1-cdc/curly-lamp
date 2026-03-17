# Functions that are used by the "weekly_variant_report_nowcast.R" script.
# Formerly just the "svycipopkg" function was in its own file. Philip Shirk
# moved the rest of the functions here on 2022-01-07.

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
                        level = level,
                        ...)
    )
  } else {
    # use a modified version of svyciprop to limit effective sample size so that
    # it's never > observed sampled size.
    res = suppressWarnings(
      svycipropkg(~VOC,
                  design = srv_design,
                  level = level,
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
svymultinom = function(mod.dat,
                       mysvy,
                       fmla = formula("as.numeric(as.factor(K_US)) ~ model_week + as.factor(HHS)"),
                       model_vars) {
  # Arguments:
  #  ~  mod.dat: source data frame
  #  ~  mysvy:   survey design object
  #  ~  fmla:    multinomial model formula
  #  ~  model_vars: vector of variants that are included in the model (in same order as K_US).
  #                 This is only used for helping to troubleshoot model-fit issues.

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


  # get the variant names in the model
  modvars <- data.frame(
    'K_US' = 1:(length(model_vars)+1),
    'Variant' = c(model_vars, 'Other')
  )
  moddatvars <- mod.dat[,
                        .(.N,
                          sum(SIMPLE_ADJ_WT)),
                        by = 'K_US']
  names(moddatvars) = c('K_US', 'N', 'Weight')
  modvars <- merge(modvars, moddatvars, all.x = TRUE, by = 'K_US')


  # aggregate data before fitting the multinomial model to vastly improve run time
  # see: https://git.biotech.cdc.gov/sars2seq/sc2_proportion_modeling/-/issues/6#note_86544
  fmla.vars = all.vars(fmla)
  mlm.dat = data.table::data.table(cbind(data.frame(mod.dat)[, fmla.vars],
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
  num_var = length(unique(with(mod.dat,
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

  # check if the Hessian is invertible
  invinf <- tryCatch(
    {
      solve(-multinom_geoid$Hessian)
    },
    error = function(cond) {
      return(NA)
    }
  )

  # if the Hessian is not invertible, avoid errors by not running the rest of the function
  # (note: the se.multinom function will still run on the output of this function,
  #  but will assume the SE is 0)
  if( is.na(invinf[1]) ){ # could use length(invinf) == 1

    # look through data if hessian is not invertible
    if(FALSE){
      # if there are problems with the Hessian, investigate counts by region
      model_counts <- mod.dat[,sum(count),by = c('HHS', 'K_US')]
      model_counts$count = model_counts$V1
      model_counts$V1 = NULL
      model_counts = model_counts[order(model_counts$count, decreasing = FALSE), ]
      model_counts$Variant = modvars$Variant[model_counts$K_US]
      model_counts[,c('HHS', 'K_US', 'Variant', 'count')]


      # also investigate small values in diagonal of Hessian
      head(sort(diag(multinom_geoid$Hessian), decreasing = FALSE))
      head(sort(diag(multinom_geoid$Hessian), decreasing = TRUE))
      modvars

      # add on some more rows to see if that helps
      # mod.dat0 = mod.dat
      extrarows <- mod.dat[ (mod.dat$HHS == 8 & mod.dat$K_US == 20),]
      extrarows$model_week = 9
      extrarows$SIMPLE_ADJ_WT = 10
      mod.dat = rbind(mod.dat,
                      extrarows)

      # update mysvy and rerun the model
      mysvy = svydesign(ids     = ~SOURCE,
                        strata  = ~STUSAB + yr_wk,
                        weights = ~SIMPLE_ADJ_WT,
                        nest = TRUE, # TRUE = disable checking for duplicate cluster ID's across strata
                        data = mod.dat)
    }

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
    # print(warning_message)

    # print out regions with the fewest observations
    # if there are problems with the Hessian, investigate counts by region
    model_counts <- mod.dat[,sum(count),by = c('HHS', 'K_US')]
    model_counts$count = model_counts$V1
    model_counts$V1 = NULL
    model_counts = model_counts[order(model_counts$count, decreasing = FALSE), ]
    model_counts$Variant = modvars$Variant[model_counts$K_US]
    print('Here are counts by region:')
    print(model_counts[,c('HHS', 'K_US', 'Variant', 'count')])

    # print out the highest and lowest values in the Hessian
    print('Also investigate very small or very large values in the Hessian')
    # also investigate small values in diagonal of Hessian
    hess_headtail <- data.frame('element' = names(sort(diag(multinom_geoid$Hessian))),
                                'value' = sort(diag(multinom_geoid$Hessian)))
    rownames(hess_headtail) <- 1:nrow(hess_headtail)
    print(hess_headtail[c(1:5, (nrow(hess_headtail)-4):nrow(hess_headtail)),])

  } else {
    # if the Hessian is invertible, proceed:

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

    #I think this just gets you what the rankings are for the variants from mod.dat
    # (i.e. the response data for the regression model)
    y = eval.parent(model.frame(formula   = ff,
                                data      = mod.dat,
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

    # gets the values for the model predictors from mod.dat
    mf = model.frame(formula   = uformula,
                     data      = mod.dat,
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
    rval$invinf = invinf # solve(-multinom_geoid$Hessian)
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
  rval$variants = modvars
  return(rval)
}


# function to extract the predictions and SE of predictions from the multinomial
# nowcast model for a particular week and location
# (and optionally grouped clades/variants)
# Note: This function estimates proportions (predicted values) and their standard
#       errors using the results of "svymultinom"; the "svymultinom" function
#       fits the model and estimates the SE of model coefficients.
# If including time as a predictor variable, make it the first covariate!
se.multinom = function(mlm, 
                      newdata_1row,
                      composite_variant=NA,
                      dy_dt=data.frame(model_week=1)) { 
  # Arguments:
  #  ~  mlm:   model output, with Hessian;
  #  ~  newdata_1row:  1-row dataframe consistent with predictors in formula
  #  ~  composite_variant: (if not NA) is a matrix with one column per model
  #                        lineage and one row per composite variant matrix
  #                        element of 1 marks each component lineage (column)
  #                        for each composite variant (row).
  #                        Example: matrix(c(1, 0, 0, 1, 0, 0), nrow = 1)
  #                        designates a variant comprising the first and fourth
  #                        lineages in the model estimates for composites will
  #                        be appended at end of p_i and se.p_i
  # ~  dy_dt:  data.frame with a column sharing the same name as the time-based 
  #            predictor variable (typically "model_week") and the number of 
  #            time units to use for calculating growth rates.

  # If mlm is without Hessian, all se's are set to zero
  # mlm not geographically stratified for geoid="USA": ~ (week - current_week)
  # mlm, for HHS Regions: ~ week + HHS (1 as reference level)

  # Output: list of 14 items:
  #  ~  p_i:     predicted values for each clade/variant (numeric vector)
  #  ~  se.p_i:  SE of predicted values (numeric vector)
  #  ~  b_i:     coefficient (beta) values (numeric vector)
  #              NOTE! b_i is only interpretable as growth rates if first
  #                    predictor term is time (eg I((date - Sys.Date())/7)
  #  ~  se.b_i:  SE of coefficients (numeric vector)
  #  ~  composite_variant: list with 3 element:
  #                        matrix = composite_variant (input matrix)
  #                        p_i    = predicted proportions for composite variants
  #                        se.p_i = SE of the predicted proportions
  #  ~ y_i:      linear predictor on the link (logit) scale
  #  ~ vc_y_i:   variance-covariance matrix of the linear predictor on the link (logit) scale
  #
  #  ~ S:        S is matrix of all variants, simple and composite; one row per variant, one column per modeled variant
  #  ~ p_S:      predicted values for each clade/variant, simple and composite (numeric vector)
  #  ~ se.p_S:   SE of predicted values, simple and composite (numeric vector)
  #  ~ g_S:      growth rates (simple and composite variants) on the link (logit) scale
  #  ~ se.g_S:   SE of the growth rates (simple and composite variants) on the link (logit) scale
  #  ~ p.vcov:   variance-covariance matrix of p_S (on the link (logit) scale)
  #  ~ g.vcov:   variance-covariance matrix of g_S (on the link (logit) scale)
  #

  # modification history:
  # Modified to handle 1-row dataframe consistent with predictors in formula (similar to interface for predict) [2022-01-21]
  # Modified to handle composite variants (eg. Delta = combination of multiple lineages in multinomial model) [2021-09-29]
  # Modified to handle output from svymultinom [2021-08-15]
  # Modified to integrate variance estimation for growth rates (including for composites) [2023-10-16]
 
  # get the model coefficients
  cf <- coefficients(mlm)
  # get the variance-covariance matrix
  if ("svyvcov" %in% names(mlm)) {
    vc <- mlm$svyvcov
  } else { 
    if ("Hessian" %in% names(mlm)) {
      vc <- solve(mlm$Hessian)
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
      vc <- matrix(data = 0,
                  nrow = length(cf),
                  ncol = length(cf))
    }
  }
  # get the model matrix from the model fit & the new data
  mm <- model.matrix(as.formula(paste("~",as.character(mlm$terms)[3])),
                    newdata_1row,
                    xlev=mlm$xlevels)
  # umm is set up for growth-rates [2023-10-11]
  umm <- mm * 0
  for (nm in names(dy_dt)) umm[, nm] = dy_dt[, nm]
  # get the covariate names for each variant
  mnames = outer(X = mlm$lev[-1],
                 Y = colnames(mm),
                 FUN = paste, sep=":")
  # vc = vcov(mlm)
  # matrix of covariate values (first create empty matrix)
  cmat = matrix(data = 0,
                nrow = nrow(mnames),
                ncol = ncol(vc),
                dimnames = list(mlm$lev[-1], colnames(mlm$Hessian)))
  # fill in covariate values into the cmat
  for (rr in 1:nrow(cmat)) cmat[rr, mnames[rr,]] = c(mm)
  # ucmat is set up for growth-rates [2023-10-11]
  ucmat = cmat * 0
  for (rr in 1:nrow(ucmat)) ucmat[rr, mnames[rr,]] = c(umm)
  cmat = rbind(`1`=0, cmat) # pad with 0s for reference variant 1
  ucmat = rbind(`1`=0, ucmat) # pad with 0s for reference variant 1
  # linear predictor values
  y_i = c(0, coefficients(mlm) %*% c(mm))
  u_i = c(0, coefficients(mlm) %*% c(umm))
  # get the variance-covariance matrix of the linear predictor
  # (after adjusting for survey design)
  vc.y_i = cmat %*% vc %*% t(cmat)
  vc.z_i = rbind(cmat, ucmat) %*% vc %*% t(rbind(cmat, ucmat))
  # calculate the predicted proportions
  p_i = exp(y_i)/sum(exp(y_i))
  g_i = u_i - sum(p_i * u_i)
  # Taylor series based variance:
  dp_dy = diag(p_i) - outer(p_i, p_i, `*`)
  # calculate variance/covariance of the predicted proportions
  p.vcov = dp_dy %*% vc.y_i %*% dp_dy
  S = diag(nrow = length(y_i))
  if (all(!is.na(composite_variant))) S = rbind(S, composite_variant)
  p.vcov = S %*% p.vcov %*% t(S) # add in composites
  se.p_S = as.vector(sqrt(diag(p.vcov)))
  p_S = as.vector(S %*% p_i)
  w_S = t(t(S) * p_i)/p_S # one row per variant (includes composites), rowsums are 1
  g_S = as.vector(w_S %*% g_i)
  dg_dz = matrix(c(-p_i * g_i, -p_i), nrow=nrow(S), ncol=nrow(vc.z_i), byrow=TRUE) + cbind(w_S * outer(-g_S, g_i, `+`), w_S)
  g.vcov = (dg_dz %*% vc.z_i %*% t(dg_dz))
  se.g_S = as.vector(sqrt(diag(g.vcov)))
  # if there are composite variants, calculate calculate their proportions & SE
  if (all(!is.na(composite_variant))) {
    composite_variant = list(matrix = composite_variant, 
                             p_i = p_S[-(1:length(y_i))],
                             se.p_i = se.p_S[-(1:length(y_i))])
  }
  dim(vc) = rep(rev(dim(cf)), 2) # To pull out coefficient of time...
  res_old = list(p_i=p_i,
                  se.p_i=se.p_S[1:length(y_i)],
                  b_i=c(0, cf[, 2]),
                  se.b_i= c(0, sqrt(diag(vc[2,,2,]))),
                  composite_variant=composite_variant,
                  y_i=y_i,
                  vc.y_i=vc.y_i)
  # updated output
  # S describes all variants (including composites), p_S are proportions, g_S are growth rates, se are standard deviations, vcov are covariance matrices
  res_new = list(S=S,
                  p_S=p_S,
                  se.p_S=se.p_S,
                  g_S=g_S,
                  se.g_S=se.g_S,
                  p.vcov=p.vcov,
                  g.vcov=g.vcov)
  c(res_old, res_new)
}

# Function to get binomial confidence interval based on the point estimated
# proportion and associated SE
# (based on output from svymultinom & se.multinom functions)
svyCI = function(p, s, ...) {
  # p = point estimate
  # s = estimated standard error
  # ... optional arguments passed on to prop.test (e.g. conf.level = 0.95)

  # if se is 0, n will be Inf (possibly because of non-invertible Hessian); return CI of 0,0
  if (s == 0) {
    return(c(0, 0)) # return confidence interval of [0,0]
  } else if (p == 0) {
    # if p == 0 or p == 1, then n will be 0 & prop.test will throw an error
    return(c(0, 0)) # return confidence interval of [0,0]
  } else if (p == 1) {
    return(c(1, 1)) # return confidence interval of [1,1]
  } else {
    # calculate the sample size
    n = p * (1 - p) / s^2

    # calculate the confidence interval using a proportion test
    out = prop.test(
      x = n * p, # a vector of counts of successes
      n = n, # a vector of counts of trials
      ...
    )$conf.int

    ### UPDATE ####
    # prop.test sometimes throws warnings. We should at least produce a flag for the estimates that are suspicious!




    # Is this incorporating the Korn-Graubard method?
    #   Korn and Graubard (1998) have suggested a method for producing confidence intervals for
    #   proportions estimated from a sample based on a complex sample design where the proportions
    #   are either very small or very large, or the sample size is small. Their method uses the exact
    #   binomial confidence intervals but with the sample size modified by dividing by the estimated
    #   design effect for the proportion in question.

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

    # return the lower and upper confidence interval limits
    return(c(out[1], out[2]))
  }
}



# Below is what was formerly in the "svycipropkg.R" script

#NOTE: the svyciprop function from the survey package uses the calculated
#     effective sample size for the KG CI (which can be larger than the actual sample size)
#     The svycipropkg function caps the effective sample size at the actual sample size
#CAVEATE: the svycipropkg code hasn't gone through comprehensive testing.

# svycipropkg function - based on svyciprop function with "beta" method
# 'design' can be an age standardized object created using 'svystandardize'
# this implementation: Crescent Martin, NCHS
# last updated 6/3/2019

svycipropkg <- function  (formula, design, level = 0.95, df = degf(design), ...) {

  # based on svyciprop function in the survey package with method = "beta"
  m <- eval(bquote(svymean(~as.numeric(.(formula[[2]])), design, ...)))
  rval <- coef(m)[1]
  attr(rval, "var") <- vcov(m)
  alpha <- 1 - level

  # if design is an age-standardized survey design object - custom code not in svyciprop function
  if (!is.null(design[['postStrata']])){
    # error checking
    if (!as.character(design[['call']][[1]]) == "svystandardize"){
      stop("svycipropkg design cannot be a subset of an age-standardized survey design object")
    }

    # create temporary data frame
    design_tmp <- design[['variables']][which(design[['prob']] != Inf), ]
    # get the age-adjustment variable
    ageadjvar <- as.character(design[['call']]$by[[2]])
    # get the age-adjustment population weights
    population <- eval(design[['call']]$population)
    # calculate the weighted proportion and get the sample size for each age adjustment group
    p <- lapply(split(design_tmp, design_tmp[[ageadjvar]]),
                function(x) weighted.mean(x[[as.character(formula[[2]])]],
                                          x[[names(design[['allprob']])]] ) )
    n <- lapply(split(design_tmp, design_tmp[[ageadjvar]]), nrow )
    p <- unlist(p)
    n <- unlist(n)
    # calculate the variance of a simple random sample of the same size, for each age group
    age_var <- p*(1-p)/n
    # normalize the age-adjustment weights to total 1
    pop <- population/sum(population)
    # accumulative the SRS variance over age groups
    varsrs_adj <- sum(pop^2 * age_var)
    # design effect
    deff_adj <- ifelse(sum(p) == 0, 1, vcov(m)/varsrs_adj)
    # adjusted effective sample size
    n.eff <- ifelse(rval == 0, nrow(design_tmp),
                    min(nrow(design_tmp),
                        nrow(design_tmp)/deff_adj * (qt(alpha/2, nrow(design_tmp) - 1)/qt(alpha/2, degf(design)))^2))

  }
  else { # crude estimates (not age-adjusted)
    # effective sample size
    n.eff <- coef(m) * (1 - coef(m))/vcov(m)
    # modification to svyciprop: cap adjusted effective sample size at the actual sample size
    # Under the assumption that the true design effect is >1, though it may be estimated as <1 due to instability of the variance estimator
    # This modification produces different estimated CIs than svyciprop IF the estimated design effect <1
    #n.eff <-     n.eff * (qt(alpha/2, nrow(design) - 1)/qt(alpha/2, degf(design)))^2
    n.eff <- min( n.eff * (qt(alpha/2, nrow(design) - 1)/qt(alpha/2, degf(design)))^2, nrow(design))
  }

  ci <- c(qbeta(alpha/2, n.eff * rval, n.eff * (1 - rval) +  1), qbeta(1 - alpha/2, n.eff * rval + 1, n.eff * (1 - rval)))
  halfalpha <- (1 - level)/2
  names(ci) <- paste(round(c(halfalpha, (1 - halfalpha)) *
                             100, 1), "%", sep = "")
  names(rval) <- paste(deparse(formula[[2]]), collapse = "")
  attr(rval, "ci") <- ci
  class(rval) <- "svyciprop"
  rval
}
