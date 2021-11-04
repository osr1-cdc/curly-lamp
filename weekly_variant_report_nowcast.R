#---------------------------------------------------------------------------------------------------------------------------
#
#       Variant share estimates and nowcasts
#           created by Prabasaj Paul
#
# Updates: ~Survey design updated to have PSU=Source and strata= week and state
#          ~Less reliable CIs are now flagged based on NCHS data presentation standards for proportions
#          ~Code now generates state-level weighted estimates for rolling 4 week time bins
#          ~Multinomial model code updated to formally account for survey design in variance estimation via svyrecvar
#          ~Function to calculate binomial CI for nowcast proportion
#
# Last updated: 10/5/2021 (by Molly Steele)
# Ver: KGCI_29Sep2021                            
library(optparse)

# optparse option list
option_list <- list(
  make_option(c("-r", "--run_number"), type = "character", default= "1",
              help="User Name",
              metavar = "character")
  
)

# parseing options list
opts <- parse_args(OptionParser(option_list = option_list))

# get base directory
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
# correction for if running from base repo interactively
if(length(script.basename) == 0) {
  script.basename = "."
}

# create results dir
dir.create(paste0(script.basename,"/results"))

#capture system time
tstart=proc.time()

# source options from config.R file
source(paste0(script.basename, "/config/config.R"))

# ----------------------------------------------------------------------
library(survey) #package with survey desgin functions
library(nnet) #package with multinomial regression for nowcast
library(data.table) # package for speeding up calculation of simple adjusted weights
options(survey.adjust.domain.lonely=T, survey.lonely.psu="average") #defines what to do when you have a lonely PSU 
options(stringsAsFactors=FALSE)


# Load output from variant_surveillance_system.r ------------------------

#Source the code with svycipropkg function created by Crescent Martin (odb4).
#NOTE: the svyciprop function from the survey package uses the calculated 
#     effective sample size for the KG CI (which can be larger than the actual sample size)
#     The svycipropkg function caps the effective sample size at the actual sample size
#CAVEATE: the svycipropkg code hasn't gone through comprehensive testing.

source(paste0(script.basename, "/svycipropkg.R"))

#load the genomic surveillance data plus survey weights
#load(paste0("svydat_", Sys.Date(), ".RData")) # Works only if svy.dat has been updated on day of report
#load(paste0("svydat_", data_date, ".RData"))
load(paste0(script.basename, "/data/svydat_", data_date, custom_tag, ".RData"))

# Convert any factor to string
fac2str = sapply(svy.dat, class)
fac2str = names(fac2str[fac2str=="factor"])
for (vv in fac2str) svy.dat[, vv] = as.character(svy.dat[, vv])
# Adding new lab streams (that are not handled by the lab-dependent weighting yet) [updated from just NS3 2021-04-02]
#svy.dat$SOURCE = svy.dat$LAB # svy.dat variable redefinition resolves this
#src.dat = subset(svy.dat, SOURCE!="OTHER" & !is.na(VARIANT) & VARIANT!="None") # VARIANT exclusions added 2021-04-21
svy.dat$SOURCE = svy.dat$LAB2 # svy.dat variable redefinition resolves this
svy.dat$count=1
if(state_source=="state_tag_included"){
  src.dat = subset(svy.dat, SOURCE!="OTHER" & !is.na(VARIANT) & VARIANT!="None") # VARIANT exclusions added 2021-04-21
  #check which state-level sources have less than 100 samples
  #check_count <- aggregate(count~SOURCE,data=src.dat,FUN=sum)
  #low_lab <- check_count$SOURCE[check_count$count<100]
  #src.dat <- src.dat[src.dat$SOURCE %notin% c("CDC",low_lab),]
} else {
  src.dat = subset(svy.dat, SOURCE %in% c("UW VIROLOGY LAB","FULGENT GENETICS","HELIX","HELIX/ILLUMINA","LABORATORY CORPORATION OF AMERICA",
                                          "AEGIS SCIENCES CORPORATION","QUEST DIAGNOSTICS INCORPORATED","BROAD INSTITUTE","INFINITY BIOLOGIX",
                                          "NS3","MAKO MEDICAL"))
  src.dat=subset(src.dat,!is.na(VARIANT) & VARIANT!="None")
}

src.dat$sgtf_weights[is.na(src.dat$sgtf_weights) | !src.dat$SGTF_UPSAMPLING] = 1 # svy.dat variable redefinition resolves this
# (Weighted) count of sequences
seq.tbl = with(src.dat, xtabs((1/sgtf_weights) ~ STUSAB + yr_wk))
# for (rr in 1:nrow(src.dat)) {
#   src.dat$SIMPLE_ADJ_WT[rr] = sqrt(src.dat$state_population[rr]/src.dat$TOTAL[rr]) *
#     src.dat$POSITIVE[rr]/seq.tbl[src.dat$STUSAB[rr], src.dat$yr_wk[rr]] / src.dat$sgtf_weights[rr]
#   # Impute by HHS region for states with missing testing data
#   if (is.na(src.dat$SIMPLE_ADJ_WT[rr])) {
#     src.dat$SIMPLE_ADJ_WT[rr] =
#       src.dat$state_population[rr] * src.dat$HHS_INCIDENCE[rr] / seq.tbl[src.dat$STUSAB[rr], src.dat$yr_wk[rr]] / src.dat$sgtf_weights[rr]
#   }
# }


# Replacing for loop with faster data.table code
# see: https://git.biotech.cdc.gov/sars2seq/sc2_proportion_modeling/-/issues/6#note_86543
src.dat = data.table::data.table(src.dat)
src.dat[, "SAW" := sqrt(state_population/TOTAL) * POSITIVE/sum(1/sgtf_weights)/sgtf_weights, .(STUSAB, yr_wk)]
src.dat[, "SAW_ALT" := state_population*HHS_INCIDENCE / sum(1/sgtf_weights)/sgtf_weights, .(STUSAB, yr_wk)]
src.dat[, "SIMPLE_ADJ_WT" := ifelse(is.na(SAW), SAW_ALT, SAW)]
src.dat[, c("SAW", "SAW_ALT") := .(NULL, NULL)]


# Just to be sure:
src.dat = subset(src.dat, !is.na(SIMPLE_ADJ_WT) & SIMPLE_ADJ_WT < Inf) # this shouldn't be necessary, just precautionary
#These are normalized weights so the total weights add up to number of sequences - these weights are used in the multinomial model
src.dat$NORM_WTS = with(src.dat, SIMPLE_ADJ_WT/sum(SIMPLE_ADJ_WT, na.rm=TRUE))
src.dat$NORM_WTS = src.dat$NORM_WTS * sum(!is.na(src.dat$NORM_WTS))

src.dat$FORTNIGHT_END = as.character(week0day1 + src.dat$week%/%2 * 14 + 13)
#Create variable that indicates date ends for 4-week time bins
src.dat$FOURWEEK_END = as.character(week0day1 + src.dat$week%/%4 * 28 + 27)
src.dat$VARIANT = as.character(src.dat$VARIANT) # if saved as factor

#Identify all the AYs in the surveillance database

AY=sort(unique(src.dat$VARIANT)[grep("AY",unique(src.dat$VARIANT))]) 
AY=AY[which(AY %notin% voc)] #vector of the AYs to aggregate
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

#Added per specifc 6/11/2021 request to aggregate sublineages of P1 and B1.351 to 
# the parent lineage (for now taking the easy route of just recoding those lineages)
# [updated by Molly 10/21/21: updated to automatically define list of sublineages to aggregate]
if(P.1_agg==TRUE) {src.dat[src.dat$VARIANT %in% P1,"VARIANT"] <- "P.1"}
if(B.1.351_agg==TRUE) {src.dat[src.dat$VARIANT %in% B351,"VARIANT"] <- "B.1.351"}
if(B.1.621_agg==TRUE) {src.dat[src.dat$VARIANT %in% B621,"VARIANT"] <- 'B.1.621' }
if(Q.1_3_agg==TRUE) {src.dat[src.dat$VARIANT %in% Q,"VARIANT"] <- "B.1.1.7"}
if(AY_agg==TRUE) {src.dat[src.dat$VARIANT %in% AY ,"VARIANT"] <- "B.1.617.2"}
if(B429_7_agg==TRUE) {src.dat[src.dat$VARIANT %in% B429 ,"VARIANT"] <- "B.1.427"}
# write.csv(src.dat, paste0("srcdat_",data_date,".csv"),row.names = FALSE)

svyNOGIS_ADJ = svydesign(ids=~SOURCE, strata=~STUSAB, weights= ~ SIMPLE_ADJ_WT, nest=TRUE, data=src.dat) # 2021-03-28: cluster corrected to SOURCE
svyREG = svydesign(ids=~STUSAB+SOURCE, strata=~HHS, weights= ~ SIMPLE_ADJ_WT, nest=TRUE, data=src.dat)

# Estimated degrees of freedom = number of clusters - number of strata
svyDF = with(src.dat, length(unique(paste(STUSAB, SOURCE))) - length(unique(STUSAB))) 


# Helper wrapper for svyciprop with CI: as vector (p, l, h), formatted string (pp.p (ll.l-hh.h)), or range (ll.l-hh.h)
# ~voc: character string of the variants of concern to analyze (if multiple listed, estimated proportion will be for the group of variants)
# ~geoid: character string/numeric of the geographic resolution to analyze (either 2 letter state code, 1:10 for HHS region, or "USA")
# ~svy: character string identifying what survey design is desired
# ~dates: character string of dates (format: "YYYY-MM-DD") to analyze
# ~ftn: logical argument indicating if dates represent fortnight end or week begin dates
# ~str: logical argument indicating if output should be formatted in single character object or as vector/dataframe
# ~range: logical argument indicating whether to present CI range vs CI bounds
# ~mut: logcial argument indicating whether to analyze mutation profiles or lineages
# ~level: numeric specifying confidence level (1-alpha)
myciprop = function(voc, geoid, svy, str=TRUE, range=FALSE, mut=FALSE,level = 0.95, ...) {
  #Create a binary indicator variable based on whether S_MUT or VARIANT is in voc
  if (mut) {VOC = grepl(paste0("^(?=.*\\",paste(voc,collapse="\\b)(?=.*\\"),"\\b)"), svy$variables$S_MUT,perl=T) }else {VOC = (svy$variables$VARIANT %in% voc)}
  
  if (ci.type=="xlogit"){
    res = svyciprop(~VOC, design = subset(update(svy, VOC=VOC), (STUSAB %in% geoid) | (HHS %in% geoid) | (geoid=="USA")),method="xlogit", ...)
  } else {
    res = svycipropkg(~VOC, design = subset(update(svy, VOC=VOC), (STUSAB %in% geoid) | (HHS %in% geoid) | (geoid=="USA")), ...)
  }
  res = c(res, confint(res))
  res = c(res, degf(subset(update(svy, VOC=VOC), (STUSAB %in% geoid) | (HHS %in% geoid) | (geoid=="USA"))))
  
  #add effective sample size to res output
  m <- eval(bquote(svymean(~as.numeric(VOC), subset(update(svy, VOC=VOC), (STUSAB %in% geoid) | (HHS %in% geoid) | (geoid=="USA"))))) #extracts the mean and SE for each estimate
  CVmean <- cv(svymean(~as.numeric(VOC), subset(update(svy, VOC=VOC), (STUSAB %in% geoid) | (HHS %in% geoid) | (geoid=="USA"))))
  deffmean <- deff(svymean(~as.numeric(VOC), subset(update(svy, VOC=VOC), (STUSAB %in% geoid) | (HHS %in% geoid) | (geoid=="USA")),deff=T))
  
  n.eff <- coef(m) * (1 - coef(m))/vcov(m)
  alpha <- 1-level #define alpha based on CI level
  n.eff.capped <- min( n.eff * (qt(alpha/2, nrow(subset(update(svy, VOC=VOC), (STUSAB %in% geoid) | (HHS %in% geoid) | (geoid=="USA"))) - 1)/qt(alpha/2, degf(subset(update(svy, VOC=VOC), (STUSAB %in% geoid) | (HHS %in% geoid) | (geoid=="USA")))))^2, nrow(subset(update(svy, VOC=VOC), (STUSAB %in% geoid) | (HHS %in% geoid) | (geoid=="USA"))))
  res = c(res,n.eff.capped,CVmean,deffmean)
  if (str) { 
    res = c(round( 100 * res[1:3], 1),res[4],res[5],res[6],res[7])
    if (range) {res = paste0(res[2], "-", res[3])} else {res = paste0(res[1], " (", res[2], "-", res[3], ") DF=",res[4],", Eff.sampsize=",res[5],", CVprop=",res[6],",DEff=",res[7])}
  }
  #modified to output degrees of freedom for each survey design subset
  
  res
} 

week_label = function(week_from_current) {
  strftime(Sys.Date() - as.numeric(strftime(Sys.Date(), format="%w")) + 7 * week_from_current, format="%m-%d")
}


#Functions for multinomial nowcast model
svymultinom = function(src.dat, mysvy, fmla=formula("as.numeric(as.factor(K_US)) ~ week + as.factor(HHS)")) {
  #src.dat=src.moddat; mysvy=mysvy; fmla=formula("as.numeric(as.factor(K_US))  ~ week")
  ## Get multinomial logistic regression coefficients (and Hessian unadapted to survey design)
  # Multinomial logistic regression (without survey design, but with survey weights) using nnet::multinom:
  # multinom_geoid = nnet::multinom(fmla, data=src.dat, weights=weights(mysvy), Hess=TRUE, maxit=1000, trace=FALSE)
  
  # aggregate data before fitting the multinomial model to vastly improve run time
  # see: https://git.biotech.cdc.gov/sars2seq/sc2_proportion_modeling/-/issues/6#note_86544
  fmla.vars = all.vars(fmla)
  mlm.dat = data.table::data.table(cbind(data.frame(src.dat)[, fmla.vars], weight=weights(mysvy)))[, .(weight=sum(weight)), by=fmla.vars]
  multinom_geoid = nnet::multinom(fmla, data=mlm.dat, weights=weight, Hess=TRUE, maxit=1000, trace=FALSE)
  
  # Format results to fit into svymle-like workflow
  num_var = length(unique(with(src.dat, eval(terms(fmla)[[2]])))) #gets the number of variants being modeled (should be length of model_vars plus 1 for others)
  fmla.no.response = formula(delete.response(terms(fmla))) #generates the model formula without the response term
   formulas = rep(list(fmla.no.response), num_var - 1) #repeats the model formula w/o response term for the number variants listed in model_vars
  # Response on first formula ensures inclusion as response later on
  formulas[[1]] = fmla #sets the first formula listed in formulas to the formula that includes a response term
  names(formulas) = paste0("b", 1:length(formulas) + 1) #sets the names of the formulas to correspond to beta coefficients for each variant
  rval = list(mlm = multinom_geoid, estimates = coefficients(multinom_geoid)) #creates a list that contains the mlm object and the coefficient estimates
  rval$estimates = as.list(data.frame(t(rval$estimates))) #transforms the estimates object to be a list where each element is a vector of the coefficients for a given variant (hhs model: Intercept,week, HHS 2:10; us model: Intercept, week)
  ## End multinomial regression
  
  
  # if the Hessian is not invertible, avoid errors by not running the rest of the function
  # (note: the se.multinom function will still run, but will assume the SE is 0)
  # if(det(multinom_geoid$Hessian) == -Inf) {
  if(det(multinom_geoid$Hessian) < -9e100) {
    # might want to replace -Inf with a very large negative number: e.g. -9e100
    
    # add empty items to the results (setting to NULL in list() does create the item)
    rval <- append(rval,
                   list(
                     'scores' = NULL,
                     'invinf' = NULL,
                     'sandwich' = NULL,
                     'SE' = NULL
                   ))
    
    # remove the Hessian
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
    # also print out an alert b/c warnings can sometimes get swamped by other minor warnings
    print(warning_message)
  } else {
    
    ## Define likelihood
    # Loglikelihood:
    lmultinom = function(v, ...) {
      # vectorized loglikelihood function
      # v (positive integer vector) is position of response variable in ordered list of possible values (1=reference)
      # ... vectors of  linear predictors
      b = cbind(0, ...)
      b = matrix(b, nrow=length(v))
      sapply(seq_along(v), function(rr) b[rr, v[rr]])  - log(rowSums(exp(b)))
    }
    # Gradient of loglikelihood:
    gmultinom = function(v, ...) {
      # vectorized partial derivatives of lmultinom(v, ...) with respect to linear predictors at ...
      # vectorized loglikelihood function
      # v (positive integer vector) is position of response variable in ordered list of possible values (1=reference)
      b = cbind(0, ...)
      b = matrix(b, nrow=length(v))
      p = exp(b)/rowSums(exp(b))
      delta_ij = p * 0
      delta_ij[1:length(v) + (v-1) * length(v)] = 1
      # Column n is partial derivative with respect to ...[n]:
      (delta_ij - p)[, -1]
    }
    ## End likelihood definitions
    
    ## Build dataframe to pass to svyrecvar
    # Adapted from code in svymle
    nms = c("", names(formulas)) #modifies formula names
    has.response = sapply(formulas, length) == 3 #logical that identifies which formula has a response variable
    ff = formulas[[which(has.response)]] #creates variable equal to regression formula that has response variable
    ff[[3]] = 1 #sets predictor terms to 1 so: as.numeric(as.factor(K_US)) ~ 1
    y = eval.parent(model.frame(ff, data = src.dat, na.action = na.pass)) #I think this just gets you what the rankings are for the variants from src.dat (i.e. the response data for the regression model)
    formulas[[which(has.response)]] = formula(delete.response(terms(formulas[[which(has.response)]]))) #sets the first formula to a formula without the response term
    mf = vector("list", length(formulas)) #creates empty list the same length of the formulas object
    vnms = unique(do.call(c, lapply(formulas, all.vars))) #vector with the names of the regression predictors
    uformula = make.formula(vnms) #converts the vnms to a formula object w/o the response term
    mf = model.frame(uformula, data = src.dat, na.action = na.pass)#gets the values for the model predictors from src.dat
    mf = cbind(`(Response)` = y, mf) #adds a column for the variant rankings to mf with the column header as as.numeric(as.factor(K_US))
    mf = mf[, !duplicated(colnames(mf)), drop = FALSE]#I think this just ensures that there aren't any duplicated columns
    weights = weights(mysvy) #gets the svy weights from the survey design object
    wtotal = sum(weights) #gets the sum of the weights
    Y = mf[, 1] #defines Y as the variant
    mmFrame = lapply(formulas, model.frame, data = mf) #for each formula listed in formula, generates the underlying data for those parameters
    mm = mapply(model.matrix, object = formulas, data = mmFrame, SIMPLIFY = FALSE) #creates model design objects for each variant
    np = c(0, cumsum(sapply(mm, NCOL))) #gets the cumulative number of columns across all dataframes listed in mm
    parnms = lapply(mm, colnames) #gets the column names from all the dataframes in mm
    for (i in 1:length(parnms)) parnms[[i]] = paste(nms[i + 1], parnms[[i]], sep = ".") #formats column names in each list element to correspond to the coefficient
    parnms = unlist(parnms) #converts the parameter names from a list to a vector
    
    theta = unlist(rval$estimates) # gets the estimated coefficients from multinom as a vector
    args = vector("list", length(nms)) #creates empty list the same length as the number of variants plus other
    args[[1]] = Y #sets first list element to the variant rank data (i.e. response variable)
    names(args) = nms #assigns the names of the args list to nms
    for (i in 2:length(nms)) args[[i]] = drop(mm[[i - 1]] %*% theta[(np[i - 1] + 1):np[i]]) #gets the coefficient estimates for each of the regression coefficients for each modeled variant, and places in the appropriate model design object
    
    deta = matrix(do.call("gmultinom", args), ncol=length(args)-1) #this takes the args list and uses the gmultinom function above to estimate the gradient of the log likelihood
    rval$scores = NULL # creates an empty list element called scores in the rval list
    reorder = na.omit(match(names(formulas), nms[-1])) # creates numerical vector the same length as the names of formulas
    for (i in reorder) rval$scores = cbind(rval$scores, deta[, i] * weights * mm[[i]]) #for each variant this multiplies the gradient the loglikelihood by the survey weights
    rval$invinf = solve(-multinom_geoid$Hessian) #not 100% sure what this step does it looks like it's solving the inverse (negative?) multinomial regression object and adds as a list element to rval
    dimnames(rval$invinf) = list(parnms, parnms) #assigns names to the rval$invinf
    db = rval$scores %*% rval$invinf # this uses matrix multiplication to multiply the "scores" (i.e., the gradient loglikelihood * survey weights) by the inverse(?) multinomial regression object
    
    #Everything up until now has been formatting steps to get the estimates/data formatted in way that's compatible with the svyrecvar function
    rval$sandwich = svyrecvar(db, mysvy$cluster, mysvy$strata, mysvy$fpc, postStrata = mysvy$postStrata) #Computes the variance of a total under multistage sampling, using a recursive descent algorithm.The object that's returned ("sandwich") is the covariance matrix
    rval$SE = sqrt(diag(rval$sandwich)) #estimating the "SE" (although isn't this sd?) as the square root of the diagonal elements of the variance-covariance matrix (i.e. the variance)
    rval$estimates = unlist(rval$estimates) #converts the "estimates" object from a list to a vector
    names(rval$estimates) = names(rval$SE) #assigns the names of "estimates" in rval to be the same as the names of "SE"
    rval$mlm$svyvcov = rval$sandwich #adds the variance-covariance matrix to the mlm object (mlm being the results object you get from solving the multinomial model)
  }
  
  rval #makes sure the rval object is the output from this function
}


se.multinom = function(mlm, week, geoid="USA", composite_variant=NULL) { 
  #mlm=svymlm_us$mlm; week=85; geoid="USA"; composite_variant=agg_var_mat
  # Modified to handle composite variants (eg. Delta = combination of multiple lineages in multinomial model) [2021-09-29]
  #   composite_variant (if not NA) is a matrix with one column per model lineage and one row per composite variant
  #     matrix element of 1 marks each component lineage (column) for each composite variant (row)
  #   example: matrix(c(1, 0, 0, 1, 0, 0), nrow=1) designates a variant comprising the first and fourth lineages in the model 
  #     estimates for composites will be appended at end of p_i and se.p_i     
  # Modified to handle geoid = state [2021-08-09]
  #   example: se.multinom(temp_state, 80, "GA")
  # Modified to handle output from svymultinom [2021-08-15]
  #   example: svymlm = function(src.dat, mysvy, fmla); se.multinom(svymlm$mlm, 80, "GA")
  # mlm = model output, with Hessian; geoid can be HHS Region (1:10)
  # If mlm is without Hessian, all se's are set to zero
  # mlm not geographically stratified for geoid="USA": ~ (week - current_week)
  # mlm, for HHS Regions: ~ week + HHS (1 as reference level)
  cf = coefficients(mlm)
  # Updated to incorporate survey design based variances (2021-08-15):
  if ("svyvcov" %in% names(mlm)) {vc = mlm$svyvcov} else {
    if ("Hessian" %in% names(mlm)) {vc = solve(mlm$Hessian)} else {vc=rep(0, length(cf)^2)}
  }
  dim(vc) = rep(rev(dim(cf)), 2) # Rearrange terms to ease pulling out relevant submatrix
  ref_geoid = ((length(mlm$xlevels) == 0) || (geoid == mlm$xlevels[[1]][1])) # Boolean: geoid is reference level
  if (ref_geoid) {n_geoid = NULL} else {n_geoid = which(mlm$xlevels[[1]] == geoid)}
  if (ref_geoid) {hhs1tail = NULL} else {hhs1tail = n_geoid + 1}
  indices = c(1, 2, hhs1tail)
  if (ref_geoid) {hhs1tail = NULL} else {hhs1tail = 1}
  coeffs = c(1, week, hhs1tail)
  # b is vector of coefficients of time (week)
  sub.vc = vc[indices,,indices,] 
  y_i = c(0, cf[, indices] %*% coeffs)
  vc.y_i = outer(1:dim(sub.vc)[2], 1:dim(sub.vc)[4], Vectorize(function(i, j) c(coeffs %*% sub.vc[, i, , j] %*% coeffs)))
  vc.y_i = rbind(0, cbind(0, vc.y_i))
  p_i = exp(y_i)/sum(exp(y_i))
  # Taylor series based variance:
  # dp_i/dy_j = delta_ij p_i - p_i * p_j
  dp_dy = diag(p_i) - outer(p_i, p_i, `*`)
  # Steps separated out to enable proportion estimation for composite variables [2021-09-29]
  #   former single step: se.p_i = as.vector(sqrt(diag(dp_dy %*% vc.y_i %*% dp_dy)))
  p.vcov = dp_dy %*% vc.y_i %*% dp_dy
  se.p_i = as.vector(sqrt(diag(p.vcov)))
  if (!is.null(composite_variant)) {
    composite_variant = list(matrix = composite_variant, 
                             p_i = as.vector(composite_variant %*% p_i),
                             se.p_i = as.vector(sqrt(diag(composite_variant %*% p.vcov %*% t(composite_variant))))
    )
  }
  list(p_i=p_i, se.p_i=se.p_i, b_i=c(0, cf[, 2]), se.b_i= c(0, sqrt(diag(vc[2,,2,]))), composite_variant=composite_variant)
}  

#Function to get binomial confidence interval based on the point estimated proportion and associated SE
svyCI = function(p, s) {
  if(s == 0){
    out = c(NA,NA)
  } else {
    n=p*(1-p)/s^2
    out=prop.test(n * p, n)$conf.int
  }
  c(out[1],out[2])
}

#######################################################
# 
# 
# ```
# 
# # Background
# 
# This brief summarizes findings from the variant surveillance system for variants of concern and of interest. 
# 
# The sequence are those submitted to NS3 and by contract with lab vendors. Both raw numbers and weighted estimates are displayed. Weights allow estimates to be representative of all infected persons, aggregated by week of specimen collection and by state. Currently, the probability of selection for each sample is assumed to be independent of data source, but that may be modified as more data on sampling protocol become available. An adjustment is applied to account for oversampling of SGTF specimens in certain data streams.
# 
# # Results ------------------------------------------------------------------------------------------
# 
# ## Current top variants by share and variants of concern or interest
# 
# Based on sequences where the specimen was collected in the current or previous `r n_recent_weeks` weeks, the count of sequences, and variant share (weighted percent) are:
# 
# Select display and model variants- this basically selects which variables are the top variables to include in the model
us_var = sort(
  prop.table(xtabs(SIMPLE_ADJ_WT ~ VARIANT, subset(src.dat, week >= current_week - n_recent_weeks & VARIANT != "None"))),
  decreasing=TRUE)
us_seq = table(subset(src.dat, week >= current_week - n_recent_weeks)$VARIANT)
us_rank = names(us_var)
hhs_var = prop.table(xtabs(SIMPLE_ADJ_WT ~ HHS + VARIANT, subset(src.dat, week >= current_week - n_recent_weeks)), 1)
hhs_rank = apply(hhs_var, 1, function(rr) names(sort(rr, decreasing=TRUE)))
all_tops = us_rank[us_rank %in% c(us_rank[1:n_top], voc)] # Ordered by national rank; for display
model_vars = us_var[all_tops] # For multinomial model
# model_vars = us_var[us_var >= share_cutoff | names(us_var) %in% c(all_tops, hhs_rank[1:n_top,])] # For multinomial model
#added variables to indicate the time window specified by Clint

if (display_option=="top7") {display_vars = head(names(model_vars), 7)} else {display_vars = voc}
display_indices = which(names(model_vars) %in% display_vars)
display_vars = names(model_vars)[display_indices] # Ensuring correct irder in display


## Model-based smoothed trends in variant share: National and by HHS Region ------------------------

#add variant ranks and current week to src.dat
src.dat$K_US = sapply(1:nrow(src.dat), function(nn) which(c(names(model_vars), src.dat$VARIANT[nn]) == src.dat$VARIANT[nn])[1]) # makes sure each seq is assigned a number based on the weighted proportion in the last few weeks (1 =most common)
src.dat$current_week = current_week

#create a subset of src.dat that only contains the weeks that will be included in multinomial model
src.moddat = subset(src.dat, week >= max(week) - model_weeks) # Modeling window

#create survey design object based on subsetted src.dat
mysvy = svydesign(ids=~SOURCE, strata=~STUSAB+yr_wk, weights= ~ SIMPLE_ADJ_WT, nest=TRUE, data=src.moddat)

# RUN2
if (length(grep("Run2",tag))!=0) { #Only run the nowcast when you're looking at Run2 i.e. expanded list of lineages
  #runs the svymultinom function and saves the output
  svymlm_hhs = svymultinom(src.moddat, mysvy,fmla=formula("as.numeric(as.factor(K_US)) ~ week + as.factor(HHS)"))
  pred_hhs.df = expand.grid(week = seq(-display_lookback, 2, 0.02) + current_week, HHS=sort(unique(src.moddat$HHS)))
  pred_hhs.df = cbind(pred_hhs.df, predict(svymlm_hhs$mlm, pred_hhs.df, type="probs"))
  bp_hhs = xtabs(NORM_WTS ~ HHS + week + VARIANT, subset(src.moddat, week < current_week & week >= current_week - display_lookback))
  bp_hhs = prop.table(bp_hhs, 1:2)[,, display_vars]
  dimnames(bp_hhs)[[2]] = week_label(as.numeric(dimnames(bp_hhs)[[2]]) - current_week)
  
  svymlm_us = svymultinom(src.moddat, mysvy,fmla=formula("as.numeric(as.factor(K_US))  ~ week"))
  pred_us.df = expand.grid(week = seq(-display_lookback, 2, 0.02) + current_week) #creates smoothed vector of week values to predict over (for the model smoothed figure)
  pred_us.df = cbind(pred_us.df, predict(svymlm_us$mlm, pred_us.df, type="probs")) #generates predicted probabilities over timeframe defined above
  bp_us = xtabs(NORM_WTS ~ week + VARIANT, subset(src.moddat, week < current_week & week >= current_week - display_lookback))
  bp_us = prop.table(bp_us, 1)
  bp_us = bp_us[, display_vars]
  rownames(bp_us) = week_label(as.numeric(rownames(bp_us)) - current_week)
  
  
  # Weighted variant shares of the top variants in the past `r display_lookback` weeks (number of sequences collected weekly above each bar), and model-based smoothed estimates, nationwide:
  # 
  # 
  # ```{r graphics us, warning=FALSE, message=FALSE, fig.width = 6, fig.height=6, dpi=300}
  
  #knitr::kable(round(100*bp_us, 1), row.names=TRUE)
  
  # Set up colors
  col.dk = hcl.colors(length(display_vars), palette="TealRose", alpha=0.8)
  names(col.dk) = display_vars
  
  stub = paste0(script.basename, "/results/wtd_shares_", format(Sys.Date(), "%Y%m%d"), "_")
  
  
  if (fig_gen_run) jpeg(paste0(stub, "barplot_US", tag, custom_tag, ".jpg"), width=1500, height=1500, pointsize=40)
  bp = barplot(100*t(bp_us), xlab="Week beginning", ylab="Weighted variant share (%)", main="Nationwide",
               border=NA, ylim=110*0:1,  col=col.dk, 
               names.arg=rownames(bp_us), 
               legend.text=display_vars, args.legend=list(x="topleft", bty="n", border=NA))
  text(bp, 3 + colSums(100*t(tail(bp_us, 12))), with(subset(src.dat, week < current_week & week >= current_week - display_lookback), table(week)), cex=0.7)
  if (fig_gen_run) dev.off()
  if (fig_gen_run) jpeg(paste0(stub, "projection_US", tag, custom_tag, ".jpg"), width=1500, height=1500, pointsize=40)
  bp = barplot(100*t(pred_us.df[, 1 + display_indices]), xlab="Week beginning", ylab="Weighted variant share (%)", main="Nationwide",
               space=0, border=NA, ylim=110*0:1,  col=col.dk,
               names.arg=ifelse(pred_us.df$week %% 1 == 0, week_label(pred_us.df$week - current_week), NA),
               legend.text=display_vars, args.legend=list(x="topleft", bty="n", border=NA))
  pc = unlist(100 * subset(pred_us.df, week==current_week)[, 1+display_indices])
  y = cumsum(pc) - pc/2
  text(1.02 * tail(bp, 1), y, round(pc, 1), cex=0.7, xpd=TRUE, adj=c(0, 0.5))
  x = bp[which(pred_us.df$week %in% (current_week + c(-2, 0, 2)))] 
  rect(x[1], 0, x[2], 100, border=NA, col="#00000020")
  rect(x[2], 0, x[3], 100, border=NA, col="#00000040")
  if (fig_gen_run) dev.off()
  
  
  # 
  # The following is a "nowcast"; the vertical axis depicts the variant share growth rate (derivative of log of variant proportion with respect to time). The axis on the right shows an estimate of variant transmissibility with respect to the overall mean transmissibility.
  # 
  
  
  us.summary = se.multinom(svymlm_us$mlm, current_week,geoid="USA")
  se.gr = with(us.summary, 100 * exp(sqrt(se.b_i^2 + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2))) - 100) 
  
  # WIP: getting better CIs (2021-04-26):
  # svyDF = with(svyNOGIS_ADJ$variables, length(unique(paste(STUSAB, SOURCE))) - length(unique(STUSAB)))
  # us.summary$se.p_i = with(us.summary, sqrt(se.p_i^2 + p_i * (1 - p_i)/svyDF))
  gr = with(us.summary,  100 * exp(b_i - sum(p_i * b_i)) - 100)
  if (fig_gen_run) png(paste0(stub, "growthrate_US",tag,".png"), width=8, height=8,units="in", pointsize=16,res=1000)
  plot(100 * us.summary$p_i, gr, log="x", type="n", ylim=range(gr + 1.96 * se.gr, gr - 1.96 * se.gr),
       xaxt="n", xlim=c(0.0005,110),
       xlab = "Nowcast Estimated Proportion (%)", ylab="Week over week growth rate (%)", main="Nationwide")
  if(is.null(svymlm_us$SE)){
    mtext(text = "*No SE estimates b/c of non-invertible Hessian in multinomial model fit.", 
          side = 3,
          line = 0, 
          cex = 0.75, 
          font = 4, 
          col = 'red' )
  }
  axis(1, at=c(0.001,0.01, 0.1, 1, 10,100), labels=c(0.001,0.01, 0.1, 1, 10,100))
  abline(h=0, col="grey65")
  for (vv in seq(model_vars)) {
    lines(100 * rep(us.summary$p_i[vv], 2), gr[vv] + 1.96 * c(1,-1) * se.gr[vv], col="blue",lwd=2)
    lines(pmax(#1e-2, 
      100 * us.summary$p_i[vv] + 196 * c(1,-1) * us.summary$se.p_i[vv]), rep(gr[vv], 2), col="blue",lwd=2)
  }
  text(100 * us.summary$p_i, gr, c(names(model_vars),""), cex=0.85, col="grey25", adj=1.15)
  
  if (fig_gen_run)  dev.off()
  
  gr_tab = cbind(variant=c(names(model_vars), "OTHER"), 
                 variant_share=(100 * us.summary$p_i), 
                 growth_rate=gr)
  
  write.csv(gr_tab, paste0(script.basename, "/results/wow_growth_variant_share", Sys.Date(), tag, custom_tag, ".csv"), row.names=FALSE)
  
  # 
  # Weighted variant shares of the top variants in the past `r display_lookback` weeks (number of sequences collected weekly above each bar), and model-based smoothed estimates, for each HHS region:
  # 
  
  for (hhs in sort(unique(src.dat$HHS))) {
    if (fig_gen_run) jpeg(paste0(stub, "barplot_HHS", hhs,tag, custom_tag, ".jpg"), width=1500, height=1500, pointsize=40)
    bp = barplot(100*t(bp_hhs[hhs,,]), xlab="Week beginning", ylab="Weighted variant share (%)", main=paste("HHS Region", hhs),
                 border=NA, ylim=110*0:1,  col=col.dk, 
                 names.arg=rownames(bp_hhs[hhs,,] - current_week), 
                 legend.text=display_vars, args.legend=list(x="topleft", bty="n", border=NA))
    text(bp, 3 + colSums(100*t(tail(bp_hhs[hhs,,], 12))), with(subset(src.dat, HHS==hhs & week < current_week & week >= current_week - display_lookback), table(week)), cex=0.7)
    if (fig_gen_run) dev.off()
    if (fig_gen_run) jpeg(paste0(stub, "projection_HHS", hhs,tag, custom_tag, ".jpg"), width=1500, height=1500, pointsize=40)
    pred.df = subset(pred_hhs.df, HHS==hhs)[, -2]
    bp = barplot(100*t(pred.df[, 1 + display_indices]), xlab="Week beginning", ylab="Weighted variant share (%)", main=paste("HHS Region", hhs),
                 space=0, border=NA, ylim=110*0:1, col=col.dk,
                 names.arg=ifelse(pred.df$week %% 1 == 0, week_label(pred.df$week - current_week), NA) ,
                 legend.text=display_vars, args.legend=list(x="topleft", bty="n", border=NA))
    pc = unlist(100 * subset(pred.df, week==current_week)[, 1+display_indices])
    y = cumsum(pc) - pc/2
    text(1.02 * tail(bp, 1), y, round(pc, 1), cex=0.7, xpd=TRUE, adj=c(0, 0.5))
    x = bp[which(pred.df$week %in% (current_week + c(-2, 0, 2)))] 
    rect(x[1], 0, x[2], 100, border=NA, col="#00000020")
    rect(x[2], 0, x[3], 100, border=NA, col="#00000040")
    if (fig_gen_run) dev.off()
    
    if (fig_gen_run) jpeg(paste0(stub, "growthrate_HHS", hhs, tag, custom_tag, ".jpg"), width=1500, height=1500, pointsize=40)
    hhs.summary = se.multinom(svymlm_hhs$mlm, current_week, hhs)
    se.gr = with(hhs.summary, 100 * exp(sqrt(se.b_i^2 + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2))) - 100) 
    gr = with(hhs.summary, 100 * exp(b_i - sum(p_i * b_i)) - 100)
    plot(100 * hhs.summary$p_i, gr, log="x", type="n", ylim=c(min(range(gr + 1.96 * se.gr, gr - 1.96 * se.gr)[1], -75),
                                                              max(range(gr + 1.96 * se.gr, gr - 1.96 * se.gr)[2],100)),
         xaxt="n", xlim=c(0.0005,110),
         xlab = "Weighted share (%)", ylab="Week over week growth rate (%)", main=paste("HHS Region", hhs))
    if(is.null(svymlm_hhs$SE)){
      mtext(text = "*No SE estimates b/c of non-invertible Hessian in multinomial model fit.", 
            side = 3,
            line = 0, 
            cex = 0.75, 
            font = 4, 
            col = 'red' )
    }
    axis(1, at=c(0.001,0.01, 0.1, 1, 10,100), labels=c(0.001,0.01, 0.1, 1, 10,100))
    abline(h=0, col="grey65")
    for (vv in seq(model_vars)) {
      lines(100 * rep(hhs.summary$p_i[vv], 2), gr[vv] + 1.96 * c(1,-1) * se.gr[vv], col="blue",lwd=2)
      lines(pmax(#1e-2, 
        100 * hhs.summary$p_i[vv] + 1.96 * c(1,-1) * hhs.summary$se.p_i[vv]), rep(gr[vv], 2), col="blue",lwd=2)
      }
    #identify which variants have share>1% or growth rate>0%
    labels=unique(c(which(100 * hhs.summary$p_i>1), which(gr>1)))
    text(100 * hhs.summary$p_i, gr, c(names(model_vars),""), cex=0.85, col="grey25", adj=1.15)
    
    if (fig_gen_run) dev.off()
    }
}



# 
# ## Technical notes
# 
# * Weights are imputed from regional incidence if testing data for a state are unavailable
# 
# * Model-based smoothing is performed using a multinomial logistic model
# 
# 
# ```{r generate spreadsheet of estimates, warning=FALSE, message=FALSE, eval=FALSE}

reported_variants = voc # voc for short list, all_tops for long; (will differ in what falls under "other")
#reported_variants = moi
if(svy.type=="svyREG"){
svyDES = svydesign(ids=~STUSAB+SOURCE, strata=~HHS, weights= ~ SIMPLE_ADJ_WT, nest=TRUE, data=src.dat)
} else {
svyDES = svydesign(ids=~SOURCE, strata=~STUSAB+yr_wk, weights= ~ SIMPLE_ADJ_WT, nest=TRUE, data=src.dat)
}


src.dat$VARIANT2=as.character(src.dat$VARIANT)
src.dat[src.dat$VARIANT %notin% voc,"VARIANT2"] <- "Other"

if (length(grep("Run3",tag))==0){ #Run the fortnight and weekly estimates when 
  
#if (update_wk == 2) {ftnts = head(sort(unique(subset(src.dat, DAY >= as.numeric(as.Date("2021-01-25") - week0day1))$FORTNIGHT_END)), -1)} else{
#  ftnts = head(sort(unique(subset(src.dat, DAY >= as.numeric(as.Date("2021-01-25") - week0day1))$FORTNIGHT_END)))}
ftnts=(unique(subset(src.dat, as.Date(FORTNIGHT_END)>=as.Date("2021-05-08") & as.Date(FORTNIGHT_END)<=time_end)$FORTNIGHT_END))

all.ftnt = expand.grid(Variant=reported_variants, Fortnight_ending=ftnts, USA_or_HHSRegion=c("USA", 1:10))[, 3:1]
ests = apply(all.ftnt, 1, function(rr) myciprop(rr[3], rr[1], subset(svyDES, FORTNIGHT_END == rr[2]), FALSE))#, method="xlogit")) took this out bc specified method = beta above for KG CI
all.ftnt = cbind(all.ftnt, Share=ests[1,], Share_lo=ests[2,], Share_hi=ests[3,],DF=ests[4,],eff.size=ests[5,], cv.mean=ests[6,], deff=ests[7,])

others = expand.grid(Variant="Other", Fortnight_ending=ftnts, USA_or_HHSRegion=c("USA", 1:10))[, 3:1]
ests.others = apply(others, 1, function(rr) myciprop(reported_variants, rr[1], subset(svyDES, FORTNIGHT_END == rr[2]), FALSE))
others = cbind(others, Share=1-ests.others[1,], Share_lo=1-ests.others[3,], Share_hi=1-ests.others[2,],DF=ests.others[4,],eff.size=ests.others[5,],cv.mean=ests.others[6,], deff=ests.others[7,])

all.ftnt = rbind(all.ftnt, others)


#Getting counts of sequences by lineage, location and date


dat2 <- src.dat[src.dat$FORTNIGHT_END>=as.Date("2021-05-08","%Y-%m-%d"),]
dat2 <- subset(dat2, as.Date(dat2$FORTNIGHT_END) <= as.Date(tail(unique(all.ftnt$Fortnight_ending),1)))
raw_counts_REG <- aggregate(count~VARIANT2+FORTNIGHT_END+HHS, data=dat2, FUN=sum, drop=FALSE)
raw_counts_REG$HHS <- as.character(raw_counts_REG$HHS)# <- aggregate(count~VARIANT2+FORTNIGHT_END+HHS, data=dat2, FUN=sum)

raw_counts_US <- aggregate(count~VARIANT2+FORTNIGHT_END, data=dat2, FUN=sum, drop=FALSE)
raw_counts_US <- cbind(raw_counts_US[,1:2],HHS="USA",count=raw_counts_US[,3])
raw_counts <- rbind.data.frame(raw_counts_US,raw_counts_REG)

#merge sequence counts with weighted proportions estimates
all.ftnt2 <- merge(all.ftnt, raw_counts,by.x=c("USA_or_HHSRegion","Fortnight_ending","Variant"),by.y=c("HHS","FORTNIGHT_END","VARIANT2"),all=T) 

all.ftnt2[is.na(all.ftnt2$count)==T,"count"] <- 0

#calculate denominator counts
dss <- aggregate(count~USA_or_HHSRegion+Fortnight_ending, data=all.ftnt2,FUN=sum)
names(dss)[grep("count",names(dss))] <- "denom_count"
all.ftnt2 <- merge(all.ftnt2, dss)

#set the Share 0 and CI limits to NA when the count for a lineage is 0
all.ftnt2$Share=ifelse(all.ftnt2$Share!=0 & all.ftnt2$count==0,0,all.ftnt2$Share)
all.ftnt2$Share_lo=ifelse(is.na(all.ftnt2$Share_lo)==F & all.ftnt2$count==0,NA,all.ftnt2$Share_lo)
all.ftnt2$Share_hi=ifelse(is.na(all.ftnt2$Share_hi)==F & all.ftnt2$count==0,NA,all.ftnt2$Share_hi)

#calculate absolute CI width
all.ftnt2$CI_width=all.ftnt2$Share_hi-all.ftnt2$Share_lo

#generate NCHS flags
all.ftnt2$flag_df=as.numeric(all.ftnt2$DF<8)
all.ftnt2$flag_eff.size=ifelse(all.ftnt2$eff.size<30|is.na(all.ftnt2$eff.size)==T,1,0)
all.ftnt2$flag_dss=ifelse(all.ftnt2$denom_count<30|is.na(all.ftnt2$denom_count)==T,1,0)
all.ftnt2$flag_abs.ciw=ifelse(all.ftnt2$CI_width>0.30|is.na(all.ftnt2$CI_width)==T,1,0)
all.ftnt2$flag_rel.ciw=ifelse(((all.ftnt2$CI_width/all.ftnt2$Share)*100)>130|is.na((all.ftnt2$CI_width/all.ftnt2$Share)*100)==T,1,0)

all.ftnt2$nchs_flag=ifelse(all.ftnt2$flag_df==1|all.ftnt2$flag_eff.size==1|all.ftnt2$denom_count==1|all.ftnt2$flag_abs.ciw==1|all.ftnt2$flag_rel.ciw==1,
                           1,0)
all.ftnt2$nchs_flag_wodf=ifelse(all.ftnt2$flag_eff.size==1|all.ftnt2$denom_count==1|all.ftnt2$flag_abs.ciw==1|all.ftnt2$flag_rel.ciw==1,
                           1,0)

all.ftnt2=all.ftnt2[,c("USA_or_HHSRegion","Fortnight_ending","Variant","Share",           
                       "Share_lo","Share_hi","count","denom_count","DF","eff.size",
                       "CI_width","nchs_flag","nchs_flag_wodf")]

all.ftnt2 <- all.ftnt2[order(all.ftnt2$USA_or_HHSRegion),]

write.csv(all.ftnt2, paste0(script.basename, "/results/variant_share_weighted_", ci.type, "CI_", svy.type, "_", data_date, tag, custom_tag, ".csv"), row.names=FALSE)

# Weekly estimates --------------------------------------------------------------------
#if (update_wk == 2) {wks = head(sort(unique(subset(src.dat, DAY >= as.numeric(as.Date("2021-01-25") - week0day1))$yr_wk)), -1)} else { 
#  wks = sort(unique(subset(src.dat, DAY >= as.numeric(as.Date("2021-01-25") - week0day1))$yr_wk))-2}

wks = sort(unique(subset(src.dat, as.Date(yr_wk) >= as.Date("2021-05-02") & as.Date(yr_wk) <= time_end-6)$yr_wk)) # all in 2021

all.wkly = expand.grid(Variant=reported_variants, Week_of=wks, USA_or_HHSRegion=c("USA", 1:10))[, 3:1]
ests = apply(all.wkly, 1, function(rr) myciprop(rr[3], rr[1], subset(svyDES, yr_wk == rr[2]), FALSE))
all.wkly = cbind(all.wkly, Share=ests[1,], Share_lo=ests[2,], Share_hi=ests[3,],DF=ests[4,],eff.size=ests[5,], cv.mean=ests[6,])

others = expand.grid(Variant="Other", Week_of=wks, USA_or_HHSRegion=c("USA", 1:10))[, 3:1]
ests.others = apply(others, 1, function(rr) myciprop(reported_variants, rr[1], subset(svyDES, yr_wk == rr[2]), FALSE))
others = cbind(others, Share=1-ests.others[1,], Share_lo=1-ests.others[3,], Share_hi=1-ests.others[2,],DF=ests.others[4,],eff.size=ests.others[5,], cv.mean=ests.others[6,])

all.wkly = rbind(all.wkly, others)
all.wkly$WEEK_END = as.Date(all.wkly$Week_of) + 6

#generate sequence counts by lineage, location and date
dat2 <- subset(src.dat, as.Date(yr_wk) >= as.Date("2021-05-02"))
dat2$WEEK_END = as.Date(dat2$yr_wk) + 6

#make sure dates match dates in proportions dataframe 
dat2 <- subset(dat2, WEEK_END <= tail(unique(all.wkly$WEEK_END),1))

raw_counts_REG <- aggregate(count~VARIANT2+WEEK_END+HHS, data=dat2, FUN=sum,drop=FALSE)
raw_counts_REG$HHS <- as.character(raw_counts_REG$HHS)

raw_counts_US <- aggregate(count~VARIANT2+WEEK_END, data=dat2, FUN=sum,drop=FALSE)
raw_counts_US <- cbind(raw_counts_US[,1:2],HHS="USA",count=raw_counts_US$count)
raw_counts <- rbind.data.frame(raw_counts_US,raw_counts_REG)

#merge sequence counts with weighted proportions estimates
all.wkly2 <- merge(all.wkly, raw_counts,by.x=c("USA_or_HHSRegion","WEEK_END","Variant"),by.y=c("HHS","WEEK_END","VARIANT2"),all=T) 

all.wkly2[is.na(all.wkly2$count)==T ,"count"] <- 0


#calculate denominator counts
dss <- aggregate(count~USA_or_HHSRegion+WEEK_END, data=all.wkly2,FUN=sum)
names(dss)[grep("count",names(dss))] <- "denom_count"
all.wkly2 <- merge(all.wkly2, dss)

#calculate absolute CI width
all.wkly2$CI_width=all.wkly2$Share_hi-all.wkly2$Share_lo

#set the Share 0 and CI limits to NA when the count for a lineage is 0
all.wkly2$Share=ifelse(all.wkly2$Share!=0 & all.wkly2$count==0,0,all.wkly2$Share)
all.wkly2$Share_lo=ifelse(is.na(all.wkly2$Share_lo)==F & all.wkly2$count==0,NA,all.wkly2$Share_lo)
all.wkly2$Share_hi=ifelse(is.na(all.wkly2$Share_hi)==F & all.wkly2$count==0,NA,all.wkly2$Share_hi)

#generate NCHS flags
all.wkly2$flag_df=ifelse(all.wkly2$DF<8,1,0)
all.wkly2$flag_eff.size=ifelse(all.wkly2$eff.size<30 | is.na(all.wkly2$eff.size)==T,1,0 )
all.wkly2$flag_dss=ifelse(all.wkly2$denom_count<30| is.na(all.wkly2$denom_count)==T,1,0 )
all.wkly2$flag_abs.ciw=ifelse(all.wkly2$CI_width>0.30 | is.na(all.wkly2$CI_width)==T,1,0)
all.wkly2$flag_rel.ciw=ifelse(((all.wkly2$CI_width/all.wkly2$Share)*100)>130 | is.na((all.wkly2$CI_width/all.wkly2$Share)*100)==T,1,0)

all.wkly2$nchs_flag=ifelse(all.wkly2$flag_df==1|all.wkly2$flag_eff.size==1|all.wkly2$denom_count==1|all.wkly2$flag_abs.ciw==1|all.wkly2$flag_rel.ciw==1,
                           1,0)
all.wkly2$nchs_flag_wodf=ifelse(all.wkly2$flag_eff.size==1|all.wkly2$denom_count==1|all.wkly2$flag_abs.ciw==1|all.wkly2$flag_rel.ciw==1,
                           1,0)

all.wkly2$count_LT20=ifelse(all.wkly2$count<20,1,0)
all.wkly2$count_LT10=ifelse(all.wkly2$count<10,1,0)

all.wkly2=all.wkly2[,c("USA_or_HHSRegion","WEEK_END","Variant","Share",           
                       "Share_lo","Share_hi","count","denom_count","DF","eff.size",
                       "CI_width","nchs_flag","nchs_flag_wodf","count_LT20","count_LT10")]

all.wkly2 <- all.wkly2[order(all.wkly2$USA_or_HHSRegion),]

write.csv(all.wkly2, paste0(script.basename, "/results/variant_share_weekly_weighted_", ci.type, "CI_", svy.type, "_", data_date, tag, custom_tag, ".csv"), row.names=FALSE)

}



# State-level estimates - Rolling 4 wk bins --------------------------------------------------------------------

# RUN3
if (length(grep("Run3",tag))!=0){ #Only Run the state-level estimates when running the state list of lineages (i.e. Run 3) 
  
  #get the week number that corresponds to date defined in state_time_end
  data_week=c()
  for(i in 1:length(state_time_end)){
    data_week[i] <-unique(src.dat$week[src.dat$yr_wk==state_time_end[i]-6]) #have to subtract 6 days because the yr_wk variable defines week starting on Sunday whereas the state_time_end is defined as week ending Saturday
  }
  
  svyDES = svydesign(ids=~SOURCE, strata=~STUSAB+yr_wk, weights= ~ SIMPLE_ADJ_WT, nest=TRUE, data=src.dat) #redefine survey design to the new design for state-level estimates
  
  #
  all.state = expand.grid(Variant=reported_variants, Roll_4wk_end=data_week,State=sort(unique(src.dat$STUSAB)))[, 3:1]
  ests = apply(all.state, 1, function(rr) myciprop(rr[3], rr[1], subset(svyDES, week >= (as.numeric(rr[2])- 3) & week < (as.numeric(rr[2])+1)), FALSE))
  all.state = cbind(all.state, Share=ests[1,], Share_lo=ests[2,], Share_hi=ests[3,],DF=ests[4,],eff.size=ests[5,], cv.mean=ests[6,])
  
  others = expand.grid(Variant="Other", Roll_4wk_end=data_week,State=sort(unique(src.dat$STUSAB)))[, 3:1]
  ests.others = apply(others, 1, function(rr) myciprop(reported_variants, rr[1],subset(svyDES, week >= (as.numeric(rr[2])- 3) & week < (as.numeric(rr[2])+1)), FALSE))
  others = cbind(others, Share=1-ests.others[1,], Share_lo=1-ests.others[3,], Share_hi=1-ests.others[2,],DF=ests.others[4,],eff.size=ests.others[5,], cv.mean=ests.others[6,])
  
  all.state = rbind(all.state, others)
  
  all.state.out <-c()
  for(i in 1:length(data_week)){
    dat2 <- src.dat[src.dat$week >= (as.numeric(data_week[i])- 3) & src.dat$week < (as.numeric(data_week[i])+1),]
    
    #make sure dates match dates in proportions dataframe 
    raw_counts_state <- aggregate(count~VARIANT2+STUSAB, data=dat2, FUN=sum,drop=FALSE)
    raw_counts_state$count <- ifelse(is.na(raw_counts_state$count)==T,0,raw_counts_state$count)
    #merge sequence counts with weighted proportions estimates
    all.state2 <- merge(all.state[all.state$Roll_4wk_end==data_week[i],], raw_counts_state,by.x=c("State","Variant"),by.y=c("STUSAB","VARIANT2"),all=T) 
    
    all.state2[is.na(all.state2$count)==T,"count"] <- 0
    
    
    #calculate denominator counts
    dss <- aggregate(count~State, data=all.state2,FUN=sum)
    names(dss)[grep("count",names(dss))] <- "denom_count"
    all.state2 <- merge(all.state2, dss,all=T)
    
    all.state2$Roll_Fourweek_ending <- unique(as.Date(src.dat$yr_wk[src.dat$week==data_week[i]]) + 6)
    
    all.state.out <- rbind.data.frame(all.state.out,all.state2)
  }
  
  #set the Share 0 and CI limits to NA when the count for a lineage is 0
  all.state.out$Share=ifelse(all.state.out$Share!=0 & all.state.out$count==0,0,all.state.out$Share)
  all.state.out$Share_lo=ifelse(is.na(all.state.out$Share_lo)==F & all.state.out$count==0,NA,all.state.out$Share_lo)
  all.state.out$Share_hi=ifelse(is.na(all.state.out$Share_hi)==F & all.state.out$count==0,NA,all.state.out$Share_hi)
  
  #calculate absolute CI width
  all.state.out$CI_width=all.state.out$Share_hi-all.state.out$Share_lo
  
  #generate NCHS flags
  all.state.out$flag_df=ifelse(all.state.out$DF<8,1,0)
  all.state.out$flag_eff.size=ifelse(all.state.out$eff.size<30 | is.na(all.state.out$eff.size)==T,1,0 )
  all.state.out$flag_dss=ifelse(all.state.out$denom_count<30| is.na(all.state.out$denom_count)==T,1,0 )
  all.state.out$flag_abs.ciw=ifelse(all.state.out$CI_width>0.30 | is.na(all.state.out$CI_width)==T,1,0)
  all.state.out$flag_rel.ciw=ifelse(((all.state.out$CI_width/all.state.out$Share)*100)>130 | is.na((all.state.out$CI_width/all.state.out$Share)*100)==T,1,0)
  
  all.state.out$nchs_flag=ifelse(all.state.out$flag_df==1| all.state.out$flag_eff.size==1| all.state.out$denom_count==1| 
                                   all.state.out$flag_abs.ciw==1| all.state.out$flag_rel.ciw==1,1,0)
  
  all.state.out$nchs_flag_wodf=ifelse( all.state.out$flag_eff.size==1| all.state.out$denom_count==1| 
                                         all.state.out$flag_abs.ciw==1| all.state.out$flag_rel.ciw==1,1,0)
  
  
  all.state.out=all.state.out[,c("State","Roll_Fourweek_ending","Variant","Share",           
                                 "Share_lo","Share_hi","count","denom_count","DF","eff.size",
                                 "CI_width","nchs_flag","nchs_flag_wodf")]
  
  
  write.csv(all.state.out, paste0(script.basename, "/results/state_weighted_roll4wk_", ci.type, "CI_svyNEW_", data_date, tag, custom_tag, ".csv"), row.names=FALSE)
}

### Model-smoothed estimates; needs Hessian for regional multinomial model ----------------------------------------

# RUN2
if (length(grep("Run2",tag))!=0){ #Run the nowcast output only when doing the run 2 list of lineages
  
  #Check to see which lineages are in model_vars, add those that need to be added
  AY_agg=names(model_vars)[grep("AY",names(model_vars))]
  
  AY_agg=AY_agg[AY_agg %notin% c("AY.1","AY.2")]
  Other_agg=names(model_vars)[names(model_vars) %notin% voc]
  
  
  #generate a matrix that indicates which lineages to aggregate for the nowcast
  #Columns are the lineages in the nowcast model, so all the defined lineages plus the other lineage
  #Rows are the aggregated lineages wanted
  agg_var_mat <- matrix(data=0,nrow=2,ncol=(length(model_vars)+1))
  colnames(agg_var_mat) <- c(names(model_vars),"Other")
  
  #Fill in matrix values: if lineage needs to be aggregated for parent lineage in given row, then value = 1, else value = 0
  agg_var_mat[1,] <- ifelse(colnames(agg_var_mat) %in% c("B.1.617.2",AY_agg),1,0)
  agg_var_mat[2,] <- ifelse(colnames(agg_var_mat) %in% c(Other_agg,"Other"),1,0)
  row.names(agg_var_mat) <-c("Delta Aggregated","Other Aggregated") 
  
  
  
  
  #define fortnights and regions to get nowcasts for
  proj_ftnts = as.Date(tail(ftnts, 2))
  proj_ftnts = sort(unique(c(proj_ftnts, proj_ftnts + 14, proj_ftnts + 28)))
  
  dfs = expand.grid(USA_or_HHSRegion=c("USA",as.character(1:10)))
  
  proj.res = c()
  for (rgn in dfs$USA_or_HHSRegion) for (ftn in proj_ftnts) {
    if (rgn=="USA") {
      mlm = svymlm_us$mlm
      geoid = rgn } else {
        mlm = svymlm_hhs$mlm
        geoid = as.numeric(rgn)}
    wk = as.numeric(as.Date(ftn, origin="1970-01-01") - week0day1) %/% 7 
    #df = dfs[dfs$USA_or_HHSRegion==rgn, "x"]
    ests = se.multinom(mlm, wk, geoid, composite_variant = agg_var_mat) #gets proportion and se estimate from multinomial model for given week/region
    
    ests = data.frame(
      USA_or_HHSRegion=rgn, 
      Fortnight_ending=as.Date(ftn, origin="1970-01-01"), 
      Variant=c(names(model_vars),"Other",row.names(ests$composite_variant$matrix)), 
      Share=c(ests$p_i,ests$composite_variant$p_i), 
      se.Share=c(ests$se.p_i, ests$composite_variant$se.p_i)#,
      #Share_lo=ifelse( (ests$p_i - (1.96*ests$se.p_i))<0,0,(ests$p_i - (1.96*ests$se.p_i)) ),
      #Share_hi=ifelse( (ests$p_i + (1.96*ests$se.p_i))>1,1,(ests$p_i + (1.96*ests$se.p_i)) )
    )
    #Get binomial CI from p_i and se.p_i
    binom.ci=apply(ests, 1 , function(rr) svyCI(as.numeric(rr[4]),as.numeric(rr[5])))
    ests$Share_lo=binom.ci[1,]
    ests$Share_hi=binom.ci[2,]
    proj.res=rbind(proj.res,ests)
    
  }
  
  proj.res = proj.res[, c("USA_or_HHSRegion", "Fortnight_ending", "Variant", "Share", "Share_lo", "Share_hi")]
  
  #Format output for the run 1 lineage list
  run_1=proj.res[proj.res$Variant %notin% c(AY_agg,"B.1.617.2", Other_agg,"Other"),]
  run_1[run_1$Variant=="Other Aggregated","Variant"] <- "Other"
  
  write.csv(run_1, paste0(script.basename, "/results/updated_nowcast_fortnightly_", data_date,"_state_tag_included_Run1", custom_tag, ".csv"), row.names=FALSE)
  
  #Format output for the run2 lineage list
  drop_lin <- row.names(agg_var_mat)[row.names(agg_var_mat) %notin% "Other Aggregated"]
  run_2=proj.res[proj.res$Variant %notin% c(drop_lin, Other_agg,"Other"),]
  run_2[run_2$Variant=="Other Aggregated","Variant"] <- "Other"
  
  write.csv(run_2, paste0(script.basename, "/results/updated_nowcast_fortnightly_", data_date,"_state_tag_included_Run2", custom_tag, ".csv"), row.names=FALSE)
  
  
  ## Weekly
  
  cast_wks = (function(dd) as.Date(seq(dd[1], dd[2], 7), origin="1970-01-01"))(range(as.Date(proj_ftnts)) + c(-7, 0)) # End-of-week dates
  
  proj.res = c()
  for (rgn in dfs$USA_or_HHSRegion) for (cwk in cast_wks) {
    if (rgn=="USA") {
      mlm = svymlm_us$mlm
      geoid = rgn } else {
        mlm = svymlm_hhs$mlm
        geoid = as.numeric(rgn)}
    wk = as.numeric(as.Date(cwk, origin="1970-01-01") - week0day1) %/% 7 
    #df = dfs[dfs$USA_or_HHSRegion==rgn, "x"]
    ests = se.multinom(mlm, wk, geoid, composite_variant = agg_var_mat) #gets proportion and se estimate from multinomial model for given week/region
    
    ests = data.frame(
      USA_or_HHSRegion=rgn, 
      Week_ending=as.Date(cwk, origin="1970-01-01"), 
      Variant=c(names(model_vars),"Other",row.names(ests$composite_variant$matrix)), 
      Share=c(ests$p_i,ests$composite_variant$p_i), 
      se.Share=c(ests$se.p_i, ests$composite_variant$se.p_i)#,
      #Share_lo=ifelse( (ests$p_i - (1.96*ests$se.p_i))<0,0,(ests$p_i - (1.96*ests$se.p_i)) ),
      #Share_hi=ifelse( (ests$p_i + (1.96*ests$se.p_i))>1,1,(ests$p_i + (1.96*ests$se.p_i)) )
    )
    #Get binomial CI from p_i and se.p_i
    binom.ci=apply(ests, 1 , function(rr) svyCI(as.numeric(rr[4]),as.numeric(rr[5])))
    ests$Share_lo=binom.ci[1,]
    ests$Share_hi=binom.ci[2,]
    proj.res=rbind(proj.res,ests)
    
  }
  
  proj.res = proj.res[, c("USA_or_HHSRegion", "Week_ending", "Variant", "Share", "Share_lo", "Share_hi")]
  
  #Format output for the run 1 lineage list
  run_1=proj.res[proj.res$Variant %notin% c(AY_agg,"B.1.617.2", Other_agg,"Other"),]
  run_1[run_1$Variant=="Other Aggregated","Variant"] <- "Other"
  
  write.csv(run_1, paste0(script.basename, "/results/updated_nowcast_weekly_", data_date,"_state_tag_included_Run1", custom_tag, ".csv"), row.names=FALSE)
  
  #Format output for the run2 lineage list
  drop_lin <- row.names(agg_var_mat)[row.names(agg_var_mat) %notin% "Other Aggregated"]
  run_2=proj.res[proj.res$Variant %notin% c(drop_lin, Other_agg,"Other"),]
  run_2[run_2$Variant=="Other Aggregated","Variant"] <- "Other"
  
  write.csv(run_2, paste0(script.basename, "/results/updated_nowcast_weekly_", data_date,"_state_tag_included_Run2", custom_tag, ".csv"), row.names=FALSE)
  
}

#capture system time
tend=proc.time()

#calculate amount of time code took to run
tdiff=tend[3]-tstart[3]

sprintf(paste0("code for ",tag," run took %f minutes"),tdiff/60)


