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

    