# Helper functions used by `weekly_variant_report_nowcast.R`.
# These cover survey-design proportion CIs and the multinomial nowcast model.
myciprop = function(voc,
                    geoid,
                    svy,
                    str   = TRUE,
                    range = FALSE,
                    mut   = FALSE,
                    level = 0.95,
                    ...) {
  # Survey-design CI for one lineage group, returned either as a vector or
  # formatted string.
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

  srv_design = subset(update(svy, VOC = VOC),
                      (STUSAB %in% geoid) |
                        (HHS %in% geoid) |
                        (geoid == "USA"))

  if (ci.type == "xlogit"){
    res = suppressWarnings(
      survey::svyciprop(~VOC,
                        design = srv_design,
                        method = "xlogit",
                        level = level,
                        ...)
    )
  } else {
    res = suppressWarnings(
      svycipropkg(~VOC,
                  design = srv_design,
                  level = level,
                  ...)
    )
  }

  res = c('estimate' = unname(res[1]),
          'lcl' = survey:::confint.svyciprop(res)[1],
          'ucl' = survey:::confint.svyciprop(res)[2])

  res = c(res,
          'DF' = survey::degf(design = srv_design))

  m <- suppressWarnings(
    survey::svymean(x = ~as.numeric(VOC),
                    design = srv_design)
  )

  CVmean <- suppressWarnings(survey::cv(object = m))

  deffmean <- suppressWarnings(
    survey::deff(object = survey::svymean(x = ~as.numeric(VOC),
                                          design = srv_design,
                                          deff = TRUE))
  )

  n.eff <- coef(m) * (1 - coef(m))/vcov(m)
  alpha <- 1-level

  n.eff.capped <- min(
    n.eff * (qt(p = alpha/2,
                df = nrow(srv_design) - 1) /
               qt(p = alpha/2,
                  df = survey::degf(design = srv_design)))^2,
    nrow(srv_design)
  )

  res = c(res,
          'n.eff.capped' = n.eff.capped,
          'CVmean'= CVmean,
          'deffmean' = unname(deffmean))

  if (str) {
    res = c(round( 100 * res[1:3], 1),
            res[4],
            res[5],
            res[6],
            res[7])

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

  res
}

# Fit the multinomial nowcast model and attach a survey-adjusted covariance
# matrix for the coefficients.
svymultinom = function(mod.dat,
                       mysvy,
                       fmla = formula("as.numeric(as.factor(K_US)) ~ model_week + as.factor(HHS)"),
                       model_vars) {
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

  # Aggregate identical design rows before fitting to reduce runtime.
  fmla.vars = all.vars(fmla)
  mlm.dat = data.table::data.table(cbind(data.frame(mod.dat)[, fmla.vars],
                                         weight = weights(mysvy)))[
                                           ,
                                           .(weight = sum(weight)), # aggregate "weight" column
                                           by = fmla.vars] # by the formula

  multinom_geoid = nnet::multinom(formula = fmla,
                                  data    = mlm.dat,
                                  weights = weight,
                                  Hess    = TRUE,
                                  maxit   = 1000,
                                  trace   = FALSE)

  num_var = length(unique(with(mod.dat,
                               eval(terms(fmla)[[2]]))))
  fmla.no.response = formula(delete.response(terms(fmla)))
  formulas = rep(list(fmla.no.response),
                 num_var - 1)
  formulas[[1]] = fmla
  names(formulas) = paste0("b", 1:length(formulas) + 1)

  rval = list(mlm = multinom_geoid,
              estimates = coefficients(multinom_geoid))
  rval$estimates = as.list(data.frame(t(rval$estimates)))
  invinf <- tryCatch(
    {
      solve(-multinom_geoid$Hessian)
    },
    error = function(cond) {
      return(NA)
    }
  )

  if( is.na(invinf[1]) ){ # could use length(invinf) == 1

    rval <- append(rval,
                   list(
                     'scores' = NULL,
                     'invinf' = NULL,
                     'sandwich' = NULL,
                     'SE' = NULL
                   ))

    rval$mlm$Hessian = NULL
    rval$estimates = unlist(rval$estimates)
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
    hess_headtail <- data.frame('element' = names(sort(diag(multinom_geoid$Hessian))),
                                'value' = sort(diag(multinom_geoid$Hessian)))
    rownames(hess_headtail) <- 1:nrow(hess_headtail)
    print(hess_headtail[c(1:5, (nrow(hess_headtail)-4):nrow(hess_headtail)),])

  } else {
    lmultinom = function(v, ...) {
      b = cbind(0, ...)
      b = matrix(b,
                 nrow = length(v))

      sapply(X = seq_along(v),
             FUN = function(rr) b[rr, v[rr]])  - log(rowSums(exp(b)))
    }

    gmultinom = function(v, ...) {
      b = cbind(0, ...)
      b = matrix(b,
                 nrow=length(v))

      p = exp(b) / rowSums(exp(b))
      delta_ij = p * 0
      delta_ij[1:length(v) + (v-1) * length(v)] = 1

      (delta_ij - p)[, -1]
    }

    nms = c("", names(formulas))
    has.response = sapply(formulas, length) == 3
    ff = formulas[[which(has.response)]]
    ff[[3]] = 1
    y = eval.parent(model.frame(formula   = ff,
                                data      = mod.dat,
                                na.action = na.pass))
    formulas[[which(has.response)]] = formula(delete.response(terms(formulas[[which(has.response)]])))
    vnms = unique(do.call(c, lapply(formulas, all.vars)))
    uformula = make.formula(vnms)
    mf = model.frame(formula   = uformula,
                     data      = mod.dat,
                     na.action = na.pass)
    mf = cbind(`(Response)` = y,
               mf)
    mf = mf[, !duplicated(colnames(mf)),
            drop = FALSE]
    weights = weights(mysvy)
    Y = mf[, 1]
    mmFrame = lapply(X    = formulas,
                     FUN  = model.frame,
                     data = mf)
    mm = mapply(FUN      = model.matrix,
                object   = formulas,
                data     = mmFrame,
                SIMPLIFY = FALSE)
    np = c(0,
           cumsum(sapply(X   = mm,
                         FUN = NCOL)))
    parnms = lapply(mm, colnames)
    for (i in 1:length(parnms)){
      parnms[[i]] = paste(nms[i + 1],
                          parnms[[i]],
                          sep = ".")
    }
    parnms = unlist(parnms)
    theta = unlist(rval$estimates)
    args = vector("list", length(nms))
    args[[1]] = Y
    names(args) = nms
    for (i in 2:length(nms)){
      args[[i]] = drop(mm[[i - 1]] %*% theta[(np[i - 1] + 1):np[i]])
    }
    deta = matrix(data = do.call(what = "gmultinom",
                                 args = args),
                  ncol = length(args) - 1)
    rval <- append(rval,
                   list('scores' = NULL))
    reorder = na.omit(match(x = names(formulas),
                            table = nms[-1]))
    for (i in reorder){
      rval$scores = cbind(rval$scores,
                          deta[, i] * weights * mm[[i]])
    }
    rval$invinf = invinf
    dimnames(rval$invinf) = list(parnms, parnms)
    db = rval$scores %*% rval$invinf
    rval$sandwich = survey::svyrecvar(x          = db,
                                      clusters   = mysvy$cluster,
                                      stratas    = mysvy$strata,
                                      fpcs       = mysvy$fpc,
                                      postStrata = mysvy$postStrata)
    rval$SE = sqrt(diag(rval$sandwich))
    rval$estimates = unlist(rval$estimates)
    names(rval$estimates) = names(rval$SE)
    rval$mlm$svyvcov = rval$sandwich
  }

  rval$variants = modvars
  return(rval)
}


# Get predicted proportions and standard errors from a fitted multinomial
# nowcast model, optionally including composite variants.
se.multinom = function(mlm, 
                      newdata_1row,
                      composite_variant=NA,
                      dy_dt=data.frame(model_week=1)) { 
  cf <- coefficients(mlm)
  if ("svyvcov" %in% names(mlm)) {
    vc <- mlm$svyvcov
  } else { 
    if ("Hessian" %in% names(mlm)) {
      vc <- solve(mlm$Hessian)
    } else {
      vc <- matrix(data = 0,
                  nrow = length(cf),
                  ncol = length(cf))
    }
  }
  mm <- model.matrix(as.formula(paste("~",as.character(mlm$terms)[3])),
                    newdata_1row,
                    xlev=mlm$xlevels)
  umm <- mm * 0
  for (nm in names(dy_dt)) umm[, nm] = dy_dt[, nm]
  mnames = outer(X = mlm$lev[-1],
                 Y = colnames(mm),
                 FUN = paste, sep=":")
  cmat = matrix(data = 0,
                nrow = nrow(mnames),
                ncol = ncol(vc),
                dimnames = list(mlm$lev[-1], colnames(mlm$Hessian)))
  for (rr in 1:nrow(cmat)) cmat[rr, mnames[rr,]] = c(mm)
  ucmat = cmat * 0
  for (rr in 1:nrow(ucmat)) ucmat[rr, mnames[rr,]] = c(umm)
  cmat = rbind(`1`=0, cmat)
  ucmat = rbind(`1`=0, ucmat)
  y_i = c(0, coefficients(mlm) %*% c(mm))
  u_i = c(0, coefficients(mlm) %*% c(umm))
  vc.y_i = cmat %*% vc %*% t(cmat)
  vc.z_i = rbind(cmat, ucmat) %*% vc %*% t(rbind(cmat, ucmat))
  p_i = exp(y_i)/sum(exp(y_i))
  g_i = u_i - sum(p_i * u_i)
  dp_dy = diag(p_i) - outer(p_i, p_i, `*`)
  p.vcov = dp_dy %*% vc.y_i %*% dp_dy
  S = diag(nrow = length(y_i))
  if (all(!is.na(composite_variant))) S = rbind(S, composite_variant)
  p.vcov = S %*% p.vcov %*% t(S)
  se.p_S = as.vector(sqrt(diag(p.vcov)))
  p_S = as.vector(S %*% p_i)
  w_S = t(t(S) * p_i)/p_S
  g_S = as.vector(w_S %*% g_i)
  dg_dz = matrix(c(-p_i * g_i, -p_i), nrow=nrow(S), ncol=nrow(vc.z_i), byrow=TRUE) + cbind(w_S * outer(-g_S, g_i, `+`), w_S)
  g.vcov = (dg_dz %*% vc.z_i %*% t(dg_dz))
  se.g_S = as.vector(sqrt(diag(g.vcov)))
  if (all(!is.na(composite_variant))) {
    composite_variant = list(matrix = composite_variant, 
                             p_i = p_S[-(1:length(y_i))],
                             se.p_i = se.p_S[-(1:length(y_i))])
  }
  dim(vc) = rep(rev(dim(cf)), 2)
  res_old = list(p_i=p_i,
                  se.p_i=se.p_S[1:length(y_i)],
                  b_i=c(0, cf[, 2]),
                  se.b_i= c(0, sqrt(diag(vc[2,,2,]))),
                  composite_variant=composite_variant,
                  y_i=y_i,
                  vc.y_i=vc.y_i)
  res_new = list(S=S,
                  p_S=p_S,
                  se.p_S=se.p_S,
                  g_S=g_S,
                  se.g_S=se.g_S,
                  p.vcov=p.vcov,
                  g.vcov=g.vcov)
  c(res_old, res_new)
}

# Wilson-style binomial CI based on a point estimate and standard error.
svyCI = function(p, s, ...) {
  if (s == 0) {
    return(c(0, 0))
  } else if (p == 0) {
    return(c(0, 0))
  } else if (p == 1) {
    return(c(1, 1))
  } else {
    n = p * (1 - p) / s^2

    out = prop.test(
      x = n * p,
      n = n,
      ...
    )$conf.int
    return(c(out[1], out[2]))
  }
}


# Modified `svyciprop` beta-method CI with capped effective sample size.
svycipropkg <- function  (formula, design, level = 0.95, df = degf(design), ...) {
  m <- eval(bquote(svymean(~as.numeric(.(formula[[2]])), design, ...)))
  rval <- coef(m)[1]
  attr(rval, "var") <- vcov(m)
  alpha <- 1 - level

  if (!is.null(design[['postStrata']])){
    if (!as.character(design[['call']][[1]]) == "svystandardize"){
      stop("svycipropkg design cannot be a subset of an age-standardized survey design object")
    }

    design_tmp <- design[['variables']][which(design[['prob']] != Inf), ]
    ageadjvar <- as.character(design[['call']]$by[[2]])
    population <- eval(design[['call']]$population)
    p <- lapply(split(design_tmp, design_tmp[[ageadjvar]]),
                function(x) weighted.mean(x[[as.character(formula[[2]])]],
                                          x[[names(design[['allprob']])]] ) )
    n <- lapply(split(design_tmp, design_tmp[[ageadjvar]]), nrow )
    p <- unlist(p)
    n <- unlist(n)
    age_var <- p*(1-p)/n
    pop <- population/sum(population)
    varsrs_adj <- sum(pop^2 * age_var)
    deff_adj <- ifelse(sum(p) == 0, 1, vcov(m)/varsrs_adj)
    n.eff <- ifelse(rval == 0, nrow(design_tmp),
                    min(nrow(design_tmp),
                        nrow(design_tmp)/deff_adj * (qt(alpha/2, nrow(design_tmp) - 1)/qt(alpha/2, degf(design)))^2))

  }
  else {
    n.eff <- coef(m) * (1 - coef(m))/vcov(m)
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
