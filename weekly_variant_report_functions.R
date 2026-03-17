# Helper functions used by `weekly_variant_report_nowcast.R`.

shares_sum_to_one <- function(dt, group_cols) {
  all(dt[, .(total_share = sum(Share)), by = group_cols][, unique(round(total_share, 5))] == 1)
}

normalize_hadoop_variant_names <- function(values) {
  values <- gsub("Delta Aggregated", "B.1.617.2", values)
  values <- gsub("Omicron Aggregated", "B.1.1.529", values)
  gsub(" Aggregated", "", values)
}

build_hadoop_output <- function(results_df,
                                base_columns,
                                data_date,
                                results_tag,
                                run_label,
                                cadence_label,
                                calc_confirmed_infections,
                                confirmed_columns) {
  hadoop_df <- data.frame(results_df[, base_columns])
  hadoop_df[is.na(hadoop_df)] <- "\\N"
  hadoop_df[, 7:15] <- "\\N"
  hadoop_df[, 16] <- "smoothed"
  hadoop_df[, 17] <- cadence_label
  hadoop_df[, 18] <- data_date
  hadoop_df[, 19] <- paste0(results_tag, "_", run_label)
  hadoop_df[, 20] <- 1

  if (calc_confirmed_infections) {
    hadoop_df[, 21:24] <- results_df[, confirmed_columns]
  }

  hadoop_df[, 3] <- normalize_hadoop_variant_names(hadoop_df[, 3])
  hadoop_df[is.na(hadoop_df)] <- "\\N"
  hadoop_df
}

write_hadoop_output <- function(results_df,
                                file,
                                base_columns,
                                data_date,
                                results_tag,
                                run_label,
                                cadence_label,
                                calc_confirmed_infections,
                                confirmed_columns) {
  hadoop_df <- build_hadoop_output(
    results_df = results_df,
    base_columns = base_columns,
    data_date = data_date,
    results_tag = results_tag,
    run_label = run_label,
    cadence_label = cadence_label,
    calc_confirmed_infections = calc_confirmed_infections,
    confirmed_columns = confirmed_columns
  )

  write.table(
    x = hadoop_df,
    file = file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = ","
  )
}

build_empiric_hadoop_output <- function(results_df,
                                        base_columns,
                                        share_columns,
                                        estimate_label,
                                        cadence_label,
                                        data_date,
                                        results_tag,
                                        run_number,
                                        calc_confirmed_infections = FALSE,
                                        confirmed_columns = c("cases", "cases_hi", "cases_lo")) {
  hadoop_df <- data.frame(results_df[, base_columns])
  hadoop_df[, 4:6] <- results_df[, share_columns]
  hadoop_df[is.na(hadoop_df)] <- "\\N"
  hadoop_df[, 14:15] <- "\\N"
  hadoop_df[, 16] <- estimate_label
  hadoop_df[, 17] <- cadence_label
  hadoop_df[, 18] <- data_date
  hadoop_df[, 19] <- paste0(results_tag, "_Run", run_number)
  hadoop_df[, 20] <- 1

  if (calc_confirmed_infections) {
    hadoop_df[, 21:23] <- results_df[, confirmed_columns]
  }

  hadoop_df[, 3] <- normalize_hadoop_variant_names(hadoop_df[, 3])
  hadoop_df[is.na(hadoop_df)] <- "\\N"
  hadoop_df
}

write_empiric_output <- function(results_df,
                                 csv_file,
                                 hadoop_file,
                                 base_columns,
                                 share_columns,
                                 estimate_label,
                                 cadence_label,
                                 data_date,
                                 results_tag,
                                 run_number,
                                 calc_confirmed_infections = FALSE,
                                 confirmed_columns = c("cases", "cases_hi", "cases_lo")) {
  write.csv(x = results_df, file = csv_file, row.names = FALSE)

  hadoop_df <- build_empiric_hadoop_output(
    results_df = results_df,
    base_columns = base_columns,
    share_columns = share_columns,
    estimate_label = estimate_label,
    cadence_label = cadence_label,
    data_date = data_date,
    results_tag = results_tag,
    run_number = run_number,
    calc_confirmed_infections = calc_confirmed_infections,
    confirmed_columns = confirmed_columns
  )

  write.table(
    x = hadoop_df,
    file = hadoop_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = ","
  )
}

get_run1_lineages <- function(agg_var_mat) {
  agg_lineages <- colnames(agg_var_mat)[colSums(agg_var_mat) > 0]
  if ("Other" %notin% agg_lineages) {
    agg_lineages <- c(agg_lineages, "Other Aggregated")
  }
  agg_lineages
}

get_run2_excluded_lineages <- function(agg_var_mat) {
  c(
    row.names(agg_var_mat)[row.names(agg_var_mat) %notin% "Other Aggregated"],
    "Other",
    colnames(agg_var_mat["Other Aggregated", , drop = FALSE])[agg_var_mat["Other Aggregated", , drop = FALSE] > 0]
  )
}

rename_other_aggregated <- function(results_df) {
  results_df[results_df$Variant == "Other Aggregated", "Variant"] <- "Other"
  results_df
}

build_run1_output <- function(proj_res, agg_var_mat) {
  rename_other_aggregated(proj_res[Variant %notin% get_run1_lineages(agg_var_mat)])
}

build_run2_output <- function(proj_res, agg_var_mat) {
  rename_other_aggregated(proj_res[Variant %notin% get_run2_excluded_lineages(agg_var_mat)])
}

split_weekly_and_daily_results <- function(results_df) {
  list(
    weekly = data.table:::subset.data.table(
      x = results_df,
      subset = model_week %% 1 == 0,
      select = !names(results_df) %in% c("total_test_positives_daily", "cases_daily", "cases_lo_daily", "cases_hi_daily")
    ),
    daily = data.table:::subset.data.table(
      x = results_df,
      select = !names(results_df) %in% c("total_test_positives_weekly", "cases_weekly", "cases_lo_weekly", "cases_hi_weekly")
    )
  )
}

write_validated_output <- function(results_df,
                                   group_cols,
                                   warning_file,
                                   invalid_message,
                                   csv_file,
                                   hadoop_file = NULL,
                                   hadoop_base_columns = NULL,
                                   data_date = NULL,
                                   results_tag = NULL,
                                   run_label = NULL,
                                   cadence_label = NULL,
                                   calc_confirmed_infections = FALSE,
                                   confirmed_columns = NULL) {
  if (!shares_sum_to_one(results_df, group_cols)) {
    warning(paste(warning_file, invalid_message))
    return(invisible(FALSE))
  }

  write.csv(x = results_df, file = csv_file, row.names = FALSE)

  if (!is.null(hadoop_file)) {
    write_hadoop_output(
      results_df = results_df,
      file = hadoop_file,
      base_columns = hadoop_base_columns,
      data_date = data_date,
      results_tag = results_tag,
      run_label = run_label,
      cadence_label = cadence_label,
      calc_confirmed_infections = calc_confirmed_infections,
      confirmed_columns = confirmed_columns
    )
  }

  invisible(TRUE)
}

results_file_path <- function(script.basename, output_folder, filename) {
  paste0(script.basename, output_folder, "/", filename)
}

write_results_csv <- function(x, script.basename, output_folder, filename, row.names = FALSE, ...) {
  write.csv(
    x = x,
    file = results_file_path(
      script.basename = script.basename,
      output_folder = output_folder,
      filename = filename
    ),
    row.names = row.names,
    ...
  )
}

write_lineage_aggregation_tables <- function(src.dat,
                                             script.basename,
                                             output_folder,
                                             ci.type,
                                             svy.type,
                                             data_date,
                                             tag) {
  write.csv(
    x = src.dat[
      ,
      .(old = lineage, new = VARIANT),
      by = c("lineage", "VARIANT")
    ][
      ,
      .(old, new)
    ][
      order(old),
    ],
    file = results_file_path(
      script.basename = script.basename,
      output_folder = output_folder,
      filename = paste0(
        "lineage_aggregations_",
        ci.type,
        "CI_",
        svy.type,
        "_",
        data_date,
        tag,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    x = src.dat[
      ,
      .(original = lineage, aggregated = VARIANT2),
      by = c("lineage", "VARIANT2")
    ][
      ,
      .(original, aggregated)
    ][
      order(original),
    ],
    file = results_file_path(
      script.basename = script.basename,
      output_folder = output_folder,
      filename = paste0(
        "lineage_aggregations_summary_",
        ci.type,
        "CI_",
        svy.type,
        "_",
        data_date,
        tag,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

save_nowcast_model_artifacts <- function(script.basename,
                                         output_folder,
                                         data_date,
                                         tag,
                                         svymlm_hhs,
                                         svymlm_us,
                                         src.moddat = NULL) {
  saveRDS(
    object = svymlm_hhs,
    file = results_file_path(
      script.basename = script.basename,
      output_folder = output_folder,
      filename = paste0("svymlm_hhs_", data_date, tag, ".RDS")
    )
  )

  saveRDS(
    object = svymlm_us,
    file = results_file_path(
      script.basename = script.basename,
      output_folder = output_folder,
      filename = paste0("svymlm_us_", data_date, tag, ".RDS")
    )
  )

  if (!is.null(src.moddat)) {
    saveRDS(
      object = src.moddat,
      file = results_file_path(
        script.basename = script.basename,
        output_folder = output_folder,
        filename = paste0("src.moddat_", data_date, tag, ".RDS")
      )
    )
  }
}

add_share_ci_columns <- function(results_df, calc_99 = FALSE) {
  binom_ci <- apply(
    X = results_df,
    MARGIN = 1,
    FUN = function(rr) svyCI(p = as.numeric(rr[4]), s = as.numeric(rr[5]))
  )

  results_df$Share_lo <- binom_ci[1, ]
  results_df$Share_hi <- binom_ci[2, ]

  if (calc_99) {
    binom_ci_99 <- array(data = numeric(), dim = c(2, nrow(results_df)))
    for (rr in seq_len(nrow(results_df))) {
      binom_ci_99[, rr] <- svyCI(
        p = as.numeric(results_df[rr, 4]),
        s = as.numeric(results_df[rr, 5]),
        conf.level = 0.99
      )
    }

    results_df$Share_lo_99 <- binom_ci_99[1, ]
    results_df$Share_hi_99 <- binom_ci_99[2, ]
  }

  results_df
}

compute_growth_metrics_from_link <- function(link, se_link) {
  growth_rate <- 100 * exp(link) - 100
  growth_rate_lo <- 100 * exp(link - 1.96 * se_link) - 100
  growth_rate_hi <- 100 * exp(link + 1.96 * se_link) - 100

  list(
    se_growth = 100 * exp(se_link) - 100,
    growth_rate = growth_rate,
    growth_rate_lo = growth_rate_lo,
    growth_rate_hi = growth_rate_hi,
    doubling_time = log(2) / link * 7,
    doubling_time_lo = log(2) / (link - 1.96 * se_link) * 7,
    doubling_time_hi = log(2) / (link + 1.96 * se_link) * 7
  )
}

compute_variant_growth_metrics <- function(ests) {
  se_link <- with(
    data = ests,
    expr = sqrt(se.b_i^2 * (1 - 2 * p_i) + sum(se.p_i^2 * b_i^2 + p_i^2 * se.b_i^2))
  )
  link <- with(data = ests, expr = (b_i - sum(p_i * b_i)))
  metrics <- compute_growth_metrics_from_link(link = link, se_link = se_link)
  metrics$growth_rate_link <- link
  metrics$se_growth_link <- se_link
  metrics
}

aggregate_growth_metrics <- function(agg_var_mat, base_share, metrics) {
  gr_agg <- data.frame(
    variant = rownames(agg_var_mat),
    gr = NA,
    se.gr = NA,
    gr_lo = NA,
    gr_hi = NA,
    dt = NA,
    dt_lo = NA,
    dt_hi = NA
  )

  for (r in seq_len(nrow(agg_var_mat))) {
    col_ind <- which(agg_var_mat[r, ] > 0)

    if (length(col_ind) == 1) {
      gr_agg[r, "gr"] <- metrics$growth_rate[col_ind]
      gr_agg[r, "se.gr"] <- metrics$se_growth[col_ind]
      gr_agg[r, "gr_lo"] <- metrics$growth_rate_lo[col_ind]
      gr_agg[r, "gr_hi"] <- metrics$growth_rate_hi[col_ind]
      gr_agg[r, "dt"] <- metrics$doubling_time[col_ind]
      gr_agg[r, "dt_lo"] <- metrics$doubling_time_lo[col_ind]
      gr_agg[r, "dt_hi"] <- metrics$doubling_time_hi[col_ind]
    } else {
      weights <- base_share[col_ind]
      gr_agg[r, "gr"] <- sum(metrics$growth_rate[col_ind] * weights) / sum(weights)
      gr_agg[r, "gr_lo"] <- sum(metrics$growth_rate_lo[col_ind] * weights) / sum(weights)
      gr_agg[r, "gr_hi"] <- sum(metrics$growth_rate_hi[col_ind] * weights) / sum(weights)
      gr_agg[r, "dt"] <- sum(metrics$doubling_time[col_ind] * weights) / sum(weights)
      gr_agg[r, "dt_lo"] <- sum(metrics$doubling_time_lo[col_ind] * weights) / sum(weights)
      gr_agg[r, "dt_hi"] <- sum(metrics$doubling_time_hi[col_ind] * weights) / sum(weights)
    }
  }

  gr_agg
}

confirmed_infections_filepath <- function(script.basename, data_date, custom_tag) {
  paste0(
    script.basename,
    "/data/backup_",
    data_date,
    custom_tag,
    "/",
    data_date,
    "_tests_aggregated",
    custom_tag,
    ".RDS"
  )
}

attach_monthly_confirmed_infections <- function(proj.res, script.basename, data_date, custom_tag) {
  test_filepath <- confirmed_infections_filepath(script.basename, data_date, custom_tag)

  if (!file.exists(test_filepath)) {
    print(paste0(
      "File ",
      test_filepath,
      " not found. Not calculating number of infections attributable to each variant for months."
    ))
    return(proj.res)
  }

  test_list <- readRDS(file = test_filepath)
  tests_monthly <- rbind(
    test_list$tests_weekly[, .("total_test_positives" = sum(POSITIVE, na.rm = TRUE)), by = "month_end"][, "HHS" := "USA"],
    test_list$tests_weekly[, .("total_test_positives" = sum(POSITIVE, na.rm = TRUE)), by = c("month_end", "HHS")]
  )

  proj.res <- merge(
    x = proj.res,
    y = tests_monthly,
    by.x = c("USA_or_HHSRegion", "Month_ending"),
    by.y = c("HHS", "month_end"),
    all.x = TRUE
  )

  proj.res[, cases := total_test_positives * Share]
  proj.res[, cases_lo := total_test_positives * Share_lo]
  proj.res[, cases_hi := total_test_positives * Share_hi]
  proj.res
}

attach_weekly_confirmed_infections <- function(proj.res, script.basename, data_date, custom_tag) {
  test_filepath <- confirmed_infections_filepath(script.basename, data_date, custom_tag)

  if (!file.exists(test_filepath)) {
    print(paste0(
      "File ",
      test_filepath,
      " not found. Not calculating number of infections attributable to each variant for weeks."
    ))
    return(proj.res)
  }

  test_list <- readRDS(file = test_filepath)

  tests_dy_us <- test_list$tests_daily[, .("total_test_positives_daily" = sum(POSITIVE_daily, na.rm = TRUE)), by = "date"][, "HHS" := "USA"]
  tests_dy_hhs <- test_list$tests_daily[, .("total_test_positives_daily" = sum(POSITIVE_daily, na.rm = TRUE)), by = c("date", "HHS")]
  tests_wk_us <- test_list$tests_weekly[, .("total_test_positives_weekly" = sum(POSITIVE, na.rm = TRUE)), by = "yr_wk"][, "HHS" := "USA"]
  tests_wk_hhs <- test_list$tests_weekly[, .("total_test_positives_weekly" = sum(POSITIVE, na.rm = TRUE)), by = c("yr_wk", "HHS")]

  proj.res <- merge(
    x = proj.res,
    y = rbind(tests_dy_us, tests_dy_hhs),
    by.x = c("USA_or_HHSRegion", "date"),
    by.y = c("HHS", "date"),
    all.x = TRUE
  )

  proj.res <- merge(
    x = proj.res,
    y = rbind(tests_wk_us, tests_wk_hhs)[, "Week_ending" := as.Date(yr_wk) + 6][, "yr_wk" := NULL],
    by.x = c("USA_or_HHSRegion", "Week_ending"),
    by.y = c("HHS", "Week_ending"),
    all.x = TRUE
  )

  proj.res[, cases_daily := total_test_positives_daily * Share]
  proj.res[, cases_lo_daily := total_test_positives_daily * Share_lo]
  proj.res[, cases_hi_daily := total_test_positives_daily * Share_hi]
  proj.res[, cases_weekly := total_test_positives_weekly * Share]
  proj.res[, cases_lo_weekly := total_test_positives_weekly * Share_lo]
  proj.res[, cases_hi_weekly := total_test_positives_weekly * Share_hi]
  proj.res
}

add_doubling_time_axis <- function(side = 4, line = 3) {
  plot_growth_rates <- axTicks(2)
  plot_doubling_times <- (log(2) / log((100 + plot_growth_rates) / 100)) * 7
  axis(side = side, at = plot_growth_rates, labels = round(plot_doubling_times, 1))
  mtext("Doubling time (days)", side = side, line = line)
}

build_aggregation_matrix <- function(voc_aggregation_method,
                                     model_vars,
                                     voc1,
                                     voc,
                                     voc_lut,
                                     script.basename,
                                     output_folder,
                                     ci.type,
                                     svy.type,
                                     data_date,
                                     tag) {
  if (tolower(voc_aggregation_method) %in% c("updated", "lineage_expanded")) {
    agg_var_mat <- matrix(data = 0, nrow = 0, ncol = (length(model_vars) + 1))
    colnames(agg_var_mat) <- c(model_vars, "Other")

    model_vars_expanded <- setNames(voc_lut$lineage_expanded, voc_lut$variant)[model_vars[model_vars != "Other"]]
    voc1_expanded <- setNames(voc_lut$lineage_expanded, voc_lut$variant)[voc1]
    voc1_expanded <- subset(voc1_expanded, !is.na(voc1_expanded))
    model_var_parents_expanded <- nearest_parent(model_vars_expanded, voc1_expanded)
    model_var_parents <- setNames(voc_lut$variant, voc_lut$lineage_expanded)[model_var_parents_expanded]

    model_var_lut <- data.frame(voc2_model_vars = model_vars, voc1_model_vars = model_var_parents)
    mvl_sub <- model_var_lut[model_var_lut$voc2_model_vars != model_var_lut$voc1_model_vars, ]
    unique_mvl_sub <- unique(mvl_sub$voc1_model_vars)

    for (i in unique_mvl_sub) {
      extra_row <- ifelse(
        colnames(agg_var_mat) %in% c(i, model_var_lut[model_var_lut$voc1_model_vars == i, "voc2_model_vars"]),
        1,
        0
      )
      agg_var_mat <- rbind(agg_var_mat, extra_row)
      row.names(agg_var_mat)[nrow(agg_var_mat)] <- paste(i, "Aggregated")
    }

    other_agg <- base::setdiff(colnames(agg_var_mat)[colSums(agg_var_mat) == 0], voc1[voc1 %in% model_vars])
    agg_var_mat <- rbind(agg_var_mat, ifelse(colnames(agg_var_mat) %in% other_agg, 1, 0))
    row.names(agg_var_mat)[nrow(agg_var_mat)] <- "Other Aggregated"

    if (max(colSums(agg_var_mat)) > 1) {
      warning(message = paste(
        "Aggregated results are invalid! These variants are being aggregated multiple times:",
        names(agg_var_mat)[colSums(agg_var_mat) > 1],
        ". Fix the aggregation matrix."
      ))
    }

    sub_mat <- agg_var_mat[row.names(agg_var_mat) != "Other Aggregated", setdiff(model_vars, voc1), drop = FALSE]
    if (!all(colSums(sub_mat) > 0)) {
      problem_vocs <- colnames(sub_mat)[colSums(sub_mat) == 0]
      warning(message = paste0(
        'Not all variants in "voc" are aggregated into something in "voc1". ',
        paste(problem_vocs, collapse = ", "),
        ' is in voc2 but is not aggregated into anything in voc1. Therefore it will be aggregated into "Other Aggregated" for both run_1 and run_2 output. Either change voc1 and/or voc2 or change the way agg_var_mat works.'
      ))
    }
  } else {
    AY_vars <- model_vars[grep("^AY\\.", model_vars, perl = TRUE)]
    BA_vars <- model_vars[grep("(^B[AC-HJ-NP-VYZ]\\.)|(^C[A-HJ-NP-WYZ]\\.)|(^D[A-HJ-NP-WYZ]\\.)|(^E[A-FHJNPQRSTVWYZ]\\.)|(^F[ABCFJKMN]\\.)", model_vars, perl = TRUE)]
    XBB_vars <- model_vars[grep("(^XBB\\.)|^E[GKLMU]\\.|^F[DEGHL]\\.", model_vars, perl = TRUE)]
    run1_lineages <- voc1
    AY_agg <- AY_vars[AY_vars %notin% run1_lineages]
    BA_agg <- BA_vars[BA_vars %notin% run1_lineages]
    XBB_agg <- XBB_vars[XBB_vars %notin% run1_lineages]
    Other_agg <- model_vars[model_vars %notin% c(voc, "B.1.617.2", "B.1.1.529")]

    agg_var_mat <- matrix(data = 0, nrow = 1, ncol = (length(model_vars) + 1))
    colnames(agg_var_mat) <- c(model_vars, "Other")
    agg_var_mat[1, ] <- ifelse(colnames(agg_var_mat) %in% c(Other_agg, "Other"), 1, 0)
    row.names(agg_var_mat) <- c("Other Aggregated")

    if (length(AY_agg) != 0) {
      extra_row <- ifelse(colnames(agg_var_mat) %in% c("B.1.617.2", AY_agg), 1, 0)
      agg_var_mat <- rbind(agg_var_mat, extra_row)
      row.names(agg_var_mat)[nrow(agg_var_mat)] <- c("Delta Aggregated")
    }
    if (length(BA_agg) != 0) {
      extra_row <- ifelse(colnames(agg_var_mat) %in% c("B.1.1.529", BA_agg), 1, 0)
      agg_var_mat <- rbind(agg_var_mat, extra_row)
      row.names(agg_var_mat)[nrow(agg_var_mat)] <- c("Omicron Aggregated")
    }
    if (length(XBB_agg) != 0) {
      extra_row <- ifelse(colnames(agg_var_mat) %in% c("XBB", XBB_agg), 1, 0)
      agg_var_mat <- rbind(agg_var_mat, extra_row)
      row.names(agg_var_mat)[nrow(agg_var_mat)] <- c("XBB Aggregated")
    }

    XBBs_in_r1l <- XBB_vars[XBB_vars %in% run1_lineages]
    for (ll in XBBs_in_r1l) {
      ll_agg <- NULL
      if (ll == "XBB.1.5") {
        ll_agg <- grep("^XBB\\.1\\.5(?![0-9])|^FD\\.", XBB_vars, perl = TRUE, value = TRUE)
        ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
        ll_agg <- c(ll_agg, "XBB.1.5")
      }
      if (ll == "XBB.1.16") {
        ll_agg <- grep("^XBB\\.1\\.16(?![0-9])", XBB_vars, perl = TRUE, value = TRUE)
        ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
        ll_agg <- c(ll_agg, "XBB.1.16")
      }
      if (ll == "XBB.1.9.2") {
        ll_agg <- grep("^XBB\\.1\\.9\\.2(?![0-9])|^EG\\.", XBB_vars, perl = TRUE, value = TRUE)
        ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
        ll_agg <- c(ll_agg, "XBB.1.9.2")
      }
      if (ll == "XBB.2.3") {
        ll_agg <- grep("^XBB\\.2\\.3(?![0-9])", XBB_vars, perl = TRUE, value = TRUE)
        ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
        ll_agg <- c(ll_agg, "XBB.2.3")
      }
      if (!is.null(ll_agg) && any(ll_agg != ll)) {
        extra_row <- ifelse(colnames(agg_var_mat) %in% ll_agg, 1, 0)
        agg_var_mat <- rbind(agg_var_mat, extra_row)
        row.names(agg_var_mat)[nrow(agg_var_mat)] <- paste(ll, "Aggregated")
        om_row <- which(row.names(agg_var_mat) == "XBB Aggregated")
        agg_var_mat[om_row, which(agg_var_mat[om_row, ] == 1 & extra_row == 1)] <- 0
      }
    }

    BAs_in_r1l <- BA_vars[BA_vars %in% run1_lineages]
    for (ll in BAs_in_r1l) {
      ll_agg <- NULL
      if (ll == "BA.1") ll_agg <- grep("(^BA\\.1)(?!(\\.1$))", BA_vars, perl = TRUE, value = TRUE)
      if (ll == "BA.1.1") ll_agg <- grep("(^BA\\.1\\.1)(?![0-9])", BA_vars, perl = TRUE, value = TRUE)
      if (ll == "BA.2") {
        ll_agg <- grep("(^BA\\.2)(?![0-9])", BA_vars, perl = TRUE, value = TRUE)
        ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
        ll_agg <- c(ll_agg, "BA.2")
      }
      if (ll == "BA.2.12.1") ll_agg <- grep("(^BA\\.2\\.12\\.1)|(^BG\\.)", BA_vars, perl = TRUE, value = TRUE)
      if (ll == "BA.2.75") {
        ll_agg <- grep("(^BA\\.2\\.75)(?![0-9])|(^BN\\.)|(^CH\\.)", BA_vars, perl = TRUE, value = TRUE)
        if (length(grep("(^CH\\.1\\.1)", run1_lineages, perl = TRUE, value = TRUE)) == 1) {
          ll_agg_not <- grep("(^CH\\.1\\.1)(?![0-9])", ll_agg, perl = TRUE, value = TRUE)
          ll_agg <- setdiff(ll_agg, ll_agg_not)
        }
        ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
        ll_agg <- c(ll_agg, "BA.2.75")
      }
      if (ll == "CH.1.1") {
        ll_agg <- grep("(^CH\\.1\\.1)(?![0-9])", BA_vars, perl = TRUE, value = TRUE)
        ll_agg <- c(ll_agg, "CH.1.1")
      }
      if (ll == "BA.3") ll_agg <- grep("(^BA\\.3)(?![0-9])", BA_vars, perl = TRUE, value = TRUE)
      if (ll == "BA.4") {
        ll_agg <- grep("(^BA\\.4)(?![0-9])", BA_vars, perl = TRUE, value = TRUE)
        if (length(grep("(^BA\\.4\\.6)", run1_lineages, perl = TRUE, value = TRUE)) == 1) {
          ll_agg <- grep("(^BA\\.4)(?!([0-9])|(\\.6))", BA_vars, perl = TRUE, value = TRUE)
        }
      }
      if (ll == "BA.4.6") ll_agg <- grep("(^BA\\.4\\.6)(?![0-9])", BA_vars, perl = TRUE, value = TRUE)
      if (ll == "BQ.1.1") ll_agg <- grep("(^BQ\\.1\\.1)(?![0-9])|(^DK\\.)", BA_vars, perl = TRUE, value = TRUE)
      if (ll == "BQ.1") {
        if (length(grep("^BQ\\.1\\.1(?![0-9])", run1_lineages, perl = TRUE, value = TRUE)) > 0) {
          ll_agg <- grep("(^BQ\\.1)(?![0-9]|(\\.1(?![0-9])))", BA_vars, perl = TRUE, value = TRUE)
        } else {
          ll_agg <- grep("(^BQ\\.1)(?![0-9])", BA_vars, perl = TRUE, value = TRUE)
        }
        ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
        ll_agg <- c(ll_agg, "BQ.1")
      }
      if (ll == "BF.11") ll_agg <- grep("(^BF\\.11)(?![0-9])", BA_vars, perl = TRUE, value = TRUE)
      if (ll == "BF.7") ll_agg <- grep("(^BF\\.7)(?![0-9])", BA_vars, perl = TRUE, value = TRUE)
      if (ll == "BA.5") {
        if (length(grep("(^BF\\.7)(?![0-9])", run1_lineages, perl = TRUE, value = TRUE)) > 0) {
          ll_agg <- grep("(^CQ\\.)|(^BA\\.5)(?![0-9])|(^BE\\.)|(^BF\\.(?!7(?![0-9])))|(^C[KR]\\.)|(^DF\\.)", BA_vars, perl = TRUE, value = TRUE)
        } else {
          ll_agg <- grep("(^CQ\\.)|(^BA\\.5)(?![0-9])|(^B[EF]\\.)|(^C[KR]\\.)|(^DF\\.)", BA_vars, perl = TRUE, value = TRUE)
        }
        ll_agg <- setdiff(ll_agg, ll_agg[ll_agg %in% run1_lineages])
        ll_agg <- c(ll_agg, "BA.5")
      }
      if (!is.null(ll_agg) && any(ll_agg != ll)) {
        extra_row <- ifelse(colnames(agg_var_mat) %in% ll_agg, 1, 0)
        agg_var_mat <- rbind(agg_var_mat, extra_row)
        row.names(agg_var_mat)[nrow(agg_var_mat)] <- paste(ll, "Aggregated")
        om_row <- which(row.names(agg_var_mat) == "Omicron Aggregated")
        agg_var_mat[om_row, which(agg_var_mat[om_row, ] == 1 & extra_row == 1)] <- 0
      }
    }

    if (any(colSums(agg_var_mat) > 1)) {
      warning(paste0("agg_var_mat not correctly specified. Some variants are aggregated more than once.", agg_var_mat))
    }

    agg_var_names <- gsub(" Aggregated", "", sub("Delta", "B.1.617.2", sub("Omicron", "B.1.1.529", row.names(agg_var_mat))))
    if (!all(diag(agg_var_mat[, agg_var_names]) == 1)) {
      stop(paste0(
        "In 'agg_var_mat', ",
        agg_var_names[diag(agg_var_mat[, agg_var_names]) != 1],
        " is not included in ",
        row.names(agg_var_mat)[diag(agg_var_mat[, agg_var_names]) != 1]
      ))
    }
  }

  write.csv(
    x = replace(agg_var_mat, agg_var_mat == 0, NA),
    file = paste0(script.basename, output_folder, "/agg_var_mat_", ci.type, "CI_", svy.type, "_", data_date, tag, ".csv"),
    row.names = TRUE,
    na = ""
  )

  agg_var_mat
}

build_growth_rate_table <- function(summary_obj, growth_metrics, model_vars, data_week_df) {
  gr_tab <- data.frame(
    variant = c(model_vars, "OTHER"),
    variant_share = 100 * summary_obj$p_i,
    variant_share_lo = 100 * (summary_obj$p_i - 1.96 * summary_obj$se.p_i),
    variant_share_hi = 100 * (summary_obj$p_i + 1.96 * summary_obj$se.p_i),
    growth_rate = growth_metrics$growth_rate,
    growth_rate_lo = growth_metrics$growth_rate_lo,
    growth_rate_hi = growth_metrics$growth_rate_hi,
    doubling_time = growth_metrics$doubling_time,
    doubling_time_lo = growth_metrics$doubling_time_lo,
    doubling_time_hi = growth_metrics$doubling_time_hi,
    model_week = data_week_df$model_week
  )

  merge(gr_tab, data_week_df, by = "model_week")
}

plot_weighted_share_barplot_us <- function(bp_us,
                                           display_vars,
                                           src.dat,
                                           current_week,
                                           display_lookback,
                                           fig_gen_run,
                                           filename) {
  if (fig_gen_run) {
    jpeg(filename = filename, width = 1500, height = 1500, pointsize = 40)
    on.exit(dev.off(), add = TRUE)
  }

  bp <- barplot(
    height = 100 * t(bp_us),
    xlab = "Week beginning",
    ylab = "Weighted variant share (%)",
    main = "Nationwide",
    border = NA,
    ylim = 110 * 0:1,
    col = col.dk,
    names.arg = rownames(bp_us),
    legend.text = display_vars,
    args.legend = list(x = "topleft", bty = "n", border = NA)
  )

  text(
    x = bp,
    y = 3 + colSums(100 * t(tail(bp_us, 12))),
    labels = with(
      subset(src.dat, week < current_week - 1 & week >= current_week - display_lookback),
      table(week)
    ),
    cex = 0.7
  )
}

plot_projection_barplot_us <- function(pred_us.df,
                                       display_indices,
                                       display_vars,
                                       data_week_df,
                                       model_week_min,
                                       week0day1,
                                       fig_gen_run,
                                       filename) {
  if (fig_gen_run) {
    jpeg(filename = filename, width = 1500, height = 1500, pointsize = 40)
    on.exit(dev.off(), add = TRUE)
  }

  bp <- barplot(
    height = 100 * t(pred_us.df[, 1 + display_indices]),
    xlab = "Week beginning",
    ylab = "Weighted variant share (%)",
    main = "Nationwide",
    space = 0,
    border = NA,
    ylim = 110 * 0:1,
    col = col.dk,
    names.arg = ifelse(
      test = pred_us.df$model_week %% 1 == 0,
      yes = format(((pred_us.df$model_week + (data_week_df$model_week - 2)) + model_week_min) * 7 + as.Date(week0day1), format = "%m-%d"),
      no = NA
    ),
    legend.text = display_vars,
    args.legend = list(x = "topleft", bty = "n", border = NA)
  )

  pc <- unlist(100 * subset(pred_us.df, model_week == data_week_df$model_week)[, 1 + display_indices])
  y <- cumsum(pc) - pc/2
  x <- bp[which(pred_us.df$model_week %in% (data_week_df$model_week + c(-2, 0, 2)))]

  text(x = x[2], y = y, labels = round(pc, 1), cex = 0.7, xpd = TRUE, adj = c(0.5, 0.5))
  rect(xleft = x[1] + (x[2] - x[1]) * .25, ybottom = 0, xright = x[2] + (x[3] - x[2]) * .25, ytop = 100, border = NA, col = "#00000020")
  text(x = mean(c(x[1], x[2])) + (x[2] - x[1]) * .25, y = 101, labels = "Nowcast", cex = 0.7, xpd = TRUE, adj = c(0.5, 0))
  rect(xleft = x[2] + (x[3] - x[2]) * .25, ybottom = 0, xright = x[3], ytop = 100, border = NA, col = "#00000040")
  text(x = mean(c(x[2], x[3])) + (x[2] - x[1]) * .125, y = 101, labels = "Future", cex = 0.7, xpd = TRUE, adj = c(0.5, 0))
}

plot_growth_rate_us <- function(gr_tab, has_se, fig_gen_run, filename) {
  if (fig_gen_run) {
    png(filename = filename, width = 8, height = 8, units = "in", pointsize = 16, res = 1000)
    on.exit(dev.off(), add = TRUE)
  }

  orpar <- par()
  on.exit(par(orpar), add = TRUE)
  par(mar = c(5.1, 4.1, 4.1, 4.1))

  gtp <- subset(gr_tab, variant_share >= 0.01)
  wow_x_scale <- floor(log(min(gtp$variant_share), 10))
  wow_x_min <- 5 * (10 ^ (wow_x_scale - 1))

  plot(
    x = gtp$variant_share,
    y = gtp$growth_rate,
    log = "x",
    type = "n",
    ylim = range(gtp$growth_rate_lo, gtp$growth_rate_hi),
    xaxt = "n",
    xlim = c(wow_x_min, 110),
    xlab = "Nowcast Estimated Proportion (%)",
    ylab = "Week over week growth rate (%)",
    main = "Nationwide"
  )

  if (!has_se) {
    mtext(
      text = "*No SE estimates b/c of non-invertible Hessian in multinomial model fit.",
      side = 3,
      line = 0,
      cex = 0.75,
      font = 4,
      col = "red"
    )
  }

  axis(side = 1, at = 10 ^ seq(wow_x_scale, 2, by = 1), labels = 10 ^ seq(wow_x_scale, 2, by = 1))
  abline(h = 0, col = "grey65")

  for (vv in seq_len(nrow(gtp))) {
    lines(x = rep(gtp$variant_share[vv], 2), y = gtp[vv, c("growth_rate_lo", "growth_rate_hi")], col = "blue", lwd = 2)
    lines(x = pmax(0.0001, gtp[vv, c("variant_share_lo", "variant_share_hi")]), y = rep(gtp$growth_rate[vv], 2), col = "blue", lwd = 2)
  }

  text(x = gtp$variant_share, y = gtp$growth_rate, labels = gtp$variant, cex = 0.85, col = "grey25", adj = 1.15)
  add_doubling_time_axis()
}

myciprop = function(voc,
                    geoid,
                    svy,
                    str   = TRUE,
                    range = FALSE,
                    mut   = FALSE,
                    level = 0.95,
                    ...) {
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

   
    warning_message = 'Hessian is non-invertible. Nowcast estimates will not have confidence intervals. Check for geographic regions with very few samples of a variant.'
    warning(warning_message)

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
