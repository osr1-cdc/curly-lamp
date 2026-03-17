#!/usr/bin/env Rscript
# Consolidated workflow entrypoint for the SC2 proportion modeling pipeline.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "")) y else x
}

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(script_path) == 0) "." else dirname(normalizePath(script_path))
})

parse_args <- function(args) {
  parsed <- list(
    command = NULL,
    user = NULL,
    password = NULL,
    dry_run = FALSE
  )

  if (length(args) == 0) {
    return(parsed)
  }

  parsed$command <- args[[1]]
  i <- 2

  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("--dry-run")) {
      parsed$dry_run <- TRUE
      i <- i + 1
      next
    }

    if (arg %in% c("-u", "--user")) {
      if (i == length(args)) stop("Missing value for ", arg, call. = FALSE)
      parsed$user <- args[[i + 1]]
      i <- i + 2
      next
    }

    if (arg %in% c("-p", "--password")) {
      if (i == length(args)) stop("Missing value for ", arg, call. = FALSE)
      parsed$password <- args[[i + 1]]
      i <- i + 2
      next
    }

    stop("Unknown argument: ", arg, call. = FALSE)
  }

  parsed
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript pipeline.R prepare-data -u <username> -p <password>",
      "  Rscript pipeline.R run1",
      "  Rscript pipeline.R run2",
      "  Rscript pipeline.R run-all -u <username> -p <password>",
      "  Rscript pipeline.R submit-all -u <username> -p <password>",
      "  Rscript pipeline.R <command> [--dry-run]",
      sep = "\n"
    ),
    "\n"
  )
}

require_credentials <- function(parsed) {
  if (is.null(parsed$user) || is.null(parsed$password)) {
    stop("Both --user and --password are required for this command.", call. = FALSE)
  }
}

run_command <- function(command, args, dry_run = FALSE) {
  pretty <- paste(shQuote(command), paste(shQuote(args), collapse = " "))
  cat(pretty, "\n")

  if (dry_run) {
    return(invisible(0))
  }

  status <- system2(command = command, args = args)
  if (!identical(status, 0L)) {
    stop("Command failed with status ", status, ": ", pretty, call. = FALSE)
  }

  invisible(status)
}

run_command_with_env <- function(command, args, env = character(), dry_run = FALSE) {
  pretty <- paste(shQuote(command), paste(shQuote(args), collapse = " "))
  cat(pretty, "\n")

  if (dry_run) {
    if (length(env) > 0) {
      cat("env:", paste(env, collapse = " "), "\n")
    }
    return(invisible(0))
  }

  status <- system2(command = command, args = args, env = env)
  if (!identical(status, 0L)) {
    stop("Command failed with status ", status, ": ", pretty, call. = FALSE)
  }

  invisible(status)
}

with_envvars <- function(env, code) {
  stopifnot(is.character(env))

  if (length(env) == 0) {
    force(code)
    return(invisible(NULL))
  }

  split_env <- strsplit(env, "=", fixed = TRUE)
  keys <- vapply(split_env, `[[`, character(1), 1)
  values <- vapply(
    split_env,
    function(parts) paste(parts[-1], collapse = "="),
    character(1)
  )

  old_values <- Sys.getenv(keys, unset = NA_character_)
  on.exit({
    for (i in seq_along(keys)) {
      if (is.na(old_values[[i]])) {
        Sys.unsetenv(keys[[i]])
      } else {
        do.call(Sys.setenv, setNames(list(old_values[[i]]), keys[[i]]))
      }
    }
  }, add = TRUE)

  do.call(Sys.setenv, stats::setNames(as.list(values), keys))
  force(code)
  invisible(NULL)
}

run_script_in_process <- function(script_path, main_function, env = character(), dry_run = FALSE) {
  pretty <- paste(shQuote("sys.source"), shQuote(script_path))
  cat(pretty, "\n")

  if (dry_run) {
    if (length(env) > 0) {
      cat("env:", paste(env, collapse = " "), "\n")
    }
    return(invisible(0))
  }

  script_env <- new.env(parent = globalenv())
  with_envvars(env, {
    sys.source(script_path, envir = script_env)
    script_env[[main_function]]()
  })

  invisible(0)
}

log_path <- file.path(script_dir, "log")

log_message <- function(message, dry_run = FALSE) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message)
  cat(line, "\n")

  if (!dry_run) {
    cat(line, "\n", file = log_path, append = TRUE)
  }
}

pipeline_env <- function(extra = character()) {
  c(
    paste0("SC2_PIPELINE_ROOT=", script_dir),
    paste0("SC2_CONFIG_PATH=", file.path(script_dir, "config", "config.R")),
    paste0("SC2_FUNCTIONS_PATH=", file.path(script_dir, "weekly_variant_report_functions.R")),
    paste0("SC2_DATA_DIR=", file.path(script_dir, "data")),
    paste0("SC2_RESULTS_DIR=", file.path(script_dir, "results")),
    paste0("SC2_RESOURCES_DIR=", file.path(script_dir, "resources")),
    paste0("SC2_CUSTOM_LINEAGES=FALSE"),
    paste0("SC2_CUSTOM_TAG="),
    extra
  )
}

run_prepare_data <- function(parsed) {
  require_credentials(parsed)

  run_script_in_process(
    script_path = file.path(script_dir, "variant_surveillance_system.R"),
    main_function = "variant_surveillance_system_main",
    env = pipeline_env(c(
      paste0("SC2_CDP_USER=", parsed$user),
      paste0("SC2_CDP_PASSWORD=", parsed$password),
      "SC2_NEXTCLADE_PANGO=FALSE"
    )),
    dry_run = parsed$dry_run
  )
}

run_model <- function(run_number, dry_run = FALSE) {
  run_script_in_process(
    script_path = file.path(script_dir, "weekly_variant_report_nowcast.R"),
    main_function = "weekly_variant_report_nowcast_main",
    env = pipeline_env(c(
      paste0("SC2_RUN_NUMBER=", run_number),
      "SC2_REDUCED_VOCS=F",
      "SC2_TRIM_WEIGHTS=quantile_99",
      "SC2_SAVE_DATASETS_TO_FILE=T",
      "SC2_PARALLEL_CORES=4",
      "SC2_WEIGHTED_METHODS=unweighted",
      "SC2_WEIGHT_TYPE=population",
      "SC2_VOC2_EXTRA_PREAGGREGATION=FALSE",
      "SC2_VOC_AGGREGATION_METHOD=updated"
    )),
    dry_run = dry_run
  )
}

run_all <- function(parsed) {
  require_credentials(parsed)
  run_prepare_data(parsed)
  run_model(run_number = 1, dry_run = parsed$dry_run)
  run_model(run_number = 2, dry_run = parsed$dry_run)
}

submit_qsub <- function(args, dry_run = FALSE) {
  run_command(command = "qsub", args = args, dry_run = dry_run)
}

submit_pipeline_job <- function(job_name, hold_jid = NULL, sync = FALSE, command_args, dry_run = FALSE) {
  qsub_args <- c("-N", job_name)

  if (!is.null(hold_jid)) {
    qsub_args <- c(qsub_args, "-hold_jid", hold_jid)
  }

  if (isTRUE(sync)) {
    qsub_args <- c(qsub_args, "-sync", "y")
  }

  qsub_args <- c(
    qsub_args,
    file.path(script_dir, "pipeline_job.sh"),
    command_args
  )

  submit_qsub(qsub_args, dry_run = dry_run)
}

submit_all <- function(parsed) {
  require_credentials(parsed)

  log_message("Modeling run started.", dry_run = parsed$dry_run)

  prep_status <- tryCatch(
    {
      submit_pipeline_job(
        job_name = "variant_surveillance",
        sync = TRUE,
        command_args = c(
          "prepare-data",
          "-u", parsed$user,
          "-p", parsed$password
        ),
        dry_run = parsed$dry_run
      )
      0L
    },
    error = function(e) {
      log_message("Data pull failed. Check Run_var_sys.err", dry_run = parsed$dry_run)
      stop(e$message, call. = FALSE)
    }
  )

  if (!identical(prep_status, 0L)) {
    log_message("Data pull failed. Check Run_var_sys.err", dry_run = parsed$dry_run)
    stop("Variant surveillance submission failed.", call. = FALSE)
  }

  log_message("Data has been pulled in.", dry_run = parsed$dry_run)

  submit_pipeline_job(
    job_name = "run1_CDT",
    hold_jid = "variant_surveillance",
    command_args = c("run1"),
    dry_run = parsed$dry_run
  )

  run2_ok <- tryCatch(
    {
      submit_pipeline_job(
        job_name = "run2_CDT",
        hold_jid = "variant_surveillance",
        sync = TRUE,
        command_args = c("run2"),
        dry_run = parsed$dry_run
      )
      TRUE
    },
    error = function(e) {
      log_message(
        "weekly_variant_report_nowcast.R errored out. Check the .err files for Run 1 and Run 2.",
        dry_run = parsed$dry_run
      )
      stop(e$message, call. = FALSE)
    }
  )

  if (isTRUE(run2_ok)) {
    log_message("Modeling run finished.", dry_run = parsed$dry_run)
  }
}

main <- function() {
  parsed <- parse_args(commandArgs(trailingOnly = TRUE))
  command <- parsed$command %||% ""

  if (command %in% c("", "-h", "--help", "help")) {
    usage()
    return(invisible(0))
  }

  switch(
    command,
    "prepare-data" = run_prepare_data(parsed),
    "run1" = run_model(run_number = 1, dry_run = parsed$dry_run),
    "run2" = run_model(run_number = 2, dry_run = parsed$dry_run),
    "run-all" = run_all(parsed),
    "submit-all" = submit_all(parsed),
    stop("Unknown command: ", command, call. = FALSE)
  )
}

main()
