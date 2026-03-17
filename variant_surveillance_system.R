load_variant_surveillance_packages <- function() {
  library(implyr)
  library(odbc)
  library(dplyr)
  library(RJDBC)
  library(optparse)
  library(data.table)
}

configure_variant_surveillance_options <- function() {
  options(
    java.parameters = "-Xmx8000m",
    stringsAsFactors = FALSE
  )
}

get_variant_pipeline_context <- function() {
  list(
    pipeline_root = Sys.getenv("SC2_PIPELINE_ROOT", unset = ""),
    pipeline_mode = nzchar(Sys.getenv("SC2_PIPELINE_ROOT", unset = "")),
    config_path = Sys.getenv("SC2_CONFIG_PATH", unset = ""),
    data_dir = Sys.getenv("SC2_DATA_DIR", unset = ""),
    resources_dir = Sys.getenv("SC2_RESOURCES_DIR", unset = "")
  )
}

parse_variant_surveillance_options <- function(pipeline_mode) {
  if (pipeline_mode) {
    return(list(
      user = Sys.getenv("SC2_CDP_USER", unset = ""),
      password = Sys.getenv("SC2_CDP_PASSWORD", unset = ""),
      nextclade_pango = Sys.getenv("SC2_NEXTCLADE_PANGO", unset = "FALSE")
    ))
  }

  option_list <- list(
    optparse::make_option(opt_str = c("-u", "--user"), type = "character", default = "", help = "User Name", metavar = "character"),
    optparse::make_option(opt_str = c("-p", "--password"), type = "character", default = "", help = "Password", metavar = "character"),
    optparse::make_option(opt_str = c("-n", "--nextclade_pango"), type = "character", default = "FALSE", help = "Whether to use nextclade_pango instead of pangolin lineage", metavar = "character")
  )

  optparse::parse_args(optparse::OptionParser(option_list = option_list))
}

resolve_variant_script_basename <- function(pipeline_mode, pipeline_root) {
  if (pipeline_mode) {
    return(pipeline_root)
  }

  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(pattern = file.arg.name, replacement = "", x = initial.options[grep(pattern = file.arg.name, x = initial.options)])
  script.basename <- dirname(script.name)
  if (length(script.basename) == 0) {
    script.basename <- "."
  }
  script.basename
}

resolve_variant_paths <- function(script.basename, context) {
  list(
    config_path = if (nzchar(context$config_path)) context$config_path else file.path(script.basename, "config", "config.R"),
    data_dir = if (nzchar(context$data_dir)) context$data_dir else file.path(script.basename, "data"),
    resources_dir = if (nzchar(context$resources_dir)) context$resources_dir else file.path(script.basename, "resources")
  )
}

get_variant_custom_settings <- function() {
  list(
    custom_lineages = toupper(Sys.getenv("SC2_CUSTOM_LINEAGES", unset = "FALSE")) %in% c("T", "TRUE", "Y", "YES"),
    custom_tag = Sys.getenv("SC2_CUSTOM_TAG", unset = "")
  )
}

initialize_variant_data_dir <- function(data_dir) {
  dir.create(path = data_dir, showWarnings = FALSE)
}

backup_dir_path <- function(data_dir, backup_date, custom_tag) {
  file.path(data_dir, paste0("backup_", backup_date, custom_tag))
}

backup_file_path <- function(data_dir, backup_date, custom_tag, filename) {
  file.path(backup_dir_path(data_dir, backup_date, custom_tag), filename)
}

backup_rds_path <- function(data_dir, backup_date, object_date, object_name, custom_tag) {
  file.path(
    backup_dir_path(data_dir, backup_date, custom_tag),
    paste0(object_date, "_", object_name, custom_tag, ".RDS")
  )
}

backup_files_exist <- function(data_dir, backup_date, object_dates, object_names, custom_tag) {
  all(vapply(
    object_names,
    FUN.VALUE = logical(1),
    FUN = function(object_name) {
      file.exists(backup_rds_path(
        data_dir = data_dir,
        backup_date = backup_date,
        object_date = object_dates[[object_name]],
        object_name = object_name,
        custom_tag = custom_tag
      ))
    }
  ))
}

read_backup_rds <- function(data_dir, backup_date, object_date, object_name, custom_tag) {
  readRDS(file = backup_rds_path(
    data_dir = data_dir,
    backup_date = backup_date,
    object_date = object_date,
    object_name = object_name,
    custom_tag = custom_tag
  ))
}

save_backup_rds <- function(object, data_dir, backup_date, object_date, object_name, custom_tag) {
  saveRDS(
    object = object,
    file = backup_rds_path(
      data_dir = data_dir,
      backup_date = backup_date,
      object_date = object_date,
      object_name = object_name,
      custom_tag = custom_tag
    )
  )
}

write_backup_csv <- function(x, data_dir, backup_date, custom_tag, filename, row.names = FALSE) {
  write.csv(
    x = x,
    file = backup_file_path(data_dir, backup_date, custom_tag, filename),
    row.names = row.names
  )
}

save_backup_rds_file <- function(object, data_dir, backup_date, custom_tag, filename) {
  saveRDS(
    object = object,
    file = backup_file_path(data_dir, backup_date, custom_tag, filename)
  )
}

find_previous_backup_rds <- function(data_dir, current_date, custom_tag, suffix, previous_dates = 6:24) {
  for (offset in previous_dates) {
    previous_date <- current_date - offset
    candidate <- backup_file_path(
      data_dir = data_dir,
      backup_date = previous_date,
      custom_tag = custom_tag,
      filename = paste0(previous_date, suffix, custom_tag, ".RDS")
    )

    if (file.exists(candidate)) {
      return(candidate)
    }
  }

  NA_character_
}

regex_lab_selector <- function(pattern) {
  function(labs) {
    grep(pattern = pattern, x = labs, ignore.case = TRUE, value = TRUE)
  }
}

rename_labs_by_selector <- function(svy_dat, unique_labs, selector, new_name) {
  matched_labs <- selector(unique_labs)
  if (!length(matched_labs)) {
    return(list(
      svy_dat = svy_dat,
      labnames_df = NULL
    ))
  }

  svy_dat[LAB %in% matched_labs, LAB2 := new_name]

  list(
    svy_dat = svy_dat,
    labnames_df = data.frame(
      old_name = matched_labs,
      new_name = new_name,
      stringsAsFactors = FALSE
    )
  )
}

apply_lab_rename_rules <- function(svy_dat, unique_labs, rules) {
  rename_logs <- vector(mode = "list", length = length(rules))

  for (ii in seq_along(rules)) {
    rename_result <- rename_labs_by_selector(
      svy_dat = svy_dat,
      unique_labs = unique_labs,
      selector = rules[[ii]]$selector,
      new_name = rules[[ii]]$new_name
    )
    svy_dat <- rename_result$svy_dat
    rename_logs[[ii]] <- rename_result$labnames_df
  }

  rename_logs <- Filter(Negate(is.null), rename_logs)
  if (!length(rename_logs)) {
    rename_logs <- list(data.frame(
      old_name = character(),
      new_name = character(),
      stringsAsFactors = FALSE
    ))
  }

  list(
    svy_dat = svy_dat,
    labnames_df = data.table::rbindlist(rename_logs, use.names = TRUE)
  )
}

variant_surveillance_system_main <- function() {

# Setup ------------------------------------------------------------------------
load_variant_surveillance_packages()
configure_variant_surveillance_options()

context <- get_variant_pipeline_context()
opts <- parse_variant_surveillance_options(context$pipeline_mode)

# If running code interactively, fill in opts values here
# opts$user = ""
# opts$password = ""


custom_settings <- get_variant_custom_settings()
custom_lineages <- custom_settings$custom_lineages
custom_tag <- custom_settings$custom_tag

## get base directory and driver location
script.basename <- resolve_variant_script_basename(context$pipeline_mode, context$pipeline_root)

# create data directory
paths <- resolve_variant_paths(script.basename, context)
config_path <- paths$config_path
data_dir <- paths$data_dir
resources_dir <- paths$resources_dir
initialize_variant_data_dir(data_dir)

# load in config/config.R variables
source(config_path)

# path to the jdbc driver --> now from the config.R file
# jdbc_driver = paste0(script.basename, "/jdbc/ClouderaImpalaJDBC-2.6.20.1024/ClouderaImpalaJDBC41-2.6.20.1024/")

# variable identifying which Hadoop table to read in. Options are:
# ~ "deduplication_cdcncbigisaid_auto": this data table is updated regularly but may be subject to cleaning issues
# ~  "analytics_metadata_frozen": this data table is updated less frequently, but is the cleanest version of the sequence data
# NOTE - if running the official Friday analysis use the "analytics_metadata_frozen" table
seq_table = "sc2_archive.analytics_metadata_frozen"
# CDP cluster node to use
# (you can use any node from 08 - 13. Sometimes nodes fail so if you get an error you can try another node)
node = "10"

# Set pangolineage definition source based on nextclade_pango flag
if(toupper(opts$nextclade_pango) %in% c('T', 'TRUE', 'Y', 'YES')){
  lineage_table = "sc2_src.nextclade"
  lineage_field = "nextclade_pango"
  expanded_lineage_field = "" # no expanded lineage for nextclade
  current_data  = TRUE   # If using nextclade_pango, there is no archived frozen data, only current data from nextclade_pango can be used
} else {
  lineage_table = "sc2_src.pangolin"
  lineage_field = "lineage"
  expanded_lineage_field = ", expanded_lineage" # include the comma for the query
}

#ref_lineage = opts$reference_lineage

# Import Data ------------------------------------------------------------------

# If the data was already pulled and you want to just use that data instead of re-pulling it, set here.
# This is useful if you aggregate some lab names at the end of this code and then want to re-run the
# script after changing which labs get aggregated.

# use previously pulled data if it exists
same_day_backup_dates <- c(
  data = data_date,
  pangolin = data_date,
  pops = data_date,
  tests_nrevss = data_date
)

frozen_backup_dates <- c(
  data = date_frozen_toread,
  pangolin = date_frozen_toread,
  tests = data_date,
  tests_nrevss = data_date,
  pops = data_date
)

if(use_previously_imported_data & date_frozen_toread == data_date &
    backup_files_exist(
      data_dir = data_dir,
      backup_date = data_date,
      object_dates = same_day_backup_dates,
      object_names = c("data", "pangolin", "pops", "tests_nrevss"),
      custom_tag = custom_tag
    )){

  print('Reading in previously pulled data')
  dat <- read_backup_rds(data_dir, data_date, data_date, "data", custom_tag)
  dat <- subset(dat, !(dat$covv_accession_id %in% c("EPI_ISL_19791260",
                                                   "EPI_ISL_19791262",
                                                   "EPI_ISL_19791263",
                                                   "EPI_ISL_19791264",
                                                   "EPI_ISL_19791265",
                                                   "EPI_ISL_19791266",
                                                   "EPI_ISL_19791267",
                                                   "EPI_ISL_19791268")))
  dat <- as.data.table(dat)                                                 

  pangolin <- read_backup_rds(data_dir, data_date, data_date, "pangolin", custom_tag)
  tests_nrevss <- read_backup_rds(data_dir, data_date, data_date, "tests_nrevss", custom_tag)
  pops <- read_backup_rds(data_dir, data_date, data_date, "pops", custom_tag)
  print('Finished reading in data.')
} else if (use_previously_imported_data & date_frozen_toread != data_date &
    backup_files_exist(
      data_dir = data_dir,
      backup_date = data_date,
      object_dates = frozen_backup_dates,
      object_names = c("data", "pangolin", "tests", "pops", "tests_nrevss"),
      custom_tag = custom_tag
    )
    ){
  print('Reading in previously pulled data. analytics_metadata frozen date different from test date')
  dat <- read_backup_rds(data_dir, data_date, date_frozen_toread, "data", custom_tag)
  pangolin <- read_backup_rds(data_dir, data_date, date_frozen_toread, "pangolin", custom_tag)
  tests <- read_backup_rds(data_dir, data_date, data_date, "tests", custom_tag)
  tests_nrevss <- read_backup_rds(data_dir, data_date, data_date, "tests_nrevss", custom_tag)
  pops <- read_backup_rds(data_dir, data_date, data_date, "pops", custom_tag)
  print('Finished reading in data.')

} else {

print('Pulling data from CDP')
## Get the sequence data
impala_classpath <- list.files(path = jdbc_driver,
                               pattern = "\\.jar$",
                               full.names = TRUE)

# Initialize the Java virtual machine
# (this function utilizes the "java.parameters" set in "options" above)
rJava::.jinit(classpath = impala_classpath)

# specify the driver for connecting to CDP database
drv <- RJDBC::JDBC(
  driverClass = "com.cloudera.impala.jdbc41.Driver",
  classPath = impala_classpath,
  identifier.quote = "`"
)

# connect to the CDP database
impala <- DBI::dbConnect(
  drv  = drv,
  url = paste0(
    "jdbc:impala://cdp-",
    node,
    ".biotech.cdc.gov:21050/sc2_air;AuthMech=3;useSasl=1;SSL=1;AllowSelfSignedCerts=1;CAIssuedCertNamesMismatch=1"),
  user     = opts$user,
  password = opts$password
)
# NOTE: The above error: "ERROR StatusLogger No Log4j 2 configuration file found." is expected. You can ignore it.

# Validate that data_date exists in the frozen source tables used by the pipeline.
valid_data_dates <- DBI::dbGetQuery(
  conn = impala,
  statement = '
	SELECT DISTINCT to_date(date_frozen)
	FROM sc2_archive.analytics_metadata_frozen
	    ')
valid_tests_dates <- DBI::dbGetQuery(
  conn = impala,
  statement = '
	SELECT DISTINCT to_date(date_frozen)
FROM sc2_archive.analytics_metadata_frozen
    ')

#throw an error if the data_date isn't in the
#sc2_archive.analytics_metadata_frozen table.
if( !((as.character(data_date) %in% unlist(valid_data_dates)) &
      (as.character(data_date) %in% unlist(valid_tests_dates))) ){
  errorCondition(paste0(
    'The "data_date" provided (',
    data_date,
    ') is not in valid.\n
    Valid options from sc2_archive.analytics_metadata_frozen.date_frozen include:\n',
    paste(sort(valid_data_dates[,1]), collapse = '\n'),
    '\nValid options from sc2_archive.hhs_protect_testing_frozen.date_frozen include:\n',
    paste(sort(valid_tests_dates[,1]), collapse = '\n'),
    '.'
  ))
}


# Get all field/column names from the table
all.vars = DBI::dbGetQuery(impala, "DESCRIBE sc2_archive.analytics_metadata_frozen")[, 1]

# specify the fields/columns that we want to get
if(seq_table == "sc2_archive.analytics_metadata_frozen"){
  # get variables that:
  #   - start with "csid", "primary", "covv"
  #   - contain    "targeted", "vendor"
  #   - equal to   "spike_mutations", "eventid_all"
  get.vars = grep(pattern = "(^csid)|(^primary)|(^covv)|(targeted)|(vendor)|(^spike_mutations$)|(^eventid_all$)|(^received_date)",
                  x = all.vars,
                  value = TRUE)
} else {
  # specify the variable to get manually
  get.vars= c(
    "primary_virus_name",
    "primary_nt_id",
    "primary_country",
    "primary_state",
    "primary_state_abv",
    "primary_host",
    "primary_specimen",
    "primary_collection_date_dt",
    "csid",
    #"contractor_vendor_id",
    "contractor_vendor_name",
    "contractor_targeted_sequencing",
    "covv_accession_id",
    "covv_orig_lab",
    "covv_add_host_info",
    "primary_collection_date",
    "primary_submitter",
    "spike_mutations",
    "lineage",
    "received_date"
  )#,"patient_age_cl","covv_patient_age")
}


# specify the database query to run
# (Subset to US sequences for faster download)
# add in geni include or exclude info and S1_key for S1_id based proportion analysis - 2022-10-13
query = paste(
  "SELECT
    CASE
      WHEN (( NOT A.spike_aln LIKE '%XX%' AND NOT A.spike_aln LIKE '%..%' ) AND NOT regexp_like( A.spike_aln, '~' ) AND LENGTH( A.spike_aln ) = 1274) THEN 'TRUE'
      ELSE 'FALSE'
      END as geni_included,
    SUBSTRING(A.spike_aln, 14, 677) as s1_key,",
  paste0('A.', get.vars, collapse=', '),
  #', coalesce(A.contractor_vendor_name, A.primary_submitter) source',
  ', COALESCE(IF(eventid_all LIKE "%1771%", "NS3", NULL), contractor_vendor_name, primary_submitter) source',
  ', lineage as pangolineage',
  paste0(
    ' FROM sc2_archive.analytics_metadata_frozen as A
    INNER JOIN
    (SELECT max(date_frozen) as max_frozen
    FROM sc2_archive.analytics_metadata_frozen ma
    WHERE to_date(ma.date_frozen) = ', date_frozen, ') as M
    ON A.date_frozen = M.max_frozen'
  ),
  ' WHERE A.primary_country in ("United States", "USA")
  AND A.primary_host = "Human"
  AND ( A.contractor_vendor_name IS NOT NULL OR A.eventid_all LIKE "%1771%" OR A.primary_sampling_strategy = "Baseline_Surveillance" )
  AND A.lineage != "Unassigned" AND A.lineage IS NOT NULL'
  # AND A.primary_submitter IS NOT NULL'
) # if unavailable, use test_deduplication_cdcncbigisaid_auto for testing

# pull the data from the CDP database to an R data.frame
dat = DBI::dbGetQuery(conn = impala,
                      statement = query)
dat <- subset(dat, !(dat$covv_accession_id %in% c("EPI_ISL_19791260",
                                                   "EPI_ISL_19791262",
                                                   "EPI_ISL_19791263",
                                                   "EPI_ISL_19791264",
                                                   "EPI_ISL_19791265",
                                                   "EPI_ISL_19791266",
                                                   "EPI_ISL_19791267",
                                                   "EPI_ISL_19791268")))
dat <- as.data.table(dat)

# Get the Pango lineages (standard Pangolin lineages only)
if(current_data){
  pangolin = DBI::dbGetQuery(
    conn = impala,
    statement = paste0(
    'SELECT DISTINCT nt_id, lineage, expanded_lineage ',
    ' FROM ', lineage_table
    ))
} else {
  pangolin = DBI::dbGetQuery(
    conn = impala,
    statement = paste0(
      '
      SELECT DISTINCT
      P.primary_nt_id as nt_id,
      P.lineage,
      P.expanded_lineage
      FROM sc2_archive.analytics_metadata_frozen as P
      WHERE to_date(P.date_frozen) = ', date_frozen
    ))
}

# download S gene mutation lists, and source (for NS3 and labs)


# download the state testing data
tests = DBI::dbGetQuery(
conn = impala,
statement = paste0(
    'SELECT to_date(H.collection_date) as collection_date,
H.reporting_state,
H.INDETERMINATE,
H.INVALID,
H.NEGATIVE,
H.POSITIVE
FROM sc2_archive.hhs_protect_testing_frozen as H
INNER JOIN
(SELECT max(date_frozen) as max_frozen
    FROM sc2_archive.hhs_protect_testing_frozen hf
    WHERE to_date(hf.date_frozen) = ', date_frozen, '
) as F
ON H.date_frozen = F.max_frozen
WHERE H.collection_date is NOT NULL'
))

# rename the columns of the testing data
colnames(tests) = c("collection_date",
                    "STUSAB", # "reporting_state",
                    "INDETERMINATE",
                    "INVALID",
                    "NEGATIVE",
                    "POSITIVE")

# add code for testing data portion exclusion
if(exclude_testing_data_portion) {

# subset unaffected states
tests_valid <- tests %>%
    filter(STUSAB %notin% exclusion_states)

# filter out exclusion state data past a given date
tests_exclusion <- tests %>%
    filter(STUSAB %in% exclusion_states &
    (as.Date(collection_date) < as.Date(testing_exclusion_cutoff) |
    as.Date(collection_date) > as.Date(testing_exclusion_cutoff_end)))

tests <- bind_rows(tests_valid, tests_exclusion)
}

# download NREVSS testing data
tests_nrevss = DBI::dbGetQuery(
  conn = impala,
  statement = paste0( # 'SELECT * FROM sc2_archive.nrevss_frozen'
'SELECT
  H.*,
  to_date(H.mmwrweek_end) as collection_week
FROM sc2_archive.nrevss_frozen H
INNER JOIN
(SELECT max(date_frozen) as max_frozen
    FROM sc2_archive.nrevss_frozen hf
    WHERE to_date(hf.date_frozen) = ', date_frozen, '  OR to_date(date_add(hf.date_frozen,-1)) = ', date_frozen, ' OR to_date(date_add(hf.date_frozen,1)) = ', date_frozen, '
) as F
ON H.date_frozen = F.max_frozen'
)) #Added condition to accomandate situation that nrevss data is frozen on a Wednesday, so that would be one day before or after the frozen date for all other data

# Get the vocs included in run 2
# only if voc2_manual is not set
if(is.na(voc2_manual)){
# SQL code from: https://cdc.sharepoint.com/teams/NCEZID-OD_CAWG/Shared%20Documents/Forms/AllItems.aspx?FolderCTID=0x01200085513C439903124B82D8FF6016BB819B&id=%2Fteams%2FNCEZID%2DOD%5FCAWG%2FShared%20Documents%2FDAV%2DActivity%2FSOPs%2Fcurrent%5Flineages%5F1%5Fpercent%2Esql&parent=%2Fteams%2FNCEZID%2DOD%5FCAWG%2FShared%20Documents%2FDAV%2DActivity%2FSOPs

# This is a heavily modified version of "current_lineages_1_percent.sql", that
# groups by regions & weeks (instead of just weeks). Resulting lineages include
# any lineage that's above 1% in any region-week combination.

# this query was formerly defined in "current_lineages_1_percent.sql",
# but has been updated to accept "data_date" as an argument.
# "data_date" must be a date on which archive data was created.

# for the moment just use the old query that does not group by HHS region:
voc2_df = DBI::dbGetQuery(
    conn = impala,
    statement = paste0(
"SELECT QQ.*,
    cor.*
FROM
(SELECT Q.lineage,
        to_date(max(biweek_ending)) AS most_recent,
        max(fraction) AS max_fraction,
        max(lineage_count) AS max_virus_count,
        to_date('", data_date, "') AS pull_date
FROM
    (SELECT l.", lineage_field, " as lineage,
            c.variant_type,
            date_add(primary_collection_date,datediff('2039-12-31',primary_collection_date)%14) AS biweek_ending,
            count(a.primary_virus_name) AS lineage_count,
            z.region_total,
            count(a.primary_virus_name) / z.region_total AS fraction,
            if(count(a.primary_virus_name) / z.region_total >= 0.01, TRUE, FALSE) AS is_one_percent,
            if(count(a.primary_virus_name) / z.region_total >= 0.005, TRUE, FALSE) AS is_zerofive_percent
    FROM sc2_archive.analytics_metadata_frozen a
    LEFT JOIN ", lineage_table, " l on a.primary_nt_id = l.nt_id
    LEFT JOIN
        (SELECT date_add(primary_collection_date, datediff('2039-12-31',primary_collection_date)%14) AS biweek_ending,
                count(za.primary_virus_name) AS region_total
        FROM sc2_archive.analytics_metadata_frozen za
        WHERE
            ( contractor_vendor_name IS NOT NULL OR eventid_all LIKE '%1771%' OR primary_sampling_strategy = 'Baseline_Surveillance' )
            AND to_date(date_frozen) = ", date_frozen, "
            AND primary_country = 'United States'
        GROUP BY biweek_ending) z ON date_add(primary_collection_date,datediff('2039-12-31',primary_collection_date)%14) = z.biweek_ending
    AND 1=1
    LEFT JOIN sc2_src.variant_definitions c ON a.lineage = c.lineage
    WHERE -- THis is generally the weeks. For the publishing week and the next off week, -9 is -1st 2-week period  -106 is the -7th 2-week period.
        biweek_ending <= date_add(date_trunc('week', date_add(to_date('", data_date, "'), 1)), -9) 
        AND biweek_ending >= date_add(date_trunc('week', date_add(to_date('", data_date, "'), 1)), -106)
        AND ( contractor_vendor_name IS NOT NULL OR eventid_all LIKE '%1771%' OR primary_sampling_strategy = 'Baseline_Surveillance' )
        AND primary_country = 'United States'
        AND to_date(date_frozen) = ", date_frozen, "
    GROUP BY l.", lineage_field, ",
            c.variant_type,
            biweek_ending,
            z.region_total
) Q
WHERE Q.is_one_percent IS TRUE --OR Q.variant_type is not null
    OR (Q.is_zerofive_percent IS TRUE
    -- For the publishing week and the next off week, this choose the -2nd 2-week period (last weighted period)
        AND Q.biweek_ending <= date_add(date_trunc('week', date_add(to_date('", data_date, "'), 1)), -23)
        AND Q.biweek_ending >= date_add(date_trunc('week', date_add(to_date('", data_date, "'), 1)), -30)
        AND lineage_count > 1)
GROUP BY lineage) QQ
LEFT JOIN sc2_air.analytics_lineage_corr cor ON QQ.lineage = cor.lineage
WHERE cor.date_range_of_calc LIKE '%US:3mo'
Order by max_fraction desc"
))

# get the variant names
voc2_auto = sort(voc2_df$lineage)

# save the results to file
saveRDS(object = voc2_auto,
        file = paste0(script.basename,
                        "/data/voc2_auto_", data_date, custom_tag, ".RDS"))
}

# end the database connection
DBI::dbDisconnect(conn = impala)

# read in static table of state populations
# (population data may need to be updated)
pops = read.delim(file = file.path(resources_dir, "ACStable_B01001_40_2018_5.txt"))[, c("STUSAB", "Total.")]

# create a folder for the backup data (data that do not get used in "weekly_variant_report_nowcast.R")
dir.create(
  path = backup_dir_path(data_dir, data_date, custom_tag),
  showWarnings = FALSE
)

save_backup_rds(dat, data_dir, data_date, data_date, "data", custom_tag)
save_backup_rds(pangolin, data_dir, data_date, data_date, "pangolin", custom_tag)
save_backup_rds(tests, data_dir, data_date, data_date, "tests", custom_tag)
save_backup_rds(tests_nrevss, data_dir, data_date, data_date, "tests_nrevss", custom_tag)
save_backup_rds(pops, data_dir, data_date, data_date, "pops", custom_tag)
print('Finished reading in data.')
}
# Data Cleaning ----------------------------------------------------------------

# Apply the core sequence and metadata cleaning rules described in the README.

# Keep track of why sequences are dropped from the analysis
{
   # get number of duplicate rows by week from dat
   dat.dt <- as.data.table(dat) # much faster as a data.table
   # add in a column for week
   dat.dt[, 'date' := as.Date(primary_collection_date)]
   dat.dt[is.na(date), 'date' := as.Date(primary_collection_date_dt)] # fill in missing values
   dat.dt[, 'week' := (date - as.numeric(strftime(date, format="%w")))]

   # get counts of original sequences by week
   cnt_by_wk <- dat.dt[, .('count' = .N), by = 'week']

   # get counts of duplicated entries by week
   dups_by_wk <- dat.dt[
      duplicated(dat.dt, by = names(dat.dt)), # filter only duplicated rows
      ][,
        .('count' = .N), # count number of rows by week
        by = 'week']

   # weeks in the analysis
   wks <- seq(from = week0day1,
              to = data_date - as.numeric(format(data_date, '%w')),
              by = 7)
   # weeks that are in the data
   datwks <- unique(dat.dt$week)
   # combine all the weeks
   allwks <- sort(unique(c(wks, datwks)))

   # long data.table to hold results
   dropped_sequences <- as.data.table(
      expand.grid(
         'week'   = allwks,
         'reason' = c('n_original',
                      #'n_dropped_duplicates',
                      #'n_dropped_not_human',
                      #'n_dropped_not_US',
                      'n_dropped_invalid_state',
                      'n_dropped_invalid_date',
                      'n_dropped_invalid_virus_name',
                      #'n_dropped_CDC_lab',
                      'n_dropped_labs_100_seqs',
                      'n_dropped_invalid_lab_name',
                      'n_dropped_invalid_variant_name',
                      'n_dropped_invalid_weight'
         )
      )
   )

   # merge in the total sequence counts
   dropped_sequences <- merge(
      x = dropped_sequences,
      y = cnt_by_wk[, .('week' = week, 'count' = count, 'reason' = 'n_original')],
      by = c('week', 'reason'),
      all.x = TRUE
   )
   # fill in 0's for rows that didn't have any sequences dropped
   dropped_sequences[is.na(count) & reason == 'n_original', 'count' := 0]
   # merge in the counts of dropped sequences
  #  dropped_sequences <- rbind(
  #     dropped_sequences[!(week %in% dups_by_wk$week & reason == 'n_dropped_duplicates')],
  #     dups_by_wk[, .('week' = week, 'reason' = 'n_dropped_duplicates', 'count' = count)]
  #  )
   # fill in 0's for rows that didn't have any sequences dropped
   #dropped_sequences[is.na(count) & reason == 'n_dropped_duplicates', 'count' := 0]
}

## exclude duplicates from each table
dat          = distinct(dat.dt) # use the data.table version to speed up calculation of why sequences are being dropped
dat <- subset(dat, !(dat$covv_accession_id %in% c("EPI_ISL_19791260",
                                                   "EPI_ISL_19791262",
                                                   "EPI_ISL_19791263",
                                                   "EPI_ISL_19791264",
                                                   "EPI_ISL_19791265",
                                                   "EPI_ISL_19791266",
                                                   "EPI_ISL_19791267",
                                                   "EPI_ISL_19791268")))
dat <- as.data.table(dat)

pangolin     = distinct(pangolin)
tests_nrevss = distinct(tests_nrevss)
pops         = distinct(pops)

## Some general parameters:
# start day is the first Sunday of 2020
#week0day1 = as.Date("2020-01-05") Get parameter from config/config.R

# HHS regions
HHS_reg = list(HHS1 = c("CT", "ME", "MA", "NH", "RI", "VT"),
               HHS2 = c("NJ", "NY", "PR", "VI"),
               HHS3 = c("DE", "DC", "MD", "PA", "VA", "WV"),
               HHS4 = c("AL", "FL", "GA", "KY", "MS", "NC", "SC", "TN"),
               HHS5 = c("IL", "IN", "MI", "MN", "OH", "WI"),
               HHS6 = c("AR", "LA", "NM", "OK", "TX"),
               HHS7 = c("IA", "KS", "MO", "NE"),
               HHS8 = c("CO", "MT", "ND", "SD", "UT", "WY"),
               HHS9 = c("AZ", "CA", "HI", "NV", "AS", "MP", "GU", "MH", "PW"),
               HHS10= c("AK", "ID", "OR", "WA")
)

# Add in populations of US territories
# (b/c they're not in the "ACStable_B01001_40_2018_5.txt" file)
# Territory 2020-07-01 populations from https://en.wikipedia.org/wiki/List_of_states_and_territories_of_the_United_States_by_population (2021-04-19)
# Marshall Islands 2018 estimate https://en.wikipedia.org/wiki/Marshall_Islands (2021-04-19)
# Palau 2018 estimate https://en.wikipedia.org/wiki/Palau#:~:text=%E2%80%A2-,2018%20estimate,18%2C024,-%5B3%5D%5B4 (2023-05-09)
pops = rbind(pops,
             data.frame(STUSAB = c("AS", "GU", "MP", "VI", "MH", "PW"),
                        `Total.`=c(49437, 168485, 51433, 106235, 58413, 18024)))
# create a lookup table with states & their HHS regions
hhs = data.frame(STUSAB = toupper(pops$STUSAB))
# get the HHS region for each state
hhs$HHS = sapply(X = hhs$STUSAB, FUN = grep, HHS_reg)

## Subset by time and place:
# - USA,
# - states with two letter abbreviations,
# - human host,
# - reasonable collection date

# Subset to exclude faulty state abbreviations
{
   # get counts of faulty state abbreviations by week
   valid_state_abb <- nchar(dat$primary_state_abv) == 2

   # counts of faulty state abbreviations by week
   fs_by_wk <- dat[ (!valid_state_abb) | is.na(valid_state_abb), .('count' = .N), by = 'week']
   # merge in the counts of faulty state abbreviations
   if(nrow(fs_by_wk) > 0)
      dropped_sequences <- rbind(
         dropped_sequences[!(week %in% fs_by_wk$week & reason == 'n_dropped_invalid_state')],
         fs_by_wk[, .('week' = week, 'reason' = 'n_dropped_invalid_state', 'count' = count)]
      )
   # fill in 0's for rows that didn't have any sequences dropped
   dropped_sequences[is.na(count) & reason == 'n_dropped_invalid_state', 'count' := 0]
}
us.dat = subset(x = dat,
                valid_state_abb)

## Disambiguate and remove unreasonable dates
# convert collection date to "date" format
us.dat$collection_date = as.Date(us.dat$primary_collection_date)
# fill in some missing dates
us.dat$collection_date = with(
  data = us.dat,
  expr = as.Date(ifelse(test = is.na(collection_date),
                        yes = as.Date(primary_collection_date_dt),
                        no = collection_date),
                 origin = "1970-01-01")
)

# Subset to exclude NULL primary_virus_name entry
{
   # counts of invalid virus name by week
   novirusname_by_wk <- us.dat[ is.na(us.dat$primary_virus_name), .('count' = .N), by = 'week']
   # merge in the counts of invalid virus name
   if(nrow(novirusname_by_wk) > 0)
      dropped_sequences <- rbind(
         dropped_sequences[!(week %in% novirusname_by_wk$week & reason == 'n_dropped_invalid_virus_name')],
         novirusname_by_wk[, .('week' = week, 'reason' = 'n_dropped_invalid_virus_name', 'count' = count)]
      )
   # fill in 0's for rows that didn't have any sequences dropped
   dropped_sequences[is.na(count) & reason == 'n_dropped_invalid_virus_name', 'count' := 0]
}
us.dat = subset(x = us.dat,
                !is.na(primary_virus_name))

# convert to date format
if ("covv_subm_date" %in% names(us.dat)) us.dat$covv_subm_date = as.Date(us.dat$covv_subm_date)
if ("contractor_receive_date_to_cdc" %in% names(us.dat)) {
  us.dat$contractor_receive_date_to_cdc = as.Date(us.dat$contractor_receive_date_to_cdc)
}

# Add a column for the date of the first day of the week [should be the same as "week" calculated above]
us.dat$yr_wk = as.character(us.dat$collection_date - as.numeric(strftime(us.dat$collection_date, format="%w")))
# Add a column for the # of days since the start date (first Sunday of 2020)
us.dat$DAY = as.numeric(us.dat$collection_date - week0day1)

# exclude unreasonably early dates
{
   valid_dates <- us.dat$collection_date >= as.Date("2019-10-01")
   # counts of invalid dates by week
   id_by_wk <- us.dat[ (!valid_dates) | is.na(valid_dates), .('count' = .N), by = 'week']
   # merge in the counts of faulty state abbreviations
   if(nrow(id_by_wk) > 0)
      dropped_sequences <- rbind(
         dropped_sequences[!(week %in% id_by_wk$week & reason == 'n_dropped_invalid_date')],
         id_by_wk[, .('week' = week, 'reason' = 'n_dropped_invalid_date', 'count' = count)]
      )
   # fill in 0's for rows that didn't have any sequences dropped
   dropped_sequences[is.na(count) & reason == 'n_dropped_invalid_date', 'count' := 0]
}
us.dat = subset(x = us.dat,
                valid_dates)

# merge Pangolin lineage into us.dat
us.dat = merge(x = us.dat,
              y = pangolin, # [, c("nt_id", "lineage", "expanded_lineage")],
              by.x = "primary_nt_id",
              by.y = "nt_id",
              all.x = TRUE)
# # fill in missing lineage data with "covv_lineage"
# # (this line of code only works if using the dedupe table)
#if ("covv_lineage" %in% names(us.dat)) us.dat$lineage = with(us.dat, ifelse(is.na(lineage), covv_lineage, lineage))

# print out counts of omicron sequences
# print('Table of all omicron sequences:')
# table(us.dat$lineage[grep(pattern = '(B\\.1\\.1\\.529)|(BA\\.[0-9])', x = us.dat$lineage)])

# Merge S-gene mutation (baseline) data into the us.dat
# NS3 identifier in baseline$source:
# us.dat = merge(x = us.dat,
#                y = baseline[, c("nt_id", "source", "primary_virus_name", "s_mut")],
#                by.x = c("primary_nt_id", "primary_virus_name"),
#                by.y = c("nt_id", "primary_virus_name"),
#                all.x = TRUE)

# Add in HHS regions
us.dat = merge(x = us.dat,
               y = hhs,
               by.x = "primary_state_abv",
               by.y = "STUSAB",
               all.x = TRUE)


## NREVSS testing data (used in place of CLERS data for test positivity only in updated weighting scheme)
tests_nrevss <- data.table::as.data.table(tests_nrevss)

# double-check that we're only using NREVSS testing aggregated to the state level (I don't think this will ever exclude anything)
tests_nrevss <- tests_nrevss[ level == 'State' & source == 'NREVSS']

# convert to date format
if(!('collection_week' %in% names(tests_nrevss))) tests_nrevss[, 'collection_week' := as.Date( as.POSIXct(mmwrweek_end, format = '%Y-%m-%d:%H:%M:%S') , format = '%Y-%m-%d')]
if( class(tests_nrevss$collection_week) == 'character' ) tests_nrevss[, 'collection_week' := as.Date(x = collection_week, format = "%Y-%m-%d")] # UPDATE format

# exclude unreasonable dates
tests_nrevss <- tests_nrevss[collection_week >= as.Date("2019-10-01") &
                               collection_week <= data_date,] # this will exclude the current week (which should be ok)

# Add a column for the date of the first day of the week
tests_nrevss[,'yr_wk' := as.character(collection_week - as.numeric(strftime(collection_week, format="%w")))]
tests_nrevss[, 'week' := as.numeric(as.Date(yr_wk) - week0day1)/7]
tests_nrevss[, 'fortnight_end' := as.character(week0day1 + week %/% 2 * 14 + 13)]

# rename columns to match "tests"
tests_nrevss[, 'TOTAL'    := as.integer(tests)]
tests_nrevss[, 'POSITIVE' := as.integer(detections)]
tests_nrevss[, 'STUSAB'   := state]

# Replace NAs with 0
tests_nrevss[is.na(POSITIVE), 'POSITIVE' := 0]

# Aggregate tests by state & week for weighting
#   this shouldn't ever be necessary b/c the data are already aggregated by state-week
#   there can be missing state-week combinations in tests_nrevss, but there shouldn't ever be multiple rows per state-week.
#     SELECT DISTINCT A.week_count FROM (SELECT count(mmwrweek_end) as week_count FROM sc2_src.nrevss GROUP BY mmwrweek_end, state) as A;
tests_nrevss_wk <- tests_nrevss[,
                  .(
                    'POSITIVE' = unique(na.omit(POSITIVE)), # I'll assume that if there are multiple rows, that the data were duplicated.
                    'TOTAL'    = unique(na.omit(TOTAL))
                  ),
                  by = c('STUSAB', 'yr_wk')]
# aggregate tests_nrevss by state & fortnight (currently no plan to use this)
tests_nrevss_fn <- tests_nrevss[,
                                .(
                                  'POSITIVE' = unique(na.omit(POSITIVE)), # I'll assume that if there are multiple rows, that the data were duplicated.
                                  'TOTAL'    = unique(na.omit(TOTAL))
                                ),
                                by = c('STUSAB', 'fortnight_end')]

# Add in HHS region & population
tests_nrevss = merge(x = merge(x = tests_nrevss,
                              y = hhs,
                              all.x=TRUE),
                    y = data.frame(STUSAB = toupper(pops$STUSAB),
                                   state_population = pops$Total.),
                    all.x = TRUE)

# summarize NREVSS testing data by HHS region & week (and hhs-fortnight)
tests_nrevss_hhs_wk <- tests_nrevss[,
                  .(
                    'POSITIVE' = sum(POSITIVE, na.rm = T), # I'll assume that if there are multiple rows, that the data were duplicated.
                    'TOTAL'    = sum(TOTAL, na.rm = T),
                    'hhs_population_reporting.nrevss' = sum(state_population[!is.na(tests)], na.rm = T) # total population of the states reporting NREVSS testing data
                  ),
                  by = c('HHS', 'yr_wk')]
# aggregate tests_nrevss by hhs & fortnight (currently no plan to use this)
tests_nrevss_hhs_fn <- tests_nrevss[,
                  .(
                    'POSITIVE' = sum(POSITIVE, na.rm = T), # I'll assume that if there are multiple rows, that the data were duplicated.
                    'TOTAL'    = sum(TOTAL, na.rm = T),
                    'hhs_population_reporting.nrevss' = sum(state_population[!is.na(tests)], na.rm = T) # total population of the states reporting NREVSS testing data
                  ),
                  by = c('HHS', 'fortnight_end')]
# tests_nrevss_hhs_wk is the info that we'll actually use in weighting (for now)


# create a column to group source by type (NS3, Contractor, or baseline tagged)
us.dat[ source == 'NS3' , source_type := 'NS3']
# us.dat[, prop.table(table(source_type))]
us.dat[ (toupper(contractor_vendor_name) != "NULL") &
          (toupper(contractor_vendor_name) != "NS3" &
             source != 'NS3'),
        source_type := 'Contractor']
# us.dat[, prop.table(table(source_type))]
us.dat[ is.na(source_type) & primary_sampling_strategy == 'Baseline_Surveillance', source_type := 'Tagged']
# us.dat[, prop.table(table(source_type))]

# Transform the received_date field as date in R
#us.dat$received_date_dt = as.Date(sapply(strsplit(x = us.dat$received_date, split = ' '), FUN = function(x) x[1]))
us.dat$received_date_dt = as.Date(us.dat$received_date)
## Creating a trimmed down survey dataset
# Removing submission/receive dates for now, for consistency with frozen dataset
svy.dat = data.table::data.table(
  STUSAB = us.dat$primary_state_abv,
  HHS    = us.dat$HHS,
  yr_wk  = us.dat$yr_wk,
  DAY    = us.dat$DAY,
  date   = us.dat$collection_date,
  nt_id  = us.dat$primary_nt_id,
  csid   = us.dat$csid,
  covv_accession_id = us.dat$covv_accession_id,
  received_date = us.dat$received_date_dt,
  # ID number just to help keep track of individual sequences
  myID = 1:nrow(us.dat),
  # SUBM_DT = us.dat$covv_subm_date,
  # CDC_DT  = us.dat$contractor_receive_date_to_cdc,
  # AGE     = us.dat$age_group,
  LAB     = toupper(us.dat$source),
  source_type     = us.dat$source_type,
  SGTF_UPSAMPLING = (us.dat$contractor_targeted_sequencing %in% "Screened for S dropout"),
  # SOURCE  = us.dat$source,
  VARIANT = us.dat$lineage,
  S_MUT   = us.dat$s_mut,
  s1_key  = us.dat$s1_key,
  geni_included = us.dat$geni_included
  #S1_GROUP = us.dat$group_name,   # s1_group names, only the groups > 0.05% in the analysis weeks will have names based on comparison to chosen reference lineage
  #S1_GROUP_KEY = us.dat$group_key # s1_group_keys
)
# add on expanded lineage (if present)
if('expanded_lineage' %in% names(us.dat)) svy.dat$expanded_lineage <- us.dat$expanded_lineage

# replace NAs with "other
# svy.dat$LAB[is.na(svy.dat$LAB)] = "OTHER"
svy.dat[ is.na(LAB), LAB := "OTHER"]

# redefine "week" field
# add a column for number of weeks since start of 2020
# svy.dat$week = as.numeric(as.Date(svy.dat$yr_wk) - week0day1)/7
svy.dat[, 'week' := as.numeric(as.Date(svy.dat$yr_wk) - week0day1)/7]
# merge in state population data
svy.dat = merge(x = svy.dat,
                y = data.table::data.table(
                  STUSAB = toupper(pops$STUSAB),
                  state_population = pops$Total.
                ),
                all.x = TRUE)#merge in state population

# Cleaning state tagged sequences-----------------------------------------------
# unique(svy.dat$LAB)

# add in a column so count occurrences using aggregate (could also use table or xtab)
# svy.dat$count = 1
# print("Unique labs include:")
# aggregate(count ~ LAB, data = svy.dat, FUN = sum)

# print('table of omicron variants before removing any labs:')
# table(svy.dat$VARIANT[grep(pattern = '(B\\.1\\.1\\.529)|(BA\\.[0-9])', x = svy.dat$VARIANT)])

# Clean lab names.
svy.dat[, 'LAB2' := as.character(LAB)]
unique_labs <- unique(svy.dat$LAB)

lab_rename_rules <- list(
  list(selector = regex_lab_selector("(MARYLAND DEPARTMENT OF HEALTH)|(MD PHL)"), new_name = "MD-DPH"),
  list(selector = regex_lab_selector("(New Jersey.+Public Health)|(NJ.PHEL)|(NJ.+Public Health)"), new_name = "NJ-DPH"),
  list(selector = regex_lab_selector("(Texas Department of state health)|(TXDSHS)"), new_name = "TX-DPH"),
  list(selector = regex_lab_selector("LSU"), new_name = "LSU LAB"),
  list(selector = regex_lab_selector("Orange County"), new_name = "Orange County PHL"),
  list(selector = regex_lab_selector("(LAURING LAB.+Michigan)|(Michigan.+Lauring Lab)"), new_name = "LAURING LAB, UNIVERSITY OF MICHIGAN"),
  list(selector = regex_lab_selector("HELIX"), new_name = "HELIX"),
  list(selector = regex_lab_selector("SOUTH DAKOTA public health"), new_name = "SD-DPH"),
  list(selector = regex_lab_selector("(omega.+MOUNES)|(oemga.+mounes)"), new_name = "OMEGA DIAGNOSTICS AT MOUNES"),
  list(selector = regex_lab_selector("Minnesota department of health"), new_name = "MN-DPH"),
  list(selector = regex_lab_selector("Sonoma county public health"), new_name = "SONOMA COUNTY PHL"),
  list(
    selector = regex_lab_selector(
      "(INFECTIOUS DISEASES,  NC SLPH COVID-19 RESPONSE TEAM)|(INFECTIOUS DISEASES,  NORTH CAROLINA STATE LABORATORY OF PUBLIC HEALTH COVID-19 RESPONSE TEAM)"
    ),
    new_name = "INFECTIOUS DISEASES, NC SLPH COVID-19 RESPONSE TEAM"
  ),
  list(
    selector = function(labs) {
      labs[
        grepl(
          pattern = "NORTH CAROLINA STATE LABORATORY OF PUBLIC HEALTH",
          x = labs,
          ignore.case = TRUE
        ) &
        !grepl(
          pattern = "(INFECTIOUS DISEASES)|(COVID-19 RESPONSE TEAM)",
          x = labs,
          ignore.case = TRUE
        )
      ]
    },
    new_name = "NCSLPH"
  ),
  list(selector = regex_lab_selector("^UNMC COVID.*RESPONSE TEAM"), new_name = "UNMC COVID-19 RESPONSE TEAM"),
  list(selector = regex_lab_selector("COUNTY OF SAN LUIS OBISPO"), new_name = "COUNTY OF SAN LUIS OBISPO PHL"),
  list(selector = regex_lab_selector("SAN JOAQUIN COUNTY PUBLIC HEALTH"), new_name = "SAN JOAQUIN COUNTY PHL"),
  list(selector = regex_lab_selector("SUMAN DAS LAB"), new_name = "DR. SUMAN DAS LAB - VANDERBILT UNIVERSITY MEDICAL CENTER"),
  list(selector = regex_lab_selector("BOISE VA MEDICAL CENTER"), new_name = "BOISE VA MEDICAL CENTER")
)

lab_rename_results <- apply_lab_rename_rules(
  svy_dat = svy.dat,
  unique_labs = unique_labs,
  rules = lab_rename_rules
)
svy.dat <- lab_rename_results$svy_dat
labnames_df <- lab_rename_results$labnames_df

# Remove GCID sequences before any broader Broad aggregation is introduced.
print(paste('Removing',
    svy.dat[ (LAB == toupper("Broad Institute's Genomic Center For Infectious Diseases (GCID)") & date >= '2023-01-01'), length(date) ],
    'sequences from "Broad Institute\'s Genomic Center For Infectious Diseases (GCID)" collected since 2023-01-01'))
svy.dat <- svy.dat[ !(LAB == toupper("Broad Institute's Genomic Center For Infectious Diseases (GCID)") & date >= '2023-01-01') ]

# save the list of lab names that were changed to file
write_backup_csv(
  x = labnames_df,
  data_dir = data_dir,
  backup_date = data_date,
  custom_tag = custom_tag,
  filename = paste0("lab_name_updates_", data_date, ".csv")
)

# Get counts of samples by lab
# check_count <- aggregate(formula = count ~ LAB2,
#                          data = svy.dat,
#                          FUN = sum)
check_count <- as.data.frame(table(svy.dat$LAB2))
names(check_count) <- c('LAB2', 'count')
save_backup_rds_file(
  object = check_count,
  data_dir = data_dir,
  backup_date = data_date,
  custom_tag = custom_tag,
  filename = paste0(data_date, "_sequence_counts_by_lab", custom_tag, ".RDS")
)
# print a list of cleaned lab names to the console
# check_count

# Use counts of sequences by lab-week to look for issues with submissions
# - 1. new labs that have been added (helpful for finding typos)
# - 2. weeks where a labs samples have gone DOWN (I'm assuming that samples aren't frequently removed)
# - 3. weeks > 8 weeks in the past where labs samples have increased (I'm not sure how how often this happens)
# - 4. labs that normally submit a lot of sequences submitting a different number of sequences than expected
# NOTE: 2&3 only look at net/total number of samples. So if 1 sample is removed & 1 is added, it will slip through the current checks. To do a more accurate check, it would require aggregating by lab-week-nt_id.
{
   # number of samples per lab per week
   counts_by_week <- svy.dat[, .(count = .N), by = c('LAB2', 'yr_wk')]

   # compare to previous week(s)
   previous_file <- find_previous_backup_rds(
     data_dir = data_dir,
     current_date = data_date,
     custom_tag = custom_tag,
     suffix = "_sequence_counts_by_lab_week"
   )

   # if a previous file was found, compare it to the new one to find:
   # - 1. new labs that have been added
   # - 2. weeks where a labs samples have gone DOWN
   # - 3. weeks > 8 weeks in the past where labs samples have increased
   # - 4. labs that normally submit a lot of sequences submitting a different number of sequences than expected
   if(!is.na(previous_file)){
      # read in the older counts
      previous_counts_by_week <- readRDS(previous_file)

      # 1. new labs that have been added: (to aid in checking for typos)
      newly_added_labs <- setdiff(unique(counts_by_week$LAB2), unique(previous_counts_by_week$LAB2))
      if(length(newly_added_labs) == 0) print('\nNo newly added labs.') else {
         print('Compare newly added labs to all labs to look for typos:')
         print('\nnewly added lab names:')
         print(sort(newly_added_labs))
         print("")
         print("All lab names")
         print(sort(unique(counts_by_week$LAB2)))
         print("")
      }

      # - 2. weeks where a labs samples have gone DOWN
      cbw <- merge(x = counts_by_week[, .(LAB2  = LAB2,
                                          yr_wk = yr_wk,
                                          N_new = count)],
                   y = previous_counts_by_week[, .(LAB2  = LAB2,
                                                   yr_wk = yr_wk,
                                                   N_old = count)],
                   by = c('LAB2', 'yr_wk'),
                   all = T)
      # calculate the number of sequences added or subtracted
      cbw[, 'N_net' := N_new - N_old]

      n_lost_print_threshold <- 100
      if(nrow(cbw[N_net < 0])){
         temp <- cbw[N_net <= 0]
         temp[, 'new_date':= data_date]
         temp[, 'old_date':= regmatches(previous_file,regexpr('[0-9]{4}-[0-9]{2}-[0-9]{2}',previous_file))]
         write_backup_csv(
           x = temp,
           data_dir = data_dir,
           backup_date = data_date,
           custom_tag = custom_tag,
           filename = paste0(data_date, "_lost_sequences", custom_tag, ".csv")
         )
         if(nrow(cbw[N_net <= -n_lost_print_threshold])){
          print(paste0('At least ', n_lost_print_threshold, ' sequences were lost from these lab-week combinations:'))
          print(temp[N_net <= -n_lost_print_threshold][order(N_net),], nrow=500) # nrow is a data.table-specific argument
         }
      }

      # - 3. weeks > 8 weeks in the past where labs samples have increased
      n_added_print_threshold <- 100
      temp <- cbw[ yr_wk <= (data_date - 8*7) & N_net > 0]
      if(nrow(temp)> 0){
         write_backup_csv(
           x = temp,
           data_dir = data_dir,
           backup_date = data_date,
           custom_tag = custom_tag,
           filename = paste0(data_date, "_old_sequence_additions", custom_tag, ".csv")
         )
        if(nrow(temp[N_net >= n_added_print_threshold]) > 0){
          print(paste0('At least ', n_added_print_threshold, ' sequences were added to these old lab-week combinations:'))
          print(temp[N_net >= n_added_print_threshold][order(N_net, decreasing = TRUE)], nrow=500) # nrow is a data.table-specific argument
        }
      }

      # - 4. labs that normally submit a lot of sequences submitting a different number of sequences than expected
      # this is not implemented yet b/c ideally we'd set up a table of sequence additions by lab, collection week, and submission week, and then look back through time over time to look for unusual submissions, but that would require more than just reading in a single file here...
   }

   # save the counts by lab & week to file
   save_backup_rds_file(
     object = counts_by_week,
     data_dir = data_dir,
     backup_date = data_date,
     custom_tag = custom_tag,
     filename = paste0(data_date, "_sequence_counts_by_lab_week", custom_tag, ".RDS")
   )
}



# NOTE! This doesn't work if new labs add data for past weeks in addition to current weeks. We'll need to load in a prior dataset & compare names that way.
# Steps to improve: 1) find previous version of "sequence_counts_by_lab" file; 2) load it in; 3) compare current lab names to previous lab names (setdiff()); 4) print new lab names to console.
# print('\nnewly added lab names:')
# recent_lab_names <- sort(unique(svy.dat$LAB2[svy.dat$week == max(svy.dat$week)]))
# older_lab_names  <- sort(unique(svy.dat$LAB2[svy.dat$week != max(svy.dat$week)]))
# new_lab_names <- recent_lab_names[ recent_lab_names %notin% older_lab_names ]
# if( length(new_lab_names) > 0) print(new_lab_names) else print('No new lab names in the most recent week of data')


# check which state-level sources have fewer than 100 samples
low_lab <- check_count$LAB2[ check_count$count < 100 ]

# Drop sequences from sources with less than 100 samples or were listed as CDC sequenced
#svy.dat <- svy.dat[svy.dat$LAB2 %notin% c("CDC",low_lab), ]
{
   # sequences that were from "CDC"
   #cdc_seqs <- svy.dat$LAB2 == 'CDC'
   # sequences from labs with few samples
   low_lab_seqs <- svy.dat$LAB2 %in% low_lab

  #  # counts of CDC sequences by week
  #  cdc_by_wk <- us.dat[ cdc_seqs, .('count' = .N), by = 'yr_wk']
  #  cdc_by_wk[, 'yr_wk' := as.Date(yr_wk)]
  #  # merge in the counts of CDC sequences by week
  #  if(nrow(cdc_by_wk) > 0)
  #     dropped_sequences <- rbind(
  #        dropped_sequences[!(week %in% cdc_by_wk$yr_wk & reason == 'n_dropped_CDC_lab')],
  #        cdc_by_wk[, .('week' = yr_wk, 'reason' = 'n_dropped_CDC_lab', 'count' = count)]
  #     )
  #  # fill in 0's for rows that didn't have any sequences dropped
  #  dropped_sequences[is.na(count) & reason == 'n_dropped_CDC_lab', 'count' := 0]

   # counts of low lab sequences by week
   low_by_wk <- svy.dat[ low_lab_seqs, .('count' = .N), by = 'yr_wk']
   low_by_wk[, 'yr_wk' := as.Date(yr_wk)]
   # merge in the counts of low lab sequences by week
   if(nrow(low_by_wk) > 0)
      dropped_sequences <- rbind(
         dropped_sequences[!(week %in% low_by_wk$yr_wk & reason == 'n_dropped_labs_100_seqs')],
         low_by_wk[, .('week' = yr_wk, 'reason' = 'n_dropped_labs_100_seqs', 'count' = count)]
      )
   # fill in 0's for rows that didn't have any sequences dropped
   dropped_sequences[is.na(count) & reason == 'n_dropped_labs_100_seqs', 'count' := 0]
}
svy.dat <- svy.dat[LAB2 %notin% c("CDC",low_lab), ]

# save the dropped_sequences counts b/c that's all that's calculated in this script
write_backup_csv(
  x = dropped_sequences,
  data_dir = data_dir,
  backup_date = data_date,
  custom_tag = custom_tag,
  filename = paste0("dropped_sequence_counts_", data_date, custom_tag, "_v1.csv")
)

save_backup_rds_file(
  object = list(
    tests_weekly = tests_nrevss_wk,
    tests_fortnight = tests_nrevss_fn,
    tests_hhs_weekly = tests_nrevss_hhs_wk,
    tests_hhs_fortnight = tests_nrevss_hhs_fn
  ),
  data_dir = data_dir,
  backup_date = data_date,
  custom_tag = custom_tag,
  filename = paste0(data_date, "_tests_nrevss_aggregated", custom_tag, ".RDS")
)

# Add in NREVSS testing data (grouped by state & week)
# currently no plan to use this
svy.dat = merge(x = svy.dat,
                y = tests_nrevss_wk[, .(STUSAB, yr_wk,
                                       # number of positive tests
                                       'POSITIVE.state.nrevss' = POSITIVE,
                                       # number of total tests
                                       'TOTAL.state.nrevss' = TOTAL)],
                      all.x = TRUE,
                      by = c("STUSAB", "yr_wk"))
# add in NREVSS tests data (grouped by HHS & week)
# this is what we use in updated weights (May 2023)
svy.dat = merge(x = svy.dat,
                y = tests_nrevss_hhs_wk[, .(HHS, yr_wk,
                                           # number of positive tests
                                           'POSITIVE.HHS.nrevss' = POSITIVE,
                                           # number of total tests
                                           'TOTAL.HHS.nrevss' = TOTAL,
                                           # population of the states reporting NREVSS testing
                                           'population_reporting.HHS.nrevss' = hhs_population_reporting.nrevss)],
                      all.x = TRUE,
                      by = c("HHS", "yr_wk"))



# Add in columns that were previously calculated in "weekly_variant_report_nowcast.R" -----
# so that they're calculated once instead of many times

# Define "source" as "LAB2"; (used for survey design weights)
#svy.dat$SOURCE = svy.dat$LAB2
svy.dat[, 'SOURCE' := LAB2]

# create a variable for aggregating counts of different groupings of samples
# Note: it would be better to update the counting methods in weekly_variant_report_nowcast.R, but for now it uses the "count" column & "aggregate()"
svy.dat[, 'count' := 1]

# calculate final date of 2-week bins for each observation
# svy.dat$FORTNIGHT_END = as.character(week0day1 + svy.dat$week%/%2 * 14 + 13)
svy.dat[, 'FORTNIGHT_END' := as.character(week0day1 + week %/% 2 * 14 + 13)]
#svy.dat[, 'MONTH_END' := as.character(week0day1 + week %/% 4 * 28 + 13)]

svy.dat[, MONTH_END := as.character(
  week0day1 + ((week %/% 2L) + ((week %/% 2L) %% 2L)) * 14L + 13L
)]


# add the current week to the source data
#svy.dat$current_week = current_week
svy.dat[, 'current_week' := current_week]


# save the data to file
if(date_frozen_toread == data_date) {
    save(svy.dat,
     file = paste0(script.basename, "/data/", "svydat_", data_date, custom_tag, ".RData"))
} else {
    save(svy.dat,
        file = paste0(script.basename, "/data/", "svydat_", data_date, custom_tag, "_", date_frozen_toread, "_frozendata",".RData"))
}

}

if (sys.nframe() == 0) {
  variant_surveillance_system_main()
}
