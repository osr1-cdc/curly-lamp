# Background -------------------------------------------------------------------

#Emerging variants of SARS-CoV-2 with increasing share among cases are of potential public health concern.
# A survey design to track variants, based on the sequencing data residing at CDC and weighted to adjust
# for multiplicity of sources as well as variations in sequencing rates, is discussed here.

# Data and methods

#An event recorded in the surveillance system is the collection of a sample for (PCR) testing (labelled by
# date and place of collection) that results in a positive test and is selected for sequencing and
# submission to the system. Prevalence of variants is estimated from sequence data from multiple sources,
# weighted to adjust for variations in sequencing and testing rates over time, source and jurisdiction,
# using a design-based approach. In this analysis:
# * The event time is date of collection of sample for testing
# * The event location is currently set to the state reporting the test result
# * The protocol for selection of positive test samples for sequencing is assumed to be random, as the
#   commercial laboratories are unable to report their selection protocol, and the NS3 instructions neglect
#   to specify any
# * The strata are U.S. states and territories.
# * The clusters are the data streams: CDC NS3 surveillance, and (currently) three contractors -- Laboratory
#   Corporation of America, Quest Diagnostics, and Helix.
# * Inverse-probability-of-selection weights are generated for estimates of prevalence among those sequenced
#   (unweighted), among test positives, as well as among all infections.

## Current data streams and workflow

#The current workflow is a provisional implementation of the methods. It begins with tapping the current
# data streams:
# * Survey sample data: Genomic data (lineage and S gene mutation list), by collection date, reporting
#   state, data source (lab, NS3, etc.) and any known sampling bias (targeted oversampling of "S gene
#   target failure" specimens, etc.), from a server internal to CDC. Other variables available that may
#   be used for stratified analysis or more careful weighting, but not used currently, are demographic
#   information (limited completeness) on the individual tested, and ZIP code of specimen collection.
# * Supplementary sample data: Pangolin lineage data table, from the CDC internal server, for current
#   lineage information
# * Weighting data: Two tables, updated on or close to the date of analysis, are downloaded from HHS Protect:
#   count of PCR tests results (positive, negative, etc.), aggregated by collection date and reporting state,
#   and count of PCR tests results (positive, negative, etc.) from selected labs (currently lab names that
#   contain substrings "labcorp", "laboratory corporation", "helix", "illumina", "quest"), aggregated by
#   collection date, reporting state and lab name.
# * Supplementary weighting data: state population estimates (ACS 2018 5 year set)

# Setup ------------------------------------------------------------------------
library(implyr)   # for connecting to impala database
library(odbc)     # for connecting to databases
library(dplyr)    # for general data wrangling
library(RJDBC)    # for connecting to a database via JDBC driver
library(optparse) # for importing command-line arguments into code
library(data.table) # for speeding up some of the calculations

# set global options
options(
  # set the java virtual machine memory to 8000 MB
  java.parameters = "-Xmx8000m",
  # do not treat strings as factors by default
  stringsAsFactors = FALSE
)

# optparse option list (get command-line arguments)
option_list <- list(
  # Get lineage definition from sc2_dev.nextclade instead of the default sc2_src.pangolin
  optparse::make_option(
    opt_str = c("-n", "--nextclade_pango"),
    type    = "character",
    default = "F",
    help    = "Whether or not to swith to get lineage defintion from sc2_dev.nextclade instead of sc2_src.pangolin (character value of T or F)",
    metavar = "character"
  ),
  # whether or not to use Custom Lineages
  optparse::make_option(
    opt_str = c("-c", "--custom_lineages"),
    type    = "character",
    default = "F",
    help    = "Whether or not to use custom lineages (character value of T or F)",
    metavar = "character"
  ),
  # CDP user ID
  optparse::make_option(
    opt_str = c("-u", "--user"),
    type    = "character",
    default = '',
    help    = "User Name",
    metavar = "character"
  ),
  # CDP password
  optparse::make_option(
    opt_str = c("-p", "--password"),
    type    = "character",
    default = '',
    help    = "Password",
    metavar = "character"
  )
)

# parse options list (put the command line arguments into an R list)
opts <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# If running code interactively, fill in opts values here
# opts$custom_lineages = "F"
# opts$user = ""
# opts$password = ""


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

## get base directory and driver location
# get command-line arguments
initial.options = commandArgs(trailingOnly = FALSE)
# get the name of this R script from the command-line arguments
file.arg.name = "--file="
script.name = sub(pattern = file.arg.name,
                  replacement = "",
                  x = initial.options[grep(pattern = file.arg.name, x = initial.options)])
# get the directory name of this R script
script.basename = dirname(script.name)
# correction for if running from base repo interactively
if(length(script.basename) == 0) {
  script.basename = "."
}

# create data directory
dir.create(path = paste0(script.basename,"/data"),
           showWarnings = FALSE)

# load in config/config.R variables
source(paste0(script.basename, "/config/config.R"))

# path to the jdbc driver
jdbc_driver = paste0(script.basename, "/jdbc/ClouderaImpalaJDBC-2.6.20.1024/ClouderaImpalaJDBC41-2.6.20.1024/")

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
  lineage_table = "sc2_dev.nextclade"
  lineage_field = "nextclade_pango"
  current_data  = TRUE   # If using nextclade_pango, there is no archived frozen data, only current data from nextclade_pango can be used
} else {
  lineage_table = "sc2_src.pangolin"
  lineage_field = "lineage"
}


# Import Data ------------------------------------------------------------------

# If the data was already pulled and you want to just use that data instead of re-pulling it, set here. 
# This is useful if you aggregate some lab names at the end of this code and then want to re-run the
# script after changing which labs get aggregated. 
use_previously_imported_data <- FALSE

# use previously pulled data if it exists
if(use_previously_imported_data &
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_data", custom_tag, ".RDS")) & 
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_pangolin", custom_tag, ".RDS")) & 
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_baseline", custom_tag, ".RDS")) & 
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_tests", custom_tag, ".RDS")) & 
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_pops", custom_tag, ".RDS"))){
      
  print('Reading in previously pulled data')
  dat <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_data", custom_tag, ".RDS"))
  pangolin <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_pangolin", custom_tag, ".RDS"))
  baseline <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_baseline", custom_tag, ".RDS"))
  tests <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_tests", custom_tag, ".RDS"))
  pops <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_pops", custom_tag, ".RDS"))
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

# Ensure that the "data_date" value is valid.
# In order to get:
# 1) CDP data
# 2) Pangolin lineages
# 3) the list of vocs for Run2,
# "data_date" must be one of the "date_frozen" dates in the
# sc2_archive.analytics_metadata_frozen table.
# valid "data_date" values include:
valid_data_dates <- DBI::dbGetQuery(
  conn = impala,
  statement = '
SELECT DISTINCT to_date(date_frozen)
FROM sc2_archive.analytics_metadata_frozen
    ')
# In order to get the tests data, "data_date" must be in the
# sc2_archive.hhs_protect_testing_frozen table.
# valid "data_date" values include:
valid_tests_dates <- DBI::dbGetQuery(
  conn = impala,
  statement = '
SELECT DISTINCT to_date(date_frozen)
FROM sc2_archive.hhs_protect_testing_frozen
    ')

# throw an error if the data_date isn't in the
# sc2_archive.analytics_metadata_frozen table.
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
  #   - start with "csid", "primary", "covv", "ct"
  #   - end with   "zip"
  #   - contain    "targeted", "vendor"
  #   - equal to   "age"
  get.vars = grep(pattern = "(^csid)|(^primary)|(^covv)|(^ct)|(zip$)|(targeted)|(vendor)|(^age$)",
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
    "primary_collection_date"
  )#,"patient_age_cl","covv_patient_age")
}


# specify the database query to run
# (Subset to US sequences for faster download)
query = paste(
  'SELECT',
  paste0('A.', get.vars, collapse=', '),
  paste0(
    ' FROM sc2_archive.analytics_metadata_frozen as A
    INNER JOIN
    (SELECT max(date_frozen) as max_frozen
    FROM sc2_archive.analytics_metadata_frozen ma
    WHERE to_date(ma.date_frozen) = ', date_frozen, ') as M
    ON A.date_frozen = M.max_frozen'
  ),
  'WHERE A.primary_country in ("United States", "USA")'
) # if unavailable, use test_deduplication_cdcncbigisaid_auto for testing

# pull the data from the CDP database to an R data.frame
dat = DBI::dbGetQuery(conn = impala,
                      statement = query)

# Get the Pango lineages (at the time of the data) from the choosen source
if(custom_lineages == TRUE) {
  custom_lineages_sql = paste0('
    CASE
      -- WHEN P.', lineage_field, ' = "AY.4.2" AND udx.substr_range(A.aa_aln, "145;222") = "HV"
      -- THEN "AY.4.2+"
      -- WHEN P.', lineage_field, ' = "AY.35" AND udx.substr_range(A.aa_aln, "484") = "Q"
      -- THEN "AY.35+"
      -- WHEN P.', lineage_field, ' = "BA.1" AND udx.substr_range(A.aa_aln, "346") = "K"
      -- THEN "BA.1+"
      WHEN regexp_like(P.', lineage_field, ', "^BA.1.{0,1}|^BC.{0,1}|^BD.{0,1}")  AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_BA1"
      WHEN regexp_like(P.', lineage_field, ', "^BA.2.75.{0,1}|^B[LMNRY].{0,1}|^C[AB].{0,1}")  AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_BA275"
      WHEN regexp_like(P.', lineage_field, ', "^BA.2.12.1.{0,1}|^BG.{0,1}")  AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_BA2121"
      WHEN regexp_like(P.', lineage_field, ', "^BA.2.{0,1}|^B[HJPS].{0,1}")  AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_BA2"
      WHEN regexp_like(P.', lineage_field, ', "^BA.4.6.{0,1}")  AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_BA46"
      WHEN regexp_like(P.', lineage_field, ', "^BA.4.{0,1}")  AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_BA4"
      WHEN regexp_like(P.', lineage_field, ', "^BF.7.{0,1}")  AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_BF7"
      WHEN regexp_like(P.', lineage_field, ', "^BQ.1.1.{0,1}") AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_BQ11"
      WHEN regexp_like(P.', lineage_field, ', "^BQ.1.{0,1}") AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_BQ1"
      WHEN regexp_like(P.', lineage_field, ', "^BE.1.1.{0,1}|^BQ.{0,1}|^CC.{0,1}") AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_BE11"
      WHEN regexp_like(P.', lineage_field, ', "^BA.5.{0,1}|^B[EFKTUVWZ].{0,1}|^C[DEFG].{0,1}")  AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_BA5"
      WHEN regexp_like(P.', lineage_field, ', "^BA.3.{0,1}|^B.1.1.529") AND udx.substr_range(A.aa_aln, "346") = "T"
      THEN "R346T_B11529"
      ELSE P.', lineage_field, '
  END as lineage')
  if(current_data){
    # if custom_lineages == TRUE & current_data == TRUE
    # define custom AY.35+ and AY.4.2+ lineages based on amino acid positions of interest
    pangolin = DBI::dbGetQuery(
      conn = impala,
      statement = paste0('
  SELECT DISTINCT
  P.nt_id,',
  custom_lineages_sql,
  ' FROM ', lineage_table, ' as P
  LEFT JOIN sc2_src.alignments as A
  ON P.nt_id = A.nt_id
  AND A.protein = "S"
  '))
  } else {
    # if custom_lineages = T & current_data = F
    pangolin = DBI::dbGetQuery(
      conn = impala,
      statement = paste0('
        SELECT DISTINCT
        P.primary_nt_id as nt_id,',
        custom_lineages_sql,
        ' FROM sc2_archive.analytics_metadata_frozen as P
        LEFT JOIN sc2_src.alignments as A
        ON P.primary_nt_id = A.nt_id
        AND A.protein = "S"
        WHERE to_date(P.date_frozen) = ', date_frozen
      ))
  }
} else {
  if(current_data){
    # if custom_lineages = F & current_data = T
    pangolin = DBI::dbGetQuery(
      conn = impala,
      statement = paste0(
      'SELECT DISTINCT nt_id, ', lineage_field, ' as lineage
      FROM ', lineage_table
      ))
  } else {
    # if custom_lineages = F & current_data = F
    pangolin = DBI::dbGetQuery(
      conn = impala,
      statement = paste0(
        '
        SELECT DISTINCT
        P.primary_nt_id as nt_id,
        P.lineage
        FROM sc2_archive.analytics_metadata_frozen as P
        WHERE to_date(P.date_frozen) = ', date_frozen
      ))
  }
}

# download S gene mutation lists, and source (for NS3 and labs)
baseline = DBI::dbGetQuery(
  conn = impala,
  statement = 'SELECT nt_id, source, primary_virus_name, s_mut, collection_date FROM sc2_dev.baselineseq'
)
# state, zip available, also in dedup; S1 slower, built at query


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

# Get the vocs included in run 2
# only if voc2_manual is not set
if(is.na(voc2_manual)){
  # SQL code from: https://cdc.sharepoint.com/teams/NCEZID-OD_CAWG/Shared%20Documents/Forms/AllItems.aspx?FolderCTID=0x01200085513C439903124B82D8FF6016BB819B&id=%2Fteams%2FNCEZID%2DOD%5FCAWG%2FShared%20Documents%2FDAV%2DActivity%2FSOPs%2Fcurrent%5Flineages%5F1%5Fpercent%2Esql&parent=%2Fteams%2FNCEZID%2DOD%5FCAWG%2FShared%20Documents%2FDAV%2DActivity%2FSOPs

  # This is a heavily modified version of "current_lineages_1_percent.sql", that
  # groups by regions & weeks (instead of just weeks). Resulting lineages include
  # any lineage that's above 1% in any region-week combination.

  #   voc2_df = DBI::dbGetQuery(
  #     conn = impala,
  #     statement = paste0(
  #       "SELECT DISTINCT t.lineage
  # FROM
  # (SELECT a.lineage,
  #         -- c.variant_type,
  #         a.week_ending,
  #         HHS.hhs_region,
  #         count(a.primary_virus_name) AS region_week_lineage_total,
  #         rwt.region_week_total,
  #         count(a.primary_virus_name) / rwt.region_week_total AS share
  #  FROM
  #     (SELECT aa.primary_virus_name,
  #             aa.lineage,
  #             aa.primary_state,
  #             aa.primary_collection_date,
  #             date_add(date_trunc('week', date_add(aa.primary_collection_date, 1)), 5) AS week_ending
  #     FROM sc2_archive.analytics_metadata_frozen AS aa
  #     WHERE to_date(aa.date_frozen) = '", data_date, "'
  #     AND aa.primary_country = 'United States'
  #     AND aa.primary_host = 'Human'
  #     AND datediff(date_add(date_trunc('week', date_add( to_timestamp('", data_date, "', 'yyyy-MM-dd') , 1)), 5), date_add(date_trunc('week', date_add(aa.primary_collection_date, 1)), 5)) >= 14
  #     AND datediff(date_add(date_trunc('week', date_add( to_timestamp('", data_date, "', 'yyyy-MM-dd') , 1)), 5), date_add(date_trunc('week', date_add(aa.primary_collection_date, 1)), 5)) < 98
  #     AND (aa.contractor_vendor_id IS NOT NULL OR aa.cdceventid = '1771')
  #     -- end frozen metadata from a particular date
  #     ) as a
  #  LEFT JOIN sc2_src.hhs_regions AS HHS
  #  ON a.primary_state = HHS.state
  #  LEFT JOIN sc2_src.variant_definitions c
  #  ON a.lineage = c.lineage
  #
  #  -- join in weekly region totals
  #  LEFT JOIN
  #         (SELECT -- z.lineage,
  #                 -- cc.variant_type,
  #                 z.week_ending,
  #                 HHS.hhs_region,
  #                 count(z.primary_virus_name) AS region_week_total
  #         FROM
  #             (SELECT az.primary_virus_name,
  #                     az.lineage,
  #                     az.primary_state,
  #                     az.primary_collection_date,
  #                     date_add(date_trunc('week', date_add(az.primary_collection_date, 1)), 5) AS week_ending
  #             FROM sc2_archive.analytics_metadata_frozen AS az
  #             WHERE to_date(az.date_frozen) = '", data_date, "'
  #             AND az.primary_country = 'United States'
  #             AND az.primary_host = 'Human'
  #             AND datediff(date_add(date_trunc('week', date_add( to_timestamp('", data_date, "', 'yyyy-MM-dd') , 1)), 5), date_add(date_trunc('week', date_add(az.primary_collection_date, 1)), 5)) >= 14
  #             AND datediff(date_add(date_trunc('week', date_add( to_timestamp('", data_date, "', 'yyyy-MM-dd') , 1)), 5), date_add(date_trunc('week', date_add(az.primary_collection_date, 1)), 5)) < 98
  #             AND (az.contractor_vendor_id IS NOT NULL OR az.cdceventid = '1771')
  #             -- end frozen metadata from a particular date
  #             ) as z
  #         LEFT JOIN sc2_src.hhs_regions AS HHS
  #         ON z.primary_state = HHS.state
  #         LEFT JOIN sc2_src.variant_definitions cc
  #         ON z.lineage = cc.lineage
  #         WHERE HHS.hhs_region IS NOT NULL
  #         GROUP BY -- z.lineage,
  #                 -- cc.variant_type, -- CDC's VOC, VOI, etc. designations
  #                 z.week_ending,
  #                 HHS.hhs_region
  #         ) as rwt
  #  ON a.week_ending = rwt.week_ending AND HHS.hhs_region = rwt.hhs_region
  #  WHERE HHS.hhs_region IS NOT NULL
  #  AND a.lineage != 'None'
  #  GROUP BY a.lineage,
  #           -- c.variant_type, -- CDC's VOC, VOI, etc. designations
  #           a.week_ending,
  #           HHS.hhs_region,
  #           rwt.region_week_total
  #  -- order by a.week_ending, HHS.hhs_region, a.lineage
  #  ) AS t
  #  WHERE t.share >= 0.01
  #  ORDER BY t.lineage"
  #     ))


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
          max(week_ending) AS most_recent,
          max(fraction) AS max_fraction,
          max(lineage_count) AS max_virus_count,
          to_timestamp('", data_date, "', 'yyyy-MM-dd') AS pull_date
   FROM
     (SELECT l.", lineage_field, " as lineage,
             c.variant_type,
             date_add(date_trunc('week', date_add(primary_collection_date, 1)), 5) AS week_ending,
             count(a.primary_virus_name) AS lineage_count,
             z.region_total,
             count(a.primary_virus_name) / z.region_total AS fraction,
             if(count(a.primary_virus_name) / z.region_total >= 0.01, TRUE, FALSE) AS is_one_percent,
             if(count(a.primary_virus_name) / z.region_total >= 0.005, TRUE, FALSE) AS is_zerofive_percent
      FROM sc2_air.analytics_metadata a
      LEFT JOIN ", lineage_table, " l on a.primary_nt_id = l.nt_id
      LEFT JOIN
        (SELECT date_add(date_trunc('week', date_add(primary_collection_date, 1)), 5) AS week_ending,
                count(za.primary_virus_name) AS region_total
         FROM sc2_air.analytics_metadata za
         WHERE (za.contractor_vendor_id IS NOT NULL OR za.cdceventid = '1771')
           AND za.primary_country = 'United States'
         GROUP BY week_ending) z ON date_add(date_trunc('week', date_add(primary_collection_date, 1)), 5) = z.week_ending
      AND 1=1
      LEFT JOIN sc2_src.variant_definitions c ON a.lineage = c.lineage
      WHERE -- THis is generally the weeks
        datediff(date_add(date_trunc('week', date_add(to_timestamp('", data_date, "', 'yyyy-MM-dd'), 1)), 5), date_add(date_trunc('week', date_add(primary_collection_date, 1)), 5))>=14
        AND datediff(date_add(date_trunc('week', date_add(to_timestamp('", data_date, "', 'yyyy-MM-dd'), 1)), 5), date_add(date_trunc('week', date_add(primary_collection_date, 1)), 5))< 98
        AND (a.contractor_vendor_id IS NOT NULL OR a.cdceventid = '1771')
        AND a.primary_country = 'United States'
      GROUP BY l.", lineage_field, ",
               c.variant_type,
               week_ending,
               z.region_total --order by week_ending, b.hhs_region, lineage
) Q
   WHERE Q.is_one_percent IS TRUE --OR Q.variant_type is not null
    OR (Q.is_zerofive_percent IS TRUE AND Q.week_ending = date_add(date_trunc('week', date_add(now(), 1)), -16))
GROUP BY lineage) QQ
LEFT JOIN sc2_air.analytics_lineage_corr cor ON QQ.lineage = cor.lineage
WHERE cor.date_range_of_calc LIKE '%US:3mo'"
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
pops = read.delim(file = paste0(script.basename, "/resources/ACStable_B01001_40_2018_5.txt"))[, c("STUSAB", "Total.")]

# save the various data objects to file as a backup
# save(dat, pangolin, baseline, tests, pops, voc2_auto,
#      file = paste0(script.basename, "/data/", "variant_survey_dat", data_date, custom_tag, ".RData"))

# create a folder for the backup data (data that do not get used in "weekly_variant_report_nowcast.R")
dir.create(
  path = paste0(script.basename, "/data/backup_", data_date, custom_tag),
  showWarnings = FALSE
)

saveRDS(dat,
  file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_data", custom_tag, ".RDS")
)
saveRDS(pangolin,
  file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_pangolin", custom_tag, ".RDS")
)
saveRDS(baseline,
  file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_baseline", custom_tag, ".RDS")
)
saveRDS(tests,
  file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_tests", custom_tag, ".RDS")
)
saveRDS(pops,
  file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_pops", custom_tag, ".RDS")
)
print('Finished reading in data.')
}
# Data Cleaning ----------------------------------------------------------------

# Things that this code filters on: [as of 2022-01-19]
# 1) only human hosts
# 2) only US sequences
# 3) drop invalid state abbreviations
# 4) drop duplicates
# 5) drop invalid dates (collection_date older than 2019-10-01 or newer than the final day in the analysis)
# 6) drop sequences from labs with < 100 sequences total
# 7) drop invalid lab names (NA's)
# 8) drop invalid variant names (NA's or "None")
# 9) drop invalid weights (NA or infinite values, which can happen if the testing data for a state is really sparse in a given week)

# The genomic data are subset to U.S. specimens among human hosts. Some light
# cleaning of the sequence data includes dropping records where the state
# abbreviation is incorrect, or the collection data is prior to October 2019.
# An HHS Region variable is merged in.

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
                      'n_dropped_duplicates',
                      'n_dropped_not_human',
                      'n_dropped_not_US',
                      'n_dropped_invalid_state',
                      'n_dropped_invalid_date',
                      'n_dropped_CDC_lab',
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
   dropped_sequences <- rbind(
      dropped_sequences[!(week %in% dups_by_wk$week & reason == 'n_dropped_duplicates')],
      dups_by_wk[, .('week' = week, 'reason' = 'n_dropped_duplicates', 'count' = count)]
   )
   # fill in 0's for rows that didn't have any sequences dropped
   dropped_sequences[is.na(count) & reason == 'n_dropped_duplicates', 'count' := 0]
}

## exclude duplicates from each table
dat      = distinct(dat.dt) # use the data.table version to speed up calculation of why sequences are being dropped
pangolin = distinct(pangolin)
baseline = distinct(baseline)
tests    = distinct(tests)
pops     = distinct(pops)

## Some general parameters:
# start day is the first Sunday of 2020
week0day1 = as.Date("2020-01-05")

# HHS regions
HHS_reg = list(HHS1 = c("CT", "ME", "MA", "NH", "RI", "VT"),
               HHS2 = c("NJ", "NY", "PR", "VI"),
               HHS3 = c("DE", "DC", "MD", "PA", "VA", "WV"),
               HHS4 = c("AL", "FL", "GA", "KY", "MS", "NC", "SC", "TN"),
               HHS5 = c("IL", "IN", "MI", "MN", "OH", "WI"),
               HHS6 = c("AR", "LA", "NM", "OK", "TX"),
               HHS7 = c("IA", "KS", "MO", "NE"),
               HHS8 = c("CO", "MT", "ND", "SD", "UT", "WY"),
               HHS9 = c("AZ", "CA", "HI", "NV", "AS", "MP", "GU", "MH"),
               HHS10= c("AK", "ID", "OR", "WA")
)

# Add in populations of US territories
# (b/c they're not in the "ACStable_B01001_40_2018_5.txt" file)
# Territory 2020-07-01 populations from https://en.wikipedia.org/wiki/List_of_states_and_territories_of_the_United_States_by_population (2021-04-19)
# Marshall Islands 2018 estimate https://en.wikipedia.org/wiki/Marshall_Islands (2021-04-19)
pops = rbind(pops,
             data.frame(STUSAB = c("AS", "GU", "MP", "VI", "MH"),
                        `Total.`=c(49437, 168485, 51433, 106235, 58413)))
# create a lookup table with states & their HHS regions
hhs = data.frame(STUSAB = toupper(pops$STUSAB))
# get the HHS region for each state
hhs$HHS = sapply(X = hhs$STUSAB, FUN = grep, HHS_reg)

## Subset by time and place:
# - USA,
# - states with two letter abbreviations,
# - human host,
# - reasonable collection date

# U.S., legitimately coded states, human host
# (dat may already be subset to US, but retained for backward compatibility)
{
   # count sequences being dropped by week
   prim_country_US <- dat$primary_country %in% c("United States", "USA")
   prim_host_human <- dat$primary_host == "Human"

   # counts of non-US sequences by week
   nonUS_by_wk <- dat[ (!prim_country_US) | is.na(prim_country_US), .('count' = .N), by = 'week']
   # merge in the counts of non-US sequences
   if(nrow(nonUS_by_wk) > 0)
      dropped_sequences <- rbind(
         dropped_sequences[!(week %in% nonUS_by_wk$week & reason == 'n_dropped_not_US')],
         nonUS_by_wk[, .('week' = week, 'reason' = 'n_dropped_not_US', 'count' = count)]
      )
   # fill in 0's for rows that didn't have any sequences dropped
   dropped_sequences[is.na(count) & reason == 'n_dropped_not_US', 'count' := 0]

   # counts of non-human sequences by week
   nonhuman_by_wk <- dat[ (!prim_host_human) | is.na(prim_host_human), .('count' = .N), by = 'week']
   # merge in the counts of non-US sequences
   if(nrow(nonhuman_by_wk) > 0)
      dropped_sequences <- rbind(
         dropped_sequences[!(week %in% nonhuman_by_wk$week & reason == 'n_dropped_not_human')],
         nonhuman_by_wk[, .('week' = week, 'reason' = 'n_dropped_not_human', 'count' = count)]
      )
   # fill in 0's for rows that didn't have any sequences dropped
   dropped_sequences[is.na(count) & reason == 'n_dropped_not_human', 'count' := 0]
}
us.dat = subset(x = dat,
                prim_country_US & prim_host_human)

# Subset to exclude faulty state abbreviations
{
   # get counts of faulty state abbreviations by week
   valid_state_abb <- nchar(us.dat$primary_state_abv) == 2

   # counts of faulty state abbreviations by week
   fs_by_wk <- us.dat[ (!valid_state_abb) | is.na(valid_state_abb), .('count' = .N), by = 'week']
   # merge in the counts of faulty state abbreviations
   if(nrow(fs_by_wk) > 0)
      dropped_sequences <- rbind(
         dropped_sequences[!(week %in% fs_by_wk$week & reason == 'n_dropped_invalid_state')],
         fs_by_wk[, .('week' = week, 'reason' = 'n_dropped_invalid_state', 'count' = count)]
      )
   # fill in 0's for rows that didn't have any sequences dropped
   dropped_sequences[is.na(count) & reason == 'n_dropped_invalid_state', 'count' := 0]
}
us.dat = subset(x = us.dat,
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
               y = pangolin[, c("nt_id", "lineage")],
               by.x = "primary_nt_id",
               by.y = "nt_id",
               all.x = TRUE)

# fill in missing lineage data with "covv_lineage"
# (this line of code only works if using the dedupe table)
if ("covv_lineage" %in% names(us.dat)) us.dat$lineage = with(us.dat, ifelse(is.na(lineage), covv_lineage, lineage))

# print out counts of omicron sequences
# print('Table of all omicron sequences:')
# table(us.dat$lineage[grep(pattern = '(B\\.1\\.1\\.529)|(BA\\.[0-9])', x = us.dat$lineage)])

# Merge S-gene mutation (baseline) data into the us.dat
# NS3 identifier in baseline$source:
us.dat = merge(x = us.dat,
               y = baseline[, c("nt_id", "source", "primary_virus_name", "s_mut")],
               by.x = c("primary_nt_id", "primary_virus_name"),
               by.y = c("nt_id", "primary_virus_name"),
               all.x = TRUE)

# Add in HHS regions
us.dat = merge(x = us.dat,
               y = hhs,
               by.x = "primary_state_abv",
               by.y = "STUSAB",
               all.x = TRUE)

# Aggregate state testing data -------------------------------------------------
# Test data are aggregated by week for weighting.

## Testing data by state and date (HHS Protect) for weighting
tests <- data.table::as.data.table(tests)

# convert to date format
# tests$collection_date = as.Date(x = tests$collection_date,
#                                 format = "%Y-%m-%d")
tests[, 'collection_date' := as.Date(x = collection_date,
                                     format = "%Y-%m-%d")]

# exclude unreasonable dates
# tests = subset(x = tests,
#                collection_date >= as.Date("2019-10-01") &
#                  collection_date <= data_date) # Sys.Date()
tests <- tests[collection_date >= as.Date("2019-10-01") &
                 collection_date <= data_date,]

# Add a column for the date of the first day of the week
# tests$yr_wk = as.character(tests$collection_date - as.numeric(strftime(tests$collection_date, format="%w")))
tests[,'yr_wk' := as.character(collection_date - as.numeric(strftime(collection_date, format="%w")))]

# calculate the total number of tests
tests$TOTAL = rowSums(tests[, c("INDETERMINATE", "INVALID", "NEGATIVE", "POSITIVE")],
                      na.rm=TRUE)

# Replace NAs with 0
# tests$POSITIVE = ifelse(test = is.na(tests$POSITIVE),
#                         yes = 0,
#                         no = tests$POSITIVE)
tests[is.na(POSITIVE), 'POSITIVE' := 0]

# Aggregate tests by state & week for weighting
# tests_wk = expand.grid(STUSAB = unique(tests$STUSAB),
#                        yr_wk = unique(tests$yr_wk),
#                        stringsAsFactors = FALSE)
#
# # calculate the number of positive tests & total tests for each week & state
# for (cc in c("POSITIVE", "TOTAL")) {
#   # calculate number of tests
#   tests_wk = merge(x = tests_wk,
#                    y = data.frame(xtabs(formula = tests[, cc] ~ STUSAB + yr_wk,
#                                         data = tests,
#                                         subset = TRUE)), # "subset" totally unnecessary, but it keeps Rstudio from complaining
#                    all.x = TRUE)
#
#   # rename the newly created column
#   names(tests_wk) = gsub(pattern = "[Ff]req",
#                          replacement = cc,
#                          x = names(tests_wk))
# }
tests_wk <- tests[,
                  .(
                    'POSITIVE' = sum(POSITIVE, na.rm = T),
                    'TOTAL'    = sum(INDETERMINATE, INVALID, NEGATIVE, POSITIVE, na.rm = T)
                  ),
                  by = c('STUSAB', 'yr_wk')]
# update some of the column names
# data.table::setnames(x = tests_wk,
#                      old = c('TOTAL', 'POSITIVE'),
#                      new = c('TOTAL_weekly', 'POSITIVE_weekly'))


# Aggregate tests by state & group for alternate weighting
{
  # Aggregate most recent 3 weeks together
  # aggregegate 2 weeks before that
  # aggregate weekly before that

  # most recent 3 weeks:
  final_3_wks <- as.character((as.Date(time_end) - 20) + (0:2)*7)
  previous_2_wks <- as.character( min(as.Date(final_3_wks)) - (2:1)*7 )

  tests$group = tests$yr_wk
  tests$group[ tests$yr_wk %in% final_3_wks] <- final_3_wks[1]
  tests$group[ tests$yr_wk %in% previous_2_wks] <- previous_2_wks[1]

  # tests_group = expand.grid(STUSAB = unique(tests$STUSAB),
  #                           stringsAsFactors = FALSE)
  #
  # # calculate the number of positive tests & total tests for each fortnight & state
  # for (cc in c("POSITIVE", "TOTAL")) {
  #   # calculate number of tests
  #   tests_group = merge(x = tests_group,
  #                       y = data.frame(xtabs(formula = tests[, cc] ~ STUSAB + group,
  #                                            data = tests,
  #                                            subset = TRUE)), # "subset" totally unnecessary, but it keeps Rstudio from complaining
  #                       all.x = TRUE)
  #
  #   # rename the newly created column
  #   names(tests_group) = gsub(pattern = "[Ff]req",
  #                             replacement = cc,
  #                             x = names(tests_group))
  # }
  tests_group <- tests[,
                    .(
                      'POSITIVE' = sum(POSITIVE, na.rm = T),
                      'TOTAL'    = sum(INDETERMINATE, INVALID, NEGATIVE, POSITIVE, na.rm = T)
                    ),
                    by = c('STUSAB', 'group')]
  # update some of the column names
  # data.table::setnames(x = tests_group,
  #                      old = c('TOTAL', 'POSITIVE'),
  #                      new = c('TOTAL_group', 'POSITIVE_group'))
}

# aggregate tests by state & day
# Note: daily aggregation is not for weighting. It's just for calculating the number
#       of PCR confirmed infections are the results of each variant (by multiplying
#       estimated proportions by the number of positive test results).
tests_dy <- tests[,
                  .(
                    'POSITIVE' = sum(POSITIVE, na.rm = T),
                    'TOTAL'    = sum(INDETERMINATE, INVALID, NEGATIVE, POSITIVE, na.rm = T)
                  ),
                  by = c('STUSAB', 'collection_date')]
# change some column names
data.table::setnames(x = tests_dy,
                     old = c('collection_date', 'TOTAL', 'POSITIVE'),
                     new = c('date', 'TOTAL_daily', 'POSITIVE_daily'))

# aggregate tests by state & fortnight
# Note: fortnightly aggregation is not for weighting. It's just for calculating the number
#       of PCR confirmed infections are the results of each variant (by multiplying
#       estimated proportions by the number of positive test results).
tests[, 'week' := as.numeric(as.Date(yr_wk) - week0day1)/7]
tests[, 'fortnight_end' := as.character(week0day1 + week %/% 2 * 14 + 13)]
tests_fn <- tests[,
                  .(
                    'POSITIVE' = sum(POSITIVE, na.rm = T),
                    'TOTAL'    = sum(INDETERMINATE, INVALID, NEGATIVE, POSITIVE, na.rm = T)
                  ),
                  by = c('STUSAB', 'fortnight_end')]

# aggregate tests by state & 4-week period
# Note: this aggregation is not for weighting. It's just for calculating the number
#       of PCR confirmed infections are the results of each variant (by multiplying
#       estimated state-level proportions (calculated in rolling 4-week bins) by
#       the number of positive test results).

tests_state_bins_list <- lapply(X = state_time_end, FUN = function(dd){
  # filter the tests dataset to only include tests in the relevant bin
  tests[ collection_date %in% (dd - 7*4):dd,
         .(
           'POSITIVE' = sum(POSITIVE, na.rm = T),
           'TOTAL'    = sum(INDETERMINATE, INVALID, NEGATIVE, POSITIVE, na.rm = T)
         ),
         by = c('STUSAB')]
})

names(tests_state_bins_list) <- state_time_end

# Weighting---------------------------------------------------------------------
# (weights are calculated here & recalculated in weekly_variant_report_nowcast.R
#  after data are subset.)

#Weights are estimated by treating each source (contracting lab, or CDC surveillance) as a cluster
#   nested with strata (states), calculated for each week (i.e., to ensure representation weekly).
#Two sets of weights are estimated: the first, $w_p$, for representation among (PCR) test positive
#   indivduals; and, second, $w_i$, for representation among all prevalent infections. Assumptions
#   (modifiable as data become available) in estimating weights are:

#  * Each positive (PCR) test is in the sampling frame of one of the source streams, so that for
#    each state, week and source
#$$w_p = \frac{\mbox{number of positive PCR test results}}{\mbox{number of sequences submitted}}$$
#  * Oversampling of SGTF samples by one source results in a reduction in weights of SGTF sequences
#   from that source by a factor that is estimated using a logistic regression model relating the
#   odds of finding an "SGTF variant" by source, state and week.
#   There is no reliable and precise method for this yet, so these weights are subject to
#   considerable uncertainties. Here, a [strategy](https://www.medrxiv.org/content/10.1101/2020.10.07.20208504v2.full-text)
#   based on test positivity is used (other [strategies](https://covid19-projections.com/estimating-true-infections-revisited/)
#   are available, too):
#  $$\frac{\mbox{number of prevalent infections}}{\mbox{number of test positives}} = \sqrt{\frac{\mbox{population of jurisdiction}}{\mbox{number of tests}}}$$
#  If each source (lab) stream is assumed to sample from a base population with the same prevalence
#   of infection as the jurisdiction (state), it can be shown that the weight specific to each source,
#   based on the test positivity of each source, is:
#  $$w_i = \frac{\mbox{number of prevalent infections in source}}{\mbox{number of source positives}}
#= \frac{\mbox{number of source positives}}{\mbox{number of source tests}}
#{\frac{\sqrt{\mbox{population of jurisdiction}\times\mbox{number of tests}}}{\mbox{number of positives}}}$$

#  The "infection" weight for each sequence is $w_p\times w_i$, and depends on stratum (state), time (week of collection) and cluster (source).

## Test tally denominator streams
# merge state populations & HHS region into the tests per week dataset
test_tallies_wk = merge(x = merge(x = tests_wk,
                                  y = hhs,
                                  all.x=TRUE),
                        y = data.frame(STUSAB = toupper(pops$STUSAB),
                                       state_population = pops$Total.),
                        all.x = TRUE)

# by grouping
test_tallies_gp = merge(x = merge(x = tests_group,
                                  y = hhs,
                                  all.x=TRUE),
                        y = data.frame(STUSAB = toupper(pops$STUSAB),
                                       state_population = pops$Total.),
                        all.x = TRUE)

# by day
test_tallies_dly = merge(x = merge(x = tests_dy,
                                   y = hhs,
                                   all.x = TRUE),
                         y = data.frame(STUSAB = toupper(pops$STUSAB),
                                        state_population = pops$Total.),
                         all.x = TRUE)
# by fortnight
test_tallies_fn = merge(x = merge(x = tests_fn,
                                   y = hhs,
                                   all.x = TRUE),
                         y = data.frame(STUSAB = toupper(pops$STUSAB),
                                        state_population = pops$Total.),
                         all.x = TRUE)


# adjust populations based on the number of weeks in each group
test_tallies_gp$group_population = test_tallies_gp$state_population
test_tallies_gp[ test_tallies_gp$group %in% final_3_wks, 'group_population'] = test_tallies_gp[ test_tallies_gp$group %in% final_3_wks, 'state_population'] * length(final_3_wks)
test_tallies_gp[ test_tallies_gp$group %in% previous_2_wks, 'group_population'] = test_tallies_gp[ test_tallies_gp$group %in% previous_2_wks, 'state_population'] * length(previous_2_wks)


# estimates of infections based on: https://www.medrxiv.org/content/10.1101/2020.10.07.20208504v2.full-text
# sqrt(state_population/TOTAL) this estimates the geometric mean of the total population versus the state population
# to try to get around the bias of assuming either total or state population

# Estimate the total number of infections based on test positivity rate
test_tallies_wk = within(data = test_tallies_wk,
                         expr = INFECTIONS <- ifelse(test = TOTAL > 0,
                                                     yes = POSITIVE * sqrt(state_population/TOTAL),
                                                     no = 0))

test_tallies_gp = within(data = test_tallies_gp,
                         expr = INFECTIONS <- ifelse(test = TOTAL > 0,
                                                     yes = POSITIVE * sqrt(group_population/TOTAL),
                                                     no = 0))
# the daily test tallies aren't used for weights, so no need to calculate infections

# Number of weeks since start date
test_tallies_wk$week = as.numeric(as.Date(test_tallies_wk$yr_wk) - week0day1)%/%7
test_tallies_gp$week = as.numeric(as.Date(test_tallies_gp$group) - week0day1)%/%7

# Aggregate infections & populations by HHS region
incidence_by_region = merge(
  x = aggregate(formula = INFECTIONS ~ yr_wk + HHS,
                data = test_tallies_wk,
                FUN = sum),
  y = aggregate(formula = state_population ~ yr_wk + HHS,
                data = test_tallies_wk,
                FUN = sum)
)
incidence_by_region_gp = merge(
  x = aggregate(formula = INFECTIONS ~ group + HHS,
                data = test_tallies_gp,
                FUN = sum),
  y = aggregate(formula = group_population ~ group + HHS,
                data = test_tallies_gp,
                FUN = sum)
)

# calculate infection rate
incidence_by_region$HHS_INCIDENCE = incidence_by_region$INFECTIONS / incidence_by_region$state_population
incidence_by_region_gp$HHS_INCIDENCE_gp = incidence_by_region_gp$INFECTIONS / incidence_by_region_gp$group_population

## SGTF over-sampling weights
# (only correcting for the Helix upsampling; other labs were not specifically SGTF up-sampling)

# # counts of lineages by contractor name
# sgtf.1 = table(us.dat$lineage,
#                paste(us.dat$contractor_vendor_name,
#                      us.dat$contractor_targeted_sequencing))
#
# # Get the labs/columns with s-gene target failure oversampling
# # (identified by having more "dropout" tests than "Illumina" tests)
# sgtf.vars = rownames(sgtf.1)[
#   sgtf.1[, grep(pattern = "dropout$", x = colnames(sgtf.1))] > sgtf.1[, grep(pattern = "Illumina $", x = colnames(sgtf.1))] ]
# # might need to update this line based on the updated lab names? "Illumina $" doesn't match any, but "Illumina" does...
# # but doesn't really matter since there's no known SGTF over-sampling being done in fall 2021.
#
# # Smoothed weights: use set of states and weeks where targeted samples were sequenced
# sgtf.sub = unique(subset(x = us.dat,
#                          contractor_targeted_sequencing %in% "Screened for S dropout")[, c("primary_state_abv", "yr_wk")])
#
# # If there are data with SGTF flag, run a binomial model to estimate probability
# # of oversampling
# if (nrow(sgtf.sub) > 0) {
#
#   #linear model that fits sgtf lineage weights
#   sgtf.glm = glm(
#     formula = I(lineage %in% sgtf.vars) ~ I(contractor_vendor_name %in% "Helix/Illumina") + primary_state_abv + yr_wk,
#     family = "binomial",
#     data = subset(x = us.dat,
#                   yr_wk %in% sgtf.sub$yr_wk &
#                     primary_state_abv %in% sgtf.sub$primary_state_abv))
#
#   # Define sgtf_weights as the relative probability of being in one of the
#   # oversampled lineages if from the lab that did oversampling vs any other lab
#   # (assumes that the samples tested by each lab had the same proportions of
#   #  different variants)
#   sgtf.sub$sgtf_weights =
#     predict(object = sgtf.glm,
#             newdata = cbind(contractor_vendor_name = c("Helix/Illumina"),
#                             sgtf.sub),
#             type = "response") /
#     predict(object = sgtf.glm,
#             newdata = cbind(contractor_vendor_name = c("Other"),
#                             sgtf.sub),
#             type = "response")
# } else {
#   # if there are no rows that had oversampling, then just set all the weights to 0
#   sgtf.sub$sgtf_weights = numeric(0)
# }
#
# # add the SGTF oversampling weights into the test_tallies
# test_tallies_wk = merge(x = test_tallies_wk,
#                         y = sgtf.sub,
#                         by.x = c("STUSAB", "yr_wk"),
#                         by.y = c("primary_state_abv", "yr_wk"),
#                         all.x = TRUE)
#
# # FIX THIS IF WE EVER START USING SGTF WEIGHTING AGAIN!!
# sgtf.sub$group = sgtf.sub$yr_wk
# # sgtf.sub$group[ sgtf.sub$yr_wk %in% final_3_wks] <- final_3_wks[1]
# # sgtf.sub$group[ sgtf.sub$yr_wk %in% previous_2_wks] <- previous_2_wks[1]
#
# test_tallies_gp = merge(x = test_tallies_gp,
#                         y = sgtf.sub[,c('primary_state_abv', 'group', 'sgtf_weights')],
#                         by.x = c("STUSAB", "group"),
#                         by.y = c("primary_state_abv", "group"),
#                         all.x = TRUE)
#
# # Set weights to 1 if weights are NA
# test_tallies_wk$sgtf_weights = ifelse(test = is.na(test_tallies_wk$sgtf_weights),
#                                       yes = 1,
#                                       no = test_tallies_wk$sgtf_weights)
# test_tallies_gp$sgtf_weights = ifelse(test = is.na(test_tallies_gp$sgtf_weights),
#                                       yes = 1,
#                                       no = test_tallies_gp$sgtf_weights)

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
  # ID number just to help keep track of individual sequences
  myID = 1:nrow(us.dat),
  # SUBM_DT = us.dat$covv_subm_date,
  # CDC_DT  = us.dat$contractor_receive_date_to_cdc,
  # AGE     = us.dat$age_group,
  LAB     = toupper(us.dat$source),
  SGTF_UPSAMPLING = (us.dat$contractor_targeted_sequencing %in% "Screened for S dropout"),
  # SOURCE  = us.dat$source,
  VARIANT = us.dat$lineage,
  S_MUT   = us.dat$s_mut
)

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

#clean LAB names
# svy.dat$LAB2 = as.character(svy.dat$LAB)
svy.dat[, 'LAB2' := as.character(LAB)]

# get a vector of unique lab names
unique_labs <- unique(svy.dat$LAB)

# Aggregate Maryland lab names
MD_labs_to_agg <- grep(pattern = "(MARYLAND DEPARTMENT OF HEALTH)|(MD PHL)",
                       x = unique_labs,
                       ignore.case = T,
                       value = TRUE)
labnames_df_md <- data.frame(old_name = MD_labs_to_agg,
                             new_name = "MD-DPH")
# svy.dat[svy.dat$LAB %in% MD_labs_to_agg,"LAB2"] <- "MD-DPH"
svy.dat[ LAB %in% MD_labs_to_agg, "LAB2" := "MD-DPH"]

# Aggregate New Jersey lab names
NJ_labs_to_agg <- grep(pattern = "(New Jersey.+Public Health)|(NJ.PHEL)|(NJ.+Public Health)",
                       x = unique_labs,
                       ignore.case = T,
                       value = TRUE)
labnames_df_nj <- data.frame(old_name = NJ_labs_to_agg,
                             new_name = "NJ-DPH")
#svy.dat[svy.dat$LAB %in% NJ_labs_to_agg,"LAB2"] <- "NJ-DPH"
svy.dat[ LAB %in% NJ_labs_to_agg, "LAB2" := "NJ-DPH"]

# Aggregate Texas lab names
TX_labs_to_agg <- grep(pattern = "(Texas Department of state health)|(TXDSHS)",
                       x = unique_labs,
                       value = TRUE,
                       ignore.case = T)
labnames_df_tx <- data.frame(old_name = TX_labs_to_agg,
                             new_name = "TX-DPH")
# svy.dat[svy.dat$LAB %in% TX_labs_to_agg,"LAB2"] <- "TX-DPH"
svy.dat[ LAB %in% TX_labs_to_agg, "LAB2" := "TX-DPH"]

# Aggregate CDC lab names
CDC_labs_to_agg <- grep(pattern = "Centers for Disease Control and Prevention",
                        x = unique_labs,
                        ignore.case = T,
                        value = TRUE)
labnames_df_cdc <- data.frame(old_name = CDC_labs_to_agg,
                              new_name = "CDC")
#svy.dat[svy.dat$LAB %in% CDC_labs_to_agg,"LAB2"] <- "CDC"
svy.dat[ LAB %in% CDC_labs_to_agg, "LAB2" := "CDC"]

# Aggregate LSU lab names
# As of 2022-05-26 this only returns "LSUHS EMERGING VIRAL THREAT LABORATORY", so I [Philip Shirk] removed it.
# LSU_labs_to_agg <- grep(pattern = "LSU",
#                         x = unique_labs,
#                         ignore.case = T,
#                         value = T)
# labnames_df_lsu <- data.frame(old_name = LSU_labs_to_agg,
#                               new_name = "LSU LAB")
# #svy.dat[svy.dat$LAB %in% LSU_labs_to_agg,"LAB2"] <- "LSU LAB"
# svy.dat[LAB %in% LSU_labs_to_agg, "LAB2" := "LSU LAB"]

# Aggregate Orange County lab names
OC_labs_to_agg <- grep(pattern = "Orange County",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_oc <- data.frame(old_name = OC_labs_to_agg,
                             new_name = "Orange County PHL")
#svy.dat[svy.dat$LAB %in% OC_labs_to_agg,"LAB2"] <- "Orange County PHL"
svy.dat[LAB %in% OC_labs_to_agg, "LAB2" := "Orange County PHL"]

# Aggregate Lauring Lab names
LL_labs_to_agg <- grep(pattern = "(LAURING LAB.+Michigan)|(Michigan.+Lauring Lab)",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_ll <- data.frame(old_name = LL_labs_to_agg,
                             new_name = "LAURING LAB, UNIVERSITY OF MICHIGAN")
#svy.dat[svy.dat$LAB %in% LL_labs_to_agg,"LAB2"] <- "LAURING LAB, UNIVERSITY OF MICHIGAN"
svy.dat[LAB %in% LL_labs_to_agg, "LAB2" := "LAURING LAB, UNIVERSITY OF MICHIGAN"]

# Aggregate Helix lab names
HE_labs_to_agg <- grep(pattern = "HELIX",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_he <- data.frame(old_name = HE_labs_to_agg,
                             new_name = "HELIX")
#svy.dat[svy.dat$LAB %in% HE_labs_to_agg,"LAB2"] <- "HELIX"
svy.dat[LAB %in% HE_labs_to_agg, "LAB2" := "HELIX"]

# Aggregate South Dakota lab names
SD_labs_to_agg <- grep(pattern = "SOUTH DAKOTA public health",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_sd <- data.frame(old_name = SD_labs_to_agg,
                             new_name = "SD-DPH")
#svy.dat[svy.dat$LAB %in% SD_labs_to_agg,"LAB2"] <- "SD-DPH"
svy.dat[LAB %in% SD_labs_to_agg, "LAB2" := "SD-DPH"]

# Aggregate Mounes lab names
OM_labs_to_agg <- grep(pattern = "(omega.+MOUNES)|(oemga.+mounes)",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_om <- data.frame(old_name = OM_labs_to_agg,
                             new_name = "OMEGA DIAGNOSTICS AT MOUNES")
# svy.dat[svy.dat$LAB %in% OM_labs_to_agg,"LAB2"] <- "OMEGA DIAGNOSTICS AT MOUNES"
svy.dat[LAB %in% OM_labs_to_agg, "LAB2" := "OMEGA DIAGNOSTICS AT MOUNES"]

# Aggregate Minnesota lab names
MN_labs_to_agg <- grep(pattern = "Minnesota department of health",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_mn <- data.frame(old_name = MN_labs_to_agg,
                             new_name = "MN-DPH")
#svy.dat[svy.dat$LAB %in% MN_labs_to_agg,"LAB2"] <- "MN-DPH"
svy.dat[LAB %in% MN_labs_to_agg, "LAB2" := "MN-DPH"]

# Aggregate Sonoma lab names
SO_labs_to_agg <- grep(pattern = "Sonoma county public health",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_so <- data.frame(old_name = SO_labs_to_agg,
                             new_name = "SONOMA COUNTY PHL")
#svy.dat[svy.dat$LAB %in% SO_labs_to_agg,"LAB2"] <- "SONOMA COUNTY PHL"
svy.dat[LAB %in% SO_labs_to_agg, "LAB2" := "SONOMA COUNTY PHL"]

# Aggregate NC labs
NC_labs_to_agg <- grep(pattern = "(INFECTIOUS DISEASES,  NC SLPH COVID-19 RESPONSE TEAM)|(INFECTIOUS DISEASES,  NORTH CAROLINA STATE LABORATORY OF PUBLIC HEALTH COVID-19 RESPONSE TEAM)",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_nc <- data.frame(old_name = NC_labs_to_agg,
                             new_name = "INFECTIOUS DISEASES, NC SLPH COVID-19 RESPONSE TEAM")
# svy.dat[svy.dat$LAB %in% NC_labs_to_agg,"LAB2"] <- labnames_df_nc$new_name[1]
svy.dat[LAB %in% NC_labs_to_agg, "LAB2" := labnames_df_nc$new_name[1]]

# Aggregate other NC labs
NC2_labs_to_agg <- unique_labs[grepl(pattern = "NORTH CAROLINA STATE LABORATORY OF PUBLIC HEALTH",
                                     x = unique_labs,
                                     ignore.case = T) &
                                 !grepl(pattern = "(INFECTIOUS DISEASES)|(COVID-19 RESPONSE TEAM)",
                                        x = unique_labs,
                                        ignore.case = T)]
labnames_df_nc2 <- data.frame(old_name = NC2_labs_to_agg,
                              new_name = "NCSLPH")
#svy.dat[svy.dat$LAB %in% NC2_labs_to_agg,"LAB2"] <- labnames_df_nc2$new_name[1]
svy.dat[LAB %in% NC2_labs_to_agg, "LAB2" := labnames_df_nc2$new_name[1]]

# Aggregate UNMC labs
UNMC_labs_to_agg <- grep(pattern = '^UNMC COVID.*RESPONSE TEAM',
                          x = unique_labs,
                          ignore.case = T,
                          value = T)
labnames_df_unmc <- data.frame(old_name = UNMC_labs_to_agg,
                               new_name = 'UNMC COVID-19 RESPONSE TEAM')
#svy.dat[svy.dat$LAB %in% UNMC_labs_to_agg, 'LAB2'] <- labnames_df_unmc$new_name[1]
svy.dat[LAB %in% UNMC_labs_to_agg, 'LAB2' := labnames_df_unmc$new_name[1]]

# Aggregate COUNTY OF SAN LUIS OBISPO labs
SLO_labs_to_agg <- grep(pattern = 'COUNTY OF SAN LUIS OBISPO',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_slo <- data.frame(old_name = SLO_labs_to_agg,
                               new_name = 'COUNTY OF SAN LUIS OBISPO PHL')
#svy.dat[svy.dat$LAB %in% SLO_labs_to_agg, 'LAB2'] <- labnames_df_slo$new_name[1]
svy.dat[LAB %in% SLO_labs_to_agg, 'LAB2' := labnames_df_slo$new_name[1]]

# Aggregate SAN JOAQUIN COUNTY PUBLIC HEALTH labs
SJ_labs_to_agg <- grep(pattern = 'SAN JOAQUIN COUNTY PUBLIC HEALTH',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_sj <- data.frame(old_name = SJ_labs_to_agg,
                               new_name = 'SAN JOAQUIN COUNTY PHL')
#svy.dat[svy.dat$LAB %in% SJ_labs_to_agg, 'LAB2'] <- labnames_df_sj$new_name[1]
svy.dat[LAB %in% SJ_labs_to_agg, 'LAB2' := labnames_df_sj$new_name[1]]

# added 2022-02-17
# Aggregate Suman Das labs
SDL_labs_to_agg <- grep(pattern = 'SUMAN DAS LAB',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_sdl <- data.frame(old_name = SDL_labs_to_agg,
                               new_name = 'DR. SUMAN DAS LAB - VANDERBILT UNIVERSITY MEDICAL CENTER')
#svy.dat[svy.dat$LAB %in% SDL_labs_to_agg, 'LAB2'] <- labnames_df_sdl$new_name[1]
svy.dat[LAB %in% SDL_labs_to_agg, 'LAB2' := labnames_df_sdl$new_name[1]]

# added 2022-04-07
# Aggregate Boise VA
BVA_labs_to_agg <- grep(pattern = 'BOISE VA MEDICAL CENTER',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_bva <- data.frame(old_name = BVA_labs_to_agg,
                               new_name = 'BOISE VA MEDICAL CENTER')
svy.dat[LAB %in% BVA_labs_to_agg, 'LAB2' := labnames_df_bva$new_name[1]]

# Aggregate Indiana State DoH
IN_labs_to_agg <- grep(pattern = 'IN.*STATE DEPARTMENT OF HEALTH LABORATORY SERVICES',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_in <- data.frame(old_name = IN_labs_to_agg,
                             new_name = 'INDIANA STATE DEPARTMENT OF HEALTH LABORATORY SERVICES')
svy.dat[LAB %in% IN_labs_to_agg, 'LAB2' := labnames_df_in$new_name[1]]

# added 2022-04-28 (Wake Forest & Michigan)
# Aggregate Wake Forest names
WF_labs_to_agg <- grep(pattern = 'WAKE FOREST school of medicine',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_wf <- data.frame(old_name = WF_labs_to_agg,
                             new_name = "WAKE FOREST SCHOOL OF MEDICINE, INTERNAL MEDICINE")
svy.dat[LAB %in% WF_labs_to_agg, 'LAB2' := labnames_df_wf$new_name[1]]

# aggregate Michigan DoHHS
MI_labs_to_agg <- grep(pattern = 'MICHIGAN DEPARTMENT of health',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_mi <- data.frame(old_name = MI_labs_to_agg,
                             new_name = "MICHIGAN DEPARTMENT OF HEALTH AND HUMAN SERVICES")
svy.dat[LAB %in% MI_labs_to_agg, 'LAB2' := labnames_df_mi$new_name[1]]

# added 2022-05-12 (Connecticut)
# Aggregate Connecticut names
CT_labs_to_agg <- grep(pattern = '(CT department of public health)|(Connecticut department of public health)',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_ct <- data.frame(old_name = CT_labs_to_agg,
                             new_name = "CT DPH")
svy.dat[LAB %in% CT_labs_to_agg, 'LAB2' := labnames_df_ct$new_name[1]]

# added 2022-05-19 (Bushman)
# Aggregate Bushman lab 
# "THE BUSHMAN LAB, SCHOOL OF MEDICINE, UNIVERSITY OF PENNSYLVANIA"
# "BUSHMAN"
BU_labs_to_agg <- grep(pattern = '(^BUSHMAN$)|(THE BUSHMAN LAB)',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_bu <- data.frame(old_name = BU_labs_to_agg,
                             new_name = "THE BUSHMAN LAB, UPENN")
svy.dat[LAB %in% BU_labs_to_agg, 'LAB2' := labnames_df_bu$new_name[1]]

# added 2022-05-26 (Delaware PHL)
# DELAWARE PUBLIC HEALTH LABORATORY
# DELAWARE PUBLIC HEALTH LAB
DE_labs_to_agg <- grep(pattern = 'DELAWARE PUBLIC HEALTH LAB',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_de <- data.frame(old_name = DE_labs_to_agg,
                             new_name = "DELAWARE PHL")
svy.dat[LAB %in% DE_labs_to_agg, 'LAB2' := labnames_df_de$new_name[1]]

# "KANSAS HEALTH AND ENVIRONMENTAL LAB"
# "KANSAS HEALTH AND ENVIRONMENTAL LABORATORIES"
KS_labs_to_agg <- grep(pattern = 'KANSAS HEALTH AND ENVIRONMENTAL',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_ks <- data.frame(old_name = KS_labs_to_agg,
                             new_name = "KANSAS HEALTH AND ENVIRONMENTAL LAB")
svy.dat[LAB %in% KS_labs_to_agg, 'LAB2' := labnames_df_ks$new_name[1]]

# "MISSISSIPPI PUBLIC HEALTH LABORATORY"
# MS PHL
MS_labs_to_agg <- grep(pattern = '(MISSISSIPPI PUBLIC HEALTH LABORATORY)|(MS PHL)',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_ms <- data.frame(old_name = MS_labs_to_agg,
                             new_name = "MISSISSIPPI PHL")
svy.dat[LAB %in% MS_labs_to_agg, 'LAB2' := labnames_df_ms$new_name[1]]

# "HOUSTON HEALTH DEPARTMENT"
# "HOUSTON HEALTH DEPT."
HHD_labs_to_agg <- grep(pattern = 'HOUSTON HEALTH DEP',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_hhd <- data.frame(old_name = HHD_labs_to_agg,
                             new_name = "HOUSTON HEALTH DEPT")
svy.dat[LAB %in% HHD_labs_to_agg, 'LAB2' := labnames_df_hhd$new_name[1]]

# GEORGIA DEPARTMENT OF PUBLIC HEALTH
# GA DEPARTMENT OF PUBLIC HEALTH
GA_labs_to_agg <- grep(pattern = '(^GEORGIA DEPARTMENT OF PUBLIC HEALTH$)|(^GA DEPARTMENT OF PUBLIC HEALTH$)',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_ga <- data.frame(old_name = GA_labs_to_agg,
                             new_name = "GA DEPARTMENT OF PUBLIC HEALTH")
svy.dat[LAB %in% GA_labs_to_agg, 'LAB2' := labnames_df_ga$new_name[1]]

# added 2022-06-30
# "BUREAU OF LABORATORIES, PENNSYLVANIA DEPARTMENT OF HEALTH"
# "PENNSYLVANIA DEPARTMENT OF HEALTH BUREAU OF LABORATORIES"
PA_labs_to_agg <- grep(pattern = 'PENNSYLVANIA DEPARTMENT OF HEALTH',
                         x = unique_labs,
                         ignore.case = T,
                         value = T)
labnames_df_pa <- data.frame(old_name = PA_labs_to_agg,
                             new_name = "PENNSYLVANIA DEPARTMENT OF HEALTH BUREAU OF LABORATORIES")
svy.dat[LAB %in% PA_labs_to_agg, 'LAB2' := labnames_df_pa$new_name[1]]


# Steps to add more lab aggregations 
# 1. find lab names that almost assuredly refer to the same lab 
# 2. copy-and-paste one of the blocks of code above 
# 3. change "XX_labs_to_agg" and "labnames_df_xx" to new AND UNIQUE names (do a control-F for the new name to make sure it's unique)
# 4. change the regex pattern to something that will return ONLY your set of labs
# 5. add "labnames_df_xx" to "labnames_df" below
# 6. after running the new code, look at ./data/backup_YYYY-MM-DD/lab_name_updates_YYYY-MM-DD.csv to make sure that ONLY the intended labs are being renamed. 

# Other labs that might be duplicates, but that I have not combined:
# 1. "INFECTIOUS DISEASE PROGRAM, BROAD INSTITUTE OF HARVARD AND MIT"
# 1. "BROAD INSTITUTE"

# 2. "RIPHL AT RUSH UNIVERSITY MEDICAL CENTER"
# 2. "RUSH UNIVERSITY MEDICAL CENTER"
# 2. "RHODE ISLAND STATE HEALTH LABORATORY"

# 3. "MASS GENERAL BRIGHAM"
# 3. "MASSACHUSETTS GENERAL HOSPITAL"

# 4. "OREGON SARS-COV-2 GENOME SEQUENCING CENTER"
# 4. "OREGON STATE PUBLIC HEALTH LABORATORY"


# These definitely are not the same lab
# "WADSWORTH CENTER, NEW YORK STATE DEPARTMENT OF HEALTH"
# "NEW YORK CITY PUBLIC HEALTH LABORATORY"

# create a dataframe of all the lab names that were changed
labnames_df <- rbind(
  labnames_df_md,
  labnames_df_nj,
  labnames_df_tx,
  labnames_df_cdc,
  # labnames_df_lsu, # removed 2022-05-26
  labnames_df_oc,
  labnames_df_ll,
  labnames_df_he,
  labnames_df_sd,
  labnames_df_om,
  labnames_df_mn,
  labnames_df_so,
  labnames_df_nc,
  labnames_df_nc2,
  labnames_df_unmc,
  labnames_df_slo,
  labnames_df_sj,
  labnames_df_sdl,
  labnames_df_bva,
  labnames_df_in,
  labnames_df_wf,
  labnames_df_mi,
  labnames_df_ct,
  labnames_df_bu,
  labnames_df_de,
  labnames_df_ks,
  labnames_df_ms,
  labnames_df_hhd,
  labnames_df_ga,
  labnames_df_pa
)

# print the list of lab names that were changed to the console
# print("Labs that were re-named:")
# for (x in seq(nrow(labnames_df))) {
#   print(
#     paste(
#       labnames_df[x,'old_name'], "     to     ", labnames_df[x,'new_name']
#     )
#   )
# }

# print the names of labs that were NOT renamed
# print("")
# print("Labs that were NOT re-named:")
# print(sort(unique_labs[unique_labs %notin% labnames_df$old_name]))

# find similar lab names that remain
# (This doesn't work very well, but it can be helpful to see things a little differently)
{
  # ul2 <- unique(svy.dat$LAB2)
  # sim_mat <- stringdist::stringdistmatrix(a = ul2,
  #                                         b = ul2,
  #                                         method = "lv") # use Levenshtein-distance
  # # remove diagonal
  # diag(sim_mat) <- NA
  # # convert matrix to long data.frame
  # names <- expand.grid(a = ul2,
  #                      b = ul2,
  #                      dist = NA)
  # names$a <- as.character(names$a)
  # names$b <- as.character(names$b)
  # # add in string distances
  # for(a in seq(ul2)){
  #   for(b in seq(ul2)){
  #     names$dist[names$a == ul2[a] &
  #                  names$b == ul2[b]] <- sim_mat[a,b]
  #   }
  # }
  #
  # # remove duplicates
  # names$c = sapply(X = seq(nrow(names)), FUN = function(x) paste(sort(c(names$a[x], names$b[x])), collapse = ' '))
  # names <- names[!duplicated(names[,c('c', 'dist')]),]
  # # return the shortest string distances
  # write.csv(x = names[order(names$dist),c('a', 'b', 'dist')],
  #           file = './data/lab_name_distances.csv',
  #           row.names = F)
}

# save the list of lab names that were changed to file
write.csv(x = labnames_df,
          file = paste0(script.basename, '/data/backup_', data_date, '/lab_name_updates_', data_date, '.csv'),
          row.names = F)

# Get counts of samples by lab
# check_count <- aggregate(formula = count ~ LAB2,
#                          data = svy.dat,
#                          FUN = sum)
check_count <- as.data.frame(table(svy.dat$LAB2))
names(check_count) <- c('LAB2', 'count')
saveRDS(object = check_count,
        file = paste0(script.basename, "/data/backup_", data_date, "/", data_date, "_sequence_counts_by_lab", custom_tag, ".RDS"))
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
   previous_file <- NA
   previous_dates <- 6:24 # search dates from 6 to 24 days ago
   previous_date_counter <- 1
   while(is.na(previous_file) & previous_date_counter <= length(previous_dates)){
      # previous file to look for
      file_to_look_for <- paste0(script.basename, "/data/backup_", data_date - previous_dates[previous_date_counter], "/", data_date - previous_dates[previous_date_counter], "_sequence_counts_by_lab_week", custom_tag, ".RDS")

      if(file.exists(file_to_look_for)) previous_file <- file_to_look_for
      previous_date_counter <- previous_date_counter + 1
   }

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
         write.csv(x = temp,
                   file = paste0(script.basename,
                                 "/data/backup_",
                                 data_date, custom_tag, "/", data_date,
                                 "_lost_sequences",
                                 custom_tag, ".csv"),
                   row.names = F)
         if(nrow(cbw[N_net <= -n_lost_print_threshold])){
          print(paste0('At least ', n_lost_print_threshold, ' sequences were lost from these lab-week combinations:'))
          print(temp[N_net <= -n_lost_print_threshold][order(N_net),], nrow=500) # nrow is a data.table-specific argument
         }
      }

      # - 3. weeks > 8 weeks in the past where labs samples have increased
      n_added_print_threshold <- 100
      temp <- cbw[ yr_wk <= (data_date - 8*7) & N_net > 0]
      if(nrow(temp)> 0){
         write.csv(x = temp,
                   file = paste0(script.basename,
                                 "/data/backup_",
                                 data_date, custom_tag, "/", data_date,
                                 "_old_sequence_additions",
                                 custom_tag, ".csv"),
                   row.names = F)
        if(nrow(temp[N_net >= n_added_print_threshold]) > 0){
          print(paste0('At least ', n_added_print_threshold, ' sequences were added to these old lab-week combinations:'))
          print(temp[N_net >= n_added_print_threshold][order(N_net, decreasing = TRUE)], nrow=500) # nrow is a data.table-specific argument
        }
      }

      # - 4. labs that normally submit a lot of sequences submitting a different number of sequences than expected
      # this is not implemented yet b/c ideally we'd set up a table of sequence additions by lab, collection week, and submission week, and then look back through time over time to look for unusual submissions, but that would require more than just reading in a single file here...
   }

   # save the counts by lab & week to file
   saveRDS(object = counts_by_week,
           file = paste0(script.basename,
                         "/data/backup_",
                         data_date, custom_tag,
                         "/", data_date,
                         "_sequence_counts_by_lab_week",
                         custom_tag,
                         ".RDS"))
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
   cdc_seqs <- svy.dat$LAB2 == 'CDC'
   # sequences from labs with few samples
   low_lab_seqs <- svy.dat$LAB2 %in% low_lab

   # counts of CDC sequences by week
   cdc_by_wk <- us.dat[ cdc_seqs, .('count' = .N), by = 'yr_wk']
   cdc_by_wk[, 'yr_wk' := as.Date(yr_wk)]
   # merge in the counts of CDC sequences by week
   if(nrow(cdc_by_wk) > 0)
      dropped_sequences <- rbind(
         dropped_sequences[!(week %in% cdc_by_wk$yr_wk & reason == 'n_dropped_CDC_lab')],
         cdc_by_wk[, .('week' = yr_wk, 'reason' = 'n_dropped_CDC_lab', 'count' = count)]
      )
   # fill in 0's for rows that didn't have any sequences dropped
   dropped_sequences[is.na(count) & reason == 'n_dropped_CDC_lab', 'count' := 0]

   # counts of low lab sequences by week
   low_by_wk <- us.dat[ low_lab_seqs, .('count' = .N), by = 'yr_wk']
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
write.csv(x = dropped_sequences,
          file = paste0(script.basename, "/data/backup_",data_date, custom_tag, "/dropped_sequence_counts_", data_date, custom_tag, "_v1.csv"),
          row.names = F)

# print('table of omicron sequences from labs with > 100 sequences:')
# table(svy.dat$VARIANT[grep('(B\\.1\\.1\\.529)|(BA\\.[0-9])', x = svy.dat$VARIANT)])
print('table of omicron sequences that will be included in analysis:')
table(svy.dat$VARIANT[ svy.dat$LAB2 != 'OTHER'][grep('(B\\.1\\.1\\.529)|(BA\\.[0-9])|(BC\\.[0-9])|(BD\\.[0-9])|(BE\\.[0-9])|(BF\\.[0-9])|(BG\\.[0-9])', x = svy.dat$VARIANT[ svy.dat$LAB2 != 'OTHER'])])

# Create survey weights --------------------------------------------------------

#Three separate survey designs are used in this analysis:

# * Unweighted, to estimate variant prevalence among sequenced samples
# * Weighted for estimation among test positives
# * Weighted for estimation among infections, but using strategies of unproven reliability in this context


## Weights
# Infections to test-positive ratio
# Geometric mean strategy: https://www.medrxiv.org/content/10.1101/2020.10.07.20208504v2.full-text
#   infections/test positives = sqrt(population/number tested)
# svy.dat = merge(x = svy.dat,
#                 y = cbind(test_tallies_wk[, c("STUSAB", "yr_wk")],
#                           I_over_POSITIVE_unadj = sqrt(test_tallies_wk$state_population/test_tallies_wk$TOTAL)),
#                 all.x = TRUE) # using same bias correcting method as above get ratio of infections to positives

# State-week totals of total and positive test counts, and population, added in to enable alternate weighting (2021-03-18)
# sgtf_weights merged in to enable alternate weighting (2021-03-19)
svy.dat = merge(x = svy.dat,
                y = test_tallies_wk[, c("STUSAB", "yr_wk", "POSITIVE", "TOTAL")],
                all.x = TRUE,
                by = c('STUSAB', 'yr_wk'))

# add in grouped test tallies
{
  # test_tallies_gp$POSITIVE_gp     = test_tallies_gp$POSITIVE
  # test_tallies_gp$TOTAL_gp        = test_tallies_gp$TOTAL
  data.table::setnames(x = test_tallies_gp,
                       old = c('TOTAL', 'POSITIVE'),
                       new = c('TOTAL_gp', 'POSITIVE_gp'))

  # svy.dat$group = svy.dat$yr_wk
  # svy.dat$group[ svy.dat$yr_wk %in% final_3_wks] <- final_3_wks[1]
  # svy.dat$group[ svy.dat$yr_wk %in% previous_2_wks] <- previous_2_wks[1]
  svy.dat[, group := yr_wk]
  svy.dat[ yr_wk %in% final_3_wks, group := final_3_wks[1]]
  svy.dat[ yr_wk %in% previous_2_wks, group := previous_2_wks[1]]

  svy.dat = merge(x = svy.dat,
                  y = test_tallies_gp[, c("STUSAB",
                                          "group",
                                          "group_population",
                                          "POSITIVE_gp",
                                          "TOTAL_gp")],
                  all.x = TRUE,
                  by = c('STUSAB', 'group'))
}

# save aggregated testing data to file
if(TRUE){
  saveRDS(list('tests_daily'  = test_tallies_dly,
               'tests_group'  = test_tallies_gp,
               'tests_weekly' = test_tallies_wk,
               'tests_fortnight' = test_tallies_fn,
               'tests_4weeks' = tests_state_bins_list),
          file = paste0(script.basename,
                        "/data/backup_",
                        data_date, custom_tag, "/",
                        data_date, "_tests_aggregated",
                        custom_tag, ".RDS"))
}

# Add in HHS region data
# (used for states missing testing data (OH), using the region level incidence)
svy.dat = merge(x = svy.dat,
                y = incidence_by_region[, c("HHS", "yr_wk", "HHS_INCIDENCE")],
                all.x = TRUE,
                by = c("HHS", "yr_wk"))
incidence_by_region_gp$yr_wk = incidence_by_region_gp$group
svy.dat = merge(x = svy.dat,
                y = incidence_by_region_gp[, c("HHS", "yr_wk", "HHS_INCIDENCE_gp")],
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

# add the current week to the source data
#svy.dat$current_week = current_week
svy.dat[, 'current_week' := current_week]


# save the data to file
save(svy.dat,
     file = paste0(script.basename, "/data/", "svydat_", data_date, custom_tag, ".RData"))



#--------------- Limitations, etc.---------------

#Merging multiple disparate data streams into a single unified survey design is achieved
# here by assuming each stream is a cluster in a PPS design. The analysis is limited by
# the limitations of this assumption, and by incomplete information about the sampling
# design within each stream. Owing to incomplete and imprecise test data, sources with
# substantial contribution to the sequence database may be underweighted.

# * Weights are derived using count of test results (including negatives) by state and
#   collection date, currently from HHS Protect. There are a few states missing, and it
#   has proved to be difficult to assign tests unambiguoisly to their source.
# * The protocols used by different sources to select for sequences are unknown; an
#   assumption of random selection is made here (with modifications for oversampling of
#   SGTF samples, where noted), leading to the possibility of biased estimates.
# * Current [NS3 guidelines for sampling](https://www.aphl.org/programs/preparedness/Crisis-Management/Documents/FEB%202021%20Revised%20NS3%20Submission%20Guidance_02052021%20FNL.pdf) are inadequate:

#  > Ideally, specimens should represent geographic, demographic (e.g., age), and clinical
#   (e.g., disease severity or outcome) diversity from across the jurisdiction. This can be
#   achieved through random selection of specimens collected within the last 7 days.


# Changelog

# * 2021-03-07: updated lineage with pangolin as master-list
# * 2021-03-07: updated test stats by collection date (some states missing!)
# * 2021-03-07: updated designation of year-week to the date of first Sunday of week
# * 2021-03-08: calibrated infections-over-positive weights to test-positivity by cluster
# * 2021-03-08: switched to a PPS design to more accurately reflect subsampling within clusters
# * 2021-03-09: added HHS region
# * 2021-04-XX: added populations for territories; removed lab based weights; infection rates
#   imputed using HHS Regional rates for states with missing testing data; removed
#   submission/receive dates from svy.dat; sgtf_weights are now weights even for SGTF_UPSAMPLING = FALSE

# Known issues/next steps

# * Smooth out positve-over-sequenced weights, if needed
# * Find test counts by collection date for states consistently missing at HHS Protect
# * Adjust multinomial logistic growth rate models for R_t
# * Locate all SGTF sequences using in silico PCR [row 5 (TCAACTCAGGACTTGTTCTTACCT) https://www.protocols.io/view/multiplexed-rt-qpcr-to-screen-for-sars-cov-2-b-1-1-br9vm966/materials]? (Probably not needed any longer)
