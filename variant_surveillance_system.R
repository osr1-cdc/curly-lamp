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
# opts$user = ""
# opts$password = ""


# Standard production: always use standard Pangolin lineages
custom_lineages = FALSE
custom_tag = ""

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
if(use_previously_imported_data & date_frozen_toread == data_date &
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_data", custom_tag, ".RDS")) &
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_pangolin", custom_tag, ".RDS")) &
    #file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_baseline", custom_tag, ".RDS")) &
    # file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_tests", custom_tag, ".RDS")) &
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_pops", custom_tag, ".RDS")) &
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_tests_nrevss", custom_tag, ".RDS"))){

  print('Reading in previously pulled data')
  dat          <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_data",         custom_tag, ".RDS"))
  dat <- subset(dat, !(dat$covv_accession_id %in% c("EPI_ISL_19791260",
                                                   "EPI_ISL_19791262",
                                                   "EPI_ISL_19791263",
                                                   "EPI_ISL_19791264",
                                                   "EPI_ISL_19791265",
                                                   "EPI_ISL_19791266",
                                                   "EPI_ISL_19791267",
                                                   "EPI_ISL_19791268")))
  dat <- as.data.table(dat)                                                 

  pangolin     <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_pangolin",     custom_tag, ".RDS"))
  #baseline    <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_baseline",     custom_tag, ".RDS"))
  # tests        <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_tests",        custom_tag, ".RDS"))
  tests_nrevss <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_tests_nrevss", custom_tag, ".RDS"))
  pops         <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_pops",         custom_tag, ".RDS"))
  #s1_groups <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_s1_groups", custom_tag, ".RDS"))
  print('Finished reading in data.')
} else if (use_previously_imported_data & date_frozen_toread != data_date &
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", date_frozen_toread, "_data",         custom_tag, ".RDS")) &
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", date_frozen_toread, "_pangolin",     custom_tag, ".RDS")) &
    #file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date,         "_baseline",     custom_tag, ".RDS")) &
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date,          "_tests",        custom_tag, ".RDS")) &
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date,          "_pops",         custom_tag, ".RDS")) &
    file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date,          "_tests_nrevss", custom_tag, ".RDS"))
    #file.exists(paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_s1_groups", custom_tag, ".RDS"))
    ){
  print('Reading in previously pulled data. analytics_metadata frozen date different from test date')
  dat          <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", date_frozen_toread, "_data",         custom_tag, ".RDS"))
  pangolin     <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", date_frozen_toread, "_pangolin",     custom_tag, ".RDS"))
  #baseline    <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date,          "_baseline",     custom_tag, ".RDS"))
  tests        <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date,          "_tests",        custom_tag, ".RDS"))
  tests_nrevss <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date,          "_tests_nrevss", custom_tag, ".RDS"))
  pops         <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date,          "_pops",         custom_tag, ".RDS"))
  #s1_groups <- readRDS(file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_s1_groups", custom_tag, ".RDS"))
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
pops = read.delim(file = paste0(script.basename, "/resources/ACStable_B01001_40_2018_5.txt"))[, c("STUSAB", "Total.")]

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
saveRDS(tests,
  file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_tests", custom_tag, ".RDS")
)
saveRDS(tests_nrevss,
  file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_tests_nrevss", custom_tag, ".RDS")
)
saveRDS(pops,
  file = paste0(script.basename, "/data/backup_", data_date, custom_tag, "/", data_date, "_pops", custom_tag, ".RDS")
)
print('Finished reading in data.')
}
# Data Cleaning ----------------------------------------------------------------

# Things that this code filters on: [as of 2022-09-08]
# 1) only human hosts  # done in sql query
# 2) only US sequences  # done in sql query
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
#baseline    = distinct(baseline)
# tests        = distinct(tests)
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
# CDC_labs_to_agg <- grep(pattern = "Centers for Disease Control and Prevention",
#                         x = unique_labs,
#                         ignore.case = T,
#                         value = TRUE)
# labnames_df_cdc <- data.frame(old_name = CDC_labs_to_agg,
#                               new_name = "CDC")
# #svy.dat[svy.dat$LAB %in% CDC_labs_to_agg,"LAB2"] <- "CDC"
# svy.dat[ LAB %in% CDC_labs_to_agg, "LAB2" := "CDC"]

# Aggregate LSU lab names
# As of 2022-05-26 this only returns "LSUHS EMERGING VIRAL THREAT LABORATORY", so I [Philip Shirk] removed it.
LSU_labs_to_agg <- grep(pattern = "LSU",
                        x = unique_labs,
                        ignore.case = T,
                        value = T)
labnames_df_lsu <- data.frame(old_name = LSU_labs_to_agg,
                              new_name = "LSU LAB")
#svy.dat[svy.dat$LAB %in% LSU_labs_to_agg,"LAB2"] <- "LSU LAB"
svy.dat[LAB %in% LSU_labs_to_agg, "LAB2" := "LSU LAB"]

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
if (length(NC_labs_to_agg)>0){
labnames_df_nc <- data.frame(old_name = NC_labs_to_agg,
                             new_name = "INFECTIOUS DISEASES, NC SLPH COVID-19 RESPONSE TEAM")
# svy.dat[svy.dat$LAB %in% NC_labs_to_agg,"LAB2"] <- labnames_df_nc$new_name[1]
svy.dat[LAB %in% NC_labs_to_agg, "LAB2" := labnames_df_nc$new_name[1]]
} else labnames_df_nc <- NULL

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


# added 2023-05-24: remove sequences from "Broad Institute's Genomic Center For Infectious Diseases (GCID)" that were collected in 2023
# added here before they get aggregated into "BROAD INSTITUTE"
print(paste('Removing',
    svy.dat[ (LAB == toupper("Broad Institute's Genomic Center For Infectious Diseases (GCID)") & date >= '2023-01-01'), length(date) ],
    'sequences from "Broad Institute\'s Genomic Center For Infectious Diseases (GCID)" collected since 2023-01-01'))
svy.dat <- svy.dat[ !(LAB == toupper("Broad Institute's Genomic Center For Infectious Diseases (GCID)") & date >= '2023-01-01') ]




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
  # labnames_df_cdc,
  labnames_df_lsu, # removed 2022-05-26
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
  labnames_df_bva
  # labnames_df_in,
  # labnames_df_wf,
  # # labnames_df_mi,
  # labnames_df_ct,
  # labnames_df_bu,
  # labnames_df_de,
  # labnames_df_ks,
  # labnames_df_ms,
  # labnames_df_hhd,
  # labnames_df_ga,
  # labnames_df_pa,
  # labnames_df_fs,
  # labnames_df_la,
  # labnames_df_sp,
  # labnames_df_pi,
  # labnames_df_broad,
  # labnames_df_um,
  # labnames_df_sflu,
  # labnames_df_prcvsi,
  # labnames_df_unm,
  # labnames_df_ummc
)

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
write.csv(x = dropped_sequences,
          file = paste0(script.basename, "/data/backup_",data_date, custom_tag, "/dropped_sequence_counts_", data_date, custom_tag, "_v1.csv"),
          row.names = F)

# save aggregated NREVSS testing data to file
if(TRUE){
  saveRDS(list(
    # aggregated by state and:
    'tests_weekly'    = tests_nrevss_wk,
    'tests_fortnight' = tests_nrevss_fn,
    # 'tests_4weeks'    = tests_nrevss_state_bins_list, # UPDATE! CREATE THIS
    # aggregated by hhs region and:
    'tests_hhs_weekly'    = tests_nrevss_hhs_wk,
    'tests_hhs_fortnight' = tests_nrevss_hhs_fn
    ),
    file = paste0(script.basename,
                  "/data/backup_",
                  data_date, custom_tag, "/",
                  data_date, "_tests_nrevss_aggregated",
                  custom_tag, ".RDS"))
}

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
