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

# set global options
options(
  # set the java virtual machine memory to 8000 MB
  java.parameters = "-Xmx8000m",
  # do not treat strings as factors by default
  stringsAsFactors = FALSE
)

# optparse option list (get command-line arguments)
option_list <- list(
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
    default = NA,
    help    = "User Name",
    metavar = "character"
  ),
  # CDP password
  optparse::make_option(
    opt_str = c("-p", "--password"),
    type    = "character",
    default = NA,
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


# Import Data ------------------------------------------------------------------

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
    Valid options from sc2_archive.analytics_metadata_frozen.date_frozen include: ',
    paste(valid_data_dates, collapse = ', '),
    '\nValid options from sc2_archive.hhs_protect_testing_frozen.date_frozen include: ',
    paste(valid_tests_dates, collapse = ', '),
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
    (SELECT max(', date_frozen, ') as max_frozen
    FROM sc2_archive.analytics_metadata_frozen) as M
    ON to_date(A.date_frozen) = M.max_frozen'
  ),
  'WHERE A.primary_country in ("United States", "USA")'
) # if unavailable, use test_deduplication_cdcncbigisaid_auto for testing

# pull the data from the CDP database to an R data.frame
dat = DBI::dbGetQuery(conn = impala,
                      statement = query)

# Get the Pangolin lineages (at the time of the data)
if(custom_lineages == TRUE) {
  if(current_data){
    # if custom_lineages == TRUE & current_data == TRUE
    # define custom AY.35+ and AY.4.2+ lineages based on amino acid positions of interest
    pangolin = DBI::dbGetQuery(
      conn = impala,
      statement = '
  SELECT DISTINCT
  P.nt_id,
  CASE
      WHEN P.lineage = "AY.4.2" AND udx.substr_range(A.aa_aln, "145;222") = "HV"
      THEN "AY.4.2+"
      WHEN P.lineage = "AY.35" AND udx.substr_range(A.aa_aln, "484") = "Q"
      THEN "AY.35+"
      ELSE P.lineage
  END as lineage,
  P.conflict,
  P.ambiguity_score,
  P.scorpio_call,
  P.scorpio_support,
  P.scorpio_conflict,
  P.version,
  P.pangolin_version,
  P.class_date,
  P.class_qc,
  P.note,
  P.scorpio_version
  FROM sc2_src.pangolin as P
  LEFT JOIN sc2_src.alignments as A
  ON P.nt_id = A.nt_id
  AND A.protein = "S"
  ')
  } else {
    # if custom_lineages = T & current_data = F
    pangolin = DBI::dbGetQuery(
      conn = impala,
      statement = paste0('
        SELECT DISTINCT
        P.primary_nt_id as nt_id,
        CASE
            WHEN P.lineage = "AY.4.2" AND udx.substr_range(A.aa_aln, "145;222") = "HV"
            THEN "AY.4.2+"
            WHEN P.lineage = "AY.35" AND udx.substr_range(A.aa_aln, "484") = "Q"
            THEN "AY.35+"
            ELSE P.lineage
        END as lineage
        FROM sc2_archive.analytics_metadata_frozen as P
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
      statement = '
      SELECT DISTINCT *
      FROM sc2_src.pangolin
      ')
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
  (SELECT max(', date_frozen, ') as max_frozen
    FROM sc2_archive.hhs_protect_testing_frozen
  ) as F
  ON to_date(H.date_frozen) = F.max_frozen
  WHERE H.collection_date is NOT NULL'
))

# rename the columns of the testing data
colnames(tests) = c("collection_date",
                    "STUSAB", # "reporting_state",
                    "INDETERMINATE",
                    "INVALID",
                    "NEGATIVE",
                    "POSITIVE")

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
     (SELECT a.lineage,
             c.variant_type,
             date_add(date_trunc('week', date_add(primary_collection_date, 1)), 5) AS week_ending,
             count(a.primary_virus_name) AS lineage_count,
             z.region_total,
             count(a.primary_virus_name) / z.region_total AS fraction,
             if(count(a.primary_virus_name) / z.region_total >= 0.01, TRUE, FALSE) AS is_one_percent
      FROM sc2_air.analytics_metadata a
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
      GROUP BY a.lineage,
               c.variant_type,
               week_ending,
               z.region_total --order by week_ending, b.hhs_region, lineage
) Q
   WHERE Q.is_one_percent IS TRUE --OR Q.variant_type is not null
GROUP BY lineage) QQ
LEFT JOIN sc2_air.analytics_lineage_corr cor ON QQ.lineage = cor.lineage
WHERE cor.date_range_of_calc LIKE '%US:3mo'"
  ))

  # get the variant names
  voc2_auto = sort(voc2_df$lineage)

  # save the results to file
  saveRDS(object = voc2_auto,
          file = paste0(script.basename,
                        "/data/voc2_auto", data_date, custom_tag, ".RDS"))
}

# end the database connection
DBI::dbDisconnect(conn = impala)

# read in static table of state populations
# (population data may need to be updated)
pops = read.delim(file = paste0(script.basename,"/resources/ACStable_B01001_40_2018_5.txt"))[, c("STUSAB", "Total.")]

# save the various data objects to file as a backup
save(dat, pangolin, baseline, tests, pops, voc2_auto,
     file = paste0(script.basename, "/data/", "variant_survey_dat", data_date, custom_tag, ".RData"))

# Data Cleaning ----------------------------------------------------------------

# The genomic data are subset to U.S. specimens among human hosts. Some light
# cleaning of the sequence data includes dropping records where the state
# abbreviation is incorrect, or the collection data is prior to October 2019.
# An HHS Region variable is merged in.

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
us.dat = subset(x = dat,
                primary_country %in% c("United States", "USA") &
                  primary_host == "Human")
# Subset to exclude faulty state abbreviations
us.dat = subset(x = us.dat,
                nchar(primary_state_abv) == 2)

## exclude duplicates from each table
us.dat   = distinct(us.dat)
pangolin = distinct(pangolin)
baseline = distinct(baseline)

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

# Add a column for the date of the first day of the week
us.dat$yr_wk = as.character(us.dat$collection_date - as.numeric(strftime(us.dat$collection_date, format="%w")))
# Add a column for the # of days since the start date (first Sunday of 2020)
us.dat$DAY = as.numeric(us.dat$collection_date - week0day1)

# exclude unreasonably early dates
us.dat = subset(x = us.dat,
                collection_date >= as.Date("2019-10-01"))

# merge Pangolin lineage into us.dat
us.dat = merge(x = us.dat,
               y = pangolin[, c("nt_id", "lineage")],
               by.x = "primary_nt_id",
               by.y = "nt_id",
               all.x = TRUE)

# fill in missing lineage data with "covv_lineage"
# (this line of code only works if using the dedupe table)
if ("covv_lineage" %in% names(us.dat)) us.dat$lineage = with(us.dat, ifelse(is.na(lineage), covv_lineage, lineage))

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

# convert to date format
tests$collection_date = as.Date(x = tests$collection_date,
                                format = "%Y-%m-%d")

# exclude unreasonable dates
tests = subset(x = tests,
               collection_date >= as.Date("2019-10-01") &
                 collection_date <= data_date) # Sys.Date()

# Add a column for the date of the first day of the week
tests$yr_wk = as.character(tests$collection_date - as.numeric(strftime(tests$collection_date, format="%w")))

# calculate the total number of tests
tests$TOTAL = rowSums(tests[, c("INDETERMINATE", "INVALID", "NEGATIVE", "POSITIVE")],
                      na.rm=TRUE)

# Replace NAs with 0
tests$POSITIVE = ifelse(test = is.na(tests$POSITIVE),
                        yes = 0,
                        no = tests$POSITIVE)

# Aggregate tests by state & week for weighting
tests_wk = expand.grid(STUSAB = unique(tests$STUSAB),
                       yr_wk = unique(tests$yr_wk),
                       stringsAsFactors = FALSE)

# calculate the number of positive tests & total tests for each week & state
for (cc in c("POSITIVE", "TOTAL")) {
  # calculate number of tests
  tests_wk = merge(x = tests_wk,
                   y = data.frame(xtabs(formula = tests[, cc] ~ STUSAB + yr_wk,
                                        data = tests,
                                        subset = TRUE)), # "subset" totally unnecessary, but it keeps Rstudio from complaining
                   all.x = TRUE)

  # rename the newly created column
  names(tests_wk) = gsub(pattern = "[Ff]req",
                         replacement = cc,
                         x = names(tests_wk))
}


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

# estimates of infections based on: https://www.medrxiv.org/content/10.1101/2020.10.07.20208504v2.full-text
# sqrt(state_population/TOTAL) this estimates the geometric mean of the total population versus the state population
# to try to get around the bias of assuming either total or state population

# Estimate the total number of infections based on test positivity rate
test_tallies_wk = within(data = test_tallies_wk,
                         expr = INFECTIONS <- ifelse(test = TOTAL > 0,
                                                     yes = POSITIVE * sqrt(state_population/TOTAL),
                                                     no = 0))
# Number of weeks since start date
test_tallies_wk$week = as.numeric(as.Date(test_tallies_wk$yr_wk) - week0day1)%/%7

# Aggregate infections & populations by HHS region
incidence_by_region = merge(
  x = aggregate(formula = INFECTIONS ~ yr_wk + HHS,
                data = test_tallies_wk,
                FUN = sum),
  y = aggregate(formula = state_population ~ yr_wk + HHS,
                data = test_tallies_wk,
                FUN = sum)
)

# calculate infection rate
incidence_by_region$HHS_INCIDENCE = incidence_by_region$INFECTIONS / incidence_by_region$state_population

## SGTF over-sampling weights
# (only correcting for the Helix upsampling; other labs were not specifically SGTF up-sampling)

# counts of lineages by contractor name
sgtf.1 = table(us.dat$lineage,
               paste(us.dat$contractor_vendor_name,
                     us.dat$contractor_targeted_sequencing))

# Get the labs/columns with s-gene target failure oversampling
# (identified by having more "dropout" tests than "Illumina" tests)
sgtf.vars = rownames(sgtf.1)[
  sgtf.1[, grep(pattern = "dropout$", x = colnames(sgtf.1))] > sgtf.1[, grep(pattern = "Illumina $", x = colnames(sgtf.1))] ]
# might need to update this line based on the updated lab names? "Illumina $" doesn't match any, but "Illumina" does...
# but doesn't really matter since there's no known SGTF over-sampling being done in fall 2021.

# Smoothed weights: use set of states and weeks where targeted samples were sequenced
sgtf.sub = unique(subset(x = us.dat,
                         contractor_targeted_sequencing %in% "Screened for S dropout")[, c("primary_state_abv", "yr_wk")])

# If there are data with SGTF flag, run a binomial model to estimate probability
# of oversampling
if (nrow(sgtf.sub) > 0) {

  #linear model that fits sgtf lineage weights
  sgtf.glm = glm(
    formula = I(lineage %in% sgtf.vars) ~ I(contractor_vendor_name %in% "Helix/Illumina") + primary_state_abv + yr_wk,
    family = "binomial",
    data = subset(x = us.dat,
                  yr_wk %in% sgtf.sub$yr_wk &
                    primary_state_abv %in% sgtf.sub$primary_state_abv))

  # Define sgtf_weights as the relative probability of being in one of the
  # oversampled lineages if from the lab that did oversampling vs any other lab
  # (assumes that the samples tested by each lab had the same proportions of
  #  different variants)
  sgtf.sub$sgtf_weights =
    predict(object = sgtf.glm,
            newdata = cbind(contractor_vendor_name = c("Helix/Illumina"),
                            sgtf.sub),
            type = "response") /
    predict(object = sgtf.glm,
            newdata = cbind(contractor_vendor_name = c("Other"),
                            sgtf.sub),
            type = "response")
} else {
  # if there are no rows that had oversampling, then just set all the weights to 0
  sgtf.sub$sgtf_weights = numeric(0)
}

# add the SGTF oversampling weights into the test_tallies
test_tallies_wk = merge(x = test_tallies_wk,
                        y = sgtf.sub,
                        by.x = c("STUSAB", "yr_wk"),
                        by.y = c("primary_state_abv", "yr_wk"),
                        all.x = TRUE)

# Set weights to 1 if weights are NA
test_tallies_wk$sgtf_weights = ifelse(test = is.na(test_tallies_wk$sgtf_weights),
                                      yes = 1,
                                      no = test_tallies_wk$sgtf_weights)

## Creating a trimmed down survey dataset
# Removing submission/receive dates for now, for consistency with frozen dataset
svy.dat = data.frame(STUSAB = us.dat$primary_state_abv,
                     HHS    = us.dat$HHS,
                     yr_wk  = us.dat$yr_wk,
                     DAY    = us.dat$DAY,
                     # SUBM_DT = us.dat$covv_subm_date,
                     # CDC_DT = us.dat$contractor_receive_date_to_cdc,
                     # AGE = us.dat$age_group,
                     LAB     = toupper(us.dat$source),
                     SGTF_UPSAMPLING = (us.dat$contractor_targeted_sequencing %in% "Screened for S dropout"),
                     SOURCE  = us.dat$source,
                     VARIANT = us.dat$lineage,
                     S_MUT   = us.dat$s_mut)

# replace NAs with "other
svy.dat$LAB[is.na(svy.dat$LAB)] = "OTHER"

# add a column for number of weeks since start of 2020
svy.dat$week = as.numeric(as.Date(svy.dat$yr_wk) - week0day1)/7

# merge in state population data
svy.dat = merge(x = svy.dat,
                y = data.frame(STUSAB = toupper(pops$STUSAB),
                               state_population = pops$Total.),
                all.x = TRUE)#merge in state population

# Cleaning state tagged sequences-----------------------------------------------
# unique(svy.dat$LAB)

# add in a column so count occurances using aggregate (could also use table or xtab)
svy.dat$count=1
print("Unique labs include:")
aggregate(count~LAB,data=svy.dat,FUN=sum)

#clean LAB names
svy.dat$LAB2 = as.character(svy.dat$LAB)

# get a vector of unique lab names
unique_labs <- unique(svy.dat$LAB)

# Aggregate Maryland lab names
MD_labs_to_agg <- grep(pattern = "(Maryland)|(MD)",
                       x = unique_labs,
                       ignore.case = T,
                       value = TRUE)
labnames_df_md <- data.frame(old_name = MD_labs_to_agg,
                             new_name = "MD-DPH")
svy.dat[svy.dat$LAB %in% MD_labs_to_agg,"LAB2"] <- "MD-DPH"
# Aggregate New Jersey lab names
NJ_labs_to_agg <- grep(pattern = "(New Jersey)|(NJ)",
                       x = unique_labs,
                       ignore.case = T,
                       value = TRUE)
labnames_df_nj <- data.frame(old_name = NJ_labs_to_agg,
                             new_name = "NJ-DPH")
svy.dat[svy.dat$LAB %in% NJ_labs_to_agg,"LAB2"] <- "NJ-DPH"
# Aggregate Texas lab names
TX_labs_to_agg <- grep(pattern = "(Texas)|(TX)",
                       x = unique_labs,
                       value = TRUE,
                       ignore.case = T)
labnames_df_tx <- data.frame(old_name = TX_labs_to_agg,
                             new_name = "TX-DPH")
svy.dat[svy.dat$LAB %in% TX_labs_to_agg,"LAB2"] <- "TX-DPH"
# Aggregate CDC lab names
CDC_labs_to_agg <- grep(pattern = "Centers for Disease Control and Prevention",
                        x = unique_labs,
                        ignore.case = T,
                        value = TRUE)
labnames_df_cdc <- data.frame(old_name = CDC_labs_to_agg,
                              new_name = "CDC")
svy.dat[svy.dat$LAB %in% CDC_labs_to_agg,"LAB2"] <- "CDC"
# Aggregate LSU lab names
LSU_labs_to_agg <- grep(pattern = "LSU",
                        x = unique_labs,
                        ignore.case = T,
                        value = T)
labnames_df_lsu <- data.frame(old_name = LSU_labs_to_agg,
                              new_name = "LSU LAB")
svy.dat[svy.dat$LAB %in% LSU_labs_to_agg,"LAB2"] <- "LSU LAB"
# Aggregate Orange County lab names
OC_labs_to_agg <- grep(pattern = "Orange",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_oc <- data.frame(old_name = OC_labs_to_agg,
                             new_name = "Orange County PHL")
svy.dat[svy.dat$LAB %in% OC_labs_to_agg,"LAB2"] <- "Orange County PHL"
# Aggregate Lauring Lab names
LL_labs_to_agg <- grep(pattern = "LAURING LAB",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_ll <- data.frame(old_name = LL_labs_to_agg,
                             new_name = "LAURING LAB, UNIVERSITY OF MICHIGAN")
svy.dat[svy.dat$LAB %in% LL_labs_to_agg,"LAB2"] <- "LAURING LAB, UNIVERSITY OF MICHIGAN"
# Aggregate Helix lab names
HE_labs_to_agg <- grep(pattern = "HELIX",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_he <- data.frame(old_name = HE_labs_to_agg,
                             new_name = "HELIX")
svy.dat[svy.dat$LAB %in% HE_labs_to_agg,"LAB2"] <- "HELIX"
# Aggregate South Dakota lab names
SD_labs_to_agg <- grep(pattern = "SOUTH DAKOTA",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_sd <- data.frame(old_name = SD_labs_to_agg,
                             new_name = "SD-DPH")
svy.dat[svy.dat$LAB %in% SD_labs_to_agg,"LAB2"] <- "SD-DPH"
# Aggregate Mounes lab names
OM_labs_to_agg <- grep(pattern = "MOUNES",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_om <- data.frame(old_name = OM_labs_to_agg,
                             new_name = "OMEGA DIAGNOSTICS AT MOUNES")
svy.dat[svy.dat$LAB %in% OM_labs_to_agg,"LAB2"] <- "OMEGA DIAGNOSTICS AT MOUNES"
# Aggregate Minnesota lab names
MN_labs_to_agg <- grep(pattern = "Minnesota department of health",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_mn <- data.frame(old_name = MN_labs_to_agg,
                             new_name = "MN-DPH")
svy.dat[svy.dat$LAB %in% MN_labs_to_agg,"LAB2"] <- "MN-DPH"
# Aggregate Sonoma lab names
SO_labs_to_agg <- grep(pattern = "Sonoma",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_so <- data.frame(old_name = SO_labs_to_agg,
                             new_name = "SONOMA COUNTY PHL")
svy.dat[svy.dat$LAB %in% SO_labs_to_agg,"LAB2"] <- "SONOMA COUNTY PHL"

# Aggregate NC labs 
NC_labs_to_agg <- grep(pattern = "(INFECTIOUS DISEASES,  NC SLPH COVID-19 RESPONSE TEAM)|(INFECTIOUS DISEASES,  NORTH CAROLINA STATE LABORATORY OF PUBLIC HEALTH COVID-19 RESPONSE TEAM)",
                       x = unique_labs,
                       ignore.case = T,
                       value = T)
labnames_df_nc <- data.frame(old_name = NC_labs_to_agg,
                             new_name = "INFECTIOUS DISEASES, NC SLPH COVID-19 RESPONSE TEAM")
svy.dat[svy.dat$LAB %in% NC_labs_to_agg,"LAB2"] <- labnames_df_nc$new_name[1]

# Aggregate other NC labs 
NC2_labs_to_agg <- unique_labs[grepl(pattern = "NORTH CAROLINA STATE LABORATORY OF PUBLIC HEALTH",
                       x = unique_labs,
                       ignore.case = T) & 
                         !grepl(pattern = "(INFECTIOUS DISEASES)|(COVID-19 RESPONSE TEAM)",
                               x = unique_labs,
                               ignore.case = T)]
labnames_df_nc2 <- data.frame(old_name = NC2_labs_to_agg,
                             new_name = "NCSLPH")
svy.dat[svy.dat$LAB %in% NC2_labs_to_agg,"LAB2"] <- labnames_df_nc2$new_name[1]



# create a dataframe of all the lab names that were changed
labnames_df <- rbind(
  labnames_df_md,
  labnames_df_nj,
  labnames_df_tx,
  labnames_df_cdc,
  labnames_df_lsu,
  labnames_df_oc,
  labnames_df_ll,
  labnames_df_he,
  labnames_df_sd,
  labnames_df_om,
  labnames_df_mn,
  labnames_df_so,
  labnames_df_nc,
  labnames_df_nc2
)
# print the list of lab names that were changed to the console
labnames_df
# save the list of lab names that were changed to file
write.csv(x = labnames_df,
          file = paste0('./data/lab_name_updates_', data_date, '.csv'),
          row.names = F)

# Get counts of samples by lab
check_count <- aggregate(formula = count ~ LAB2,
                         data = svy.dat,
                         FUN = sum)

# print a list of cleaned lab names to the console
check_count

# check which state-level sources have fewer than 100 samples
low_lab <- check_count$LAB2[ check_count$count < 100 ]

# Drop sequences from sources with less than 100 samples or were listed as CDC sequenced
svy.dat <- svy.dat[svy.dat$LAB2 %notin% c("CDC",low_lab), ]


# Create survey weights --------------------------------------------------------

#Three separate survey designs are used in this analysis:

# * Unweighted, to estimate variant prevalence among sequenced samples
# * Weighted for estimation among test positives
# * Weighted for estimation among infections, but using strategies of unproven reliability in this context


## Weights
# Infections to test-positive ratio
# Geometric mean strategy: https://www.medrxiv.org/content/10.1101/2020.10.07.20208504v2.full-text
#   infections/test positives = sqrt(population/number tested)
svy.dat = merge(x = svy.dat,
                y = cbind(test_tallies_wk[, c("STUSAB", "yr_wk")],
                          I_over_POSITIVE_unadj = sqrt(test_tallies_wk$state_population/test_tallies_wk$TOTAL)),
                all.x = TRUE) # using same bias correcting method as above get ratio of infections to positives

# State-week totals of total and positive test counts, and population, added in to enable alternate weighting (2021-03-18)
# sgtf_weights merged in to enable alternate weighting (2021-03-19)
svy.dat = merge(x = svy.dat,
                y = test_tallies_wk[, c("STUSAB", "yr_wk", "POSITIVE", "TOTAL", "INFECTIONS", "sgtf_weights")],
                all.x = TRUE)

# replace NAs with 1
svy.dat$sgtf_weights[is.na(svy.dat$sgtf_weights) | !svy.dat$SGTF_UPSAMPLING] = 1

# Add in HHS region data
# (used for states missing testing data (OH), using the region level incidence)
svy.dat = merge(x = svy.dat,
                y = incidence_by_region[, c("HHS", "yr_wk", "HHS_INCIDENCE")],
                all.x = TRUE)

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