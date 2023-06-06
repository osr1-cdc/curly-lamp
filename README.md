# SARS-CoV-2 Proportion Modeling and Nowcast

Code for applying a weights and confidence intervals to SARS-CoV-2 proportion data using a survey design-based approach.

 **Written by** Prabasaj Paul <vig5@cdc.gov>, Molly Steele <xzn9@cdc.gov>
 **with updates from** Norman Hassell <ncy6@cdc.gov>, Philip Shirk <rsv4@cdc.gov>, Xiaoyu Sherry Zheng <qiu5@cdc.gov>

**Reference**:
Paul P, France AM, Aoki Y, et al. Genomic Surveillance for SARS-CoV-2 Variants Circulating in the United States, December 2020–May 2021. MMWR Morb Mortal Wkly Rep 2021;70:846–850. DOI: [http://dx.doi.org/10.15585/mmwr.mm7023a3](http://dx.doi.org/10.15585/mmwr.mm7023a3)


## Code
- **config/config.R** - Specify various configuration settings, such as data dates, VOCs to include in particular runs, and figure settings.
- **variant_surveillance_system.R** - Generates the analytic dataset with the survey weights
- **variant_surveillance_system.sh** - Wrapper script for submitting variant_surveillance_system.R to run in HPC.
   - qsub variant_surveillance_system.sh <`username`> <`password`> <T/F for custom lineage> <T/F for using nextclade pango calls>
- **weekly_variant_report_nowcast.R** - Creates variant proportion estimates using the dataset created in `variant_surveillance_system.R`. 
   - `weekly_variant_report_nowcast.R` accomodates 3 different "runs", each of which has its own set of "vocs" (i.e. variants for which to calculate proportions)
      1) Run 1 calculates variant share/proportion and confidence intervals estimated using survey design for both fortnights (HHS regions & nationally) and weeks (HHS regions & nationally) for vocs that are shown on CDC COVID data tracker website (voc1). 
      2) Run 2 generates the same 2 files as Run 1, but with vocs specified for Run 2 (voc2). Run 2 voc2s are usually automatically generated based on preset selection criteria. Currently, Run 2 autoselected vocs include all variants that occur at a frequency >= 1% of the unweighted sequence data in any of the 2-week periods from -1 to -7 2-week periods and variants at a frequency >= 0.5% of the unweighted sequence data in the -1 2-week period (See table below). All variants from Run 1 are also added to Run 2 voc2. Additionally, any variant of special interest can also be added to voc2 using the `voc2_additional` variable. Run 2 fits model-based smoothed trends in variant share to both national and HHS regional estimates (i.e. "_Nowcast_"). Run 2 also generates nowcast estimates for voc list for Run1 (voc1) by aggregating estimates for lineages from Run 2 to their corresponding parental lineages.
         ### Run 2 vocs auto selection criteria

         | fortnight |      -7     |      -6     |      -5     |     -4    |     -3    |     -2    |      -1     |  Current Fortnight  |
         |:---------:|:-----------:|:-----------:|:-----------:|:---------:|:---------:|:---------:|:-----------:|      :-------------------:|
         |    week   | -15 and -14 | -13 and -12 | -11 and -10 | -9 and -8 | -7 and -6 | -5 and -4 |  -3 and -2  | -1 and Current Week |
         | Nowcast   |  >=1% (unw)  |  >=1% (unw)  |  >=1% (unw)  | >=1% (unw) | >=1% (unw) | >=1% (unw) | >=0.5% (unw) |                     |
      3) Run3 generates state-level estimates in rolling 4 wk bins using survey design (same as Run 1). NO LONGER NEEDED. 
- **proportion_modeling_run.sh** - Wrapper script for running variant_surveillance_system.sh and weekly_variant_report_nowcast.R.
- **weekly_s1_variant_report_nowcast.R** - Creates s1_species proportion estimates using the dataset created in `variant_surveillance_system.R`. 
   - qsub s1_run.sh <`username`> <`password`> <`reference lineage`>


## Data Requirements
- Sequence data (compiled by SSEV Bioinformatics infrastructure team) 
   - Source: CDP database `sc2_archive.analytics_metadata_frozen`
- Pangolin lineages 
   - Source: default CDP database tables `sc2_src.pangolin`, same information is also preserved in the `lineage` field of `sc2_archive.analytics_metadata_frozen`. The `nextclade_pango` field of `sc2_src.nextclade` can also be used instead; however when using `nextclade_pango`, the 4th argument for variant_surveillance_system.sh should be 'T', and aggregation codes need to be manually updated
- NREVSS testing data
   - Source: CDP database `sc2_archive.nrevss_frozen`
- HHS protect RT-PCR testing data
   - Source: CDP database `sc2_archive.hhs_protect_testing_frozen`
- For auto-voc2 selection
   - Source: CDP database `sc2_air.analytics_metadata`, `sc2_src.variant_definitions`, `sc2_air.analytics_lineage_corr`
- Population Data - State population data (as of 2018) 
   - Source: `./resources/ACStable_B01001_40_2018_5.txt`

## Other requirements
- Cloudera Impala JDBC driver: provided in `./jdbc/ClouderaImpalaJDBC-2.6.20.1024/ClouderaImpalaJDBC41-2.6.20.1024`

## Output Files
Result folder: `paste0("/scicomp/groups/Projects/SARS2Seq/repos/sc2_proportion_modeling/results_", data_date, '_', results_tag)`
   1) Run 1
      - Pre-modeling Lineage aggregation results: `paste0("lineage_aggregataion_summary_KGCI_svyNEW_", data_date, "state_tag_included_Run1.csv")`
      - Fornightly weighted estimates: `paste0("variant_share_weighted_KGCI_svyNEW_", data_date, "state_tag_included_Run1", results_tag,".csv")`
      - Weekly weighted estimates: `paste0("variant_share_weekly_weighted_KGCI_svyNEW_", data_date, "state_tag_included_Run1", results_tag,".csv")`      
      - Fornightly weighted estimates formatted for hadoop table `sc2_archive.proportion_modeling`: `paste0("variant_share_weighted_KGCI_svyNEW_", data_date, "state_tag_included_Run1", results_tag,"_hadoop.csv")`
      - Weekly weighted estimates formatted for hadoop table `sc2_archive.proportion_modeling`: `paste0("variant_share_weekly_weighted_KGCI_svyNEW_", data_date, "state_tag_included_Run1", results_tag,"_hadoop.csv")`

   2) Run 2 
   - *NOTE! Run 2 produces results for both Run 2 AND Run 1!*
      - Pre-modeling Lineage aggregation results: `paste0("lineage_aggregataion_summary_KGCI_svyNEW_", data_date, "_state_tag_included_Run2.csv")`
      - Fornightly weighted estimates: `paste0("variant_share_weighted_KGCI_svyNEW_", data_date, "_state_tag_included_Run2", results_tag,".csv")`
      - Weekly weighted estimates: `paste0("variant_share_weekly_weighted_KGCI_svyNEW_", data_date, "_state_tag_included_Run2", results_tag,".csv")`
      - Fornightly nowcast estimates: `paste0("updated_nowcast_fornightly_", weighted_methods, "_", data_date, "_state_tag_included_Run2", results_tag,".csv")`
      - Weekly nowcast estimates: `paste0("updated_nowcast_weekly_", weighted_methods, "_", data_date, "_state_tag_included_Run2", results_tag,".csv")`
      - Fornightly weighted estimates formatted for hadoop table `sc2_archive.proportion_modeling`: `paste0("variant_share_weighted_KGCI_svyNEW_", data_date, "_state_tag_included_Run2", results_tag,"_hadoop.csv")`
      - Weekly weighted estimates formatted for hadoop table `sc2_archive.proportion_modeling`: `paste0("variant_share_weekly_weighted_KGCI_svyNEW_", data_date, "_state_tag_included_Run2", results_tag,"_hadoop.csv")`
      - Fornightly nowcast estimates formatted for hadoop table `sc2_archive.proportion_modeling`: `paste0("updated_nowcast_fornightly_", data_date, "_state_tag_included_Run2", results_tag,"_hadoop.csv")`
      - Weekly nowcast estimates formatted for hadoop table `sc2_archive.proportion_modeling`: `paste0("updated_nowcast_weekly_", data_date, "_state_tag_included_Run2", results_tag,"_hadoop.csv")`\
   - The following files are nowcast estimates for the Run1 vocs
      - Fornightly nowcast estimates: `paste0("updated_nowcast_fornightly_", weighted_methods, "_", data_date, "state_tag_included_Run1", results_tag,".csv")`
      - Weekly nowcast estimates: `paste0("updated_nowcast_weekly_", weighted_methods, "_", data_date, "_state_tag_included_Run1", results_tag,".csv")`
      - Fornightly nowcast estimates formatted for hadoop table `sc2_archive.proportion_modeling`: `paste0("updated_nowcast_fornightly_", data_date, "_state_tag_included_Run1", results_tag,"_hadoop.csv")`
      - Weekly nowcast estimates formatted for hadoop table `sc2_archive.proportion_modeling`: `paste0("updated_nowcast_weekly_", data_date, "_state_tag_included_Run1", results_tag,"_hadoop.csv")`
      - Post-modeling lineage aggregation results (aggregateing voc2 variants to their parental lineages in voc1): `paste0(agg_var_mat_KGCI_svyNEW_", data_date, "_state_tag_included_Run2.csv")`
   - Image files
      - National data
         - `paste0("/results/wtd_shares_",data_date,"_","barplot_US",tag,".jpg")`
         - `paste0("/results/wtd_shares_",data_date,"_","projection_US",tag,".jpg")`
         - `paste0("/results/wtd_shares_",data_date,"_","growthrate_US",tag,".png")`
      - For each HHS region
         - `paste0("/results/wtd_shares_",data_date,"_","barplot_HHS",hhs,tag,".jpg")`
         - `paste0("/results/wtd_shares_",data_date,"_","projection_HHS",hhs,tag,".jpg")`
         - `paste0("/results/wtd_shares_",data_date,"_","growthrate_HHS", hhs,tag,".jpg")`
   3) Run3 generates state-level estimates in rolling 4 wk bins using survey design (same as Run 1). NO LONGER NEEDED. 
    - State level weighted estimates: `paste0("/state_weighted_roll4wk_KGCI_svyNEW_", data_date, "_state_tag_included_Run3.csv")`
    - State level weighted estimates formatted for hadoop table `sc2_archive.state_proportion_modeling`: `paste0(/state_weighted_roll4wk_KGCI_svyNEW_", data_date, "_state_tag_included_Run3_hadoop.csv")`


## How to run from a scicomp location
It is best to run this from rosalind (`rosalind.biotech.cdc.gov`). Option A and Option B are two ways to run the codes. Choose either one of them to run.

### Option A
This is the simplest way to run when variant list for CDC COVID data tracker has been determined and no other special modifications/testings are needed.

1. Navigate to the folder from which you will run the analyses. Typically this is the shared project folder. 
   ```bash
   cd /scicomp/groups/Projects/SARS2Seq/repos/sc2_proportion_modeling
   ```

2. Make sure the files are up to date with the git repository
   ```bash 
	git status 
	# make sure it is on the appropriate branch, master branch should generally be used
	git fetch --all
	git status
    # if not up-to-date, pull updates from remote repository
    git pull
   ```
3. Get the list of autoselected voc2s by running the [voc2_autoselection sql](https://cdp-01.biotech.cdc.gov:8889/hue/editor?editor=74361) in Hue. Decide the list of variants to be shown on CDC COVID data tracker website (voc1) and any additional variants to add to voc2 (that have not receached the selection criteria).

   - To change voc1 list (Variants shown in CDT): change `voc1` definition in `config\config.R`
   - Make sure all voc1s are included in `voc2_additional` definition in `config\config.R`
   - To add additional variants to voc2: change `voc2_additionl` definition in `config\config.R`

4. Make sure the configuation settings are all correct for the current runs

   - `data_date`: Should be one of the `date_frozen` from `sc2_archive.analytics_metadata_frozen`. This value would be used in the result folder, result file names and the `analysis_date` field when the data is ingested to `sc2_archive.proportion_modeling`.\
   `date_frozen` can be checked in Hue using the following sql
      ```sql
      select distinct date_frozen from sc2_archive.analytics_metadata_frozen 
      order by date_frozen desc
      ```
   
   - `results_tag`: This value would be used in the result folder, result file names and the `note` field when the data is ingested to `sc2_archive.proportion_modeling`.

5. Make sure the email notification settings are all correct in the bash scripts below by changing the value in the following line
   ```bash
   #$ -M <user email address>
   ```

   - proportion_modeling_run.sh
   - variant_surveillance_system.sh
   - run1_trim.sh
   - run2_trim.sh

6. Run the wrapper script by submitting it as a HPC job to the server.
   ```bash
   qsub proportion_modeling_run.sh -u <username> -p <cdp password> -c <'T' or 'F' for custom lineage>
   ```
   
   This wrapper script submits the job for running variant_surveillance_system.sh first. When that is finished, it submits the jobs for run1_trim.sh and run2_trim.sh. So user will receive individual notification emails for the above three jobs and finally a notification for finishing this wrapper job.

7. Check the wrapper job logfile and .err file to make sure all steps run successfully.\
   The log file is `/scicomp/groups/Projects/SARS2Seq/repos/sc2_proportion_modeling/log`

8. Check results in the result folder: 
   ```paste0("/scicomp/groups/Projects/SARS2Seq/repos/sc2_proportion_modeling/results_", data_date, '_', results_tag)```

   a. Check .err and .out files to make sure no running error occured.\
      The following error message can be ignored:\
          ```WARNING: sun.reflect.Reflection.getCallerClass is not supported. This will impact performance. ERROR StatusLogger No Log4j 2 configuration file found. Using default configuration (logging only errors to the console), or user programmatically provided configurations. Set system property 'log4j2.debug' to show Log4j 2 internal initialization logging. See https://logging.apache.org/log4j/2.x/manual/configuration.html for instructions on how to configure Log4j 2```

   b. Move the corresponding .err and .out files to result folder.\
      ```bash
      mv Run1_<results_tag>.[eo]* results/results_<data_date>_<results_tag>/
      ```
   
   c. automated aggregation is implemented based on the `extended_lineage` field in `sc2_src.pangolin`. This only works when custom lineage option is NOT used and nextclade_pango option is NOT used. lineage aggregation results can be checked in the following files:
      - `lineage_aggregation_summary...` list the aggregation done before the actuall weighting and modeling processes.
      - `agg_var_mat...` lists the aggregation done after the nowcast modeling processes are finished in Run2, to aggregate variants in voc2 but not in voc1 to their parental lineages in voc1, and generate the nowcast results for Run 1.
   
9. Backup the code that was used for a production run to the git repository
   ```bash
   git status
   git add config/config.R # and other files modified for the run
	git commit -m 'Production runs: YYYY-MM-DD' # and any other brief notes for changes made for the week
	git push
   ```
10. Input run information in [SC2_Proportion_Modeling_Run_Records](https://cdc.sharepoint.com/:x:/r/teams/NCEZID-OD_CAWG/_layouts/15/Doc.aspx?sourcedoc=%7B44002A5C-B5B1-49ED-AA6B-27F34EEAC8CC%7D&file=SC2_Proportion_Modeling_Run_Records.xlsx&action=default&mobileredirect=true&cid=13859ed2-691c-4865-b2a7-c548e4b1a585)

### Option B
This is the expanded way to run the whole process step by step without using the wrapper script. Use this option when debugging, or when special modifications/requests/tests are needed for the modeling run.

1. Navigate to the folder from which you will run the analyses. Typically this is the shared project folder. 
   ```bash
   cd /scicomp/groups/Projects/SARS2Seq/repos/sc2_proportion_modeling
   ```
2. Make sure the files are up to date with the git repository
   ```bash 
	git status 
	# make sure it is on the appropriate branch
	git fetch --all
	git status
   # if not up-to-date, pull updates from remote repository
   git pull
   ```
3. Update configuration settings in `config/config.R` if required. The variables you will most likely need to update include `data_date`, `voc1`, `voc2_additional`. Details about these variables can check [code](#code) and [Option A](#option-a).
4. Run the `variant_surveillance_system.sh` script to pull and clean the data. Currently two other different options can be chosen.   
   - One is to pull data with custom lineages. This requires defining the custom lineage in sql in` variant_suerveillance_system.R` and update the lineage aggregation code blocks in `weekly_variant_report_nowcast.R`. 
   - The other is to pull lineage definition from the `nexclade_pango` field from `sc2_src.nextclade` table. If using this option, lineage aggregation code blocks in `weekly_variant_report_nowcast.R` need to be manually updated.
   ```bash
   qsub variant_surveillance_system.sh <username> <password> <'T' or 'F' for custom lineage> <'T' or 'F' for nextclade_pango>
   ```
   This  will generate the survey dataset. Data results will be output in a folder titled `data/` and will be dated using `data_date` variable from `config/config.R`.

5. Wait until `variant_surveillance_system.sh` finishes which creates the databaset. Check `Run_var_sys.err` and `Run_var_sys.out`.

6. Run the `weekly_variant_report_nowcast.R` script with the desired specified runs. If special modifications are needed, for example using different weighting method, variable inputs in the run1_trim.sh or run2_trim.sh need to be modified. See details in the [code](#code) section.
   ```bash
   qsub run1_trim.sh 
   qsub run2_trim.sh
   ```
7. Check results in the result folder: 
   ```paste0("/scicomp/groups/Projects/SARS2Seq/repos/sc2_proportion_modeling/results_", data_date, '_', results_tag)```

   a. Check .err and .out files to make sure no running error occured.\
      The following error message can be ignored:\
          ```WARNING: sun.reflect.Reflection.getCallerClass is not supported. This will impact performance. ERROR StatusLogger No Log4j 2 configuration file found. Using default configuration (logging only errors to the console), or user programmatically provided configurations. Set system property 'log4j2.debug' to show Log4j 2 internal initialization logging. See https://logging.apache.org/log4j/2.x/manual/configuration.html for instructions on how to configure Log4j 2```

   b. Move the corresponding .err and .out files to result folder.\
      ```bash
      mv Run1_<results_tag>.[eo]* results/results_<data_date>_<results_tag>/
      ```
   
   c. automated aggregation is implemented based on the `extended_lineage` field in `sc2_src.pangolin`. This only works when custom lineage option is NOT used and nextclade_pango option is NOT used. lineage aggregation results can be checked in the following files:
      - `lineage_aggregation_summary...` list the aggregation done before the actuall weighting and modeling processes.
      - `agg_var_mat...` lists the aggregation done after the nowcast modeling processes are finished in Run2, to aggregate variants in voc2 but not in voc1 to their parental lineages in voc1, and generate the nowcast results for Run 1.

8. Backup the code that was used for a production run to the git repository
   ```bash
   git status
   git add config/config.R # and other files modified for the run
	git commit -m 'Production runs: YYYY-MM-DD' # and any other brief notes for changes made for the week
	git push
   ```

# Other Important Notes
## Lab Aggregation
Some submitting labs might input slightly different names for their sequences, so they should be aggregated/normalized to be the name in order to be calculated correctly in the survey design method. This is done in `variant_surveillance_system.R` around Line 1548 - 1095.

- Check the labs that are in the `Run_var_sys.out` after "\nnewly added lab names:"\
Compare that to the list of all lab names immediately after to see if any new labs should be combined with any previous labs.\
For example, the below two labs should be combined.\
	- "BUREAU OF LABORATORIES, PENNSYLVANIA DEPARTMENT OF HEALTH"
	- "PENNSYLVANIA DEPARTMENT OF HEALTH BUREAU OF LABORATORIES"
- every now and again it's good to look at the "data/backup_YYYY-MM-DD/lab_name_updates_YYYY-MM-DD.csv" file to see if any labs were combined that should NOT be combined. 
- Steps to add more lab aggregations
   1. Find lab names that almost assuredly refer to the same lab
   2. Copy-and-paste one of the blocks of code above
   3. Change "XX_labs_to_agg" and "labnames_df_xx" to new AND UNIQUE names (do a control-F for the new name to make sure it's unique)
   4. Change the regex pattern to something that will return ONLY your set of labs
   5. Add "labnames_df_xx" to "labnames_df" below
   6. Run variant_surveillance_system.sh again
   7. After running the new code, look at ./data/backup_YYYY-MM-DD/lab_name_updates_YYYY-MM-DD.csv to make sure that ONLY the intended labs are being renamed.

## Lineage Aggregation
## Weight Scalling

# Color Codes used on CDT
|     Variant            |     Hex Value    |   |
|------------------------|------------------|---|
|     B.1.1.529          |     #E26028      |   |
|     BA.1.1             |     #FF824C      |   |
|     BA.2               |     #9CCC65      |   |
|     BA.2.12.1          |     #7CB342      |   |
|     BA.2.75            |     #D4E157      |   |
|     BA.2.75.2          |     #C0CA33      |   |
|     CH.1.1             |     #827717      |   |
|     BN.1               |     #9e9d24      |   |
|     BA.4               |     #FFD54F      |   |
|     BA.4.6             |     #FFB300      |   |
|     BA.5               |     #80CBC4      |   |
|     BF.7               |     #81D4FA      |   |
|     BF.11              |     #29B6F6      |   |
|     BA.5.2.6           |     #009688      |   |
|     BQ.1               |     #006064      |   |
|     BQ.1.1             |     #00838F      |   |
|     XBB                |     #9FA8DA      |   |
|     XBB.1.5            |     #3F51B5      |   |
|     XBB.1.5.1          |     #1a237e      |   |
|     XBB.1.9.1          |     #304ffe      |   |
|     XBB.1.9.2          |     #536DFE      |   |
|     FD.2               |     #8C9EFF      |   |
|     XBB.1.16           |     #4527A0      |   |
|     XBB.1.16.1         |     #f89b9b      |   |
|     XBB.2.3            |     #f7bcdf      |   |
|     B.1.617.2          |     #B39DDB      |   |
|     Other              |     #797979      |   |

## Column Descriptions for Run1 and Run2 proportion modeling output (sc2_archive.state_proportion_modeling)
| name                 | type      | comment                                                                                                                                                                                 |
|----------------------|-----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| usa_or_hhsregion     | string    | Specifies the geographic area for the estimate; either national   (USA), or HHS regions 1 to 10.                                                                                        |
| week_ending          | timestamp | The date of the last day (Saturday) of the week (or 2-week period) for   which an estimate is made.                                                                                     |
| variant              | string    | The name of the variant that the estimate is for. (Note that some   variants are aggregates of all that variantâ€™s subvariants.)                                                       |
| share                | float     | The proportion of the variant.                                                                                                                                                          |
| share_lo             | float     | The upper bound of the 95% confidence interval for the estimated variant   share (prediction intervals for smoothed estimates).                                                         |
| share_hi             | float     | The lower bound of the 95% confidence interval for the estimated variant   share (prediction intervals for smoothed estimates).                                                         |
| count                | int       | The number of observed samples of the given variant in the given week and   region.                                                                                                     |
| denom_count          | int       | The total number of samples (of all variants) in the given week and   region.                                                                                                           |
| df                   | int       | The approximate degrees of freedom.                                                                                                                                                     |
| eff_size             | float     | The effective sample size.                                                                                                                                                              |
| ci_width             | float     | An NCHS QA for proportions: the absolute width of the confidence interval   (i.e. share_hi - share_lo)                                                                                  |
| nchs_flag            | tinyint   | Whether or not any of the NCHS QA checks identified the estimate as less   reliable.                                                                                                    |
| nchs_flag_wodf       | tinyint   | Whether or not any of the NCHS QA checks (other than the Degrees of   Freedom check) identified the estimate as less reliable.                                                          |
| count_lt20           | tinyint   | An NCHS QA for proportions: If the estimate is based on less than 20   samples.                                                                                                         |
| count_lt10           | tinyint   | An NCHS QA for proportions: If the estimate is based on less than 10   samples.                                                                                                         |
| modeltype            | string    | Estimate is based on weighted data ("weighted") or Nowcasted   data ("smoothed").                                                                                                       |
| interval             | string    | Whether data is base on one week interval ("weekly") or two   week intervals ("biweekly").                                                                                              |
| creation_date        | timestamp | Date when the estimate was created.                                                                                                                                                     |
| notes                | string    |                                                                                                                                                                                         |
| tagged_included      | tinyint   | Whether or not tagged samples were included in the estimate. Tagged   samples are those that labs other than the CDC contracting labs   self-identified as being surveillance quality.  |
| total_test_positives | int       | The total number of PCR positive tests from a given region in a given   week.                                                                                                           |
| cases                | float     | The estimated number of infections attributable to a given variant in a   given week. (cases = share * total_test_positives)                                                            |
| cases_lo             | float     | The lower bound of the 95% confidence interval for the estimated number   of infections attributable to a given variant in a given week. (cases_lo =   share_lo * total_test_positives) |
| cases_hi             | float     | The upper bound of the 95% confidence interval for the estimated number   of infections attributable to a given variant in a given week. (cases_hi =   share_hi * total_test_positives) |


## Column Descriptions for Run3 State proportionmodeling output (sc2_archive.state_proportion_modeling)
| name                 | type      | comment                                                                                                                                                                                 |   |   |
|----------------------|-----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---|---|
| state                | string    | Specifies the geographic area (abbreviation for state or jurisdiction)   for the estimate                                                                                               |   |   |
| roll_fourweek_ending | timestamp | The date of the last day (Saturday) of the 4-week time period for which   an estimate is made.                                                                                          |   |   |
| variant              | string    | The name of the variant that the estimate is for. (Note that some   variants are aggregates of all that variantâ€™s subvariants.)                                                       |   |   |
| share                | float     | The proportion of the variant.                                                                                                                                                          |   |   |
| share_lo             | float     | The upper bound of the 95% confidence interval for the estimated variant   share (prediction intervals for smoothed estimates).                                                         |   |   |
| share_hi             | float     | The lower bound of the 95% confidence interval for the estimated variant   share (prediction intervals for smoothed estimates).                                                         |   |   |
| count                | int       | The number of observed samples of the given variant in the given week and   region.                                                                                                     |   |   |
| denom_count          | int       | The total number of samples (of all variants) in the given week and   region.                                                                                                           |   |   |
| df                   | int       | The approximate degrees of freedom.                                                                                                                                                     |   |   |
| eff_size             | float     | The effective sample size.                                                                                                                                                              |   |   |
| ci_width             | float     | An NCHS QA for proportions: the absolute width of the confidence interval   (i.e. share_hi - share_lo)                                                                                  |   |   |
| nchs_flag            | int       | Whether or not any of the NCHS QA checks identified the estimate as less   reliable.                                                                                                    |   |   |
| nchs_flag_wodf       | int       | Whether or not any of the NCHS QA checks (other than the Degrees of   Freedom check) identified the estimate as less reliable.                                                          |   |   |
| creation_date        | timestamp |                                                                                                                                                                                         |   |   |
| note                 | string    | Date when the estimate was created.                                                                                                                                                     |   |   |
| total_test_positives | int       | The total number of PCR positive tests from a given region in a given   week.                                                                                                           |   |   |
| cases                | float     | The estimated number of infections attributable to a given variant in a   given week. (cases = share * total_test_positives)                                                            |   |   |
| cases_lo             | float     | The lower bound of the 95% confidence interval for the estimated number   of infections attributable to a given variant in a given week. (cases_lo =   share_lo * total_test_positives) |   |   |
| cases_hi             | float     | The upper bound of the 95% confidence interval for the estimated number   of infections attributable to a given variant in a given week. (cases_hi =   share_hi * total_test_positives) |   |   |


## Multinomial Logistic Regression growth rates
Multinomial logistic regression is essentially a series of logistic regressions (comparing each outcome category to the same reference category) with an added [normalization factor](https://en.wikipedia.org/wiki/Multinomial_logistic_regression#As_a_log-linear_model:~:text=as%20well%20as%20an%20additional%20normalization%20factor) to ensure that the probabilities always sum to 1. The normalization factor changes depending on the values of predictor variables and coefficient values. Because the relationship between coefficient value(s) and the predicted proportion vary with the normalization factor, one cannot just use multinomial regression coefficient value(s) as growth rates (as one can do with logistic regression). Rather, we calculate growth rate as the exponent of the derivative (i.e. instantaneous rate of change) of the log of the model-estimated proportion at time $`t`$:  $`p_i(t)`$. The estimated proportion of variant $`i`$ at time $`t`$ is:  


```math 
p_i(t) = \frac{e^{b_{0i}+b_{1i}*t}}{ \sum_je^{b_{0j}+b_{1j}*t} } 
```

If we take the log of $p_i(t)$ (so that we’re on the [linear predictor scale](https://en.wikipedia.org/wiki/Multinomial_logistic_regression#As_a_log-linear_model:~:text=we%20model%20the%20logarithm%20of%20the%20probability) instead of the response scale): 

```math
log⁡(p_i(t))=log⁡(\frac{e^{b_{0i}+b_{1i}*t}}{ \sum_je^{b_{0j}+b_{1j}*t} })
```

we can then take the derivative of $log(p_i(t))$ to get the instantaneous rate of change:

```math
\frac{d log⁡(p_i(t))}{dt}=b_{1i}-\sum_j{p_j*b_{1j}}
```

and exponentiate to get our growth rate (then multiply to 100 so that it's a percent instead of a decimal and subtract 100 so that "no change" is 0 instead of (100% of the current value).)

```math 
growth\ rate=100*e^{b_{1i}-\sum_j{p_j*b_{1j}}}-100
```

The resulting growth rate is the relative amount that the proportion is expected to change in a given time period. For example, a growth rate of 100 means that (if the growth rate were maintained) the proportion would grow by 100% over the given time period (i.e. double). However, these growth rates are _not_ constant over time (they decline with time).
