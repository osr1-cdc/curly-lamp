# SARS-CoV-2 Proportion Modeling and Nowcast

Code for applying a weights and confidence intervals to SARS-CoV-2 proportion data using a survey design-based approach.

 **Written by** Prabasaj Paul <vig5@cdc.gov>, Molly Steele <xzn9@cdc.gov>
 **with updates from** Norman Hassell <ncy6@cdc.gov>, Philip Shirk <rsv4@cdc.gov>

**Reference**:
Paul P, France AM, Aoki Y, et al. Genomic Surveillance for SARS-CoV-2 Variants Circulating in the United States, December 2020–May 2021. MMWR Morb Mortal Wkly Rep 2021;70:846–850. DOI: [http://dx.doi.org/10.15585/mmwr.mm7023a3](http://dx.doi.org/10.15585/mmwr.mm7023a3)


## Code
- **config/config.R** - Specify various configuration settings, such as data dates, VOCs to include in particular runs, and figure settings.
- **variant_surveillance_system.R** - Generates the analytic dataset with the survey weights
- **variant_surveillance_system.sh** - Wrapper script for submitting variant_surveillance_system.R to run in HPC.
   - qsub variant_surveillance_system.sh <`username`> <`password`> <T/F for custom lineage> <T/F for using nextclade pango calls> <ref_lineage from geni.sc2_lineage_rep>
- **weekly_variant_report_nowcast.R** - Creates variant proportion estimates using the dataset created in `variant_surveillance_system.R`. 
   - `weekly_variant_report_nowcast.R` accomodates 3 different "runs", each of which has its own set of "vocs" (i.e. variants for which to calculate proportions)
      1) Run 1 calculates variant share/proportion and confidence intervals estimated using survey design for both fortnights (HHS regions & nationally) and weeks (HHS regions & nationally). 
         - Output includes 2 files (1 for fortnightly estimates; 1 for weekly estimates): 
            - `paste0("results/variant_share_weighted_",       ci.type,"CI_",svy.type,"_",data_date,tag,".csv")`
            - `paste0("results/variant_share_weekly_weighted_",ci.type,"CI_",svy.type,"_",data_date,tag,".csv")`
      2) Run 2 generates the same 2 files as Run 1, but with vocs specified for Run 2. Run 2 vocs generally include all variants that occur at a frequency >= 1% of the unweighted sequence data in recent weeks (i.e. 2 to 14 weeks prior to date of running code), but can also be manually altered. Additionally, Run 2 fits model-based smoothed trends in variant share to both national and HHS regional estimates (i.e. "_Nowcast_").
         - Output includes 7 csv & 33 image files. *NOTE! Run 2 produces results for both Run 2 AND Run 1!*
            - CSV files
               - Files matching those produced by Run 1: 
                  - `paste0("results/variant_share_weighted_",       ci.type,"CI_",svy.type,"_",data_date,tag,".csv")`
                  - `paste0("results/variant_share_weekly_weighted_",ci.type,"CI_",svy.type,"_",data_date,tag,".csv")`
               - "_Nowcast_" output:
                  - `paste0("/results/wow_growth_variant_share",data_date,tag,".csv")`
                  - `paste0("/results/updated_nowcast_fortnightly_",data_date,"_state_tag_included_Run1.csv")`
                  - `paste0("/results/updated_nowcast_fortnightly_",data_date,"_state_tag_included_Run2.csv")`
                  - `paste0("/results/updated_nowcast_weekly_",data_date,"_state_tag_included_Run1.csv")`
                  - `paste0("/results/updated_nowcast_weekly_",data_date,"_state_tag_included_Run2.csv")`
            - Image files
               - National data
                  - `paste0("/results/wtd_shares_",data_date,"_","barplot_US",tag,".jpg")`
                  - `paste0("/results/wtd_shares_",data_date,"_","projection_US",tag,".jpg")`
                  - `paste0("/results/wtd_shares_",data_date,"_","growthrate_US",tag,".png")`
               - For each HHS region
                  - `paste0("/results/wtd_shares_",data_date,"_","barplot_HHS",hhs,tag,".jpg")`
                  - `paste0("/results/wtd_shares_",data_date,"_","projection_HHS",hhs,tag,".jpg")`
                  - `paste0("/results/wtd_shares_",data_date,"_","growthrate_HHS", hhs,tag,".jpg")`
      3) Run3 generates state-level estimates in rolling 4 wk bins using survey design (same as Run 1). 
         - Output includes 1 file
            - `paste0("/results/state_weighted_roll4wk_",ci.type,"CI_svyNEW_",data_date,tag,".csv")`
- **weekly_s1_variant_report_nowcast.R** - Creates s1_species proportion estimates using the dataset created in `variant_surveillance_system.R`. 




## Data Requirements
- Sequence data (compiled by SSEV Bioinformatics) 
   - Source: CDP database `sc2_archive.analytics_metadata_frozen`
- Pangolin lineages 
   - Source: depends on dates & custom lineages, but may include CDP database tables `sc2_src.pangolin`, `sc2_archive.analytics_metadata_frozen`, and/or `sc2_src.alignments`
- S-gene mutations
   - Source: CDP database `sc2_dev.baselineseq`
- HHS protect RT-PCR testing data
   - Source: CDP database `sc2_archive.hhs_protect_testing_frozen`
- Common variants in recent weeks (for Run 2)
   - Source: CDP database `sc2_air.analytics_metadata`, `sc2_src.variant_definitions`, `sc2_air.analytics_lineage_corr`
   - (Sources may change to: `sc2_archive.analytics_metadata_frozen`, `sc2_src.hhs_regions`, `sc2_src.variant_definitions`)
- Population Data - State population data (as of 2018) 
   - Source: `./resources/ACStable_B01001_40_2018_5.txt`

## Other requirements
- Cloudera Impala JDBC driver: provided in `./jdbc/ClouderaImpalaJDBC-2.6.20.1024/ClouderaImpalaJDBC41-2.6.20.1024`

## Important Note
The current branch has the following custom defined lineages:
- AY.4.2+ - AY.4.2 with Y145H and A222V
- AY.35+ - AY.35 with E484Q
All other lineages (including AY.4.2 and AY.35) are from default pangolin calls.
In order to change the definitions of custom lineages, one much change both the vector `custom_lineage_names` in `config/config.R` and also the definitions of each custom lineage in the Pangolin table database query in `variant_surveillance_system.R`. 

## How to run from a scicomp location
It is best to run this from DVD-VM (`dvd-sars2seq-dev01.biotech.cdc.gov`), but can also be run from other scicomp locations. 

1. Navigate to the folder from which you will run the analyses. Typically this is the shared project folder. 
   ```bash
   cd /scicomp/groups/Projects/SARS2Seq/repos/sc2_proportion_modeling
   ```
2. Save the results of previous runs
   ```bash
   # do things interactively via GUI
   xdg-open .
	# Create a new subfolder for results files: results/YYY-MM-DD
   # Move previous week’s results files into subfolder as a backup
   ```
3. Make sure the files are up to date with the git repository
   ```bash 
	git status 
	# make sure it's on the appropriate branch
	git fetch --all
	git status
   # if not up-to-date, pull updates from remote repository
   git pull
   ```
4. Activate the conda environment (requires conda to be installed for your user):
   ```bash
   conda activate /scicomp/groups/Projects/SARS2Seq/bin/miniconda/envs/prop_model
   ```
5. Update configuration settings in `config/config.R` if required. The variables you will most likely need to update include `data_date`, `voc1`, `voc2_additional`, `voc3`. However, updating configuration settings from one week to another is not necessarily required (as it was in the past).
6. Run the `variant_surveillance_system.R` script submitting your username and password (FreeIPA CDP password) as well as specifying whether or not to include custom lineages. 
   ```bash
   Rscript variant_surveillance_system.R -u <username> -p <password> -c F
   # optionally get custom dataset simultaneously from another terminal window or screen session
   Rscript variant_surveillance_system.R -u <username> -p <password> -c T
   ```
   This  will generate the survey dataset. Data results will be output in a folder titled `data/` and will be dated using `data_date` variable from `config/config.R`.
7. Wait until databasets are created.
   - As the code runs, it will print out lab names to the terminal window. `LAB` includes lab source names as present in the dataset. `LAB2` includes cleaned lab source names. Look for typos in `LAB2` that result in individual labs being listed more than once (e.g. `Montana Public Health Lab` and `Montana PHL`). If `LAB2` contains duplicates, edit `variant_surveillance_system.R` to combine them. 
8. Run the `weekly_variant_report_nowcast.R` script with the desired specified runs (1-3) and specifying whether or not to include custom lineages. 
   ```bash
   qsub run1.sh 
   qsub run2.sh
	qsub run3.sh
   qsub run1_custom.sh
	qsub run2_custom.sh
   # or run on current node
   # Rscript weekly_variant_report_nowcast.R -r 1 -c F
   ```
9. Backup the code that was used for a production run to the git repository
   ```bash
   git add *
	git commit -m 'Production runs: YYYY-MM-DD'
	git push
   ```


## Color Codes used on CDT
| Variant	| Hex Code|
| -------|-------- |
| B.1.621	| #cee1ec |
| B.1.617.2	| #f28e2b |
| B.1.1.7	| #1e80a0 |
| P.1	| #83962a |
| AY.2	| #c14f22 |
| AY.1	| #a13703 |
| B.1.351	| #8e49b5 |
| B.1.526	| #eec9e5 |
| B.1.525	| #e78cc7 |
| B.1.617.1	| #55aaff |
| B.1.617.3	| #ff55ff |
| Other	| #152d44 |
| C.37	| #9eacb4 |
| B.1.1.529 / BA.1	| #88419d |
| BA.1.1	| #4d004b |
| BA.2	| #f1b6da |
| BA.2.12.1	| #e15759 |
| BA.4	| #86bcb6 |
| BA.5	| #499894 |

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
