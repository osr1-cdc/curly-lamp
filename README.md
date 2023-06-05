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
      2) Run 2 generates the same 2 files as Run 1, but with vocs specified for Run 2 (voc2). Run 2 voc2s are usually automatically generated based on preset selection criteria. Currently, Run 2 autoselected vocs include all variants that occur at a frequency >= 1% of the unweighted sequence data in any of the 2-week periods from -1 to -7 2-week periods and variants at a frequency >= 0.5% of the unweighted sequence data in the -1 2-week period (See table below). All variants from Run 1 are also added to Run 2 voc2. Additionally, any variant of special interest can also be added to voc2 using the `voc2_additional` parameter. Run 2 fits model-based smoothed trends in variant share to both national and HHS regional estimates (i.e. "_Nowcast_"). Run 2 also generates nowcast estimates for voc list for Run1 (voc1) by aggregating estimates for lineages from Run 2 to their corresponding parental lineages.
         ### Run 2 vocs auto selection criteria

         | fortnight |      -7     |      -6     |      -5     |     -4    |     -3    |     -2    |      -1     |  Current Fortnight  |
         |:---------:|:-----------:|:-----------:|:-----------:|:---------:|:---------:|:---------:|:-----------:|      :-------------------:|
         |    week   | -15 and -14 | -13 and -12 | -11 and -10 | -9 and -8 | -7 and -6 | -5 and -4 |  -3 and -2  | -1 and Current Week |
         | Nowcast   |  >=1% (unw)  |  >=1% (unw)  |  >=1% (unw)  | >=1% (unw) | >=1% (unw) | >=1% (unw) | >=0.5% (unw) |                     |
      3) Run3 generates state-level estimates in rolling 4 wk bins using survey design (same as Run 1). NO LONGER NEEDED. 
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
   1) Run 1
      - Fornightly weighted estimates: `paste0("results/variant_share_weighted_",       ci.type,"CI_",svy.type,"_",data_date,tag,".csv")`
      - Weekly weighted estimates: `paste0("results/variant_share_weekly_weighted_",ci.type,"CI_",svy.type,"_",data_date,tag,".csv")`
   2) Run 2 
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
   3) Run3 generates state-level estimates in rolling 4 wk bins using survey design (same as Run 1). NO LONGER NEEDED. 
      - Output includes 1 file
         - `paste0("/results/state_weighted_roll4wk_",ci.type,"CI_svyNEW_",data_date,tag,".csv")`


## How to run from a scicomp location
It is best to run this from rosalind (`rosalind.biotech.cdc.gov`). 

# Option A
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
   qsub proportion_modeling_run.sh -u <username> -p <cdp password> -c <T/F for custom lineage>
   ```
   
   This wrapper script submits the job for running variant_surveillance_system.sh first. When that is finished, it submits the jobs for run1_trim.sh and run2_trim.sh. So user will receive individual notification emails for the above three jobs and finally a notification for finishing this wrapper job.

7. Check the wrapper job logfile to make sure all steps run successfully. The log file is `/scicomp/groups/Projects/SARS2Seq/repos/sc2_proportion_modeling/log`

8. Check results in the result folder: 
   ```paste0("/scicomp/groups/Projects/SARS2Seq/repos/sc2_proportion_modeling/results_", data_date, '_', results_tag)```

   a. Check .err and .out files to make sure no running error occured.\
      The following error message can be ignored:\
          ```WARNING: sun.reflect.Reflection.getCallerClass is not supported. This will impact performance. ERROR StatusLogger No Log4j 2 configuration file found. Using default configuration (logging only errors to the console), or user programmatically provided configurations. Set system property 'log4j2.debug' to show Log4j 2 internal initialization logging. See https://logging.apache.org/log4j/2.x/manual/configuration.html for instructions on how to configure Log4j 2```

   b. automated aggregation is implemented based on the `extended_lineage` field in `sc2_src.pangolin`. This only works when custom lineage option is NOT used and nextclade_pango option is NOT used. lineage aggregation results can be checked in the following files:
      - `lineage_aggregation_summary...` list the aggregation done before the actuall weighting and modeling processes.
      - `agg_var_mat...` lists the aggregation done after the nowcast modeling processes are finished in Run2, to aggregate variants in voc2 but not in voc1 to their parental lineages in voc1, and generate the nowcast results for Run 1.
   
9. Backup the code that was used for a production run to the git repository
   ```bash
   git add *
	git commit -m 'Production runs: YYYY-MM-DD'
	git push
   ```
10. Input run information in [SC2_Proportion_Modeling_Run_Records](https://cdc.sharepoint.com/:x:/r/teams/NCEZID-OD_CAWG/_layouts/15/Doc.aspx?sourcedoc=%7B44002A5C-B5B1-49ED-AA6B-27F34EEAC8CC%7D&file=SC2_Proportion_Modeling_Run_Records.xlsx&action=default&mobileredirect=true&cid=13859ed2-691c-4865-b2a7-c548e4b1a585)

# Option B
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
3. Update configuration settings in `config/config.R` if required. The variables you will most likely need to update include `data_date`, `voc1`, `voc2_additional`, `voc3`. However, updating configuration settings from one week to another is not necessarily required (as it was in the past).
4. Run the `variant_surveillance_system.sh` script submitting your username and password (FreeIPA CDP password) as well as specifying whether or not to include custom lineages. 
   ```bash
   qsub variant_surveillance_system.sh <username> <password> 'F' 'F'
   # optionally get custom dataset simultaneously from another terminal window or screen session
   qsub variant_surveillance_system.sh <username> <password> 'T' 'F'
   ```
   This  will generate the survey dataset. Data results will be output in a folder titled `data/` and will be dated using `data_date` variable from `config/config.R`.
5. Wait until databasets are created.
   - As the code runs, it will print out lab names to the terminal window. `LAB` includes lab source names as present in the dataset. `LAB2` includes cleaned lab source names. Look for typos in `LAB2` that result in individual labs being listed more than once (e.g. `Montana Public Health Lab` and `Montana PHL`). If `LAB2` contains duplicates, edit `variant_surveillance_system.R` to combine them. 
6. Run the `weekly_variant_report_nowcast.R` script with the desired specified runs (1-3) and specifying whether or not to include custom lineages. 
   ```bash
   qsub run1_trim.sh 
   qsub run2_trim.sh
   ```
7. Backup the code that was used for a production run to the git repository
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
