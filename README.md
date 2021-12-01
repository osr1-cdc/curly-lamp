# SARS-CoV-2 Proportion Modeling and Nowcast

Code for applying a weights and confidence intervals to SARS-CoV-2 proportion data using a survey design-based approach.

 **Written by** Prabasaj Paul <vig5@cdc.gov>, Molly Steele <xzn9@cdc.gov>
 **with updates from** Norman Hassell <ncy6@cdc.gov>, Philip Shirk <rsv4@cdc.gov>

**Reference**:
Paul P, France AM, Aoki Y, et al. Genomic Surveillance for SARS-CoV-2 Variants Circulating in the United States, December 2020–May 2021. MMWR Morb Mortal Wkly Rep 2021;70:846–850. DOI: [http://dx.doi.org/10.15585/mmwr.mm7023a3](http://dx.doi.org/10.15585/mmwr.mm7023a3)


## Code
- **config/config.R** - Specify various configuration settings, such as data dates, VOCs to include in particular runs, and figure settings.
- **variant_surveillance_system.R** - Generates the analytic dataset with the survey weights
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
   # optionally get custom dataset simultaneously from another terminal window
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