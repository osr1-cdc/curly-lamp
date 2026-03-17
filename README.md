# SARS-CoV-2 Proportion Modeling and Nowcast

**Standard Production Code for CDC Monthly Reporting**

Code for applying weights and confidence intervals to SARS-CoV-2 proportion data using a survey design-based approach.

 **Written by** Prabasaj Paul <vig5@cdc.gov>, Molly Steele <xzn9@cdc.gov>
 **with updates from** Norman Hassell <ncy6@cdc.gov>, Philip Shirk <rsv4@cdc.gov>, Xiaoyu Sherry Zheng <qiu5@cdc.gov>
 **recent contributors (past year)** Clint Paden <fep2@cdc.gov>, Juan Castro <jcb0@cdc.gov>, Peter Cook <ooj4@cdc.gov>, Casey Smith osr1 <osr1@cdc.gov>, Megha Aggarwal <tlz5@cdc.gov>

**Reference**:
Paul P, France AM, Aoki Y, et al. Genomic Surveillance for SARS-CoV-2 Variants Circulating in the United States, December 2020–May 2021. MMWR Morb Mortal Wkly Rep 2021;70:846–850. DOI: [http://dx.doi.org/10.15585/mmwr.mm7023a3](http://dx.doi.org/10.15585/mmwr.mm7023a3)

- [Code](#code)
- [Data Requirements](#data-requirements)
- [Other Requirements](#other-requirements)
- [Output Files](#output-files)
- [How to Run Monthly Reporting](#how-to-run-monthly-reporting)
- [Lab Aggregation](#lab-aggregation)
- [Column Descriptions](#column-descriptions)
- [Load Hadoop](#load-hadoop)
- [Implementation Notes](#implementation-notes)
- [Methods](#methods)

## Code

- **config/config.R** - Specify configuration settings: data dates, VOCs (variants shown on CDC COVID Data Tracker), and figure settings.
- **config/environment.yml** - Conda environment definition for the R and system dependencies used by the pipeline.
- **config/conda_env.sh** - Shared shell activation settings used by the HPC wrappers.

- **pipeline.R** - Main entrypoint for orchestration.
  - Commands: `prepare-data`, `run1`, `run2`, `run-all`, `submit-all`
  - Usage: `Rscript pipeline.R --help`

- **pipeline_job.sh** - Generic Grid Engine wrapper that runs a single `pipeline.R` subcommand inside the configured conda environment.

- **proportion_modeling_run.sh** - Thin submission wrapper for the full production workflow.
  - Usage: `qsub proportion_modeling_run.sh -u <username> -p <password>`
  - This wrapper calls `pipeline.R submit-all`

- **variant_surveillance_system.R** - Generates the analytic dataset with survey weights from CDC databases.
  - Input: CDC database queries for sequences, lineages, and testing data
  - Output: Processed data with survey design weights

- **weekly_variant_report_nowcast.R** - Creates variant proportion estimates using the dataset created in `variant_surveillance_system.R`.
  - Run 1: Calculates variant share/proportion and confidence intervals for fortnights and weeks (HHS regions & nationally)
  - Run 2: Generates the same output as Run 1, plus model-based smoothed trends (nowcasts) for both national and HHS regional estimates

- **weekly_variant_report_functions.R** - Library of helper functions used by `weekly_variant_report_nowcast.R`

- **pango_tree_diagram.R** - Optional utility for generating Pango lineage tree visualizations.
- **Lineage_Colors_on_CDT.csv** - Optional color map used by `pango_tree_diagram.R`.

## Data Requirements

- Sequence data (compiled by SSEV Bioinformatics infrastructure team)
  - Source: CDP database `sc2_archive.analytics_metadata_frozen`
- Pangolin lineages
  - Source: CDP database `sc2_src.pangolin`
- NREVSS testing data
  - Source: CDP database `sc2_archive.nrevss_frozen`
- HHS Protect RT-PCR testing data
  - Source: CDP database `sc2_archive.hhs_protect_testing_frozen`
- Population Data - State population data (as of 2018)
  - Source: `./resources/ACStable_B01001_40_2018_5.txt`

## Other Requirements

- Cloudera Impala JDBC driver: configured in [config.R](/scicomp/home-pure/osr1/repos/sc2_proportion_modeling/config/config.R)
- Conda environment: [environment.yml](/scicomp/home-pure/osr1/repos/sc2_proportion_modeling/config/environment.yml)
- Shared shell activation settings: [conda_env.sh](/scicomp/home-pure/osr1/repos/sc2_proportion_modeling/config/conda_env.sh)

## Output Files

Result folder: `results_<data_date>_<results_tag>`

### Run 1 Outputs
- Pre-modeling lineage aggregation results: `voc_aggregation_table_<date><tag>.csv`
- Lineage aggregation summaries: `lineage_aggregations_*.csv`, `lineage_aggregations_summary_*.csv`
- Fortnightly empirical estimates: `variant_share_<weighted_method>_<ci.type>CI_<svy.type>_<date><tag>_<results_tag>.csv`
- Weekly empirical estimates: `variant_share_weekly_<weighted_method>_<ci.type>CI_<svy.type>_<date><tag>_<results_tag>.csv`
- Hadoop-ready versions of the empirical outputs: matching `*_hadoop.csv` files

### Run 2 Outputs
Run 2 produces both empirical outputs and smoothed nowcast outputs.

- Pre-modeling lineage aggregation results: `voc_aggregation_table_<date><tag>.csv`
- Post-modeling lineage aggregation results: `agg_var_mat_<ci.type>CI_<svy.type>_<date><tag>.csv`
- Fortnightly nowcast estimates: `updated_nowcast_fortnightly_<weighted_method>_<date><tag>_<results_tag>.csv`
- Weekly nowcast estimates: `updated_nowcast_weekly_<weighted_method>_<date><tag>_<results_tag>.csv`
- Daily nowcast estimates: `updated_nowcast_weekly_<weighted_method>_<date><tag>_<results_tag>_daily.csv`
- Hadoop-ready versions of the nowcast outputs: matching `*_hadoop.csv` files
- Diagnostic/model artifacts such as `svymlm_*.RDS`, `src.moddat_*.RDS`, `src.dat_*.RDS`, and `wow_growth_variant_share_*.csv`

### Visualization Files (Run 2 only)
National data:
- `wtd_shares_<data_date>_barplot_US<tag>.png`
- `wtd_shares_<data_date>_projection_US<tag>.jpg`
- `wtd_shares_<data_date>_growthrate_US<tag>.png`

For each HHS region (1-10):
- `wtd_shares_<data_date>_barplot_HHS<region><tag>.jpg`
- `wtd_shares_<data_date>_projection_HHS<region><tag>.jpg`
- `wtd_shares_<data_date>_growthrate_HHS<region><tag>.jpg`

## How to Run Monthly Reporting

Best run from rosalind (`rosalind.biotech.cdc.gov`).

### Step 1: Prepare Configuration
Navigate to the repository:
```bash
cd /path/to/sc2_proportion_modeling
```

Create or update the environment if needed:
```bash
conda env create -f config/environment.yml
```

Make sure files are up to date:
```bash
git status
git fetch --all
git status
git pull  # if not up-to-date
```

### Step 2: Update Configuration Settings
Edit `config/config.R` with:
- **`data_date`**: Should be one of the `date_frozen` from `sc2_archive.analytics_metadata_frozen`. Check with:
  ```sql
  select distinct date_frozen from sc2_archive.analytics_metadata_frozen
  order by date_frozen desc
  ```
- **`voc1`**: Variants to display on CDC COVID Data Tracker
- **`voc2_additional`**: Any additional variants to include in Run 2 beyond auto-selected ones
- **`results_tag`**: Label for this run's results (used in filenames and database notes)

### Step 3: Get Auto-Selected VOC2 List
Run the voc2_autoselection SQL in Hue to identify variants automatically selected for Run 2.
- Ensure all voc1 variants are included in voc2_additional
- Add any special interest variants to voc2_additional

### Step 4: Update Email Notifications
In `proportion_modeling_run.sh`, update the email address:
```bash
#$ -M <your email address>
```

If the conda installation or env prefix differs on your host, update `CONDA_ACTIVATE` and `CONDA_ENV_PREFIX` in `config/conda_env.sh`. The shell entrypoints now fail fast if either path is invalid.

### Step 5: Submit the Master Wrapper Script
```bash
qsub proportion_modeling_run.sh -u <username> -p <cdp_password>
```

This automatically executes:
1. `pipeline.R submit-all` - HPC orchestration
2. `pipeline_job.sh prepare-data` - Data preparation job
3. `pipeline_job.sh run1` - Run 1 job
4. `pipeline_job.sh run2` - Run 2 job

### Step 6: Monitor Progress
Check the repo log file: `./log`

Each script generates `.err` and `.out` files. The following warning can be ignored:
```
WARNING: sun.reflect.Reflection.getCallerClass is not supported. This will impact performance.
ERROR StatusLogger No Log4j 2 configuration file found. Using default configuration (logging only errors to the console)
```

### Step 7: Verify Results
Check results folder: `results_<data_date>_<results_tag>`
- Verify all expected output files are present
- Check .err and .out files for errors
- Move any relevant scheduler `.err` and `.out` files into the matching results folder if you want them archived with the run.

### Step 8: Archive Code
Commit configuration changes to git for record-keeping:
```bash
git status
git add config/config.R
git commit -m "Production run: YYYY-MM-DD"
git push
```

### Step 9: Record Run Information
Input run information in [SC2_Proportion_Modeling_Run_Records](https://cdc.sharepoint.com/:x:/r/teams/NCEZID-OD_CAWG/_layouts/15/Doc.aspx?sourcedoc=%7B44002A5C-B5B1-49ED-AA6B-27F34EEAC8CC%7D&file=SC2_Proportion_Modeling_Run_Records.xlsx&action=default&mobileredirect=true&cid=13859ed2-691c-4865-b2a7-c548e4b1a585)

## Lab Aggregation

Some submitting labs use slightly different names for their sequences, so they are aggregated/normalized in `variant_surveillance_system.R` during data preparation.

- Check the labs in the Run_var_sys.out file after "\nnewly added lab names:"
- Compare to the list of all lab names immediately after
- Combine similar lab names if they refer to the same lab
- Check `data/backup_YYYY-MM-DD/lab_name_updates_YYYY-MM-DD.csv` to verify aggregations

**To add more lab aggregations:**
1. Find lab names that should be combined
2. Add a new rule to the `lab_rename_rules` list in `variant_surveillance_system.R`
3. Use either `regex_lab_selector(...)` or a small custom selector function
4. Set the replacement `new_name`
5. Re-run `Rscript pipeline.R prepare-data -u <username> -p <password>`
6. Verify aggregations in `data/backup_YYYY-MM-DD/lab_name_updates_YYYY-MM-DD.csv`

## Column Descriptions

### Table: sc2_archive.proportion_modeling

Output from both Run 1 and Run 2 (rows where modeltype = 'weighted' and 'smoothed')

| name                 | type      | comment                                                                                                                                                                                 |
|----------------------|-----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| usa_or_hhsregion     | string    | Specifies the geographic area for the estimate; either national (USA), or HHS regions 1 to 10                                                                                        |
| week_ending          | timestamp | The date of the last day (Saturday) of the week (or 2-week period) for which an estimate is made                                                                                     |
| variant              | string    | The name of the variant that the estimate is for. (Note that some variants are aggregates of all that variant's subvariants.)                                                       |
| share                | float     | The proportion of the variant                                                                                                                                                         |
| share_lo             | float     | The lower bound of the 95% confidence interval for the estimated variant share (prediction intervals for smoothed estimates)                                                         |
| share_hi             | float     | The upper bound of the 95% confidence interval for the estimated variant share (prediction intervals for smoothed estimates)                                                         |
| count                | int       | The number of observed samples of the given variant in the given week and region                                                                                                     |
| denom_count          | int       | The total number of samples (of all variants) in the given week and region                                                                                                           |
| df                   | int       | The approximate degrees of freedom                                                                                                                                                    |
| eff_size             | float     | The effective sample size                                                                                                                                                             |
| ci_width             | float     | An NCHS QA for proportions: the absolute width of the confidence interval (share_hi - share_lo)                                                                                      |
| nchs_flag            | tinyint   | Whether or not any of the NCHS QA checks identified the estimate as less reliable                                                                                                    |
| nchs_flag_wodf       | tinyint   | Whether or not any of the NCHS QA checks (other than the Degrees of Freedom check) identified the estimate as less reliable                                                          |
| count_lt20           | tinyint   | An NCHS QA for proportions: If the estimate is based on less than 20 samples                                                                                                         |
| count_lt10           | tinyint   | An NCHS QA for proportions: If the estimate is based on less than 10 samples                                                                                                         |
| modeltype            | string    | Estimate type: "weighted" (Run 1/2 weighted data) or "smoothed" (Run 2 nowcasted data)                                                                                              |
| interval             | string    | Time interval: "weekly" (1 week) or "biweekly" (2 week)                                                                                                                               |
| creation_date        | timestamp | Date when the estimate was created                                                                                                                                                    |
| notes                | string    | Additional notes                                                                                                                                                                       |
| tagged_included      | tinyint   | Whether tagged samples were included. Tagged samples are those from labs other than CDC contracting labs that self-identified as surveillance quality                               |
| total_test_positives | int       | The total number of PCR positive tests from a given region in a given week                                                                                                           |
| cases                | float     | Estimated number of infections from a given variant: cases = share × total_test_positives                                                                                           |
| cases_lo             | float     | Lower bound of 95% CI for estimated infections: cases_lo = share_lo × total_test_positives                                                                                          |
| cases_hi             | float     | Upper bound of 95% CI for estimated infections: cases_hi = share_hi × total_test_positives                                                                                          |

## Load Hadoop

Once data is available, load the CSV files into [sc2_archive.proportion_modeling](https://cdp-01.biotech.cdc.gov:8889/hue/filebrowser/view=%2Fwarehouse%2Ftablespace%2Fexternal%2Fhive%2Fsc2_archive.db%2Fproportion_modeling) at `hdfs:///warehouse/tablespace/external/hive/sc2_archive.db/proportion_modeling` on CDP.

Load the following files using Hue or hput:
- `updated_nowcast_fortnightly_*_hadoop.csv`
- `updated_nowcast_weekly_*_hadoop.csv`
- `variant_share_*_hadoop.csv`
- `variant_share_weekly_*_hadoop.csv`
- (And all other `*_hadoop.csv` files)

**Important:** Refresh the table after loading:
```sql
REFRESH sc2_archive.proportion_modeling
```

Extract data by creation_date into Tableau workbooks.

**Downstream tables:** See https://git.biotech.cdc.gov/sars2seq/sc2_variant_proportion_socrata_update/-/blob/main/README.md for updates to publishing pipelines.

## Implementation Notes

### Pipeline Structure

- `pipeline.R` is the main entrypoint for orchestration and shared runtime setup.
- `variant_surveillance_system.R` prepares the analytic dataset.
- `weekly_variant_report_nowcast.R` produces weighted share outputs and Run 2 nowcasts.
- `weekly_variant_report_functions.R` contains the survey-design and multinomial helper functions used by the nowcast step.

### Data Preparation Model

The preparation step treats each retained record as a sequenced positive specimen tied to a collection date and reporting state.

- Event time is specimen collection date.
- Event location is the reporting state.
- Strata are U.S. states and territories.
- Cluster structure follows the sequencing source stream, such as NS3 or contractor labs.
- Inputs include sequence metadata, Pangolin lineage data, testing data, and state population data.

### Weighting Notes

The modeling step supports multiple weight constructions.

- `population` weights represent the population covered by each sequence within a state-week.
- `updated` weights use regional NREVSS positivity and CELR testing totals to proxy infections, then allocate that burden back to states by population.
- SGTF upsampling weights currently default to `1` in the nowcast script; any future SGTF-specific weighting logic should be documented here rather than as inline TODO blocks.
- Invalid weights are removed before analysis, and optional trimming can cap extreme weights.

### Run Structure

- Run 1 produces weighted share estimates and confidence intervals for weekly and fortnightly outputs.
- Run 2 produces the same weighted outputs plus smoothed nowcast estimates and figures.
- Run 2 also writes both Run 1-style and Run 2-style downstream output tables.

### Data Cleaning Rules

Before the survey design object is created, the preparation step drops or repairs records that would make the analysis unreliable.

- Keep U.S. human-host sequences only.
- Drop invalid state abbreviations.
- Drop duplicates.
- Drop unreasonable specimen collection dates.
- Drop invalid lab names and invalid variant names.
- Drop records with invalid derived weights.

### Comment Policy

High-level workflow notes and assumptions should live in this README rather than as large inline comment blocks in the R scripts.

# Methods

## Weights

We weight sequences based on their inverse probability of selection from all SARS-CoV-2 infections in the United States. We estimate total infections in each state-week combination using the methods of [Chiu & Ndeffo-Mbah](https://doi.org/10.1371/journal.pcbi.1009374), based on observed correlation between testing and total infections.

```math
I_{rw} \approx \sqrt{ \frac{p_{rw}}{t_{rw}} } * t_{rw+}
```

>$I_{rw}$ = number of infections in HHS region $r$ in week $w$
>$p_{rw}$ = population of states in region $r$ that reported testing in week $w$
>$t_{rw}$ = total tests conducted in region $r$ in week $w$
>$t_{rw+}$ = number of positive tests in region $r$ in week $w$

We handle multiple test data sources (NREVSS and CELR) by separating test positivity:

```math
I_{rw} \approx \sqrt{ \frac{t_{rw+,N}}{t_{rw,N}} * p_{rwC} * t_{rw+C}}
```

>$t_{rw,N}$ = total tests in region $r$ week $w$ from NREVSS
>$t_{rw+,N}$ = positive tests in region $r$ week $w$ from NREVSS
>$t_{rw+,C}$ = positive tests in region $r$ week $w$ from CELR
>$p_{rwC}$ = population of states in region $r$ reporting to CELR in week $w$

Regional infections are distributed to states by population:

```math
I_{sw} = \frac{p_{s}}{p_{r}} * I_{rw}
```

Final weights are infections per sequence:

```math
w_{sw} = \frac{I_{sw}}{s_{sw}}
```

>$s_{sw}$ = number of sequences from state $s$ in week $w$

All sequences from the same state and week have the same weight.

## Lineage Aggregation

Lineages are aggregated to their nearest parent lineage listed in the VOC list. The nearest parent is determined using the `expanded_lineage` field and finding the longest lineage in the VOC that is a perfect subset.

For example, `XBB.1.5.1` aggregates to `XBB.1.5` if `XBB.1.5` is in VOC but not `XBB.1.5.1`.

We validate parent-child relationships by checking that the parent lineage (with `.` appended) is a subset of the child. This ensures `BA.1` is not considered a parent of `BA.11`.

Aggregation is validated in output files:
- `voc_aggregation_table...` - Shows aggregation before weighting and modeling
- `agg_var_mat...` (Run 2 only) - Shows aggregation after nowcast modeling (cells with value 1 indicate aggregation)

## Multinomial Logistic Regression Prediction Intervals

Multinomial logistic regression models are fit using the `nnet` package. The Hessian matrix is needed for confidence intervals. When numerical under/overflow occurs during Hessian estimation, we rescale weights by factors specified in `config/config.R` (`rescale_model_weights_by`) until successful estimation occurs.

Scaling changes intermediate values used in fitting without changing the model fit itself.

Model uncertainty is adjusted using survey design, as samples from the same state, week, and lab are correlated and contribute less information than random samples. The adjusted uncertainty is based on effective sample size rather than actual sequence count, and is typically larger than unadjusted model uncertainty.

## Multinomial Logistic Regression Growth Rates

Growth rate is calculated as the exponent of the derivative of log-proportions at time $t$:

```math
p_i(t) = \frac{e^{b_{0i}+b_{1i}*t}}{ \sum_je^{b_{0j}+b_{1j}*t} }
```

Taking the derivative of the log-proportion:

```math
\frac{d log⁡(p_i(t))}{dt}=b_{1i}-\sum_j{p_j*b_{1j}}
```

And exponentiating (multiplying by 100) for percent change:

```math
growth\ rate=100*(e^{b_{1i}-\sum_j{p_j*b_{1j}}}-1)
```

This represents relative proportional change over the time period. A growth rate of 100 means the proportion would double over the period if the growth rate were maintained. These growth rates decline over time and are not constant.
