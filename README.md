# SARS-CoV-2 Proportion Modeling and Nowcast

Code for applying a weights and confidence intervals to SARS-CoV-2 proportion data using a survey design-based approach.

 **Written by** Prabasaj Paul <vig5@cdc.gov>, Molly Steele <xzn9@cdc.gov>

**Reference**:
Paul P, France AM, Aoki Y, et al. Genomic Surveillance for SARS-CoV-2 Variants Circulating in the United States, December 2020–May 2021. MMWR Morb Mortal Wkly Rep 2021;70:846–850. DOI: [http://dx.doi.org/10.15585/mmwr.mm7023a3](http://dx.doi.org/10.15585/mmwr.mm7023a3)


## Code
- **variant_surveillance_system.R** - Generates the analytic dataset with the survey weights
-  **weekly_variant_report_nowcast.R** - Reads in the analytic dataset, drops and observations where the sequencing source or the lineage/sequence information is unknown, and subsets the data depending on whether runs include/exclude the state tagged data

## Required Data
- analytics_metadata - Genomics data compiled by SSEV Bioinformatics
- HHS protect RT-PCR data - Testing dataset downloaded from HHS Protect. Example provided in Resources : *Tests_ByCol_Date_Weekly.csv*.
- Population Data - State population data. 2018 data in Resources: *ACStable_B01001_40_2012_5.txt*

## How to run from a scicomp location

Before running the scripts, variables for running must be updated in the `config/config.R` file.

1. Activate the conda environment (requires conda to be installed for your user):
   ```bash
   conda activate /scicomp/groups/Projects/SARS2Seq/bin/miniconda/envs/prop_model
   ```
2. Run the `variant_surveillance_system.R` script submitting your username and password (FreeIPA CDP password):
   ```bash
   Rscript variant_surveillance_system.R -u <username> -p <password>
   ```
   This  will generate the survey dataset. Data results will be output in a folder titled `data/` and will be dated.
3. Run the `weekly_variant_report_nowcast.R` script with the desired specified run type (run 1-3):
   ```bash
   Rscript weekly_variant_report_nowcast.R -r 1
   ```

If multiple runs need to be submitted (which is always), script runs can be submitted through cluster job 
submissions or background jobs.

