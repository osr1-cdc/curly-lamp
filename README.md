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

## How to run
*This is a work in progress...*

