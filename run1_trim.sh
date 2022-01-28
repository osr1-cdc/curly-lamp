#!/bin/bash -l
#
# Embedded Grid Engine parameters
# All Embedded parameters must start with "#$"
# This is the same as adding these lines to the actual qsub line
#
# save the standard output text to this file instead of the the default jobID.o file
#$ -o Run1_trim.out
#
# save the standard error text to this file instead of the the default jobID.e file
#$ -e Run1_trim.err
# 
# Rename the job to be this string instead of the default which is the name of the script
# -N Run1_trim_proportion_modeling
#$ -N s3_r1
# 
# Refer all file reference to work the current working directory which is
# the directory from which the script was qsubbed
#$ -cwd
#
# Set up mail address for script
# -M ncy6@cdc.gov,fep2@cdc.gov,nyy7@cdc.gov,qiu5@cdc.gov,rsv4@cdc.gov,oow9@cdc.gov
#$ -M rsv4@cdc.gov
# ncy6 = Norman Hassell
# fep2 = Clinton Paden
# nyy7 = Sandra Seby
# qiu5 = Sherry (Xhiaoyu) Zheng
# rsv4 = Philip Shirk
# oow9 = Roopa Nagilla
#
# Mail options (b = beginning; e = end; a = aborted)
#$ -m ea
# 
# Always run in the default queue
#$ -q all.q
#
# Set the parallel_environment to "smp" and use 4 cores (smp = Symmetric multiprocessing or shared-memory multiprocessing)
#$ -pe smp 4

source /scicomp/groups-pure/Projects/SARS2Seq/bin/miniconda/bin/activate /scicomp/groups-pure/Projects/SARS2Seq/bin/miniconda/envs/prop_model-pure

Rscript weekly_variant_report_nowcast.R -r 1 -c F -v T -t quantile_99 -s T
# -r = run number
# -c = include custom lineages
# -v = use reduced vocs
# -t = trim weights
# -s = save datasets (results are always saved)