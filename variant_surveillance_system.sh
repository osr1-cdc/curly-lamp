#!/bin/bash
#
# Embedded Grid Engine parameters
# All Embedded parameters must start with "#$"
# This is the same as adding these lines to the actual qsub line
#
# save the standard output text to this file instead of the the default jobID.o file
#$ -o Run_var_sys.out
#
# save the standard error text to this file instead of the the default jobID.e file
#$ -e Run_var_sys.err
# 
# Rename the job to be this string instead of the default which is the name of the script
# -N Run1_trim_proportion_modeling
#$ -N Run_var_sys
# 
# Refer all file reference to work the current working directory which is
# the directory from which the script was qsubbed
#$ -cwd
#
# Set up mail address for script
# -M ncy6@cdc.gov,fep2@cdc.gov,nyy7@cdc.gov,qiu5@cdc.gov,rsv4@cdc.gov,oow9@cdc.gov
#$ -M ukc2@cdc.gov,qiu5@cdc.gov
# ncy6 = Norman Hassell
# fep2 = Clinton Paden
# nyy7 = Sandra Seby
# qiu5 = Sherry (Xiaoyu) Zheng
# rsv4 = Philip Shirk
# oow9 = Roopa Nagilla
#
# Mail options (b = beginning; e = end; a = aborted)
#$ -m ea
# 
# Choose queue
#$ -q all.q
# -q covid.q
#
# Set the parallel_environment to "smp" and use xx cores (smp = Symmetric multiprocessing or shared-memory multiprocessing); MAKE SURE THIS IS >= p CORES!
#$ -pe smp 2
# Set the amount of RAM (per processor) to use (default is 32 GB)
#$ -l h_vmem=64G
# set the run-time <hh:mm:ss> (default is 72 hrs)
#$ -l h_rt=01:00:00
source /etc/profile

source /scicomp/groups-pure/Projects/SARS2Seq/bin/miniconda/bin/activate /scicomp/groups-pure/Projects/SARS2Seq/bin/miniconda/envs/prop_model-pure












Rscript variant_surveillance_system.R -u ${1} -p ${2} -c ${3} -n ${4}
# -u = username
# -p = password
# -c = custom_lineages
# -n = nextclade_pango
