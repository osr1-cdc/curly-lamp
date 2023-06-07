#!/bin/bash -l

# Title: proportion_modeling_run.sh
# Usage: proportion_modeling_run.sh <username> <password>
# Version: v1.1 by Xiaoyu Sherry Zheng (qiu5@cdc.gov)
# Date: 05/16/2023
# Description: This wrapper script is created to automatically execute the regular proportion_modeling runs.
#              This script currently automatically start a regular modeling runs (run1, run2 and run3) using settings that are already set in config/config.R.


# Embedded Grid Engine parameters
# All Embedded parameters must start with "#$"
# This is the same as adding these lines to the actual qsub line
#
# save the standard output text to this file instead of the the default jobID.o file
#$ -o modeling_run.out
#
# save the standard error text to this file instead of the the default jobID.e file
#$ -e modeling_run.err
# 
# Rename the job to be this string instead of the default which is the name of the script
#$ -N modeling_run
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
#$ -pe smp 1
# Set the amount of RAM (per processor) to use (default is 32 GB)
#$ -l h_vmem=32G
# set the run-time <hh:mm:ss> (default is 72 hrs)
#$ -l h_rt=08:00:00
##############################################################
# HELPER FUNCTIONS
##############################################################
#sourced from Norman Hassell's nettree script
function time_stamp() {
	local t=$(date +"%Y-%m-%d %k:%M:%S")
	echo -e "[$t] $1"
}

##############################################################
# MAIN
##############################################################
# Run all scripts in the sc2_proportion_modeling folder
# cd /scicomp/groups-pure/Projects/SARS2Seq/repos/sc2_proportion_modeling

#---------------------
# Defining variables
#---------------------
# Identify running directory and set up log file.
bpath=$(pwd)
LOGFILE=${bpath}/log

# Get user-specified argument values
while getopts "u:p:c:h" opt; do
    case $opt in
        u)  username="$OPTARG"
            ;;
        p)  password="$OPTARG"
            ;;
        c)  custom="$OPTARG"
            ;;
        h)  echo "$usage"
            exit
            ;;
        \?) echo "Invalid option: -$OPTARG" >> ${LOGFILE}
            exit
            ;;
    esac
done

# If no user-supplied username and password, print error msg and usage, exit
if [[ (-z ${username}) || (-z ${password}) ]]; then
    echo "Data pull failed because username and password have to be provided." >> ${LOGFILE}
    exit
fi

#-------------------------------------------
# Step 1: Pull data in
#-------------------------------------------
time_stamp "Modeling run started." >> ${LOGFILE}
# submit variant_surveillance_system.sh as hpc job
qsub -N variant_surveillance -sync y ./variant_surveillance_system.sh "${username}" "${password}" "${custom}" 'F'

#-------------------------------------------
# Step 2: Start modeling runs
#-------------------------------------------
if [[ $? -eq 0 ]] && [[ ${custom} == 'F' ]]; then
    time_stamp "Data has been pulled in." >> ${LOGFILE}
    qsub -hold_jid variant_surveillance ./run1_trim.sh
    qsub -hold_jid variant_surveillance -sync y ./run2_trim.sh
#    qsub -hold_jid variant_surveillance ./run3_trim.sh
elif [[ $? -eq 0 ]] && [[ ${custom} == 'T' ]]; then
    time_stamp "Data has been pulled in." >> ${LOGFILE}
    qsub -hold_jid variant_surveillance ./run1_trim_custom.sh
    qsub -hold_jid variant_surveillance -sync y ./run2_trim_custom.sh
#    qsub -hold_jid variant_surveillance ./run3_trim_custom.sh
else
    time_stamp "Data pull failed. Check Run_var_sys.err" >> ${LOGFILE}
    exit 1
fi
#-------------------------------------------
# Step 3: Send finishing notifications
#-------------------------------------------
if [[ $? -eq 1 ]] ; then
    time_stamp "Weekly_variant_report_nowcast.R errored out. Check the .err files for Run 1 and Run 2." >> ${LOGFILE}
    exti 1
else
    time_stamp "Modeling run finished." >> ${LOGFILE}
fi
