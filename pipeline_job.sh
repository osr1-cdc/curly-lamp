#!/bin/bash
#
# Generic Grid Engine wrapper for pipeline.R subcommands.

# save the standard output text to this file instead of the the default jobID.o file
#$ -o pipeline_job.out

# save the standard error text to this file instead of the the default jobID.e file
#$ -e pipeline_job.err

# Rename the job to be this string instead of the default which is the name of the script
#$ -N pipeline_job

# Refer all file reference to work the current working directory which is the directory from which the script was qsubbed
#$ -cwd

# Set up mail address for script
#$ -M ukc2@cdc.gov

# Mail options (b = beginning; e = end; a = aborted)
#$ -m ea

# Choose queue
#$ -q all.q

source /etc/profile

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/config/conda_env.sh"
activate_prop_model_env || exit 1

Rscript pipeline.R "$@"
