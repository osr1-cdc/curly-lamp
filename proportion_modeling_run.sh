#!/bin/bash
# SGE params
#$ -o modeling_run.out
#$ -e modeling_run.err
#$ -N modeling_run
#$ -cwd
#$ -M ukc2@cdc.gov
#$ -m ea
#$ -q all.q
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=08:00:00

source /etc/profile
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/config/conda_env.sh"
activate_prop_model_env || exit 1

Rscript pipeline.R submit-all "$@"
