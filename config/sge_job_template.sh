#!/bin/bash
# SC2 Proportion Modeling - SGE Job Template
#
# This script submits the SC2 proportion modeling pipeline to the SGE scheduler.
# Modify paths and settings as needed for your HPC environment.
#
# Usage:
#   qsub sc2_monthly_job.sh
#   qsub -v CONTAINER_PATH=/path/to/sc2.sif sc2_monthly_job.sh
#

#$ -o /path/to/logs/sc2-modeling.out
#$ -e /path/to/logs/sc2-modeling.err
#$ -N sc2-proportion-modeling
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=08:00:00
#$ -M your_email@cdc.gov
#$ -m eas

# Set error handling
set -euo pipefail

# Email notification (optional)
echo "Starting SC2 Proportion Modeling Pipeline at $(date)"
echo "Job ID: $JOB_ID"
echo "Host: $(hostname)"

# Configuration
CONTAINER_PATH="${CONTAINER_PATH:-/path/to/sc2-proportion-modeling.sif}"
DATA_DIR="${DATA_DIR:-/path/to/data}"
RESULTS_DIR="${RESULTS_DIR:-$DATA_DIR/results}"
CACHE_DIR="${CACHE_DIR:-$DATA_DIR/cache}"
CONFIG_FILE="${CONFIG_DIR:-$DATA_DIR/config.yml}"

# Verify container exists
if [ ! -f "$CONTAINER_PATH" ]; then
    echo "ERROR: Container not found at $CONTAINER_PATH" >&2
    exit 1
fi

# Create output directories
mkdir -p "$RESULTS_DIR" "$CACHE_DIR"

# Load Singularity module (or Apptainer if available)
module load singularity || module load apptainer || true

# Verify Singularity/Apptainer is available
command -v singularity > /dev/null || command -v apptainer > /dev/null || {
    echo "ERROR: Singularity or Apptainer not found in PATH" >&2
    exit 1
}

# Run the pipeline
echo "Running: singularity run"
echo "  --bind $DATA_DIR:/app/data"
echo "  --bind $RESULTS_DIR:/app/results"
echo "  --bind $CACHE_DIR:/app/cache"
echo "  $CONTAINER_PATH"
echo "  --config /app/config/config.yml"
echo ""

singularity run \
    --bind "$DATA_DIR:/app/data" \
    --bind "$RESULTS_DIR:/app/results" \
    --bind "$CACHE_DIR:/app/cache" \
    "$CONTAINER_PATH" \
    --config /app/data/config.yml

EXITCODE=$?

if [ $EXITCODE -eq 0 ]; then
    echo "✅ Pipeline completed successfully"
    echo "Results available in: $RESULTS_DIR"
else
    echo "❌ Pipeline failed with exit code: $EXITCODE"
fi

exit $EXITCODE
