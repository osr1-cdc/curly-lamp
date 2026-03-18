# HPC Deployment Guide: SC2 Proportion Modeling with Singularity

## Overview

This guide provides detailed instructions for deploying and running the SC2 Proportion Modeling pipeline on HPC environments using Singularity/Apptainer containers and SGE (Sun Grid Engine) job scheduling.

**Target Environment:**
- HPC Cluster with Singularity/Apptainer runtime
- SGE job scheduler
- Shared filesystem for data/results
- Module system (environment modules)

---

## Quick Start (5 minutes)

### Step 1: Download the Container

```bash
# Get the latest release
wget https://github.com/anthropics/sc2-proportion-modeling/releases/download/v2.0.0/sc2-proportion-modeling.sif

# Verify the download
ls -lh sc2-proportion-modeling.sif
# Expected: ~1.2GB
```

### Step 2: Prepare Data Directory

```bash
mkdir -p sc2-run-{DATE}
cd sc2-run-{DATE}

# Create subdirectories
mkdir -p config results cache logs

# Copy your configuration file
cp /path/to/your/config.yml config/config.yml

# Verify config
cat config/config.yml
```

### Step 3: Run the Pipeline

```bash
# Test with dry-run first
singularity run \
  --bind $(pwd)/config:/app/config \
  --bind $(pwd)/results:/app/results \
  --bind $(pwd)/cache:/app/cache \
  ./sc2-proportion-modeling.sif \
  --config /app/config/config.yml \
  --dry-run

# If successful, run full pipeline
singularity run \
  --bind $(pwd)/config:/app/config \
  --bind $(pwd)/results:/app/results \
  --bind $(pwd)/cache:/app/cache \
  ./sc2-proportion-modeling.sif \
  --config /app/config/config.yml
```

### Step 4: Check Results

```bash
ls -la results/

# View outputs
head -20 results/nowcasts.csv
cat results/metadata.json | python -m json.tool
```

---

## Detailed Setup

### Module Loading

On CDC/NIAID HPC systems, load required modules:

```bash
# Load Singularity/Apptainer
module load singularity
# or: module load apptainer

# Optional: Load R if you need interactive testing (not required for container)
module load R/4.2  # If available

# Verify container runtime
singularity --version
# or: apptainer --version
```

### Directory Structure

Create a standard working directory structure:

```
/path/to/sc2-analysis/
├── sc2-proportion-modeling.sif    (container image)
├── data/
│   ├── config.yml                 (your configuration)
│   └── cache/                     (intermediate results, auto-created)
├── results/                       (output directory)
├── logs/                          (execution logs)
└── scripts/                       (job submission scripts)
```

### Configuration File

Create/modify `config/config.yml`:

```yaml
# SC2 Proportion Modeling Configuration

run:
  data_date: "2026-03-18"           # Data cutoff date
  results_tag: "CDT"                # Results tag (e.g., CDT, Production)
  use_previously_imported_data: false

impala:
  host: "localhost"                 # Or your Impala server hostname
  port: 21050
  database: "covid"
  timeout_seconds: 300

database:
  connection_string: ""             # Optional: JDBC connection for direct query

variants:
  voc1:
    - "BA.1.1"
    - "BA.2"
    - "XBB.1.5"
    - "JN.1"

model:
  weeks_lookback: 8
  weeks_predict: 4
  nowcast_method: "multinomial_regression"
  credible_interval: 0.95
  ci_type: "KG"

output:
  results_dir: "/path/to/results"   # Absolute path or relative to working dir
  cache_dir: "/path/to/cache"       # For parquet caching
  export_formats:
    - "csv"
    - "json"
    - "parquet"
    - "png"
    - "pdf"
```

---

## SGE Job Submission

### Basic Job Script Template

Create `sc2_monthly_job.sh`:

```bash
#!/bin/bash
#$ -o /path/to/logs/modeling_%Y%m%d_%H%M%S.out
#$ -e /path/to/logs/modeling_%Y%m%d_%H%M%S.err
#$ -N sc2-proportion-modeling
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=08:00:00
#$ -M your_email@cdc.gov
#$ -m ea

# Load required modules
module load singularity

# Set working directory
WORK_DIR="/path/to/sc2-analysis"
cd $WORK_DIR

# Container path
CONTAINER="$WORK_DIR/sc2-proportion-modeling.sif"

# Run the pipeline
singularity run \
  --bind $WORK_DIR/config:/app/config \
  --bind $WORK_DIR/results:/app/results \
  --bind $WORK_DIR/cache:/app/cache \
  $CONTAINER \
  --config /app/config/config.yml \
  -v

# Check exit code
if [ $? -eq 0 ]; then
  echo "Pipeline completed successfully"
  # Optional: Send notification
  echo "SC2 pipeline completed" | mail -s "SC2 Complete" $USER@cdc.gov
else
  echo "Pipeline failed with status $?"
  echo "Check logs: $WORK_DIR/logs/"
  exit 1
fi
```

### Submit Job

```bash
# Make script executable
chmod +x sc2_monthly_job.sh

# Submit to queue
qsub sc2_monthly_job.sh

# Monitor job
qstat -u $USER

# View job details
qstat -j <job_id>

# Cancel job if needed
qdel <job_id>
```

### Advanced Job Script with Data Staging

```bash
#!/bin/bash
#$ -o /path/to/logs/modeling.out
#$ -e /path/to/logs/modeling.err
#$ -N sc2-modeling
#$ -pe smp 2
#$ -l h_vmem=48G
#$ -l h_rt=10:00:00

set -e  # Exit on error

module load singularity

WORK_DIR="/path/to/sc2-analysis"
SCRATCH_DIR="/scratch/$USER/sc2-$$"  # Temporary local storage

# Create scratch directory
mkdir -p $SCRATCH_DIR
trap "rm -rf $SCRATCH_DIR" EXIT  # Cleanup on exit

# Stage data to local scratch (faster I/O)
cp -r $WORK_DIR/{config,cache} $SCRATCH_DIR/

# Run container with scratch directories
singularity run \
  --bind $SCRATCH_DIR/config:/app/config \
  --bind $SCRATCH_DIR/results:/app/results \
  --bind $SCRATCH_DIR/cache:/app/cache \
  $WORK_DIR/sc2-proportion-modeling.sif \
  --config /app/config/config.yml

# Stage results back to shared storage
cp -r $SCRATCH_DIR/results/* $WORK_DIR/results/
cp -r $SCRATCH_DIR/cache/* $WORK_DIR/cache/

echo "Job completed successfully"
```

---

## Volume Binding (--bind Flags)

Singularity `--bind` maps host directories into the container:

```bash
singularity run \
  --bind /host/path:/container/path \
  --bind /another/host/path:/container/path2 \
  image.sif
```

### Standard Mounts for SC2

| Host Path | Container Path | Purpose | Notes |
|-----------|----------------|---------|-------|
| `./config` | `/app/config` | Configuration files | Must contain config.yml |
| `./results` | `/app/results` | Output directory | Results written here |
| `./cache` | `/app/cache` | Intermediate cache | Parquet files, Stan models |
| `./logs` | `/app/logs` | Execution logs | Optional, for loguru output |
| `/data/impala` | `/app/data` | External data (optional) | For mounted Impala tables |

### Example with Multiple Mounts

```bash
singularity run \
  --bind $(pwd)/config:/app/config:ro \
  --bind $(pwd)/results:/app/results:rw \
  --bind $(pwd)/cache:/app/cache:rw \
  --bind /mnt/shared/impala_data:/app/data:ro \
  sc2-proportion-modeling.sif \
  --config /app/config/config.yml
```

**Mount Options:**
- `:ro` - Read-only (recommended for config)
- `:rw` - Read-write (default for results/cache)

---

## Runtime Options

### Common Command-Line Flags

```bash
# Show help
singularity run sc2-proportion-modeling.sif --help

# Dry-run (validate config without executing)
singularity run ... sc2-proportion-modeling.sif --dry-run

# Verbose logging
singularity run ... sc2-proportion-modeling.sif -vv

# Use cached intermediate data
singularity run ... sc2-proportion-modeling.sif --use-cache

# Override output directory
singularity run ... sc2-proportion-modeling.sif --output /app/results

# Override config
singularity run --env CONFIG_PATH=/app/config/alternate.yml ... sc2-proportion-modeling.sif
```

### Environment Variables

```bash
# Set via command line
singularity run \
  --env PYTHONUNBUFFERED=1 \
  --env LOG_LEVEL=DEBUG \
  ... sc2-proportion-modeling.sif

# Or in job script
export PYTHONUNBUFFERED=1
export LOG_LEVEL=DEBUG
singularity run ...
```

---

## Troubleshooting

### Container Won't Start

**Error:** `FATAL: kernel too old`
```bash
# Solution: Use older Singularity container build or upgrade system
# Check kernel version
uname -r
# Need kernel 3.8+, usually available on modern HPCs
```

**Error:** `permission denied` on bind mount
```bash
# Solution: Ensure path exists and is readable
[ -d ./config ] || mkdir -p config
chmod 755 config
```

### Pipeline Fails Inside Container

**Error:** `ModuleNotFoundError: No module named 'polars'`
```bash
# This shouldn't happen - container includes all dependencies
# Check container integrity:
singularity exec sc2-proportion-modeling.sif python -c "import polars; print(polars.__version__)"
```

**Error:** `Impala connection timeout`
```bash
# Check Impala server availability
singularity exec sc2-proportion-modeling.sif \
  python -c "from sc2.config import ImpalaConfig; \
            config = ImpalaConfig(host='your-host'); \
            print(f'Impala endpoint: {config.host}:{config.port}')"

# Test connectivity from host
nc -zv impala-hostname 21050
```

**Error:** `R error: package 'survey' not found`
```bash
# R packages are inside container - shouldn't happen
# Verify R installation
singularity exec sc2-proportion-modeling.sif R --version

# Check R packages
singularity exec sc2-proportion-modeling.sif \
  R -e "library('survey'); print('OK')"
```

### Performance Issues

**Slow data access:**
```bash
# Solution: Use local scratch for temporary files
# See "Advanced Job Script with Data Staging" above
```

**Out of memory (OOM):**
```bash
# Increase SGE resource limit
#$ -l h_vmem=48G  # Increase from 32G

# Or reduce workload:
# - Limit date range (shorter lookback period)
# - Subset variants to track
# - Reduce prediction horizon
```

**Slow R startup:**
```bash
# R startup takes 20-30 seconds (normal for Stan/ggplot2)
# Not a failure - just wait longer
# Total pipeline typically runs in 30-60 seconds
```

---

## Monitoring and Debugging

### Check Job Status

```bash
# List all your jobs
qstat -u $USER

# Detailed job info
qstat -j <job_id>

# Watch job progress
watch -n 5 'qstat | grep sc2'
```

### Review Output Logs

```bash
# SGE stdout/stderr (first 100 lines)
head -100 modeling_20260318_143000.out
head -100 modeling_20260318_143000.err

# Python execution logs
# (Created inside container if LOG_DIR is bound)
ls -la logs/
tail -50 logs/sc2-run.log
```

### Interactive Debugging

```bash
# Start interactive session
singularity shell \
  --bind $(pwd)/config:/app/config \
  --bind $(pwd)/results:/app/results \
  sc2-proportion-modeling.sif

# Inside container, test components
python
>>> from sc2.config import SC2Config
>>> config = SC2Config.from_yaml('/app/config/config.yml')
>>> print(config)
```

### Manual Step-by-Step Execution

```bash
# Run individual pipeline stages
singularity exec ./sc2-proportion-modeling.sif \
  python -c "
from sc2.pipeline.fetch import VariantDataFetcher
from sc2.config import SC2Config
config = SC2Config.from_yaml('/app/config/config.yml')
fetcher = VariantDataFetcher(config.database)
print('Fetcher initialized successfully')
"
```

---

## Performance Characteristics

### Typical Runtime

| Stage | Time | Notes |
|-------|------|-------|
| Fetch | 10-30 sec | Data download from Impala |
| Aggregate | 5-10 sec | Lineage mapping |
| Weight | 15-30 sec | Survey weighting + clustering |
| Model | 20-40 sec | Stan MCMC sampling |
| Export | 30-60 sec | **R startup dominates (~20 sec)** |
| **Total** | **80-170 sec** | **Typically 2-3 minutes** |

### Resource Requirements

**Standard Job:**
```
- CPUs: 1 core (no parallelization needed)
- Memory: 32GB (for Stan posterior samples)
- Runtime: 8 hours (generous upper bound, ~2-3 min typical)
- Disk: 10GB (results + cache)
```

**High-Volume Job (many variants/regions):**
```
- CPUs: 2-4 cores (for parallel Stan chains)
- Memory: 48GB
- Runtime: 12 hours
- Disk: 20GB
```

---

## Production Checklist

- [ ] Singularity/Apptainer module available on HPC
- [ ] Container image downloaded and verified
- [ ] Configuration file created and validated
- [ ] Working directory structure set up
- [ ] Data connectivity tested (Impala connection works)
- [ ] SGE job script created and executable
- [ ] Dry-run successful (`--dry-run` flag passes)
- [ ] Full pipeline executed successfully
- [ ] Results validated (outputs present, reasonable values)
- [ ] Output directory permissions correct
- [ ] Notifications/monitoring set up (optional)
- [ ] Backup strategy for results defined

---

## Advanced: Building Custom Container

If you need to customize the container (e.g., additional R packages):

```bash
# Get the Dockerfile from GitHub
curl -O https://raw.githubusercontent.com/anthropics/sc2-proportion-modeling/master/Dockerfile

# Modify as needed (e.g., add R packages in RUN instruction)

# Build locally (requires Docker)
docker build -t sc2-proportion-modeling:custom .

# Convert to Singularity
singularity build --fakeroot sc2-proportion-modeling-custom.sif docker-daemon://sc2-proportion-modeling:custom

# Use in job script
singularity run ./sc2-proportion-modeling-custom.sif ...
```

---

## Support and Issues

For problems or questions:

1. **Check logs first** - Most issues visible in stdout/stderr
2. **Try dry-run** - Validate config without full execution
3. **Report to:** sc2-team@cdc.gov
4. **Include:**
   - Job ID (qstat output)
   - Error messages and logs
   - Configuration (sanitized)
   - SC2 version (from container)

---

## See Also

- **CI_CD_GUIDE.md** - GitHub Actions and development workflows
- **PHASE6_IMPLEMENTATION.md** - Technical container architecture
- **DEPLOYMENT.md** - General deployment overview
- **README.md** - Main project documentation

Generated: 2026-03-18
