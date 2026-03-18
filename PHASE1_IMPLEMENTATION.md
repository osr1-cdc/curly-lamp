# Phase 1: Python Skeleton + YAML Configuration

## Overview

Phase 1 establishes the foundational Python package structure and configuration system for SC2 Proportion Modeling v2.0.

## Completed Work

### 1. Directory Structure

```
src/sc2/
├── __init__.py                 # Package root
├── config.py                   # Pydantic configuration schema
├── pipeline/
│   ├── __init__.py
│   ├── exceptions.py           # Custom pipeline exceptions
│   ├── fetch.py                # Data fetching module (stub)
│   ├── aggregate.py            # Lineage aggregation module (stub)
│   ├── weight.py               # Survey weighting module (stub)
│   ├── model.py                # Statistical modeling module (stub)
│   └── export.py               # Results export module (stub)
├── models/
│   ├── __init__.py
│   └── nowcast.stan            # Stan nowcast model
└── scripts/
    ├── __init__.py
    └── run.py                  # CLI entry point

data/
├── config.yml                  # Example production configuration

tests/
├── __init__.py
├── conftest.py                 # Pytest fixtures
├── test_config.py              # Configuration tests
└── test_pipeline.py            # Pipeline module tests
```

### 2. Configuration System (src/sc2/config.py)

**Features:**
- Pydantic v2 schema for type-safe validation
- Hierarchical config structure (RunConfig, DatabaseConfig, VariantConfig, ModelConfig, OutputConfig)
- YAML loading/exporting with `from_yaml()` / `to_yaml()` methods
- Environment variable support via pydantic-settings
- Example config dictionary for documentation

**Validation:**
- Data date parsing (string -> date object)
- Port validation for Impala connection
- Credible interval bounds (0-1)
- Variant list non-empty enforcement

**Usage:**
```python
# From YAML
config = SC2Config.from_yaml(Path("data/config.yml"))

# From environment
config = SC2Config.from_env()

# Programmatic
config = SC2Config(
    run=RunConfig(data_date="2026-03-18", results_tag="CDT"),
    variants=VariantConfig(voc1=["BA.1.1", "BA.2"])
)

# Export to YAML
config.to_yaml(Path("config_backup.yml"))
```

### 3. Core Modules (Pipeline Stubs)

All modules follow a consistent pattern:
- **Type annotations** for clarity and IDE support
- **Docstrings** with full documentation
- **TODO comments** marking implementation points
- **Exception handling** with custom exceptions

**Modules:**
- `fetch.py`: VariantDataFetcher class (Impala query architecture)
- `aggregate.py`: LineageAggregator class (hierarchical lineage mapping)
- `weight.py`: SurveyWeighter class (rpy2 integration for R survey package)
- `model.py`: NowcastModel class (Stan model orchestration)
- `export.py`: ResultsExporter class (multi-format output: CSV, JSON, Parquet, PNG/PDF)
- `exceptions.py`: Custom exception hierarchy for error handling

### 4. Stan Model (src/sc2/models/nowcast.stan)

**Features:**
- Multinomial logistic regression with time-series structure
- Random walk priors on log-odds for temporal smoothing
- Per-variant hyperpriors on trend SD
- Posterior smoothing for historical periods
- Prediction intervals for future periods (nowcasting)
- Identifiability via (K-1) log-odds parameterization

**Key Components:**
```stan
parameters {
  matrix theta[N_weeks, N_variants-1]    // Log-odds (latent)
  real tau_variant[N_variants]           // Per-variant trend SD
}

generated quantities {
  // Posterior predictions for extrapolation
  matrix p_pred[N_weeks, N_variants]     // Predicted proportions
}
```

### 5. CLI Entry Point (src/sc2/scripts/run.py)

**Commands:**
- `sc2-run`: Main pipeline execution
  - `--config`: YAML configuration file path
  - `--output`: Override output directory
  - `--use-cache`: Use previously imported data
  - `--dry-run`: Validate configuration only
  - `--verbose` / `-v`: Increase logging verbosity

- `sc2-run validate-config`: Validate configuration file syntax

**Features:**
- Structured logging via loguru
- Graceful error handling with stack traces (verbose mode)
- Configuration loading from YAML or environment variables
- Placeholder for full pipeline execution (Phase 2+)

**Usage Examples:**
```bash
# Standard monthly run
python -m sc2.scripts.run --config data/config.yml

# Dry run validation
python -m sc2.scripts.run --config data/config.yml --dry-run

# Custom output path
python -m sc2.scripts.run --config data/config.yml --output results_2026-03-18_CDT

# Verbose debugging
python -m sc2.scripts.run --config data/config.yml -vvv

# Configuration validation
python -m sc2.scripts.run validate-config data/config.yml
```

### 6. Dependency Management (pyproject.toml)

**Build System:**
- Modern setuptools-based configuration
- Python 3.11+ requirement

**Core Dependencies:**
- `polars>=0.20.0`: Efficient DataFrame processing (replaces pandas for speed)
- `pydantic>=2.0.0`: Configuration schema and validation
- `pyyaml>=6.0`: YAML file handling
- `cmdstanpy>=1.1.0`: Stan model compilation and sampling
- `rpy2>=3.16.0`: R integration for survey package
- `loguru>=0.7.0`: Structured logging
- `typer>=0.9.0`: CLI framework for subcommands

**Development Dependencies:**
- `pytest>=7.4.0`: Testing framework
- `pytest-cov>=4.1.0`: Code coverage tracking
- `black>=23.0.0`: Code formatting
- `isort>=5.12.0`: Import sorting
- `mypy>=1.4.0`: Type checking
- `ruff>=0.0.280`: Fast linting

**Tool Configuration:**
- `black`: Line length 100, Python 3.11 target
- `isort`: Black-compatible import sorting
- `ruff`: Common lint rules (E, F, W, I, UP)
- `pytest`: Test discovery in tests/ directory, -v output, coverage reporting
- `mypy`: Basic type checking with missing import allowance

**Entry Points:**
- `sc2-run`: Maps to `sc2.scripts.run:main`
- `sc2-monthly`: Maps to `sc2.scripts.monthly_report:main` (future)

### 7. Example Configuration (data/config.yml)

**Provides:**
- Impala connection parameters
- Run date and results tagging
- Variant tracking list (voc1: 28 variants)
- Reduced comparison set (voc1_reduced: 2 variants)
- Model parameters (weeks_lookback=8, weeks_predict=4)
- Output configuration (formats: csv, json, parquet, png)

### 8. Testing Framework

**Test Files:**
- `test_config.py`: Configuration validation tests
  - Date parsing (string and date object inputs)
  - YAML loading/exporting
  - Default values
  - Optional field handling

- `test_pipeline.py`: Pipeline module scaffolding
  - Placeholder tests for Phase 2+ implementation
  - Marked with `@pytest.skip` for implementation

- `conftest.py`: Shared pytest fixtures
  - `sample_config`: Pre-configured SC2Config instance
  - `tmp_results_dir`: Temporary results directory

**Running Tests:**
```bash
# Install with dev dependencies
pip install -e ".[dev]"

# Run all tests with coverage
pytest --cov=src/sc2

# Run specific test file
pytest tests/test_config.py -v

# Run only integration tests
pytest -m integration

# Run with verbose output
pytest -vv
```

### 9. Updated .gitignore

Added comprehensive Python entries:
- `__pycache__/`, `.pyc` files
- Virtual environment directories (venv/, env/, .venv)
- Build artifacts (dist/, build/, *.egg-info/)
- pytest cache and coverage (.pytest_cache/, .coverage)
- IDE files (.vscode/, .idea/, *.swp)
- Poetry lock file (poetry.lock)
- Jupyter notebooks (.ipynb, .ipynb_checkpoints/)

### 10. Package Documentation

**Package Docstring** (src/sc2/__init__.py):
```
SC2 Proportion Modeling v2.0 - Modern Python Implementation

A modular, testable, and maintainable Python pipeline for CDC SARS-CoV-2
proportion modeling with statistical nowcasting using survey design-based
weighting and multinomial regression analysis.
```

## What's Implemented

✅ Configuration schema with pydantic validation
✅ YAML config file loading and exporting
✅ CLI entry point with subcommands
✅ Pipeline module stubs with full docstrings and type hints
✅ Stan nowcast model (mathematical formulation)
✅ Dependency management (pyproject.toml)
✅ Basic test framework with fixtures
✅ Custom exception hierarchy
✅ Logging infrastructure (loguru)

## What's Not Yet Implemented

⏳ Impala database connection logic (Phase 2)
⏳ Lineage aggregation implementation (Phase 2)
⏳ Survey weighting calculations (Phase 3)
⏳ Stan model compilation and fitting (Phase 4)
⏳ Results export to CSV/JSON/Parquet (Phase 5)
⏳ Figure generation (Phase 5)
⏳ Docker containerization (Phase 6)
⏳ CI/CD GitHub Actions (Phase 6)

## Next Steps (Phase 2)

1. **Implement Fetch Module**
   - Impala JDBC connection setup
   - Parameterized SQL query builder
   - Data validation and logging

2. **Implement Aggregate Module**
   - Load Pangolin lineage hierarchy from CSV
   - Implement hierarchical prefix matching
   - Test against real lineage examples

3. **Integration Tests**
   - Mock database queries with test data
   - Verify dataframe transformations
   - Validate output shapes and types

## Installation & Development

```bash
# Clone repository and switch to branch
git clone https://path/to/repo
cd sc2_proportion_modeling
git checkout feature/python-rewrite

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install in development mode with all dependencies
pip install -e ".[dev]"

# Run tests
pytest --cov=src/sc2

# Format code
black src/ tests/
isort src/ tests/

# Type check
mypy src/

# Lint
ruff check src/
```

## Design Rationale

This Phase 1 implementation establishes:

1. **Type Safety**: Pydantic schemas catch configuration errors at load time
2. **Modularity**: Each pipeline stage as independent class for testability
3. **Extensibility**: Stan model and rpy2 bridge enable algorithm evolution
4. **Maintainability**: Clear docstrings, comprehensive type hints, fixture-based tests
5. **Reproducibility**: YAML config versioning, dependency locking via pyproject.toml
6. **Professionalism**: Modern Python patterns (dataclasses, context managers, logging)
