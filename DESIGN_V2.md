# SC2 Proportion Modeling – Rewrite Design Document

## 1. Language Choice: Python with R for Statistics

### Why Python Primary + R Secondary

**Data Pipeline Orchestration**: Python excels at ETL workflows, multiprocessing, and cloud integration for future-proofing
**Database Querying & Data Manipulation**: Pandas/Polars are far more ergonomic than R for Impala queries and large datasets
**Configuration & CLI**: Python's click/typer provide better CLI tooling than R's optparse/argparse
**Packaging & Deployment**: Poetry/uv solve dependency issues that plague R renv
**Statistical Components**: Call R via rpy2 or embed Stan models directly for nowcasting regression
**Reproducibility**: Docker + Python environment superior to shell scripts + qsub workflows

### Technical Stack

```
Primary Layer (Python 3.12+):
  • Data: polars (faster than pandas for this use case)
  • Config: pydantic + YAML
  • Database: sqlalchemy + Impala connector
  • Orchestration: Prefect or Dask for HPC parity
  • Testing: pytest + hypothesis

Secondary Layer (R 4.4+):
  • Survey analysis: survey package (already in use)
  • Nowcasting model: cmdstanr or brms for Bayesian multinomial
  • Visualization: ggplot2/plotly for reproducible outputs
  • Bridging: rpy2 or subprocess for calling R functions
```

## 2. Architecture: Modular Data Pipeline

```
sc2-proportion-modeling/
├── src/
│   ├── pipeline/
│   │   ├── fetch.py         # Query Impala, apply date filters
│   │   ├── aggregate.py     # Lineage aggregation logic
│   │   ├── weight.py        # Survey weighting calculations
│   │   ├── model.py         # Nowcast model execution
│   │   └── export.py        # Output generation (CSV, JSON, figures)
│   ├── models/
│   │   └── nowcast.stan     # Stan model for multinomial regression
│   ├── config.py            # Config schema + loading
│   └── __init__.py
├── data/
│   ├── cache/               # Intermediate results (polars parquet)
│   ├── config.yml           # Production configuration
│   └── lineage_hierarchy.csv # Pangolin aggregation rules
├── scripts/
│   ├── run.py              # Main entry point (replaces proportion_modeling_run.sh)
│   └── monthly_report.py   # CLI for standard monthly workflow
├── analysis/
│   └── nowcast_analysis.Rmd # Exploratory R notebook for model development
├── tests/
│   ├── test_aggregation.py
│   ├── test_weighting.py
│   └── test_model.py
├── docker/
│   └── Dockerfile          # Reproducible environment
├── pyproject.toml          # Dependencies + metadata
└── README.md
```

## 3. Core Implementation Strategy

### Fetch Module (replaces variant_surveillance_system.R)

```python
class VariantDataFetcher:
    def __init__(self, impala_config):
        self.conn = ImpalaConnection(impala_config)

    def fetch_sequences(self, date_range, regions):
        # Single query, no custom lineage conditionals
        query = f"""
        SELECT collection_date, hhs_region, pango_lineage, ...
        FROM covid.sequencing_data
        WHERE collection_date BETWEEN {date_range[0]} AND {date_range[1]}
        """
        return pl.from_arrow(self.conn.execute(query))
```

### Aggregate Module

```python
class LineageAggregator:
    def __init__(self, hierarchy_df):
        self.hierarchy = hierarchy_df

    def aggregate_to_voc(self, lineages: list[str]) -> list[str]:
        # Deterministic aggregation based on Pangolin hierarchy
        # No custom lineage branching logic
        return self._apply_aggregation_rules(lineages)
```

### Weight Module (calls R for survey package)

```python
class SurveyWeighter:
    def __init__(self, r_session):
        self.r = r_session

    def calculate_weighted_proportions(self, data_df, design_spec):
        # Call R's survey package, return clean DataFrame
        proportions = self.r.calculate_survey_estimates(data_df, design_spec)
        return pl.from_dict(proportions)
```

### Model Module (Stan for nowcasting)

```python
class NowcastModel:
    def __init__(self):
        self.model = stan.load_model(Path("models/nowcast.stan"))

    def fit_and_predict(self, proportions_df, prediction_weeks=4):
        # Multinomial logistic regression via Stan
        # Returns: point estimates + credible intervals
        fit = self.model.sample(data=prepare_stan_data(proportions_df))
        return fit.draw_predictions(weeks_ahead=prediction_weeks)
```

## 4. Configuration as Code

```yaml
# config.yml
run:
  data_date: 2026-03-18
  results_tag: CDT

database:
  impala:
    host: impala-host
    port: 21050

variants:
  voc1:
    - BA.1.1
    - BA.2
    # ... (explicit list, no custom lineage branching)

model:
  weeks_lookback: 8
  weeks_predict: 4
  nowcast_method: multinomial_regression
  credible_interval: 0.95
```

## 5. Key Improvements Over Current Code

| Current | Improved |
|---------|----------|
| Shell orchestration + R scripts | Python pipeline + R components |
| Conditional branching for Run 3, custom lineages, Nextclade | Single canonical workflow |
| Hard-coded paths, environment variables | Config file + dependency injection |
| qsub jobs scattered across scripts | Unified job orchestration (Dask/Prefect) |
| Manual test runs, debugging via comments | Pytest + integration tests |
| Docker: None | Reproducible containerized environment |
| CI/CD: None | GitHub Actions for monthly runs & PR validation |

## 6. Workflow Execution

```bash
# Single intuitive entry point
python scripts/run.py \
  --config data/config.yml \
  --output results_2026-03-18_CDT/ \
  --run-nowcast

# Or for monthly automation:
python scripts/monthly_report.py $(date +%Y-%m-%d)
```

## 7. Why This Approach Works

✅ **Maintainability**: Clear separation of concerns (fetch → aggregate → weight → model → export)
✅ **Testing**: Each module testable in isolation with pytest
✅ **Reproducibility**: Docker + Poetry lock file + version control
✅ **Scalability**: Dask can handle 10x larger datasets without rework
✅ **Modern practices**: Type hints, pydantic validation, logging
✅ **Documentation**: Self-documenting code structure
✅ **Collaboration**: Python is more widely known than R shell scripting

## 8. Phased Migration Path

### Phase 1: Python Skeleton + YAML Config (1-2 weeks)
- Set up directory structure and pyproject.toml
- Create config schema with pydantic
- Implement config loading from YAML
- Create base classes and type definitions
- Set up pytest structure

### Phase 2: Fetch & Aggregate Modules (1-2 weeks)
- Implement VariantDataFetcher for Impala queries
- Implement LineageAggregator with Pangolin hierarchy
- Create integration tests against test dataset

### Phase 3: Weight Module (1 week)
- Implement rpy2 bridge to R's survey package
- Call existing R weighting functions from Python
- Validate outputs match current code

### Phase 4: Modernize Nowcast Model (2-3 weeks)
- Design Stan model for multinomial with time-series
- Implement StanModelWrapper class
- Port existing R nowcast logic to Python workflow

### Phase 5: Export Module + Testing (1 week)
- Implement output CSV/JSON generation
- Implement figure generation (ggplot2 via rpy2)
- Full test coverage

### Phase 6: Docker + CI/CD Setup (1 week)
- Create Dockerfile for reproducible environment
- Set up GitHub Actions for monthly runs
- Add pre-commit hooks and linting

### Phase 7: Validation (2-3 weeks)
- Run side-by-side with current code
- Verify outputs match to numerical precision
- Document any methodology fixes

## 9. Backward Compatibility & Validation

During Phase 7, both pipelines will run in parallel:
- **Old code**: shell scripts + qsub on current branch
- **New code**: Python pipeline with all outputs logged
- **Comparison**: CSV diffs on all outputs (proportions, confidence intervals, nowcast predictions)
- **Tolerance**: Allow 1e-10 relative error for floating-point rounding

Once validated, old scripts can be retired and existing monthly reports switch to Python pipeline.

## 10. Future Extensions Enabled by This Architecture

- **Cloud deployment**: Easily migrate to AWS Lambda or Google Cloud Run
- **Interactive dashboards**: Connect Plotly Dash for real-time monitoring
- **Model improvements**: Swap Stan for PyMC3 or probabilistic programming framework without pipeline changes
- **Parallelization**: Process regions/timepoints independently with Dask
- **Continuous retraining**: Automate updates as new sequences arrive
- **Feature engineering**: Add covariates (vaccination rates, testing volume, etc.)
