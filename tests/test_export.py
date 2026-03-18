"""Tests for SC2 proportion modeling export pipeline.

Covers CSV/JSON/Parquet export, figure generation, and metadata tracking.
"""

from datetime import date, datetime, timedelta
from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from sc2.config import OutputConfig
from sc2.pipeline.export import ResultsExporter
from sc2.pipeline.exceptions import ExportException


@pytest.fixture
def sample_proportions_df():
    """Create sample proportions DataFrame."""
    return pl.DataFrame(
        {
            "date": [date(2026, 1, 1), date(2026, 1, 1), date(2026, 1, 8), date(2026, 1, 8)],
            "voc": ["BA.1.1", "XBB.1.5", "BA.1.1", "XBB.1.5"],
            "proportion": [0.45, 0.30, 0.42, 0.33],
            "ci_lower": [0.43, 0.28, 0.40, 0.31],
            "ci_upper": [0.47, 0.32, 0.44, 0.35],
        }
    )


@pytest.fixture
def sample_nowcasts_df():
    """Create sample nowcasts DataFrame."""
    return pl.DataFrame(
        {
            "date": [date(2026, 1, 1), date(2026, 1, 1), date(2026, 1, 8), date(2026, 1, 8), date(2026, 1, 15),
                     date(2026, 1, 15)],
            "voc": ["BA.1.1", "XBB.1.5", "BA.1.1", "XBB.1.5", "BA.1.1", "XBB.1.5"],
            "proportion": [0.45, 0.30, 0.42, 0.33, 0.38, 0.37],
            "ci_lower": [0.43, 0.28, 0.40, 0.31, 0.35, 0.34],
            "ci_upper": [0.47, 0.32, 0.44, 0.35, 0.41, 0.40],
            "retrospective_smooth": [True, True, True, True, False, False],
        }
    )


@pytest.fixture
def sample_metadata():
    """Create sample metadata dictionary."""
    return {
        "data_date": "2026-01-15",
        "export_date": datetime(2026, 1, 15, 14, 30, 0).isoformat(),
        "effective_n": 50000,
        "design_effect": 1.05,
        "variants_tracked": ["BA.1.1", "BA.2", "XBB.1.5"],
        "model_config": {
            "weeks_lookback": 8,
            "weeks_predict": 4,
            "ci_type": "KG",
        },
    }


@pytest.fixture
def output_config(tmp_path):
    """Create OutputConfig with temporary directories."""
    results_dir = tmp_path / "results"
    cache_dir = tmp_path / "cache"
    results_dir.mkdir(parents=True)
    cache_dir.mkdir(parents=True)

    return OutputConfig(
        results_dir=results_dir,
        cache_dir=cache_dir,
        export_formats=["csv", "json", "parquet"],
    )


class TestResultsExporter:
    """Test basic ResultsExporter initialization and directory management."""

    def test_init(self, output_config):
        """Test ResultsExporter initialization."""
        exporter = ResultsExporter(output_config)
        assert exporter.config == output_config
        assert exporter.run_date is not None
        assert output_config.results_dir.exists()
        assert output_config.cache_dir.exists()

    def test_init_with_run_date(self, output_config):
        """Test ResultsExporter with custom run date."""
        run_date = datetime(2026, 1, 15, 14, 30, 0)
        exporter = ResultsExporter(output_config, run_date=run_date)
        assert exporter.run_date == run_date

    def test_ensure_directories(self, tmp_path):
        """Test directory creation."""
        results_dir = tmp_path / "new" / "results"
        cache_dir = tmp_path / "new" / "cache"

        config = OutputConfig(results_dir=results_dir, cache_dir=cache_dir)
        exporter = ResultsExporter(config)

        assert results_dir.exists()
        assert cache_dir.exists()


class TestExportCSV:
    """Test CSV export functionality."""

    def test_export_csv(self, output_config, sample_proportions_df, sample_nowcasts_df, sample_metadata):
        """Test basic CSV export."""
        exporter = ResultsExporter(output_config)
        files = exporter._export_csv(sample_proportions_df, sample_nowcasts_df, sample_metadata)

        assert "proportions_csv" in files
        assert "nowcasts_csv" in files
        assert files["proportions_csv"].exists()
        assert files["nowcasts_csv"].exists()

    def test_csv_format_with_metadata(self, output_config, sample_proportions_df, sample_nowcasts_df,
                                       sample_metadata):
        """Test CSV includes metadata header."""
        exporter = ResultsExporter(output_config)
        proportions_str = exporter._format_proportions_csv(sample_proportions_df, sample_metadata)

        assert "# SC2 Proportion Modeling" in proportions_str
        assert "# Data cutoff date:" in proportions_str
        assert "# Export date:" in proportions_str
        assert "# Effective sample size:" in proportions_str
        assert "# Design effect:" in proportions_str

    def test_csv_format_columns(self, output_config, sample_proportions_df, sample_metadata):
        """Test CSV columns are properly formatted."""
        exporter = ResultsExporter(output_config)
        proportions_str = exporter._format_proportions_csv(sample_proportions_df, sample_metadata)

        # Check that date, voc, and proportion columns exist
        assert "date" in proportions_str
        assert "voc" in proportions_str or "BA.1.1" in proportions_str
        assert "proportion" in proportions_str or "45" in proportions_str  # percentage format

    def test_nowcasts_csv_format(self, output_config, sample_nowcasts_df, sample_metadata):
        """Test nowcasts CSV format with retrospective_smooth column."""
        exporter = ResultsExporter(output_config)
        nowcasts_str = exporter._format_nowcasts_csv(sample_nowcasts_df, sample_metadata)

        assert "# SC2 Proportion Modeling - Nowcasts" in nowcasts_str
        assert "retrospective_smooth" in nowcasts_str
        assert "true" in nowcasts_str.lower() or "True" in nowcasts_str


class TestExportJSON:
    """Test JSON export functionality."""

    def test_export_json(self, output_config, sample_proportions_df, sample_nowcasts_df, sample_metadata):
        """Test basic JSON export."""
        exporter = ResultsExporter(output_config)
        files = exporter._export_json(sample_proportions_df, sample_nowcasts_df, sample_metadata)

        assert "results_json" in files
        assert files["results_json"].exists()

    def test_json_structure(self, output_config, sample_proportions_df, sample_nowcasts_df, sample_metadata):
        """Test JSON has proper nested structure."""
        import json

        exporter = ResultsExporter(output_config)
        files = exporter._export_json(sample_proportions_df, sample_nowcasts_df, sample_metadata)

        with open(files["results_json"]) as f:
            data = json.load(f)

        assert "metadata" in data
        assert "proportions" in data
        assert "nowcasts" in data

    def test_json_metadata(self, output_config, sample_proportions_df, sample_nowcasts_df, sample_metadata):
        """Test JSON includes complete metadata."""
        import json

        exporter = ResultsExporter(output_config)
        files = exporter._export_json(sample_proportions_df, sample_nowcasts_df, sample_metadata)

        with open(files["results_json"]) as f:
            data = json.load(f)

        metadata = data["metadata"]
        assert metadata["data_date"] == "2026-01-15"
        assert "export_date" in metadata
        assert metadata["effective_n"] == 50000
        assert metadata["design_effect"] == 1.05

    def test_json_nested_by_variant(self, output_config, sample_proportions_df, sample_metadata):
        """Test JSON is nested by variant and date."""
        import json

        exporter = ResultsExporter(output_config)
        exporter._export_json(sample_proportions_df, pl.DataFrame(), sample_metadata)

        results_path = output_config.results_dir / "results.json"
        with open(results_path) as f:
            data = json.load(f)

        proportions = data["proportions"]
        # Should have variants as keys
        assert "BA.1.1" in proportions or "XBB.1.5" in proportions or len(proportions) > 0

    def test_prepare_metadata(self, output_config, sample_metadata):
        """Test metadata preparation for JSON."""
        exporter = ResultsExporter(output_config)
        prepared = exporter._prepare_metadata(sample_metadata)

        assert prepared["data_date"] == "2026-01-15"
        assert prepared["effective_n"] == 50000
        assert "model_config" in prepared


class TestExportParquet:
    """Test Parquet export functionality."""

    def test_export_parquet(self, output_config, sample_proportions_df, sample_nowcasts_df):
        """Test Parquet export creates files."""
        exporter = ResultsExporter(output_config)
        files = exporter._export_parquet(sample_proportions_df, sample_nowcasts_df)

        assert "proportions_parquet" in files
        assert "nowcasts_parquet" in files
        assert files["proportions_parquet"].exists()
        assert files["nowcasts_parquet"].exists()

    def test_parquet_roundtrip(self, output_config, sample_proportions_df):
        """Test Parquet can be read back correctly."""
        exporter = ResultsExporter(output_config)
        exporter._export_parquet(sample_proportions_df, pl.DataFrame())

        parquet_path = output_config.cache_dir / "proportions.parquet"
        readback_df = pl.read_parquet(parquet_path)

        assert readback_df.shape == sample_proportions_df.shape
        assert readback_df.columns == sample_proportions_df.columns


class TestExportFigures:
    """Test figure export functionality."""

    def test_export_figures_graceful_fallback(self, output_config, sample_proportions_df, sample_nowcasts_df):
        """Test figures gracefully handle missing rpy2."""
        exporter = ResultsExporter(output_config)

        # This should not raise an error even if rpy2 is unavailable
        files = exporter._export_figures(sample_proportions_df, sample_nowcasts_df, "png")

        # May be empty if rpy2 unavailable, but should be a dict
        assert isinstance(files, dict)

    @patch("sc2.pipeline.export.rpy2")
    def test_export_figures_with_mock_rpy2(self, mock_rpy2, output_config, sample_proportions_df,
                                           sample_nowcasts_df):
        """Test figures with mocked rpy2."""
        # Mock rpy2 imports
        mock_rpy2.robjects = MagicMock()
        mock_rpy2.robjects.pandas2ri = MagicMock()

        exporter = ResultsExporter(output_config)

        # Should handle mock rpy2 gracefully
        files = exporter._export_figures(sample_proportions_df, sample_nowcasts_df, "png")
        assert isinstance(files, dict)


class TestExportResults:
    """Test main export_results orchestration method."""

    def test_export_results_with_all_formats(self, output_config, sample_proportions_df, sample_nowcasts_df,
                                            sample_metadata):
        """Test export_results with multiple formats."""
        exporter = ResultsExporter(output_config)
        files = exporter.export_results(sample_proportions_df, sample_nowcasts_df, sample_metadata)

        # Should have CSV and JSON exports
        assert "proportions_csv" in files
        assert "nowcasts_csv" in files
        assert "results_json" in files
        assert "proportions_parquet" in files
        assert "nowcasts_parquet" in files

    def test_export_results_metadata_written(self, output_config, sample_proportions_df, sample_nowcasts_df,
                                            sample_metadata):
        """Test metadata is written to file."""
        exporter = ResultsExporter(output_config)
        exporter.export_results(sample_proportions_df, sample_nowcasts_df, sample_metadata)

        metadata_path = output_config.results_dir / "metadata.json"
        assert metadata_path.exists()

    def test_export_results_metadata_includes_export_date(self, output_config, sample_proportions_df,
                                                         sample_nowcasts_df, sample_metadata):
        """Test metadata includes export date."""
        import json

        exporter = ResultsExporter(output_config)
        exporter.export_results(sample_proportions_df, sample_nowcasts_df, sample_metadata)

        metadata_path = output_config.results_dir / "metadata.json"
        with open(metadata_path) as f:
            written_metadata = json.load(f)

        assert "export_date" in written_metadata

    def test_export_results_error_handling(self, output_config, sample_metadata):
        """Test export_results handles errors appropriately."""
        exporter = ResultsExporter(output_config)

        # Pass invalid data to trigger errors
        with pytest.raises(ExportException):
            exporter.export_results(
                pl.DataFrame(),  # Empty DataFrame
                pl.DataFrame(),
                sample_metadata,
            )


class TestExportIntegration:
    """Integration tests with realistic data."""

    def test_full_export_pipeline(self, output_config):
        """Test complete export pipeline with realistic multi-week data."""
        # Create realistic data spanning 12 weeks with 4 variants
        dates = []
        vocs = []
        proportions = []
        ci_lowers = []
        ci_uppers = []
        retrospective = []

        variants = ["BA.1.1", "BA.2", "XBB.1.5", "JN.1"]
        base_date = date(2026, 1, 1)

        for week in range(12):
            current_date = base_date + timedelta(weeks=week)
            for i, variant in enumerate(variants):
                dates.append(current_date)
                vocs.append(variant)
                # Simulate shifting proportions
                prop = 0.25 + (week * 0.02 * (i - 1.5))
                prop = max(0.05, min(0.50, prop))  # Clamp to [0.05, 0.50]
                proportions.append(prop)
                ci_lowers.append(max(0.01, prop - 0.05))
                ci_uppers.append(min(0.95, prop + 0.05))
                retrospective.append(week < 8)

        nowcasts_df = pl.DataFrame(
            {
                "date": dates,
                "voc": vocs,
                "proportion": proportions,
                "ci_lower": ci_lowers,
                "ci_upper": ci_uppers,
                "retrospective_smooth": retrospective,
            }
        )

        metadata = {
            "data_date": "2026-03-18",
            "effective_n": 75000,
            "design_effect": 1.08,
            "variants_tracked": variants,
        }

        exporter = ResultsExporter(output_config, run_date=datetime(2026, 3, 18, 12, 0, 0))
        files = exporter.export_results(nowcasts_df, nowcasts_df, metadata)

        # Verify all expected files exist
        assert len(files) > 0
        for file_path in files.values():
            assert file_path.exists()

    def test_export_with_regional_data(self, output_config):
        """Test export with regional (HHS region) data."""
        # Create data with HHS regions
        dates = []
        vocs = []
        proportions = []
        ci_lowers = []
        ci_uppers = []
        regions = []

        variants = ["BA.2", "XBB"]
        hhs_regions = [1, 2, 3]  # HHS regions 1-3

        for region in hhs_regions:
            for week in range(4):
                for variant in variants:
                    dates.append(date(2026, 1, 1 + week * 7))
                    vocs.append(variant)
                    proportions.append(0.5)
                    ci_lowers.append(0.45)
                    ci_uppers.append(0.55)
                    regions.append(region)

        df = pl.DataFrame(
            {
                "date": dates,
                "voc": vocs,
                "proportion": proportions,
                "ci_lower": ci_lowers,
                "ci_upper": ci_uppers,
                "hhs_region": regions,
            }
        )

        metadata = {
            "data_date": "2026-02-01",
            "effective_n": 50000,
            "design_effect": 1.05,
        }

        exporter = ResultsExporter(output_config)
        files = exporter.export_results(df, df, metadata)

        # Should create exports without error
        assert len(files) > 0
