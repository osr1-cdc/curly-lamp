"""Output export stage for SC2 proportion modeling.

Responsible for generating and exporting results in multiple formats:
- CSV tables for public reporting
- JSON for downstream systems
- Parquet for efficient storage/caching
- PNG/PDF figures for dashboards and reports
"""

from datetime import datetime
from pathlib import Path
from typing import Optional

import polars as pl

from sc2.config import OutputConfig
from sc2.pipeline.exceptions import ExportException


class ResultsExporter:
    """Exports SC2 proportion modeling results in multiple formats.

    Handles:
    - Standardized output directory structure
    - Format-specific serialization (CSV, JSON, Parquet)
    - Figure generation and styling
    - Metadata tracking (run date, data cutoff, parameters)
    """

    def __init__(self, config: OutputConfig, run_date: Optional[datetime] = None):
        """Initialize exporter.

        Args:
            config: OutputConfig specifying output paths and formats
            run_date: Timestamp of analysis run (defaults to now)
        """
        self.config = config
        self.run_date = run_date or datetime.now()
        self._ensure_directories()

    def export_results(
        self,
        proportions: pl.DataFrame,
        nowcasts: pl.DataFrame,
        metadata: dict,
    ) -> dict[str, Path]:
        """Export all results in configured formats.

        Args:
            proportions: Variant proportions with confidence intervals
            nowcasts: Smoothed/predicted proportions from nowcast model
            metadata: Analysis metadata (data_date, results_tag, etc.)

        Returns:
            Dictionary mapping format -> output file paths

        Raises:
            ExportException: If export fails
        """
        try:
            output_files = {}

            for fmt in self.config.export_formats:
                if fmt == "csv":
                    output_files.update(self._export_csv(proportions, nowcasts, metadata))
                elif fmt == "json":
                    output_files.update(self._export_json(proportions, nowcasts, metadata))
                elif fmt == "parquet":
                    output_files.update(self._export_parquet(proportions, nowcasts))
                elif fmt in ("png", "pdf"):
                    output_files.update(self._export_figures(proportions, nowcasts, fmt))

            self._write_metadata(metadata)
            return output_files
        except Exception as e:
            raise ExportException(f"Failed to export results: {e}") from e

    def _export_csv(
        self,
        proportions: pl.DataFrame,
        nowcasts: pl.DataFrame,
        metadata: dict,
    ) -> dict[str, Path]:
        """Export results as CSV files.

        Args:
            proportions: Proportions DataFrame
            nowcasts: Nowcasts DataFrame
            metadata: Metadata dictionary

        Returns:
            Mapping of description -> output file path
        """
        # TODO: Export with human-readable formatting
        # - proportions.csv: Main results table
        # - nowcasts.csv: Smoothed/predicted proportions
        # - metadata.json: Run parameters and data cutoff date
        files = {}
        prop_path = self.config.results_dir / "proportions.csv"
        proportions.write_csv(prop_path)
        files["proportions_csv"] = prop_path

        nowcast_path = self.config.results_dir / "nowcasts.csv"
        nowcasts.write_csv(nowcast_path)
        files["nowcasts_csv"] = nowcast_path

        return files

    def _export_json(
        self,
        proportions: pl.DataFrame,
        nowcasts: pl.DataFrame,
        metadata: dict,
    ) -> dict[str, Path]:
        """Export results as JSON files.

        Args:
            proportions: Proportions DataFrame
            nowcasts: Nowcasts DataFrame
            metadata: Metadata dictionary

        Returns:
            Mapping of description -> output file path
        """
        # TODO: Export with nested structure for easy API consumption
        # - Include metadata in root
        # - Organize by variant and date
        files = {}
        output_path = self.config.results_dir / "results.json"
        # Combine data with metadata
        combined = {
            "metadata": metadata,
            "proportions": proportions.to_dicts(),
            "nowcasts": nowcasts.to_dicts(),
        }
        # Write JSON
        import json

        with open(output_path, "w") as f:
            json.dump(combined, f, indent=2, default=str)
        files["results_json"] = output_path

        return files

    def _export_parquet(
        self,
        proportions: pl.DataFrame,
        nowcasts: pl.DataFrame,
    ) -> dict[str, Path]:
        """Export results as Parquet files (efficient storage).

        Args:
            proportions: Proportions DataFrame
            nowcasts: Nowcasts DataFrame

        Returns:
            Mapping of description -> output file path
        """
        # TODO: Export for fast reload and caching
        files = {}
        prop_path = self.config.cache_dir / "proportions.parquet"
        proportions.write_parquet(prop_path)
        files["proportions_parquet"] = prop_path

        nowcast_path = self.config.cache_dir / "nowcasts.parquet"
        nowcasts.write_parquet(nowcast_path)
        files["nowcasts_parquet"] = nowcast_path

        return files

    def _export_figures(
        self,
        proportions: pl.DataFrame,
        nowcasts: pl.DataFrame,
        fmt: str,
    ) -> dict[str, Path]:
        """Export visualization figures.

        Args:
            proportions: Proportions DataFrame
            nowcasts: Nowcasts DataFrame
            fmt: Output format ("png" or "pdf")

        Returns:
            Mapping of description -> output file path
        """
        # TODO: Generate figures using R's ggplot2 via rpy2
        # - Variant proportion time series with CI bands
        # - Nowcast smoothing + prediction visualization
        # - Multiple regions on subplots
        files = {}
        # Placeholder implementation
        return files

    def _write_metadata(self, metadata: dict) -> None:
        """Write analysis metadata to JSON file.

        Args:
            metadata: Metadata dictionary
        """
        import json

        metadata_path = self.config.results_dir / "metadata.json"
        metadata["export_date"] = self.run_date.isoformat()
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2, default=str)

    def _ensure_directories(self) -> None:
        """Create output directories if they don't exist."""
        self.config.results_dir.mkdir(parents=True, exist_ok=True)
        self.config.cache_dir.mkdir(parents=True, exist_ok=True)
