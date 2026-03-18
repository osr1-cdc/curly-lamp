"""Data fetching stage for SC2 proportion modeling.

Responsible for querying Impala database and retrieving SARS-CoV-2 sequencing data
within specified date ranges and regions.
"""

from dataclasses import dataclass
from datetime import date, timedelta
from typing import Optional

import polars as pl

from sc2.config import DatabaseConfig
from sc2.pipeline.exceptions import DataFetchException


@dataclass
class DateRange:
    """Represents a date range for data queries."""

    start: date
    end: date

    def __post_init__(self):
        if self.start > self.end:
            raise ValueError(f"Invalid date range: {self.start} > {self.end}")

    def days(self) -> int:
        """Get number of days in range (inclusive)."""
        return (self.end - self.start).days + 1


class VariantDataFetcher:
    """Fetches SARS-CoV-2 sequencing data from Impala database.

    Handles:
    - Connection management
    - Query construction with date/region filters
    - Data validation
    - Error handling
    """

    def __init__(self, config: DatabaseConfig):
        """Initialize fetcher with database configuration.

        Args:
            config: DatabaseConfig instance with Impala connection parameters

        Raises:
            ConfigurationException: If config is invalid
        """
        self.config = config
        self.conn = None
        # TODO: Initialize Impala connection in __enter__ for context manager support

    def fetch_sequences(
        self,
        date_range: DateRange,
        regions: Optional[list[str]] = None,
    ) -> pl.DataFrame:
        """Fetch sequences from database within date range and regions.

        Args:
            date_range: DateRange object specifying query time window
            regions: Optional list of HHS regions. If None, fetch all regions.

        Returns:
            Polars DataFrame with columns:
            - collection_date: Date sequence was submitted
            - hhs_region: HHS region assignment
            - pango_lineage: Pangolin lineage classification
            - ... (other sequence metadata)

        Raises:
            DataFetchException: If query fails or connection issues occur
        """
        try:
            query = self._build_query(date_range, regions)
            df = pl.read_database(query, connection=self.conn)
            return df
        except Exception as e:
            raise DataFetchException(f"Failed to fetch sequences: {e}") from e

    def _build_query(self, date_range: DateRange, regions: Optional[list[str]] = None) -> str:
        """Build SQL query for sequence data.

        Args:
            date_range: Date range for filtering
            regions: Optional list of regions

        Returns:
            SQL query string (parameterized to prevent injection)
        """
        # TODO: Implement SQL query construction
        # This will query the covid.sequencing_data table with:
        # - Date filters on collection_date
        # - Region filters on hhs_region (if specified)
        # - All columns needed for downstream processing
        pass

    def __enter__(self):
        """Context manager entry."""
        # TODO: Establish database connection
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        # TODO: Close database connection
        pass
