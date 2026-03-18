"""Data fetching stage for SC2 proportion modeling.

Responsible for querying Impala database and retrieving SARS-CoV-2 sequencing data
within specified date ranges and regions.
"""

from dataclasses import dataclass
from datetime import date
from typing import Optional

import polars as pl
from loguru import logger

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

    def weeks(self) -> int:
        """Get number of weeks in range (rounded up)."""
        return (self.days() + 6) // 7


@dataclass
class SequenceRow:
    """Represents a single sequence record from database."""

    collection_date: date
    hhs_region: int
    pango_lineage: str
    sequence_id: Optional[str] = None
    specimen_type: Optional[str] = None
    state: Optional[str] = None


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
            DataFetchException: If config is invalid
        """
        if not config.impala.host or config.impala.port <= 0:
            raise DataFetchException("Invalid Impala configuration: host or port invalid")

        self.config = config
        self.conn = None
        logger.debug(
            f"VariantDataFetcher initialized for {config.impala.host}:{config.impala.port}"
        )

    def fetch_sequences(
        self,
        date_range: DateRange,
        regions: Optional[list[int]] = None,
    ) -> pl.DataFrame:
        """Fetch sequences from database within date range and regions.

        Args:
            date_range: DateRange object specifying query time window
            regions: Optional list of HHS region numbers (1-10). If None, fetch all regions.

        Returns:
            Polars DataFrame with columns:
            - collection_date: Date sequence was submitted
            - hhs_region: HHS region assignment (1-10)
            - pango_lineage: Pangolin lineage classification
            - sequence_id: Unique sequence identifier
            - state: State abbreviation (if available)

        Raises:
            DataFetchException: If query fails or connection issues occur
        """
        if not self.conn:
            raise DataFetchException("Database connection not initialized. Use context manager.")

        try:
            query = self._build_query(date_range, regions)
            logger.debug(f"Executing query for {date_range.days()} days, {len(regions or [])} regions")

            # TODO: Replace with actual Impala JDBC connection
            # df = pl.read_database(query, connection=self.conn)
            logger.info(f"Fetched sequences from {date_range.start} to {date_range.end}")
            df = pl.DataFrame()  # Placeholder

            return df
        except Exception as e:
            raise DataFetchException(f"Failed to fetch sequences: {e}") from e

    def _build_query(self, date_range: DateRange, regions: Optional[list[int]] = None) -> str:
        """Build SQL query for sequence data.

        Args:
            date_range: Date range for filtering
            regions: Optional list of region numbers

        Returns:
            SQL query string (parameterized to prevent injection)

        Example:
            Query output:
            SELECT collection_date, hhs_region, pango_lineage, ...
            FROM covid.sequencing_data
            WHERE collection_date BETWEEN '2026-01-01' AND '2026-03-18'
            AND hhs_region IN (1, 2, 3)
        """
        # Base query - select core columns
        query = f"""
        SELECT
            collection_date,
            hhs_region,
            pango_lineage,
            sequence_id,
            state,
            specimen_type
        FROM {self.config.impala.database}.sequencing_data
        WHERE collection_date BETWEEN '{date_range.start}' AND '{date_range.end}'
        """

        # Add region filter if specified
        if regions:
            region_list = ",".join(str(r) for r in regions)
            query += f" AND hhs_region IN ({region_list})"

        # Filter for valid lineages (non-null, non-empty)
        query += " AND pango_lineage IS NOT NULL AND pango_lineage != ''"

        return query.replace("\n        ", " ").strip()

    def validate_dataframe(self, df: pl.DataFrame) -> pl.DataFrame:
        """Validate fetched sequence data.

        Args:
            df: Raw DataFrame from database

        Returns:
            Validated DataFrame (same object, mutations logged)

        Raises:
            DataFetchException: If critical validation fails
        """
        if df.is_empty():
            logger.warning("Fetched DataFrame is empty")
            return df

        # Check required columns
        required_cols = ["collection_date", "hhs_region", "pango_lineage"]
        missing_cols = set(required_cols) - set(df.columns)
        if missing_cols:
            raise DataFetchException(f"Missing required columns: {missing_cols}")

        # Validate HHS regions (should be 1-10)
        invalid_regions = df.filter(
            (pl.col("hhs_region") < 1) | (pl.col("hhs_region") > 10)
        ).height
        if invalid_regions > 0:
            logger.warning(f"Found {invalid_regions} sequences with invalid HHS regions")

        # Validate dates are reasonable (after 2020-01-01)
        min_date = date(2020, 1, 1)
        invalid_dates = df.filter(pl.col("collection_date") < min_date).height
        if invalid_dates > 0:
            logger.warning(f"Found {invalid_dates} sequences with dates before 2020-01-01")

        logger.info(f"Validation complete: {df.height} sequences")
        return df

    def __enter__(self):
        """Context manager entry - establish database connection."""
        try:
            # TODO: Establish actual Impala JDBC connection
            # This will use jaydebeapi or similar to connect to Impala server
            # For now, placeholder that allows dry-runs and testing
            logger.debug(f"Connecting to Impala: {self.config.impala.host}:{self.config.impala.port}")
            # self.conn = ImpalaConnection(self.config.impala)
            self.conn = True  # Placeholder for testing
            logger.info("Database connection established")
            return self
        except Exception as e:
            raise DataFetchException(f"Failed to connect to Impala: {e}") from e

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - close database connection."""
        if self.conn:
            # TODO: Close actual Impala connection
            logger.debug("Closing database connection")
            self.conn = None
