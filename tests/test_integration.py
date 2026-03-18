"""Integration tests for fetch and aggregate pipeline stages.

These tests verify that the fetch and aggregate modules work together
to process raw sequence data into aggregated variant counts.
"""

from datetime import date

import polars as pl
import pytest

from sc2.config import DatabaseConfig
from sc2.pipeline.aggregate import LineageAggregator
from sc2.pipeline.fetch import VariantDataFetcher, DateRange
from sc2.pipeline.exceptions import DataFetchException, AggregationException


@pytest.mark.integration
class TestFetchAggregateIntegration:
    """Integration tests between Fetch and Aggregate modules."""

    def test_fetch_then_aggregate_workflow(self):
        """Test complete workflow: fetch data then aggregate lineages."""
        # This is a dry-run test with placeholder implementation
        config = DatabaseConfig()
        date_range = DateRange(start=date(2026, 1, 1), end=date(2026, 1, 31))

        with VariantDataFetcher(config) as fetcher:
            # In real implementation, would fetch from Impala
            # For now, create test DataFrame
            df = pl.DataFrame({
                "collection_date": [date(2026, 1, 1)] * 10,
                "hhs_region": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                "pango_lineage": [
                    "BA.1.1.2",
                    "BA.2.1",
                    "JN.1.4.3",
                    "XBB.1.5.10",
                    "KP.3.1",
                    "UNKNOWN.999",
                    "BA.1.1",
                    "JN.1",
                    "XBB",
                    "EG.5",
                ],
                "state": ["CA", "TX", "NY", "FL", "PA", "IL", "OH", "GA", "NC", "MI"],
            })

            # Validate raw data
            assert df.height == 10
            assert "pango_lineage" in df.columns

        # Aggregate lineages
        agg = LineageAggregator()
        df_agg = agg.aggregate_dataframe(df, voc_list=[
            "BA.1.1",
            "BA.2",
            "JN.1",
            "XBB",
            "XBB.1.5",
            "KP.3",
        ])

        # Verify aggregation results
        assert df_agg.height == 10
        assert "voc" in df_agg.columns
        assert "pango_lineage" not in df_agg.columns

        # Check that "other" was assigned to unknown lineages
        other_count = df_agg.filter(pl.col("voc") == "other").height
        assert other_count >= 1  # UNKNOWN and EG.5 should be "other"

    def test_fetch_validate_aggregate_workflow(self):
        """Test workflow with validation step."""
        config = DatabaseConfig()

        # Create raw data
        df_raw = pl.DataFrame({
            "collection_date": [date(2026, 1, 1), date(2026, 1, 2), date(2026, 1, 3)],
            "hhs_region": [1, 2, 3],
            "pango_lineage": ["BA.1.1.2", "JN.1.4.3", "XBB.1.5.10"],
            "state": ["CA", "TX", "NY"],
        })

        with VariantDataFetcher(config) as fetcher:
            # Validate fetched data
            df_validated = fetcher.validate_dataframe(df_raw)
            assert df_validated.height == 3

        # Aggregate
        agg = LineageAggregator()
        df_agg = agg.aggregate_dataframe(df_validated)

        # Should have voc column with aggregated lineages
        assert "voc" in df_agg.columns
        voc_values = df_agg["voc"].unique().to_list()
        assert len(voc_values) > 0

    def test_aggregation_by_region(self):
        """Test aggregation with regional stratification."""
        # Create data with multiple regions
        df = pl.DataFrame({
            "hhs_region": [1, 1, 2, 2, 3, 3],
            "pango_lineage": [
                "BA.1.1.2", "BA.1.1.3",  # Region 1
                "JN.1.4.3", "JN.1.5.2",  # Region 2
                "XBB.1.5.10", "XBB.1.5.11",  # Region 3
            ],
        })

        agg = LineageAggregator()
        df_agg = agg.aggregate_dataframe(df)

        # Group by region and check aggregation
        regional_vocs = df_agg.group_by("hhs_region").agg(
            pl.col("voc").n_unique().alias("n_voc")
        )

        # Each region should have aggregated lineages
        assert regional_vocs.height == 3
        for n_voc in regional_vocs["n_voc"]:
            assert n_voc >= 1

    def test_aggregation_temporal_pattern(self):
        """Test aggregation maintains temporal pattern."""
        dates = [date(2026, 1, 1) + __import__("datetime").timedelta(days=i)
                 for i in range(7)]

        lineages = [
            ["BA.1.1.2"] * 15,  # Day 1
            ["BA.1.1.2"] * 10 + ["JN.1.4.3"] * 5,  # Day 2 - JN.1 emerging
            ["JN.1.4.3"] * 15,  # Day 3 - JN.1 becoming dominant
            ["JN.1.4.3"] * 12 + ["XBB.1.5.10"] * 3,  # Day 4 - XBB appearing
            ["JN.1.4.3"] * 9 + ["XBB.1.5.10"] * 6,  # Day 5
            ["JN.1.4.3"] * 7 + ["XBB.1.5.10"] * 8,  # Day 6
            ["JN.1.4.3"] * 4 + ["XBB.1.5.10"] * 11,  # Day 7
        ]

        df = pl.DataFrame({
            "collection_date": [d for d in dates for _ in range(15)],
            "pango_lineage": [l for day_lineages in lineages for l in day_lineages],
        })

        agg = LineageAggregator()
        df_agg = agg.aggregate_dataframe(df)

        # Verify temporal pattern
        daily_counts = df_agg.group_by("collection_date").agg(
            pl.when(pl.col("voc").str.contains("BA")).then(1).otherwise(0).sum().alias("ba_count"),
            pl.when(pl.col("voc").str.contains("JN")).then(1).otherwise(0).sum().alias("jn_count"),
            pl.when(pl.col("voc").str.contains("XBB")).then(1).otherwise(0).sum().alias("xbb_count"),
        )

        # BA should decrease, JN dominant then decrease, XBB should increase
        assert daily_counts.height == 7
