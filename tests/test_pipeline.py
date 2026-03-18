"""Tests for pipeline modules."""

from datetime import date

import polars as pl
import pytest

from sc2.config import DatabaseConfig, ImpalaConfig
from sc2.pipeline.aggregate import LineageAggregator, AggregationRule
from sc2.pipeline.fetch import VariantDataFetcher, DateRange
from sc2.pipeline.exceptions import AggregationException, DataFetchException


class TestDateRange:
    """Tests for DateRange helper class."""

    def test_date_range_valid(self):
        """Test creating valid date range."""
        dr = DateRange(start=date(2026, 1, 1), end=date(2026, 1, 31))
        assert dr.days() == 31

    def test_date_range_single_day(self):
        """Test date range for single day."""
        dr = DateRange(start=date(2026, 1, 1), end=date(2026, 1, 1))
        assert dr.days() == 1

    def test_date_range_invalid_order(self):
        """Test date range with reversed dates."""
        with pytest.raises(ValueError):
            DateRange(start=date(2026, 1, 31), end=date(2026, 1, 1))

    def test_date_range_weeks(self):
        """Test weeks calculation."""
        dr = DateRange(start=date(2026, 1, 1), end=date(2026, 1, 8))  # 8 days
        assert dr.weeks() == 2  # Rounded up: (8+6)//7 = 2


class TestVariantDataFetcher:
    """Tests for VariantDataFetcher."""

    def test_fetcher_initialization(self):
        """Test fetcher can be initialized."""
        config = DatabaseConfig()
        fetcher = VariantDataFetcher(config)
        assert fetcher is not None

    def test_fetcher_invalid_config(self):
        """Test fetcher initialization with invalid config."""
        config = DatabaseConfig(
            impala=ImpalaConfig(host="", port=0)
        )
        with pytest.raises(DataFetchException):
            VariantDataFetcher(config)

    def test_fetcher_context_manager(self):
        """Test fetcher can be used as context manager."""
        config = DatabaseConfig()
        with VariantDataFetcher(config) as fetcher:
            assert fetcher.conn is not None

    def test_build_query_basic(self):
        """Test SQL query building with basic date range."""
        config = DatabaseConfig()
        fetcher = VariantDataFetcher(config)

        date_range = DateRange(start=date(2026, 1, 1), end=date(2026, 1, 31))
        query = fetcher._build_query(date_range)

        assert "SELECT" in query
        assert "collection_date BETWEEN" in query
        assert "2026-01-01" in query
        assert "2026-01-31" in query
        assert "pango_lineage IS NOT NULL" in query

    def test_build_query_with_regions(self):
        """Test SQL query building with region filter."""
        config = DatabaseConfig()
        fetcher = VariantDataFetcher(config)

        date_range = DateRange(start=date(2026, 1, 1), end=date(2026, 1, 31))
        query = fetcher._build_query(date_range, regions=[1, 2, 3])

        assert "hhs_region IN (1,2,3)" in query

    def test_validate_dataframe_empty(self):
        """Test validation of empty DataFrame."""
        config = DatabaseConfig()
        fetcher = VariantDataFetcher(config)

        df = pl.DataFrame()
        result = fetcher.validate_dataframe(df)
        assert result.is_empty()

    def test_validate_dataframe_missing_columns(self):
        """Test validation fails with missing required columns."""
        config = DatabaseConfig()
        fetcher = VariantDataFetcher(config)

        df = pl.DataFrame({"collection_date": [date(2026, 1, 1)]})

        with pytest.raises(DataFetchException):
            fetcher.validate_dataframe(df)


class TestLineageAggregator:
    """Tests for LineageAggregator."""

    def test_aggregator_initialization(self):
        """Test aggregator initializes with default rules."""
        agg = LineageAggregator()
        assert agg is not None
        assert len(agg.rules) > 0

    def test_aggregate_lineage_exact_match(self):
        """Test aggregating lineage with exact match."""
        agg = LineageAggregator()
        # BA.1.1 is in the default rules
        result = agg.aggregate_lineage("BA.1.1")
        assert result == "BA.1.1"

    def test_aggregate_lineage_prefix_match(self):
        """Test aggregating lineage with prefix matching."""
        agg = LineageAggregator()
        # JN.1.4.3 should match JN.1 prefix
        result = agg.aggregate_lineage("JN.1.4.3")
        # Should match one of the JN.1 sublineages or JN.1 prefix
        assert "JN.1" in result

    def test_aggregate_lineage_xbb_family(self):
        """Test aggregating XBB family lineages."""
        agg = LineageAggregator()
        # XBB.1.5.10 should map to XBB.1.5 hierarchy
        result = agg.aggregate_lineage("XBB.1.5.10.3")
        # Should be in XBB family
        assert "XBB" in result

    def test_aggregate_lineage_unknown(self):
        """Test aggregating unknown lineage returns original."""
        agg = LineageAggregator()
        # Unknown lineage should be returned as-is
        result = agg.aggregate_lineage("UNKNOWN.999")
        assert result == "UNKNOWN.999"

    def test_aggregate_dataframe(self):
        """Test aggregating lineages in a DataFrame."""
        agg = LineageAggregator()

        df = pl.DataFrame({
            "pango_lineage": ["BA.1.1", "BA.1.1.2", "JN.1.4.3", "UNKNOWN.999"],
            "collection_date": [
                date(2026, 1, 1),
                date(2026, 1, 1),
                date(2026, 1, 1),
                date(2026, 1, 1),
            ],
        })

        result = agg.aggregate_dataframe(df)

        assert "voc" in result.columns
        assert "pango_lineage" not in result.columns
        assert result.height() == 4

    def test_aggregate_dataframe_with_voc_list(self):
        """Test aggregation enforces VOC list."""
        agg = LineageAggregator()

        df = pl.DataFrame({
            "pango_lineage": ["BA.1.1", "JN.1.4.3", "UNKNOWN.999"],
        })

        # Only allow BA.1.1 and variants
        result = agg.aggregate_dataframe(df, voc_list=["BA.1.1"])

        # Other variants should become "other"
        assert "other" in result["voc"].to_list()

    def test_aggregate_with_custom_rules(self):
        """Test aggregator with custom aggregation rules."""
        custom_rules = [
            AggregationRule(parent="TEST.1", children=["TEST.1.1", "TEST.1.2"]),
        ]
        agg = LineageAggregator(rules=custom_rules)

        # Exact match should work
        assert agg.aggregate_lineage("TEST.1.1") == "TEST.1"

        # Prefix should work
        assert agg.aggregate_lineage("TEST.1.1.5") == "TEST.1"

    def test_build_lookup_tables(self):
        """Test lookup table construction."""
        agg = LineageAggregator()

        # parent_map should contain exact lineages
        assert len(agg.parent_map) > 0

        # prefix_map should contain parent lineages
        assert len(agg.prefix_map) > 0

        # Verify some known mappings
        assert "BA.1.1" in agg.prefix_map
        assert "XBB.1.5" in agg.prefix_map

