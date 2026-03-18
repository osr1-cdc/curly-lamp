"""Tests for pipeline modules."""

import pytest

from sc2.pipeline.aggregate import LineageAggregator
from sc2.pipeline.exceptions import AggregationException


class TestLineageAggregator:
    """Tests for LineageAggregator."""

    def test_aggregator_initialization(self):
        """Test aggregator can be initialized."""
        agg = LineageAggregator()
        assert agg is not None

    def test_aggregate_lineage_single(self):
        """Test aggregating a single lineage."""
        agg = LineageAggregator()
        # TODO: Once aggregation rules are implemented, test with real lineages
        # result = agg.aggregate_lineage("JN.1.4.3")
        # assert result == "JN.1"
        pass

    @pytest.mark.skip(reason="Implementation pending")
    def test_aggregate_dataframe(self):
        """Test aggregating lineages in DataFrame."""
        pass

    @pytest.mark.skip(reason="Implementation pending")
    def test_aggregate_with_voc_list(self):
        """Test aggregation respects VOC list."""
        pass
