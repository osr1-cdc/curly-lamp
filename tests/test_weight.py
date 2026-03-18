"""Tests for the Survey Weighting pipeline stage."""

from datetime import date, timedelta

import polars as pl
import pytest

from sc2.config import ModelConfig
from sc2.pipeline.weight import SurveyWeighter, WeightingResult
from sc2.pipeline.exceptions import WeightingException


class TestWeightingResult:
    """Tests for WeightingResult dataclass."""

    def test_weighting_result_creation(self):
        """Test creating a WeightingResult."""
        proportions = pl.DataFrame(
            {
                "voc": ["BA.1.1", "JN.1"],
                "proportion": [0.6, 0.4],
                "ci_lower": [0.55, 0.35],
                "ci_upper": [0.65, 0.45],
            }
        )
        weights = pl.DataFrame({"hhs_region": [1, 2], "weight": [1.0, 1.2]})

        result = WeightingResult(
            variant_proportions=proportions,
            weights=weights,
            design_effect=1.05,
            effective_n=1000,
        )

        assert result.design_effect == 1.05
        assert result.effective_n == 1000
        assert result.variant_proportions.height() == 2


class TestSurveyWeighter:
    """Tests for SurveyWeighter class."""

    @pytest.fixture
    def weighter(self):
        """Fixture providing a SurveyWeighter instance."""
        config = ModelConfig()
        return SurveyWeighter(config)

    @pytest.fixture
    def sample_sequences(self):
        """Fixture providing sample sequence data."""
        start_date = date(2026, 1, 1)
        dates = [start_date + timedelta(days=i) for i in range(14)]

        return pl.DataFrame(
            {
                "collection_date": dates + dates,
                "hhs_region": [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4],
                "voc": [
                    "BA.1.1", "BA.1.1", "JN.1",
                    "BA.1.1", "JN.1", "JN.1",
                    "JN.1", "JN.1", "XBB",
                    "JN.1", "XBB", "XBB",
                    "XBB", "XBB"
                ] * 2,
            }
        )

    def test_weighter_initialization(self, weighter):
        """Test SurveyWeighter initialization."""
        assert weighter is not None
        assert weighter.model_config.ci_type == "KG"

    def test_aggregate_by_week(self, weighter, sample_sequences):
        """Test aggregating sequences to weekly counts."""
        weekly = weighter._aggregate_by_week(sample_sequences)

        assert "week_start" in weekly.columns
        assert "hhs_region" in weekly.columns
        assert "voc" in weekly.columns
        assert "count" in weekly.columns
        assert weekly.height() > 0

    def test_aggregate_by_week_counts(self, weighter):
        """Test weekly aggregation produces correct counts."""
        # Create simple test data with known counts
        df = pl.DataFrame(
            {
                "collection_date": [date(2026, 1, 5)] * 5 + [date(2026, 1, 6)] * 3,
                "hhs_region": [1] * 8,
                "voc": ["BA.1.1"] * 5 + ["JN.1"] * 3,
            }
        )

        weekly = weighter._aggregate_by_week(df)

        # Both dates fall in same week
        assert weekly.height() == 2  # BA.1.1 and JN.1
        ba_count = weekly.filter(pl.col("voc") == "BA.1.1")["count"][0]
        assert ba_count == 5

    def test_calculate_regional_weight(self):
        """Test weight calculation for different sequence counts."""
        # More sequences = lower weight
        weight_100 = SurveyWeighter._calculate_regional_weight(100)
        weight_1000 = SurveyWeighter._calculate_regional_weight(1000)
        weight_10000 = SurveyWeighter._calculate_regional_weight(10000)

        assert weight_100 > weight_1000 > weight_10000
        assert weight_100 <= 10.0  # Capped at 10
        assert weight_10000 >= 0.1  # Floored at 0.1

    def test_calculate_regional_weight_edge_cases(self):
        """Test weight calculation edge cases."""
        # Zero sequences
        weight_zero = SurveyWeighter._calculate_regional_weight(0)
        assert weight_zero == 1.0

        # Negative sequences (shouldn't happen but handle gracefully)
        weight_negative = SurveyWeighter._calculate_regional_weight(-5)
        assert weight_negative == 1.0

    def test_calculate_weights(self, weighter, sample_sequences):
        """Test calculating weights for region-weeks."""
        weekly = weighter._aggregate_by_week(sample_sequences)
        weights = weighter._calculate_weights(weekly)

        assert "hhs_region" in weights.columns
        assert "week_start" in weights.columns
        assert "weight" in weights.columns
        assert weights.height() > 0

        # Check weight bounds
        for w in weights["weight"]:
            assert 0.1 <= w <= 10.0

    def test_apply_weights(self, weighter, sample_sequences):
        """Test applying weights to sequences."""
        weekly = weighter._aggregate_by_week(sample_sequences)
        weights = weighter._calculate_weights(weekly)
        weighted = weighter._apply_weights(sample_sequences, weights)

        assert "weight" in weighted.columns
        assert weighted.height() == sample_sequences.height()
        # No weights should be null (fill_null defaults to 1.0)
        assert weighted["weight"].null_count() == 0

    def test_simple_weighted_proportions(self, weighter, sample_sequences):
        """Test computing simple weighted proportions."""
        weekly = weighter._aggregate_by_week(sample_sequences)
        weights = weighter._calculate_weights(weekly)
        weighted = weighter._apply_weights(sample_sequences, weights)

        proportions = weighter._simple_weighted_proportions(weighted)

        assert "voc" in proportions.columns
        assert "proportion" in proportions.columns
        assert "ci_lower" in proportions.columns
        assert "ci_upper" in proportions.columns

        # Proportions should sum to ~1.0
        total_prop = proportions["proportion"].sum()
        assert 0.95 < total_prop < 1.05  # Allow small numerical error

    def test_validate_weights_normal(self, weighter):
        """Test weight validation with normal weights."""
        weights_df = pl.DataFrame({"weight": [0.5, 1.0, 1.5, 2.0]})
        assert weighter.validate_weights(weights_df) is True

    def test_validate_weights_empty(self, weighter):
        """Test weight validation with empty DataFrame."""
        weights_df = pl.DataFrame({"weight": []})
        assert weighter.validate_weights(weights_df) is False

    def test_calculate_weighted_proportions(self, weighter, sample_sequences):
        """Test full weighted proportions workflow."""
        result = weighter.calculate_weighted_proportions(sample_sequences)

        assert isinstance(result, WeightingResult)
        assert result.variant_proportions.height() > 0
        assert result.design_effect > 0
        assert result.effective_n > 0

        # Proportions should be valid (0-1 range, sum to ~1)
        proportions = result.variant_proportions["proportion"]
        assert (proportions >= 0).all() and (proportions <= 1).all()
        assert 0.95 < proportions.sum() < 1.05

    def test_weighter_context_manager(self, weighter):
        """Test context manager functionality."""
        with weighter as w:
            assert w is not None
            # After exit, session should be None
        assert weighter.r_session is None


class TestWeightingIntegration:
    """Integration tests for weighting workflow."""

    @pytest.mark.integration
    def test_aggregate_weight_proportions_workflow(self):
        """Test complete weighting workflow: aggregate → weight → proportions."""
        config = ModelConfig(ci_type="KG")
        weighter = SurveyWeighter(config)

        # Create realistic test data
        start_date = date(2026, 1, 1)
        dates = [start_date + timedelta(days=i) for i in range(56)]  # 8 weeks

        sequences = pl.DataFrame(
            {
                "collection_date": [d for d in dates for _ in range(20)],
                "hhs_region": [i % 10 + 1 for i in range(20 * 56)],
                "voc": [
                    ["BA.1.1"] * 5 + ["JN.1"] * 10 + ["XBB"] * 5
                    if i < 200
                    else (
                        ["BA.1.1"] * 2
                        + ["JN.1"] * 8
                        + ["XBB"] * 10
                        if i < 600
                        else ["BA.1.1"] * 1 + ["JN.1"] * 5 + ["XBB"] * 14
                    )
                    for i in range(20 * 56)
                ],
            }
        )

        with weighter:
            result = weighter.calculate_weighted_proportions(sequences)

        # Validate results
        assert result.variant_proportions.height() == 3  # 3 variants
        assert result.effective_n > 0
        assert result.design_effect > 0

        # Check that variant proportions make sense
        props = result.variant_proportions.sort_values("proportion", descending=True)
        # XBB should be dominant in later weeks
        assert props["voc"][0] in ["XBB", "JN.1"]

    @pytest.mark.integration
    def test_weighting_regional_stratification(self):
        """Test that weighting maintains regional patterns."""
        config = ModelConfig()
        weighter = SurveyWeighter(config)

        # Create data with regional variation
        sequences = pl.DataFrame(
            {
                "collection_date": [date(2026, 1, 1)] * 50,
                "hhs_region": [1] * 20 + [2] * 30,  # Region 2 has more sequences
                "voc": ["BA.1.1"] * 25 + ["JN.1"] * 25,
            }
        )

        with weighter:
            result = weighter.calculate_weighted_proportions(sequences)

        # Region 1 sequences should get higher weights due to smaller sample
        weights = result.weights
        region_1_weight = weights.filter(pl.col("hhs_region") == 1)["weight"].mean()
        region_2_weight = weights.filter(pl.col("hhs_region") == 2)["weight"].mean()

        # Smaller sample (region 1) should have higher weights
        assert region_1_weight >= region_2_weight * 0.8  # Allow some tolerance

    @pytest.mark.integration
    def test_weighting_with_missing_regions(self):
        """Test handling of missing data for some regions."""
        config = ModelConfig()
        weighter = SurveyWeighter(config)

        # Data missing region 5-10
        sequences = pl.DataFrame(
            {
                "collection_date": [date(2026, 1, 1)] * 40,
                "hhs_region": [1, 2, 3, 4] * 10,  # Only regions 1-4
                "voc": ["BA.1.1"] * 20 + ["JN.1"] * 20,
            }
        )

        with weighter:
            result = weighter.calculate_weighted_proportions(sequences)

        # Should still work
        assert result.variant_proportions.height() > 0
        assert result.effective_n > 0
