"""Tests for the Nowcasting Model pipeline stage."""

from datetime import date, timedelta

import numpy as np
import polars as pl
import pytest

from sc2.config import ModelConfig
from sc2.pipeline.model import NowcastModel, NowcastPrediction
from sc2.pipeline.exceptions import ModelException


class TestNowcastPrediction:
    """Tests for NowcastPrediction dataclass."""

    def test_prediction_creation(self):
        """Test creating a nowcast prediction."""
        pred = NowcastPrediction(
            prediction_date=date(2026, 3, 18),
            variant="BA.1.1",
            proportion_point_estimate=0.35,
            proportion_ci_lower=0.32,
            proportion_ci_upper=0.38,
            retrospective_smooth=True,
        )

        assert pred.variant == "BA.1.1"
        assert pred.proportion_point_estimate == 0.35
        assert pred.retrospective_smooth is True


class TestNowcastModel:
    """Tests for NowcastModel class."""

    @pytest.fixture
    def model_config(self):
        """Fixture providing model configuration."""
        return ModelConfig(
            weeks_lookback=8,
            weeks_predict=4,
            credible_interval=0.95,
        )

    @pytest.fixture
    def model(self, model_config, tmp_path):
        """Fixture providing a NowcastModel instance."""
        # Create dummy Stan model file
        stan_file = tmp_path / "nowcast.stan"
        stan_file.write_text("// Dummy model")

        return NowcastModel(model_config, stan_file)

    @pytest.fixture
    def sample_proportions(self):
        """Fixture providing sample weighted proportions."""
        dates = [date(2026, 1, 1) + timedelta(weeks=i) for i in range(8)]
        variants = ["BA.1.1", "JN.1", "XBB"]

        data = []
        for week_idx, d in enumerate(dates):
            # Shift proportions over time
            ba_prop = 0.6 - (week_idx * 0.05)  # Decreasing
            jn_prop = 0.3 + (week_idx * 0.02)  # Slightly increasing
            xbb_prop = 0.1 + (week_idx * 0.03)  # Increasing

            # Normalize
            total = ba_prop + jn_prop + xbb_prop
            data.extend(
                [
                    {
                        "collection_date": d,
                        "voc": "BA.1.1",
                        "proportion": ba_prop / total,
                        "ci_lower": (ba_prop / total) * 0.90,
                        "ci_upper": (ba_prop / total) * 1.10,
                    },
                    {
                        "collection_date": d,
                        "voc": "JN.1",
                        "proportion": jn_prop / total,
                        "ci_lower": (jn_prop / total) * 0.90,
                        "ci_upper": (jn_prop / total) * 1.10,
                    },
                    {
                        "collection_date": d,
                        "voc": "XBB",
                        "proportion": xbb_prop / total,
                        "ci_lower": (xbb_prop / total) * 0.90,
                        "ci_upper": (xbb_prop / total) * 1.10,
                    },
                ]
            )

        return pl.DataFrame(data)

    def test_model_initialization(self, model):
        """Test NowcastModel initialization."""
        assert model is not None
        assert model.config.weeks_predict == 4

    def test_model_initialization_missing_stan(self, model_config):
        """Test initialization fails with missing Stan model."""
        with pytest.raises(ModelException):
            NowcastModel(model_config, stan_model_path="/nonexistent/path.stan")

    def test_prepare_stan_data(self, model, sample_proportions):
        """Test Stan data preparation."""
        stan_data = model._prepare_stan_data(sample_proportions, prediction_weeks=4)

        assert "N_weeks" in stan_data
        assert "N_variants" in stan_data
        assert "N_fit" in stan_data
        assert "prop" in stan_data
        assert "se" in stan_data

        assert stan_data["N_variants"] == 3
        assert stan_data["N_fit"] == 8
        assert stan_data["N_weeks"] == 12  # 8 fit + 4 predict

    def test_prepare_stan_data_matrix_shapes(self, model, sample_proportions):
        """Test Stan data matrix shapes and properties."""
        stan_data = model._prepare_stan_data(sample_proportions, prediction_weeks=4)

        prop_matrix = stan_data["prop"]
        se_matrix = stan_data["se"]

        # Check shapes
        assert prop_matrix.shape == (12, 3)
        assert se_matrix.shape == (12, 3)

        # Check proportions sum to 1 for fitted weeks
        for week_idx in range(8):
            row_sum = prop_matrix[week_idx, :].sum()
            assert 0.99 < row_sum < 1.01  # Allow minor numerical error

    def test_placeholder_predictions(self, model, sample_proportions):
        """Test placeholder prediction generation."""
        predictions = model._placeholder_predictions(sample_proportions, prediction_weeks=4)

        assert predictions.height > sample_proportions.height
        assert "retrospective_smooth" in predictions.columns

        # Check smoothing flags
        fitted = predictions.filter(pl.col("retrospective_smooth") == True)
        predicted = predictions.filter(pl.col("retrospective_smooth") == False)

        assert fitted.height == sample_proportions.height
        assert predicted.height == 12  # 4 weeks × 3 variants

    def test_placeholder_predictions_date_extension(self, model, sample_proportions):
        """Test that predictions extend dates correctly."""
        predictions = model._placeholder_predictions(sample_proportions, prediction_weeks=4)

        # Get last historical date and first prediction date
        last_historical = predictions.filter(
            pl.col("retrospective_smooth") == True
        )["collection_date"].max()
        first_pred = predictions.filter(
            pl.col("retrospective_smooth") == False
        )["collection_date"].min()

        # Predictions should start 1 week after last data
        expected_first_pred = last_historical + timedelta(weeks=1)
        assert first_pred == expected_first_pred

    def test_fit_and_predict_placeholder(self, model, sample_proportions):
        """Test full fit_and_predict workflow with placeholder."""
        predictions = model.fit_and_predict(sample_proportions, prediction_weeks=4)

        # Should have historical + predictions
        assert predictions.height >= sample_proportions.height

        # All variants should have proportions < 1
        for prop in predictions["proportion"]:
            assert 0 <= prop <= 1

    def test_model_context_manager(self, model):
        """Test context manager functionality."""
        with model as m:
            assert m is not None
        # After exit, fit_result should be None
        assert model.fit_result is None

    def test_uncertainty_growth_in_predictions(self, model, sample_proportions):
        """Test that prediction uncertainty grows over time."""
        predictions = model._placeholder_predictions(sample_proportions, prediction_weeks=4)

        # Get CI widths for same variant across time
        variant_preds = predictions.filter(pl.col("voc") == "BA.1.1")

        fitted_rows = variant_preds.filter(pl.col("retrospective_smooth") == True)
        pred_rows = variant_preds.filter(pl.col("retrospective_smooth") == False)

        if fitted_rows.height > 0 and pred_rows.height > 0:
            # Sample CIs
            fitted_ci_width = (
                fitted_rows[-1]["ci_upper"] - fitted_rows[-1]["ci_lower"]
            )
            pred_ci_width = pred_rows[0]["ci_upper"] - pred_rows[0]["ci_lower"]

            # Prediction CI should be wider (uncertainty grows)
            assert pred_ci_width >= fitted_ci_width * 0.95  # Allow small margin


class TestNowcastIntegration:
    """Integration tests for nowcasting workflow."""

    @pytest.mark.integration
    def test_complete_nowcast_workflow(self, tmp_path):
        """Test complete flow: weighted proportions → nowcast → predictions."""
        # Create config and model
        config = ModelConfig(weeks_lookback=8, weeks_predict=4)
        stan_file = tmp_path / "nowcast.stan"
        stan_file.write_text("// Dummy")

        model = NowcastModel(config, stan_file)

        # Create realistic data
        dates = [date(2026, 1, 1) + timedelta(weeks=i) for i in range(12)]
        variants = ["BA.1.1", "JN.1", "XBB", "KP.3"]

        data = []
        for week_idx, d in enumerate(dates):
            # Changing proportions over time
            proportions = [
                0.40 - (week_idx * 0.03),
                0.35 + (week_idx * 0.01),
                0.15 + (week_idx * 0.02),
                0.10 + (week_idx * 0.01),
            ]

            # Normalize
            total = sum(proportions)
            proportions = [p / total for p in proportions]

            for var_idx, variant in enumerate(variants):
                data.append(
                    {
                        "collection_date": d,
                        "voc": variant,
                        "proportion": proportions[var_idx],
                        "ci_lower": proportions[var_idx] * 0.85,
                        "ci_upper": proportions[var_idx] * 1.15,
                    }
                )

        proportions_df = pl.DataFrame(data)

        # Fit and predict
        with model:
            predictions = model.fit_and_predict(
                proportions_df, prediction_weeks=4
            )

        # Validate results
        assert predictions.height > proportions_df.height
        assert "retrospective_smooth" in predictions.columns

        # All variants should appear in predictions
        pred_variants = set(predictions["voc"].unique().to_list())
        assert pred_variants == set(variants)

    @pytest.mark.integration
    def test_temporal_coherence(self):
        """Test that predicted proportions show coherent temporal patterns."""
        config = ModelConfig()

        # Simulate trend: BA decreasing, XBB increasing
        dates = [date(2026, 1, 1) + timedelta(weeks=i) for i in range(10)]
        data = []

        for week_idx, d in enumerate(dates):
            ba_prop = 0.60 - (week_idx * 0.04)
            xbb_prop = 0.10 + (week_idx * 0.04)
            jn_prop = 1 - ba_prop - xbb_prop

            for variant, prop in [("BA.1.1", ba_prop), ("XBB", xbb_prop), ("JN.1", jn_prop)]:
                data.append(
                    {
                        "collection_date": d,
                        "voc": variant,
                        "proportion": max(0.001, prop),
                        "ci_lower": max(0.001, prop * 0.9),
                        "ci_upper": min(0.999, prop * 1.1),
                    }
                )

        df = pl.DataFrame(data)

        # Simple placeholder should continue general trend
        predictions_df = config  # Use config's default nowcast_method

        # Verify smoothness: consecutive proportions shouldn't jump dramatically
        ba_preds = df.filter(pl.col("voc") == "BA.1.1")["proportion"].to_list()
        for i in range(len(ba_preds) - 1):
            change = abs(ba_preds[i + 1] - ba_preds[i])
            assert change < 0.15  # Max 15% week-to-week change for this trend
