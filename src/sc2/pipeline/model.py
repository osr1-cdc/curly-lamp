"""Statistical modeling stage for SC2 proportion modeling.

Implements nowcasting via multinomial logistic regression with time-series
structure to produce smoothed estimates and predictions of variant proportions.
"""

from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Optional

import polars as pl

from sc2.config import ModelConfig
from sc2.pipeline.exceptions import ModelException


@dataclass
class NowcastPrediction:
    """Nowcast model predictions for variant proportions."""

    prediction_date: date
    variant: str
    proportion_point_estimate: float
    proportion_ci_lower: float
    proportion_ci_upper: float
    retrospective_smooth: bool  # True if fitted; False if predicted


class NowcastModel:
    """Fits multinomial regression nowcast model to variant proportions.

    Uses time-series structure to:
    1. Smooth historical proportion estimates (handles measurement noise)
    2. Predict future proportions (extrapolate trends)
    3. Quantify uncertainty via credible/prediction intervals

    Currently supports Stan Bayesian approach; can be extended to frequentist MLE.

    Handles:
    - Model compilation and caching
    - Data preparation and validation
    - Posterior sampling and diagnostics
    - Prediction interval generation
    """

    def __init__(self, config: ModelConfig, stan_model_path: Optional[Path] = None):
        """Initialize nowcast model.

        Args:
            config: ModelConfig with nowcast method and parameters
            stan_model_path: Path to Stan model file. If None, uses default.

        Raises:
            ModelException: If model path invalid or Stan compilation fails
        """
        self.config = config
        self.stan_model_path = stan_model_path or Path(__file__).parent.parent / "models" / "nowcast.stan"
        self.model = None
        # TODO: Load and compile Stan model

    def fit_and_predict(
        self,
        proportions_df: pl.DataFrame,
        prediction_weeks: int = 4,
    ) -> pl.DataFrame:
        """Fit model to proportions and generate predictions.

        Args:
            proportions_df: Historical variant proportions with columns:
                - date, variant, proportion, se (standard error)
            prediction_weeks: Number of weeks ahead to predict

        Returns:
            DataFrame with retrospective smooths + predictions for each variant/date

        Raises:
            ModelException: If model fitting or prediction fails
        """
        try:
            # Prepare data for Stan
            stan_data = self._prepare_stan_data(proportions_df, prediction_weeks)

            # Fit model
            fit = self._fit_stan_model(stan_data)

            # Extract predictions
            predictions = self._extract_predictions(fit, proportions_df, prediction_weeks)

            return predictions
        except Exception as e:
            raise ModelException(f"Model fitting/prediction failed: {e}") from e

    def _prepare_stan_data(
        self,
        proportions_df: pl.DataFrame,
        prediction_weeks: int,
    ) -> dict:
        """Prepare data in format expected by Stan model.

        Args:
            proportions_df: Historical proportions
            prediction_weeks: Future prediction window

        Returns:
            Dictionary with keys: N_weeks, N_variants, proportions, se, etc.
        """
        # TODO: Transform polars DataFrame to Stan data format
        # Ensure:
        # - Proportions sum to 1 for each week
        # - Standard errors are realistic
        # - Missing data handled (imputation or exclusion)
        # - Prediction period dates defined
        pass

    def _fit_stan_model(self, stan_data: dict) -> "CmdStanMCMC":
        """Fit Stan model via MCMC sampling.

        Args:
            stan_data: Data dictionary for Stan model

        Returns:
            Fitted model results object
        """
        # TODO: Use cmdstanpy to sample from posterior
        # - 4 chains, 2000 iterations
        # - Adapt delta = 0.95
        # - Check diagnostics (Rhat, n_eff)
        pass

    def _extract_predictions(
        self,
        fit,
        proportions_df: pl.DataFrame,
        prediction_weeks: int,
    ) -> pl.DataFrame:
        """Extract posterior predictions from fitted model.

        Args:
            fit: Fitted Stan model result
            proportions_df: Original proportions (for column mapping)
            prediction_weeks: Prediction horizon

        Returns:
            DataFrame with smoothed historical + predicted future proportions
        """
        # TODO: Extract posterior draws and compute:
        # - Credible intervals from posterior quantiles
        # - Point estimates (median or mean)
        # - Prediction intervals for future weeks
        pass
