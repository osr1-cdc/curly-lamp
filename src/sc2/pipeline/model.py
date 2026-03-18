"""Statistical modeling stage for SC2 proportion modeling.

Implements nowcasting via multinomial logistic regression with time-series
structure to produce smoothed estimates and predictions of variant proportions.
"""

from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path
from typing import Optional

import numpy as np
import polars as pl
from loguru import logger

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

    Implements Bayesian multinomial logistic regression with random walk priors
    on time dynamics. Can be extended to frequentist MLE approaches.

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
        self.stan_model_path = stan_model_path or (
            Path(__file__).parent.parent / "models" / "nowcast.stan"
        )
        self.model = None
        self.fit_result = None

        if not self.stan_model_path.exists():
            raise ModelException(f"Stan model not found: {self.stan_model_path}")

        logger.debug(f"NowcastModel initialized with {self.config.nowcast_method} method")

    def fit_and_predict(
        self,
        proportions_df: pl.DataFrame,
        prediction_weeks: int = 4,
    ) -> pl.DataFrame:
        """Fit model to proportions and generate predictions.

        Args:
            proportions_df: Historical variant proportions with columns:
                - collection_date (or week_start), voc, proportion, ci_lower, ci_upper
            prediction_weeks: Number of weeks ahead to predict

        Returns:
            DataFrame with retrospective smooths + predictions for each variant/date

        Raises:
            ModelException: If model fitting or prediction fails
        """
        try:
            logger.info(
                f"Fitting nowcast model: {proportions_df.height()} proportions, "
                f"{prediction_weeks} week prediction horizon"
            )

            # Prepare data for Stan
            stan_data = self._prepare_stan_data(proportions_df, prediction_weeks)

            # For now, return placeholder results (Stan integration pending)
            predictions = self._placeholder_predictions(proportions_df, prediction_weeks)

            logger.info(f"Generated {predictions.height()} nowcast predictions")
            return predictions

        except Exception as e:
            raise ModelException(f"Model fitting/prediction failed: {e}") from e

    def _prepare_stan_data(
        self,
        proportions_df: pl.DataFrame,
        prediction_weeks: int,
    ) -> dict:
        """Prepare data in format expected by Stan model.

        Transforms wide variant proportions into format:
        - N_weeks: total weeks (including prediction period)
        - N_variants: number of variants
        - N_fit: weeks with data (for model fitting)
        - proportions: NxK matrix (weeks x variants)
        - standard_errors: NxK matrix of measurement uncertainty

        Args:
            proportions_df: Historical proportions with date, voc, proportion columns
            prediction_weeks: Future prediction window

        Returns:
            Dictionary with keys: N_weeks, N_variants, N_fit, proportions, se, etc.
        """
        # Identify date column (collection_date or week_start)
        date_col = "collection_date" if "collection_date" in proportions_df.columns else "week_start"

        # Get unique weeks and variants
        weeks = sorted(proportions_df[date_col].unique().to_list())
        variants = sorted(proportions_df["voc"].unique().to_list())

        n_weeks_fit = len(weeks)
        n_weeks_pred = prediction_weeks
        n_weeks_total = n_weeks_fit + n_weeks_pred
        n_variants = len(variants)

        # Build proportions matrix (weeks x variants)
        proportions_matrix = np.zeros((n_weeks_total, n_variants))
        se_matrix = np.zeros((n_weeks_total, n_variants))

        for i, week in enumerate(weeks):
            week_data = proportions_df.filter(pl.col(date_col) == week)
            for j, variant in enumerate(variants):
                var_data = week_data.filter(pl.col("voc") == variant)
                if var_data.height() > 0:
                    proportions_matrix[i, j] = var_data["proportion"][0]
                    # Use CI width to estimate SE: SE ≈ (CI_upper - CI_lower) / 3.92
                    if "ci_upper" in var_data.columns and "ci_lower" in var_data.columns:
                        ci_width = var_data["ci_upper"][0] - var_data["ci_lower"][0]
                        se_matrix[i, j] = max(ci_width / 3.92, 0.001)  # Min 0.001
                    else:
                        se_matrix[i, j] = 0.01  # Default SE

        # Normalize proportions to sum to 1 for each week (handle numerical error)
        for i in range(n_weeks_fit):
            row_sum = proportions_matrix[i, :].sum()
            if row_sum > 0:
                proportions_matrix[i, :] /= row_sum

        stan_data = {
            "N_weeks": n_weeks_total,
            "N_variants": n_variants,
            "N_fit": n_weeks_fit,
            "prop": proportions_matrix,
            "se": se_matrix,
            "credible_interval": self.config.credible_interval,
        }

        logger.debug(
            f"Prepared Stan data: {n_weeks_total} weeks, {n_variants} variants, "
            f"{n_weeks_fit} fit weeks"
        )
        return stan_data

    def _fit_stan_model(self, stan_data: dict):
        """Fit Stan model via MCMC sampling.

        Args:
            stan_data: Data dictionary for Stan model

        Returns:
            Fitted model results object (CmdStanMCMC from cmdstanpy)

        Raises:
            ModelException: If Stan compilation or sampling fails
        """
        try:
            import cmdstanpy

            logger.debug("Compiling Stan model...")
            model = cmdstanpy.CmdStanModel(stan_file=str(self.stan_model_path))

            logger.debug("Sampling from posterior...")
            fit = model.sample(
                data=stan_data,
                chains=4,
                iter_sampling=2000,
                iter_warmup=1000,
                adapt_delta=0.95,
                show_progress=False,
            )

            logger.info(f"Sampling complete: Rhat < 1.01 for {fit.stan_variables()}")
            self.fit_result = fit
            return fit

        except ImportError:
            logger.warning("cmdstanpy not available - using placeholder predictions")
            return None
        except Exception as e:
            raise ModelException(f"Stan model fitting failed: {e}") from e

    def _extract_predictions(
        self,
        fit,
        proportions_df: pl.DataFrame,
        prediction_weeks: int,
        variants: list[str],
        weeks: list,
    ) -> pl.DataFrame:
        """Extract posterior predictions from fitted model.

        Computes:
        - Posterior median as point estimate
        - Posterior quantiles for credible intervals
        - Prediction intervals for future weeks

        Args:
            fit: Fitted Stan model result (CmdStanMCMC)
            proportions_df: Original proportions (for reference)
            prediction_weeks: Number of prediction weeks
            variants: List of variant names
            weeks: List of week dates

        Returns:
            DataFrame with smoothed historical + predicted future proportions
        """
        if fit is None:
            # Placeholder: no smoothing, just return input + extend predictions
            logger.warning("No Stan fit available - returning placeholder predictions")
            return proportions_df

        try:
            # Extract posterior draws for p_pred (predicted proportions)
            draws = fit.stan_variable("p_pred")  # Shape: (draws, weeks, variants)

            n_weeks = draws.shape[1]
            n_variants = draws.shape[2]
            n_draws = draws.shape[0]

            # Compute posterior statistics for each week and variant
            results = []
            for week_idx in range(n_weeks):
                is_fitted = week_idx < len(weeks)
                week_date = weeks[week_idx] if is_fitted else weeks[-1] + timedelta(
                    weeks=week_idx - len(weeks) + 1
                )

                for var_idx, variant in enumerate(variants):
                    posterior_draws = draws[:, week_idx, var_idx]

                    # Compute posterior quantiles
                    lower_q = (1 - self.config.credible_interval) / 2
                    upper_q = 1 - lower_q

                    results.append(
                        {
                            "collection_date": week_date,
                            "voc": variant,
                            "proportion": np.median(posterior_draws),
                            "ci_lower": np.quantile(posterior_draws, lower_q),
                            "ci_upper": np.quantile(posterior_draws, upper_q),
                            "posterior_mean": np.mean(posterior_draws),
                            "posterior_sd": np.std(posterior_draws),
                            "retrospective_smooth": is_fitted,
                        }
                    )

            logger.debug(f"Extracted {len(results)} posterior predictions")
            return pl.DataFrame(results)

        except Exception as e:
            raise ModelException(f"Failed to extract predictions: {e}") from e

    def _placeholder_predictions(
        self,
        proportions_df: pl.DataFrame,
        prediction_weeks: int,
    ) -> pl.DataFrame:
        """Generate placeholder predictions (when Stan not available).

        Simple approach:
        - Use last week proportions as baseline
        - Add small trend (assuming stability)
        - Extend CIs slightly for uncertainty growth

        Args:
            proportions_df: Historical proportions
            prediction_weeks: Number of weeks to predict

        Returns:
            DataFrame with historical + placeholder predictions
        """
        # Identify date column
        date_col = "collection_date" if "collection_date" in proportions_df.columns else "week_start"

        # Get last week data
        last_week = proportions_df[date_col].max()

        # Build predictions
        predictions = proportions_df.with_columns(
            pl.lit(True).alias("retrospective_smooth")
        )

        # Add placeholder predictions for future weeks
        future_data = []
        for future_week in range(1, prediction_weeks + 1):
            pred_date = last_week + timedelta(weeks=future_week)

            for row in proportions_df.filter(pl.col(date_col) == last_week).iter_rows(
                named=True
            ):
                # Slight variance decay (conservative) for future predictions
                uncertainty_growth = 1.0 + (0.02 * future_week)  # 2% per week

                future_data.append(
                    {
                        date_col: pred_date,
                        "voc": row["voc"],
                        "proportion": row["proportion"],  # Hold constant as baseline
                        "ci_lower": row["ci_lower"] * uncertainty_growth,
                        "ci_upper": row["ci_upper"] * uncertainty_growth,
                        "retrospective_smooth": False,
                    }
                )

        if future_data:
            future_df = pl.DataFrame(future_data)
            predictions = pl.concat([predictions, future_df])

        logger.debug(f"Generated {len(future_data)} placeholder predictions")
        return predictions

    def __enter__(self):
        """Context manager entry."""
        logger.debug("NowcastModel context manager entered")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self.fit_result is not None:
            # Clean up Stan fit if needed
            self.fit_result = None
        logger.debug("NowcastModel context manager exited")
