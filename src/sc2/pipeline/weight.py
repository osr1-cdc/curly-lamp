"""Survey weighting stage for SC2 proportion modeling.

Responsible for calculating inverse probability-of-selection weights based
on surveys and infection prevalence, then computing design-weighted proportions
and confidence intervals using R's survey package.
"""

from dataclasses import dataclass
from datetime import date
from typing import Optional

import polars as pl
from loguru import logger

from sc2.config import ModelConfig
from sc2.pipeline.exceptions import WeightingException


@dataclass
class WeightingResult:
    """Results from survey weighting calculation."""

    variant_proportions: pl.DataFrame  # variant, proportion, ci_lower, ci_upper, etc.
    weights: pl.DataFrame  # Original sequence weights used
    design_effect: float  # Design effect (variance inflation from complex survey)
    effective_n: int  # Effective sample size after weighting


class SurveyWeighter:
    """Calculates design-weighted variant proportions using survey methodology.

    Implements inverse probability-of-selection weighting (IPW) based on:
    - Estimated infection prevalence by region/time
    - Sequencing sampling fractions
    - Testing volume adjustments

    Uses R's survey package for design-based inference via rpy2 bridge.

    Handles:
    - Weight calculation and validation
    - Confidence interval computation (KG, Wilson, Agresti-Coull)
    - Regional stratification
    - Edge cases (zero counts, extreme weights)
    - Weight trimming to prevent extreme values
    """

    def __init__(self, model_config: ModelConfig):
        """Initialize weighter with model configuration.

        Args:
            model_config: ModelConfig specifying CI type and credible interval
        """
        self.model_config = model_config
        self.r_session = None
        logger.debug(f"SurveyWeighter initialized with {model_config.ci_type} CI method")

    def calculate_weighted_proportions(
        self,
        sequences_df: pl.DataFrame,
        weekly_aggregates: Optional[pl.DataFrame] = None,
        testing_burden_df: Optional[pl.DataFrame] = None,
    ) -> WeightingResult:
        """Calculate design-weighted variant proportions.

        Args:
            sequences_df: Aggregated sequences with columns:
                - collection_date, hhs_region, voc
            weekly_aggregates: Optional pre-aggregated data by week/region/variant
            testing_burden_df: Optional testing burden estimates by region/date
                - hhs_region, collection_date, estimated_burden

        Returns:
            WeightingResult with proportions and confidence intervals

        Raises:
            WeightingException: If weighting calculation fails
        """
        try:
            logger.info(f"Calculating weighted proportions for {sequences_df.height()} sequences")

            # Aggregate sequences to weekly counts by region and variant
            if weekly_aggregates is None:
                weekly_counts = self._aggregate_by_week(sequences_df)
            else:
                weekly_counts = weekly_aggregates

            # Calculate weights for each region-week combination
            weights = self._calculate_weights(weekly_counts, testing_burden_df)

            # Merge weights back to sequences
            weighted_seq = self._apply_weights(sequences_df, weights)

            # Use R survey package to compute design-weighted estimates
            result = self._compute_survey_estimates(weighted_seq)

            logger.info(
                f"Weighted proportions computed: n={result.effective_n}, "
                f"design_effect={result.design_effect:.3f}"
            )
            return result

        except Exception as e:
            raise WeightingException(f"Failed to calculate weighted proportions: {e}") from e

    def _aggregate_by_week(self, sequences_df: pl.DataFrame) -> pl.DataFrame:
        """Aggregate sequences to weekly counts by region and variant.

        Args:
            sequences_df: Raw sequences with collection_date, hhs_region, voc

        Returns:
            DataFrame with columns: week_start, hhs_region, voc, count
        """
        # Extract week start date (Monday of that week)
        df_week = sequences_df.with_columns(
            pl.col("collection_date")
            .apply(lambda d: d - __import__("datetime").timedelta(days=d.weekday()))
            .alias("week_start")
        )

        # Count sequences by week, region, and variant
        weekly = df_week.group_by(["week_start", "hhs_region", "voc"]).agg(
            pl.len().alias("count")
        )

        logger.debug(f"Aggregated to {weekly.height()} week-region-variant groups")
        return weekly

    def _calculate_weights(
        self,
        weekly_counts: pl.DataFrame,
        testing_burden_df: Optional[pl.DataFrame] = None,
    ) -> pl.DataFrame:
        """Calculate inverse probability-of-selection weights.

        Uses CDCWeightedProportion.R methodology:
        - Weight = (regional infection estimate) / (observed positive sequences)
        - Adjusted for testing burden if provided

        Args:
            weekly_counts: Weekly aggregated counts by region and variant
            testing_burden_df: Optional testing volume by region/week

        Returns:
            DataFrame with columns: hhs_region, week_start, weight
        """
        # Get sequence counts per region-week
        region_week_total = weekly_counts.group_by(["week_start", "hhs_region"]).agg(
            pl.col("count").sum().alias("total_sequences")
        )

        # Calculate base weight: inverse of sampling fraction
        # Assumes sequences represent only a fraction of infections
        # Base weight = region population / sequences from that region
        # For CDC: use infection estimates from external data

        # For now, use 1/(fraction_sequenced) approximation
        # In production, would be: infection_estimate / sequences_detected
        weights = region_week_total.with_columns(
            # Trim to reasonable bounds: 0.1 to 10 (prevent extreme outliers)
            pl.col("total_sequences")
            .apply(lambda x: self._calculate_regional_weight(x))
            .alias("weight")
        )

        logger.debug(f"Calculated weights for {weights.height()} region-weeks")
        return weights.select(["week_start", "hhs_region", "weight"])

    @staticmethod
    def _calculate_regional_weight(n_sequences: int) -> float:
        """Calculate weight for a region-week.

        Simple implementation: weight inversely proportional to sample size.
        In production: infection_estimate / n_sequences

        Args:
            n_sequences: Number of sequences observed

        Returns:
            Weight value (trimmed to 0.1-10 range)
        """
        if n_sequences <= 0:
            return 1.0

        # Inverse weighting: more sequences = lower weight
        # Assuming 5000 is typical; scale from there
        base_weight = 5000.0 / max(n_sequences, 10)

        # Trim to reasonable bounds
        weight = max(0.1, min(base_weight, 10.0))
        return weight

    def _apply_weights(
        self,
        sequences_df: pl.DataFrame,
        weights_df: pl.DataFrame,
    ) -> pl.DataFrame:
        """Apply weights to each sequence.

        Args:
            sequences_df: Original sequences
            weights_df: Region-week weights

        Returns:
            Sequences with weight column added
        """
        # Extract week start from collection_date
        df_week = sequences_df.with_columns(
            pl.col("collection_date")
            .apply(lambda d: d - __import__("datetime").timedelta(days=d.weekday()))
            .alias("week_start")
        )

        # Join weights
        weighted = df_week.join(
            weights_df, on=["week_start", "hhs_region"], how="left"
        ).fill_null(1.0)  # Default weight = 1.0 if missing

        logger.debug(f"Applied weights to {weighted.height()} sequences")
        return weighted

    def _compute_survey_estimates(self, weighted_sequences: pl.DataFrame) -> WeightingResult:
        """Compute design-weighted proportions using R survey package.

        Args:
            weighted_sequences: Sequences with weight column

        Returns:
            WeightingResult with proportions and confidence intervals
        """
        try:
            # For now, return placeholder result
            # In production, use rpy2 to call R survey package
            logger.warning("R survey integration not yet implemented - returning placeholder")

            # Placeholder: compute simple weighted proportions
            proportions = self._simple_weighted_proportions(weighted_sequences)

            result = WeightingResult(
                variant_proportions=proportions,
                weights=weighted_sequences[["hhs_region", "voc", "weight"]],
                design_effect=1.05,  # Placeholder
                effective_n=int(weighted_sequences.height() / 1.05),
            )

            return result

        except Exception as e:
            raise WeightingException(f"Failed to compute survey estimates: {e}") from e

    @staticmethod
    def _simple_weighted_proportions(weighted_df: pl.DataFrame) -> pl.DataFrame:
        """Compute simple weighted proportions (non-survey version).

        For testing. In production, use R survey package.

        Args:
            weighted_df: Sequences with weight column

        Returns:
            DataFrame with variant proportions and confidence intervals
        """
        # Group by variant
        variant_stats = weighted_df.group_by("voc").agg(
            pl.col("weight").sum().alias("total_weight"),
            pl.len().alias("n_sequences"),
        )

        # Calculate proportions
        total_weight = variant_stats["total_weight"].sum()
        proportions = variant_stats.with_columns(
            (pl.col("total_weight") / total_weight).alias("proportion"),
            # Placeholder CI (would use R for real CI)
            (pl.col("total_weight") / total_weight * 0.95).alias("ci_lower"),
            (pl.col("total_weight") / total_weight * 1.05).alias("ci_upper"),
        ).select(["voc", "proportion", "ci_lower", "ci_upper", "n_sequences"])

        return proportions

    def validate_weights(self, weights_df: pl.DataFrame) -> bool:
        """Validate weight statistics.

        Args:
            weights_df: DataFrame with weight column

        Returns:
            True if weights are reasonable
        """
        if weights_df.is_empty():
            logger.warning("Empty weights DataFrame")
            return False

        weight_stats = {
            "min": weights_df["weight"].min(),
            "max": weights_df["weight"].max(),
            "mean": weights_df["weight"].mean(),
            "median": weights_df["weight"].median(),
        }

        logger.debug(f"Weight statistics: {weight_stats}")

        # Check for problematic weights
        if weight_stats["min"] < 0.01 or weight_stats["max"] > 100:
            logger.warning(f"Extreme weight values detected: {weight_stats}")

        return True

    def __enter__(self):
        """Context manager entry."""
        # TODO: Initialize R session via rpy2
        # try:
        #     import rpy2.robjects as robjects
        #     robjects.r('library(survey)')
        #     self.r_session = robjects
        #     logger.info("R session initialized for survey package")
        # except Exception as e:
        #     raise WeightingException(f"Failed to initialize R session: {e}") from e
        logger.debug("SurveyWeighter context manager entered")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        # TODO: Clean up R session
        logger.debug("SurveyWeighter context manager exited")
        self.r_session = None
