"""Survey weighting stage for SC2 proportion modeling.

Responsible for calculating inverse probability-of-selection weights based
on surveys and infection prevalence, then computing design-weighted proportions
and confidence intervals using R's survey package.
"""

from dataclasses import dataclass
from typing import Optional

import polars as pl

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
    - Estimated testing volume by region/time
    - Infection prevalence from external sources
    - Sequencing sampling fractions

    Uses R's survey package for design-based inference via rpy2 bridge.

    Handles:
    - Weight calculation and validation
    - Confidence interval computation (KG, Wilson, Agresti-Coull)
    - Regional stratification
    - Edge cases (zero counts, extreme weights)
    """

    def __init__(self, model_config: ModelConfig):
        """Initialize weighter with model configuration.

        Args:
            model_config: ModelConfig specifying CI type and method
        """
        self.model_config = model_config
        self.r_bridge = None
        # TODO: Initialize rpy2 connection to R

    def calculate_weighted_proportions(
        self,
        sequences_df: pl.DataFrame,
        infection_estimates_df: pl.DataFrame,
        testing_volume_df: Optional[pl.DataFrame] = None,
    ) -> WeightingResult:
        """Calculate design-weighted variant proportions.

        Args:
            sequences_df: Aggregated sequences with columns:
                - collection_date, hhs_region, voc, weight (if pre-weighted)
            infection_estimates_df: External infection prevalence data for IPW
            testing_volume_df: Optional testing volume by region/date

        Returns:
            WeightingResult with proportions and confidence intervals

        Raises:
            WeightingException: If weighting calculation fails
        """
        try:
            # Calculate IPW weights
            weights = self._calculate_ipw_weights(
                sequences_df, infection_estimates_df, testing_volume_df
            )

            # Add weights to sequences
            weighted_seq = sequences_df.join(weights, on="hhs_region")

            # Call R survey package for design-based estimates
            result = self._call_r_survey_estimates(weighted_seq)

            return result
        except Exception as e:
            raise WeightingException(f"Failed to calculate weighted proportions: {e}") from e

    def _calculate_ipw_weights(
        self,
        sequences_df: pl.DataFrame,
        infection_estimates_df: pl.DataFrame,
        testing_volume_df: Optional[pl.DataFrame] = None,
    ) -> pl.DataFrame:
        """Calculate inverse probability-of-selection weights.

        Args:
            sequences_df: Sequence data
            infection_estimates_df: External infection estimates by region/date
            testing_volume_df: Optional testing volume adjustments

        Returns:
            DataFrame with weights by hhs_region and time period
        """
        # TODO: Implement IPW calculation following current R implementation
        # Weight_i = 1 / (probability of being sequenced)
        #          = 1 / (testing_burden_i * fraction_sequenced_i)
        pass

    def _call_r_survey_estimates(self, weighted_sequences: pl.DataFrame) -> WeightingResult:
        """Call R survey package to compute design-weighted estimates.

        Args:
            weighted_sequences: Sequences with weights column

        Returns:
            WeightingResult with proportions and confidence intervals
        """
        # TODO: Use rpy2 to call R survey package functions:
        # - svydesign(): Define survey design
        # - svymean(): Calculate proportions
        # - confint(): Generate confidence intervals (type specified in config)
        pass

    def __enter__(self):
        """Context manager entry."""
        # TODO: Initialize R session via rpy2
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        # TODO: Clean up R session
        pass
