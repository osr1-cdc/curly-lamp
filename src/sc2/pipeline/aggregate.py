"""Lineage aggregation stage for SC2 proportion modeling.

Responsible for aggregating fine-grained Pangolin lineages to higher-level
variants of concern/interest (VOC/VOI) based on taxonomic hierarchy rules.
"""

from dataclasses import dataclass
from typing import Optional

import polars as pl

from sc2.pipeline.exceptions import AggregationException


@dataclass
class AggregationRule:
    """Single aggregation rule mapping child lineages to parent."""

    parent: str
    children: list[str]
    description: Optional[str] = None


class LineageAggregator:
    """Aggregates Pangolin lineages to variants of concern.

    Uses Pangolin extended lineage hierarchy to deterministically map
    fine-grained lineages (e.g., JN.1.4.3) to VOC/VOI groups (e.g., JN.1).

    Handles:
    - Hierarchical lineage rules
    - Un-aggregated lineages (passed through as-is)
    - Validation that all lineages in VOC list are covered
    """

    def __init__(self, rules: Optional[list[AggregationRule]] = None):
        """Initialize aggregator with rules.

        Args:
            rules: List of AggregationRule objects. If None, uses default Pangolin hierarchy.
        """
        self.rules = rules or self._load_default_rules()
        self._build_lookup_tables()

    def aggregate_dataframe(
        self,
        df: pl.DataFrame,
        lineage_column: str = "pango_lineage",
        voc_list: Optional[list[str]] = None,
    ) -> pl.DataFrame:
        """Aggregate lineages in a DataFrame.

        Args:
            df: Input DataFrame containing lineage column
            lineage_column: Name of column containing Pangolin lineages
            voc_list: Optional list of VOC to enforce. Lineages not mapping to VOC
                     become "other".

        Returns:
            DataFrame with original lineage_column replaced by aggregated "voc" column

        Raises:
            AggregationException: If aggregation rules are invalid
        """
        try:
            df_agg = df.with_columns(
                pl.col(lineage_column)
                .map_elements(self.aggregate_lineage, return_dtype=pl.Utf8)
                .alias("voc")
            )

            if voc_list:
                df_agg = df_agg.with_columns(
                    pl.when(pl.col("voc").is_in(voc_list))
                    .then(pl.col("voc"))
                    .otherwise(pl.lit("other"))
                    .alias("voc")
                )

            return df_agg.drop(lineage_column)
        except Exception as e:
            raise AggregationException(f"Failed to aggregate lineages: {e}") from e

    def aggregate_lineage(self, lineage: str) -> str:
        """Aggregate single lineage to VOC.

        Args:
            lineage: Pangolin lineage string (e.g., "JN.1.4.3")

        Returns:
            Aggregated lineage (e.g., "JN.1") or original if not in hierarchy
        """
        # TODO: Implement hierarchical matching
        # Try to find longest matching prefix in aggregation rules
        # Fall back to original lineage if no match
        pass

    def _build_lookup_tables(self):
        """Build efficient lookup tables from aggregation rules."""
        # TODO: Create:
        # - parent_map: dict mapping each lineage to its parent
        # - prefix_tree: trie-like structure for efficient matching
        pass

    @staticmethod
    def _load_default_rules() -> list[AggregationRule]:
        """Load default Pangolin lineage aggregation rules.

        Returns:
            List of AggregationRule objects representing standard Omicron
            and other aggregation hierarchies.
        """
        # TODO: Define rules based on Pangolin lineage hierarchy
        # Examples:
        # - BA.1.1.* -> BA.1.1
        # - BA.2.* -> BA.2 (except specific sublineages like BA.2.75 that stay distinct)
        # - XBB.1.5.* -> XBB.1.5
        # - JN.1.* -> JN.1
        return []
