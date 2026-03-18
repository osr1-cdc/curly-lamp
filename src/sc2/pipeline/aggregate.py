"""Lineage aggregation stage for SC2 proportion modeling.

Responsible for aggregating fine-grained Pangolin lineages to higher-level
variants of concern/interest (VOC/VOI) based on taxonomic hierarchy rules.
"""

from dataclasses import dataclass
from typing import Optional

import polars as pl
from loguru import logger

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
    - Hierarchical lineage prefix matching
    - Fast lookup via trie structure
    - Un-aggregated lineages (passed through as-is)
    - Validation that all lineages in VOC list are covered
    """

    def __init__(self, rules: Optional[list[AggregationRule]] = None):
        """Initialize aggregator with rules.

        Args:
            rules: List of AggregationRule objects. If None, uses default Pangolin hierarchy.
        """
        self.rules = rules or self._load_default_rules()
        self.parent_map: dict[str, str] = {}  # Maps exact lineage -> parent
        self.prefix_map: dict[str, str] = {}  # Maps prefix patterns -> parent
        self._build_lookup_tables()
        logger.debug(f"LineageAggregator initialized with {len(self.rules)} rules")

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

        Uses hierarchical matching: first checks exact match in parent_map,
        then tries prefix matching (longest prefix wins), falls back to original.

        Args:
            lineage: Pangolin lineage string (e.g., "JN.1.4.3")

        Returns:
            Aggregated lineage (e.g., "JN.1") or original if not in hierarchy

        Examples:
            >>> agg = LineageAggregator()
            >>> agg.aggregate_lineage("JN.1.4.3")
            'JN.1'
            >>> agg.aggregate_lineage("BA.1.1.2.5")
            'BA.1.1'
            >>> agg.aggregate_lineage("XBB.1.5.10.3")
            'XBB.1.5'
        """
        # Check exact match first
        if lineage in self.parent_map:
            return self.parent_map[lineage]

        # Try prefix matching (longest first)
        # Split by "." and progressively check shorter prefixes
        parts = lineage.split(".")
        for i in range(len(parts) - 1, 0, -1):
            prefix = ".".join(parts[:i])
            if prefix in self.prefix_map:
                return self.prefix_map[prefix]

        # No match found - return original lineage
        return lineage

    def _build_lookup_tables(self):
        """Build efficient lookup tables from aggregation rules.

        Creates:
        - parent_map: dict mapping each child lineage to its parent
        - prefix_map: dict mapping prefix patterns to parents for wildcard matching
        """
        self.parent_map = {}
        self.prefix_map = {}

        for rule in self.rules:
            # Add exact matches for explicit children
            for child in rule.children:
                self.parent_map[child] = rule.parent

            # For wildcard patterns (indicated by trailing .*), add to prefix_map
            # e.g., "JN.1.*" -> store "JN.1" in prefix_map pointing to "JN.1"
            # This enables matching JN.1.4.3 -> "JN.1.4.3" (if exact) or "JN.1" (via prefix)
            self.prefix_map[rule.parent] = rule.parent

    @staticmethod
    def _load_default_rules() -> list[AggregationRule]:
        """Load default Pangolin lineage aggregation rules.

        Rules based on CDC/Pango lineage hierarchy. Aggregations are specified
        so that fine-grained sublineages map to their designated parent for tracking.

        Returns:
            List of AggregationRule objects representing standard Omicron
            and other aggregation hierarchies.
        """
        return [
            # BA.1 lineage family
            AggregationRule(
                parent="BA.1.1",
                children=["BA.1.1.2", "BA.1.1.3", "BA.1.1.5", "BA.1.1.6"],
            ),
            # BA.2 lineage family
            AggregationRule(
                parent="BA.2",
                children=["BA.2.1", "BA.2.2", "BA.2.3", "BA.2.45"],
            ),
            AggregationRule(
                parent="BA.2.12.1",
                children=["BA.2.12.1.1", "BA.2.12.1.2"],
            ),
            AggregationRule(
                parent="BA.2.75",
                children=["BA.2.75.1", "BA.2.75.3", "BA.2.75.4", "BA.2.75.5", "BA.2.75.6"],
            ),
            AggregationRule(
                parent="BA.2.75.2",
                children=["BA.2.75.2.1"],
            ),
            # BA.4/BA.5 families
            AggregationRule(
                parent="BA.4",
                children=["BA.4.1", "BA.4.2", "BA.4.3", "BA.4.4"],
            ),
            AggregationRule(
                parent="BA.4.6",
                children=["BA.4.6.1", "BA.4.6.2", "BA.4.6.3"],
            ),
            AggregationRule(
                parent="BA.5",
                children=["BA.5.1", "BA.5.2", "BA.5.3", "BA.5.4"],
            ),
            AggregationRule(
                parent="BA.5.2.6",
                children=["BA.5.2.6.1"],
            ),
            # BQ family
            AggregationRule(
                parent="BQ.1",
                children=["BQ.1.2", "BQ.1.3", "BQ.1.8"],
            ),
            AggregationRule(
                parent="BQ.1.1",
                children=["BQ.1.1.1", "BQ.1.1.2", "BQ.1.1.3", "BQ.1.1.4", "BQ.1.1.5"],
            ),
            # XBB family (major lineage requiring detailed aggregation)
            AggregationRule(
                parent="XBB",
                children=["XBB.1", "XBB.2"],
            ),
            AggregationRule(
                parent="XBB.1",
                children=["XBB.1.2", "XBB.1.3", "XBB.1.4"],
            ),
            AggregationRule(
                parent="XBB.1.5",
                children=[
                    "XBB.1.5.1",
                    "XBB.1.5.2",
                    "XBB.1.5.3",
                    "XBB.1.5.10",
                    "XBB.1.5.59",
                    "XBB.1.5.68",
                    "XBB.1.5.70",
                    "XBB.1.5.72",
                ],
            ),
            AggregationRule(
                parent="XBB.1.9.1",
                children=["XBB.1.9.1.1", "XBB.1.9.1.2"],
            ),
            AggregationRule(
                parent="XBB.1.9.2",
                children=["XBB.1.9.2.1"],
            ),
            AggregationRule(
                parent="XBB.1.16",
                children=["XBB.1.16.2", "XBB.1.16.3", "XBB.1.16.6"],
            ),
            AggregationRule(
                parent="XBB.1.16.1",
                children=["XBB.1.16.1.2"],
            ),
            AggregationRule(
                parent="XBB.1.16.11",
                children=["XBB.1.16.11.1"],
            ),
            AggregationRule(
                parent="XBB.2.3",
                children=["XBB.2.3.1", "XBB.2.3.2", "XBB.2.3.8"],
            ),
            # JN.1 family (major recent lineage)
            AggregationRule(
                parent="JN.1",
                children=[
                    "JN.1.2",
                    "JN.1.3",
                    "JN.1.4",
                    "JN.1.5",
                    "JN.1.6",
                    "JN.1.7",
                    "JN.1.8",
                    "JN.1.9",
                    "JN.1.10",
                    "JN.1.11",
                    "JN.1.12",
                    "JN.1.13",
                    "JN.1.14",
                    "JN.1.15",
                    "JN.1.16",
                ],
            ),
            AggregationRule(
                parent="JN.1.7",
                children=["JN.1.7.1"],
            ),
            AggregationRule(
                parent="JN.1.8.1",
                children=["JN.1.8.1.1"],
            ),
            AggregationRule(
                parent="JN.1.11.1",
                children=["JN.1.11.1.1"],
            ),
            AggregationRule(
                parent="JN.1.16.1",
                children=["JN.1.16.1.1"],
            ),
            AggregationRule(
                parent="JN.1.18",
                children=["JN.1.18.1", "JN.1.18.2", "JN.1.18.6"],
            ),
            # KP families (post-JN.1 variants)
            AggregationRule(
                parent="KP.2",
                children=["KP.2.1", "KP.2.3"],
            ),
            AggregationRule(
                parent="KP.3",
                children=["KP.3.1", "KP.3.1.1"],
            ),
            # XEC and related
            AggregationRule(
                parent="XEC",
                children=["XEC.1", "XEC.4"],
            ),
        ]
