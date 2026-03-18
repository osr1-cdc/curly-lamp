"""Pipeline module for SC2 proportion modeling.

Contains data processing stages: fetch, aggregate, weight, model, export.
Each stage is designed to be independently testable and composable.
"""

__all__ = ["fetch", "aggregate", "weight", "model", "export", "exceptions"]
