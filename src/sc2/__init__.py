"""SC2 Proportion Modeling v2.0 - Modern Python Implementation

A modular, testable, and maintainable Python pipeline for CDC SARS-CoV-2
proportion modeling with statistical nowcasting using survey design-based
weighting and multinomial regression analysis.

Modules:
    config: Configuration schema and loading
    pipeline: Data processing pipeline modules
    models: Statistical models (Stan, etc.)
    scripts: Entry point CLI commands
"""

__version__ = "2.0.0"
__author__ = "CDC Emerging Infectious Diseases Branch"
__all__ = ["config", "pipeline", "models"]
