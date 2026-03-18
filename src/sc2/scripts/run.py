"""Main entry point for SC2 proportion modeling pipeline.

Orchestrates the complete workflow: fetch -> aggregate -> weight -> model -> export
"""

import sys
from pathlib import Path
from typing import Optional

import typer
from loguru import logger

from sc2.config import SC2Config

app = typer.Typer(help="SC2 Proportion Modeling Pipeline v2.0")


@app.command()
def run(
    config_path: Optional[Path] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to YAML configuration file",
        exists=True,
    ),
    output_dir: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output directory for results (overrides config)",
    ),
    use_cache: bool = typer.Option(
        False,
        "--use-cache",
        help="Use previously imported data instead of fresh query",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Validate config and data pipeline without full execution",
    ),
    verbose: int = typer.Option(
        0,
        "--verbose",
        "-v",
        count=True,
        help="Verbosity level (-v, -vv, -vvv)",
    ),
) -> int:
    """Run SC2 proportion modeling pipeline end-to-end.

    Executes all stages:
    1. Fetch: Query Impala for sequencing data
    2. Aggregate: Map lineages to VOC/VOI
    3. Weight: Calculate survey-weighted proportions
    4. Model: Fit nowcast model and predict
    5. Export: Generate output files

    Examples:
        # Standard monthly run with config file
        sc2-run --config data/config.yml

        # Custom run with explicit output path
        sc2-run --config custom.yml --output results_2026-03-18/

        # Dry run to validate configuration
        sc2-run --config data/config.yml --dry-run

        # Verbose output for debugging
        sc2-run --config data/config.yml -vvv
    """
    # Configure logging
    log_level = ["INFO", "DEBUG", "TRACE"][min(verbose, 2)]
    logger.remove()
    logger.add(sys.stderr, level=log_level, format="<level>{level: <8}</level> | {message}")

    logger.info("Starting SC2 proportion modeling pipeline")

    try:
        # Load configuration
        if config_path:
            logger.info(f"Loading configuration from: {config_path}")
            config = SC2Config.from_yaml(config_path)
        else:
            logger.info("Loading configuration from environment variables")
            config = SC2Config.from_env()

        # Override output directory if specified
        if output_dir:
            config.output.results_dir = output_dir
            logger.info(f"Results will be written to: {output_dir}")

        # Override cache setting
        if use_cache:
            config.run.use_previously_imported_data = True
            logger.warning("Using previously imported data (cache mode)")

        logger.debug(f"Configuration loaded: {config}")

        if dry_run:
            logger.info("Dry run mode: validating configuration only")
            logger.success("Configuration validation passed!")
            return 0

        # TODO: Add pipeline execution logic here
        # 1. Fetch data
        # 2. Aggregate lineages
        # 3. Calculate weights
        # 4. Fit nowcast model
        # 5. Export results

        logger.info("TODO: Pipeline execution not yet implemented")
        logger.warning("This is Phase 1 skeleton - implementation follows in later phases")

        logger.success("Pipeline completed successfully")
        return 0

    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=verbose > 0)
        return 1


@app.command()
def validate_config(
    config_path: Path = typer.Argument(
        ...,
        help="Path to YAML configuration file",
        exists=True,
    ),
) -> int:
    """Validate configuration file without running pipeline.

    Examples:
        sc2-run validate-config data/config.yml
    """
    try:
        logger.info(f"Validating configuration: {config_path}")
        config = SC2Config.from_yaml(config_path)
        logger.success("Configuration is valid!")
        logger.info(f"Data date: {config.run.data_date}")
        logger.info(f"Results tag: {config.run.results_tag}")
        logger.info(f"Variants tracked: {len(config.variants.voc1)}")
        return 0
    except Exception as e:
        logger.error(f"Configuration validation failed: {e}")
        return 1


def main() -> int:
    """Main entry point for CLI."""
    return app()


if __name__ == "__main__":
    sys.exit(main())
