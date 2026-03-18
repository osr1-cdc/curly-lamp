"""Pytest configuration and shared fixtures."""

import pytest
from datetime import date
from pathlib import Path

from sc2.config import SC2Config, RunConfig, VariantConfig


@pytest.fixture
def sample_config() -> SC2Config:
    """Fixture providing a sample SC2Config for testing."""
    return SC2Config(
        run=RunConfig(
            data_date=date(2026, 3, 18),
            results_tag="TEST",
            use_previously_imported_data=False,
        ),
        variants=VariantConfig(
            voc1=[
                "BA.1.1",
                "BA.2",
                "BA.2.75",
                "XBB",
                "XBB.1.5",
                "JN.1",
                "KP.3",
                "XEC",
            ],
            voc1_reduced=["B.1.1.529", "B.1.617.2"],
        ),
    )


@pytest.fixture
def tmp_results_dir(tmp_path) -> Path:
    """Fixture providing temporary results directory."""
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    return results_dir
