"""Tests for configuration schema and loading."""

from datetime import date
from pathlib import Path

import pytest

from sc2.config import SC2Config, ImpalaConfig, RunConfig, VariantConfig


def test_impala_config_defaults():
    """Test ImpalaConfig with default values."""
    config = ImpalaConfig()
    assert config.host == "localhost"
    assert config.port == 21050
    assert config.database == "covid"


def test_run_config_date_parsing():
    """Test RunConfig parses string dates correctly."""
    config = RunConfig(
        data_date="2026-03-18",
        results_tag="TEST",
    )
    assert isinstance(config.data_date, date)
    assert config.data_date == date(2026, 3, 18)


def test_run_config_date_object():
    """Test RunConfig accepts date objects."""
    test_date = date(2026, 3, 18)
    config = RunConfig(
        data_date=test_date,
        results_tag="TEST",
    )
    assert config.data_date == test_date


def test_sc2_config_from_dict():
    """Test SC2Config can be instantiated from dictionary."""
    config_dict = {
        "run": {
            "data_date": "2026-03-18",
            "results_tag": "CDT",
        },
        "variants": {
            "voc1": ["BA.1.1", "BA.2", "XBB"],
        },
    }
    config = SC2Config(**config_dict)
    assert config.run.data_date == date(2026, 3, 18)
    assert len(config.variants.voc1) == 3


@pytest.mark.integration
def test_sc2_config_from_yaml(tmp_path):
    """Test SC2Config loads from YAML file."""
    config_file = tmp_path / "test_config.yml"
    config_file.write_text(
        """
run:
  data_date: "2026-03-18"
  results_tag: "TEST"

variants:
  voc1:
    - BA.1.1
    - BA.2
"""
    )

    config = SC2Config.from_yaml(config_file)
    assert config.run.data_date == date(2026, 3, 18)
    assert config.run.results_tag == "TEST"
    assert len(config.variants.voc1) == 2


def test_sc2_config_to_yaml(tmp_path):
    """Test SC2Config exports to YAML file."""
    config = SC2Config(
        run=RunConfig(
            data_date="2026-03-18",
            results_tag="CDT",
        ),
        variants=VariantConfig(
            voc1=["BA.1.1", "BA.2"],
        ),
    )

    output_file = tmp_path / "output.yml"
    config.to_yaml(output_file)
    assert output_file.exists()

    # Reload and verify
    config2 = SC2Config.from_yaml(output_file)
    assert config2.run.data_date == config.run.data_date
    assert config2.run.results_tag == config.run.results_tag


def test_variant_config_with_reduced():
    """Test VariantConfig with both voc1 and voc1_reduced."""
    config = VariantConfig(
        voc1=["BA.1.1", "BA.2", "XBB"],
        voc1_reduced=["B.1.1.529", "B.1.617.2"],
    )
    assert len(config.voc1) == 3
    assert len(config.voc1_reduced) == 2
