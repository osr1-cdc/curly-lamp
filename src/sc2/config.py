"""Configuration schema and loading for SC2 proportion modeling.

Uses pydantic v2 for validation and yaml for configuration files.
"""

from datetime import date, datetime
from pathlib import Path
from typing import Literal, Optional

from pydantic import BaseModel, Field, validator
from pydantic_settings import BaseSettings


class ImpalaConfig(BaseModel):
    """Impala database connection configuration."""

    host: str = Field(default="localhost", description="Impala server hostname")
    port: int = Field(default=21050, description="Impala server port")
    database: str = Field(default="covid", description="Default database name")
    timeout_seconds: int = Field(default=300, description="Query timeout in seconds")


class DatabaseConfig(BaseModel):
    """Database configuration."""

    impala: ImpalaConfig = Field(default_factory=ImpalaConfig)


class VariantConfig(BaseModel):
    """Variant tracking configuration."""

    voc1: list[str] = Field(
        ...,
        description="List of variants of concern/interest to track",
    )
    voc1_reduced: Optional[list[str]] = Field(
        default=None,
        description="Reduced set of variants for comparison",
    )


class ModelConfig(BaseModel):
    """Statistical model configuration."""

    weeks_lookback: int = Field(default=8, description="Weeks of historical data for fitting")
    weeks_predict: int = Field(default=4, description="Weeks ahead to predict")
    nowcast_method: Literal["multinomial_regression", "stan_bayesian"] = Field(
        default="multinomial_regression",
        description="Nowcasting methodology",
    )
    credible_interval: float = Field(
        default=0.95,
        description="Credible interval for uncertainty quantification",
    )
    ci_type: Literal["KG", "wilson", "agresti-coull"] = Field(
        default="KG",
        description="Confidence interval type for proportions",
    )


class RunConfig(BaseModel):
    """Run-specific settings."""

    data_date: date = Field(..., description="Data cutoff date (YYYY-MM-DD format)")
    results_tag: str = Field(default="CDT", description="Results folder tag")
    use_previously_imported_data: bool = Field(
        default=False,
        description="Use cached data instead of fresh database query",
    )

    @validator("data_date", pre=True)
    def parse_date(cls, v):
        if isinstance(v, str):
            return datetime.strptime(v, "%Y-%m-%d").date()
        return v


class OutputConfig(BaseModel):
    """Output and caching configuration."""

    results_dir: Path = Field(default=Path("results"), description="Results output directory")
    cache_dir: Path = Field(default=Path("data/cache"), description="Intermediate data cache")
    export_formats: list[Literal["csv", "json", "parquet", "png", "pdf"]] = Field(
        default=["csv", "json", "png"],
        description="Output formats for results",
    )


class SC2Config(BaseSettings):
    """Root configuration object combining all settings."""

    run: RunConfig = Field(..., description="Run-specific settings")
    database: DatabaseConfig = Field(default_factory=DatabaseConfig)
    variants: VariantConfig = Field(..., description="Variant definitions")
    model: ModelConfig = Field(default_factory=ModelConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)

    class Config:
        """Pydantic settings configuration."""

        env_file = ".env"
        case_sensitive = False
        yaml_file = None

    @classmethod
    def from_yaml(cls, path: Path) -> "SC2Config":
        """Load configuration from YAML file.

        Args:
            path: Path to YAML configuration file

        Returns:
            SC2Config instance

        Raises:
            FileNotFoundError: If YAML file does not exist
            ValueError: If YAML is invalid
        """
        import yaml

        if not path.exists():
            raise FileNotFoundError(f"Configuration file not found: {path}")

        with open(path) as f:
            data = yaml.safe_load(f)

        if not isinstance(data, dict):
            raise ValueError(f"Invalid YAML configuration: {path}")

        return cls(**data)

    @classmethod
    def from_env(cls) -> "SC2Config":
        """Load configuration from environment variables and .env file.

        Returns:
            SC2Config instance
        """
        return cls()

    def to_yaml(self, path: Path) -> None:
        """Export configuration to YAML file.

        Args:
            path: Output path for YAML file
        """
        import yaml

        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            yaml.dump(self.model_dump(exclude_none=True), f, default_flow_style=False)


# Example configuration for development and testing
EXAMPLE_CONFIG = {
    "run": {
        "data_date": "2026-03-18",
        "results_tag": "CDT",
        "use_previously_imported_data": False,
    },
    "database": {
        "impala": {
            "host": "impala-host",
            "port": 21050,
            "database": "covid",
        }
    },
    "variants": {
        "voc1": [
            "BA.1.1",
            "BA.2",
            "BA.2.12.1",
            "BA.2.75",
            "BA.2.86",
            "XBB",
            "XBB.1.5",
            "JN.1",
            "KP.3",
            "XEC",
        ],
        "voc1_reduced": ["B.1.1.529", "B.1.617.2"],
    },
    "model": {
        "weeks_lookback": 8,
        "weeks_predict": 4,
        "nowcast_method": "multinomial_regression",
        "credible_interval": 0.95,
        "ci_type": "KG",
    },
    "output": {
        "results_dir": "results",
        "cache_dir": "data/cache",
        "export_formats": ["csv", "json", "png"],
    },
}
