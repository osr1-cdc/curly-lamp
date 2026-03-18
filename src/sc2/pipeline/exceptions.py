"""Custom exceptions for SC2 proportion modeling pipeline."""


class SC2PipelineException(Exception):
    """Base exception for all SC2 pipeline errors."""

    pass


class DataFetchException(SC2PipelineException):
    """Raised when data fetching from Impala fails."""

    pass


class ConfigurationException(SC2PipelineException):
    """Raised when configuration is invalid or incomplete."""

    pass


class AggregationException(SC2PipelineException):
    """Raised when lineage aggregation fails."""

    pass


class WeightingException(SC2PipelineException):
    """Raised when survey weighting calculation fails."""

    pass


class ModelException(SC2PipelineException):
    """Raised when statistical model fitting or prediction fails."""

    pass


class ExportException(SC2PipelineException):
    """Raised when exporting results fails."""

    pass
