"""
Experimental Data Models for CryoProtect.

This package provides models for experimental data in the CryoProtect system,
including protocols, experiments, results, and time series data.
"""

from .base_model import (
    BaseModel, 
    Uncertainty, 
    Provenance, 
    ValidationError, 
    DataIntegrityError,
    ModelError
)

from .experimental_models import (
    Protocol,
    ProtocolStep,
    ExperimentType,
    TissueType,
    Experiment,
    ExperimentResult,
    TimeSeries,
    TimeSeriesDataPoint
)

__all__ = [
    'BaseModel',
    'Uncertainty',
    'Provenance',
    'ValidationError',
    'DataIntegrityError',
    'ModelError',
    'Protocol',
    'ProtocolStep',
    'ExperimentType',
    'TissueType',
    'Experiment',
    'ExperimentResult',
    'TimeSeries',
    'TimeSeriesDataPoint'
]