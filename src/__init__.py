"""
Phi Modulation Analysis Framework

A scientifically defensible framework for testing φ-modulation
as an empirical pattern in cosmic structure.

This cosmological research work is independent and separate from any other projects.
It originated as a thought experiment and hypothesis exploring whether the Golden Ratio
might be fundamental to cosmic structure, given self-similarity patterns observed
from plants to galaxies.
"""

__version__ = "0.1.0"

# Core modules
from .phi_modulation import PhiModulationModel

# Optional imports (may not be available in all environments)
try:
    from .systematics import SystematicErrorBudget
except ImportError:
    SystematicErrorBudget = None

try:
    from .bayesian_tools import BayesianEvidence, compute_bic, interpret_bic
except ImportError:
    BayesianEvidence = None
    compute_bic = None
    interpret_bic = None

__all__ = [
    'PhiModulationModel',
    'SystematicErrorBudget',
    'BayesianEvidence',
    'compute_bic',
    'interpret_bic',
]

