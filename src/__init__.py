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

try:
    from .hz_analysis import run_hz_analysis, fit_forward_branch, fit_reverse_branch
    from .cmb_analysis import run_cmb_analysis, fit_dual_oscillations
    from .lss_analysis import run_lss_analysis, compute_phi_scales
except ImportError:
    run_hz_analysis = None
    run_cmb_analysis = None
    run_lss_analysis = None

try:
    from .g_phi_analysis import run_g_phi_analysis
except ImportError:
    run_g_phi_analysis = None

try:
    from .desi_bao import run_desi_bao_analysis, fetch_desi_bao
except ImportError:
    run_desi_bao_analysis = None
    fetch_desi_bao = None

__all__ = [
    'PhiModulationModel',
    'SystematicErrorBudget',
    'BayesianEvidence',
    'compute_bic',
    'interpret_bic',
    'run_hz_analysis',
    'run_cmb_analysis',
    'run_lss_analysis',
    'compute_phi_scales',
    'run_g_phi_analysis',
]

