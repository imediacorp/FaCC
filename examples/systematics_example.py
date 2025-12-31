"""
Example: Using systematic error analysis with φ-modulation forecasts

This example demonstrates how to use the SystematicErrorBudget class
to include systematic errors in DESI forecast calculations.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib.pyplot as plt
from src.phi_modulation import PhiModulationModel
from src.systematics import SystematicErrorBudget

def main():
    """Run example systematic error analysis"""
    
    print("=" * 60)
    print("Systematic Error Analysis Example")
    print("=" * 60)
    
    # Initialize model
    model = PhiModulationModel()
    
    # Run forecast with systematics
    print("\n1. Running DESI forecast with systematic errors...")
    A_phi = 0.01
    result = model.forecast_desi_sensitivity_with_systematics(
        A_phi_true=A_phi,
        k_min=0.01,
        k_max=0.3,
        n_k=50,
        include_systematics=True
    )
    
    # Print results
    print("\n2. Forecast Results:")
    print(f"   Statistical uncertainty: σ_Aφ = {result['sigma_Aphi_stat']:.6f}")
    print(f"   Systematic uncertainty:  σ_Aφ_sys = {result['sigma_Aphi_sys']:.6f}")
    print(f"   Total uncertainty:       σ_Aφ_total = {result['sigma_Aphi_total']:.6f}")
    print(f"   Signal-to-noise ratio:   SNR = {result['SNR']:.2f}σ")
    
    # Systematic error breakdown
    sys_budget = result['systematic_budget']
    print("\n3. Systematic Error Breakdown:")
    print(f"   Photo-z error fraction:     {np.mean(sys_budget['sigma_P_photoz'] / sys_budget['sigma_P_total']):.2%}")
    print(f"   Bias error fraction:        {np.mean(sys_budget['sigma_P_bias'] / sys_budget['sigma_P_total']):.2%}")
    print(f"   Geometry error fraction:    {np.mean(sys_budget['sigma_P_geometry'] / sys_budget['sigma_P_total']):.2%}")
    print(f"   Total systematic fraction:  {np.mean(sys_budget['fraction_sys']):.2%}")
    
    # Compare with statistical-only forecast
    print("\n4. Comparison with statistical-only forecast:")
    result_stat = model.forecast_desi_sensitivity(
        A_phi_true=A_phi,
        k_min=0.01,
        k_max=0.3,
        n_k=50
    )
    sigma_increase = (result['sigma_Aphi_total'] / result_stat['sigma_Aphi'] - 1) * 100
    print(f"   Uncertainty increase: {sigma_increase:.1f}%")
    print(f"   SNR reduction: {(1 - result['SNR'] / result_stat['SNR']) * 100:.1f}%")
    
    # Create visualization
    print("\n5. Creating visualization...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    k = result['k']
    
    # Panel 1: Power spectrum modulation
    ax = axes[0, 0]
    ratio = result['Pk_mod'] / result['Pk_base'] - 1
    ax.semilogx(k, ratio * 100, 'b-', linewidth=2, label='Modulation')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('k [h/Mpc]')
    ax.set_ylabel('ΔP/P [%]')
    ax.set_title('φ-modulation in Power Spectrum')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Panel 2: Error breakdown
    ax = axes[0, 1]
    ax.semilogx(k, result['sigma_P'] * 100 / result['Pk_base'], 
                'k-', linewidth=2, label='Statistical')
    ax.semilogx(k, sys_budget['sigma_P_sys'] * 100 / result['Pk_base'], 
                'r--', linewidth=2, label='Systematic')
    ax.semilogx(k, sys_budget['sigma_P_total'] * 100 / result['Pk_base'], 
                'b:', linewidth=2, label='Total')
    ax.set_xlabel('k [h/Mpc]')
    ax.set_ylabel('Relative Error [%]')
    ax.set_title('Error Breakdown')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Panel 3: Systematic error components
    ax = axes[1, 0]
    ax.semilogx(k, sys_budget['sigma_P_photoz'] * 100 / result['Pk_base'], 
                'r-', linewidth=2, label='Photo-z')
    ax.semilogx(k, sys_budget['sigma_P_bias'] * 100 / result['Pk_base'], 
                'g--', linewidth=2, label='Bias')
    ax.semilogx(k, sys_budget['sigma_P_geometry'] * 100 / result['Pk_base'], 
                'b:', linewidth=2, label='Geometry')
    ax.set_xlabel('k [h/Mpc]')
    ax.set_ylabel('Relative Error [%]')
    ax.set_title('Systematic Error Components')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Panel 4: Systematic error fraction
    ax = axes[1, 1]
    ax.semilogx(k, sys_budget['fraction_sys'] * 100, 'purple', linewidth=2)
    ax.axhline(50, color='k', linestyle='--', alpha=0.3, label='50%')
    ax.set_xlabel('k [h/Mpc]')
    ax.set_ylabel('Systematic Fraction [%]')
    ax.set_title('Fraction of Total Error from Systematics')
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_ylim([0, 100])
    
    plt.tight_layout()
    
    # Save figure
    output_dir = '../figures'
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'systematics_analysis.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"   Figure saved to: {output_path}")
    
    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)

if __name__ == '__main__':
    main()

