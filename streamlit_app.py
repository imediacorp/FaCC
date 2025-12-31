"""
Streamlit Dashboard for φ-Modulation Analysis

Interactive dashboard for exploring φ-modulated power spectrum forecasts
and systematic error analysis for DESI surveys.

Run with: streamlit run streamlit_app.py
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

try:
    from phi_modulation import PhiModulationModel
    HAS_PHI_MODULATION = True
except ImportError:
    HAS_PHI_MODULATION = False
    st.error("⚠️ PhiModulationModel not found. Please ensure src/phi_modulation.py is available.")

try:
    from systematics import SystematicErrorBudget
    HAS_SYSTEMATICS = True
except ImportError:
    HAS_SYSTEMATICS = False
    st.warning("⚠️ SystematicErrorBudget not available. Some features will be disabled.")

# Page configuration
st.set_page_config(
    page_title="φ-Modulation Analysis Dashboard",
    page_icon="🌌",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.5rem;
        color: #666;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-container {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
</style>
""", unsafe_allow_html=True)

# Title
st.markdown('<p class="main-header">φ-Modulation Analysis Dashboard</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Interactive DESI Forecast Explorer</p>', unsafe_allow_html=True)

# Sidebar
st.sidebar.header("⚙️ Parameters")

# Initialize model
if HAS_PHI_MODULATION:
    model = PhiModulationModel()
    
    st.sidebar.markdown("### Golden Ratio")
    st.sidebar.info(f"φ = {model.phi:.8f}")
    st.sidebar.info(f"ln(φ) = {model.lnphi:.8f}")
    
    st.sidebar.markdown("---")
    
    # Parameter controls
    st.sidebar.markdown("### Forecast Parameters")
    
    A_phi = st.sidebar.slider(
        "A_φ (Modulation Amplitude)",
        min_value=0.001,
        max_value=0.05,
        value=0.01,
        step=0.001,
        help="Amplitude of φ-modulation"
    )
    
    z_eff = st.sidebar.slider(
        "Effective Redshift (z_eff)",
        min_value=0.0,
        max_value=2.0,
        value=0.8,
        step=0.1,
        help="Effective redshift for power spectrum"
    )
    
    k_min = st.sidebar.slider(
        "k_min [h/Mpc]",
        min_value=0.001,
        max_value=0.05,
        value=0.01,
        step=0.001,
        format="%.3f"
    )
    
    k_max = st.sidebar.slider(
        "k_max [h/Mpc]",
        min_value=0.1,
        max_value=1.0,
        value=0.3,
        step=0.01
    )
    
    include_systematics = st.sidebar.checkbox(
        "Include Systematic Errors",
        value=True,
        disabled=not HAS_SYSTEMATICS,
        help="Include systematic error budget in forecasts"
    )
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("### Cosmological Parameters")
    
    H0 = st.sidebar.slider("H₀ [km/s/Mpc]", 60.0, 75.0, float(model.params['H0']), 0.1)
    ombh2 = st.sidebar.slider("Ω_b h²", 0.020, 0.025, float(model.params['ombh2']), 0.0001, format="%.4f")
    omch2 = st.sidebar.slider("Ω_c h²", 0.10, 0.14, float(model.params['omch2']), 0.001)
    ns = st.sidebar.slider("n_s", 0.90, 1.00, float(model.params['ns']), 0.001)
    
    # Update model parameters if changed
    if (H0 != model.params['H0'] or ombh2 != model.params['ombh2'] or 
        omch2 != model.params['omch2'] or ns != model.params['ns']):
        custom_params = model.params.copy()
        custom_params.update({'H0': H0, 'ombh2': ombh2, 'omch2': omch2, 'ns': ns})
        model = PhiModulationModel(params=custom_params)
    
    # Main content tabs
    tab1, tab2, tab3, tab4 = st.tabs(["📊 Forecast Analysis", "🔬 Systematic Errors", "📈 Power Spectrum", "📋 Summary"])
    
    with tab1:
        st.header("DESI Forecast Analysis")
        
        # Run forecast
        try:
            if include_systematics and HAS_SYSTEMATICS:
                forecast = model.forecast_desi_sensitivity_with_systematics(
                    A_phi_true=A_phi,
                    k_min=k_min,
                    k_max=k_max,
                    n_k=100,
                    include_systematics=True
                )
                sigma_Aphi_stat = forecast['sigma_Aphi_stat']
                sigma_Aphi_sys = forecast['sigma_Aphi_sys']
                sigma_Aphi_total = forecast['sigma_Aphi_total']
                has_systematics = True
            else:
                forecast = model.forecast_desi_sensitivity(
                    A_phi_true=A_phi,
                    k_min=k_min,
                    k_max=k_max,
                    n_k=100
                )
                sigma_Aphi_stat = forecast['sigma_Aphi']
                sigma_Aphi_sys = 0.0
                sigma_Aphi_total = forecast['sigma_Aphi']
                has_systematics = False
            
            SNR = forecast['SNR']
            k = forecast['k']
            Pk_base = forecast['Pk_base']
            Pk_mod = forecast['Pk_mod']
            
            # Metrics
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric(
                    "A_φ (True)",
                    f"{A_phi:.4f}",
                    help="True modulation amplitude"
                )
            
            with col2:
                st.metric(
                    "σ_Aφ (Statistical)",
                    f"{sigma_Aphi_stat:.6f}",
                    help="Statistical uncertainty"
                )
            
            with col3:
                if has_systematics:
                    st.metric(
                        "σ_Aφ (Systematic)",
                        f"{sigma_Aphi_sys:.6f}",
                        help="Systematic uncertainty"
                    )
                else:
                    st.metric("σ_Aφ (Systematic)", "N/A", help="Systematics module not available")
            
            with col4:
                st.metric(
                    "Signal-to-Noise",
                    f"{SNR:.2f}σ",
                    delta=f"{SNR - 3:.2f}σ" if SNR > 3 else None,
                    help="Detection significance"
                )
            
            st.markdown("---")
            
            # Forecast visualization
            fig, axes = plt.subplots(2, 2, figsize=(14, 10))
            
            # Panel 1: Power spectrum ratio
            ax = axes[0, 0]
            ratio = (Pk_mod / Pk_base - 1) * 100
            ax.semilogx(k, ratio, 'b-', linewidth=2, label=f'A_φ = {A_phi:.3f}')
            ax.axhline(0, color='k', linestyle='--', alpha=0.3)
            ax.set_xlabel('k [h/Mpc]', fontsize=11)
            ax.set_ylabel('ΔP/P [%]', fontsize=11)
            ax.set_title('φ-modulation in Power Spectrum', fontsize=12)
            ax.grid(True, alpha=0.3)
            ax.legend()
            
            # Panel 2: Power spectrum comparison
            ax = axes[0, 1]
            ax.loglog(k, Pk_base, 'k-', linewidth=2, label='ΛCDM base', alpha=0.7)
            ax.loglog(k, Pk_mod, 'b--', linewidth=2, label='φ-modulated', alpha=0.7)
            ax.set_xlabel('k [h/Mpc]', fontsize=11)
            ax.set_ylabel('P(k) [(Mpc/h)³]', fontsize=11)
            ax.set_title('Power Spectrum Comparison', fontsize=12)
            ax.grid(True, alpha=0.3)
            ax.legend()
            
            # Panel 3: Error comparison
            ax = axes[1, 0]
            if has_systematics and 'sigma_P' in forecast:
                sigma_P_stat = forecast.get('sigma_P', Pk_base * 0.1)
                if 'systematic_budget' in forecast:
                    sys_budget = forecast['systematic_budget']
                    sigma_P_total = sys_budget['sigma_P_total']
                    ax.semilogx(k, sigma_P_stat * 100 / Pk_base, 'k-', linewidth=2, label='Statistical')
                    ax.semilogx(k, sys_budget['sigma_P_sys'] * 100 / Pk_base, 'r--', linewidth=2, label='Systematic')
                    ax.semilogx(k, sigma_P_total * 100 / Pk_base, 'b:', linewidth=2, label='Total')
                    ax.set_xlabel('k [h/Mpc]', fontsize=11)
                    ax.set_ylabel('Relative Error [%]', fontsize=11)
                    ax.set_title('Error Breakdown', fontsize=12)
                    ax.grid(True, alpha=0.3)
                    ax.legend()
            
            # Panel 4: SNR visualization
            ax = axes[1, 1]
            # Create SNR bar chart
            categories = ['Statistical', 'Systematic', 'Total']
            if has_systematics:
                values = [
                    A_phi / sigma_Aphi_stat,
                    A_phi / sigma_Aphi_sys if sigma_Aphi_sys > 0 else 0,
                    SNR
                ]
            else:
                values = [SNR, 0, SNR]
            
            colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
            bars = ax.bar(categories, values, color=colors, alpha=0.7)
            ax.axhline(3, color='r', linestyle='--', linewidth=2, label='3σ threshold')
            ax.set_ylabel('Signal-to-Noise Ratio', fontsize=11)
            ax.set_title('Detection Significance', fontsize=12)
            ax.grid(True, alpha=0.3, axis='y')
            ax.legend()
            
            # Add value labels on bars
            for bar, val in zip(bars, values):
                if val > 0:
                    height = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width()/2., height,
                           f'{val:.2f}σ',
                           ha='center', va='bottom', fontsize=10)
            
            plt.tight_layout()
            st.pyplot(fig)
            
        except Exception as e:
            st.error(f"Error running forecast: {str(e)}")
            st.exception(e)
    
    with tab2:
        st.header("Systematic Error Analysis")
        
        if not HAS_SYSTEMATICS:
            st.warning("⚠️ SystematicErrorBudget module not available. Please ensure src/systematics.py is in the path.")
        else:
            try:
                # Get power spectrum for systematics analysis
                k, z_arr, Pk_full = model.get_base_power_spectrum(
                    k_min=k_min*0.5, k_max=k_max*2, npoints=200, z=z_eff
                )
                Pk = Pk_full[0] if len(Pk_full.shape) > 1 else Pk_full
                
                # Interpolate to desired k range
                from scipy.interpolate import interp1d
                k_target = np.logspace(np.log10(k_min), np.log10(k_max), 100)
                Pk_interp = interp1d(k, Pk, kind='linear', bounds_error=False, fill_value='extrapolate')
                Pk_base = Pk_interp(k_target)
                
                # Statistical error (10% for demonstration)
                sigma_P_stat = Pk_base * 0.1
                
                # Create systematic error budget
                sys_budget = SystematicErrorBudget(z_eff=z_eff)
                
                # Compute systematic budget
                sys_result = sys_budget.compute_systematic_budget(
                    k_target, Pk_base, sigma_P_stat
                )
                
                # Display metrics
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    avg_photoz = np.mean(sys_result['sigma_P_photoz'] / sys_result['sigma_P_total']) * 100
                    st.metric("Photo-z Contribution", f"{avg_photoz:.1f}%")
                
                with col2:
                    avg_bias = np.mean(sys_result['sigma_P_bias'] / sys_result['sigma_P_total']) * 100
                    st.metric("Bias Contribution", f"{avg_bias:.1f}%")
                
                with col3:
                    avg_geometry = np.mean(sys_result['sigma_P_geometry'] / sys_result['sigma_P_total']) * 100
                    st.metric("Geometry Contribution", f"{avg_geometry:.1f}%")
                
                st.markdown("---")
                
                # Visualization
                fig, axes = plt.subplots(2, 2, figsize=(14, 10))
                
                # Panel 1: Systematic error components
                ax = axes[0, 0]
                ax.semilogx(k_target, sys_result['sigma_P_photoz'] * 100 / Pk_base, 
                           'r-', linewidth=2, label='Photo-z')
                ax.semilogx(k_target, sys_result['sigma_P_bias'] * 100 / Pk_base, 
                           'g--', linewidth=2, label='Bias')
                ax.semilogx(k_target, sys_result['sigma_P_geometry'] * 100 / Pk_base, 
                           'b:', linewidth=2, label='Geometry')
                ax.set_xlabel('k [h/Mpc]', fontsize=11)
                ax.set_ylabel('Relative Error [%]', fontsize=11)
                ax.set_title('Systematic Error Components', fontsize=12)
                ax.grid(True, alpha=0.3)
                ax.legend()
                
                # Panel 2: Total error breakdown
                ax = axes[0, 1]
                ax.semilogx(k_target, sigma_P_stat * 100 / Pk_base, 
                           'k-', linewidth=2, label='Statistical')
                ax.semilogx(k_target, sys_result['sigma_P_sys'] * 100 / Pk_base, 
                           'r--', linewidth=2, label='Systematic')
                ax.semilogx(k_target, sys_result['sigma_P_total'] * 100 / Pk_base, 
                           'b:', linewidth=2, label='Total')
                ax.set_xlabel('k [h/Mpc]', fontsize=11)
                ax.set_ylabel('Relative Error [%]', fontsize=11)
                ax.set_title('Error Breakdown', fontsize=12)
                ax.grid(True, alpha=0.3)
                ax.legend()
                
                # Panel 3: Systematic fraction
                ax = axes[1, 0]
                ax.semilogx(k_target, sys_result['fraction_sys'] * 100, 
                           'purple', linewidth=2)
                ax.axhline(50, color='k', linestyle='--', alpha=0.3, label='50%')
                ax.set_xlabel('k [h/Mpc]', fontsize=11)
                ax.set_ylabel('Systematic Fraction [%]', fontsize=11)
                ax.set_title('Fraction of Total Error from Systematics', fontsize=12)
                ax.set_ylim([0, 100])
                ax.grid(True, alpha=0.3)
                ax.legend()
                
                # Panel 4: Pie chart of error contributions
                ax = axes[1, 1]
                mean_photoz = np.mean(sys_result['sigma_P_photoz']**2)
                mean_bias = np.mean(sys_result['sigma_P_bias']**2)
                mean_geometry = np.mean(sys_result['sigma_P_geometry']**2)
                total_sys_sq = mean_photoz + mean_bias + mean_geometry
                
                if total_sys_sq > 0:
                    sizes = [mean_photoz/total_sys_sq*100, 
                            mean_bias/total_sys_sq*100, 
                            mean_geometry/total_sys_sq*100]
                    labels = ['Photo-z', 'Bias', 'Geometry']
                    colors = ['#ff7f0e', '#2ca02c', '#1f77b4']
                    ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
                    ax.set_title('Systematic Error Contributions\n(Mean Square)', fontsize=12)
                
                plt.tight_layout()
                st.pyplot(fig)
                
            except Exception as e:
                st.error(f"Error in systematic error analysis: {str(e)}")
                st.exception(e)
    
    with tab3:
        st.header("Power Spectrum Visualization")
        
        try:
            # Get power spectrum
            k_full, z_arr, Pk_full = model.get_base_power_spectrum(
                k_min=k_min*0.5, k_max=k_max*2, npoints=500, z=z_eff
            )
            Pk_base_full = Pk_full[0] if len(Pk_full.shape) > 1 else Pk_full
            
            # Apply modulation
            Pk_mod_full, mod_factor_full = model.apply_phi_modulation(
                k_full, Pk_base_full, A_phi=A_phi
            )
            
            # Filter to desired range
            mask = (k_full >= k_min) & (k_full <= k_max)
            k = k_full[mask]
            Pk_base = Pk_base_full[mask]
            Pk_mod = Pk_mod_full[mask]
            mod_factor = mod_factor_full[mask]
            
            # Create visualization
            fig, axes = plt.subplots(3, 1, figsize=(12, 12))
            
            # Panel 1: Power spectrum
            ax = axes[0]
            ax.loglog(k, Pk_base, 'k-', linewidth=2, label='ΛCDM base', alpha=0.7)
            ax.loglog(k, Pk_mod, 'b--', linewidth=2, label='φ-modulated', alpha=0.7)
            ax.set_xlabel('k [h/Mpc]', fontsize=11)
            ax.set_ylabel('P(k) [(Mpc/h)³]', fontsize=11)
            ax.set_title(f'Power Spectrum (z = {z_eff})', fontsize=12)
            ax.grid(True, alpha=0.3)
            ax.legend()
            
            # Panel 2: Modulation factor
            ax = axes[1]
            ax.semilogx(k, mod_factor, 'g-', linewidth=2)
            ax.axhline(1, color='k', linestyle='--', alpha=0.3)
            ax.set_xlabel('k [h/Mpc]', fontsize=11)
            ax.set_ylabel('Modulation Factor', fontsize=11)
            ax.set_title('φ-Modulation Factor', fontsize=12)
            ax.grid(True, alpha=0.3)
            
            # Panel 3: Relative difference
            ax = axes[2]
            ratio = (Pk_mod / Pk_base - 1) * 100
            ax.semilogx(k, ratio, 'r-', linewidth=2)
            ax.axhline(0, color='k', linestyle='--', alpha=0.3)
            ax.fill_between(k, -A_phi*100*2, A_phi*100*2, alpha=0.1, color='gray')
            ax.set_xlabel('k [h/Mpc]', fontsize=11)
            ax.set_ylabel('(P_mod / P_base - 1) [%]', fontsize=11)
            ax.set_title('Relative Difference', fontsize=12)
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            st.pyplot(fig)
            
        except Exception as e:
            st.error(f"Error generating power spectrum: {str(e)}")
            st.exception(e)
    
    with tab4:
        st.header("Analysis Summary")
        
        try:
            # Run forecast for summary
            if include_systematics and HAS_SYSTEMATICS:
                forecast = model.forecast_desi_sensitivity_with_systematics(
                    A_phi_true=A_phi,
                    k_min=k_min,
                    k_max=k_max,
                    include_systematics=True
                )
            else:
                forecast = model.forecast_desi_sensitivity(
                    A_phi_true=A_phi,
                    k_min=k_min,
                    k_max=k_max
                )
            
            st.markdown("### Forecast Parameters")
            col1, col2 = st.columns(2)
            
            with col1:
                st.write(f"**Modulation Amplitude (A_φ):** {A_phi:.4f}")
                st.write(f"**Effective Redshift:** {z_eff:.2f}")
                st.write(f"**k range:** {k_min:.3f} - {k_max:.2f} h/Mpc")
            
            with col2:
                st.write(f"**H₀:** {model.params['H0']:.2f} km/s/Mpc")
                st.write(f"**Ω_b h²:** {model.params['ombh2']:.4f}")
                st.write(f"**Ω_c h²:** {model.params['omch2']:.3f}")
                st.write(f"**n_s:** {model.params['ns']:.4f}")
            
            st.markdown("---")
            st.markdown("### Forecast Results")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.metric("Statistical σ_Aφ", f"{forecast.get('sigma_Aphi_stat', forecast['sigma_Aphi']):.6f}")
            
            with col2:
                if include_systematics and HAS_SYSTEMATICS:
                    st.metric("Systematic σ_Aφ", f"{forecast.get('sigma_Aphi_sys', 0):.6f}")
                else:
                    st.metric("Systematic σ_Aφ", "N/A")
            
            with col3:
                st.metric("Total σ_Aφ", f"{forecast.get('sigma_Aphi_total', forecast['sigma_Aphi']):.6f}")
            
            st.markdown("---")
            st.markdown("### Detection Significance")
            
            SNR = forecast['SNR']
            
            if SNR >= 5:
                st.success(f"✅ **Very Strong Detection:** {SNR:.2f}σ (≥5σ threshold)")
            elif SNR >= 3:
                st.success(f"✅ **Strong Detection:** {SNR:.2f}σ (≥3σ threshold)")
            elif SNR >= 2:
                st.warning(f"⚠️ **Marginal Detection:** {SNR:.2f}σ (≥2σ but <3σ)")
            else:
                st.info(f"ℹ️ **Below Detection Threshold:** {SNR:.2f}σ (<2σ)")
            
            st.markdown("---")
            st.markdown("### Interpretation")
            
            interpretation = f"""
            For A_φ = {A_phi:.4f}:
            - DESI Year 5 can detect this modulation with **{SNR:.2f}σ** significance
            - The statistical uncertainty is σ_Aφ = {forecast.get('sigma_Aphi_stat', forecast['sigma_Aphi']):.6f}
            """
            
            if include_systematics and HAS_SYSTEMATICS:
                interpretation += f"""
            - Including systematic errors, the total uncertainty is σ_Aφ = {forecast.get('sigma_Aphi_total', forecast['sigma_Aphi']):.6f}
            - Systematic errors increase the uncertainty by {((forecast.get('sigma_Aphi_total', forecast['sigma_Aphi']) / forecast.get('sigma_Aphi_stat', forecast['sigma_Aphi']) - 1) * 100):.1f}%
                """
            
            st.info(interpretation)
            
        except Exception as e:
            st.error(f"Error generating summary: {str(e)}")
            st.exception(e)

else:
    st.error("""
    # ⚠️ Module Not Found
    
    The `PhiModulationModel` class could not be imported. 
    
    Please ensure:
    1. The `src/phi_modulation.py` file exists
    2. All dependencies are installed (especially CAMB)
    3. The Python path is configured correctly
    
    To install dependencies:
    ```bash
    pip install -r requirements.txt
    ```
    """)

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666;'>
    <p>φ-Modulation Analysis Dashboard | FaCC Framework</p>
    <p>For more information, see the <a href='https://github.com/imediacorp/FaCC'>repository</a></p>
</div>
""", unsafe_allow_html=True)

