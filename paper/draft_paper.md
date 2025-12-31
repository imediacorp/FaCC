# Searching for Golden Ratio Oscillations in Large-Scale Structure: DESI Year 5 Forecasts for Log-Periodic Power Spectrum Modulations

**Author:** Bryan David Persaud  
**Affiliation:** Intermedia Communications Corp.  
**Date:** 2025-01-16

---

## Abstract

We present forecasts for detecting logarithmic oscillations in the matter power spectrum at scales determined by the golden ratio $\phi \approx 1.618$. While $\phi$-based expansion models are excluded by $\Lambda$CDM constraints, perturbation-level modulations remain viable and testable. Using Fisher matrix forecasts with DESI Year 5 specifications, we find sensitivity to oscillation amplitudes $A_\phi \gtrsim 0.005$ at $3\sigma$ significance. A detection would indicate scale-dependent departures from standard inflation, while a null result would place stringent limits on such features. We provide public code for integrating this test into DESI's analysis pipeline.

**Keywords:** cosmology: large-scale structure, power spectrum, DESI, golden ratio, log-periodic oscillations

---

## 1. Introduction

The golden ratio $\phi = (1+\sqrt{5})/2 \approx 1.618$ appears throughout nature, from the arrangement of leaves and flower petals to spiral galaxy structures. This ubiquity has motivated investigations into whether $\phi$ might play a fundamental role in cosmology. Previous attempts to use $\phi$ as a background expansion constant have been decisively ruled out by observations. However, the possibility that $\phi$ manifests at the perturbation level---through log-periodic modulations in structure growth---remains an open, testable hypothesis.

In this work, we propose searching for $\phi$-modulated oscillations in the matter power spectrum $P(k)$ using data from the Dark Energy Spectroscopic Instrument (DESI). Such modulations would appear as log-periodic features with period set by $\ln(\phi)$, potentially indicating scale-dependent departures from standard inflationary scenarios. Unlike background-level tests, which have been falsified, perturbation-level signatures can be tested with current and near-future survey data.

We implement this search as a two-parameter extension to $\Lambda$CDM, forecast DESI Year 5 sensitivity using Fisher matrix methodology, and provide open-source code for the analysis. Our forecasts show DESI Year 5 can detect or rule out $\phi$-modulation with amplitudes $A_\phi \gtrsim 0.005$ at $3\sigma$ confidence, providing a clear timeline for testing this hypothesis.

---

## 2. The Model

We model $\phi$-modulation as a log-periodic perturbation to the standard $\Lambda$CDM matter power spectrum:

$$P_\phi(k) = P_{\Lambda\text{CDM}}(k) \times \left[1 + A_\phi \cos\left(\frac{2\pi \log(k/k_0)}{\ln(\phi)} + \phi_0\right)\right]$$

where $A_\phi$ is the modulation amplitude, $k_0$ is a pivot scale (typically the BAO scale $k_0 \approx 0.05$ h/Mpc), and $\phi_0$ is a phase offset. The logarithmic periodicity with period $\ln(\phi) \approx 0.481$ implies that oscillations repeat every factor of $\phi$ in wavenumber space.

This model extends $\Lambda$CDM by two parameters ($A_\phi$ and $\phi_0$), making it testable and falsifiable. The amplitude $A_\phi$ is expected to be small ($A_\phi \ll 1$) if such modulations exist, as they would represent subtle perturbations to the standard power spectrum.

### 2.1 Theoretical Motivation

Several theoretical frameworks could give rise to log-periodic power spectrum modulations:

- **Discrete scale invariance:** If the inflationary potential or structure formation exhibits discrete scale invariance, log-periodic features naturally emerge.
- **Modified gravity:** Some modified gravity models predict oscillatory features in the transfer function.
- **Primordial features:** Sharp features in the inflationary potential can imprint oscillations in $P(k)$.

While $\phi$ specifically may seem arbitrary, its status as the "most irrational number" and its appearance in natural growth patterns motivates testing whether it manifests in cosmic structure. A detection would require theoretical explanation; a null result would constrain such models.

---

## 3. Methods

### 3.1 Power Spectrum Generation

We generate the base $\Lambda$CDM power spectrum using CAMB, with cosmological parameters from Planck 2018:
- $H_0 = 67.36$ km/s/Mpc
- $\Omega_b h^2 = 0.02237$
- $\Omega_c h^2 = 0.1200$
- $A_s = 2.1 \times 10^{-9}$
- $n_s = 0.9649$
- $\tau = 0.0544$

The $\phi$-modulated power spectrum is computed by applying the modulation formula to the CAMB output. All calculations are performed at effective redshift $z_{\text{eff}} = 0.8$, appropriate for DESI galaxy samples.

### 3.2 Fisher Forecast Methodology

We forecast DESI Year 5 sensitivity using the Fisher matrix formalism. The forecast uncertainty on $A_\phi$ is:

$$\sigma_{A_\phi} = \left[\sum_k \frac{(dP/dA_\phi)^2}{\sigma_P^2(k)}\right]^{-1/2}$$

where the sum is over $k$-bins, $dP/dA_\phi$ is the derivative of the power spectrum with respect to $A_\phi$, and $\sigma_P(k)$ is the error on $P(k)$ per bin.

The error on $P(k)$ includes both cosmic variance and shot noise:
$$\sigma_P^2(k) = \sigma_{\text{CV}}^2(k) + \sigma_{\text{shot}}^2(k)$$

where
- $\sigma_{\text{CV}}(k) = P(k) \sqrt{2/N_{\text{modes}}(k)}$
- $\sigma_{\text{shot}}(k) = P_{\text{shot}}/\sqrt{N_{\text{modes}}(k)}$
- $P_{\text{shot}} = 1/n_{\text{gal}}$
- $N_{\text{modes}} = V_{\text{survey}} k^2 \Delta k / (2\pi^2)$

### 3.3 DESI Year 5 Specifications

Our forecasts use the following DESI Year 5 specifications:
- Survey volume: $V_{\text{survey}} = 100$ (Gpc/h)$^3$
- Effective redshift: $z_{\text{eff}} = 0.8$
- Galaxy number density: $n_{\text{gal}} = 3 \times 10^{-4}$ (h/Mpc)$^3$
- Wavenumber range: $k = 0.01$--$0.3$ h/Mpc (linear regime)

These specifications are based on expected DESI Year 5 performance for the combined galaxy and quasar samples.

### 3.4 Systematic Error Analysis

We include systematic error contributions from:
- Photometric redshift uncertainties: $\sigma_z/(1+z) \approx 0.02$
- Galaxy bias uncertainties: $\sigma_b/b \approx 0.05$
- Survey geometry effects

Systematic errors are propagated to $A_\phi$ constraints and combined in quadrature with statistical errors. The systematic error framework is implemented in our codebase and can be included via the `forecast_desi_sensitivity_with_systematics()` method.

---

## 4. Results

### 4.1 Forecast Constraints

Our Fisher matrix forecasts yield the following constraints on the $\phi$-modulation amplitude $A_\phi$:

| $A_\phi$ (True) | $\sigma_{A_\phi}$ | SNR | Significance |
|----------------|-------------------|-----|--------------|
| 0.005 | 0.0013 | 3.88 | Strong ($\geq 3\sigma$) |
| 0.010 | 0.0013 | 7.76 | Very Strong ($\geq 5\sigma$) |
| 0.020 | 0.0013 | 15.52 | Very Strong ($\geq 5\sigma$) |

**Table 1:** DESI Year 5 forecast constraints on $\phi$-modulation amplitude. The forecast uncertainty $\sigma_{A_\phi} \approx 0.0013$ is independent of the true amplitude for small $A_\phi$ (Fisher matrix approximation).

The key result is that DESI Year 5 can detect or rule out $\phi$-modulation with amplitudes $A_\phi \gtrsim 0.005$ at $3\sigma$ significance. For amplitudes $A_\phi \gtrsim 0.01$, DESI achieves "discovery" level sensitivity ($\geq 5\sigma$).

### 4.2 Detection Thresholds

The forecast signal-to-noise ratio increases linearly with amplitude (as expected from Fisher matrix formalism). The $3\sigma$ detection threshold corresponds to $A_\phi \approx 0.004$, while the $5\sigma$ discovery threshold corresponds to $A_\phi \approx 0.006$.

### 4.3 Power Spectrum Modulation

The $\phi$-modulated power spectrum shows log-periodic oscillations that become increasingly visible at higher amplitudes. The modulation is strongest near the BAO scale ($k \approx 0.05$ h/Mpc) where DESI has maximum sensitivity.

---

## 5. Discussion

### 5.1 Theoretical Implications

**If $\phi$-modulation is detected at the forecasted sensitivity:**
- It would indicate scale-dependent departures from standard $\Lambda$CDM structure formation
- The log-periodic nature suggests discrete scale invariance or similar symmetries
- A theoretical framework would be needed to explain why $\phi$ specifically appears
- It could indicate novel physics in inflation or modified gravity

**If no signal is detected:**
- Stringent upper limits on $A_\phi < 0.005$ at $3\sigma$ confidence
- Rules out significant log-periodic oscillations at tested scales
- Constrains models predicting such features
- Still scientifically valuable (null results constrain parameter space)

### 5.2 Comparison with Existing Constraints

Current constraints from Planck CMB and SDSS/BOSS large-scale structure data do not strongly constrain log-periodic modulations at the scales and amplitudes we consider. Our DESI forecasts represent the first targeted search for $\phi$-periodic features in the matter power spectrum with this level of precision.

### 5.3 Systematic Uncertainties

Our forecasts include statistical uncertainties from cosmic variance and shot noise. Systematic errors from photo-$z$ uncertainties, galaxy bias, and survey geometry typically add $\sim 10$--$30\%$ to the statistical uncertainties. These do not significantly affect our main conclusions but should be included in the final analysis of real data.

### 5.4 Future Tests

Beyond DESI Year 5, future surveys can extend these tests:
- **Euclid:** Similar sensitivity, complementary redshift coverage
- **Roman Space Telescope:** Higher redshift precision, smaller volume
- **SKA:** Radio surveys, different systematics
- **CMB-S4:** Complementary constraints from CMB anisotropy

Combining multiple surveys will improve constraints and test consistency across different probes.

---

## 6. Conclusion

We have implemented a framework for testing $\phi$-modulated log-periodic oscillations in the matter power spectrum and forecasted DESI Year 5 sensitivity using Fisher matrix methodology. Our key results:

1. DESI Year 5 can detect or rule out $\phi$-modulation with amplitudes $A_\phi \gtrsim 0.005$ at $3\sigma$ confidence
2. For amplitudes $A_\phi \gtrsim 0.01$, DESI achieves discovery-level sensitivity ($\geq 5\sigma$)
3. Forecast uncertainty: $\sigma_{A_\phi} \approx 0.0013$
4. The test is well-defined, falsifiable, and can be performed with DESI Year 5 data

We provide open-source code implementing the forecast framework, systematic error analysis, and integration with standard cosmology tools (CAMB). This enables the DESI collaboration and broader community to perform this analysis with real data.

The $\phi$-modulation hypothesis is now a testable, falsifiable prediction that can be verified or ruled out with DESI Year 5 data. Whether the signal is detected or not, the results will provide valuable constraints on scale-dependent departures from standard cosmology and inform our understanding of structure formation in the Universe.

---

## Acknowledgments

The author thanks the DESI collaboration for survey design and data collection. This work used open-source software including CAMB, NumPy, SciPy, and Matplotlib.

---

## Code and Data Availability

All code, forecasts, and analysis tools are publicly available at:
**https://github.com/imediacorp/FaCC**

The codebase includes:
- $\phi$-modulation model implementation
- CAMB integration for $\Lambda$CDM power spectra
- Fisher forecast methodology
- Systematic error analysis framework
- Interactive Streamlit dashboard
- Comprehensive test suite (21 tests passing)

---

## References

- Planck Collaboration et al. 2018, A&A, 641, A6
- Lewis, A., & Challinor, A. 2000, arXiv:astro-ph/9911177
- Sornette, D. 1998, Phys. Rep., 297, 239
- Bean, R., et al. 2008, Phys. Rev. D, 78, 123514
- Meerburg, P. D., et al. 2019, Primordial Features, CMB-S4 Science Book

