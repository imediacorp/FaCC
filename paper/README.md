# Paper Draft

This directory contains the draft paper for the φ-modulation forecast analysis.

## Files

- `draft_paper.tex` - LaTeX version of the paper
- `draft_paper.md` - Markdown version (easier to read/edit)

## Paper Status

**Current Status:** Draft v1.0  
**Date:** 2025-01-16  
**Ready for:** Internal review, collaboration feedback

## Abstract

We present forecasts for detecting logarithmic oscillations in the matter power spectrum at scales determined by the golden ratio $\phi \approx 1.618$. Using Fisher matrix forecasts with DESI Year 5 specifications, we find sensitivity to oscillation amplitudes $A_\phi \gtrsim 0.005$ at $3\sigma$ significance.

## Key Results

- **Forecast uncertainty:** $\sigma_{A_\phi} \approx 0.0013$
- **Detection threshold:** $A_\phi \gtrsim 0.005$ at $3\sigma$
- **Discovery threshold:** $A_\phi \gtrsim 0.01$ at $5\sigma$
- **Methodology:** Fisher matrix forecasts with CAMB integration
- **Code:** Publicly available at https://github.com/imediacorp/FaCC

## Paper Structure

1. **Introduction** - Motivation, background, shift from falsified background model to testable perturbations
2. **The Model** - φ-modulation formalism, theoretical motivation
3. **Methods** - CAMB integration, Fisher forecasts, DESI specifications, systematic errors
4. **Results** - Forecast constraints, detection thresholds, power spectrum visualizations
5. **Discussion** - Theoretical implications, comparison with existing constraints, future tests
6. **Conclusion** - Summary, timeline for testing with real data

## Next Steps

1. **Review draft** - Internal review and refinement
2. **Generate figures** - Create publication-quality figures from notebook outputs
3. **Add references** - Complete bibliography with relevant literature
4. **Collaboration feedback** - Share with DESI collaboration members
5. **Submission** - Prepare for arXiv submission (forecast paper)
6. **Journal submission** - Target: JCAP or PRD after DESI data analysis

## Figure Requirements

The paper references several figures that should be generated:

1. **Figure 1:** Power spectrum modulation $P_\phi(k)/P_{\Lambda\text{CDM}}(k) - 1$ vs $k$ for different $A_\phi$
2. **Figure 2:** SNR vs amplitude plot (from notebook Panel 3)
3. **Figure 3:** BAO signature comparison (from notebook Panel 2)
4. **Figure 4:** Forecast summary visualization (optional)

These can be generated from the notebook `notebooks/01_desi_forecasts.ipynb`.

## Notes

- The paper emphasizes the **testable, falsifiable** nature of the hypothesis
- Honest about background model falsification
- Clear forecast methodology
- Ready for DESI data analysis
- Code availability statement included
