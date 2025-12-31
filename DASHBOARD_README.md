# Streamlit Dashboard

Interactive web dashboard for exploring φ-modulation analysis and DESI forecasts.

## Quick Start

### Install Dependencies

```bash
pip install streamlit
# or
pip install -r requirements.txt  # includes streamlit
```

### Run the Dashboard

```bash
streamlit run streamlit_app.py
```

The dashboard will open in your default web browser at `http://localhost:8501`

## Features

### 📊 Forecast Analysis Tab
- Interactive parameter sliders (A_φ, redshift, k-range)
- Real-time power spectrum visualization
- Forecast sensitivity calculations
- Signal-to-noise ratio (SNR) visualization
- 4-panel forecast figure

### 🔬 Systematic Errors Tab
- Photo-z error analysis
- Galaxy bias uncertainty propagation
- Survey geometry effects
- Systematic error breakdown visualization
- Error contribution pie chart

### 📈 Power Spectrum Tab
- Power spectrum comparison (ΛCDM vs. φ-modulated)
- Modulation factor visualization
- Relative difference plots
- Multi-panel visualization

### 📋 Summary Tab
- Complete forecast summary
- Parameter overview
- Detection significance interpretation
- Key metrics and statistics

## Dashboard Controls

### Sidebar Parameters

- **A_φ (Modulation Amplitude):** 0.001 - 0.05
- **Effective Redshift (z_eff):** 0.0 - 2.0
- **k_min, k_max:** Wavenumber range [h/Mpc]
- **Include Systematic Errors:** Toggle systematic error analysis
- **Cosmological Parameters:** H₀, Ω_b h², Ω_c h², n_s

### Interactive Features

- **Real-time Updates:** All visualizations update automatically when parameters change
- **Multiple Views:** Switch between tabs for different aspects of the analysis
- **Export Ready:** Figures are publication-quality and can be saved

## Requirements

- Python 3.9+
- Streamlit >= 1.28.0
- All dependencies from `requirements.txt`
- CAMB (for power spectrum generation)
- Optional: Systematics module (for systematic error analysis)

## Troubleshooting

### Dashboard won't start

1. Check that Streamlit is installed:
   ```bash
   pip install streamlit
   ```

2. Verify modules are accessible:
   ```bash
   python -c "import sys; sys.path.insert(0, 'src'); from phi_modulation import PhiModulationModel; print('OK')"
   ```

3. Check for CAMB:
   ```bash
   python -c "import camb; print('CAMB version:', camb.__version__)"
   ```

### Module not found errors

- Ensure you're running from the repository root directory
- The dashboard automatically adds `src/` to the Python path
- If issues persist, try: `export PYTHONPATH="${PYTHONPATH}:$(pwd)/src"`

### Slow performance

- Reduce `n_k` parameter in forecast calls (default is 100)
- Use smaller k ranges for faster calculations
- Close other resource-intensive applications

## Deployment

### Local Network Access

To access the dashboard from other devices on your network:

```bash
streamlit run streamlit_app.py --server.address 0.0.0.0
```

Then access via: `http://YOUR_IP_ADDRESS:8501`

### Cloud Deployment

#### Streamlit Cloud (Recommended)

1. Push your repository to GitHub
2. Go to https://streamlit.io/cloud
3. Connect your GitHub account
4. Select your repository
5. Deploy!

The dashboard will be automatically deployed and accessible via a public URL.

#### Heroku

See Streamlit's deployment guide: https://docs.streamlit.io/streamlit-community-cloud/deploy-your-app

#### Docker

Create a `Dockerfile`:

```dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

EXPOSE 8501

CMD ["streamlit", "run", "streamlit_app.py", "--server.port=8501", "--server.address=0.0.0.0"]
```

Build and run:
```bash
docker build -t facc-dashboard .
docker run -p 8501:8501 facc-dashboard
```

## Usage Tips

1. **Start with defaults:** Use the default parameters to get familiar with the dashboard
2. **Explore systematically:** Change one parameter at a time to understand its effect
3. **Check all tabs:** Each tab provides different insights into the analysis
4. **Use tooltips:** Hover over parameter controls for helpful descriptions
5. **Compare scenarios:** Open multiple browser windows to compare different parameter sets

## Screenshots

The dashboard includes:
- Interactive parameter controls
- Real-time visualization updates
- Multi-panel forecast figures
- Systematic error breakdowns
- Summary statistics and interpretations

## Support

For issues or questions:
- Check the main README.md
- Open a GitHub issue
- Review the documentation in `docs/`

## License

Same license as the main FaCC project (see LICENSE file)

