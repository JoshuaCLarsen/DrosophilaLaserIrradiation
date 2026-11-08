# Drosophila Laser Irradiation

## Overview

This repository is being curated to accompany the study *Laser Induced Mortality on Vinegar Flies in Free Fall*. It will contain the code, data, and supporting files used to characterize the laser and model its effects on *Drosophila hydei*.

## Repository contents

### `LaserCharacterization/`

Files related to measuring and modeling the laser beam profile:

- `rawLaserIntensityPattern.png` - captured image of the blue laser beam on a diffuse target. This is the image input to the profile-fitting script.
- `LaserDiodeCharacterization.txt` - measured laser output power (W) versus laser setting / PWM duty cycle (%), used for the manuscript's laser power characterization.
- `laserIntensityFitter.m` - MATLAB function that reads the beam image, converts image pixels to millimeter coordinates, fits a 2D super-Gaussian with fixed exponents 2 and 6, and writes the fit parameters and diagnostic plots.

The script generates `LaserFittingResults.csv`, `laserImage.pdf`, `LaserFit.pdf`, and `LaserResiduals.pdf` in the current working directory. The fitting step requires MATLAB's Curve Fitting Toolbox; image reading and graphics also use MATLAB functions.

### `MonteCarloSimulations/`

Files for the free-fall dose and mortality model described in the manuscript:

- `freefall_dose_model_monte6.py` - runs the manuscript's fixed-position model for uniform and wall-band fly entry positions, fits the dose-response parameters, and writes lethality curves, fit summaries, uncertainty data, and simulated fly energies.
- `freefall_dose_model_monte.py` - shared beam, trajectory, dose, and fitting implementation used by the run and plotting scripts; it also contains an earlier moving-fly trajectory model.
- `make_dose_figure_monte6.py` - builds the dose-distribution and example-path figure from the Monte Carlo outputs.
- `beam_geometry.py` - prints beam divergence, beam widths, coverage estimates, and transit geometry used in the manuscript.
- `flyDeathRates.txt` - measured laser power (W) and mortality fraction data for the ten experimental power settings, including the laser-off control. This is the fit input for the Monte Carlo analysis.
- `monte6Energy.npz` - included NumPy archive of per-fly simulated energies, example trajectories, and dose-model fit values used by the figure script.

The Monte Carlo model integrates the fitted super-Gaussian beam irradiance over a disc representing each fly as it traverses the 50 mm interaction region. The manuscript uses 20,000 simulated flies, 150 depth planes, and 61 points per fly disc, then fits a logistic fluence-response model to measured mortality.

## Requirements

- MATLAB
- MATLAB Curve Fitting Toolbox (for `fittype` and `fit`)
- Python 3 with NumPy and SciPy for the Monte Carlo scripts
- Matplotlib for `make_dose_figure_monte6.py`

## Usage

From MATLAB, set the current folder to `LaserCharacterization` and call the function. Its defaults match the manuscript's beam-image setup: 1.0 m camera-to-target distance and 35 mm focal length. Both input values are in meters.

```matlab
laserFit = laserIntensityFitter();
```

To override the defaults, provide both inputs, for example `laserFit = laserIntensityFitter(1.0, 35e-3);`.

### Monte Carlo simulation

From the repository root, run the manuscript model with:

```bash
python MonteCarloSimulations/freefall_dose_model_monte6.py
```

Then generate the dose figure with:

```bash
python MonteCarloSimulations/make_dose_figure_monte6.py
```

The scripts read `MonteCarloSimulations/flyDeathRates.txt` and write model summaries, curves, error bars, and figures into `MonteCarloSimulations/`. They do not require a separate `Figs/` directory. `monte6Energy.npz` is a precomputed simulation output, not the experimental input.

## Data and outputs

The beam image is provided as the input measurement. The function returns a structure containing the fitted peak intensity, widths, center, and fixed super-Gaussian exponents; it also writes the fit parameters and diagnostic figures listed above. Run MATLAB with `LaserCharacterization/` as the current folder because image and output filenames are relative to the current folder.

## Reproducibility notes

The function accepts camera-to-target distance and lens focal length as inputs rather than hard-coding the focal length. Use the values recorded for the acquisition to reproduce the millimeter calibration; the arguments must be in meters.

The manuscript reports fitted 1/e half-widths of 36.8 mm and 40.6 mm at the 1 m plane, with super-Gaussian exponents 2 and 6. The script uses the same fixed exponents and processes the blue image channel. Its reported widths depend on the focal length and image distance supplied for pixel scaling, so those acquisition values should accompany any published fit results. The raw beam image and fitting script are appropriate to retain in the repository.

The Monte Carlo scripts encode the model parameters in `freefall_dose_model_monte.py` and reproduce the uniform and wall-band fits discussed in the paper. The included mortality data supplies the experimental fractions used to fit the dose-response parameters and regenerate the reported lethality curves.

## Citation

Citation details will be added when the associated paper is finalized.

## License

License information will be added.
