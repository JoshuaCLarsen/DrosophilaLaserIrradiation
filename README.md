# Drosophila Laser Irradiation

## Overview

This repository is being curated to accompany the study *Laser Induced Mortality on Vinegar Flies in Free Fall*. It will contain the code, data, and supporting files used to characterize the laser and model its effects on *Drosophila hydei*.

## Repository contents

### `LaserCharacterization/`

Files related to measuring and modeling the laser beam profile:

- `rawLaserIntensityPattern.png` - captured image of the blue laser beam on a diffuse target. This is the image input to the profile-fitting script.
- `laserIntensityFitter.m` - MATLAB function that reads the beam image, converts image pixels to millimeter coordinates, fits a 2D super-Gaussian with fixed exponents 2 and 6, and writes the fit parameters and diagnostic plots.

The script generates `LaserFittingResults.csv`, `laserImage.pdf`, `LaserFit.pdf`, and `LaserResiduals.pdf` in the current working directory. The fitting step requires MATLAB's Curve Fitting Toolbox; image reading and graphics also use MATLAB functions.

## Requirements

- MATLAB
- MATLAB Curve Fitting Toolbox (for `fittype` and `fit`)

## Usage

From MATLAB, set the current folder to `LaserCharacterization` and call the function. Its defaults match the manuscript's beam-image setup: 1.0 m camera-to-target distance and 35 mm focal length. Both input values are in meters.

```matlab
laserFit = laserIntensityFitter();
```

To override the defaults, provide both inputs, for example `laserFit = laserIntensityFitter(1.0, 35e-3);`.

## Data and outputs

The beam image is provided as the input measurement. The function returns a structure containing the fitted peak intensity, widths, center, and fixed super-Gaussian exponents; it also writes the fit parameters and diagnostic figures listed above. Run MATLAB with `LaserCharacterization/` as the current folder because image and output filenames are relative to the current folder.

## Reproducibility notes

The function accepts camera-to-target distance and lens focal length as inputs rather than hard-coding the focal length. Use the values recorded for the acquisition to reproduce the millimeter calibration; the arguments must be in meters.

The manuscript reports fitted 1/e half-widths of 36.8 mm and 40.6 mm at the 1 m plane, with super-Gaussian exponents 2 and 6. The script uses the same fixed exponents and processes the blue image channel. Its reported widths depend on the focal length and image distance supplied for pixel scaling, so those acquisition values should accompany any published fit results. The raw beam image and fitting script are appropriate to retain in the repository.

## Citation

Citation details will be added when the associated paper is finalized.

## License

License information will be added.
