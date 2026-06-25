# Sequencing-coverage-calculator


Optimal coverage calculator for somatic variant detection.

This folder contains two MATLAB scripts for modeling somatic variant detection as a function of variant allele frequency (VAF), sequencing depth, sequencing error, and the minimum number of supporting reads required to call a variant.

These scripts were used for the analyses and visualizations associated with Figure 3 and Supplementary Figures S10-S11.

## Scripts

### `modelFitting.m`

Fits and compares several analytical detection-rate models against observed recall-rate data from HapMap mixture truth sets.

The script includes hard-coded recall-rate vectors for four sequencing centers:

- `TRUE_WASHU`
- `TRUE_BCM`
- `TRUE_NYGC`
- `TRUE_UW`

Each vector contains recall rates across six VAF levels: 10%, 5%, 2%, 1%, 0.5%, and 0.25%.

The script compares six model configurations:

- Binomial
- Binomial with sequencing error
- Binomial with Poisson-distributed sequencing error

It also computes SSE/RMSE summary statistics and generates model-comparison plots, including the RMSE comparison used for Figure 3 and related supplemental figures. The current version is configured to fit the WASHU dataset by default.

### `VAF_heatmap.m`

Generates a heatmap of expected variant detection rate across VAF levels and sequencing coverage values.

This script prompts the user for:

- `k0`: minimum number of supporting reads required to confirm variant existence
- `VAF`: variant allele frequency, entered as a fraction such as `0.01` for 1%
- `n`: total read depth
- `pE`: sequencing error rate, such as `0.001` for 0.1%

If the user presses Enter without providing a value, the script uses default values:

- `k0 = 1`
- `VAF = 0.01`
- `n = 500`
- `pE = 0.001`

The script reports the expected detection rate for the user-provided VAF/depth/error setting and creates a heatmap across preset VAF and coverage ranges.

## Requirements

- MATLAB
- MATLAB functions used by the scripts: `nchoosek`, `beta`, `heatmap`, `bar`, `plot`, `input`, and standard plotting functions
- `customcolormap_preset`, used to apply the red-white-blue heatmap colormap

If `customcolormap_preset` is not already installed or on the MATLAB path, either add it to the path or replace this line:

```matlab
J = customcolormap_preset('red-white-blue');
```

with any MATLAB-supported colormap, for example:

```matlab
J = parula;
```

## How to Run

Open MATLAB, navigate to the folder containing the scripts, and run:

```matlab
modelFitting
```

or:

```matlab
VAF_heatmap
```

In `modelFitting.m`, the default minimum supporting-read threshold is:

```matlab
k0 = 1;
```

The script is currently configured to model WASHU coverage only:

```matlab
n_total = [n_WASHU];
```

To compare all four sequencing centers, use:

```matlab
n_total = [n_WASHU, n_BCM, n_NYGC, n_UW];
```

Adjust `k0` or `n_total` if a different calling threshold or sequencing-center set is required.

## Model Notes

The scripts estimate detection probability by calculating the probability of observing at least `k0` supporting reads for a variant at a given VAF and read depth.

Sequencing error is modeled using:

```matlab
VAF * (1 - pE) + (1 - VAF) * pE / 3
```

where `pE` is the sequencing error rate and the `/3` term assumes errors are distributed among the three non-reference bases.

## Citation

If you use these scripts or results generated from them, please cite:

**A Pangenomic Method for Establishing a Somatic Variant Detection Resource in HapMap Mixtures**
