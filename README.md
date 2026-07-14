# Sequencing Coverage Calculator

An optimal coverage calculator for estimating somatic variant detection probability from variant allele frequency (VAF), sequencing depth, sequencing error rate, and the minimum number of supporting reads required to call a variant.

This repository contains MATLAB modeling scripts and a browser-based version of the VAF heatmap calculator. The browser version was added so users can run the calculator without installing MATLAB.

The original MATLAB scripts were used for analyses and visualizations associated with Figure 3 and Supplementary Figures S10-S11.

## Browser Calculator

Open the self-contained HTML calculator:

```text
Sequencing_coverage_calculator.html
```

[Launch the calculator]([Sequencing_coverage_calculator.html](https://nahyunkong.github.io/Sequencing-coverage-calculator/Sequencing_coverage_calculator.html))


No installation is required. The file can be opened directly in any modern web browser.

The browser calculator reproduces the main functionality of `VAF_heatmap.m`:

- accepts user-provided `k0`, `VAF`, `n`, and `pE`
- reports the expected variant detection rate
- displays a heatmap of detection rate across preset VAF and coverage values

Example input:

```text
k0 = 1
VAF = 0.01
n = 100
pE = 0.001
```

Example output:

```text
Variant detection rate = 0.645612
```

or:

```text
64.5612%
```

## MATLAB Scripts

### `modelFitting.m`

Fits and compares analytical detection-rate models against observed recall-rate data from HapMap mixture truth sets.

The script includes hard-coded recall-rate vectors for four sequencing centers:

- `TRUE_WASHU`
- `TRUE_BCM`
- `TRUE_NYGC`
- `TRUE_Broad`

Each vector contains recall rates across VAF levels such as:

```text
10%, 5%, 2%, 1%, 0.5%, 0.25%
```

The script compares detection-rate model configurations including:

- binomial model
- binomial model with sequencing error
- binomial model with Poisson-distributed sequencing error

It also computes summary statistics such as SSE and RMSE, then generates model-comparison plots. These outputs are used to evaluate which analytical model best matches the observed recall-rate data.

### `VAF_heatmap.m`

Calculates variant detection probability for a user-provided sequencing scenario and generates a heatmap across preset VAF and sequencing coverage values.

This script is the MATLAB source for the browser calculator.

## How `VAF_heatmap.m` Works

### 1. Collect user inputs

The script asks the user to enter:

- `k0`: minimum number of supporting reads required to confirm variant existence
- `VAF`: variant allele frequency, entered as a fraction
- `n`: total read depth
- `pE`: sequencing error rate

Example:

```text
Enter how many supporting read are required (k0): 1
Enter the variant allele frequency (VAF): 0.01
Enter the total read depth (n): 100
Enter the sequencing error rate (pE): 0.001
```

Here, `VAF = 0.01` means 1%, and `pE = 0.001` means 0.1%.

If the user presses Enter without providing a value, the script uses defaults:

```text
k0 = 1
VAF = 0.01
n = 500
pE = 0.001
```

### 2. Validate inputs

The script checks that:

- `k0` is a positive integer
- `VAF` is between 0 and 1
- `n` is a positive integer
- `pE` is between 0 and 1

If an input is outside the expected range, MATLAB prints a warning.

### 3. Define heatmap axes

The script evaluates detection rate across preset VAF and coverage values:

```text
VAF:      10%, 5%, 2.5%, 1%, 0.5%, 0.25%, 0.1%
Coverage: 30X, 60X, 100X, 200X, 300X, 400X, 500X, 600X, 1000X
```

In the MATLAB code, VAF is stored as fractions:

```matlab
VAF = [0.1 0.05 0.025 0.01 0.005 0.0025 0.001];
```

Coverage is stored as read depth:

```matlab
n_seq = [30 60 100 200 300 400 500 600 1000];
```

### 4. Adjust the supporting-read probability for sequencing error

For each read, the probability that it supports the variant is modeled as:

```text
q = VAF * (1 - pE) + (1 - VAF) * pE / 3
```

This has two parts:

- `VAF * (1 - pE)`: the variant is truly present and the read is not affected by sequencing error
- `(1 - VAF) * pE / 3`: the variant is not truly present in that read, but sequencing error creates one of the three possible non-reference bases

The `/3` assumes sequencing errors are distributed equally among the three non-reference bases.

### 5. Calculate detection probability

The model treats the number of supporting reads as a binomial random variable:

```text
X ~ Binomial(n, q)
```

where:

- `n` is total read depth
- `q` is the probability that one read supports the variant

The detection rate is the probability of observing at least `k0` supporting reads:

```text
P(detection) = P(X >= k0)
```

In the MATLAB script, this is calculated by summing the binomial probabilities from `k0` to `n`:

```matlab
for k = k0:n_temp
    p_detect = p_detect + ...
        nchoosek(n_temp,k) * q.^k .* (1-q).^(n_temp-k);
end
```

The browser calculator uses the same model, but computes the equivalent value as:

```text
P(detection) = 1 - P(X < k0)
```

This is numerically convenient and produces the same result.

### 6. Generate the heatmap

For every combination of preset VAF and coverage, the script calculates detection probability and stores it in `RR`.

The values are multiplied by 100:

```matlab
RR = RR * 100;
```

Then MATLAB displays the matrix as a heatmap:

```matlab
h = heatmap(RR');
```

The transpose, `RR'`, is used so VAF appears on the y-axis and coverage appears on the x-axis.

### 7. Report the user-specific detection rate

After generating the heatmap, the script repeats the same binomial calculation using the user-provided `VAF_input` and `n_input`.

It prints the expected detection rate:

```text
======== Variant detection rate ========
Variant detection rate will be: 0.645612
========================================
```

## Requirements for MATLAB Version

- MATLAB
- MATLAB functions: `nchoosek`, `heatmap`, `bar`, `plot`, `input`, and standard plotting functions
- `customcolormap_preset.m`
- `customcolormap.m`

The heatmap script uses:

```matlab
J = customcolormap_preset('red-white-blue');
```

If `customcolormap_preset.m` is not available or not on the MATLAB path, add it to the path or replace the line with a MATLAB built-in colormap, for example:

```matlab
h.Colormap = parula;
```

## How to Run the MATLAB Scripts

Open MATLAB, navigate to the folder containing the scripts, and run:

```matlab
modelFitting
```

or:

```matlab
VAF_heatmap
```

For `VAF_heatmap-2.m`, MATLAB may display the runnable script name as `VAF_heatmap` if the file is renamed to `VAF_heatmap.m`.

## Files

```text
outputs/VAF_heatmap_calculator.html
```

Self-contained browser calculator.

```text
README.md
```

Repository documentation.

Original MATLAB files:

```text
modelFitting.m
VAF_heatmap.m or VAF_heatmap-2.m
customcolormap.m
customcolormap_preset.m
```

## Citation

If you use these scripts or results generated from them, please cite:

```text
A Pangenomic Method for Establishing a Somatic Variant Detection Resource in HapMap Mixtures
```
