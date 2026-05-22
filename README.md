# Galaxy Colour Evolution — DESI BGS

**Advisor: Dr. Barbara Ryden**

Analysis of galaxy colour bimodality across cosmic time and local environment using 6.2 million galaxies from the DESI Bright Galaxy Survey.

<p align="center">
  <img src="figures/cmd.png" width="600"/>
</p>

---

## Overview

Galaxies fall into two main populations in colour–magnitude space: blue star-forming galaxies and red quiescent galaxies, separated by the *green valley*. This project constructs a statistical green-valley divider and tracks how it evolves with both redshift (cosmic time) and local stellar mass density (environment) over a lookback time of ~5 Gyr.

The divider is defined independently in each redshift, environment, and magnitude bin using double-Gaussian fits to the (g−r) colour distribution, optimized via a completeness–reliability metric.

---

## Key Results

- **Redshift evolution:** The divider shifts ~0.17 mag blueward at M_r = −21 and ~0.05 mag at M_r = −18 over 0 < z ≤ 0.3, indicating that the blue/red boundary was redder in the past.

<p align="center">
  <img src="figures/dividerEvolution_Z.png" width="600"/>
</p>

- **Environmental dependence:** At fixed redshift, the divider is ~0.15 mag redder in high-density vs. low-density environments at M_r = −18, with weaker dependence (~0.02 mag) at M_r = −21 — lower-luminosity galaxies show stronger environmental sensitivity.

<p align="center">
  <img src="figures/dividerEvolution_Rho.png" width="600"/>
</p>

- **Joint evolution:** At fixed environment the divider evolves blueward with redshift. At fixed redshift the divider shifts subtly redward with increasing local stellar mass density, with the effect most apparent at higher redshifts.

<p align="center">
  <img src="figures/dividerEvolution_ZxRho.png"/>
</p>

---

## Dataset

| Property | Value |
|---|---|
| Survey | DESI Bright Galaxy Survey (BGS) |
| Sample size | 6,192,490 galaxies |
| Redshift range | 0 < z ≤ 0.5 (~5 Gyr lookback time) |
| Environment metric | Log total stellar mass within 2 Mpc |
| Colour diagnostic | Rest-frame (g−r) |

---

## Methods

<p align="center">
  <img src="figures/dgFit_Row3.png" width="700"/>
</p>

**Double-Gaussian fitting** — Each (g−r) colour distribution slice is modelled as a sum of two Gaussians representing the blue and red populations using `scipy.optimize.curve_fit`.

<p align="center">
  <img src="figures/GVDiv_Row3.png" width="700"/>
</p>

**Completeness–Reliability optimization** — The green-valley divider is placed at the colour that maximizes $\tau$ :

$\tau = C_B\ ×\ C_R\ ×\ R_B\ ×\ R_R$

where $C_B$ and $C_R$ are the completeness of the blue and red sequences, and $R_B$ and $R_R$ are their respective reliabilities. This is equivalent to maximizing the geometric mean of a 2×2 classification matrix.

**Environment binning** — Local environment is quantified as the log total stellar mass of neighbouring galaxies within a 2 Mpc projected radius. Galaxies are split into: no neighbours, low-density (0–33rd percentile), mid-density (33rd–67th), and high-density (67th–100th).

**HPC pipeline** — All computation was run on the Perlmutter supercomputer at NERSC using Python (NumPy, pandas, SciPy, Matplotlib) with results cached as `.npz` files for downstream visualization.

---

## Repository Structure

```
├── functions.py          # Core library: Gaussian fitting, completeness-reliability optimizer
├── cmdZ_calc.py          # CMD analysis binned by redshift (20 bins, 0 < z ≤ 0.5)
├── cmdEnv_calc.py        # CMD analysis binned by environment density (10 percentile bins)
├── cmdZxRho_calc.py      # Joint redshift × environment CMD analysis (20 × 4 bins)
├── AppendixD.py          # Diagnostic plots: double-Gaussian fits per bin
├── Intro_Methods.ipynb   # Data selection, methodology, and pipeline walkthrough
├── Results.ipynb         # Figures and results from the analysis
└── figures/              # All plots from the main text
```

---

## Dependencies

```
numpy
pandas
scipy
matplotlib
pyarrow
```

---

## Data Availability

The raw DESI BGS data used in this analysis may be available through the DESI data releases. All computation was run on the Perlmutter supercomputer at NERSC. The processed .parquet file and cached .npz outputs are not included in this repository due to file size.

---

## Citation

> Reddy, S. (2026). *The Evolution of Galaxy Colour as a Function of Time and Local Environment*. Undergraduate Honors Thesis, The Ohio State University. Advisor: Prof. Barbara Ryden.
