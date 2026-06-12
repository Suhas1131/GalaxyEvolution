# Modeling Galaxy Evolution as a Function of Cosmic Time and Local Environment

**Advisor:** Dr. Barbara Ryden  
**Project type:** Undergraduate honors thesis / large-scale data processing and statistical modeling  

This repository contains public analysis code, selected figures, and methodology notes for a project studying how galaxy populations change across cosmic time and local environment using data from the DESI Bright Galaxy Survey.

The project analyzes a cleaned sample of approximately **6.2 million galaxies** and builds a statistical boundary between two major galaxy populations: **blue, star-forming galaxies** and **red, quiescent galaxies**. In astronomy, the transition region between these populations is commonly called the **green valley**. The goal of this project is to measure how that boundary changes with redshift, luminosity, and local environment.

<p align="center">
  <img src="figures/cmd.png" width="650"/>
</p>

---

## Project Summary

Galaxies are not distributed randomly in color–magnitude space. They tend to separate into two broad groups:

- **Blue galaxies:** systems that are typically still forming stars and contain younger stellar populations
- **Red galaxies:** systems that are typically older and less actively star-forming

Between these groups is a transition region known as the **green valley**. Measuring this boundary helps quantify how galaxy populations evolve over time and how that evolution depends on environment.

This project constructs a statistical green-valley divider using double-Gaussian fits to galaxy color distributions. The divider is then tracked across:

- **Redshift:** used as a proxy for cosmic time
- **Absolute magnitude:** a measure of intrinsic galaxy brightness
- **Local environment:** total stellar mass within a 2 Mpc radius around each galaxy

The public version of this repository emphasizes the analysis workflow, statistical methodology, code structure, and selected results. Raw survey catalogs and large intermediate data products are not included.

<sub>1 Mpc, or one megaparsec, is approximately 3.26 million light-years.</sub>

---

## Technical Highlights

- Analyzed a cleaned sample of approximately **6.2 million galaxies** from a larger DESI Bright Galaxy Survey catalog.
- Built Python-based analysis workflows for large tabular astronomical data.
- Engineered features based on redshift, rest-frame color, absolute magnitude, stellar mass, and local environment.
- Modeled galaxy color distributions using double-Gaussian fits.
- Defined a statistical galaxy-population divider using a completeness–reliability optimization metric.
- Compared divider behavior across redshift, luminosity, and local-density bins.
- Used cached intermediate outputs to make repeated analysis and visualization more efficient on NERSC Perlmutter.
- Produced diagnostic plots to validate model fits and inspect population-level trends.

---

## Key Results

### 1. Evolution with Redshift

The green-valley divider changes with redshift, indicating that the boundary between blue and red galaxy populations evolves over cosmic time.

At higher redshift, the divider is generally redder; equivalently, the boundary shifts blueward toward the present day. The size of this shift depends on galaxy luminosity.

<p align="center">
  <img src="figures/dividerEvolution_Z.png" width="650"/>
</p>

---

### 2. Dependence on Local Environment

At fixed redshift, galaxies in denser environments show a redder divider than galaxies in lower-density environments. This effect is stronger for lower-luminosity galaxies and weaker for more luminous galaxies.

This suggests that local environment plays a measurable role in galaxy population evolution, especially for less luminous systems.

<p align="center">
  <img src="figures/dividerEvolution_Rho.png" width="650"/>
</p>

---

### 3. Joint Redshift and Environment Trends

The divider evolves with both redshift and local environment. At fixed environment, the divider changes over cosmic time. At fixed redshift, the divider shifts with increasing local stellar mass density.

<p align="center">
  <img src="figures/dividerEvolution_ZxRho.png" width="750"/>
</p>

---

## Dataset

| Property | Description |
|---|---|
| Survey | DESI Bright Galaxy Survey |
| Initial catalog size | Approximately 6.34 million galaxies before final cleaning |
| Cleaned analysis sample | Approximately 6.2 million galaxies |
| Full redshift range considered | 0 < z ≤ 0.5 |
| Primary results range | 0 < z ≤ 0.3 |
| Main color diagnostic | Rest-frame g−r color |
| Luminosity measure | Rest-frame r-band absolute magnitude |
| Environment metric | Log total stellar mass of neighboring galaxies within 2 Mpc |

---

## Methodology

### Color–Magnitude Diagrams

A color–magnitude diagram compares galaxy color against luminosity. In this project, galaxy color is measured using rest-frame **g−r color**, and luminosity is represented using **r-band absolute magnitude**.

The diagram reveals two broad galaxy populations: a blue population and a red population. The project focuses on measuring the statistical divider between them.

---

### Double-Gaussian Fitting

Within each redshift, magnitude, and environment bin, the galaxy color distribution is modeled as the sum of two Gaussian components:

- One Gaussian for the blue population
- One Gaussian for the red population

This provides a statistical way to identify the two populations without manually drawing a dividing line.

<p align="center">
  <img src="figures/dgFit_Row3.png" width="750"/>
</p>

---

### Completeness–Reliability Optimization

After fitting the blue and red populations, the green-valley divider is selected by maximizing a completeness–reliability metric.

The goal is to place the divider where it best separates the two modeled populations while minimizing cross-contamination between the blue and red groups.

The optimization metric is:

```text
tau = C_B × C_R × R_B × R_R
```

where:

- `C_B` is the completeness of the blue population
- `C_R` is the completeness of the red population
- `R_B` is the reliability of the blue population
- `R_R` is the reliability of the red population

In practical terms, this is a way to choose the divider that best balances classification performance for both galaxy populations.

<p align="center">
  <img src="figures/GVDiv_Row3.png" width="750"/>
</p>

---

### Environment Binning

Local environment is quantified using the total stellar mass of neighboring galaxies within a projected radius of 2 Mpc.

Galaxies are grouped by environment so that the divider can be compared across low-density and high-density regions. This helps test whether galaxy evolution depends only on cosmic time or also on a galaxy’s surroundings.

---

### HPC Workflow

The analysis was run on the Perlmutter supercomputer at the National Energy Research Scientific Computing Center using Python scientific-computing tools.

Intermediate results were cached as `.npz` files to avoid repeatedly recomputing expensive binning, fitting, and visualization steps. These cached outputs are not included in the public repository due to file size.

---

## Repository Structure

```text
├── functions.py          # Core helper functions for Gaussian fitting and divider optimization
├── cmdZ_calc.py          # Color–magnitude divider analysis binned by redshift
├── cmdEnv_calc.py        # Color–magnitude divider analysis binned by local environment
├── cmdZxRho_calc.py      # Joint redshift × environment divider analysis
├── AppendixD.py          # Diagnostic plots for double-Gaussian fits and divider behavior
├── Intro_Methods.ipynb   # Data selection, methodology, and pipeline walkthrough
├── Results.ipynb         # Result figures and analysis interpretation
└── figures/              # Selected public-facing figures from the analysis
```

---

## Dependencies

The core analysis uses:

```text
numpy
pandas
scipy
matplotlib
pyarrow
```

Depending on the execution environment, additional packages may be needed for notebook display or file-path management.

---

## Data Availability and Reproducibility

This repository does **not** include raw DESI catalogs, processed parquet files, cached `.npz` outputs, or private HPC file paths.

The public version is intended to document the analysis workflow, statistical methods, and selected results. Users with access to equivalent DESI Bright Galaxy Survey data products can adapt the scripts by updating the relevant input and output paths.

Because the full processed dataset and cached intermediate files are not included, the repository should be treated as a public research-code and methodology repository rather than a fully self-contained reproducibility package.

---

## Skills Demonstrated

This project demonstrates experience with:

- Large-scale tabular data analysis
- Scientific Python workflows
- Statistical modeling
- Feature engineering
- Model validation and diagnostic plotting
- High-performance computing workflows
- Translating domain-specific research questions into quantitative analysis pipelines
- Communicating technical results to both scientific and non-specialist audiences

---

## Citation

Reddy, S. (2026). *The Evolution of Galaxy Color as a Function of Time and Local Environment*. Undergraduate Honors Thesis, The Ohio State University. Advisor: Prof. Barbara Ryden.