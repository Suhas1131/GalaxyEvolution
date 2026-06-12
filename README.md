# Modeling Galaxy Evolution as a Function of Cosmic Time and Local Environment

**Advisor:** Dr. Barbara Ryden 
**Project type:** Undergraduate honors thesis / large-scale data processing and statistical modeling

This repository contains the public analysis code, selected figures, and methodology notes for a project studying how galaxy populations change across cosmic time and local environment using data from the DESI Bright Galaxy Survey.

The project analyzes millions of galaxies and builds a statistical boundary between two major galaxy populations: **blue, star-forming galaxies** and **red, quiescent galaxies**. In astronomy, the transition region between these populations is commonly called the **green valley**. The goal of this project is to measure how that boundary changes with redshift (cosmic time), luminosity (intrinsic brightness), and local environment.

---

## Project Summary

Galaxies are not distributed randomly in color–magnitude space. They tend to separate into two broad groups:

- **Blue galaxies**: systems that are typically still forming stars and consist of younger stars
- **Red galaxies**: systems which are typically older and less actively star-forming

Between these groups is a transition region known as the **green valley**. Measuring this boundary helps quantify how galaxy populations evolve over time and how that evolution depends on environment.

This project constructs a statistical green-valley divider using double-Gaussian fits to galaxy color distributions. The divider is then tracked across:

- **Redshift**: acts as a proxy for cosmic time
- **Absolute magnitude**: logarithmic representation of galaxy luminosity
- **Local environment**: total stellar mass within a 2 Mpc* radius around a central galaxy

The public version of this repository emphasizes the analysis workflow, statistical methodology, code structure, and selected results. Raw survey catalogs and large intermediate data products are not included.

<sub>*1 Mpc, or one megaparsec, is approximately 3.26 million light-years.</sub>

---

## Technical Highlights

- Analyzed a cleaned sample of approximately **6.2 million galaxies** from a larger DESI Bright Galaxy Survey dataset.
- Built Python-based analysis workflows for large tabular astronomical data.
- Engineered features based on redshift, rest-frame color, absolute magnitude, stellar mass, and local environment.
- Modeled galaxy color distributions using double-Gaussian fits with `scipy.optimize.curve_fit`.
- Defined a statistical galaxy-population divider using a completeness–reliability optimization metric.
- Compared divider behavior across redshift, luminosity, and local-density bins.
- Used cached intermediate outputs to make repeated analysis and visualization more efficient on NERSC Perlmutter.
- Produced diagnostic plots to validate model fits and inspect population-level trends.

---

## Key Results

### 1. Evolution with Redshift

The green-valley divider changes with redshift, indicating that the boundary between blue and red galaxy populations evolves over cosmic time.

At higher redshift, the divider is generally redder; equivalently, the boundary shifts blueward toward the present day. The size of this shift depends on galaxy luminosity.

---

### 2. Dependence on Local Environment

At fixed redshift, galaxies in denser environments show a redder divider than galaxies in lower-density environments. This effect is stronger for lower-luminosity galaxies and weaker for more luminous galaxies.

This suggests that local environment plays a measurable role in galaxy population evolution, especially for less luminous systems.

---

### 3. Joint Redshift and Environment Trends

The divider evolves with both redshift and local environment. At fixed environment, the divider changes over cosmic time. At fixed redshift, the divider shifts with increasing local stellar mass density.

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

The diagram reveals two broad galaxy populations: a blue population and a red population. The project focuses on measuring the divider between them.

---

### Double-Gaussian Fitting

Within each redshift, magnitude, and environment bin, the galaxy color distribution is modeled as the sum of two Gaussian components:

- One Gaussian for the blue population
- One Gaussian for the red population

This provides a statistical way to identify the two populations without manually drawing a dividing line.

---

### Completeness–Reliability Optimization

After fitting the blue and red populations, the green-valley divider is selected by maximizing a completeness–reliability metric.

The goal is to place the divider where it best separates the two modeled populations while minimizing cross-contamination.

The optimization metric is:

```text
tau = C_B × C_R × R_B × R_R