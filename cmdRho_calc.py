"""
Calculate CMD density grids and green-valley dividers across environment bins.

This script reads the processed BGS dataset, computes color-magnitude density histograms in local-environment bins, fits double-Gaussian color distributions across absolute-magnitude slices, and saves the cached outputs for plotting.

Run DataProcessing.py and EnvironmentFeatures.py before running this script.
"""

# ----------------------------- LIBRARIES -----------------------------

import time
import numpy as np
import pandas as pd
import functions as fn

# ----------------------------- PATHS -----------------------------

InputPath = "data/BGS_data.parquet"
OutputPath = "data/cmdEnv_cache.npz"

# ----------------------------- TIMER -----------------------------

T0 = time.time()

def Stamp(Message):
    """
    Print elapsed runtime and status.
    """

    TimeElapsed = time.time() - T0
    print(f"[{TimeElapsed:7.1f}s] {Message}", flush=True)

# ----------------------------- LOAD DATA -----------------------------

Stamp("Loading data")

Data = pd.read_parquet(InputPath)

Z = Data["Z"].to_numpy()
M_r = Data["M_r"].to_numpy()
gr = Data["gr"].to_numpy()
rho = Data["LogMstar_N_2Mpc"].to_numpy()

# ----------------------------- BIN DEFINITIONS -----------------------------

mrMin, mrMax = -24.0, -14.0
grMin, grMax = -0.5, 1.75

numBins = 500
grBins = np.linspace(grMin, grMax, numBins + 1)
mrBins = np.linspace(mrMin, mrMax, numBins + 1)

mrStep = 0.2
mrEdges = np.arange(mrMin, mrMax + mrStep, mrStep)
mrCenters = (mrEdges[:-1] + mrEdges[1:]) / 2

# ----------------------------- ENVIRONMENT BINS -----------------------------

mask0 = rho != 0

# 9 equal percentile bins for galaxies with nonzero neighboring stellar mass
percentiles = np.linspace(0, 100, 10)
pEdges = np.percentile(rho[mask0], percentiles)

nEnv = 10  # 1 no neighbor bin + 9 percentile bins

# ----------------------------- CONTAINERS -----------------------------

h2dList = []  # Stores 2D histogram data
dividerPoints = []  # Stores x-coordinates of dividers
dividerMr = []  # Stores y-coordinates of dividers
alphaVals = []  # Stores the slopes of the linear divider
counts = []  # Stores the number of galaxies in each environment bin

# ----------------------------- CALCULATIONS -----------------------------

for j in range(nEnv):
    Stamp(f"Processing environment bin {j + 1}/{nEnv}")

    # Select environment bin
    if j == 0:
        maskEnv = rho == 0
    else:
        lo = pEdges[j - 1]
        hi = pEdges[j]

        if j == 1:
            maskEnv = (rho > 0) & (rho <= hi)
        else:
            maskEnv = (rho > lo) & (rho <= hi)

    figMr = M_r[maskEnv]
    figGr = gr[maskEnv]

    counts.append(len(figGr))

    # Build CMD 2D histogram
    h2d, _, _ = np.histogram2d(
        figGr,
        figMr,
        bins=[grBins, mrBins]
    )

    h2dList.append(h2d)

    # Calculate dividers across absolute-magnitude slices
    dividers = []
    mrVals = []
    weights = []

    for mrCenter in mrCenters:
        maskMr = (figMr >= mrCenter - mrStep / 2) & (figMr < mrCenter + mrStep / 2)

        nSlice = np.count_nonzero(maskMr)
        grVals = figGr[maskMr]

        try:
            binCenters, countsHist, dgFit, gauss1, gauss2, popt = fn.fitDoubleGauss(
                grVals,
                bins=np.linspace(grMin, grMax, 101)
            )

            xDiv, _ = fn.CR_Div(binCenters, gauss1, gauss2)

            if np.isfinite(xDiv):
                dividers.append(xDiv)
                mrVals.append(mrCenter)
                weights.append(nSlice)

        except Exception:
            continue

    dividerPoints.append(np.array(dividers))
    dividerMr.append(np.array(mrVals))

    if len(dividers) >= 2:
        coefs = np.polyfit(mrVals, dividers, deg=1, w=weights)
        alphaVals.append(abs(coefs[0]))
    else:
        alphaVals.append(np.nan)

# ----------------------------- SAVE CACHE -----------------------------

Stamp("Saving cache")

np.savez(
    OutputPath,
    h2dList=np.array(h2dList),
    dividerPoints=np.array(dividerPoints, dtype=object),
    dividerMr=np.array(dividerMr, dtype=object),
    alphaVals=np.array(alphaVals),
    counts=np.array(counts),
    grBins=grBins,
    mrBins=mrBins,
    pEdges=pEdges
)

Stamp("Done")