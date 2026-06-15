"""
Calculate CMD density grids and green-valley dividers across redshift and environment bins.

This script reads the processed BGS dataset, computes color-magnitude density histograms in joint redshift/environment bins, fits double-Gaussian color distributions across absolute-magnitude slices, and saves the cached outputs for plotting.

Run DataProcessing.py and EnvironmentFeatures.py before running this script.
"""

# ----------------------------- LIBRARIES -----------------------------

import time
import numpy as np
import pandas as pd
import functions as fn

# ----------------------------- PATHS -----------------------------

InputPath = "data/BGS_data.parquet"
OutputPath = "data/cmdZxRho_cache.npz"

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

zMin, zMax, zStep = 0.0, 0.5, 0.025
zEdges = np.arange(zMin, zMax + zStep, zStep)
nZ = len(zEdges) - 1

mrMin, mrMax, mrStep = -24.0, -14.0, 0.2
mrEdges = np.arange(mrMin, mrMax + mrStep, mrStep)
mrCenters = 0.5 * (mrEdges[:-1] + mrEdges[1:])

grMin, grMax = -0.5, 1.75
grEdges = np.linspace(grMin, grMax, 101)

numBins = 500
grBins = np.linspace(grMin, grMax, numBins + 1)
mrBins = np.linspace(mrMin, mrMax, numBins + 1)

nEnv = 4

# ----------------------------- CONTAINERS -----------------------------

histCache = np.zeros((nZ, nEnv, numBins, numBins))
dividerXCache = [[None] * nEnv for _ in range(nZ)]
dividerMrCache = [[None] * nEnv for _ in range(nZ)]
alphaCache = np.zeros((nZ, nEnv))
countsCache = np.zeros((nZ, nEnv))

# ----------------------------- CALCULATIONS -----------------------------

for iz in range(nZ):

    zLo = zEdges[iz]
    zHi = zEdges[iz + 1]

    Stamp(f"z bin {iz + 1}/{nZ} : {zLo:.3f} < z <= {zHi:.3f}")

    maskZ = (zLo < Z) & (Z <= zHi)

    grZ = gr[maskZ]
    mrZ = M_r[maskZ]
    rhoZ = rho[maskZ]

    mask0Z = rhoZ != 0

    p33 = np.percentile(rhoZ[mask0Z], 33.3333333333)
    p67 = np.percentile(rhoZ[mask0Z], 66.6666666667)

    for j in range(nEnv):
        if j == 0:
            grCell = grZ[~mask0Z]
            mrCell = mrZ[~mask0Z]
        else:
            rhoNz = rhoZ[mask0Z]
            grNz = grZ[mask0Z]
            mrNz = mrZ[mask0Z]

            if j == 1:
                maskEnv = rhoNz <= p33
            elif j == 2:
                maskEnv = (p33 < rhoNz) & (rhoNz <= p67)
            else:
                maskEnv = rhoNz > p67

            grCell = grNz[maskEnv]
            mrCell = mrNz[maskEnv]

        countsCache[iz, j] = len(grCell)

        # Build CMD 2D histogram
        h2d, _, _ = np.histogram2d(
            grCell,
            mrCell,
            bins=[grBins, mrBins]
        )

        histCache[iz, j] = h2d

        # Calculate dividers across absolute-magnitude slices
        dividers = []
        mrVals = []
        weights = []
        
        for mrCenter in mrCenters:
            maskMr = (mrCell >= mrCenter - mrStep / 2) & (mrCell < mrCenter + mrStep / 2)
            nSlice = np.count_nonzero(maskMr)
        
            grVals = grCell[maskMr]
        
            try:
                binCenters, countsHist, dgFit, gauss1, gauss2, popt = fn.fitDoubleGauss(
                    grVals,
                    bins=grEdges
                )
        
                xDiv, _ = fn.CR_Div(binCenters, gauss1, gauss2)
        
                if np.isfinite(xDiv):
                    dividers.append(xDiv)
                    mrVals.append(mrCenter)
                    weights.append(nSlice)
        
            except Exception:
                continue

        dividerXCache[iz][j] = np.array(dividers)
        dividerMrCache[iz][j] = np.array(mrVals)

        if len(dividers) >= 2:
            coefs = np.polyfit(mrVals, dividers, deg=1, w=weights)
            alphaCache[iz, j] = abs(coefs[0])
        else:
            alphaCache[iz, j] = np.nan

# ----------------------------- SAVE CACHE -----------------------------

Stamp("Saving cache")

np.savez(
    OutputPath,

    histCache=histCache,
    dividerXCache=np.array(dividerXCache, dtype=object),
    dividerMrCache=np.array(dividerMrCache, dtype=object),
    alphaCache=alphaCache,
    countsCache=countsCache,
    zEdges=zEdges,
    grBins=grBins,
    mrBins=mrBins
)

Stamp("Finished")