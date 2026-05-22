import numpy as np
import pandas as pd
import functions as fn
import time

T0 = time.time()
def stamp(msg):
    dt = time.time() - T0
    print(f"[{dt:7.1f}s] {msg}", flush=True)

# -----------------------------
# Load data
# -----------------------------
stamp("Loading data")

filePath = "/pscratch/sd/s/suhas31/Data/BGS_data.parquet"
data = pd.read_parquet(filePath)

Z   = data["Z"].values
M_r = data["M_r"].values
gr  = data["gr"].values
rho = data["LogMstar_N_2Mpc"].values   # environment

# -----------------------------
# Bin definitions
# -----------------------------
mrMin, mrMax = -24.0, -14.0
grMin, grMax = -0.5, 1.75

numBins = 500
grBins = np.linspace(grMin, grMax, numBins + 1)
mrBins = np.linspace(mrMin, mrMax, numBins + 1)

mrStep = 0.2
mrEdges = np.arange(mrMin, mrMax + mrStep, mrStep)
mrCenters = (mrEdges[:-1] + mrEdges[1:]) / 2

# -----------------------------
# Environment bins
# -----------------------------
mask0 = rho != 0

# 9 equal percentile bins (nonzero only)
percentiles = np.linspace(0, 100, 10)  # 0,11.1,...,100
p_edges = np.percentile(rho[mask0], percentiles)

nEnv = 10  # 1 (zero) + 9 percentile bins

# -----------------------------
# Containers
# -----------------------------
h2d_list = []
divider_points = []
divider_mr = []
alpha_vals = []
counts = []

stamp("Starting calculations")

for j in range(nEnv):

    stamp(f"Processing environment bin {j+1}/{nEnv}")

    # -----------------------------
    # Select environment bin
    # -----------------------------
    if j == 0:
        maskEnv = (rho == 0)
    else:
        lo = p_edges[j-1]
        hi = p_edges[j]

        if j == 1:
            maskEnv = (rho > 0) & (rho <= hi)
        else:
            maskEnv = (rho > lo) & (rho <= hi)

    Fig_Mr = M_r[maskEnv]
    Fig_gr = gr[maskEnv]

    counts.append(len(Fig_gr))

    # -----------------------------
    # Histogram
    # -----------------------------
    h2d, _, _ = np.histogram2d(
        Fig_gr,
        Fig_Mr,
        bins=[grBins, mrBins]
    )

    h2d_list.append(h2d)

    # -----------------------------
    # Divider calculation
    # -----------------------------
    dividers = []
    mrVals = []
    weights = []

    for mrCenter in mrCenters:

        maskMr = (
            (Fig_Mr >= mrCenter - mrStep/2) &
            (Fig_Mr <  mrCenter + mrStep/2)
        )

        Nslice = np.count_nonzero(maskMr)
        grVals = Fig_gr[maskMr]

        try:
            binCenters, counts_hist, dgFit, gauss1, gauss2, popt = \
                fn.fit_double_gauss(grVals, bins=np.linspace(grMin, grMax, 101))

            xDiv, _ = fn.CR_div(binCenters, gauss1, gauss2)

            if xDiv is not None and np.isfinite(xDiv):
                dividers.append(xDiv)
                mrVals.append(mrCenter)
                weights.append(Nslice)

        except:
            continue

    divider_points.append(np.array(dividers))
    divider_mr.append(np.array(mrVals))

    if len(dividers) >= 2:
        coefs = np.polyfit(mrVals, dividers, deg=1, w=weights)
        alpha_vals.append(abs(coefs[0]))
    else:
        alpha_vals.append(np.nan)

# -----------------------------
# Save cache
# -----------------------------
stamp("Saving cache")

np.savez(
    "/pscratch/sd/s/suhas31/numpyFiles/cmdEnv_cache.npz",

    h2d_list = np.array(h2d_list),
    divider_points = np.array(divider_points, dtype=object),
    divider_mr = np.array(divider_mr, dtype=object),

    alpha_vals = np.array(alpha_vals),
    counts = np.array(counts),

    grBins = grBins,
    mrBins = mrBins,
    p_edges = p_edges
)

stamp("Done")