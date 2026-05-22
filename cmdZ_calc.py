import numpy as np
import pandas as pd
import functions as fn
import time

T0 = time.time()
def stamp(msg):
    dt = time.time() - T0
    print(f"[{dt:7.1f}s] {msg}", flush=True)

stamp("Loading data")

filePath = "/pscratch/sd/s/suhas31/Data/BGS_data.parquet"
data = pd.read_parquet(filePath)

Z   = data["Z"].values
M_r = data["M_r"].values
gr  = data["gr"].values

# -----------------------------
# Bin definitions
# -----------------------------

z_min, z_max, z_step = 0.0, 0.5, 0.025
zBins = np.arange(z_min, z_max + z_step, z_step)

mrMin, mrMax = -24.0, -14.0
grMin, grMax = -0.5, 1.75

numBins = 500
grBins = np.linspace(grMin, grMax, numBins + 1)
mrBins = np.linspace(mrMin, mrMax, numBins + 1)

mrStep = 0.2
mrEdges = np.arange(mrMin, mrMax + mrStep, mrStep)
mrCenters = (mrEdges[:-1] + mrEdges[1:]) / 2

# -----------------------------
# Containers
# -----------------------------

h2d_list = []
divider_points = []
divider_mr = []
alpha_vals = []
counts = []

stamp("Starting calculations")

for iz in range(len(zBins) - 1):

    z_lo = zBins[iz]
    z_hi = zBins[iz + 1]

    stamp(f"Processing {z_lo:.3f} < z ≤ {z_hi:.3f}")

    maskZ = (Z > z_lo) & (Z <= z_hi)

    FigZ_Mr = M_r[maskZ]
    FigZ_gr = gr[maskZ]

    counts.append(len(FigZ_gr))

    # -----------------------------
    # Histogram
    # -----------------------------

    h2d, _, _ = np.histogram2d(
        FigZ_gr,
        FigZ_Mr,
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
            (FigZ_Mr >= mrCenter - mrStep/2) &
            (FigZ_Mr <  mrCenter + mrStep/2)
        )

        Nslice = np.count_nonzero(maskMr)
        grVals = FigZ_gr[maskMr]

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
    "/pscratch/sd/s/suhas31/numpyFiles/cmdZ_cache.npz",

    h2d_list = np.array(h2d_list),
    divider_points = np.array(divider_points, dtype=object),
    divider_mr = np.array(divider_mr, dtype=object),

    alpha_vals = np.array(alpha_vals),
    counts = np.array(counts),

    zBins = zBins,
    grBins = grBins,
    mrBins = mrBins
)

stamp("Done")