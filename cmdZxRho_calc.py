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
LogMstar_N_2Mpc = data["LogMstar_N_2Mpc"].values


zMin, zMax, zStep = 0.0, 0.5, 0.025
zEdges = np.arange(zMin, zMax + zStep, zStep)
nZ = len(zEdges) - 1

mrMin, mrMax, mrStep = -24., -14., 0.2
mrEdges   = np.arange(mrMin, mrMax + mrStep, mrStep)
mrCenters = 0.5 * (mrEdges[:-1] + mrEdges[1:])

grMin, grMax = -0.5, 1.75
grEdges = np.linspace(grMin, grMax, 101)

numBins = 500
grBins = np.linspace(grMin, grMax, numBins + 1)
mrBins = np.linspace(mrMin, mrMax, numBins + 1)

nEnv = 4


# STORAGE
hist_cache = np.zeros((nZ, nEnv, numBins, numBins))
divider_x_cache = [[None]*nEnv for _ in range(nZ)]
divider_mr_cache = [[None]*nEnv for _ in range(nZ)]
alpha_cache = np.zeros((nZ, nEnv))
counts_cache = np.zeros((nZ, nEnv))


stamp("Starting calculations")

for iz in range(nZ):

    zlo, zhi = zEdges[iz], zEdges[iz+1]

    stamp(f"z bin {iz+1}/{nZ} : {zlo:.3f} < z ≤ {zhi:.3f}")

    mZ = (zlo < Z) & (Z <= zhi)

    gr_z  = gr[mZ]
    Mr_z  = M_r[mZ]
    rho_z = LogMstar_N_2Mpc[mZ]

    mask0_z = rho_z != 0

    p33 = np.percentile(rho_z[mask0_z], 33.3333333333)
    p67 = np.percentile(rho_z[mask0_z], 66.6666666667)


    for j in range(nEnv):

        if j == 0:

            gr_cell = gr_z[~mask0_z]
            Mr_cell = Mr_z[~mask0_z]

        else:

            rho_nz = rho_z[mask0_z]
            gr_nz  = gr_z[mask0_z]
            Mr_nz  = Mr_z[mask0_z]

            if j == 1:
                mEnv = rho_nz <= p33
            elif j == 2:
                mEnv = (p33 < rho_nz) & (rho_nz <= p67)
            else:
                mEnv = rho_nz > p67

            gr_cell = gr_nz[mEnv]
            Mr_cell = Mr_nz[mEnv]


        counts_cache[iz, j] = len(gr_cell)


        # Histogram
        h2d, _, _ = np.histogram2d(
            gr_cell,
            Mr_cell,
            bins=[grBins, mrBins]
        )

        hist_cache[iz, j] = h2d


        # Divider
        dividers = []
        Mr_vals = []
        wts = []
        
        for Mr_c in mrCenters:
        
            mMr = (Mr_cell >= Mr_c - mrStep/2) & (Mr_cell < Mr_c + mrStep/2)
            n_in = np.count_nonzero(mMr)
        
            grVals = gr_cell[mMr]
        
            try:
                binCenters, counts, dgFit, gauss1, gauss2, popt = \
                    fn.fit_double_gauss(grVals, bins=grEdges)
        
                xDiv, _ = fn.CR_div(binCenters, gauss1, gauss2)
        
                # keep ANY valid numeric result
                if np.isfinite(xDiv):
        
                    dividers.append(xDiv)
                    Mr_vals.append(Mr_c)
                    wts.append(n_in)
        
            except Exception:
                # only skip if the code actually crashes
                continue


        divider_x_cache[iz][j] = np.array(dividers)
        divider_mr_cache[iz][j] = np.array(Mr_vals)

        if len(dividers) >= 2:

            coefs = np.polyfit(Mr_vals, dividers, deg=1, w=wts)
            alpha_cache[iz, j] = abs(coefs[0])

        else:

            alpha_cache[iz, j] = np.nan

# SAVE CACHE
stamp("Saving cache")

np.savez(
    "/pscratch/sd/s/suhas31/numpyFiles/cmdZxRho_cache.npz",

    hist_cache = hist_cache,

    divider_x_cache = np.array(divider_x_cache, dtype=object),
    divider_mr_cache = np.array(divider_mr_cache, dtype=object),

    alpha_cache = alpha_cache,
    counts_cache = counts_cache,

    zEdges = zEdges,
    grBins = grBins,
    mrBins = mrBins
)

stamp("Finished")