"""
Generate diagnostic divider plots across redshift and environment bins.

This script reads the processed BGS dataset, fits double-Gaussian color
distributions across absolute-magnitude slices, selects green-valley dividers,
and saves one 50-panel diagnostic figure for each redshift/environment bin.

Run DataProcessing.py and EnvironmentFeatures.py before running this script.
"""

# ----------------------------- LIBRARIES -----------------------------

import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import functions as fn

# ----------------------------- PATHS -----------------------------

InputPath = "data/BGS_data.parquet"

# Update this path before running
OutputDir = "insert/path/to/scratch/folder/DividerFigures_ZxRho"

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
nMr = len(mrCenters)

grMin, grMax = -0.5, 1.75
grEdges = np.linspace(grMin, grMax, 101)

# ----------------------------- ENVIRONMENT LABELS -----------------------------

envFolderNames = [
    "NoNeighbours",
    "LowDensity",
    "MidDensity",
    "HighDensity",
]

envTitles = [
    "No neighbours",
    "Low density",
    "Mid density",
    "High density",
]

# Create output folders
for Folder in envFolderNames:
    os.makedirs(os.path.join(OutputDir, Folder), exist_ok=True)

# ----------------------------- FIGURE GENERATION -----------------------------

Stamp("Starting calculations")

for iz in range(nZ):

    # Select redshift slice
    zLo = zEdges[iz]
    zHi = zEdges[iz + 1]

    Stamp(f"z bin {iz + 1}/{nZ}: {zLo:.3f} < z <= {zHi:.3f}")

    maskZ = (Z > zLo) & (Z <= zHi)

    figZMr = M_r[maskZ]
    figZgr = gr[maskZ]
    figZrho = rho[maskZ]

    # Split nonzero-density galaxies into three percentile bins
    mask0Z = figZrho != 0

    if np.count_nonzero(mask0Z) > 0:
        p33 = np.percentile(figZrho[mask0Z], 33.3333333333)
        p67 = np.percentile(figZrho[mask0Z], 66.6666666667)
    else:
        p33 = np.nan
        p67 = np.nan

    for j in range(4):

        Stamp(f"  Environment: {envTitles[j]}")

        # Select environment bin
        if j == 0:
            grCell = figZgr[~mask0Z]
            mrCell = figZMr[~mask0Z]

        else:
            rhoNz = figZrho[mask0Z]
            grNz = figZgr[mask0Z]
            mrNz = figZMr[mask0Z]

            if j == 1:
                maskEnv = rhoNz <= p33
            elif j == 2:
                maskEnv = (p33 < rhoNz) & (rhoNz <= p67)
            else:
                maskEnv = rhoNz > p67

            grCell = grNz[maskEnv]
            mrCell = mrNz[maskEnv]

        # One panel for each absolute-magnitude slice
        fig, axs = plt.subplots(
            10, 5,
            figsize=(40, 50),
            sharex=True,
            sharey=True
        )

        axs = axs.ravel()

        for i in range(nMr):

            # Select absolute-magnitude slice
            mrCenter = mrCenters[i]
            mrLo = mrCenter - mrStep / 2
            mrHi = mrCenter + mrStep / 2

            maskMr = (mrCell >= mrLo) & (mrCell < mrHi)
            figData = grCell[maskMr]

            ax = axs[i]

            ax.set_title(
                f"${mrLo:.1f} \\leq M_r < {mrHi:.1f}$",
                fontsize=32
            )

            ax.grid(True, which="both", alpha=0.5, linewidth=0.5)

            # Skip empty slices
            if len(figData) == 0:
                ax.text(
                    0.02, 0.95,
                    f"N = {len(figData)}",
                    transform=ax.transAxes,
                    ha="left", va="top",
                    fontsize=32
                )
                ax.text(
                    0.5, 0.5,
                    "Fit Failed",
                    transform=ax.transAxes,
                    ha="center", va="center",
                    fontsize=28,
                    color="red",
                    weight="bold"
                )

                ax.set_xlim(grMin, grMax)
                ax.set_ylim(0, 1.25)
                ax.set_xticks([-0.5, 0.0, 0.5, 1.0, 1.5])
                ax.set_yticks([0.0, 0.5, 1.0])

                continue

            counts, _ = np.histogram(figData, bins=grEdges)
            peak = counts.max()

            if peak == 0:
                ax.text(
                    0.02, 0.95,
                    f"N = {len(figData)}",
                    transform=ax.transAxes,
                    ha="left", va="top",
                    fontsize=32
                )
                ax.text(
                    0.5, 0.5,
                    "Fit Failed",
                    transform=ax.transAxes,
                    ha="center", va="center",
                    fontsize=28,
                    color="red",
                    weight="bold"
                )

                ax.set_xlim(grMin, grMax)
                ax.set_ylim(0, 1.25)
                ax.set_xticks([-0.5, 0.0, 0.5, 1.0, 1.5])
                ax.set_yticks([0.0, 0.5, 1.0])

                continue

            ax.hist(
                figData,
                bins=grEdges,
                weights=np.ones_like(figData) / peak,
                alpha=0.3,
                color="gray"
            )

            try:
                # Fit double Gaussian model
                binCenters, countsFit, dgFit, gauss1, gauss2, popt = fn.fitDoubleGauss(
                    figData,
                    bins=grEdges
                )

                if np.sum(gauss1) > 0 and np.sum(gauss2) > 0:

                    gauss1 = gauss1 / peak
                    gauss2 = gauss2 / peak
                    dgFit = dgFit / peak

                    # Keep blue and red components consistent
                    mean1 = popt[1]
                    mean2 = popt[4]

                    if mean1 <= mean2:
                        blue = gauss1
                        red = gauss2
                    else:
                        blue = gauss2
                        red = gauss1

                    ax.plot(binCenters, blue, "--", color="blue")
                    ax.plot(binCenters, red, "--", color="red")
                    ax.plot(binCenters, dgFit, "-", color="black")

                    # Calculate and plot green-valley divider
                    xDiv, tau = fn.CR_Div(binCenters, blue, red)

                    if np.isfinite(xDiv):
                        ax.axvline(xDiv, color="green")

                    if tau is None or not np.isfinite(tau):
                        tauText = "NaN"
                    else:
                        tauText = f"{tau:.2f}"

                    ax.text(
                        0.02, 0.95,
                        f"N = {len(figData)}\n$\\tau$ = {tauText}",
                        transform=ax.transAxes,
                        ha="left", va="top",
                        fontsize=32
                    )

                else:
                    ax.text(
                        0.02, 0.95,
                        f"N = {len(figData)}",
                        transform=ax.transAxes,
                        ha="left", va="top",
                        fontsize=32
                    )
                    ax.text(
                        0.5, 0.5,
                        "Fit Failed",
                        transform=ax.transAxes,
                        ha="center", va="center",
                        fontsize=28,
                        color="red",
                        weight="bold"
                    )

            except Exception as Error:
                print(
                    f"Bin {i} failed | env={envTitles[j]} | "
                    f"z={zLo:.3f}-{zHi:.3f} | N={len(figData)} | error={repr(Error)}",
                    flush=True
                )

                ax.text(
                    0.02, 0.95,
                    f"N = {len(figData)}",
                    transform=ax.transAxes,
                    ha="left", va="top",
                    fontsize=32
                )
                ax.text(
                    0.5, 0.5,
                    "Fit Failed",
                    transform=ax.transAxes,
                    ha="center", va="center",
                    fontsize=28,
                    color="red",
                    weight="bold"
                )

            # Scale axes
            ax.set_xlim(grMin, grMax)
            ax.set_ylim(0, 1.25)

            ax.set_xticks([-0.5, 0.0, 0.5, 1.0, 1.5])
            ax.set_yticks([0.0, 0.5, 1.0])

        # Shared figure labels
        fig.supxlabel(r"$\mathbf{g-r}$", fontsize=34)
        fig.supylabel("Normalized Counts", fontsize=34, weight="bold")

        fig.tight_layout(rect=[0.02, 0.01, 0.98, 0.98])

        # Save figure
        saveDir = os.path.join(OutputDir, envFolderNames[j])
        fileName = os.path.join(
            saveDir,
            f"Z{zLo:.3f}-{zHi:.3f}_{envFolderNames[j]}.pdf"
        )

        plt.savefig(fileName)
        plt.close(fig)

        Stamp(f"    Saved {fileName}")