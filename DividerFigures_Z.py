"""
Generate diagnostic divider plots across redshift bins.

This script reads the processed BGS dataset, fits double-Gaussian color distributions across absolute-magnitude slices, selects green-valley dividers, and saves one 50-panel diagnostic figure per redshift bin.

Run DataProcessing.py and EnvironmentFeatures.py before running this script.
"""

# ----------------------------- LIBRARIES -----------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import functions as fn

# ----------------------------- PATHS -----------------------------

InputPath = "data/BGS_data.parquet"

# Update this path before running
OutputDir = "insert/path/to/scratch/folder/DividerFigures_Z/"

# ----------------------------- LOAD DATA -----------------------------

Data = pd.read_parquet(InputPath)

Z = Data["Z"].to_numpy()
M_r = Data["M_r"].to_numpy()
gr = Data["gr"].to_numpy()

# ----------------------------- BIN DEFINITIONS -----------------------------

zMin, zMax, zStep = 0.0, 0.5, 0.025
zBins = np.arange(zMin, zMax + zStep, zStep)

mrMin, mrMax, mrStep = -14.0, -24.0, -0.2
mrEdges = np.arange(mrMin, mrMax + mrStep, mrStep)

grMin, grMax = -0.5, 1.75
grEdges = np.linspace(grMin, grMax, 101)

# ----------------------------- FIGURE GENERATION -----------------------------

for iz in range(len(zBins) - 1):

    # Select redshift slice
    zLo = zBins[iz]
    zHi = zBins[iz + 1]

    maskZ = (Z > zLo) & (Z <= zHi)

    figZMr = M_r[maskZ]
    figZgr = gr[maskZ]

    # One panel for each absolute-magnitude slice
    fig, axs = plt.subplots(10, 5, figsize=(40, 50), sharex=True, sharey=True)
    axs = axs.ravel()

    for i in range(50):

        # Select one absolute-magnitude slice
        mrHi = mrEdges[i]
        mrLo = mrEdges[i + 1]

        mrMask = (figZMr <= mrHi) & (figZMr > mrLo)
        figData = figZgr[mrMask]

        axs[i].set_title(
            f"${mrHi:.1f} > M_r \\geq {mrLo:.1f}$",
            fontsize=18
        )

        axs[i].grid(True, which="both", alpha=0.4, linewidth=0.5)

        # Skip empty slices
        if len(figData) == 0:
            axs[i].text(
                0.02, 0.95,
                f"N = {len(figData)}",
                transform=axs[i].transAxes,
                ha="left", va="top",
                fontsize=12
            )
            axs[i].text(
                0.5, 0.5,
                "Fit Failed",
                transform=axs[i].transAxes,
                ha="center", va="center",
                fontsize=12,
                color="red"
            )
            continue

        counts, _ = np.histogram(figData, bins=grEdges)
        peak = counts.max()

        if peak == 0:
            continue

        axs[i].hist(
            figData,
            bins=grEdges,
            weights=np.ones_like(figData) / peak,
            alpha=0.2,
            color="gray"
        )

        try:
            # Fit color distribution with a double Gaussian model
            binCenters, counts, dgFit, gauss1, gauss2, popt = fn.fitDoubleGauss(
                figData,
                bins=grEdges
            )

            if np.sum(gauss1) > 0 and np.sum(gauss2) > 0:

                gaussScale = np.max(dgFit)

                gauss1 = gauss1 / gaussScale
                gauss2 = gauss2 / gaussScale
                dgFit = dgFit / gaussScale

                axs[i].plot(binCenters, gauss1, "--", color="blue")
                axs[i].plot(binCenters, gauss2, "--", color="red")
                axs[i].plot(binCenters, dgFit, "-", color="black")

                # Calculate and plot green-valley divider
                xDiv, tau = fn.CR_Div(binCenters, gauss1, gauss2)

                axs[i].axvline(xDiv, color="green")

                axs[i].text(
                    0.02, 0.95,
                    f"N = {len(figData)}\n$\\tau$ = {tau:.2f}",
                    transform=axs[i].transAxes,
                    ha="left", va="top",
                    fontsize=12
                )

        except Exception:
            axs[i].text(
                0.02, 0.95,
                f"N = {len(figData)}",
                transform=axs[i].transAxes,
                ha="left", va="top",
                fontsize=12
            )
            axs[i].text(
                0.5, 0.5,
                "Fit Failed",
                transform=axs[i].transAxes,
                ha="center", va="center",
                fontsize=12,
                color="red"
            )

        # Scale axes
        axs[i].set_xlim(grMin, grMax)
        axs[i].set_ylim(0, 1.25)

        axs[i].set_xticks([-0.5, 0.0, 0.5, 1.0, 1.5])
        axs[i].set_yticks([0.0, 0.5, 1.0])

    # Create shared figure labels
    fig.supxlabel("$g-r$", fontsize=28)
    fig.supylabel("Normalized Counts", fontsize=28)

    fig.tight_layout(rect=[0.02, 0.02, 0.98, 0.98])

    # Save figure
    fileName = f"{OutputDir}Z{zLo:.3f}-{zHi:.3f}.png"

    plt.savefig(fileName, dpi=100)
    plt.close(fig)

    print(f"Saved {fileName}")