# -------------------------------------------------
# Generate double-Gaussian + divider plots
# One 50-panel figure per redshift bin, saved into
# separate folders for each environment
# -------------------------------------------------

import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import functions as fn

# -----------------------------
# Timing helper
# -----------------------------
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
rho = data["LogMstar_N_2Mpc"].values

# -----------------------------
# Bins (same as cache script)
# -----------------------------
zMin, zMax, zStep = 0.0, 0.5, 0.025
zEdges = np.arange(zMin, zMax + zStep, zStep)
nZ = len(zEdges) - 1

mrMin, mrMax, mrStep = -24.0, -14.0, 0.2
mrEdges   = np.arange(mrMin, mrMax + mrStep, mrStep)
mrCenters = 0.5 * (mrEdges[:-1] + mrEdges[1:])
nMr = len(mrCenters)   # should be 50

grMin, grMax = -0.5, 1.75
grEdges = np.linspace(grMin, grMax, 101)

# -----------------------------
# Output folders
# -----------------------------
base_path = "/pscratch/sd/s/suhas31/AppendixFigures"

env_folder_names = [
    "NoNeighbours",
    "LowDensity",
    "MidDensity",
    "HighDensity",
]

env_titles = [
    "No neighbours",
    "Low Density",
    "Mid Density",
    "High Density",
]

for folder in env_folder_names:
    os.makedirs(os.path.join(base_path, folder), exist_ok=True)

stamp("Starting calculations")

# -----------------------------
# Loop over redshift bins
# -----------------------------
for iz in range(nZ):

    z_lo = zEdges[iz]
    z_hi = zEdges[iz + 1]

    stamp(f"z bin {iz+1}/{nZ}: {z_lo:.3f} < z <= {z_hi:.3f}")

    maskZ = (Z > z_lo) & (Z <= z_hi)

    FigZ_Mr  = M_r[maskZ]
    FigZ_gr  = gr[maskZ]
    FigZ_rho = rho[maskZ]

    # same environment logic as cache builder
    mask0_z = FigZ_rho != 0

    if np.count_nonzero(mask0_z) > 0:
        p33 = np.percentile(FigZ_rho[mask0_z], 33.3333333333)
        p67 = np.percentile(FigZ_rho[mask0_z], 66.6666666667)
    else:
        p33 = np.nan
        p67 = np.nan

    # -----------------------------
    # Loop over environments
    # -----------------------------
    for j in range(4):

        stamp(f"  Environment: {env_titles[j]}")

        if j == 0:
            # no neighbours
            gr_cell = FigZ_gr[~mask0_z]
            Mr_cell = FigZ_Mr[~mask0_z]

        else:
            rho_nz = FigZ_rho[mask0_z]
            gr_nz  = FigZ_gr[mask0_z]
            Mr_nz  = FigZ_Mr[mask0_z]

            if j == 1:
                mEnv = rho_nz <= p33
            elif j == 2:
                mEnv = (p33 < rho_nz) & (rho_nz <= p67)
            else:
                mEnv = rho_nz > p67

            gr_cell = gr_nz[mEnv]
            Mr_cell = Mr_nz[mEnv]

        fig, axs = plt.subplots(
            10, 5,
            figsize=(40, 50),
            sharex=True,
            sharey=True
        )

        axs = axs.ravel()

        for i in range(nMr):

            Mr_c = mrCenters[i]
            Mr_lo = Mr_c - mrStep / 2
            Mr_hi = Mr_c + mrStep / 2

            mMr = (Mr_cell >= Mr_lo) & (Mr_cell < Mr_hi)
            Fig_data = gr_cell[mMr]

            ax = axs[i]

            ax.set_title(
                f"${Mr_lo:.1f} \\leq M_r < {Mr_hi:.1f}$",
                fontsize=32
            )

            ax.grid(True, which="both", alpha=0.5, linewidth=0.5)

            # -----------------------------
            # Empty bin
            # -----------------------------
            if len(Fig_data) == 0:
                ax.text(
                    0.02, 0.95,
                    f"N = {len(Fig_data)}",
                    transform=ax.transAxes,
                    ha="left", va="top",
                    fontsize=32
                )
                ax.text(
                    0.5, 0.5, "Fit Failed",
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

            # -----------------------------
            # Histogram + normalization
            # -----------------------------
            counts, _ = np.histogram(Fig_data, bins=grEdges)
            peak = counts.max()

            if peak == 0:
                ax.text(
                    0.02, 0.95,
                    f"N = {len(Fig_data)}",
                    transform=ax.transAxes,
                    ha="left", va="top",
                    fontsize=32
                )
                ax.text(
                    0.5, 0.5, "Fit Failed",
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
                Fig_data,
                bins=grEdges,
                weights=np.ones_like(Fig_data) / peak,
                alpha=0.3,
                color="gray"
            )

            # -----------------------------
            # Fit double Gaussian + divider
            # -----------------------------
            try:
                bin_centers, counts_fit, dg_fit, gauss1, gauss2, popt = \
                    fn.fit_double_gauss(Fig_data, bins=grEdges)

                if np.sum(gauss1) > 0 and np.sum(gauss2) > 0:

                    # normalize to histogram peak, matching your old script
                    gauss1 = gauss1 / peak
                    gauss2 = gauss2 / peak
                    dg_fit = dg_fit / peak

                    # enforce blue = lower-mean component, red = higher-mean component
                    mean1 = popt[1]
                    mean2 = popt[4]

                    if mean1 <= mean2:
                        blue = gauss1
                        red = gauss2
                    else:
                        blue = gauss2
                        red = gauss1

                    ax.plot(bin_centers, blue, "--", color="blue")
                    ax.plot(bin_centers, red,  "--", color="red")
                    ax.plot(bin_centers, dg_fit, "-", color="black")

                    x_div, tau = fn.CR_div(bin_centers, blue, red)

                    if np.isfinite(x_div):
                        ax.axvline(x_div, color="green")

                    if tau is None or not np.isfinite(tau):
                        ax.text(
                            0.02, 0.95,
                            f"N = {len(Fig_data)}\n$\\tau$ = NaN",
                            transform=ax.transAxes,
                            ha="left", va="top",
                            fontsize=32
                        )
                    else:
                        ax.text(
                            0.02, 0.95,
                            f"N = {len(Fig_data)}\n$\\tau$ = {tau:.2f}",
                            transform=ax.transAxes,
                            ha="left", va="top",
                            fontsize=32
                        )

                else:
                    ax.text(
                        0.02, 0.95,
                        f"N = {len(Fig_data)}",
                        transform=ax.transAxes,
                        ha="left", va="top",
                        fontsize=32
                    )
                    ax.text(
                        0.5, 0.5, "Fit Failed",
                        transform=ax.transAxes,
                        ha="center", va="center",
                        fontsize=28,
                        color="red",
                        weight="bold"
                    )

            except Exception as e:
                print(
                    f"Bin {i} failed | env={env_titles[j]} | "
                    f"z={z_lo:.3f}-{z_hi:.3f} | N={len(Fig_data)} | error={repr(e)}",
                    flush=True
                )

                ax.text(
                    0.02, 0.95,
                    f"N = {len(Fig_data)}",
                    transform=ax.transAxes,
                    ha="left", va="top",
                    fontsize=32
                )
                ax.text(
                    0.5, 0.5, "Fit Failed",
                    transform=ax.transAxes,
                    ha="center", va="center",
                    fontsize=28,
                    color="red",
                    weight="bold"
                )

            # -----------------------------
            # Fixed axis scaling
            # -----------------------------
            ax.set_xlim(grMin, grMax)
            ax.set_ylim(0, 1.25)

            ax.set_xticks([-0.5, 0.0, 0.5, 1.0, 1.5])
            ax.set_yticks([0.0, 0.5, 1.0])

        # -----------------------------
        # Global labels + save
        # -----------------------------
        fig.supxlabel(r"$\mathbf{g-r}$", fontsize=34)
        fig.supylabel("Normalized Counts", fontsize=34, weight="bold")

        fig.tight_layout(rect=[0.02, 0.01, 0.98, 0.98])

        save_dir = os.path.join(base_path, env_folder_names[j])
        filename = os.path.join(
            save_dir,
            f"Z{z_lo:.3f}-{z_hi:.3f}_{env_folder_names[j]}.pdf"
        )

        plt.savefig(filename)
        plt.close(fig)

        stamp(f"    Saved {filename}")

stamp("Finished")