"""
Compute local-environment features for the BGS galaxy sample.

This script reads the processed BGS dataset, computes 3D neighbor counts and neighboring stellar-mass density within fixed comoving radii, and saves the catalog with added environment columns.

Run DataProcessing.py before running this script.
"""

# ----------------------------- LIBRARIES -----------------------------

import numpy as np
import pandas as pd
from astropy.cosmology import Planck18
import astropy.units as u
from astropy.coordinates import SkyCoord, search_around_3d

# ----------------------------- PATHS -----------------------------

InputPath = OutputPath = "data/BGS_data.parquet"
# OutputPath = "data/BGS_data.parquet"

# ----------------------------- RADII (Units: Mpc) -----------------------------

SearchRadii = [1, 2, 3]

# ----------------------------- LOAD DATA -----------------------------

def LoadData():
    """
    Load the processed DESI BGS catalog.
    """

    Data = pd.read_parquet(InputPath)

    print(f"Loaded {len(Data):,} galaxies")
    print(f"Loaded {len(Data.columns):,} columns")

    return Data

# ----------------------------- COORDINATES -----------------------------

def BuildSkyCoords(Data):
    """
    Build 3D sky coordinates using RA, Dec, and comoving distance.
    """

    ComovingDist = Planck18.comoving_distance(Data["Z"])

    Coords = SkyCoord(
        ra=np.array(Data["RA"]) * u.deg, 
        dec=np.array(Data["Dec"]) * u.deg, 
        distance=ComovingDist
    )

    return Coords

# ----------------------------- ENVIRONMENT FEATURES -----------------------------

def CalcEnvFeatures(Data, Coords, Radius):
    """
    Calculate neighbor count and neighboring stellar mass within a fixed radius.
    """

    print(f"Calculating environment features within {Radius} Mpc")

    idx1, idx2, sep2d, sep3d = search_around_3d(
        Coords,
        Coords,
        distlimit=Radius * u.Mpc
    )

    # Sort by central-galaxy index
    Sort = np.argsort(idx1)

    idx1 = idx1[Sort]
    idx2 = idx2[Sort]

    # Remove self-matches
    Mask = idx1 != idx2

    idx1 = idx1[Mask]
    idx2 = idx2[Mask]

    # Create arrays to store data
    NeighborCount = np.zeros(len(Data), dtype=int)
    NeighborLogMstar = np.zeros(len(Data))

    UniqueCenters = np.unique(idx1)

    for Center in UniqueCenters:
        NeighborIdx = idx2[idx1 == Center]

        NeighborCount[Center] = len(NeighborIdx)
        NeighborLogMstar[Center] = np.log10(np.sum(10 ** Data["LogMstar"].to_numpy()[NeighborIdx]))

    Data[f"N_N_{Radius}Mpc"] = NeighborCount
    Data[f"LogMstar_N_{Radius}Mpc"] = NeighborLogMstar

    return Data

def AddEnvironmentFeatures(Data):
    """
    Add neighbor-count and neighboring stellar-mass features for all radii.
    """
    
    Coords = BuildSkyCoords(Data)

    for Radius in SearchRadii:
        Data = CalcEnvFeatures(Data, Coords, Radius)

    return Data

# ----------------------------- RUNNING SCRIPT -----------------------------

if __name__ == "__main__":

    # Load processed catalog
    Data = LoadData()

    # Add local-environment features
    Data = AddEnvironmentFeatures(Data)

    # Save updated catalog
    Data.to_parquet(OutputPath, index=False)

    print(f"Saved {len(Data):,} galaxies with environment features to {OutputPath}")