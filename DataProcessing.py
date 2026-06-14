"""
Prepare the DESI BGS data used in the galaxy evolution analysis.

This script extracts galaxy records from DESI FastSpecFit catalogs, applies
basic quality filters, converts cosmology-dependent quantities to the Planck18
convention, derives color/magnitude features, and saves the processed catalog.

Raw DESI catalogs and generated parquet files are not included in this public
repository. 

** Update DesiDir and OutputPath before running. **

"""

# ----------------------------- LIBRARIES -----------------------------

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import Distance
from astropy.cosmology import Planck18
from astropy.io import fits

# ----------------------------- PATHS -----------------------------

# Update path to DESI FastSpecFit catalogs
DesiDir = "/path/to/desi/fastspecfit/catalogs/"

# Update output path for processed catalog
OutputPath = "data/BGS_data.parquet"

# ----------------------------- COLUMNS -----------------------------

MetadataCols = ["TARGETID", "RA", "DEC", "Z"]

SpecphotCols = [
    "TARGETID",
    "LOGMSTAR",
    "ABSMAG01_SDSS_U",
    "ABSMAG01_SDSS_G",
    "ABSMAG01_SDSS_R",
    "ABSMAG01_SDSS_I",
    "ABSMAG01_SDSS_Z",
    "KCORR01_SDSS_U",
    "KCORR01_SDSS_G",
    "KCORR01_SDSS_R",
    "KCORR01_SDSS_I",
    "KCORR01_SDSS_Z",
    "ABSMAG01_IVAR_SDSS_U",
    "ABSMAG01_IVAR_SDSS_G",
    "ABSMAG01_IVAR_SDSS_R",
    "ABSMAG01_IVAR_SDSS_I",
    "ABSMAG01_IVAR_SDSS_Z"
]

RenameCols = {
    "TARGETID": "TargetID",
    "DEC": "Dec",
    "LOGMSTAR": "LogMstar",
    "ABSMAG01_SDSS_U": "M_u",
    "ABSMAG01_SDSS_G": "M_g",
    "ABSMAG01_SDSS_R": "M_r",
    "ABSMAG01_SDSS_I": "M_i",
    "ABSMAG01_SDSS_Z": "M_z",
    "KCORR01_SDSS_U": "kcorr_u",
    "KCORR01_SDSS_G": "kcorr_g",
    "KCORR01_SDSS_R": "kcorr_r",
    "KCORR01_SDSS_I": "kcorr_i",
    "KCORR01_SDSS_Z": "kcorr_z",
    "ABSMAG01_IVAR_SDSS_U": "ivar_M_u",
    "ABSMAG01_IVAR_SDSS_G": "ivar_M_g",
    "ABSMAG01_IVAR_SDSS_R": "ivar_M_r",
    "ABSMAG01_IVAR_SDSS_I": "ivar_M_i",
    "ABSMAG01_IVAR_SDSS_Z": "ivar_M_z"
}

MagCols = ["M_u", "M_g", "M_r", "M_i", "M_z"]

# ----------------------------- READ FILES -----------------------------

def ReadFiles(FilePath):
    """
    Read FastSpecFit file and return galaxy records.
    """
    
    with fits.open(FilePath) as Hdul:
        
        Metadata = Hdul[1].data
        Specphot = Hdul[2].data

        GalaxyMask = Metadata["SPECTYPE"] == "GALAXY"

        MetadataDf = pd.DataFrame({Col: Metadata[Col][GalaxyMask] for Col in MetadataCols})
        SpecphotDf = pd.DataFrame({Col: Specphot[Col] for Col in SpecphotCols})

    Data = MetadataDf.merge(SpecphotDf, on="TARGETID", how="inner")
    Data = Data.rename(columns=RenameCols)

    return Data

def LoadCatalog():
    """
    Load and combine the DESI BGS FastSpecFit files.
    """
    
    DataChunks = []

    for i in range(12):
        FileName = f"fastspec-iron-main-bright-nside1-hp{i:02d}.fits"
        FilePath = DesiDir + FileName

        print(f"Reading {FileName}")
        DataChunks.append(ReadFiles(FilePath))

    Data = pd.concat(DataChunks, ignore_index=True)

    print(f"Loaded {len(Data):,} galaxies before filters")

    return Data

# ----------------------------- CLEANING -----------------------------

def ApplyFilters(Data):
    """
    Apply redshift and absolute-magnitude filters.
    """
    
    Mask = Data["Z"] > 0.004

    for Col in MagCols:
        Mask = Mask & (-30 < Data[Col]) & (Data[Col] < -10)

    Data = Data[Mask].copy()

    print(f"Retained {len(Data):,} galaxies after filtering")

    return Data

# ----------------------------- RESCALING -----------------------------

def RescaleAbsMag(Mag, h=Planck18.h):
    """
    Convert DESI absolute magnitudes from h=1.0 (default) to the Planck18 convention.

    The correction is:
        M_rescaled = M_catalog + 5 log10(h)
    """

    AbsMag = Mag + 5 * np.log10(h)
    
    return AbsMag

def RescaleCosmology(Data):
    """
    Apply h-dependent logarithmic stellar mass and absolute magnitude corrections.
    """
    
    Data["LogMstar"] = Data["LogMstar"] - 2 * np.log10(Planck18.h)

    for Col in MagCols:
        Data[Col] = RescaleAbsMag(Data[Col])

    return Data

# ----------------------------- CREATING FEATURES -----------------------------

def CalcDistMod(Data):
    """
    Calculate luminosity distance in parsecs.
    """
    
    Data["D_pc"] = Distance(z=Data["Z"], unit=u.pc, cosmology=Planck18).value

    return Data

def CalcAppMag(Data):
    """
    Calculate apparent magnitudes from absolute magnitudes and luminosity distance.
    """
    
    Data["m_u"] = Data["M_u"] + 5 * np.log10(Data["D_pc"]) - 5
    Data["m_g"] = Data["M_g"] + 5 * np.log10(Data["D_pc"]) - 5
    Data["m_r"] = Data["M_r"] + 5 * np.log10(Data["D_pc"]) - 5
    Data["m_i"] = Data["M_i"] + 5 * np.log10(Data["D_pc"]) - 5
    Data["m_z"] = Data["M_z"] + 5 * np.log10(Data["D_pc"]) - 5

    return Data

def CalcColor(Data):
    """
    Calculate color features.
    """
    
    Data["gr"] = Data["m_g"] - Data["m_r"]
    Data["ur"] = Data["m_u"] - Data["m_r"]

    return Data

# ----------------------------- RUNNING SCRIPT -----------------------------

if __name__ == "__main__":

    # Call functions in order
    Data = LoadCatalog()
    Data = ApplyFilters(Data)
    Data = RescaleCosmology(Data)
    Data = CalcDistMod(Data)
    Data = CalcAppMag(Data)
    Data = CalcColor(Data)

    # Save data
    Data.to_parquet(OutputPath, index=False)
    print(f"Saved {len(Data):,} processed galaxies to {OutputPath}")