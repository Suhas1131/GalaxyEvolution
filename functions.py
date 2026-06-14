"""
Helper functions for fitting galaxy color distributions and selecting
green-valley dividers.
"""

# ----------------------------- LIBRARIES -----------------------------

import numpy as np
from scipy.optimize import curve_fit

# ----------------------------- FUNCTIONS -----------------------------

def gaussian(x, A, mean, sigma):
    """
    Calculates a single Gaussian curve.

    Inputs:
        x (array-like): Data values.
        A (float): Gaussian amplitude.
        mean (float): Gaussian mean.
        sigma (float): Gaussian standard deviation.

    Returns:
        gauss (array-like): Gaussian curve.
    """

    gauss = A * np.exp(-0.5 * ((x - mean) / sigma) ** 2)

    return gauss

def doubleGaussian(x, A1, mean1, sigma1, A2, mean2, sigma2):
    """
    Calculates a double Gaussian curve.

    Inputs:
        x (array-like): Data values.
        A1 (float): Amplitude of the first Gaussian.
        mean1 (float): Mean of the first Gaussian.
        sigma1 (float): Standard deviation of the first Gaussian.
        A2 (float): Amplitude of the second Gaussian.
        mean2 (float): Mean of the second Gaussian.
        sigma2 (float): Standard deviation of the second Gaussian.

    Returns:
        doubleGauss (array-like): Sum of the two Gaussian components.
        gauss1 (array-like): Bluer Gaussian component.
        gauss2 (array-like): Redder Gaussian component.
    """
    
    if mean1 <= mean2:
        gauss1 = gaussian(x, A1, mean1, sigma1)
        gauss2 = gaussian(x, A2, mean2, sigma2)
    else:
        gauss2 = gaussian(x, A1, mean1, sigma1)
        gauss1 = gaussian(x, A2, mean2, sigma2)

    doubleGauss = gauss1 + gauss2

    return doubleGauss, gauss1, gauss2

def fitDoubleGauss(x, bins=50, weights=None):
    """
    Fits a histogram of galaxy colors with a double Gaussian model.

    Inputs:
        x (array-like): One-dimensional color data to histogram and fit.
        bins (int or array-like): Number of histogram bins or bin edges.
        weights (array-like): Optional histogram weights.

    Returns:
        centers (ndarray): Histogram bin centers.
        counts (ndarray): Histogram counts.
        doubleGaussianFit (ndarray): Best-fit double Gaussian curve.
        gauss1 (ndarray): Bluer Gaussian component.
        gauss2 (ndarray): Redder Gaussian component.
        popt (ndarray): Best-fit parameters.
    """

    # Build histogram
    counts, edges = np.histogram(x, bins=bins, density=False, weights=weights)
    centers = (edges[:-1] + edges[1:]) / 2

    # Model to pass through curve_fit
    def sumGauss(x, A1, mean1, sigma1, A2, mean2, sigma2):
        doubleGauss, _, _ = doubleGaussian(x, A1, mean1, sigma1, A2, mean2, sigma2)

        return doubleGauss

    # Initial guesses
    p0 = [
        counts.max(), np.percentile(x, 25), np.std(x) * 0.3,
        counts.max() / 2, np.percentile(x, 75), np.std(x) * 0.3
    ]

    lower = [0, 0, 1e-6, 0, 0, 1e-6]
    upper = [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]

    # Fit double Gaussian model
    popt, pcov = curve_fit(
        sumGauss,
        centers,
        counts,
        p0=p0,
        bounds=(lower, upper),
        maxfev=100000
    )

    # Calculate fitted curves
    doubleGaussianFit, gauss1, gauss2 = doubleGaussian(centers, *popt)

    return centers, counts, doubleGaussianFit, gauss1, gauss2, popt

def CR_Div(bin_centers, gauss1, gauss2, eps=1e-12):
    """
    Selects the green-valley divider by maximizing a completeness-reliability score.

    Inputs:
        bin_centers (array-like): Color-bin centers.
        gauss1 (array-like): Bluer Gaussian component.
        gauss2 (array-like): Redder Gaussian component.
        eps (float): Small value used to avoid division by zero.

    Returns:
        bestX (float): Selected g-r divider.
        bestTau (float): Maximum completeness-reliability score.
    """

    g1Tot = np.sum(gauss1)
    g2Tot = np.sum(gauss2)

    # Return nan if tau is undefined
    if (
        (not np.isfinite(g1Tot)) or
        (not np.isfinite(g2Tot)) or
        (g1Tot <= eps) or
        (g2Tot <= eps)
    ):
        return np.nan, np.nan

    bestTau = -np.inf
    bestX = np.nan

    for k in bin_centers[1:]:
        mask1 = bin_centers < k
        mask2 = ~mask1

        g1B = np.sum(gauss1[mask1])
        g2B = np.sum(gauss2[mask1])
        g1R = np.sum(gauss1[mask2])
        g2R = np.sum(gauss2[mask2])

        denomB = g1B + g2B
        denomR = g1R + g2R

        # Skip dividers where reliability is undefined
        if denomB <= eps or denomR <= eps:
            continue

        CB = g1B / g1Tot
        CR = g2R / g2Tot

        RB = g1B / denomB
        RR = g2R / denomR

        tauNew = CB * CR * RB * RR  # Always <=1 

        if np.isfinite(tauNew) and tauNew > bestTau:
            bestTau = tauNew
            bestX = k

    if not np.isfinite(bestTau):
        return np.nan, np.nan

    return bestX, bestTau