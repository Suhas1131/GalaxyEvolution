from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
import numpy as np
import warnings



def gaussian(x, A, mean, sigma):
    """
    Calculates the gaussian.

    Parameters:
        - x (array-like): Array of data points to be operated on.
        - A (float): Amplitude of peak the gaussian.
        - mean (float): Mean of data for the gaussian.
        - sigma (float): Standard deviation of data for the gaussian.

    Returns:
        - gauss (array-like): Gaussian curve
    
    Last Updated: 22 February 2025
    """

    # Calculating the gaussian curve
    gauss = A * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
    return gauss



def double_gaussian(x, A1, mean1, sigma1, A2, mean2, sigma2):
    """
    Calculates the double gaussian curve.

    Parameters:
        - x (array-like): Array of data points to be operated on.
        - A1 (float): Amplitude of the peak of the first gaussian.
        - mean1 (float): Mean of the data for the first gaussian.
        - sigma1 (float): Standard deviation of the data for the first gaussian.
        - A2 (float): Amplitude of the peak of the second gaussian.
        - mean2 (float): Mean of the data for of the second gaussian.
        - sigma2 (float): Standard deviation of the data for the second gaussian.
        
    Returns:
        - gauss1+gauss2 (array-like): Double gaussian curve.
        - gauss1 (array-like): First gaussian curve.
        - gauss2 (array-like): Second gaussian curve.
    
    Last Updated: 28 May 2025
    """
    # Adding 2 gaussians
    if mean1 <= mean2:
        gauss1 = gaussian(x, A1, mean1, sigma1)
        gauss2 = gaussian(x, A2, mean2, sigma2)
    else:
        gauss2 = gaussian(x, A1, mean1, sigma1)
        gauss1 = gaussian(x, A2, mean2, sigma2)
    return gauss1+gauss2, gauss1, gauss2


    
def fit_double_gauss(x, bins=50, weights=None):
    """
    Fit a histogram of x with a double Gaussian and return the curves.

    Parameters:
        - x (array-like): 1D data to be histogrammed and fitted.
        - bins (int, optional): Number of bins (or bin edges) for the histogram. Default is 50.

    Returns:
        - centers (ndarray): The bin-center values at which the Gaussians are evaluated.
        - counts (ndarray): The histogram counts for each bin.
        - double_gaussian_fit (ndarray): The best-fit sum of the two Gaussians.
        - gauss1 (ndarray): The first Gaussian component.
        - gauss2 (ndarray): The second Gaussian component.
        - popt (ndarray): Optimal fit parameters [A1, mean1, sigma1, A2, mean2, sigma2].

    Last Updated: 21 Oct 2025
    """
    # Getting histogram info
    counts, edges = np.histogram(x, bins=bins, density=False, weights=weights)
    centers = (edges[:-1] + edges[1:]) / 2

    def sum_gauss(x, A1, m1, s1, A2, m2, s2):
        gsum, _, _ = double_gaussian(x, A1, m1, s1, A2, m2, s2)
        return gsum
    
    # Initial guesses
    p0 = [
        counts.max(), np.percentile(x,25), np.std(x)*0.3,
        counts.max()/2, np.percentile(x,75), np.std(x)*0.3
    ]
    lower = [0, 0, 1e-6,  0, 0, 1e-6]
    upper = [np.inf, np.inf, np.inf,  np.inf, np.inf, np.inf]

    # Adjusting parameters to fit the data [scipy.optimize.curve_fit]
    popt, pcov = curve_fit(
        sum_gauss, centers, counts,
        p0=p0, bounds=(lower, upper),
        maxfev=100000
    )

    # Find indic=vidual gaussians and double gaussian
    double_gaussian_fit, gauss1, gauss2 = double_gaussian(centers, *popt)

    return centers, counts, double_gaussian_fit, gauss1, gauss2, popt


# Version 1 (Outdated: Feb 21, 2026)
# def CR_div(bin_centers, gauss1, gauss2):
#     """
#     Finds the best divider based on completeness and reliability.

#     Parameters:
#         - bin_centers (ndarray):
#         - gauss1 (ndarray):
#         - gauss2 (ndarray):

#     Returns:
#         - x (float): g-r divider
#         - best_Tau (float): value of best Tau
#     """
#     best_Tau, x = -np.inf, None
#     for k in bin_centers[1:]:
#         mask1 = bin_centers < k
#         mask2 = bin_centers >= k

#         denomB = np.sum(gauss1[mask1]) + np.sum(gauss2[mask1])
#         denomR = np.sum(gauss1[mask2]) + np.sum(gauss2[mask2])        
#         if denomB == 0 or denomR == 0:
#             return np.nan, np.nan

#         CB = np.sum(gauss1[mask1])/np.sum(gauss1)
#         CR = np.sum(gauss2[mask2])/np.sum(gauss2)
            
#         RB = np.sum(gauss1[mask1])/denomB
#         RR = np.sum(gauss2[mask2])/denomR
        
#         Tau_new = CB * CR * RB * RR
#         if Tau_new > best_Tau:
#             x = k
#             best_Tau = Tau_new
#     return x, best_Tau

def CR_div(bin_centers, gauss1, gauss2, eps=1e-12):
    """
    Finds the best divider based on completeness and reliability.

    Returns:
        x (float): g-r divider
        best_Tau (float): value of best Tau

    Last updated: Feb 21 2026
    """
    g1_tot = np.sum(gauss1)
    g2_tot = np.sum(gauss2)

    # If a component is effectively absent, tau isn't meaningful
    if (not np.isfinite(g1_tot)) or (not np.isfinite(g2_tot)) or (g1_tot <= eps) or (g2_tot <= eps):
        return np.nan, np.nan

    best_Tau = -np.inf
    best_x = np.nan

    for k in bin_centers[1:]:
        mask1 = bin_centers < k
        mask2 = ~mask1  # bin_centers >= k

        g1_B = np.sum(gauss1[mask1])
        g2_B = np.sum(gauss2[mask1])
        g1_R = np.sum(gauss1[mask2])
        g2_R = np.sum(gauss2[mask2])

        denomB = g1_B + g2_B
        denomR = g1_R + g2_R

        # Skip dividers where reliability is undefined (0/0 or tiny/tiny)
        if denomB <= eps or denomR <= eps:
            continue

        CB = g1_B / g1_tot
        CR = g2_R / g2_tot

        RB = g1_B / denomB
        RR = g2_R / denomR

        Tau_new = CB * CR * RB * RR

        if np.isfinite(Tau_new) and Tau_new > best_Tau:
            best_Tau = Tau_new
            best_x = k

    if not np.isfinite(best_Tau):
        return np.nan, np.nan

    return best_x, best_Tau