"""
================================================================================
All-Sky Downwelling Longwave Radiation (DLR) Estimation Algorithm
================================================================================

This module implements the all-sky DLR estimation framework described in:

    "Estimation of Longwave Sky Emissivity under All-Sky Conditions Using 
    GHI, Temperature, and Relative Humidity Measurements"
    
    Wang, J., Srivastava, A., Huang, Z., Briggs, C., & Coimbra, C.F.M.
    University of California San Diego
    Submitted to Solar Energy, 2026

Usage:
------
    from dlr_estimation import estimate_dlr, DLREstimator
    
    # Quick calculation
    dlr = estimate_dlr(T_air=298.15, RH=0.60, altitude=213, GHI=650, DHI=195)
    
    # Or use the class for batch processing
    estimator = DLREstimator()
    dlr = estimator.calculate(T_air, RH, altitude, GHI, DHI)

Author: Coimbra Research Group, UC San Diego
License: MIT
================================================================================
"""

import numpy as np
from typing import Union, Optional, Tuple
from dataclasses import dataclass

# Physical Constants
STEFAN_BOLTZMANN = 5.67e-8  # W/(m²·K⁴)
P0 = 101325.0  # Standard atmospheric pressure [Pa]
SCALE_HEIGHT = 8500.0  # Atmospheric scale height [m]


@dataclass
class DLRResult:
    """Container for DLR estimation results with intermediate values."""
    dlr: float  # Downwelling Longwave Radiation [W/m²]
    emissivity_all_sky: float  # All-sky emissivity [-]
    emissivity_clear: float  # Clear-sky emissivity [-]
    gamma: float  # Cloud fraction factor [-]
    k_d: float  # Diffuse fraction [-]
    k_d_source: str  # 'measured' or 'estimated'
    p_w: float  # Dimensionless water vapor pressure [-]


def saturation_vapor_pressure(T_air: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Calculate saturation vapor pressure using Alduchov & Eskridge (1996) formula.
    
    Parameters
    ----------
    T_air : float or array-like
        Air temperature at screening height [K]
        
    Returns
    -------
    P_s : float or array-like
        Saturation vapor pressure [Pa]
        
    Reference
    ---------
    Alduchov, O.A., Eskridge, R.E., 1996. Improved Magnus form approximation 
    of saturation vapor pressure. J. Appl. Meteorol. 35, 601-609.
    """
    T_celsius = T_air - 273.15
    return 610.94 * np.exp(17.625 * T_celsius / (T_air - 30.11))


def clear_sky_emissivity(p_w: Union[float, np.ndarray], 
                         altitude: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Calculate clear-sky emissivity using Matsunobu et al. (2024) formula.
    
    Parameters
    ----------
    p_w : float or array-like
        Dimensionless partial pressure of water vapor [-]
    altitude : float or array-like
        Site altitude [m]
        
    Returns
    -------
    epsilon_clear : float or array-like
        Clear-sky emissivity [-]
        
    Reference
    ---------
    Matsunobu, T., et al., 2024. Effective sky emissivity for longwave radiation.
    """
    altitude_correction = 0.150 * (np.exp(-altitude / SCALE_HEIGHT) - 1)
    return 0.600 + 1.652 * np.sqrt(p_w) + altitude_correction


def erbs_model(k_t: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Estimate diffuse fraction k_d from clearness index k_t using Erbs et al. (1982).
    
    Parameters
    ----------
    k_t : float or array-like
        Clearness index (GHI_measured / GHI_clear-sky) [-]
        
    Returns
    -------
    k_d : float or array-like
        Estimated diffuse fraction [-]
        
    Reference
    ---------
    Erbs, D.G., Klein, S.A., Duffie, J.A., 1982. Estimation of the diffuse 
    radiation fraction for hourly, daily and monthly-average global radiation.
    Solar Energy 28(4), 293-302.
    """
    k_t = np.asarray(k_t)
    k_d = np.zeros_like(k_t, dtype=float)
    
    # Region 1: k_t <= 0.22
    mask1 = k_t <= 0.22
    k_d[mask1] = 1.0 - 0.09 * k_t[mask1]
    
    # Region 2: 0.22 < k_t <= 0.80
    mask2 = (k_t > 0.22) & (k_t <= 0.80)
    kt2 = k_t[mask2]
    k_d[mask2] = (0.9511 - 0.1604 * kt2 + 4.388 * kt2**2 
                  - 16.638 * kt2**3 + 12.336 * kt2**4)
    
    # Region 3: k_t > 0.80
    mask3 = k_t > 0.80
    k_d[mask3] = 0.165
    
    return float(k_d) if k_d.ndim == 0 else k_d


def orgill_hollands_model(k_t: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Estimate diffuse fraction k_d from clearness index k_t using Orgill-Hollands (1977).
    
    Parameters
    ----------
    k_t : float or array-like
        Clearness index [-]
        
    Returns
    -------
    k_d : float or array-like
        Estimated diffuse fraction [-]
        
    Reference
    ---------
    Orgill, J.F., Hollands, K.G.T., 1977. Correlation equation for hourly 
    diffuse radiation on a horizontal surface. Solar Energy 19(4), 357-359.
    """
    k_t = np.asarray(k_t)
    k_d = np.zeros_like(k_t, dtype=float)
    
    mask1 = k_t <= 0.35
    k_d[mask1] = 1.0 - 0.249 * k_t[mask1]
    
    mask2 = (k_t > 0.35) & (k_t <= 0.75)
    k_d[mask2] = 1.557 - 1.84 * k_t[mask2]
    
    mask3 = k_t > 0.75
    k_d[mask3] = 0.177
    
    return float(k_d) if k_d.ndim == 0 else k_d


def reindl_model(k_t: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Estimate diffuse fraction k_d from clearness index k_t using Reindl et al. (1990).
    
    Parameters
    ----------
    k_t : float or array-like
        Clearness index [-]
        
    Returns
    -------
    k_d : float or array-like
        Estimated diffuse fraction [-]
        
    Reference
    ---------
    Reindl, D.T., Beckman, W.A., Duffie, J.A., 1990. Diffuse fraction correlations.
    Solar Energy 45(1), 1-7.
    """
    k_t = np.asarray(k_t)
    k_d = np.zeros_like(k_t, dtype=float)
    
    mask1 = k_t <= 0.30
    k_d[mask1] = 1.45 - 1.67 * k_t[mask1]
    
    mask2 = (k_t > 0.30) & (k_t <= 0.78)
    k_d[mask2] = 1.02 - 1.02 * k_t[mask2]
    
    mask3 = k_t > 0.78
    k_d[mask3] = 0.147
    
    return float(k_d) if k_d.ndim == 0 else k_d


def cloud_fraction_factor(k_d: Union[float, np.ndarray], 
                          c1: float = 0.585, 
                          c2: float = 1.748) -> Union[float, np.ndarray]:
    """
    Calculate cloud fraction factor gamma from diffuse fraction k_d.
    
    This is the core empirical relationship proposed in the paper.
    
    Parameters
    ----------
    k_d : float or array-like
        Diffuse fraction (DHI/GHI) [-]
    c1 : float, optional
        First fitting parameter (default: 0.585)
    c2 : float, optional
        Second fitting parameter (default: 1.748)
        
    Returns
    -------
    gamma : float or array-like
        Cloud fraction factor [-]
        
    Notes
    -----
    The default parameters c1=0.585 and c2=1.748 are optimized using 
    2010-2022 SURFRAD data from seven stations across diverse climates.
    
    Expected performance (with measured k_d):
        - rRMSE: ~4.6%
        - R²: ~0.83
    """
    return c1 * np.power(k_d, c2)


def all_sky_emissivity(epsilon_clear: Union[float, np.ndarray], 
                       gamma: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Calculate all-sky emissivity from clear-sky emissivity and cloud factor.
    
    Parameters
    ----------
    epsilon_clear : float or array-like
        Clear-sky emissivity [-]
    gamma : float or array-like
        Cloud fraction factor [-]
        
    Returns
    -------
    epsilon_all : float or array-like
        All-sky emissivity [-]
        
    Reference
    ---------
    Coimbra, C.F.M., 2025. Effective total emissivity model for cloudy sky.
    """
    return epsilon_clear + gamma * (1 - epsilon_clear)


def estimate_dlr(T_air: Union[float, np.ndarray],
                 RH: Union[float, np.ndarray],
                 altitude: Union[float, np.ndarray],
                 GHI: Optional[Union[float, np.ndarray]] = None,
                 DHI: Optional[Union[float, np.ndarray]] = None,
                 GHI_clearsky: Optional[Union[float, np.ndarray]] = None,
                 k_d_model: str = 'erbs',
                 return_details: bool = False) -> Union[float, np.ndarray, DLRResult]:
    """
    Estimate Downwelling Longwave Radiation (DLR) under all-sky conditions.
    
    This is the main function implementing the complete DLR estimation algorithm.
    
    Parameters
    ----------
    T_air : float or array-like
        Air temperature at screening height [K]
    RH : float or array-like
        Relative humidity [0-1]
    altitude : float or array-like
        Site altitude above sea level [m]
    GHI : float or array-like, optional
        Global Horizontal Irradiance [W/m²]. Required if DHI is provided.
    DHI : float or array-like, optional
        Diffuse Horizontal Irradiance [W/m²]. If provided, k_d is calculated directly.
    GHI_clearsky : float or array-like, optional
        Clear-sky GHI [W/m²]. Required if DHI is not provided but k_d estimation
        is needed. Can be obtained from pvlib or similar libraries.
    k_d_model : str, optional
        Model for estimating k_d from k_t when DHI is not available.
        Options: 'erbs' (default), 'orgill_hollands', 'reindl'
    return_details : bool, optional
        If True, return DLRResult object with intermediate calculations.
        Default is False (returns only DLR value).
        
    Returns
    -------
    DLR : float, array-like, or DLRResult
        Downwelling Longwave Radiation [W/m²]. If return_details=True, 
        returns DLRResult object containing all intermediate values.
        
    Raises
    ------
    ValueError
        If insufficient inputs are provided for k_d calculation.
        
    Examples
    --------
    # Scenario A: Full radiation data available (GHI + DHI)
    >>> dlr = estimate_dlr(T_air=298.15, RH=0.60, altitude=213, GHI=650, DHI=195)
    >>> print(f"DLR = {dlr:.1f} W/m²")
    DLR = 374.4 W/m²
    
    # Scenario B: Only GHI available (need clear-sky model)
    >>> dlr = estimate_dlr(T_air=298.15, RH=0.60, altitude=213, 
    ...                    GHI=650, GHI_clearsky=900, k_d_model='erbs')
    
    # Scenario C: Meteorological data only (clear-sky estimate)
    >>> dlr = estimate_dlr(T_air=298.15, RH=0.60, altitude=213)
    >>> print(f"Clear-sky DLR = {dlr:.1f} W/m²")
    
    # Get detailed results
    >>> result = estimate_dlr(T_air=298.15, RH=0.60, altitude=213, 
    ...                       GHI=650, DHI=195, return_details=True)
    >>> print(f"All-sky emissivity = {result.emissivity_all_sky:.4f}")
    """
    
    # Step 1-2: Compute water vapor pressure
    P_s = saturation_vapor_pressure(T_air)
    P_w = P_s * RH
    p_w = P_w / P0
    
    # Step 3: Clear-sky emissivity
    eps_clear = clear_sky_emissivity(p_w, altitude)
    
    # Step 4: Determine diffuse fraction
    if DHI is not None and GHI is not None:
        # Direct calculation from measurements
        k_d = DHI / GHI
        k_d_source = 'measured'
    elif GHI is not None and GHI_clearsky is not None:
        # Estimate from clearness index
        k_t = GHI / GHI_clearsky
        k_t = np.clip(k_t, 0, 1.5)  # Physical bounds
        
        if k_d_model == 'erbs':
            k_d = erbs_model(k_t)
        elif k_d_model == 'orgill_hollands':
            k_d = orgill_hollands_model(k_t)
        elif k_d_model == 'reindl':
            k_d = reindl_model(k_t)
        else:
            raise ValueError(f"Unknown k_d model: {k_d_model}")
        k_d_source = f'estimated ({k_d_model})'
    else:
        # Clear-sky only (no cloud correction)
        k_d = 0.0
        k_d_source = 'clear-sky (no clouds)'
    
    # Step 5: Cloud fraction factor
    gamma = cloud_fraction_factor(k_d)
    
    # Step 6: All-sky emissivity
    eps_all = all_sky_emissivity(eps_clear, gamma)
    
    # Step 7: DLR
    dlr = eps_all * STEFAN_BOLTZMANN * np.power(T_air, 4)
    
    if return_details:
        return DLRResult(
            dlr=float(dlr) if np.ndim(dlr) == 0 else dlr,
            emissivity_all_sky=float(eps_all) if np.ndim(eps_all) == 0 else eps_all,
            emissivity_clear=float(eps_clear) if np.ndim(eps_clear) == 0 else eps_clear,
            gamma=float(gamma) if np.ndim(gamma) == 0 else gamma,
            k_d=float(k_d) if np.ndim(k_d) == 0 else k_d,
            k_d_source=k_d_source,
            p_w=float(p_w) if np.ndim(p_w) == 0 else p_w
        )
    
    return dlr


class DLREstimator:
    """
    Class-based interface for DLR estimation with configurable parameters.
    
    Useful for batch processing or when integrating with larger pipelines.
    
    Parameters
    ----------
    c1 : float, optional
        Cloud fraction factor coefficient (default: 0.585)
    c2 : float, optional
        Cloud fraction factor exponent (default: 1.748)
    k_d_model : str, optional
        Default model for k_d estimation (default: 'erbs')
        
    Attributes
    ----------
    c1, c2 : float
        Model parameters for gamma = c1 * k_d^c2
    k_d_model : str
        Default diffuse fraction estimation model
        
    Examples
    --------
    >>> estimator = DLREstimator()
    >>> dlr = estimator.calculate(T_air=298.15, RH=0.60, altitude=213, GHI=650, DHI=195)
    
    >>> # Use custom parameters (e.g., station-specific fit)
    >>> estimator_bon = DLREstimator(c1=0.607, c2=1.744)  # Bondville fit
    """
    
    def __init__(self, c1: float = 0.585, c2: float = 1.748, k_d_model: str = 'erbs'):
        self.c1 = c1
        self.c2 = c2
        self.k_d_model = k_d_model
        
    def calculate(self, 
                  T_air: Union[float, np.ndarray],
                  RH: Union[float, np.ndarray],
                  altitude: Union[float, np.ndarray],
                  GHI: Optional[Union[float, np.ndarray]] = None,
                  DHI: Optional[Union[float, np.ndarray]] = None,
                  GHI_clearsky: Optional[Union[float, np.ndarray]] = None,
                  return_details: bool = False) -> Union[float, np.ndarray, DLRResult]:
        """
        Calculate DLR using the estimator's configured parameters.
        
        See estimate_dlr() for parameter documentation.
        """
        # Use custom gamma calculation with instance parameters
        P_s = saturation_vapor_pressure(T_air)
        P_w = P_s * RH
        p_w = P_w / P0
        
        eps_clear = clear_sky_emissivity(p_w, altitude)
        
        if DHI is not None and GHI is not None:
            k_d = DHI / GHI
            k_d_source = 'measured'
        elif GHI is not None and GHI_clearsky is not None:
            k_t = np.clip(GHI / GHI_clearsky, 0, 1.5)
            if self.k_d_model == 'erbs':
                k_d = erbs_model(k_t)
            elif self.k_d_model == 'orgill_hollands':
                k_d = orgill_hollands_model(k_t)
            elif self.k_d_model == 'reindl':
                k_d = reindl_model(k_t)
            else:
                raise ValueError(f"Unknown k_d model: {self.k_d_model}")
            k_d_source = f'estimated ({self.k_d_model})'
        else:
            k_d = 0.0
            k_d_source = 'clear-sky'
        
        # Use instance parameters for gamma
        gamma = self.c1 * np.power(k_d, self.c2)
        eps_all = all_sky_emissivity(eps_clear, gamma)
        dlr = eps_all * STEFAN_BOLTZMANN * np.power(T_air, 4)
        
        if return_details:
            return DLRResult(
                dlr=float(dlr) if np.ndim(dlr) == 0 else dlr,
                emissivity_all_sky=float(eps_all) if np.ndim(eps_all) == 0 else eps_all,
                emissivity_clear=float(eps_clear) if np.ndim(eps_clear) == 0 else eps_clear,
                gamma=float(gamma) if np.ndim(gamma) == 0 else gamma,
                k_d=float(k_d) if np.ndim(k_d) == 0 else k_d,
                k_d_source=k_d_source,
                p_w=float(p_w) if np.ndim(p_w) == 0 else p_w
            )
        return dlr


# Station-specific parameters (optional fine-tuning)
STATION_PARAMS = {
    'BON': {'c1': 0.607, 'c2': 1.744, 'altitude': 213},   # Bondville, IL
    'DRA': {'c1': 0.323, 'c2': 0.977, 'altitude': 1007},  # Desert Rock, NV
    'FPK': {'c1': 0.580, 'c2': 1.505, 'altitude': 634},   # Fort Peck, MT
    'GWN': {'c1': 0.580, 'c2': 1.505, 'altitude': 98},    # Goodwin Creek, MS
    'PSU': {'c1': 0.625, 'c2': 1.651, 'altitude': 376},   # Penn State, PA
    'SXF': {'c1': 0.623, 'c2': 1.591, 'altitude': 473},   # Sioux Falls, SD
    'TBL': {'c1': 0.500, 'c2': 1.508, 'altitude': 1689},  # Table Mountain, CO
}


def create_station_estimator(station: str) -> DLREstimator:
    """
    Create a DLREstimator with station-specific parameters.
    
    Parameters
    ----------
    station : str
        SURFRAD station code: 'BON', 'DRA', 'FPK', 'GWN', 'PSU', 'SXF', or 'TBL'
        
    Returns
    -------
    DLREstimator
        Estimator configured with station-specific c1 and c2 parameters
        
    Notes
    -----
    The global model (c1=0.585, c2=1.748) is recommended for general use.
    Station-specific parameters may provide marginal improvements for sites
    with similar climatic conditions to the reference station.
    """
    if station not in STATION_PARAMS:
        raise ValueError(f"Unknown station: {station}. Available: {list(STATION_PARAMS.keys())}")
    
    params = STATION_PARAMS[station]
    return DLREstimator(c1=params['c1'], c2=params['c2'])


# =============================================================================
# VERIFICATION: Worked Example from Paper
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("DLR Estimation Algorithm - Verification")
    print("=" * 70)
    
    # Test case: Bondville, IL - typical summer conditions
    T_air = 298.15  # K (25°C)
    RH = 0.60       # 60% relative humidity
    altitude = 213  # m
    GHI = 650       # W/m²
    DHI = 195       # W/m²
    
    print(f"\nInput Parameters:")
    print(f"  T_air    = {T_air} K ({T_air - 273.15}°C)")
    print(f"  RH       = {RH * 100}%")
    print(f"  Altitude = {altitude} m")
    print(f"  GHI      = {GHI} W/m²")
    print(f"  DHI      = {DHI} W/m²")
    
    # Calculate with details
    result = estimate_dlr(T_air, RH, altitude, GHI, DHI, return_details=True)
    
    print(f"\nIntermediate Results:")
    print(f"  p_w (dim. water vapor)   = {result.p_w:.5f}")
    print(f"  ε_clear (clear-sky)      = {result.emissivity_clear:.4f}")
    print(f"  k_d (diffuse fraction)   = {result.k_d:.4f} [{result.k_d_source}]")
    print(f"  γ (cloud factor)         = {result.gamma:.4f}")
    print(f"  ε_all-sky                = {result.emissivity_all_sky:.4f}")
    
    print(f"\nFinal Result:")
    print(f"  DLR = {result.dlr:.1f} W/m²")
    
    # Expected value for verification
    print(f"\n  Expected (from worked example): ~374.4 W/m²")
    print(f"  Difference: {abs(result.dlr - 374.4):.2f} W/m²")
    
    print("\n" + "=" * 70)
    print("Test: Clear-sky only (no radiation data)")
    print("=" * 70)
    
    result_clear = estimate_dlr(T_air, RH, altitude, return_details=True)
    print(f"  Clear-sky DLR = {result_clear.dlr:.1f} W/m²")
    print(f"  ε_clear       = {result_clear.emissivity_clear:.4f}")
    
    print("\n" + "=" * 70)
    print("Test: Estimated k_d from k_t (Erbs model)")
    print("=" * 70)
    
    GHI_clearsky = 900  # Hypothetical clear-sky GHI
    result_erbs = estimate_dlr(T_air, RH, altitude, GHI=GHI, GHI_clearsky=GHI_clearsky,
                                k_d_model='erbs', return_details=True)
    print(f"  k_t (clearness index)    = {GHI/GHI_clearsky:.3f}")
    print(f"  k_d (estimated, Erbs)    = {result_erbs.k_d:.4f}")
    print(f"  DLR (estimated k_d)      = {result_erbs.dlr:.1f} W/m²")
    
    print("\n" + "=" * 70)
