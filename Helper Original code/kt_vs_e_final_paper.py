#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script Name: kt_vs_e_three_models.py
Description:
  - Reads minute-level SURFRAD HDF5 data;
  - Applies minute-level filters;
  - Computes k_t, k_d_obs, p_w, e_sky_obs;
  - Resamples to hourly averages;
  - Defines three kd models;
  - Computes theoretical e_sky surface over 0 ≤ k_t ≤ 1.3;
  - Generates 3D and 2D projection plots per model, saving nine PNGs.
Author: Jiedong Wang
Date: 2025-05-23
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Use built-in MathText with Times New Roman-like serif (STIX)
plt.rcParams.update({
    'text.usetex': False,
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'mathtext.fontset': 'stix',
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'axes.grid': False
})

# Station altitudes (m)
station_altitudes = {
    'BON': 213, 'DRA': 1007, 'FPK': 634,
    'GWN': 98, 'PSU': 376, 'SXF': 473,
    'TBL': 1689
}


# kd model definitions
def kd_orgill_1977(kt):
    return np.where(kt <= 0.35, 1.0 - 0.249 * kt,
                    np.where(kt <= 0.75, 1.557 - 1.84 * kt, 0.177))


def kd_erbs_1982(kt):
    return np.where(kt <= 0.22, 1.0 - 0.09 * kt,
                    np.where(kt <= 0.80, 0.9511 - 0.1604 * kt + 4.388 * kt ** 2
                             - 16.638 * kt ** 3 + 12.336 * kt ** 4, 0.165))


def kd_reindl_1990a(kt):
    return np.where(kt <= 0.30, 1.020 - 0.248 * kt,
                    np.where(kt < 0.78, 1.45 - 1.67 * kt, 0.147))


models = {
    'Orgill–Hollands (1977)': kd_orgill_1977,
    'Erbs et al. (1982)': kd_erbs_1982,
    'Reindl et al. (1990a)': kd_reindl_1990a
}

# Configuration
h5_file = r'D:\Graduate study\Coimbra_research\output_test.h5'
fig_folder = r'D:\Graduate study\Coimbra_research\figures\5_30'
stations = ['BON', 'DRA', 'FPK', 'GWN', 'PSU', 'SXF', 'TBL']
os.makedirs(fig_folder, exist_ok=True)
G_sc, H = 1353, 8434


# Load, filter, compute hourly means
def load_hourly():
    all_hr = []
    with pd.HDFStore(h5_file, 'r') as store:
        for st in stations:
            if st == 'DRA': continue
            df = store[st].copy()
            df.index = pd.to_datetime(df.index)
            # assign altitude
            alt = station_altitudes.get(st, 0)
            # apply minute filters
            m1 = (df['sza'] < 72.5) & (df['ghi_m'] > 0) & (df['dhi_m'] > 0) & (df['ghi_m'] / df['ghi_c'] < 1.5)
            m2 = df['ghi_m'] < (1.2 * G_sc * np.cos(np.radians(df['sza'])) ** 1.2 + 50)
            m3 = df['dni_m'] < (0.95 * G_sc * np.cos(np.radians(df['sza'])) ** 0.2 + 10)
            m4 = (df['dhi_m'] / df['ghi_m'] <= 1)
            m5 = (df['ghi_m'] / df['ghi_c'] >= 0.1) & (df['ghi_m'] / df['ghi_c'] <= 1.5)
            df = df[m1 & m2 & m3 & m4 & m5]
            # compute variables
            df['e_s'] = 6.112 * np.exp(17.625 * df['temp'] / (df['temp'] + 273.15 - 30.11))
            df['pw_hpa'] = df['e_s'] * df['rh'] / 100
            df = df[df['pw_hpa'] > 0]
            df['kt'] = df['ghi_m'] / df['ghi_c']
            df['kd_obs'] = df['dhi_m'] / df['ghi_m']
            df['e_clear'] = 0.6 + 1.652 * np.sqrt(df['pw_hpa'] / 1013.25) + 0.15 * (np.exp(-alt / H) - 1)
            df['gamma_obs'] = 0.585 * df['kd_obs'] ** 1.75
            df['e_sky_obs'] = (1 - df['gamma_obs']) * df['e_clear'] + df['gamma_obs']
            # hourly average
            hr = df[['kt', 'pw_hpa', 'e_sky_obs']].resample('h').mean().dropna()
            all_hr.append(hr)
    return pd.concat(all_hr)


# Load data
global_df = load_hourly()
# filter kt
global_df = global_df[global_df['kt'] <= 1.3]

# create grids for model surface
pw_vals = np.linspace(global_df['pw_hpa'].min(), global_df['pw_hpa'].max(), 120)
t_vals = np.linspace(0, 1.3, 120)
K_grid, P_grid = np.meshgrid(t_vals, pw_vals)

# Pre-calculate representative altitude for theoretical surface
t_altitudes = np.array(list(station_altitudes.values()))
alt_med = np.median(t_altitudes)

# Plot for each model
for name, model in models.items():
    KD = model(K_grid)
    GAM = 0.585 * KD ** 1.75
    # Use median altitude for theoretical e_clear
    Eclr = 0.6 + 1.652 * np.sqrt(P_grid / 1013.25) + 0.15 * (np.exp(-alt_med / H) - 1)
    Esrf = (1 - GAM) * Eclr + GAM

    # 3D surface
    fig1 = plt.figure(figsize=(6, 6), dpi=300)
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.plot_surface(K_grid, P_grid, Esrf, cmap='plasma', alpha=0.6)
    ax1.scatter(global_df['kt'], global_df['pw_hpa'], global_df['e_sky_obs'],
                s=0.5, c='gray', marker='o', edgecolors='none')
    ax1.set_title(f"3D $e_{{sky}}$ Surface: {name}", fontfamily='serif', fontsize=14)
    ax1.set_xlabel('$k_t$', fontfamily='serif')
    ax1.set_ylabel('$p_w$ (hPa)', fontfamily='serif')
    ax1.set_zlabel('$e_{sky}$', fontfamily='serif')
    fig1.savefig(os.path.join(fig_folder, f'Fig1_{name.replace(" ", "")}.png'), dpi=300)
    plt.close(fig1)

    # Projection k_t–e_sky
    fig2 = plt.figure(figsize=(6, 4), dpi=300)
    ax2 = fig2.add_subplot(111)
    ax2.scatter(global_df['kt'], global_df['e_sky_obs'],
                s=0.5, c='gray', marker='o', edgecolors='none')
    cs2 = ax2.contourf(K_grid, Esrf, P_grid, levels=15, cmap='plasma', alpha=0.7)
    ax2.set_title(f"$e_{{sky}}=f(k_t,p_w)$ Projection: {name}", fontfamily='serif', fontsize=14)
    ax2.set_xlabel('$k_t$', fontfamily='serif')
    ax2.set_ylabel('$e_{sky}$', fontfamily='serif')
    fig2.savefig(os.path.join(fig_folder, f'Fig2_{name.replace(" ", "")}.png'), dpi=300)
    plt.close(fig2)

    # Projection p_w–e_sky
    fig3 = plt.figure(figsize=(6, 4), dpi=300)
    ax3 = fig3.add_subplot(111)
    ax3.scatter(global_df['pw_hpa'], global_df['e_sky_obs'],
                s=0.5, c='gray', marker='o', edgecolors='none')
    cs3 = ax3.contourf(P_grid, Esrf, K_grid, levels=15, cmap='plasma', alpha=0.7)
    ax3.set_title(f"$e_{{sky}}=f(k_t,p_w)$ Projection: {name}", fontfamily='serif', fontsize=14)
    ax3.set_xlabel('$p_w$ (hPa)', fontfamily='serif')
    ax3.set_ylabel('$e_{sky}$', fontfamily='serif')
    fig3.savefig(os.path.join(fig_folder, f'Fig3_{name.replace(" ", "")}.png'), dpi=300)
    plt.close(fig3)
