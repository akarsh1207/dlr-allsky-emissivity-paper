import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl

mpl.rc('font', family='Times New Roman')
# 1. 全局字体与尺寸设置
plt.rcParams["font.family"] = "Times New Roman"  # Times New Roman
plt.rcParams["axes.labelsize"] = 14             # 坐标轴标签字体大小
plt.rcParams["xtick.labelsize"] = 12            # x轴刻度字体大小
plt.rcParams["ytick.labelsize"] = 12            # y轴刻度字体大小
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Times New Roman"
plt.rcParams["mathtext.it"] = "Times New Roman:italic"
plt.rcParams["mathtext.bf"] = "Times New Roman:bold"
# ----------------------------
# 1) Constants and File Path Definitions
# ----------------------------
sigma = 5.670374419e-8  # Stefan-Boltzmann constant (W·m⁻²·K⁻⁴)
G_sc = 1353            # Solar constant (W/m²)
H = 8500               # Scale height (m)

# Station altitude data (unit: meters)
station_altitudes = {
    'BON': 213,
    'DRA': 1007,
    'FPK': 634,
    'GWN': 98,
    'PSU': 376,
    'SXF': 473,
    'TBL': 1689
}

# Update paths below to actual file locations
precipitation_h5_path = r'/Users/akarsh1207/Desktop/Lab/Coimbra Research Group/Jiedong_work/all_precipitation_data.h5'
file_path = r'/Users/akarsh1207/Desktop/Lab/Coimbra Research Group/Jiedong_work/output.h5'
output_dir_all = r'/Users/akarsh1207/Desktop/Lab/Coimbra Research Group/Jiedong_work/scatter_plots'

# ----------------------------
# 2) Known Global Model Parameters (no re-fitting required)
#    - Model 1: CF = a1 * k_d^b1
#    - Model 2: CF = a2 * k_t^b2
#    - Model 3: CF = a3 * k_t_dni^b3
# ----------------------------
params_model1 = [0.5854, 1.7482]   # 基于 k_d
params_model2 = [0.0902, -1.3879]  # 基于 k_t(ghi)
params_model3 = [0.0790, -0.3221]  # Based on k_t(dni)

# ----------------------------
# 3) Define Three Model Functions (calculate emissivity)
# ----------------------------
def model_function_1(X, a1, b1):
    """
    Model 1: CF = a1 * k_d^b1
    """
    k_d, sqrt_pw, altitude = X[:, 0], X[:, 1], X[:, 2]
    CF = a1 * (k_d ** b1)
    e_clear_sky = 0.6 + 1.652 * sqrt_pw + 0.15 * (np.exp(-altitude / H) - 1)
    return (1 - CF) * e_clear_sky + CF * 1  # e_cloud = 1

def model_function_2(X, a2, b2):
    """
    Model 2: CF = a2 * k_t^b2
    """
    k_t, sqrt_pw, altitude = X[:, 0], X[:, 1], X[:, 2]
    CF = a2 * (k_t ** b2)
    e_clear_sky = 0.6 + 1.652 * sqrt_pw + 0.15 * (np.exp(-altitude / H) - 1)
    return (1 - CF) * e_clear_sky + CF * 1

def model_function_3(X, a3, b3):
    """
    Model 3: CF = a3 * k_t_dni^b3
    """
    k_t_dni, sqrt_pw, altitude = X[:, 0], X[:, 1], X[:, 2]
    CF = a3 * (k_t_dni ** b3)
    e_clear_sky = 0.6 + 1.652 * sqrt_pw + 0.15 * (np.exp(-altitude / H) - 1)
    return (1 - CF) * e_clear_sky + CF * 1

# ----------------------------
# 4) Data Reading and Preprocessing
#    - Read station data from HDF5
#    - Filter conditions: precipitation=0, exclude certain years, sza < 72.5, etc.
#    - Merge into df_all
# ----------------------------
stations = ['BON', 'DRA', 'FPK', 'GWN', 'PSU', 'SXF', 'TBL']
data_list = []

# Read precipitation data
with pd.HDFStore(precipitation_h5_path, 'r') as precip_store:
    precipitation_data = {
        st: precip_store[st].set_index('datetime')['PRECIPITATIONCAL']
        for st in stations
    }

# Read and filter data for each station
with pd.HDFStore(file_path, 'r') as h5_file:
    for station in stations:
        dataset = h5_file[station].copy()
        dataset.index = pd.to_datetime(dataset.index)

        # Add station name as a column
        dataset['station'] = station

        # For DRA station, delete data from 2017-2023
        if station == 'DRA':
            dataset = dataset[~((dataset.index.year >= 2017) & (dataset.index.year <= 2023))]
        dataset = dataset[~(dataset.index.year == 2023)]

        # Only keep dates with zero precipitation
        precip_data = precipitation_data[station]
        daily_precip = precip_data.groupby(precip_data.index.date).mean()
        valid_dates = daily_precip[daily_precip == 0].index
        valid_dates_index = pd.to_datetime(valid_dates).normalize()
        dataset = dataset[dataset.index.normalize().isin(valid_dates_index)]

        # Calculate day angle and G_on (extraterrestrial solar radiation)
        day_of_year = dataset.index.dayofyear
        B = (day_of_year - 1) * 360 / 365
        G_on = G_sc * (1.000110 + 0.034221 * np.cos(np.radians(B)) +
                       0.001280 * np.sin(np.radians(B)) +
                       0.000719 * np.cos(2 * np.radians(B)) +
                       0.000077 * np.sin(2 * np.radians(B)))

        # Multiple filtering conditions (aligned with density_plot.ipynb)
        condition_1 = (dataset['sza'] < 72.5) & (dataset['ghi_m'] > 0) & \
                      (dataset['dhi_m'] > 0) & \
                      (dataset['ghi_m'] / dataset['ghi_c'] < 1.5)
        condition_2 = dataset.index.time >= pd.to_datetime("08:00").time()
        condition_4 = dataset['ghi_m'] < (1.2 * G_on * np.cos(np.radians(dataset['sza'])) ** 1.2 + 50)
        condition_6 = dataset['dni_m'] < (0.95 * G_on * np.cos(np.radians(dataset['sza'])) ** 0.2 + 10)
        condition_8 = (dataset['dhi_m'] / dataset['ghi_m'] <= 1)
        condition_9 = ((dataset['ghi_m'] / dataset['ghi_c']) > 0.1)
        dataset = dataset[condition_1 & condition_2 & condition_4 & condition_6 & condition_8 & condition_9]

        # Calculate related parameters
        dataset['k_t'] = dataset['ghi_m'] / dataset['ghi_c']
        dataset['k_t_dni'] = dataset['dni_m'] / dataset['dni_c']
        dataset['k_d'] = dataset['dhi_m'] / dataset['ghi_m']

        dataset['e_s'] = 6.112 * np.exp(17.625 * dataset['temp'] / (dataset['temp'] - 30.11 + 273.15))
        dataset['pw_hpa'] = dataset['e_s'] * dataset['rh'] / 100

        condition_3 = (dataset['pw_hpa'] >= 0) & (dataset['e_s'] >= 0) & \
                      (dataset['dlw'] > 0) & (dataset['temp'] <= 90) & (dataset['temp'] >= -80)
        condition_7 = (dataset['qc_direct_n'] == 0) & (dataset['qc_dwsolar'] == 0) & \
                      (dataset['qc_diffuse'] == 0) & (dataset['qc_dwir'] == 0) & \
                      (dataset['qc_temp'] == 0) & (dataset['qc_rh'] == 0) & (dataset['qc_pressure'] == 0)
        dataset = dataset[condition_3 & condition_7]

        dataset['sqrt_pw'] = np.sqrt(dataset['pw_hpa'] / 1013.25)
        dataset['altitude'] = station_altitudes[station]

        # Calculate e_sky and e_clear_sky
        dataset['e_sky'] = dataset['dlw'] / (sigma * ((dataset['temp'] + 273.15) ** 4))
        dataset['e_clear_sky'] = 0.6 + 1.652 * dataset['sqrt_pw'] + \
                                  0.15 * (np.exp(-dataset['altitude'] / H) - 1)
        dataset['r'] = (dataset['e_sky'] - dataset['e_clear_sky']) / (1 - dataset['e_clear_sky'])
        dataset = dataset[dataset['r'] <= 1]
        
        # Filter k_d values (aligned with density_plot.ipynb)
        dataset = dataset[dataset['k_d'] <= 1]

        dataset = dataset.replace([np.inf, -np.inf], np.nan).dropna()

        data_list.append(dataset)

# Balance data across stations (use minimum sample size)
station_sizes = {st: ds.shape[0] for st, ds in zip(stations, data_list)}
min_size = min(station_sizes.values())
for i, station in enumerate(stations):
    ds = data_list[i]
    if ds.shape[0] > min_size:
        data_list[i] = ds.sample(n=min_size, random_state=42)

# Merge all station data
df_all = pd.concat(data_list).dropna()

# ----------------------------
# 5) Randomly sample 1000 data points for plotting
# ----------------------------
sample_df = df_all.sample(n=150, random_state=30)

# Construct independent variables for three models
X_model1_sample = sample_df[['k_d', 'sqrt_pw', 'altitude']].values
X_model2_sample = sample_df[['k_t', 'sqrt_pw', 'altitude']].values
X_model3_sample = sample_df[['k_t_dni', 'sqrt_pw', 'altitude']].values

# Calculate predicted emissivity for three models
y_model1 = model_function_1(X_model1_sample, *params_model1)
y_model2 = model_function_2(X_model2_sample, *params_model2)
y_model3 = model_function_3(X_model3_sample, *params_model3)

# Measured emissivity
y_true = sample_df['e_sky'].values

# Measured DLW and temperature
dlw_true = sample_df['dlw'].values
temp_sample = sample_df['temp'].values

# Calculate predicted DLW for three models
dlw_model1 = y_model1 * sigma * ((temp_sample + 273.15) ** 4)
dlw_model2 = y_model2 * sigma * ((temp_sample + 273.15) ** 4)
dlw_model3 = y_model3 * sigma * ((temp_sample + 273.15) ** 4)

# ----------------------------
# 6) Plotting: Subplots (a) Emissivity, (b) DLR
#    Place (a), (b) labels on the top left of each subplot
# ----------------------------
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(8, 10), sharex=False)

# Adjust spacing between subplots
fig.subplots_adjust(hspace=0.25)

# ========== (a) Emissivity comparison plot ==========
min_val_e = min(y_true.min(), y_model1.min(), y_model2.min(), y_model3.min())
max_val_e = max(y_true.max(), y_model1.max(), y_model2.max(), y_model3.max())

# Add some padding to the limits
e_padding = (max_val_e - min_val_e) * 0.02
min_val_e -= e_padding
max_val_e += e_padding

ax1.plot([min_val_e, max_val_e], [min_val_e, max_val_e], 'k-', linewidth=1.5)
ax1.scatter(y_true, y_model1, marker='o', facecolors='none', edgecolors='#1E90FF', s=30,
            label='Model 1', alpha=0.9)
ax1.scatter(y_true, y_model2, marker='o', facecolors='none', edgecolors='#FFA500', s=30,
            label='Model 2', alpha=0.9)
ax1.scatter(y_true, y_model3, marker='o', facecolors='none', edgecolors='#32CD32', s=30,
            label='Model 3', alpha=0.9)

ax1.set_xlabel('Measured emissivity [-]', fontsize=14)
ax1.set_ylabel('Predicted emissivity [-]', fontsize=14)
ax1.set_xlim([min_val_e, max_val_e])
ax1.set_ylim([min_val_e, max_val_e])

# Set nice tick spacing for emissivity
ax1.xaxis.set_major_locator(plt.MultipleLocator(0.05))
ax1.yaxis.set_major_locator(plt.MultipleLocator(0.05))

ax1.legend(loc='upper left', fontsize=11, frameon=True, edgecolor='black', fancybox=False)
ax1.grid(False)

# Place (a) label outside the plot, above the axes on the left
ax1.text(-0.02, 1.02, '(a)', transform=ax1.transAxes,
         ha='left', va='bottom', fontsize=14)

# ========== (b) DLR comparison plot ==========
min_val_dlw = min(dlw_true.min(), dlw_model1.min(), dlw_model2.min(), dlw_model3.min())
max_val_dlw = max(dlw_true.max(), dlw_model1.max(), dlw_model2.max(), dlw_model3.max())

# Add some padding to the limits
dlw_padding = (max_val_dlw - min_val_dlw) * 0.02
min_val_dlw -= dlw_padding
max_val_dlw += dlw_padding

ax2.plot([min_val_dlw, max_val_dlw], [min_val_dlw, max_val_dlw], 'k-', linewidth=1.5)
ax2.scatter(dlw_true, dlw_model1, marker='o', facecolors='none', edgecolors='#1E90FF', s=30,
            label='Model 1', alpha=0.9)
ax2.scatter(dlw_true, dlw_model2, marker='o', facecolors='none', edgecolors='#FFA500', s=30,
            label='Model 2', alpha=0.9)
ax2.scatter(dlw_true, dlw_model3, marker='o', facecolors='none', edgecolors='#32CD32', s=30,
            label='Model 3', alpha=0.9)

ax2.set_xlabel(r'Measured DLR [W/m$^2$]', fontsize=14)
ax2.set_ylabel(r'Predicted DLR [W/m$^2$]', fontsize=14)
ax2.set_xlim([min_val_dlw, max_val_dlw])
ax2.set_ylim([min_val_dlw, max_val_dlw])

# Set nice tick spacing for DLR (every 50 W/m²)
ax2.xaxis.set_major_locator(plt.MultipleLocator(50))
ax2.yaxis.set_major_locator(plt.MultipleLocator(50))

ax2.legend(loc='upper left', fontsize=11, frameon=True, edgecolor='black', fancybox=False)
ax2.grid(False)

# Place (b) label outside the plot, above the axes on the left
ax2.text(-0.02, 1.02, '(b)', transform=ax2.transAxes,
         ha='left', va='bottom', fontsize=14)

plt.tight_layout()

# Save figure
out_fig_path = os.path.join(output_dir_all, 'compaer_model_plots.png')
plt.savefig(out_fig_path, dpi=300, bbox_inches='tight')
plt.show()
