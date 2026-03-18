# model1_rmse.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error
import os

# ---------------------------
# 1. 数据读取与预处理
# ---------------------------
sigma = 5.670374419e-8
G_sc = 1353
H = 8500

station_altitudes = {
    'BON': 213, 'DRA': 1007, 'FPK': 634,
    'GWN': 98, 'PSU': 376, 'SXF': 473, 'TBL': 1689
}

precipitation_h5_path = r'C:\Users\jiw181\PycharmProjects\pythonProject1\testing_focus_on_altitude_correction\precipitation_data\all_precipitation_data.h5'
file_path = r'C:\Users\jiw181\PycharmProjects\pythonProject1\output.h5'
output_dir_all = r'C:\Users\jiw181\PycharmProjects\pythonProject1\2009-2023data_process\regression_results'
os.makedirs(output_dir_all, exist_ok=True)

stations = ['BON', 'DRA', 'FPK', 'GWN', 'PSU', 'SXF', 'TBL']
data_list = []

with pd.HDFStore(precipitation_h5_path, 'r') as precip_store:
    precipitation_data = {
        station: precip_store[station].set_index('datetime')['PRECIPITATIONCAL']
        for station in stations
    }

with pd.HDFStore(file_path, 'r') as h5_file:
    for station in stations:
        dataset = h5_file[station].copy()
        dataset.index = pd.to_datetime(dataset.index)
        if station == 'DRA':
            dataset = dataset[~((dataset.index.year >= 2017) & (dataset.index.year <= 2023))]
        else:
            dataset = dataset[~(dataset.index.year == 2023)]
        precip_data = precipitation_data[station]
        daily_precip = precip_data.groupby(precip_data.index.date).mean()
        valid_dates = daily_precip[daily_precip == 0].index
        valid_dates_index = pd.Index(valid_dates)
        dataset = dataset[dataset.index.normalize().isin(valid_dates_index)]
        day_of_year = dataset.index.dayofyear
        B = (day_of_year - 1)*360/365
        G_on = G_sc * (1.000110 + 0.034221*np.cos(np.radians(B)) +
                       0.001280*np.sin(np.radians(B)) +
                       0.000719*np.cos(2*np.radians(B)) +
                       0.000077*np.sin(2*np.radians(B)))
        condition_1 = (dataset['sza'] < 72.5) & (dataset['ghi_m']>0) & (dataset['dhi_m']>0) & (dataset['dni_m']>0) & (dataset['ghi_m']/dataset['ghi_c']<1.5)
        condition_2 = dataset.index.time >= pd.to_datetime("08:00").time()
        condition_4 = dataset['ghi_m'] < (1.2 * G_on * np.cos(np.radians(dataset['sza']))**1.2 + 50)
        condition_6 = dataset['dni_m'] < (0.95 * G_on * np.cos(np.radians(dataset['sza']))**0.2 + 10)
        condition_8 = ((dataset['dhi_m']/dataset['ghi_m'])<=1)
        condition_9 = ((dataset['ghi_m']/dataset['ghi_c'])>=0.1) & ((dataset['ghi_m']/dataset['ghi_c'])<=1.5)
        dataset = dataset[condition_1 & condition_2 & condition_4 & condition_6 & condition_8 & condition_9]
        dataset['k_d'] = dataset['dhi_m']/dataset['ghi_m']
        dataset['sqrt_pw'] = np.sqrt(dataset['pw_hpa'] if 'pw_hpa' in dataset.columns else  (6.112 * np.exp(17.625*dataset['temp']/(dataset['temp']-30.11+273.15)) * dataset['rh']/100))
        if 'pw_hpa' not in dataset.columns:
            dataset['pw_hpa'] = 6.112 * np.exp(17.625*dataset['temp']/(dataset['temp']-30.11+273.15)) * dataset['rh']/100
        altitude = station_altitudes[station]
        dataset['altitude'] = altitude
        dataset['e_sky'] = dataset['dlw']/(sigma*((dataset['temp']+273.15)**4))
        dataset['e_clear_sky'] = 0.6 + 1.652*dataset['sqrt_pw'] + 0.15*(np.exp(-dataset['altitude']/H)-1)
        dataset['r'] = (dataset['e_sky']-dataset['e_clear_sky'])/(1-dataset['e_clear_sky'])
        dataset = dataset[dataset['r']<=1]
        dataset = dataset.replace([np.inf, -np.inf], np.nan).dropna()
        data_list.append(dataset)

min_size = min([ds.shape[0] for ds in data_list])
data_list = [ds.sample(n=min_size, random_state=42) if ds.shape[0]>min_size else ds for ds in data_list]
df_all = pd.concat(data_list)

# ---------------------------
# 2. Model1 拟合（自变量： [k_d, sqrt_pw, altitude]）
# ---------------------------
y_all = df_all['e_sky'].values
X_model1 = df_all[['k_d', 'sqrt_pw', 'altitude']].values

def model_function_1(X, a, b):
    k_d = X[:, 0]
    sqrt_pw = X[:, 1]
    altitude = X[:, 2]
    CF = a * (k_d ** b)
    e_clear_sky = 0.6 + 1.652*sqrt_pw + 0.15*(np.exp(-altitude/H)-1)
    return (1-CF)*e_clear_sky + CF*1.0

initial_guess = [0.1, 1.0]
params_model1, _ = curve_fit(model_function_1, X_model1, y_all, p0=initial_guess, maxfev=10000)
print(f"Model1 fitted parameters: a = {params_model1[0]:.4f}, b = {params_model1[1]:.4f}")

def predict_model1(X, a, b):
    k_d = X[:, 0]
    sqrt_pw = X[:, 1]
    altitude = X[:, 2]
    CF = a * (k_d ** b)
    e_clear_sky = 0.6 + 1.652*sqrt_pw + 0.15*(np.exp(-altitude/H)-1)
    return (1-CF)*e_clear_sky + CF*1.0

def compute_rmse(y_true, y_pred):
    return np.sqrt(mean_squared_error(y_true, y_pred))

def rmse_on_grid(X, y_true, predict_func, a_vals, b_vals):
    RMSE_map = np.zeros((len(a_vals), len(b_vals)))
    for i, a in enumerate(a_vals):
        for j, b in enumerate(b_vals):
            y_pred = predict_func(X, a, b)
            RMSE_map[i, j] = compute_rmse(y_true, y_pred)
    return RMSE_map

# ---------------------------
# 3. 扫描范围设置（Model1： a in [0,1], b in [1.5,2.5]）及 RMSE 计算
# ---------------------------
a_vals1 = np.linspace(0, 1, 30)
b_vals1 = np.linspace(1.5, 2.5, 30)
RMSE_map1 = rmse_on_grid(X_model1, y_all, predict_model1, a_vals1, b_vals1)

# ---------------------------
# 4. 绘制等值线图
# ---------------------------
plt.figure(figsize=(6,5))
contour = plt.contourf(b_vals1, a_vals1, RMSE_map1, levels=30, cmap='coolwarm')
plt.plot(params_model1[1], params_model1[0], '^', color='black', markersize=8,
         label=f'Best: (a={params_model1[0]:.4f}, b={params_model1[1]:.4f})')
plt.xlabel('b')
plt.ylabel('a')
plt.title('Model1 RMSE Contour')
plt.legend(loc='best')
plt.colorbar(contour, label='RMSE')
plt.tight_layout()
plt.show()
# 可保存： plt.savefig(os.path.join(output_dir_all, "model1_rmse_contour.png"), dpi=300, bbox_inches='tight')
