# model1_rmse.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import os

# 设置全局字体为 Times New Roman
# 1. 全局字体与尺寸设置
plt.rcParams["font.family"] = "Times New Roman"  # Times New Roman
plt.rcParams["axes.labelsize"] = 14             # 坐标轴标签字体大小
plt.rcParams["xtick.labelsize"] = 12            # x轴刻度字体大小
plt.rcParams["ytick.labelsize"] = 12            # y轴刻度字体大小
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Times New Roman"
plt.rcParams["mathtext.it"] = "Times New Roman:italic"
plt.rcParams["mathtext.bf"] = "Times New Roman:bold"
# ---------------------------
# 1. 数据读取与预处理
# ---------------------------
sigma = 5.670374419e-8
G_sc = 1353
H = 8500

station_altitudes = {
    'BON': 213, 'DRA': 1007, 'FPK': 634,
    'GWN': 98, 'PSU': 376, 'SXF': 473,
    'TBL': 1689
}

# 修改为你实际的文件路径
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
        # 针对 DRA 站点删除 2017-2023 年数据，其它站点删除 2023 年数据
        if station == 'DRA':
            dataset = dataset[~((dataset.index.year >= 2017) & (dataset.index.year <= 2023))]
        else:
            dataset = dataset[~(dataset.index.year == 2023)]

        # 只保留降雨量为 0 的日子
        precip_data = precipitation_data[station]
        daily_precip = precip_data.groupby(precip_data.index.date).mean()
        valid_dates = daily_precip[daily_precip == 0].index
        dataset = dataset[dataset.index.normalize().isin(valid_dates)]

        # 计算日序 B
        day_of_year = dataset.index.dayofyear
        B = (day_of_year - 1) * 360 / 365
        G_on = G_sc * (1.000110 + 0.034221 * np.cos(np.radians(B)) +
                       0.001280 * np.sin(np.radians(B)) +
                       0.000719 * np.cos(2 * np.radians(B)) +
                       0.000077 * np.sin(2 * np.radians(B)))

        # 筛选条件
        condition_1 = (dataset['sza'] < 72.5) & (dataset['ghi_m'] > 0) & (dataset['dhi_m'] > 0) & (
                    dataset['dni_m'] > 0) & \
                      (dataset['ghi_m'] / dataset['ghi_c'] < 1.5)
        condition_2 = dataset.index.time >= pd.to_datetime("08:00").time()
        condition_4 = dataset['ghi_m'] < (1.2 * G_on * np.cos(np.radians(dataset['sza'])) ** 1.2 + 50)
        condition_6 = dataset['dni_m'] < (0.95 * G_on * np.cos(np.radians(dataset['sza'])) ** 0.2 + 10)
        condition_8 = (dataset['dhi_m'] / dataset['ghi_m'] <= 1)
        condition_9 = ((dataset['ghi_m'] / dataset['ghi_c']) >= 0.1) & ((dataset['ghi_m'] / dataset['ghi_c']) <= 1.5)
        dataset = dataset[condition_1 & condition_2 & condition_4 & condition_6 & condition_8 & condition_9]

        # 计算 k_d, sqrt_pw 等
        dataset['k_d'] = dataset['dhi_m'] / dataset['ghi_m']
        if 'sqrt_pw' not in dataset.columns:
            dataset['e_s'] = 6.112 * np.exp(17.625 * dataset['temp'] / (dataset['temp'] - 30.11 + 273.15))
            dataset['pw_hpa'] = dataset['e_s'] * dataset['rh'] / 100
            dataset['sqrt_pw'] = np.sqrt(dataset['pw_hpa'] / 1013.25)

        altitude = station_altitudes[station]
        dataset['altitude'] = altitude

        dataset['e_sky'] = dataset['dlw'] / (sigma * ((dataset['temp'] + 273.15) ** 4))
        dataset['e_clear_sky'] = 0.6 + 1.652 * dataset['sqrt_pw'] + 0.15 * (np.exp(-dataset['altitude'] / H) - 1)
        dataset['r'] = (dataset['e_sky'] - dataset['e_clear_sky']) / (1 - dataset['e_clear_sky'])
        dataset = dataset[dataset['r'] <= 1]
        dataset = dataset.replace([np.inf, -np.inf], np.nan).dropna()
        data_list.append(dataset)

# 使所有站点数据量一致
min_size = min(ds.shape[0] for ds in data_list)
data_list = [ds.sample(n=min_size, random_state=42) if ds.shape[0] > min_size else ds for ds in data_list]
df_all = pd.concat(data_list)

# ---------------------------
# 2. 已知最佳参数 (c_1, c_2) for Model1
# ---------------------------
best_c1 = 0.5854
best_c2 = 1.7482

X_model1 = df_all[['k_d', 'sqrt_pw', 'altitude']].values
y_all = df_all['e_sky'].values


def predict_model1(X, c_1, c_2):
    k_d = X[:, 0]
    sqrt_pw = X[:, 1]
    altitude = X[:, 2]
    CF = c_1 * (k_d ** c_2)
    e_clear_sky = 0.6 + 1.652 * sqrt_pw + 0.15 * (np.exp(-altitude / H) - 1)
    return (1 - CF) * e_clear_sky + CF * 1.0


def compute_rmse(y_true, y_pred):
    return np.sqrt(mean_squared_error(y_true, y_pred))


def rmse_on_grid(X, y_true, predict_func, c1_vals, c2_vals):
    RMSE_map = np.zeros((len(c1_vals), len(c2_vals)))
    for i, c1 in enumerate(c1_vals):
        for j, c2 in enumerate(c2_vals):
            y_pred = predict_func(X, c1, c2)
            RMSE_map[i, j] = compute_rmse(y_true, y_pred)
    return RMSE_map


# 网格范围： c_1 in [0,1], c_2 in [1.5,2.5]
c1_vals1 = np.linspace(0.5, 0.7, 10)  # 少一些色带 => 10 个步长
c2_vals1 = np.linspace(1.5, 2.5, 10)

# 计算网格上 RMSE
RMSE_map1 = rmse_on_grid(X_model1, y_all, predict_model1, c1_vals1, c2_vals1)

# 计算最佳点的 RMSE
best_rmse = compute_rmse(y_all, predict_model1(X_model1, best_c1, best_c2))

# ---------------------------
# 3. 绘图
# ---------------------------
plt.figure(figsize=(5.5, 5))
# 设置 levels=10 同时可自定义区间，如想手动指定，可用 np.linspace(RMSE_map1.min(), RMSE_map1.max(), 10)
contour = plt.contourf(c2_vals1, c1_vals1, RMSE_map1, levels=10, cmap='coolwarm')

# 画出最佳点
plt.plot(best_c2, best_c1, '^', color='black', markersize=8)
# 在最佳点附近标注文字： (c_1, c_2) RMSE: ...
plt.text(best_c2, best_c1, f"({best_c1:.3f}, {best_c2:.3f})\nRMSE: {best_rmse:.4f}",
         ha='left', va='bottom', fontsize=12,
         # 适当微调位置可使用: x + 0.02, y - 0.02 等
         )

plt.xlabel('$c_2$')
plt.ylabel('$c_1$')
plt.title('Model1 RMSE Contour')
plt.colorbar(contour, label='RMSE')
plt.tight_layout()
plt.show()

# 若需要保存：
# plt.savefig(os.path.join(output_dir_all, "model1_rmse_contour.png"), dpi=300, bbox_inches='tight')
