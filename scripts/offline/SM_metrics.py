import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata, LinearNDInterpolator
from scipy.spatial import Delaunay, cKDTree
from scipy.stats import pearsonr
import warnings
import os
import csv

warnings.filterwarnings("ignore")

# ====================== 1. 核心加速：向量化 TCC 计算函数 ======================
def vectorized_pearsonr(a, b, axis=0):
    mask = ~np.isnan(a) & ~np.isnan(b)
    a_safe = np.where(mask, a, np.nan)
    b_safe = np.where(mask, b, np.nan)
    
    a_mean = np.nanmean(a_safe, axis=axis, keepdims=True)
    b_mean = np.nanmean(b_safe, axis=axis, keepdims=True)
    a_anom = a_safe - a_mean
    b_anom = b_safe - b_mean
    
    cov = np.nansum(a_anom * b_anom, axis=axis)
    std_a = np.sqrt(np.nansum(a_anom**2, axis=axis))
    std_b = np.sqrt(np.nansum(b_anom**2, axis=axis))
    
    with np.errstate(divide='ignore', invalid='ignore'):
        tcc = cov / (std_a * std_b)
        
    valid_count = np.sum(mask, axis=axis)
    tcc = np.where(valid_count >= 3, tcc, np.nan)
    return tcc

# ====================== 2. 路径与季节配置 ======================
path          = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
sim1_file     = path + "colmoff_2001-2023_monmean_nogravel.nc"
sim2_file     = path + "colmoff_2001-2023_monmean_gravel.nc"
cn05_file     = "/share/home/dq117/CN05.1/CN05.1_Tm_2020_daily_025x025.nc"

seasons_config = {
    'MAM': {'obs_suffix': '_MAM.nc', 'model_idx_start': 0, 'model_idx_end': 3, 'name': 'Spring'},
    'JJA': {'obs_suffix': '_JJA.nc', 'model_idx_start': 3, 'model_idx_end': 6, 'name': 'Summer'}
}

# ====================== 3. 自定义土壤层及权重设置 ======================
target_layers = [
    {
        'name': '0-10cm', 
        'obs_var': 'SMs',
        'obs_prefix': '/share/home/dq135/reference/soil_moisture/SMs_1991-2023_GLEAM_v4.2a',
        'wgt': np.array([0.175, 0.276, 0.455, 0.094]),
        'sim_lev_end': 4
    },
    {
        'name': '0-100cm', 
        'obs_var': 'SMrz',
        'obs_prefix': '/share/home/dq135/reference/soil_moisture/SMrz_1991-2023_GLEAM_v4.2a',
        'wgt': np.array([0.0175, 0.0276, 0.0455, 0.0749, 0.1236, 0.2038, 0.336, 0.1711]),
        'sim_lev_end': 8
    }
]

# ====================== 4. 读取公共网格与预构建插值器 ======================
print("加载网格与掩膜数据...")
f_wrf = Dataset(wrfinput_file)
lat2d = np.array(f_wrf.variables['XLAT'][0, :, :])
lon2d = np.array(f_wrf.variables['XLONG'][0, :, :])
f_wrf.close()

ff_temp = Dataset(target_layers[0]['obs_prefix'] + '_MAM.nc')
lat1d = ff_temp.variables['lat'][:]
lon1d = ff_temp.variables['lon'][:]
ff_temp.close()

print("  -> 正在构建空间插值索引 ...")
pts2d = np.column_stack((lon2d.ravel(), lat2d.ravel()))
gx, gy = np.meshgrid(lon1d, lat1d)
target_pts = np.column_stack((gx.ravel(), gy.ravel()))

tri_mesh = Delaunay(pts2d)
kd_tree = cKDTree(pts2d)
_, nearest_inds = kd_tree.query(target_pts)

def fast_regrid(var2d):
    vals = var2d.ravel()
    interp = LinearNDInterpolator(tri_mesh, vals)
    res = interp(gx, gy)
    nan_mask = np.isnan(res)
    if nan_mask.any():
        nearest_res = vals[nearest_inds].reshape(gx.shape)
        res[nan_mask] = nearest_res[nan_mask]
    return res

cn05 = Dataset(cn05_file)
obs_cn05 = np.nanmean(cn05.variables['tm'][:], axis=0)
lat_cn = cn05.variables['lat'][:]
lon_cn = cn05.variables['lon'][:]
cn05.close()

mask_cn05 = ~np.isnan(obs_cn05)
points_cn = np.column_stack((np.tile(lon_cn, len(lat_cn)), np.repeat(lat_cn, len(lon_cn))))
values_cn = mask_cn05.ravel().astype(float)
mask_target = griddata(points_cn, values_cn, (gx, gy), method='nearest') > 0.5

rad = np.pi / 180
cos_lat = np.cos(lat1d * rad)
dlon = lon1d[1] - lon1d[0]
wgt_2d = np.outer(cos_lat, np.ones_like(lon1d)) * dlon
wgt_2d = np.where(mask_target, wgt_2d, np.nan)

def spatial_weighted_avg(data, weight):
    mask = ~np.isnan(data) & ~np.isnan(weight)
    if np.sum(mask) == 0: return np.nan
    return np.average(data[mask], weights=weight[mask])

def weighted_spatial_correlation(obs, sim, wgt):
    o, s, w = obs.flatten(), sim.flatten(), wgt.flatten()
    valid = ~np.isnan(o) & ~np.isnan(s) & ~np.isnan(w)
    if np.sum(valid) < 3: return np.nan, np.nan
    o_v, s_v, w_v = o[valid], s[valid], w[valid]
    mean_o = np.average(o_v, weights=w_v)
    mean_s = np.average(s_v, weights=w_v)
    cov = np.average((o_v - mean_o) * (s_v - mean_s), weights=w_v)
    var_o = np.average((o_v - mean_o)**2, weights=w_v)
    var_s = np.average((s_v - mean_s)**2, weights=w_v)
    scc = cov / np.sqrt(var_o * var_s)
    # p_val 暂时忽略，因为核心仅需相关系数
    return scc, np.nan

# ====================== 5. 核心计算与统计 ======================
table_metrics_data = []
nyear = 17

for season in ['MAM', 'JJA']:
    cfg = seasons_config[season]
    print(f"\n================ 开始处理季节: {season} ================")
    
    for layer in target_layers:
        layer_name = layer['name']
        obs_var = layer['obs_var']
        wgt = layer['wgt']
        sim_lev_end = layer['sim_lev_end']
        obs_file = layer['obs_prefix'] + cfg['obs_suffix']
        
        print(f"  -> [{layer_name}] 开始提取与插值...")

        # 1. 读取观测数据并降维至年际
        ff = Dataset(obs_file)
        var_obs = ff.variables[obs_var]
        obs_offset = 30
        total_steps = nyear * 3
        data_obs = var_obs[obs_offset : obs_offset + total_steps, :, :] 
        if hasattr(var_obs, '_FillValue'): data_obs[data_obs == var_obs._FillValue] = np.nan
        data_obs[np.abs(data_obs) > 1e30] = np.nan
        ff.close()

        # temp_17: shape 为 (17年, lat, lon)
        temp_17 = np.nanmean(data_obs.reshape(nyear, 3, lat1d.size, lon1d.size), axis=1)
        obs_mean = np.nanmean(temp_17, axis=0) 
        mask_obs = ~np.isnan(obs_mean)

        # 2. 读取模拟数据并插值
        sim1_17 = np.zeros((nyear, lat1d.size, lon1d.size))
        sim2_17 = np.zeros((nyear, lat1d.size, lon1d.size))
        idx_s, idx_e = cfg['model_idx_start'], cfg['model_idx_end']

        for yr in range(nyear):
            ds1 = Dataset(sim1_file)
            d1 = ds1.variables['f_h2osoi'][6*yr+idx_s : 6*yr+idx_e, :, :, 0:sim_lev_end]
            if hasattr(ds1.variables['f_h2osoi'], '_FillValue'): d1[d1 == ds1.variables['f_h2osoi']._FillValue] = np.nan
            d1[np.abs(d1) > 1e30] = np.nan
            sim1_mam = np.sum(np.nanmean(d1, axis=0) * wgt, axis=-1)
            sim1_17[yr] = fast_regrid(sim1_mam)
            ds1.close()

            ds2 = Dataset(sim2_file)
            d2 = ds2.variables['f_h2osoi'][6*yr+idx_s : 6*yr+idx_e, :, :, 0:sim_lev_end]
            if hasattr(ds2.variables['f_h2osoi'], '_FillValue'): d2[d2 == ds2.variables['f_h2osoi']._FillValue] = np.nan
            d2[np.abs(d2) > 1e30] = np.nan
            sim2_mam = np.sum(np.nanmean(d2, axis=0) * wgt, axis=-1)
            sim2_17[yr] = fast_regrid(sim2_mam)
            ds2.close()

        for yr in range(nyear):
            sim1_17[yr] = np.where(~np.isnan(temp_17[yr]), sim1_17[yr], np.nan)
            sim2_17[yr] = np.where(~np.isnan(temp_17[yr]), sim2_17[yr], np.nan)

        print(f"  -> [{layer_name}] 开始计算空间平均指标...")
        
        valid_mask = mask_target & mask_obs
        wgt_valid = np.where(valid_mask, wgt_2d, np.nan)

        # ----------------------------------------------------
        # 指标 A: RMSE 的面积加权空间平均
        # ----------------------------------------------------
        rmse1_map = np.sqrt(np.nanmean((sim1_17 - temp_17)**2, axis=0))
        rmse2_map = np.sqrt(np.nanmean((sim2_17 - temp_17)**2, axis=0))
        
        rmse1_avg = spatial_weighted_avg(np.where(valid_mask, rmse1_map, np.nan), wgt_valid)
        rmse2_avg = spatial_weighted_avg(np.where(valid_mask, rmse2_map, np.nan), wgt_valid)

        # ----------------------------------------------------
        # 指标 B: TCC 的面积加权空间平均
        # ----------------------------------------------------
        TCC1_raw = vectorized_pearsonr(sim1_17, temp_17, axis=0)
        TCC2_raw = vectorized_pearsonr(sim2_17, temp_17, axis=0)
        
        TCC1_avg = spatial_weighted_avg(np.where(valid_mask, TCC1_raw, np.nan), wgt_valid)
        TCC2_avg = spatial_weighted_avg(np.where(valid_mask, TCC2_raw, np.nan), wgt_valid)

        # ----------------------------------------------------
        # 指标 C: ACC (异常相关系数) 的时间平均
        # 逻辑：算出气候态(mean)，每年减去mean获得anomaly，算每年的空间相关系数，再时间平均
        # ----------------------------------------------------
        sim1_mean = np.nanmean(sim1_17, axis=0)
        sim2_mean = np.nanmean(sim2_17, axis=0)
        
        acc1_list = []
        acc2_list = []
        
        for yr in range(nyear):
            obs_anom = temp_17[yr] - obs_mean
            sim1_anom = sim1_17[yr] - sim1_mean
            sim2_anom = sim2_17[yr] - sim2_mean
            
            acc1, _ = weighted_spatial_correlation(obs_anom, sim1_anom, wgt_valid)
            acc2, _ = weighted_spatial_correlation(obs_anom, sim2_anom, wgt_valid)
            
            acc1_list.append(acc1)
            acc2_list.append(acc2)
            
        acc1_avg = np.nanmean(acc1_list)
        acc2_avg = np.nanmean(acc2_list)

        # 将结果存入表格数据
        table_metrics_data.append([season, layer_name, "CTL OFF", rmse1_avg, TCC1_avg, acc1_avg])
        table_metrics_data.append([season, layer_name, "EXP OFF", rmse2_avg, TCC2_avg, acc2_avg])

# ====================== 6. 导出表格 ======================
csv_file = "SM_metrics_table.csv"
with open(csv_file, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Season", "Layer", "Experiment", "RMSE", "TCC", "ACC"])
    # 格式化小数位数输出以保证整洁 (保留4位小数)
    for row in table_metrics_data:
        formatted_row = [row[0], row[1], row[2], f"{row[3]:.4f}", f"{row[4]:.4f}", f"{row[5]:.4f}"]
        writer.writerow(formatted_row)
        
print(f"\n✅ 所有数据处理完成！统计指标已导出至表格: {csv_file}")