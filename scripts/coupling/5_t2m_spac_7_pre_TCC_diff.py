import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, RegularGridInterpolator, PchipInterpolator
import warnings
warnings.filterwarnings("ignore")
from scipy.stats import pearsonr, ttest_rel
from scipy.stats import t as t_dist
from scipy.signal import detrend
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.ticker as mticker
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import concurrent.futures
import matplotlib.gridspec as gridspec

import geopandas as gpd
from matplotlib.path import Path

# ====================== 1. 全局配置 ======================
TARGET_VAR = 'pre'        # 切换为 'pre' 画降水组图，切换为 't2m' 画气温组图
APPLY_GRAVEL_MASK = True  # 开关：是否在气温偏差图(FIG5)中应用砾石含量>=0.3的条件遮罩
SIG_LEVEL = 0.05          # 【统一调整为0.05】自定义显著性检验的 p 值阈值大小

PLOT_EXTENT = [80, 130, 10, 55] 
ACC_EXTENT = [80, 130, 10, 55] 

PROJ_TYPE = 'LambertConformal' 
PROJ_CENTRAL_LON = 110.0        
PROJ_CENTRAL_LAT = 40.0         
PROJ_STD_PARALLELS = (30.0, 60.0) 

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['hatch.linewidth'] = 0.4

path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
shp_dir = "../../shapefile_China/"

# 变量独立配置字典
VAR_CONFIG = {
    't2m': {
        'obs_var': 'tm', 'mod_var': 'AT2M', 'convert': lambda x: x - 273.15,
        'obs_file_mam': "/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc",
        'obs_file_jja': "/share/home/dq135/reference/CN05.1_Tm_1991_2023_JJA_025x025.nc",
        'mod_file_ctl_mam': path + "wrfout_2001-2017_MAM_nogravel.nc",
        'mod_file_exp_mam': path + "wrfout_2001-2017_MAM_gravel.nc",
        'mod_file_ctl_jja': path + "wrfout_2001-2017_JJA_nogravel.nc",
        'mod_file_exp_jja': path + "wrfout_2001-2017_JJA_gravel.nc",
        'sp_obs_levels': np.arange(-30, 31, 2),   'sp_obs_cmap': 'RdYlBu_r', 
        'sp_bias_levels': np.arange(-4, 4.1, 0.5), 'sp_bias_cmap': 'RdYlBu_r',
        'sp_diff_levels': np.arange(-1, 1.1, 0.1), 'sp_diff_cmap': 'RdYlBu_r',
        'tcc_levels': np.arange(0, 1.01, 0.1), 'tcc_cmap': 'YlOrRd',
        'tcc_diff_levels': np.arange(-0.2, 0.21, 0.05), 'tcc_diff_cmap': 'RdBu_r'
    },
    'pre': {
        'obs_var': 'pre', 'mod_var': 'PRAVG', 'convert': lambda x: x * 86400.0,
        'obs_file_mam': "/share/home/dq135/reference/CN05.1_Pre_1991_2023_MAM_025x025.nc",
        'obs_file_jja': "/share/home/dq135/reference/CN05.1_Pre_1991_2023_JJA_025x025.nc",
        'mod_file_ctl_mam': path + "wrfout_2001-2017_MAM_nogravel.nc",
        'mod_file_exp_mam': path + "wrfout_2001-2017_MAM_gravel.nc",
        'mod_file_ctl_jja': path + "wrfout_2001-2017_JJA_nogravel.nc",
        'mod_file_exp_jja': path + "wrfout_2001-2017_JJA_gravel.nc",
        'sp_obs_levels': np.arange(0, 11, 1),      'sp_obs_cmap': 'GnBu', 
        'sp_bias_levels': np.arange(-9.0, 9.1, 1.0), 'sp_bias_cmap': 'RdBu',
        'sp_diff_levels': np.arange(-1.0, 1.1, 0.1),'sp_diff_cmap': 'RdBu',
        'tcc_levels': np.arange(0, 1.01, 0.1), 'tcc_cmap': 'YlGnBu',     
        'tcc_diff_levels': np.arange(-0.2, 0.21, 0.05), 'tcc_diff_cmap': 'RdBu_r'       
    }
}
cfg = VAR_CONFIG[TARGET_VAR]
seasons = ['MAM', 'JJA']
nyears = 17 

# ====================== 2. 数据处理函数 ======================
def regrid_rcm2rgrid(var2d, lat2d, lon2d, lat1d, lon1d):
    pts = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    vals = var2d.ravel()
    gx, gy = np.meshgrid(lon1d, lat1d)
    interp = griddata(pts, vals, (gx, gy), method='linear')
    nan_mask = np.isnan(interp)
    if nan_mask.any():
        interp[nan_mask] = griddata(pts, vals, (gx, gy), method='nearest')[nan_mask]
    return interp

def get_yearly_mean(file_path, var_name, is_obs, convert_func=None):
    nc = Dataset(file_path)
    data = nc.variables[var_name][:]
    nc.close()
    if is_obs:
        subset = data[30:30+nyears*3]
        yearly = np.nanmean(subset.reshape(nyears, 3, subset.shape[1], subset.shape[2]), axis=1)
    else:
        subset = data[:nyears*92]
        yearly = np.nanmean(subset.reshape(nyears, 92, subset.shape[1], subset.shape[2]), axis=1)
        if convert_func:
            yearly = convert_func(yearly)
    return yearly

def spatial_weighted_avg(data, weight):
    mask = ~np.isnan(data) & ~np.isnan(weight)
    if np.sum(mask) == 0: return np.nan
    return np.average(data[mask], weights=weight[mask])

def weighted_spatial_correlation(obs, sim, wgt):
    o, s, w = obs.flatten(), sim.flatten(), wgt.flatten()
    valid = ~np.isnan(o) & ~np.isnan(s) & ~np.isnan(w)
    if np.sum(valid) < 3: return np.nan
    o_v, s_v, w_v = o[valid], s[valid], w[valid]
    
    mean_o = np.average(o_v, weights=w_v)
    mean_s = np.average(s_v, weights=w_v)
    cov = np.average((o_v - mean_o) * (s_v - mean_s), weights=w_v)
    var_o = np.average((o_v - mean_o)**2, weights=w_v)
    var_s = np.average((s_v - mean_s)**2, weights=w_v)
    
    if var_o == 0 or var_s == 0: return np.nan
    return cov / np.sqrt(var_o * var_s)

def spatial_acc(anom1, anom2, wgt):
    return weighted_spatial_correlation(anom1, anom2, wgt)

def calc_row_tcc(args):
    i, c_mat, e_mat, o_mat, m_row = args
    lon_size = m_row.shape[0]
    t1, pv1, t2, pv2, pv_diff = [np.full(lon_size, np.nan) for _ in range(5)]
    
    for j in range(lon_size):
        if m_row[j]:
            s1, s2, o = c_mat[:, j], e_mat[:, j], o_mat[:, j]
            v = ~np.isnan(s1) & ~np.isnan(s2) & ~np.isnan(o)
            n = v.sum()
            if n >= 3:
                s1_v, s2_v, o_v = s1[v], s2[v], o[v]
                r12, p1_2tail = pearsonr(s1_v, o_v)
                r13, p2_2tail = pearsonr(s2_v, o_v)
                r23, _  = pearsonr(s1_v, s2_v)
                
                t1[j] = r12
                t2[j] = r13
                
                pv1[j] = p1_2tail / 2.0 if r12 > 0 else 1.0
                pv2[j] = p2_2tail / 2.0 if r13 > 0 else 1.0
                
                if n > 3:
                    rmean = (r12 + r13) / 2.0
                    R = 1.0 - r12**2 - r13**2 - r23**2 + 2.0*r12*r13*r23
                    R = max(0.0, R)
                    den = np.sqrt(2.0 * ((n - 1.0) / (n - 3.0)) * R + (rmean**2) * ((1.0 - r23)**3))
                    if den > 1e-8:
                        t_val = (r13 - r12) * np.sqrt((n - 1.0) * (1.0 + r23)) / den
                        pv_diff[j] = 2.0 * t_dist.sf(np.abs(t_val), n - 3)
                    else:
                        pv_diff[j] = 1.0
                else:
                    pv_diff[j] = 1.0
    return i, t1, pv1, t2, pv2, pv_diff

def get_min_max(data, mask):
    valid_data = data[mask & ~np.isnan(data)]
    if len(valid_data) == 0: 
        return np.nan, np.nan
    return np.min(valid_data), np.max(valid_data)

# ====================== 主程序 ======================
if __name__ == '__main__':
    f_wrf = Dataset(wrfinput_file)
    lat2d = f_wrf.variables['XLAT'][0, :, :]
    lon2d = f_wrf.variables['XLONG'][0, :, :]
    f_wrf.close()

    f_obs_ref = Dataset(cfg['obs_file_mam'])
    lat1d = f_obs_ref.variables['lat'][:]
    lon1d = f_obs_ref.variables['lon'][:]
    f_obs_ref.close()
    lon_grid, lat_grid = np.meshgrid(lon1d, lat1d)

    os.makedirs("../illustration/cpl", exist_ok=True)
    proj = ccrs.LambertConformal(central_longitude=PROJ_CENTRAL_LON, central_latitude=PROJ_CENTRAL_LAT, standard_parallels=PROJ_STD_PARALLELS)

    letters = [chr(97 + i) for i in range(26)]
    
    # ------------------ 计算高分辨砾石数据及降解插值 ------------------
    gravel_cond_mask = None
    qtp_mask = None
    if TARGET_VAR == 't2m':
        print("\n---> 正在提取高分辨率砾石数据，并插值到0.25°网格生成(mean_gravel>=0.3)掩膜...")
        input_grav_file = "/share/home/dq013/zhwei/colm/data/CoLMrawdata/soil/vf_gravels_s.nc" 
        ds_grav = Dataset(input_grav_file)
        if 'longitude' in ds_grav.variables:
            lon_all = ds_grav.variables['longitude'][:]
            lat_all = ds_grav.variables['latitude'][:]
        else:
            lon_all = np.linspace(-180, 180, 86400)
            lat_all = np.linspace(90, -90, 43200)

        lon_idx = np.where((lon_all >= 70) & (lon_all <= 140))[0]
        lat_idx = np.where((lat_all >= 5) & (lat_all <= 60))[0]
        lat_start, lat_end = np.min(lat_idx), np.max(lat_idx) + 1
        lon_start, lon_end = np.min(lon_idx), np.max(lon_idx) + 1

        stride = 5  
        lon_subset = lon_all[lon_start:lon_end:stride]
        lat_subset = lat_all[lat_start:lat_end:stride]

        dz_values = np.array([0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058])
        weight_coefs = dz_values / np.sum(dz_values)
        layer_weights = {
            1: weight_coefs[0] + weight_coefs[1], 2: weight_coefs[2],
            3: weight_coefs[3], 4: weight_coefs[4], 5: weight_coefs[5],
            6: weight_coefs[6], 7: weight_coefs[7], 8: weight_coefs[8] + weight_coefs[9]
        }

        sum_gravel_weighted = np.zeros((len(lat_subset), len(lon_subset)), dtype=np.float32)
        sum_weights         = np.zeros((len(lat_subset), len(lon_subset)), dtype=np.float32)

        for i in range(1, 9):
            var_name = f'vf_gravels_s_l{i}'
            weight = layer_weights[i]
            data = ds_grav.variables[var_name][lat_start:lat_end:stride, lon_start:lon_end:stride]
            if hasattr(data, 'mask'): data = np.ma.filled(data, np.nan)
            data = np.where((data > 1000) | (data < 0), np.nan, data)
            valid_mask = ~np.isnan(data)
            sum_gravel_weighted[valid_mask] += data[valid_mask] * weight
            sum_weights[valid_mask] += weight
        ds_grav.close()

        with np.errstate(divide='ignore', invalid='ignore'):
            mean_gravel_raw = np.where(sum_weights > 0, sum_gravel_weighted / sum_weights, np.nan)

        if lat_subset[0] > lat_subset[-1]:
            lat_subset = lat_subset[::-1]
            mean_gravel_raw = mean_gravel_raw[::-1, :]
        if lon_subset[0] > lon_subset[-1]:
            lon_subset = lon_subset[::-1]
            mean_gravel_raw = mean_gravel_raw[:, ::-1]

        interp_func = RegularGridInterpolator((lat_subset, lon_subset), mean_gravel_raw, method='nearest', bounds_error=False, fill_value=np.nan)
        pts = np.stack((lat_grid, lon_grid), axis=-1)
        mean_gravel_025 = interp_func(pts)
        
        gravel_cond_mask = mean_gravel_025 >= 0.3
        print("砾石条件掩膜提取完毕！\n")
        
        print("---> 读取青藏高原 Shapefile 以生成严格边界掩膜...")
        qtp_shp_path = "../../shapefile_China/TPBoundary_new_2021/TPBoundary_new(2021).shp"
        qtp_gdf = gpd.read_file(qtp_shp_path)

        poly = qtp_gdf.geometry.union_all()
        polygons = list(poly.geoms) if poly.geom_type == 'MultiPolygon' else [poly]

        points_flat = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))
        qtp_mask_np = np.zeros(lon_grid.size, dtype=bool)

        for p in polygons:
            path_obj = Path(np.asarray(p.exterior.coords))
            qtp_mask_np |= path_obj.contains_points(points_flat, radius=0)

        qtp_mask = qtp_mask_np.reshape(lon_grid.shape)
        print("青藏高原边界掩膜提取完毕！\n")
    
    # ------------------ 画布与独立网格排版初始化 ------------------
    if TARGET_VAR == 't2m':
        fig = plt.figure(figsize=(15, 8)) 
        gs = gridspec.GridSpec(2, 3, figure=fig, left=0.06, right=0.98, top=0.94, bottom=0.18, hspace=0.15, wspace=0.14)
    else:
        fig2 = plt.figure(figsize=(15, 7))
        
        gs2_left  = gridspec.GridSpec(2, 1, figure=fig2, left=0.05, right=0.24, top=0.95, bottom=0.10, hspace=0.14, height_ratios=[0.817, 1])
        gs2_right = gridspec.GridSpec(2, 3, figure=fig2, left=0.28, right=0.98, top=0.95, bottom=0.18, hspace=0.08, wspace=0.06)

    cf_sp_obs = cf_sp_bias = cf_sp_diff = None
    cf_tcc = cf_tcc_diff = None

    panel_idx1 = 0
    panel_idx2 = 0
    
    ax_dist_list = []
    ax_map_tcc_list = [[], []]

    def apply_map_style(ax):
        ax.set_extent(PLOT_EXTENT, crs=ccrs.PlateCarree())
        for shp_file, lw in [("province.shp", 0.4), ("china.shp", 0.6), ("south_china_sea.shp", 0.8)]:
            p_shp = os.path.join(shp_dir, shp_file)
            if os.path.exists(p_shp):
                color = 'blue' if 'river' in shp_file else 'black'
                try:
                    reader = shpreader.Reader(p_shp)
                    ax.add_geometries(reader.geometries(), crs=ccrs.PlateCarree(), facecolor='none', edgecolor=color, linewidth=lw, zorder=3)
                except: pass
        
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, 
                          linewidth=0.6, color='gray', alpha=0.5, linestyle='--', zorder=2)
        gl.xlocator = mticker.FixedLocator(np.arange(60, 140, 10))
        gl.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
        
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = True
        gl.bottom_labels = True
        
        gl.xlabel_style = {'size': 9, 'rotation': 0, 'va': 'top', 'ha': 'center'}
        gl.ylabel_style = {'size': 9, 'rotation': 0, 'va': 'center', 'ha': 'right'}
        try:
            gl.padding = 3
        except Exception:
            pass

    for row_idx, season in enumerate(seasons):
        season_name = "Spring" if season == "MAM" else "Summer"
        print(f"---> 开始处理 {TARGET_VAR} 季节: {season}")
        
        obs_17 = get_yearly_mean(cfg[f'obs_file_{season.lower()}'], cfg['obs_var'], is_obs=True)
        ctl_17_raw = get_yearly_mean(cfg[f'mod_file_ctl_{season.lower()}'], cfg['mod_var'], is_obs=False, convert_func=cfg['convert'])
        exp_17_raw = get_yearly_mean(cfg[f'mod_file_exp_{season.lower()}'], cfg['mod_var'], is_obs=False, convert_func=cfg['convert'])
        
        ctl_17, exp_17 = np.zeros_like(obs_17), np.zeros_like(obs_17)
        for yr in range(nyears):
            ctl_17[yr] = regrid_rcm2rgrid(ctl_17_raw[yr], lat2d, lon2d, lat1d, lon1d)
            exp_17[yr] = regrid_rcm2rgrid(exp_17_raw[yr], lat2d, lon2d, lat1d, lon1d)
            ctl_17[yr] = np.where(~np.isnan(obs_17[yr]), ctl_17[yr], np.nan)
            exp_17[yr] = np.where(~np.isnan(obs_17[yr]), exp_17[yr], np.nan)
            
        obs_mean = np.nanmean(obs_17, axis=0)
        ctl_mean = np.nanmean(ctl_17, axis=0)
        exp_mean = np.nanmean(exp_17, axis=0)
        
        bias_ctl = ctl_mean - obs_mean
        bias_exp = exp_mean - obs_mean
        diff_sp = exp_mean - ctl_mean
        mask_obs = ~np.isnan(obs_mean)

        rad = np.pi / 180
        wgt_2d = np.where(mask_obs, np.outer(np.cos(lat1d * rad), np.ones_like(lon1d)), np.nan)

        # === 修改：将显著性检验 (T-test) 及其掩膜步骤提前执行 ===
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, p_bias_ctl = ttest_rel(ctl_17, obs_17, axis=0, nan_policy='omit')
            _, p_bias_exp = ttest_rel(exp_17, obs_17, axis=0, nan_policy='omit')
            _, p_diff_sp  = ttest_rel(exp_17, ctl_17, axis=0, nan_policy='omit')

        # 仅保留通过显著性检验的网格偏差值（未通过的设为 np.nan）
        bias_ctl_sig = np.where((p_bias_ctl < SIG_LEVEL) & mask_obs, bias_ctl, np.nan)
        bias_exp_sig = np.where((p_bias_exp < SIG_LEVEL) & mask_obs, bias_exp, np.nan)
        diff_sp_sig  = np.where((p_diff_sp < SIG_LEVEL) & mask_obs, diff_sp, np.nan)
        # ====================================================

        # === 修改：针对通过显著性检验的有效网格计算并打印均值 ===
        if TARGET_VAR == 't2m' and gravel_cond_mask is not None and qtp_mask is not None:
            # 定义排除掩膜: lat 25~29.5N, lon 91~98E
            exclude_mask = (lat_grid >= 25.0) & (lat_grid <= 29.5) & (lon_grid >= 91.0) & (lon_grid <= 98.0)
            
            # 基础条件：有效观测 + 砾石 >= 0.3 + 在高原Shapefile内部 + 排除特殊矩形区域
            final_cond_mask = mask_obs & gravel_cond_mask & qtp_mask & (~exclude_mask)
            
            wgt_2d_region = np.where(final_cond_mask, wgt_2d, np.nan)
            
            # OBS使用原始值(不需要检验)，模式偏差使用过滤后的_sig值
            sub_gravel_obs = spatial_weighted_avg(obs_mean, wgt_2d_region)
            sub_gravel_bias = spatial_weighted_avg(bias_ctl_sig, wgt_2d_region)
            sub_gravel_diff = spatial_weighted_avg(diff_sp_sig, wgt_2d_region)
            
            # 提取最值 (同样使用 _sig 变量)
            obs_min, obs_max = get_min_max(obs_mean, final_cond_mask)
            bias_min, bias_max = get_min_max(bias_ctl_sig, final_cond_mask)
            diff_min, diff_max = get_min_max(diff_sp_sig, final_cond_mask)
            
            print(f"    [Gravel>=0.3 & QTP Shapefile & 扣除(25-29.5N, 91-98E) & p<{SIG_LEVEL}] {season} 2m气温 区域特征:")
            print(f"        OBS     -> 加权平均: {sub_gravel_obs:.2f} °C | 最小值: {obs_min:.2f} °C | 最大值: {obs_max:.2f} °C")
            print(f"        CTL-OBS -> 加权平均: {sub_gravel_bias:.2f} °C | 最小值: {bias_min:.2f} °C | 最大值: {bias_max:.2f} °C")
            print(f"        EXP-CTL -> 加权平均: {sub_gravel_diff:.2f} °C | 最小值: {diff_min:.2f} °C | 最大值: {diff_max:.2f} °C")
        # ======================================================================

        # -------------------- ACC / TCC 基础计算 --------------------
        if TARGET_VAR == 't2m':
            def apply_detrend(data_3d):
                out = np.copy(data_3d)
                vy, vx = np.where(mask_obs)
                for y, x in zip(vy, vx):
                    ts = data_3d[:, y, x]
                    if (~np.isnan(ts)).sum() > 2:
                        out[~np.isnan(ts), y, x] = detrend(ts[~np.isnan(ts)]) + np.mean(ts[~np.isnan(ts)])
                return out
            obs_17_dt = apply_detrend(obs_17)
            ctl_17_dt = apply_detrend(ctl_17)
            exp_17_dt = apply_detrend(exp_17)
            
            if APPLY_GRAVEL_MASK:
                wgt_2d_metrics = np.where(gravel_cond_mask, wgt_2d, np.nan)
            else:
                wgt_2d_metrics = wgt_2d
        else:
            obs_17_dt, ctl_17_dt, exp_17_dt = obs_17, ctl_17, exp_17
            wgt_2d_metrics = wgt_2d

        acc_ctl, acc_exp = np.zeros(nyears), np.zeros(nyears)
        for yr in range(nyears):
            acc_ctl[yr] = spatial_acc(obs_17_dt[yr]-np.nanmean(obs_17_dt,axis=0), ctl_17_dt[yr]-np.nanmean(ctl_17_dt,axis=0), wgt_2d_metrics)
            acc_exp[yr] = spatial_acc(obs_17_dt[yr]-np.nanmean(obs_17_dt,axis=0), exp_17_dt[yr]-np.nanmean(exp_17_dt,axis=0), wgt_2d_metrics)

        rmse_ctl = np.sqrt(spatial_weighted_avg((ctl_mean - obs_mean)**2, wgt_2d_metrics))
        rmse_exp = np.sqrt(spatial_weighted_avg((exp_mean - obs_mean)**2, wgt_2d_metrics))
        tcc_sp_ctl = weighted_spatial_correlation(obs_mean, ctl_mean, wgt_2d_metrics)
        tcc_sp_exp = weighted_spatial_correlation(obs_mean, exp_mean, wgt_2d_metrics)
        acc_time_ctl = np.nanmean(acc_ctl)
        acc_time_exp = np.nanmean(acc_exp)
        
        metrics_str_ctl = f"RMSE: {rmse_ctl:.2f}\nTCC: {tcc_sp_ctl:.2f}\nACC: {acc_time_ctl:.2f}"
        metrics_str_exp = f"RMSE: {rmse_exp:.2f}\nTCC: {tcc_sp_exp:.2f}\nACC: {acc_time_exp:.2f}"
        
        if TARGET_VAR == 'pre':
            TCC1, P1, TCC2, P2, P_DIFF = [np.full_like(mask_obs, np.nan, dtype=float) for _ in range(5)]
            with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
                args = [(i, ctl_17_dt[:, i, :], exp_17_dt[:, i, :], obs_17_dt[:, i, :], mask_obs[i, :]) for i in range(lat1d.size)]
                for i, t1, pv1, t2, pv2, pv_diff in executor.map(calc_row_tcc, args):
                    TCC1[i, :], P1[i, :], TCC2[i, :], P2[i, :], P_DIFF[i, :] = t1, pv1, t2, pv2, pv_diff
            
            diff_tcc = TCC2 - TCC1

            sig_mask_ctl_tcc = mask_obs & ~np.isnan(TCC1) & (P1 < SIG_LEVEL)
            sig_mask_exp_tcc = mask_obs & ~np.isnan(TCC2) & (P2 < SIG_LEVEL)
            sig_mask_union = sig_mask_ctl_tcc | sig_mask_exp_tcc
        
        # ==================== 核心绘图逻辑 ====================
        if TARGET_VAR == 't2m':
            sp_titles = ["OBS", "CTL CPL - OBS", "(EXP - CTL) CPL"]
            sp_data = [obs_mean, bias_ctl_sig, diff_sp_sig]
            sp_levels = [cfg['sp_obs_levels'], cfg['sp_bias_levels'], cfg['sp_diff_levels']]
            sp_cmaps = [cfg['sp_obs_cmap'], cfg['sp_bias_cmap'], cfg['sp_diff_cmap']]

            for col in range(3):
                ax = fig.add_subplot(gs[row_idx, col], projection=proj)
                apply_map_style(ax)
                
                plot_data = sp_data[col]
                
                if col in [1, 2] and APPLY_GRAVEL_MASK:
                    ax.contourf(lon_grid, lat_grid, np.where((~gravel_cond_mask) & mask_obs, 1, np.nan), 
                                levels=[0, 2], colors=['lightgray'], transform=ccrs.PlateCarree(), zorder=1)
                    plot_data = np.where(gravel_cond_mask & mask_obs, plot_data, np.nan)
                    ax.text(0.97, 0.96, "%gravel≥0.3", transform=ax.transAxes, fontsize=10, 
                            zorder=5, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1.5))

                cf = ax.contourf(lon_grid, lat_grid, plot_data, levels=sp_levels[col], cmap=sp_cmaps[col], extend='both', transform=ccrs.PlateCarree(), zorder=1, antialiased=False)
                try:
                    for c_ in cf.collections:
                        c_.set_edgecolor("face")
                        c_.set_linewidth(1e-8)
                except AttributeError:
                    cf.set_edgecolor("face")
                    cf.set_linewidth(1e-8)
                
                if row_idx == 0: ax.set_title(sp_titles[col], loc='left', fontsize=11)
                if col == 0: ax.text(-0.2, 0.5, season_name, transform=ax.transAxes, fontsize=12, va='center', rotation=90)
                
                if row_idx == 0 and col == 2:
                    ax.text(1.0, 1.02, "℃", transform=ax.transAxes, fontsize=11, va='bottom', ha='right')

                ax.text(0.03, 0.96, f"({letters[panel_idx1]})", transform=ax.transAxes, fontsize=11, zorder=5, va='top', ha='left', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1.5))
                
                if col == 1: 
                    ax.text(0.03, 0.04, metrics_str_ctl, transform=ax.transAxes, fontsize=9, va='bottom', ha='left', zorder=5, bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', lw=0.5, pad=2.5))
                elif col == 2: 
                    ax.text(0.03, 0.04, metrics_str_exp, transform=ax.transAxes, fontsize=9, va='bottom', ha='left', zorder=5, bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', lw=0.5, pad=2.5))
                
                if col == 0: cf_sp_obs = cf
                elif col == 1: cf_sp_bias = cf
                elif col == 2: cf_sp_diff = cf
                panel_idx1 += 1

        elif TARGET_VAR == 'pre':
            ax_dist = fig2.add_subplot(gs2_left[row_idx, 0])
            ax_dist_list.append(ax_dist)
            
            valid_mask_ctl = mask_obs & ~np.isnan(TCC1)
            valid_mask_exp = mask_obs & ~np.isnan(TCC2)
            
            valid_tcc_ctl = TCC1[valid_mask_ctl]
            valid_tcc_exp = TCC2[valid_mask_exp]
            
            if len(valid_tcc_ctl) > 1 and len(valid_tcc_exp) > 1:
                global_min = min(np.nanmin(valid_tcc_ctl), np.nanmin(valid_tcc_exp))
                min_xlim = np.floor(global_min * 10) / 10.0
                max_xlim = 1.0
                
                bin_width = 0.05
                bins_edges = np.arange(min_xlim, max_xlim + bin_width, bin_width)
                bin_centers = bins_edges[:-1] + bin_width / 2.0
                
                counts_ctl, _ = np.histogram(valid_tcc_ctl, bins=bins_edges)
                counts_exp, _ = np.histogram(valid_tcc_exp, bins=bins_edges)
                
                freq_ctl = counts_ctl / len(valid_tcc_ctl)
                freq_exp = counts_exp / len(valid_tcc_exp)
                
                x_smooth = np.linspace(min_xlim, max_xlim, 300)
                spline_ctl = PchipInterpolator(bin_centers, freq_ctl)
                spline_exp = PchipInterpolator(bin_centers, freq_exp)
                
                y_smooth_ctl = np.clip(spline_ctl(x_smooth), 0, None)
                y_smooth_exp = np.clip(spline_exp(x_smooth), 0, None)

                sig_all_pos = np.concatenate([TCC1[(P1 < SIG_LEVEL) & (TCC1 > 0) & mask_obs], 
                                              TCC2[(P2 < SIG_LEVEL) & (TCC2 > 0) & mask_obs]])
                stat_threshold = np.nanmin(sig_all_pos) if len(sig_all_pos) > 0 else 0.412
                
                mean_sig_ctl = np.nanmean(TCC1[(TCC1 >= stat_threshold) & mask_obs])
                mean_sig_exp = np.nanmean(TCC2[(TCC2 >= stat_threshold) & mask_obs])

                ax_dist.plot(x_smooth, y_smooth_ctl, color='forestgreen', ls='-', lw=1.8, alpha=0.9, 
                             label=f'CTL ({mean_sig_ctl:.2f})', zorder=3)
                ax_dist.plot(x_smooth, y_smooth_exp, color='darkorange', ls='--', lw=1.8, alpha=0.9, 
                             label=f'EXP ({mean_sig_exp:.2f})', zorder=2)
                
                ax_dist.axvline(x=stat_threshold, color='red', linestyle=':', linewidth=2.0, zorder=1)
                    
                ax_dist.set_xlim(min_xlim, max_xlim)
                ax_dist.xaxis.set_major_locator(mticker.MaxNLocator(nbins=6, prune=None))
                ax_dist.set_ylim(bottom=0)
            
            ax_dist.set_xlabel('') 
            ax_dist.set_ylabel('Probability', fontsize=11, labelpad=8)
            if row_idx == 0: ax_dist.set_title("TCC", loc='left', fontsize=11)
            ax_dist.text(0.03, 0.96, f"({letters[panel_idx2]}) {season_name}", transform=ax_dist.transAxes, 
                         fontsize=11, zorder=5, va='top', ha='left', fontweight='normal')
            
            ax_dist.legend(loc='best', frameon=False, fontsize=9)
            ax_dist.grid(True, ls='--', alpha=0.4)
            panel_idx2 += 1
            
            tcc_titles = ["CTL CPL TCC", "EXP CPL TCC", "TCC Difference"]
            
            tcc_data = [TCC1, TCC2, diff_tcc]
            tcc_levels = [cfg['tcc_levels'], cfg['tcc_levels'], cfg['tcc_diff_levels']]
            tcc_cmaps = [cfg['tcc_cmap'], cfg['tcc_cmap'], cfg['tcc_diff_cmap']]
            
            for col in range(3):
                ax = fig2.add_subplot(gs2_right[row_idx, col], projection=proj)
                ax_map_tcc_list[row_idx].append(ax)
                apply_map_style(ax)
                
                plot_data = tcc_data[col]
                
                # ==================== 修改开始 ====================
                plot_data_contour = np.copy(plot_data) 
                
                # 如果是最后一列(TCC Difference)，通过显著性掩膜将非打点区域设为 nan
                if col == 2:
                    sig_diff = sig_mask_union & mask_obs & ~np.isnan(plot_data)
                    plot_data_contour = np.where(sig_diff, plot_data, np.nan)
                
                # 绘图时使用 plot_data_contour 替代原来的 plot_data
                cf = ax.contourf(lon_grid, lat_grid, plot_data_contour, levels=tcc_levels[col], cmap=tcc_cmaps[col], extend='both', transform=ccrs.PlateCarree(), zorder=1, antialiased=False)
                # ==================== 修改结束 ====================

                try:
                    for c_ in cf.collections:
                        c_.set_edgecolor("face")
                        c_.set_linewidth(1e-8)
                except AttributeError:
                    cf.set_edgecolor("face")
                    cf.set_linewidth(1e-8)
                
                if col < 2:
                    P_val = P1 if col == 0 else P2
                    sig_diff = mask_obs & ~np.isnan(plot_data) & (P_val < SIG_LEVEL)
                    
                    stride_dot = 4
                    sparse_mask = np.zeros_like(sig_diff, dtype=bool)
                    sparse_mask[::stride_dot, ::stride_dot] = True
                    final_sig_mask = sig_diff & sparse_mask

                    lon_sig = lon_grid[final_sig_mask]
                    lat_sig = lat_grid[final_sig_mask]
                    if len(lon_sig) > 0:
                        ax.scatter(lon_sig, lat_sig, s=3.5, color='black', marker='o', edgecolors='none', alpha=0.8, transform=ccrs.PlateCarree(), zorder=3)
                else:
                    # ==================== 修改开始 ====================
                    # 移除了原本的白色打点代码 (ax.scatter)，因为现在只在这些区域填色
                    pass
                    # ==================== 修改结束 ====================

                if row_idx == 0: ax.set_title(tcc_titles[col], loc='left', fontsize=11)
                ax.text(0.03, 0.96, f"({letters[panel_idx2]})", transform=ax.transAxes, fontsize=11, zorder=5, va='top', ha='left', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1.5))
                
                if col == 0:
                    ax.text(0.03, 0.04, metrics_str_ctl, transform=ax.transAxes, fontsize=9, va='bottom', ha='left', zorder=5, bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', lw=0.5, pad=2.5))
                    cf_tcc = cf
                elif col == 1:
                    ax.text(0.03, 0.04, metrics_str_exp, transform=ax.transAxes, fontsize=9, va='bottom', ha='left', zorder=5, bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', lw=0.5, pad=2.5))
                elif col == 2:
                    cf_tcc_diff = cf
                    
                    valid_diff_masked = ~np.isnan(plot_data) & mask_obs
                    total_wgt_for_fraction = np.nansum(wgt_2d[valid_diff_masked])
                    
                    area = [0, 0]
                    if total_wgt_for_fraction > 0:
                        area[0] = np.nansum(wgt_2d[(plot_data < 0) & valid_diff_masked]) / total_wgt_for_fraction * 100.0
                        area[1] = np.nansum(wgt_2d[(plot_data > 0) & valid_diff_masked]) / total_wgt_for_fraction * 100.0

                    ax_inset = inset_axes(ax, width="70%", height="4.5%", loc='lower left', 
                                          bbox_to_anchor=(0.05, 0.08, 0.9, 0.9), bbox_transform=ax.transAxes, borderpad=0)
                    ax_inset.patch.set_alpha(0.0)
                    left = 0
                    for val, c in zip(area, ['lightblue', 'lightcoral']):
                        ax_inset.barh(0, val, left=left, color=c, height=0.7, edgecolor='black', linewidth=0.6)
                        if val > 5: 
                            ax_inset.text(left + val/2, 0, f"{val:.1f}%", va='center', ha='center', color='black', fontsize=7, fontweight='normal')
                        left += val
                    for spine in ax_inset.spines.values(): spine.set_visible(False)
                    ax_inset.set(xlim=(0, 100), ylim=(-0.5, 0.5), xticks=[], yticks=[])
                    ax_inset.set_title('Areal Fraction (%)', fontsize=7, fontweight='normal', pad=2)
                    ax_inset.legend(handles=[Rectangle((0,0),1,1, facecolor=c, edgecolor='black', lw=0.6) for c in ['lightblue', 'lightcoral']],
                                    labels=['Negative', 'Positive'], loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=7, frameon=False)
                panel_idx2 += 1


    # ====================== 3. Colorbar 极简排版与保存 ======================
    if TARGET_VAR == 't2m':
        plt.figure(fig.number)
        plt.draw() 
        pos0 = fig.axes[3].get_position()
        pos1 = fig.axes[4].get_position()
        pos2 = fig.axes[5].get_position()
        
        cbar_y = 0.10
        cbar_h = 0.012
        
        cax_obs = fig.add_axes([pos0.x0, cbar_y, pos0.width, cbar_h])
        plt.colorbar(cf_sp_obs, cax=cax_obs, orientation='horizontal')
        
        cax_bias = fig.add_axes([pos1.x0, cbar_y, pos1.width, cbar_h])
        plt.colorbar(cf_sp_bias, cax=cax_bias, orientation='horizontal')

        cax_diff = fig.add_axes([pos2.x0, cbar_y, pos2.width, cbar_h])
        plt.colorbar(cf_sp_diff, cax=cax_diff, orientation='horizontal')
        
        save_path = f"../illustration/cpl/FIG5.{TARGET_VAR.upper()}_clim_bias.pdf"
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"气温气候态偏差图已保存至：{save_path}")

    elif TARGET_VAR == 'pre':
        plt.figure(fig2.number)
        plt.draw()
        
        pos_tcc_ctl = ax_map_tcc_list[1][0].get_position() 
        pos_tcc_exp = ax_map_tcc_list[1][1].get_position()
        pos_tcc_diff = ax_map_tcc_list[1][2].get_position()
        
        cbar_y2 = 0.10   
        cbar_h2 = 0.015  
        
        cax_tcc = fig2.add_axes([pos_tcc_ctl.x0, cbar_y2, pos_tcc_exp.x1 - pos_tcc_ctl.x0, cbar_h2])
        plt.colorbar(cf_tcc, cax=cax_tcc, orientation='horizontal')

        cax_tcc_diff = fig2.add_axes([pos_tcc_diff.x0, cbar_y2, pos_tcc_diff.width, cbar_h2])
        plt.colorbar(cf_tcc_diff, cax=cax_tcc_diff, orientation='horizontal')
        
        save_path2 = f"../illustration/cpl/FIG7.{TARGET_VAR.upper()}_TCC_prob.pdf"
        fig2.savefig(save_path2, dpi=300, bbox_inches='tight')
        plt.close(fig2)
        print(f"降水 TCC 综合图版已保存至：{save_path2}")