import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, RegularGridInterpolator # === 导入插值器 ===
import warnings
from scipy.stats import ttest_rel
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.ticker as mticker
import cftime
import os

warnings.filterwarnings("ignore")

# ====================== 字体全局设置 ======================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ====================== 1. 路径与配置 ======================
path          = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
sim_ctl_file  = path + "colmrun_2001-2017_nogravel.nc"
sim_exp_file  = path + "colmrun_2001-2017_gravel.nc"
cn05_file     = "/share/home/dq117/CN05.1/CN05.1_Tm_2020_daily_025x025.nc"
shp_dir       = "../../shapefile_China/"
out_dir       = "./figs/"

os.makedirs(out_dir, exist_ok=True)
os.environ['CARTOPY_OFFLINE'] = 'true'

obs_files = {
    '0cm': "/share/home/dq135/reference/soil_temperature/HOM_TS0cm_daily_025x025_1960_2017.nc",
    '15cm': "/share/home/dq135/reference/soil_temperature/HOM_TS15cm_daily_025x025_1960_2017.nc",
    '40cm': "/share/home/dq135/reference/soil_temperature/HOM_TS40cm_daily_025x025_1960_2017.nc"
}

layers_config = [
    {'name': '0cm', 'colm_depth': 0.00},
    {'name': '15cm', 'colm_depth': 0.15},
    {'name': '40cm', 'colm_depth': 0.40}
]
colm_lev = np.array([0.0071, 0.0279, 0.0623, 0.1189, 0.2122, 0.3661, 0.6198, 1.0380, 1.7276, 2.8646])
nyear = 17

# ====================== 2. 读取公共网格与 Mask ======================
print("加载公共网格与 Mask 模板...")
f_wrf = Dataset(wrfinput_file)
lat2d = f_wrf.variables['XLAT'][0, :, :]
lon2d = f_wrf.variables['XLONG'][0, :, :]
f_wrf.close()

f_obs_ref = Dataset(obs_files['0cm'])
lat1d = f_obs_ref.variables['latitude'][:]
lon1d = f_obs_ref.variables['longitude'][:]
f_obs_ref.close()

cn05 = Dataset(cn05_file)
obs_cn05 = np.nanmean(cn05.variables['tm'][:], axis=0)
lat_cn = cn05.variables['lat'][:]
lon_cn = cn05.variables['lon'][:]
cn05.close()

mask_cn05 = ~np.isnan(obs_cn05)
grid_x, grid_y = np.meshgrid(lon1d, lat1d)
points = np.column_stack((np.tile(lon_cn, len(lat_cn)), np.repeat(lat_cn, len(lon_cn))))
values = mask_cn05.ravel().astype(float)
mask_target = griddata(points, values, (grid_x, grid_y), method='nearest') > 0.5

# ====================== 2.5 提取高分辨率砾石数据并生成条件掩膜 ======================
print("提取高分辨率砾石数据，并生成 (砾石>=0.3 + 经度<=110°E) 的条件掩膜...")
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
pts = np.stack((grid_y, grid_x), axis=-1)
mean_gravel_025 = interp_func(pts)

# 条件1: 砾石含量 >= 0.3
gravel_cond_mask = mean_gravel_025 >= 0.3
# 条件2: 经度范围小于等于 110°E
lon_cond_mask = grid_x <= 110.0

# 最终复合掩膜
final_cond_mask = mask_target & gravel_cond_mask & lon_cond_mask

# 纬度面积加权矩阵
rad = np.pi / 180
wgt_2d = np.cos(grid_y * rad)
print("掩膜提取与加权矩阵准备完毕！\n")

# ====================== 3. 核心计算与功能函数 ======================
def regrid_rcm2rgrid(var2d, lat_in, lon_in, lat_out, lon_out):
    pts = np.column_stack((lon_in.ravel(), lat_in.ravel()))
    vals = var2d.ravel()
    gx, gy = np.meshgrid(lon_out, lat_out)
    interp = griddata(pts, vals, (gx, gy), method='linear')
    nan_mask = np.isnan(interp)
    if nan_mask.any():
        interp[nan_mask] = griddata(pts, vals, (gx, gy), method='nearest')[nan_mask]
    return interp

def extract_obs_yearly(file_path, target_months, start_yr=2001, end_yr=2017):
    ds = Dataset(file_path)
    times = ds.variables['time']
    dates = cftime.num2date(times[:], units=times.units, calendar=times.calendar if hasattr(times, 'calendar') else 'standard')
    
    years = np.array([d.year for d in dates])
    months = np.array([d.month for d in dates])
    
    mask_time = (years >= start_yr) & (years <= end_yr) & np.isin(months, target_months)
    filtered_years = years[mask_time] 
    
    data = ds.variables['TS'][mask_time, :, :]
    if hasattr(ds.variables['TS'], '_FillValue'):
        data = np.where(data == ds.variables['TS']._FillValue, np.nan, data)
    data = np.where(np.abs(data) > 1e30, np.nan, data)
    ds.close()
    
    unique_years = np.unique(filtered_years)
    yearly_mean_data = np.empty((len(unique_years), data.shape[1], data.shape[2]))
    
    for i, yr in enumerate(unique_years):
        yr_mask = (filtered_years == yr)
        yearly_mean_data[i] = np.nanmean(data[yr_mask], axis=0)
        
    return yearly_mean_data

def get_model_yearly_all_layers(filename, time_indices, target_depths):
    ds = Dataset(filename)
    var = ds.variables['f_t_soisno']
    data = var[time_indices, :, :, 5:15]
    
    if hasattr(var, '_FillValue'):
        data = np.where(data == var._FillValue, np.nan, data)
    data = np.where(np.abs(data) > 1e30, np.nan, data)
    data = np.where(data == -9999, np.nan, data)
    data = data - 273.15
    
    data_yr = np.nanmean(data.reshape(nyear, 3, data.shape[1], data.shape[2], data.shape[3]), axis=1)
    
    ip_data = np.empty((nyear, data.shape[1], data.shape[2], len(target_depths)))
    for yr in range(nyear):
        ip_data[yr] = np.apply_along_axis(lambda x: np.interp(target_depths, colm_lev, x), axis=-1, arr=data_yr[yr])
    ds.close()
    return ip_data

def spatial_weighted_avg(data, mask, wgt):
    valid = mask & ~np.isnan(data) & ~np.isnan(wgt)
    if np.sum(valid) == 0: return np.nan
    return np.average(data[valid], weights=wgt[valid])

def get_min_max(data, mask): # === 新增：安全提取指定区域内最值的函数 ===
    valid_data = data[mask & ~np.isnan(data)]
    if len(valid_data) == 0: 
        return np.nan, np.nan
    return np.min(valid_data), np.max(valid_data)

# ====================== 4. 预先加载并处理模式数据 ======================
print("提取 17 年逐年模式结果并插值...")
mam_idx = [6 * i + j for i in range(nyear) for j in range(3)]
jja_idx = [6 * i + 3 + j for i in range(nyear) for j in range(3)]
target_depth_vals = [layer['colm_depth'] for layer in layers_config]

ctl_mam_yr = get_model_yearly_all_layers(sim_ctl_file, mam_idx, target_depth_vals)
ctl_jja_yr = get_model_yearly_all_layers(sim_ctl_file, jja_idx, target_depth_vals)
exp_mam_yr = get_model_yearly_all_layers(sim_exp_file, mam_idx, target_depth_vals)
exp_jja_yr = get_model_yearly_all_layers(sim_exp_file, jja_idx, target_depth_vals)

# ====================== 5. 大图绘制与终端值打印 ======================
print("\n================ 开始绘制极致紧凑的大面板 ================")
seasons = {'MAM': [3, 4, 5], 'JJA': [6, 7, 8]}
proj = ccrs.LambertConformal(central_longitude=110, central_latitude=40, standard_parallels=(30, 60))

fig, axes = plt.subplots(6, 3, figsize=(12, 22), subplot_kw={'projection': proj})

# 极致压缩 wspace 和 hspace
plt.subplots_adjust(wspace=0.03, hspace=0.03, left=0.08, right=0.98, bottom=0.08, top=0.95)

cf_obs = cf_bias = cf_diff = None

for s_idx, (season, months) in enumerate(seasons.items()):
    for l_idx, layer in enumerate(layers_config):
        row_idx = s_idx * 3 + l_idx
        layer_name = layer['name']
        print(f" -> 处理并检验 {season} {layer_name}")
        
        # 1. 提取 17 年观测
        obs_yr = extract_obs_yearly(obs_files[layer_name], months)
        obs_mean = np.nanmean(obs_yr, axis=0)
        obs_mean_masked = np.where(mask_target, obs_mean, np.nan)
        
        # 2. 提取并插值 17 年模式
        if season == 'MAM':
            ctl_layer_lambert = ctl_mam_yr[:, :, :, l_idx]
            exp_layer_lambert = exp_mam_yr[:, :, :, l_idx]
        else:
            ctl_layer_lambert = ctl_jja_yr[:, :, :, l_idx]
            exp_layer_lambert = exp_jja_yr[:, :, :, l_idx]
            
        ctl_yr_reg = np.empty((nyear, lat1d.size, lon1d.size))
        exp_yr_reg = np.empty((nyear, lat1d.size, lon1d.size))
        
        for yr in range(nyear):
            ctl_yr_reg[yr] = regrid_rcm2rgrid(ctl_layer_lambert[yr], lat2d, lon2d, lat1d, lon1d)
            exp_yr_reg[yr] = regrid_rcm2rgrid(exp_layer_lambert[yr], lat2d, lon2d, lat1d, lon1d)
            
        ctl_mean = np.nanmean(ctl_yr_reg, axis=0)
        exp_mean = np.nanmean(exp_yr_reg, axis=0)
        
        # 3. 统计学 T 检验 (p < 0.05)
        _, p_bias = ttest_rel(ctl_yr_reg, obs_yr, axis=0, nan_policy='omit')
        _, p_diff = ttest_rel(exp_yr_reg, ctl_yr_reg, axis=0, nan_policy='omit')
        
        # 4. 生成显著性填色 Mask (仅 p<0.05 保留)
        bias = ctl_mean - obs_mean
        diff = exp_mean - ctl_mean
        
        # ================== 新增：终端打印特定区域的均值与最值 ==================
        # 计算在砾石>=0.3且经度<=110°E区域的面积加权平均值
        sub_gravel_obs = spatial_weighted_avg(obs_mean, final_cond_mask, wgt_2d)
        sub_gravel_bias = spatial_weighted_avg(bias, final_cond_mask, wgt_2d)
        sub_gravel_diff = spatial_weighted_avg(diff, final_cond_mask, wgt_2d)
        
        # 计算区域内的最大、最小值
        obs_min, obs_max = get_min_max(obs_mean, final_cond_mask)
        bias_min, bias_max = get_min_max(bias, final_cond_mask)
        diff_min, diff_max = get_min_max(diff, final_cond_mask)
        
        print(f"    [Gravel>=0.3 & Lon<=110°E] {season} {layer_name} 区域特征:")
        print(f"        OBS     -> 加权平均: {sub_gravel_obs:.2f} °C | 最小值: {obs_min:.2f} °C | 最大值: {obs_max:.2f} °C")
        print(f"        CTL-OBS -> 加权平均: {sub_gravel_bias:.2f} °C | 最小值: {bias_min:.2f} °C | 最大值: {bias_max:.2f} °C")
        print(f"        EXP-CTL -> 加权平均: {sub_gravel_diff:.2f} °C | 最小值: {diff_min:.2f} °C | 最大值: {diff_max:.2f} °C")
        # ======================================================================

        bias_masked = np.where(mask_target & (p_bias < 0.05), bias, np.nan)
        diff_masked = np.where(mask_target & (p_diff < 0.05), diff, np.nan)
        
        # 5. 子图绘制
        ax_row = axes[row_idx]
        for col_idx, ax in enumerate(ax_row):
            ax.set_extent([80, 130, 10, 55], crs=ccrs.PlateCarree())
            for shp_file, lw in [("province.shp", 0.4), ("china.shp", 0.6), ("south_china_sea.shp", 0.8)]:
                p_shp = os.path.join(shp_dir, shp_file)
                if os.path.exists(p_shp):
                    color = 'blue' if 'river' in shp_file else 'black'
                    try:
                        reader = shpreader.Reader(p_shp)
                        ax.add_geometries(reader.geometries(), crs=ccrs.PlateCarree(), facecolor='none', edgecolor=color, linewidth=lw, zorder=3)
                    except Exception:
                        pass
            
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                              x_inline=False, y_inline=False,
                              linewidth=0.6, color='gray', alpha=0.5, linestyle='--', zorder=2)
            gl.xlocator = mticker.FixedLocator(np.arange(60, 140, 10))
            gl.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
            
            # 精准控制外围标签显示
            gl.top_labels = False
            gl.right_labels = False
            gl.left_labels = (col_idx == 0)   # 仅最左列保留纬度
            gl.bottom_labels = (row_idx == 5) # 仅最底行保留经度
            
            gl.xlabel_style = {'size': 10, 'rotation': 0, 'ha': 'center', 'va': 'top'}
            gl.ylabel_style = {'size': 10, 'rotation': 0, 'ha': 'right', 'va': 'center'}
            try:
                gl.padding = 6
            except Exception:
                pass

        if row_idx == 0:
            ax_row[0].set_title("OBS", loc='left', fontsize=11, pad=5)
            ax_row[1].set_title("CTL CPL - OBS", loc='left', fontsize=11, pad=5)
            ax_row[2].set_title("(EXP - CTL) CPL", loc='left', fontsize=11, pad=5)
            ax_row[2].set_title("°C", loc='right', fontsize=11, pad=5)
            
        season_str = 'Spring' if season == 'MAM' else 'Summer'
        ax_row[0].text(-0.20, 0.5, season_str, transform=ax_row[0].transAxes, fontsize=14, fontweight='normal', va='center', ha='center', rotation=90)

        # 填色范围设置
        levels_obs = np.linspace(-30, 30, 61)
        levels_bias = np.linspace(-6, 6, 25)
        levels_diff = np.linspace(-1.5, 1.5, 31)
        
        # 设置 antialiased=False
        kwargs = dict(cmap='RdBu_r', extend='both', transform=ccrs.PlateCarree(), zorder=1, antialiased=False)
        cf_obs = ax_row[0].contourf(grid_x, grid_y, obs_mean_masked, levels=levels_obs, **kwargs)
        cf_bias = ax_row[1].contourf(grid_x, grid_y, bias_masked, levels=levels_bias, **kwargs)
        cf_diff = ax_row[2].contourf(grid_x, grid_y, diff_masked, levels=levels_diff, **kwargs)

        for cf in [cf_obs, cf_bias, cf_diff]:
            try:
                for c in cf.collections:
                    c.set_edgecolor("face")
                    c.set_linewidth(1e-8)
            except AttributeError:
                cf.set_edgecolor("face")
                cf.set_linewidth(1e-8)

        for i in range(3):
            ax_row[i].text(0.03, 0.96, f"({chr(97 + row_idx * 3 + i)}) {layer_name}", 
                           transform=ax_row[i].transAxes, fontsize=11, va='top', ha='left', zorder=5, 
                           bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

# ================= 6. Colorbars 完美对齐 =================
plt.draw() 

pos0 = axes[5, 0].get_position()
pos1 = axes[5, 1].get_position()
pos2 = axes[5, 2].get_position()

cbar_y = pos0.y0 - 0.02
cbar_h = 0.008

cax_obs = fig.add_axes([pos0.x0, cbar_y, pos0.width, cbar_h])
plt.colorbar(cf_obs, cax=cax_obs, orientation='horizontal', ticks=np.arange(-30, 31, 10))

cax_bias = fig.add_axes([pos1.x0, cbar_y, pos1.width, cbar_h])
plt.colorbar(cf_bias, cax=cax_bias, orientation='horizontal', ticks=np.arange(-6, 7, 2))

cax_diff = fig.add_axes([pos2.x0, cbar_y, pos2.width, cbar_h])
plt.colorbar(cf_diff, cax=cax_diff, orientation='horizontal', ticks=np.arange(-1.5, 1.6, 0.5))

save_name = os.path.join(out_dir, "add_cpl_soilt_diff_raw.pdf")
plt.savefig(save_name, dpi=300, bbox_inches='tight')
print(f"  -> 已成功生成: {save_name}")
plt.close()