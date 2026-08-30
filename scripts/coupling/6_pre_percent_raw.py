import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import geopandas as gpd
from scipy.interpolate import griddata
from scipy.stats import ttest_rel
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import os
import warnings
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator

# ========================= 1. 环境与路径设置 =========================
warnings.filterwarnings('ignore')
os.environ['CARTOPY_OFFLINE'] = 'true' 
os.environ['PROJ_NETWORK'] = 'OFF'

path_data = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
path_ref = "/share/home/dq135/reference/"
shp_dir = "../../shapefile_China/"
n_years = 17 

# 投影参数
proj = ccrs.LambertConformal(
    central_longitude=110, 
    central_latitude=40,   
    standard_parallels=(30, 60)
)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ========================= 2. 数据处理函数 =========================

def regrid_rcm2rgrid(var2d, lat2d, lon2d, lon_mesh, lat_mesh):
    pts = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    vals = var2d.ravel()
    interp = griddata(pts, vals, (lon_mesh, lat_mesh), method='linear')
    nan_mask = np.isnan(interp)
    if nan_mask.any():
        interp[nan_mask] = griddata(pts, vals, (lon_mesh, lat_mesh), method='nearest')[nan_mask]
    return interp

def get_annual_data(season, case, is_obs=False):
    if is_obs:
        obs_file = f"{path_ref}CN05.1_Pre_1991_2023_{season}_025x025.nc"
        ds = xr.open_dataset(obs_file)
        var = ds['pre'].isel(time=slice(30, 30 + n_years * 3))
        annual = var.values.reshape(n_years, 3, var.shape[1], var.shape[2]).mean(axis=1)
        lon, lat = ds.lon.values, ds.lat.values
        return lon, lat, annual
    else:
        sim_file = f"{path_data}wrfout_2001-2017_{season}_{case}.nc"
        ds = xr.open_dataset(sim_file)
        s_time_dim = 'Times' if 'Times' in ds.dims else 'Time'
        raw = ds['PRAVG'].isel({s_time_dim: slice(0, n_years * 92)})
        annual = raw.values.reshape(n_years, 92, raw.shape[1], raw.shape[2]).mean(axis=1) * 86400.0
        return annual

f_in_path = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
f_in = xr.open_dataset(f_in_path)
i_time_dim = 'Time' if 'Time' in f_in['XLAT'].dims else 'Times'
lat2d = f_in['XLAT'].isel({i_time_dim: 0}).values
lon2d = f_in['XLONG'].isel({i_time_dim: 0}).values

# ========================= 2.5 提取高分辨率砾石数据及辅助函数 =========================
print("提取高分辨率砾石数据，并生成条件掩膜...")
lon_obs_temp, lat_obs_temp, _ = get_annual_data("MAM", None, is_obs=True)
lon_mesh_temp, lat_mesh_temp = np.meshgrid(lon_obs_temp, lat_obs_temp)
wgt_2d = np.cos(np.radians(lat_mesh_temp))

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
    valid_mask_g = ~np.isnan(data)
    sum_gravel_weighted[valid_mask_g] += data[valid_mask_g] * weight
    sum_weights[valid_mask_g] += weight
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
pts = np.stack((lat_mesh_temp, lon_mesh_temp), axis=-1)
mean_gravel_025 = interp_func(pts)

def get_min_max(data, mask):
    valid_data = data[mask & ~np.isnan(data)]
    if len(valid_data) == 0: 
        return np.nan, np.nan
    return np.min(valid_data), np.max(valid_data)

def spatial_weighted_avg(data, mask, wgt):
    valid = mask & ~np.isnan(data) & ~np.isnan(wgt)
    if np.sum(valid) == 0: return np.nan
    return np.average(data[valid], weights=wgt[valid])

# ========================= 3. 地图样式 =========================

def apply_reference_style(ax, letter, title_str="", is_left=False, is_bottom=False):
    ax.set_extent([80, 130, 10, 55], crs=ccrs.PlateCarree())
    
    def add_local_shp(ax, name, lw=0.6, color='black', zorder=3):
        path = os.path.join(shp_dir, name)
        if os.path.exists(path):
            try:
                gdf = gpd.read_file(path)
                if gdf.crs is None:
                    gdf = gdf.set_crs("EPSG:4326")
                gdf.to_crs(ax.projection).plot(ax=ax, facecolor='none', edgecolor=color, linewidth=lw, zorder=zorder)
            except:
                pass

    add_local_shp(ax, "china.shp", lw=1.2, color='black', zorder=5) 
    add_local_shp(ax, "province.shp", lw=0.5, color='black', zorder=4) 
    add_local_shp(ax, "south_china_sea.shp", lw=1.0, color='black', zorder=5) 

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False,
                      linewidth=0.8, color='grey', alpha=0.5, linestyle='--', zorder=2)
    gl.xlocator = mticker.FixedLocator(np.arange(60, 140, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
    
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = is_left      
    gl.bottom_labels = is_bottom  
    
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    gl.xlabel_style = {'size': 11, 'rotation': 0, 'va': 'top', 'ha': 'center'}
    gl.ylabel_style = {'size': 11, 'rotation': 0, 'va': 'center', 'ha': 'right'}

    if title_str:
       ax.set_title(title_str, loc='left', fontsize=12, fontweight='normal', pad=5)

    ax.text(0.02, 0.96, f"({letter})", transform=ax.transAxes, fontsize=12, 
            va='top', ha='left', zorder=10, 
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

# ========================= 4. 图表核心组装 =========================

seasons_code = ["MAM", "JJA"]
seasons_name = ["Spring", "Summer"]
letters = [['a', 'b', 'c', 'd'], ['e', 'f', 'g', 'h']]

levels_obs = np.arange(0, 11, 1)  
cmap_obs = plt.get_cmap('GnBu') 

levels_bias = np.arange(-100, 101, 20)  
cmap_bias = plt.get_cmap('RdBu')      

levels_diff = np.arange(-20, 21, 4)    
cmap_diff = plt.get_cmap('RdBu')        

alpha_val = 0.05

fig, axes = plt.subplots(2, 4, figsize=(16.5, 8.5), subplot_kw={'projection': proj})
plt.subplots_adjust(wspace=0.01, hspace=0.01, left=0.05, right=0.98, bottom=0.15, top=0.95)

cf0 = cf1 = cf3 = None

for r, season in enumerate(seasons_code):
    s_name = seasons_name[r]
    print(f"=========================================")
    print(f"正在处理 {season} 季节数据...")
    
    lon_obs, lat_obs, obs_annual = get_annual_data(season, None, is_obs=True)
    lon_mesh, lat_mesh = np.meshgrid(lon_obs, lat_obs)
    obs_mean = obs_annual.mean(axis=0)
    
    ctl_annual_raw = get_annual_data(season, "nogravel")
    exp_annual_raw = get_annual_data(season, "gravel")
    
    ctl_annual = np.zeros_like(obs_annual)
    exp_annual = np.zeros_like(obs_annual)
    for y in range(n_years):
        ctl_annual[y] = regrid_rcm2rgrid(ctl_annual_raw[y], lat2d, lon2d, lon_mesh, lat_mesh)
        exp_annual[y] = regrid_rcm2rgrid(exp_annual_raw[y], lat2d, lon2d, lon_mesh, lat_mesh)
        
    ctl_mean = ctl_annual.mean(axis=0)
    exp_mean = exp_annual.mean(axis=0)
    
    _, p_val = ttest_rel(exp_annual, ctl_annual, axis=0)

    with np.errstate(divide='ignore', invalid='ignore'):
        ctl_bias_raw = (ctl_mean - obs_mean) / np.where(obs_mean > 0.01, obs_mean, np.nan) * 100.0
        exp_bias_raw = (exp_mean - obs_mean) / np.where(obs_mean > 0.01, obs_mean, np.nan) * 100.0
    
    diff_mean = exp_mean - ctl_mean            # 绝对降水差值
    diff_bias_raw = exp_bias_raw - ctl_bias_raw # 百分比降水差值

    # 画图用的变量，继续应用截断防止 colorbar 爆炸
    ctl_bias = np.clip(ctl_bias_raw, -100, 100)
    exp_bias = np.clip(exp_bias_raw, -100, 100)
    diff_bias = exp_bias - ctl_bias
    
    valid_mask = ~np.isnan(obs_mean)
    
    # 用有效掩膜过滤画图数据
    ctl_bias = np.where(valid_mask, ctl_bias, np.nan)
    exp_bias = np.where(valid_mask, exp_bias, np.nan)
    diff_bias = np.where(valid_mask, diff_bias, np.nan)
    
    # 用有效掩膜过滤 raw 统计数据
    ctl_bias_raw = np.where(valid_mask, ctl_bias_raw, np.nan)
    exp_bias_raw = np.where(valid_mask, exp_bias_raw, np.nan)
    diff_bias_raw = np.where(valid_mask, diff_bias_raw, np.nan)

    # ================= 【修改】添加剔除特定区域的逻辑与均值计算 =================
    exclude_mask = (lat_mesh >= 25.0) & (lat_mesh <= 29.5) & (lon_mesh >= 91.0) & (lon_mesh <= 98.0)
    
    mask_all_calc = valid_mask & (~exclude_mask)
    mask_gravel_calc = (mean_gravel_025 >= 0.3) & (lon_mesh < 110.0) & valid_mask & (~exclude_mask)
    mask_sc_calc = (lon_mesh >= 105.0) & (lon_mesh <= 115.0) & (lat_mesh >= 20.0) & (lat_mesh <= 25.0) & valid_mask & (~exclude_mask)
    
    obs_min_all, obs_max_all = get_min_max(obs_mean, mask_all_calc)
    obs_min_g, obs_max_g = get_min_max(obs_mean, mask_gravel_calc)
    obs_mean_g = spatial_weighted_avg(obs_mean, mask_gravel_calc, wgt_2d)
    obs_min_sc, obs_max_sc = get_min_max(obs_mean, mask_sc_calc)
    obs_mean_sc = spatial_weighted_avg(obs_mean, mask_sc_calc, wgt_2d)
    
    ctl_min_g, ctl_max_g = get_min_max(ctl_mean, mask_gravel_calc)
    ctl_mean_g = spatial_weighted_avg(ctl_mean, mask_gravel_calc, wgt_2d)
    ctl_min_sc, ctl_max_sc = get_min_max(ctl_mean, mask_sc_calc)
    ctl_mean_sc = spatial_weighted_avg(ctl_mean, mask_sc_calc, wgt_2d)
    
    ctl_b_min_g, ctl_b_max_g = get_min_max(ctl_bias_raw, mask_gravel_calc)
    ctl_b_mean_g = spatial_weighted_avg(ctl_bias_raw, mask_gravel_calc, wgt_2d)
    ctl_b_min_sc, ctl_b_max_sc = get_min_max(ctl_bias_raw, mask_sc_calc)
    ctl_b_mean_sc = spatial_weighted_avg(ctl_bias_raw, mask_sc_calc, wgt_2d)
    
    exp_min_g, exp_max_g = get_min_max(exp_mean, mask_gravel_calc)
    exp_mean_g = spatial_weighted_avg(exp_mean, mask_gravel_calc, wgt_2d)
    exp_min_sc, exp_max_sc = get_min_max(exp_mean, mask_sc_calc)
    exp_mean_sc = spatial_weighted_avg(exp_mean, mask_sc_calc, wgt_2d)
    
    exp_b_min_g, exp_b_max_g = get_min_max(exp_bias_raw, mask_gravel_calc)
    exp_b_mean_g = spatial_weighted_avg(exp_bias_raw, mask_gravel_calc, wgt_2d)
    exp_b_min_sc, exp_b_max_sc = get_min_max(exp_bias_raw, mask_sc_calc)
    exp_b_mean_sc = spatial_weighted_avg(exp_bias_raw, mask_sc_calc, wgt_2d)
    
    diff_min_g, diff_max_g = get_min_max(diff_mean, mask_gravel_calc)
    diff_mean_g = spatial_weighted_avg(diff_mean, mask_gravel_calc, wgt_2d)
    diff_min_sc, diff_max_sc = get_min_max(diff_mean, mask_sc_calc)
    diff_mean_sc = spatial_weighted_avg(diff_mean, mask_sc_calc, wgt_2d)
    
    diff_b_min_g, diff_b_max_g = get_min_max(diff_bias_raw, mask_gravel_calc)
    diff_b_mean_g = spatial_weighted_avg(diff_bias_raw, mask_gravel_calc, wgt_2d)
    diff_b_min_sc, diff_b_max_sc = get_min_max(diff_bias_raw, mask_sc_calc)
    diff_b_mean_sc = spatial_weighted_avg(diff_bias_raw, mask_sc_calc, wgt_2d)

    # ================= 提取显著性区域散点坐标 =================
    sig_mask = (p_val < alpha_val) & valid_mask
    
    # 【修改 1】将步长从 6 提高到 8，稍微降低打点的密度
    step = 10  
    lon_sig = lon_mesh[sig_mask][::step]
    lat_sig = lat_mesh[sig_mask][::step]

    is_bottom_row = (r == 1)

    # ================= 第一列：OBS =================
    ax_obs = axes[r, 0]
    cf0 = ax_obs.contourf(lon_mesh, lat_mesh, obs_mean, levels=levels_obs, cmap=cmap_obs, extend='max', transform=ccrs.PlateCarree(), zorder=1)
    title_obs = "OBS" if r == 0 else ""
    apply_reference_style(ax_obs, letters[r][0], title_obs, is_left=True, is_bottom=is_bottom_row)
    
    ax_obs.text(-0.16, 0.5, s_name, va='center', ha='center', rotation=90, transform=ax_obs.transAxes, 
                fontsize=12, fontweight='normal', color='black')

    # ================= 第二列：CTL CPL =================
    ax_ctl = axes[r, 1]
    cf1 = ax_ctl.contourf(lon_mesh, lat_mesh, ctl_bias, levels=levels_bias, cmap=cmap_bias, extend='both', transform=ccrs.PlateCarree(), zorder=1)
    title_ctl = "CTL CPL" if r == 0 else ""
    apply_reference_style(ax_ctl, letters[r][1], title_ctl, is_left=False, is_bottom=is_bottom_row)
    
    # ================= 第三列：EXP CPL =================
    ax_exp = axes[r, 2]
    cf2 = ax_exp.contourf(lon_mesh, lat_mesh, exp_bias, levels=levels_bias, cmap=cmap_bias, extend='both', transform=ccrs.PlateCarree(), zorder=1)
    title_exp = "EXP CPL" if r == 0 else ""
    apply_reference_style(ax_exp, letters[r][2], title_exp, is_left=False, is_bottom=is_bottom_row)

    # ================= 第四列：(EXP - CTL) =================
    ax_comp = axes[r, 3]
    cf3 = ax_comp.contourf(lon_mesh, lat_mesh, diff_bias, levels=levels_diff, cmap=cmap_diff, extend='both', transform=ccrs.PlateCarree(), zorder=1)
    
    # ----------------------------------------------------
    # 【修改 2】增加星号透明度：将 alpha 从 1.0 降低为 0.7
    # ----------------------------------------------------
    if len(lon_sig) > 0:
         ax_comp.scatter(lon_sig, lat_sig, marker='x', facecolors='black', edgecolors='black', 
                         s=15, linewidths=1.3, alpha=0.7, transform=ccrs.PlateCarree(), zorder=6)

    title_comp = "(EXP - CTL) CPL" if r == 0 else ""
    apply_reference_style(ax_comp, letters[r][3], title_comp, is_left=False, is_bottom=is_bottom_row)

# ========================= 5. 底部色标动态布局 =========================

fig.canvas.draw()

pos_ax0 = axes[1, 0].get_position()
pos_ax1 = axes[1, 1].get_position()
pos_ax2 = axes[1, 2].get_position()
pos_ax3 = axes[1, 3].get_position()

cb_height = 0.015
cb_y = pos_ax0.y0 - 0.075

# Colorbar 0：观测绝对值 (第1列)
cb0_x = pos_ax0.x0
cb0_w = pos_ax0.x1 - pos_ax0.x0
cbar_ax_obs = fig.add_axes([cb0_x, cb_y, cb0_w, cb_height]) 
cb_obs = plt.colorbar(cf0, cax=cbar_ax_obs, orientation='horizontal', spacing='uniform', extend='max')
cb_obs.set_ticks(levels_obs)
cb_obs.ax.tick_params(labelsize=11)
cb_obs.ax.set_title("mm/day", fontsize=11, pad=4)

# Colorbar 1：相对偏差 (第2、3列)
cb1_x = pos_ax1.x0
cb1_w = pos_ax2.x1 - pos_ax1.x0
cbar_ax_bias = fig.add_axes([cb1_x, cb_y, cb1_w, cb_height]) 
cb_bias = plt.colorbar(cf1, cax=cbar_ax_bias, orientation='horizontal', spacing='uniform', extend='both')
cb_bias.set_ticks(levels_bias)
cb_bias.ax.tick_params(labelsize=11)
cb_bias.ax.set_title("%", fontsize=11, pad=4)

# Colorbar 2：差值 (第4列)
cb2_x = pos_ax3.x0
cb2_w = pos_ax3.x1 - pos_ax3.x0
cbar_ax_diff = fig.add_axes([cb2_x, cb_y, cb2_w, cb_height]) 
cb_diff = plt.colorbar(cf3, cax=cbar_ax_diff, orientation='horizontal', spacing='uniform', extend='both')
cb_diff.set_ticks(levels_diff)
cb_diff.ax.tick_params(labelsize=11)
cb_diff.ax.set_title("%", fontsize=11, pad=4)

# ========================= 6. 保存 =========================
save_path = "../illustration/cpl/FIG6.precip_bias_percent_raw.png"
os.makedirs(os.path.dirname(save_path), exist_ok=True)
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"绘图完成，合并图像已保存为：{save_path}")