import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from scipy import stats
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
from netCDF4 import Dataset
import os
import warnings

warnings.filterwarnings("ignore")

# ==============================================================================
# 0. 自定义绘图配置
# ==============================================================================
RIDGE_SOURCE = 'CTL' 
HGT200_LEVEL = 12500 if RIDGE_SOURCE == 'ERA5' else 12550
HGT500_LEVEL = 5860  if RIDGE_SOURCE == 'ERA5' else 5870  

# ==============================================================================
# 1. 全局字体与环境设置
# ==============================================================================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

os.environ['CARTOPY_OFFLINE'] = 'true'
os.environ['PROJ_NETWORK'] = 'OFF'

# ==============================================================================
# 2. 基础设置与文件路径 
# ==============================================================================
path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
shp_dir = "../../shapefile_China/"
season = "JJA" 
nyears = 17

files = {
    'ctl_200_hgt': path + f"wrfpost_gh200_2001_2017_{season}_nogravel.nc",
    'exp_200_hgt': path + f"wrfpost_gh200_2001_2017_{season}_gravel.nc",
    'ctl_200_u': path + f"wrfpost_u200_2001_2017_{season}_nogravel.nc",
    'exp_200_u': path + f"wrfpost_u200_2001_2017_{season}_gravel.nc",
    'ctl_200_v': path + f"wrfpost_v200_2001_2017_{season}_nogravel.nc",
    'exp_200_v': path + f"wrfpost_v200_2001_2017_{season}_gravel.nc",
    'ctl_500_hgt': path + f"wrfpost_gh500_2001_2017_{season}_nogravel.nc",
    'exp_500_hgt': path + f"wrfpost_gh500_2001_2017_{season}_gravel.nc",
    'ctl_850_q': path + f"wrfpost_qx850_2001_2017_{season}_nogravel.nc",
    'exp_850_q': path + f"wrfpost_qx850_2001_2017_{season}_gravel.nc",
    'ctl_850_hgt': path + f"wrfpost_gh850_2001_2017_{season}_nogravel.nc", 
    'exp_850_hgt': path + f"wrfpost_gh850_2001_2017_{season}_gravel.nc",
}

geo_file = "/share/home/dq135/wrfpost/wrfinput_d01"

era5_file = "/share/home/dq135/reference/era5_1991_2023_monthly_025x025.nc"
sim_ctl_file = path + "colmoff_2001-2023_monmean_nogravel.nc"
sim_exp_file = path + "colmoff_2001-2023_monmean_gravel.nc"

# ==============================================================================
# 3. 数据处理与计算函数
# ==============================================================================
def get_seasonal_mean_per_year(ds, var_name, is_geopt=False, scale=1.0):
    data = ds[var_name][:, 0, :, :].values if len(ds[var_name].shape) == 4 else ds[var_name].values
    if is_geopt:
        data = data / 9.8
    data = data * scale
    years, days_per_season = nyears, 92
    lat, lon = data.shape[1], data.shape[2]
    data_reshaped = data[:years*days_per_season].reshape(years, days_per_season, lat, lon)
    return np.mean(data_reshaped, axis=1)

def calc_diff_and_ttest(ctl_data, exp_data):
    ctl_mean = np.nanmean(ctl_data, axis=0)
    exp_mean = np.nanmean(exp_data, axis=0)
    diff = exp_mean - ctl_mean
    t_stat, p_value = stats.ttest_rel(exp_data, ctl_data, axis=0, nan_policy='omit')
    return ctl_mean, exp_mean, diff, p_value

def calc_area_mean_idx(data_yr, mask):
    """提取特定掩膜内的区域平均并进行标准化计算异常"""
    idx = np.array([np.nanmean(year_data[mask]) for year_data in data_yr])
    return (idx - np.mean(idx)) / np.std(idx)

def calc_correlation_and_ttest(idx, data):
    """计算相关系数场并返回真实的 P 值"""
    idx_anom = idx - np.nanmean(idx)
    data_anom = data - np.nanmean(data, axis=0)
    
    cov = np.nanmean(idx_anom[:, None, None] * data_anom, axis=0)
    std_idx = np.nanstd(idx_anom)
    std_data = np.nanstd(data_anom, axis=0)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        corr = cov / (std_idx * std_data)
        n = len(idx)
        corr_clipped = np.clip(corr, -0.99999, 0.99999) 
        t_stat = corr_clipped * np.sqrt(n - 2) / np.sqrt(1.0 - corr_clipped**2)
        p_val = stats.t.sf(np.abs(t_stat), n - 2) * 2
        
    return corr, p_val

def get_colm_soil_temp_L1L2_yearly(filename, time_indices, nyear=17):
    """提取 CoLM 第一、二层土壤温度求平均 (切片 5:7)"""
    ds = Dataset(filename)
    var = ds.variables['f_t_soisno']
    data = var[time_indices, :, :, 5:7] 
    
    if hasattr(var, '_FillValue'):
        data = np.where(data == var._FillValue, np.nan, data)
    data = np.where(np.abs(data) > 1e30, np.nan, data)
    data = np.where(data == -9999, np.nan, data)
    data = data - 273.15
    
    data_mean_layers = np.nanmean(data, axis=3)
    data_yr = np.nanmean(data_mean_layers.reshape(nyear, 3, data.shape[1], data.shape[2]), axis=1)
    ds.close()
    return data_yr

print("开始加载并计算 WRF 与相关数据...")
ds_geo = xr.open_dataset(geo_file)
lons = ds_geo['XLONG'][0].values
lats = ds_geo['XLAT'][0].values
cosalpha = ds_geo['COSALPHA'][0].values
sinalpha = ds_geo['SINALPHA'][0].values

# ----------- [步骤 A] 计算原图气候态变量 (a, b, c) -----------
ctl_200_u_yr = get_seasonal_mean_per_year(xr.open_dataset(files['ctl_200_u']), 'U')
exp_200_u_yr = get_seasonal_mean_per_year(xr.open_dataset(files['exp_200_u']), 'U')
ctl_200_v_yr = get_seasonal_mean_per_year(xr.open_dataset(files['ctl_200_v']), 'V')
exp_200_v_yr = get_seasonal_mean_per_year(xr.open_dataset(files['exp_200_v']), 'V')

exp_200_u_true_yr = exp_200_u_yr * cosalpha - exp_200_v_yr * sinalpha
ctl_200_u_true_yr = ctl_200_u_yr * cosalpha - ctl_200_v_yr * sinalpha

u200_c_grid, u200_e_grid, _, u200_pval = calc_diff_and_ttest(ctl_200_u_yr, exp_200_u_yr)
v200_c_grid, v200_e_grid, _, v200_pval = calc_diff_and_ttest(ctl_200_v_yr, exp_200_v_yr)

u200_c = u200_c_grid * cosalpha - v200_c_grid * sinalpha
u200_e = u200_e_grid * cosalpha - v200_e_grid * sinalpha
u200_d = u200_e - u200_c

ctl_200_z_yr = get_seasonal_mean_per_year(xr.open_dataset(files['ctl_200_hgt']), 'geopt', is_geopt=True)
exp_200_z_yr = get_seasonal_mean_per_year(xr.open_dataset(files['exp_200_hgt']), 'geopt', is_geopt=True)
z200_c, z200_e, z200_d, z200_pval = calc_diff_and_ttest(ctl_200_z_yr, exp_200_z_yr)
z200_c_smooth = gaussian_filter(z200_c, sigma=3.0)

ctl_500_z_yr = get_seasonal_mean_per_year(xr.open_dataset(files['ctl_500_hgt']), 'geopt', is_geopt=True)
exp_500_z_yr = get_seasonal_mean_per_year(xr.open_dataset(files['exp_500_hgt']), 'geopt', is_geopt=True)
z500_c, z500_e, z500_d, z500_pval = calc_diff_and_ttest(ctl_500_z_yr, exp_500_z_yr)

ctl_850_uq_yr = get_seasonal_mean_per_year(xr.open_dataset(files['ctl_850_q']), 'uq', scale=1000.0)
exp_850_uq_yr = get_seasonal_mean_per_year(xr.open_dataset(files['exp_850_q']), 'uq', scale=1000.0)
uq850_c, _, uq850_d, uq850_pval = calc_diff_and_ttest(ctl_850_uq_yr, exp_850_uq_yr)

ctl_850_vq_yr = get_seasonal_mean_per_year(xr.open_dataset(files['ctl_850_q']), 'vq', scale=1000.0)
exp_850_vq_yr = get_seasonal_mean_per_year(xr.open_dataset(files['exp_850_q']), 'vq', scale=1000.0)
vq850_c, _, vq850_d, vq850_pval = calc_diff_and_ttest(ctl_850_vq_yr, exp_850_vq_yr)

ctl_850_z_yr = get_seasonal_mean_per_year(xr.open_dataset(files['ctl_850_hgt']), 'geopt', is_geopt=True)
exp_850_z_yr = get_seasonal_mean_per_year(xr.open_dataset(files['exp_850_hgt']), 'geopt', is_geopt=True)
_, _, z850_d, z850_pval = calc_diff_and_ttest(ctl_850_z_yr, exp_850_z_yr)


# ----------- [步骤 B] 提取 ERA5 真实观测及计算模型土壤温度差值 -----------
print("开始加载 ERA5 再分析数据进行统计关系计算...")
ds_era5 = xr.open_dataset(era5_file)

time_mask = (ds_era5['valid_time'].dt.year >= 2001) & \
            (ds_era5['valid_time'].dt.year <= 2017) & \
            (ds_era5['valid_time'].dt.month.isin([6, 7, 8]))
ds_era5_jja = ds_era5.isel(valid_time=time_mask)

ds_era5_jja.coords['year'] = ds_era5_jja['valid_time'].dt.year
ds_era5_yr = ds_era5_jja.groupby('year').mean(dim='valid_time')

z200_era5_yr = ds_era5_yr['z'].sel(pressure_level=200).values / 9.80665
u200_era5_yr = ds_era5_yr['u'].sel(pressure_level=200).values
z500_era5_yr = ds_era5_yr['z'].sel(pressure_level=500).values / 9.80665
v850_era5_yr = ds_era5_yr['v'].sel(pressure_level=850).values

lons_era5 = ds_era5_yr['longitude'].values
lats_era5 = ds_era5_yr['latitude'].values
lon_era5_g, lat_era5_g = np.meshgrid(lons_era5, lats_era5)

era5_mean_clim = ds_era5_yr.mean(dim='year')
z200_era5 = era5_mean_clim['z'].sel(pressure_level=200).values / 9.80665
z500_era5 = era5_mean_clim['z'].sel(pressure_level=500).values / 9.80665

# 提取离线土壤温度
jja_indices = [6 * i + 3 + j for i in range(nyears) for j in range(3)]
ctl_st_raw_yr = get_colm_soil_temp_L1L2_yearly(sim_ctl_file, jja_indices, nyear=nyears)
exp_st_raw_yr = get_colm_soil_temp_L1L2_yearly(sim_exp_file, jja_indices, nyear=nyears)


# ----------- [步骤 C] 提取 [EXP-CTL] 异常指数并生成相关场 -----------

# (d) TP_ST_idx: 29-33N, 80-92E (Exp-Ctl Soil Temp)
tp_mask = (lats >= 29.0) & (lats <= 33.0) & (lons >= 80.0) & (lons <= 92.0)
diff_st_raw_yr = exp_st_raw_yr - ctl_st_raw_yr
tp_st_idx = calc_area_mean_idx(diff_st_raw_yr, tp_mask)
u200_tp_corr, u200_tp_pval = calc_correlation_and_ttest(tp_st_idx, u200_era5_yr)

# (e) EAJI_idx: 30-40N, 80-120E (Exp-Ctl U200)
eaji_mask_wrf = (lats >= 30.0) & (lats <= 40.0) & (lons >= 80.0) & (lons <= 120.0)
u200_diff_yr = exp_200_u_true_yr - ctl_200_u_true_yr
eaji_idx = calc_area_mean_idx(u200_diff_yr, eaji_mask_wrf)
h500_eaji_corr, h500_eaji_pval = calc_correlation_and_ttest(eaji_idx, z500_era5_yr)

# ⭐ 修改处: (f) H500_idx: 15-30N, 110-145E (Exp-Ctl Z500)
h500_mask_wrf = (lats >= 15.0) & (lats <= 30.0) & (lons >= 110.0) & (lons <= 145.0)
z500_diff_yr = exp_500_z_yr - ctl_500_z_yr
h500_idx = calc_area_mean_idx(z500_diff_yr, h500_mask_wrf)
v850_h500_corr, v850_h500_pval = calc_correlation_and_ttest(h500_idx, v850_era5_yr)

# ----------- [步骤 D] 彻底修复四角的插值重构函数 -----------
def regrid_era5_to_wrf(era5_data_2d):
    """将基于矩形的 ERA5 投影转换并填满曲线网格的 WRF lons/lats 空间"""
    pts = np.column_stack((lon_era5_g.ravel(), lat_era5_g.ravel()))
    vals = era5_data_2d.ravel()
    
    interp = griddata(pts, vals, (lons, lats), method='linear')
    nan_mask = np.isnan(interp)
    # 利用 nearest 填补四个边角外的空缺
    if nan_mask.any():
        interp[nan_mask] = griddata(pts, vals, (lons[nan_mask], lats[nan_mask]), method='nearest')
    return interp

print("插值重构：解决地图边缘填充问题...")
# 对 相关系数 及 P值 进行统一网格对齐
u200_tp_corr_wrf = regrid_era5_to_wrf(u200_tp_corr)
u200_tp_pval_wrf = regrid_era5_to_wrf(u200_tp_pval)

h500_eaji_corr_wrf = regrid_era5_to_wrf(h500_eaji_corr)
h500_eaji_pval_wrf = regrid_era5_to_wrf(h500_eaji_pval)

v850_h500_corr_wrf = regrid_era5_to_wrf(v850_h500_corr)
v850_h500_pval_wrf = regrid_era5_to_wrf(v850_h500_pval)


# ==============================================================================
# 4. 绘图底图及双重打点函数
# ==============================================================================
proj = ccrs.LambertConformal(central_longitude=110, central_latitude=40, standard_parallels=(30, 60))
data_crs = ccrs.PlateCarree()

def apply_style_and_shp(ax):
    ax.set_extent([83, 130, 15, 55], crs=data_crs)
    for shp_file, lw in [("province.shp", 0.4), ("china.shp", 0.6), ("south_china_sea.shp", 0.8)]:
        p_shp = os.path.join(shp_dir, shp_file)
        if os.path.exists(p_shp):
            color = 'gray' if 'province' in shp_file else 'black'
            try:
                reader = shpreader.Reader(p_shp)
                ax.add_geometries(reader.geometries(), crs=data_crs, facecolor='none', edgecolor=color, linewidth=lw, zorder=3)
            except: pass
            
    gl = ax.gridlines(crs=data_crs, draw_labels=True, x_inline=False, y_inline=False, linewidth=0.6, color='gray', alpha=0.5, linestyle='--', zorder=2)
    gl.xlocator = mticker.FixedLocator(np.arange(70, 135, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(10, 55, 10))
    gl.top_labels = gl.right_labels = False
    gl.bottom_labels = gl.left_labels = True
    
    # 放大字体，拉开边界
    gl.xlabel_style = {'size': 11, 'rotation': 0, 'ha': 'center', 'va': 'top'}
    gl.ylabel_style = {'size': 11, 'rotation': 0, 'ha': 'right', 'va': 'center'}
    gl.padding = 10 

def add_stippling_wrf(ax, p_val):
    """统一采用 WRF 坐标系统的双重显著性打点"""
    stride = 4
    lons_sub = lons[::stride, ::stride]
    lats_sub = lats[::stride, ::stride]
    pval_sub = p_val[::stride, ::stride]
    
    # 0.05 <= p < 0.2: 灰色实心圆点
    mask_02 = (pval_sub < 0.2) & (pval_sub >= 0.05) & ~np.isnan(pval_sub)
    if np.any(mask_02):
        ax.scatter(lons_sub[mask_02], lats_sub[mask_02], s=6, color='gray', marker='o', edgecolors='none', alpha=0.9, transform=data_crs, zorder=5)
        
    # p < 0.05: 白色实心圆点
    mask_05 = (pval_sub < 0.05) & ~np.isnan(pval_sub)
    if np.any(mask_05):
        ax.scatter(lons_sub[mask_05], lats_sub[mask_05], s=6, color='white', marker='o', edgecolors='none', alpha=0.9, transform=data_crs, zorder=6)


# ==============================================================================
# 5. 开始多图版绘制 (2x3 排版布局)
# ==============================================================================
fig = plt.figure(figsize=(18, 12)) 
plt.subplots_adjust(wspace=0.12, hspace=0.15, bottom=0.10) 

lvl_corr = np.arange(-0.7, 0.75, 0.1)

# ----------------- 上排: (a, b, c) 原气候态偏差 -----------------
ax1 = fig.add_subplot(2, 3, 1, projection=proj)
apply_style_and_shp(ax1)
ax1.set_title("(a) 200 hPa HGT & Westerly Jet Difference", loc='left', fontsize=11, pad=5)
cf1 = ax1.contourf(lons, lats, z200_d, levels=np.arange(-3.2, 3.3, 0.2), cmap='RdBu_r', extend='both', transform=data_crs, zorder=1)
add_stippling_wrf(ax1, u200_pval)

if RIDGE_SOURCE == 'ERA5':
    ax1.contour(lon_era5_g, lat_era5_g, z200_era5, levels=[HGT200_LEVEL], colors='darkgreen', linewidths=2.0, transform=data_crs, zorder=5)
else:
    ax1.contour(lons, lats, z200_c_smooth, levels=[HGT200_LEVEL], colors='darkgreen', linewidths=2.0, transform=data_crs, zorder=5)

u200_d_smooth = gaussian_filter(u200_d, sigma=3.0)
ax1.contour(lons, lats, u200_d_smooth, levels=np.arange(0.1, 0.4, 0.1), colors='tomato', linewidths=1.8, linestyles='solid', transform=data_crs, zorder=6)
ax1.contour(lons, lats, u200_d_smooth, levels=np.arange(-0.3, 0.0, 0.1), colors='lightskyblue', linewidths=1.8, linestyles='dashed', transform=data_crs, zorder=6)
plt.colorbar(cf1, ax=ax1, orientation='horizontal', fraction=0.046, pad=0.08)

ax2 = fig.add_subplot(2, 3, 2, projection=proj)
apply_style_and_shp(ax2)
ax2.set_title("(b) 500 hPa HGT Difference", loc='left', fontsize=11, pad=5)
cf2 = ax2.contourf(lons, lats, z500_d, levels=np.arange(-1.8, 1.9, 0.1), cmap='RdBu_r', extend='both', transform=data_crs, zorder=1)
add_stippling_wrf(ax2, z500_pval)

if RIDGE_SOURCE == 'ERA5':
    ax2.contour(lon_era5_g, lat_era5_g, z500_era5, levels=[HGT500_LEVEL], colors='darkgreen', linewidths=2.0, transform=data_crs, zorder=5)
else:
    ax2.contour(lons, lats, gaussian_filter(z500_c, 3.0), levels=[HGT500_LEVEL], colors='darkgreen', linewidths=2.0, transform=data_crs, zorder=5)
plt.colorbar(cf2, ax=ax2, orientation='horizontal', fraction=0.046, pad=0.08)

ax3 = fig.add_subplot(2, 3, 3, projection=proj)
apply_style_and_shp(ax3)
ax3.set_title("(c) 850 hPa HGT & Moisture Flux Diff", loc='left', fontsize=11, pad=5)
cf3 = ax3.contourf(lons, lats, z850_d, levels=np.arange(-1.8, 1.9, 0.1), cmap='RdBu_r', extend='both', transform=data_crs, zorder=1)
add_stippling_wrf(ax3, z850_pval)

stride3 = 12
q3 = ax3.quiver(lons[::stride3, ::stride3], lats[::stride3, ::stride3], uq850_d[::stride3, ::stride3], vq850_d[::stride3, ::stride3],
                transform=data_crs, color='blue', width=0.003, scale=80, zorder=5)
ax3.quiverkey(q3, 0.76, 1.03, 5, "5 g/(kg·m/s)", labelpos='E', coordinates='axes', fontproperties={'size': 11})
plt.colorbar(cf3, ax=ax3, orientation='horizontal', fraction=0.046, pad=0.08)


# ----------------- 下排: (d, e, f) 基于 ERA5 观测 的 模式异常相关分析 -----------------
# (d) TP ST
ax4 = fig.add_subplot(2, 3, 4, projection=proj)
apply_style_and_shp(ax4)
ax4.set_title("(d) ERA5: U200 Corr with TP ST (EXP-CTL)", loc='left', fontsize=11, pad=5)
cf4 = ax4.contourf(lons, lats, u200_tp_corr_wrf, levels=lvl_corr, cmap='RdBu_r', extend='both', transform=data_crs, zorder=1)
add_stippling_wrf(ax4, u200_tp_pval_wrf)
rect_d = mpatches.Rectangle((80, 29), 12, 4, transform=data_crs, fill=False, edgecolor='darkgreen', linewidth=1.5, zorder=6)
ax4.add_patch(rect_d)
plt.colorbar(cf4, ax=ax4, orientation='horizontal', fraction=0.046, pad=0.08)

# (e) EAJI
ax5 = fig.add_subplot(2, 3, 5, projection=proj)
apply_style_and_shp(ax5)
ax5.set_title("(e) ERA5: H500 Corr with EAJI (EXP-CTL)", loc='left', fontsize=11, pad=5)
cf5 = ax5.contourf(lons, lats, h500_eaji_corr_wrf, levels=lvl_corr, cmap='RdBu_r', extend='both', transform=data_crs, zorder=1)
add_stippling_wrf(ax5, h500_eaji_pval_wrf)
rect_e = mpatches.Rectangle((80, 30), 40, 10, transform=data_crs, fill=False, edgecolor='darkgreen', linewidth=1.5, zorder=6)
ax5.add_patch(rect_e)
plt.colorbar(cf5, ax=ax5, orientation='horizontal', fraction=0.046, pad=0.08)

# ⭐ 修改处: (f) H500 idx (15-30N, 110-145E)
ax6 = fig.add_subplot(2, 3, 6, projection=proj)
apply_style_and_shp(ax6)
ax6.set_title("(f) ERA5: V850 Corr with WPSHI (EXP-CTL)", loc='left', fontsize=11, pad=5)
cf6 = ax6.contourf(lons, lats, v850_h500_corr_wrf, levels=lvl_corr, cmap='RdBu_r', extend='both', transform=data_crs, zorder=1)
add_stippling_wrf(ax6, v850_h500_pval_wrf)
# 修改为最新的框线坐标：起始点 (110E, 15N)，宽度 35，高度 15
rect_f = mpatches.Rectangle((110, 15), 35, 15, transform=data_crs, fill=False, edgecolor='darkgreen', linewidth=1.5, zorder=6)
ax6.add_patch(rect_f)
plt.colorbar(cf6, ax=ax6, orientation='horizontal', fraction=0.046, pad=0.08)

# ==============================================================================
# 6. 保存出图
# ==============================================================================
save_name = "../illustration/cpl/comprehensive_circle_correlation_ERA5.pdf"
out_dir = os.path.dirname(save_name)
os.makedirs(out_dir, exist_ok=True) 

plt.savefig(save_name, dpi=600, bbox_inches='tight')
plt.close()

print(f"绘图成功！图表无缺角、双重打点、全景居中完美版已保存至：{save_name}")