import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
import warnings
import os
import string

# 屏蔽空切片求均值的警告
warnings.filterwarnings('ignore', category=RuntimeWarning, message='Mean of empty slice')

# ========================= 1. 物理常数与土壤深度 =========================
denice = 917.0      # density of ice [kg/m3]
cpliq  = 4188.0     # Specific heat of water [J/(kg K)]
cpice  = 2117.27    # Specific heat of ice [J/(kg K)]

# CoLM 10层土壤厚度 [m]
dz = np.array([0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058])

# ========================= 2. 数据读取与处理 =========================
print("正在读取和处理数据...")
path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
f1 = xr.open_dataset(path + "colmoff_2001-2023_monmean_nogravel.nc")
f2 = xr.open_dataset(path + "colmoff_2001-2023_monmean_gravel.nc")

# 读取冬季DJF数据 (51个维度)
f1_djf = xr.open_dataset(path + "colmoff_2001-2017_DJF_nogravel.nc")
f2_djf = xr.open_dataset(path + "colmoff_2001-2017_DJF_gravel.nc")

# 提取土壤液态水和冰 (春、夏)
w_ctl_mon = f1['f_wliq_soisno'].isel(soilsnow=slice(5, 15)).values
w_exp_mon = f2['f_wliq_soisno'].isel(soilsnow=slice(5, 15)).values
i_ctl_mon = f1['f_wice_soisno'].isel(soilsnow=slice(5, 15)).values
i_exp_mon = f2['f_wice_soisno'].isel(soilsnow=slice(5, 15)).values

# 提取土壤液态水和冰 (冬)
w_ctl_win_mon = f1_djf['f_wliq_soisno'].isel(soilsnow=slice(5, 15)).values
w_exp_win_mon = f2_djf['f_wliq_soisno'].isel(soilsnow=slice(5, 15)).values
i_ctl_win_mon = f1_djf['f_wice_soisno'].isel(soilsnow=slice(5, 15)).values
i_exp_win_mon = f2_djf['f_wice_soisno'].isel(soilsnow=slice(5, 15)).values

n_years = 17
n_lat, n_lon, n_layers = w_ctl_mon.shape[1], w_ctl_mon.shape[2], 10

# 初始化季节平均数组
w_ctl_win, w_exp_win = np.zeros((n_layers, n_lat, n_lon)), np.zeros((n_layers, n_lat, n_lon))
w_ctl_spr, w_exp_spr = np.zeros((n_layers, n_lat, n_lon)), np.zeros((n_layers, n_lat, n_lon))
w_ctl_sum, w_exp_sum = np.zeros((n_layers, n_lat, n_lon)), np.zeros((n_layers, n_lat, n_lon))
i_ctl_win, i_exp_win = np.zeros((n_layers, n_lat, n_lon)), np.zeros((n_layers, n_lat, n_lon))
i_ctl_spr, i_exp_spr = np.zeros((n_layers, n_lat, n_lon)), np.zeros((n_layers, n_lat, n_lon))
i_ctl_sum, i_exp_sum = np.zeros((n_layers, n_lat, n_lon)), np.zeros((n_layers, n_lat, n_lon))

# 统一计算季节(DJF, MAM, JJA)多年平均
for k in range(n_layers):
    for i in range(n_years):
        # 冬季DJF数据是连续的51个月，每年占3个月
        start_win, end_win = 3 * i, 3 * i + 3
        # 春夏季数据在一个文件中，每年占6个月 (前3月MAM，后3月JJA)
        start_spr, end_spr = 6 * i, 6 * i + 3
        start_sum, end_sum = 6 * i + 3, 6 * i + 6
        
        # 冬季处理
        w_ctl_win[k] += np.nanmean(w_ctl_win_mon[start_win:end_win, :, :, k], axis=0) / n_years
        w_exp_win[k] += np.nanmean(w_exp_win_mon[start_win:end_win, :, :, k], axis=0) / n_years
        i_ctl_win[k] += np.nanmean(i_ctl_win_mon[start_win:end_win, :, :, k], axis=0) / n_years
        i_exp_win[k] += np.nanmean(i_exp_win_mon[start_win:end_win, :, :, k], axis=0) / n_years

        # 春季处理
        w_ctl_spr[k] += np.nanmean(w_ctl_mon[start_spr:end_spr, :, :, k], axis=0) / n_years
        w_exp_spr[k] += np.nanmean(w_exp_mon[start_spr:end_spr, :, :, k], axis=0) / n_years
        i_ctl_spr[k] += np.nanmean(i_ctl_mon[start_spr:end_spr, :, :, k], axis=0) / n_years
        i_exp_spr[k] += np.nanmean(i_exp_mon[start_spr:end_spr, :, :, k], axis=0) / n_years
        
        # 夏季处理
        w_ctl_sum[k] += np.nanmean(w_ctl_mon[start_sum:end_sum, :, :, k], axis=0) / n_years
        w_exp_sum[k] += np.nanmean(w_exp_mon[start_sum:end_sum, :, :, k], axis=0) / n_years
        i_ctl_sum[k] += np.nanmean(i_ctl_mon[start_sum:end_sum, :, :, k], axis=0) / n_years
        i_exp_sum[k] += np.nanmean(i_exp_mon[start_sum:end_sum, :, :, k], axis=0) / n_years

# 读取经纬度网格和土壤参数
f_input = xr.open_dataset("/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01")
lat2d = f_input['XLAT'].isel(Time=0).values
lon2d = f_input['XLONG'].isel(Time=0).values

f3 = xr.open_dataset("interp_soil_params.nc")
csol_ctl_interp = f3['csol_ctl_interp'].values
csol_exp_interp = f3['csol_exp_interp'].values

diff_win = np.zeros((n_layers, n_lat, n_lon))
diff_spr = np.zeros((n_layers, n_lat, n_lon))
diff_sum = np.zeros((n_layers, n_lat, n_lon))

# ========================= 3. 计算热容差 =========================
print("正在计算土壤热容...")
for k in range(10):
    p = 0 if k <= 1 else (7 if k >= 8 else k - 1)
        
    c_sol_ctl = csol_ctl_interp[p, :, :]
    c_sol_exp = csol_exp_interp[p, :, :]

    # 将 NaN 填充为 0.0 用于计算
    w_c_win = np.where(np.isnan(w_ctl_win[k]), 0.0, w_ctl_win[k])
    i_c_win = np.where(np.isnan(i_ctl_win[k]), 0.0, i_ctl_win[k])
    w_e_win = np.where(np.isnan(w_exp_win[k]), 0.0, w_exp_win[k])
    i_e_win = np.where(np.isnan(i_exp_win[k]), 0.0, i_exp_win[k])

    w_c_spr = np.where(np.isnan(w_ctl_spr[k]), 0.0, w_ctl_spr[k])
    i_c_spr = np.where(np.isnan(i_ctl_spr[k]), 0.0, i_ctl_spr[k])
    w_e_spr = np.where(np.isnan(w_exp_spr[k]), 0.0, w_exp_spr[k])
    i_e_spr = np.where(np.isnan(i_exp_spr[k]), 0.0, i_exp_spr[k])

    w_c_sum = np.where(np.isnan(w_ctl_sum[k]), 0.0, w_ctl_sum[k])
    i_c_sum = np.where(np.isnan(i_ctl_sum[k]), 0.0, i_ctl_sum[k])
    w_e_sum = np.where(np.isnan(w_exp_sum[k]), 0.0, w_exp_sum[k])
    i_e_sum = np.where(np.isnan(i_exp_sum[k]), 0.0, i_exp_sum[k])

    # 热容计算
    c_total_ctl_win = c_sol_ctl + (w_c_win / dz[k]) * cpliq + (i_c_win / dz[k]) * cpice
    c_total_exp_win = c_sol_exp + (w_e_win / dz[k]) * cpliq + (i_e_win / dz[k]) * cpice
    
    c_total_ctl_spr = c_sol_ctl + (w_c_spr / dz[k]) * cpliq + (i_c_spr / dz[k]) * cpice
    c_total_exp_spr = c_sol_exp + (w_e_spr / dz[k]) * cpliq + (i_e_spr / dz[k]) * cpice
    
    c_total_ctl_sum = c_sol_ctl + (w_c_sum / dz[k]) * cpliq + (i_c_sum / dz[k]) * cpice
    c_total_exp_sum = c_sol_exp + (w_e_sum / dz[k]) * cpliq + (i_e_sum / dz[k]) * cpice

    # 差值 (EXP - CTL)
    diff_win[k, :, :] = (c_total_exp_win - c_total_ctl_win) * 1.0e-6
    diff_spr[k, :, :] = (c_total_exp_spr - c_total_ctl_spr) * 1.0e-6
    diff_sum[k, :, :] = (c_total_exp_sum - c_total_ctl_sum) * 1.0e-6

# ========================= 4. 观测掩膜处理 =========================
f_obs = xr.open_dataset("/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc")
t_obs_mon = f_obs['tm'].isel(time=0).values
lat1d = f_obs['lat'].values
lon1d = f_obs['lon'].values

lon1d_mesh, lat1d_mesh = np.meshgrid(lon1d, lat1d)
t_mask_2d = griddata(
    (lon1d_mesh.flatten(), lat1d_mesh.flatten()), 
    t_obs_mon.flatten(), 
    (lon2d, lat2d), 
    method='nearest'
)

t_mask_3d = np.broadcast_to(t_mask_2d, diff_spr.shape)

diff_win_ip = np.where(np.isnan(t_mask_3d), np.nan, diff_win)
diff_spr_ip = np.where(np.isnan(t_mask_3d), np.nan, diff_spr)
diff_sum_ip = np.where(np.isnan(t_mask_3d), np.nan, diff_sum)

diff_win_ip = np.where(np.abs(diff_win_ip) < 1e-6, np.nan, diff_win_ip)
diff_spr_ip = np.where(np.abs(diff_spr_ip) < 1e-6, np.nan, diff_spr_ip)
diff_sum_ip = np.where(np.abs(diff_sum_ip) < 1e-6, np.nan, diff_sum_ip)

# ========================= 5. 绘图  =========================
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker

print("正在绘图...")

os.environ['CARTOPY_OFFLINE'] = 'true'
os.environ['PROJ_NETWORK'] = 'OFF'

COASTLINE_50M = cfeature.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='gray', facecolor='none')
LAKES_50M = cfeature.NaturalEarthFeature('physical', 'lakes', '50m', edgecolor='gray', facecolor='none')

proj = ccrs.LambertConformal(
    central_longitude=110,
    central_latitude=40,
    standard_parallels=(30, 60)
)

fig, axes = plt.subplots(10, 3, figsize=(18, 40), subplot_kw={'projection': proj})
plt.subplots_adjust(wspace=0.05, hspace=0.15, left=0.05, right=0.95, bottom=0.05, top=0.98)

shp_dir = "../../shapefile_China/"

def apply_style_and_shp(ax, title_str, is_left, is_bottom):
    ax.set_extent([80, 130, 15, 55], crs=ccrs.PlateCarree())
    
    ax.add_feature(COASTLINE_50M, linewidth=0.6, zorder=2)
    ax.add_feature(LAKES_50M, linewidth=0.5, zorder=2)
    
    def add_shp(ax, name, lw=0.6, color='black'):
        path = f"{shp_dir}{name}"
        if not os.path.exists(path):
            return
        try:
            import geopandas as gpd
            gdf = gpd.read_file(path)
            if gdf.crs is None:
                gdf = gdf.set_crs("EPSG:4326")
            gdf.to_crs(proj).plot(ax=ax, facecolor='none', edgecolor=color, linewidth=lw, zorder=3)
        except Exception:
            pass

    add_shp(ax, "province.shp", lw=0.4, color='black')
    add_shp(ax, "china.shp", lw=0.6, color='black')
    add_shp(ax, "south_china_sea.shp", lw=0.6, color='black')
    add_shp(ax, "river.nc", lw=0.4, color='blue')

    gl = ax.gridlines(
        crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False,
        linewidth=0.6, color='gray', alpha=0.5, linestyle='--', zorder=2
    )
    gl.xlocator = mticker.FixedLocator(np.arange(80, 140, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
    
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = is_left
    gl.bottom_labels = is_bottom
    gl.xlabel_style = {'size': 10, 'rotation': 0, 'va': 'top', 'ha': 'center'}
    gl.ylabel_style = {'size': 10, 'rotation': 0, 'va': 'center', 'ha': 'right'}
    
    ax.set_title(title_str, loc='left', fontsize=12, fontweight='bold', pad=10)

levels_diff = np.arange(-0.5, 0.52, 0.02)
cmap_diff = plt.get_cmap('RdBu_r')

letters = list(string.ascii_lowercase) + [f"a{l}" for l in string.ascii_lowercase]

for j in range(10):
    for i, season_data in enumerate([diff_win_ip[j], diff_spr_ip[j], diff_sum_ip[j]]):
        ax = axes[j, i]
        
        season_str = ["DJF", "MAM", "JJA"][i]
        letter = letters[j * 3 + i]
        title = f"({letter}) {season_str} (Layer n={j+1})"
        
        is_left_col = (i == 0)
        is_bottom_row = (j == 9)
        
        apply_style_and_shp(ax, title, is_left_col, is_bottom_row)
        
        cf = ax.contourf(
            lon2d, lat2d, season_data,
            levels=levels_diff, cmap=cmap_diff, extend='both',
            transform=ccrs.PlateCarree(), zorder=1
        )
        
        ax.text(0.98, 0.98, "EXP - CTL", transform=ax.transAxes, fontsize=9, 
                verticalalignment='top', horizontalalignment='right', 
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

# ====================== 6. 统一水平色标 ======================
cbar_ax = fig.add_axes([0.2, 0.02, 0.6, 0.008]) 
cb = plt.colorbar(cf, cax=cbar_ax, orientation='horizontal')
cb.set_label('Heat Capacity Difference ($10^6$ J m$^{-3}$ K$^{-1}$)', fontsize=14, fontweight='bold')
cb.set_ticks(np.arange(-0.6, 0.7, 0.1))
cb.ax.tick_params(labelsize=12)

# ====================== 7. 保存 ======================
save_name = "./figs/heat_capacity_diff.pdf"
os.makedirs("./figs", exist_ok=True)
plt.savefig(save_name, dpi=1200, bbox_inches='tight')
plt.close()

print(f"绘图完成保存至：{save_name}")