import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import geopandas as gpd
from scipy.interpolate import griddata
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
start_year = 2001 # 设定起始年份

# 投影参数
proj = ccrs.LambertConformal(
    central_longitude=110, 
    central_latitude=40,   
    standard_parallels=(30, 60)
)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ========================= 2. 数据处理辅助函数 =========================

def regrid_rcm2rgrid(var2d, lat2d, lon2d, lon_mesh, lat_mesh):
    pts = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    vals = var2d.ravel()
    interp = griddata(pts, vals, (lon_mesh, lat_mesh), method='linear')
    nan_mask = np.isnan(interp)
    if nan_mask.any():
        interp[nan_mask] = griddata(pts, vals, (lon_mesh, lat_mesh), method='nearest')[nan_mask]
    return interp

def get_min_max(data, mask):
    valid_data = data[mask & ~np.isnan(data)]
    if len(valid_data) == 0: 
        return np.nan, np.nan
    return np.min(valid_data), np.max(valid_data)

def spatial_weighted_avg(data, mask, wgt):
    valid = mask & ~np.isnan(data) & ~np.isnan(wgt)
    if np.sum(valid) == 0: return np.nan
    return np.average(data[valid], weights=wgt[valid])

# 读取模式网格
f_in_path = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
f_in = xr.open_dataset(f_in_path)
i_time_dim = 'Time' if 'Time' in f_in['XLAT'].dims else 'Times'
lat2d = f_in['XLAT'].isel({i_time_dim: 0}).values
lon2d = f_in['XLONG'].isel({i_time_dim: 0}).values
f_in.close()

# ========================= 2.5 提取高分辨率砾石数据及辅助函数 =========================
print("正在读取观测网格底图并生成砾石条件掩膜...")
# 临时读取一个观测文件来获取标准 0.25 经纬度网格
obs_temp_file = f"{path_ref}CN05.1_Pre_1991_2023_MAM_025x025.nc"
ds_obs_temp = xr.open_dataset(obs_temp_file)
lon_obs, lat_obs = ds_obs_temp.lon.values, ds_obs_temp.lat.values
ds_obs_temp.close()

lon_mesh, lat_mesh = np.meshgrid(lon_obs, lat_obs)
wgt_2d = np.cos(np.radians(lat_mesh))

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
pts = np.stack((lat_mesh, lon_mesh), axis=-1)
mean_gravel_025 = interp_func(pts)

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


# ========================= 4. 数据预读取与处理 (精细到月) =========================
print("=========================================")
print("开始预读取并严格按照 3-8 月份进行切片与插值...")
data_all = {}
months_list = [3, 4, 5, 6, 7, 8]

for season in ["MAM", "JJA"]:
    # 1. 读取观测数据并转换为 (年份, 3个月, lat, lon)
    obs_file = f"{path_ref}CN05.1_Pre_1991_2023_{season}_025x025.nc"
    ds_obs = xr.open_dataset(obs_file)
    var_obs = ds_obs['pre'].isel(time=slice(30, 30 + n_years * 3)).values
    var_obs_reshaped = var_obs.reshape(n_years, 3, var_obs.shape[1], var_obs.shape[2])
    ds_obs.close()

    # 2. 读取模拟数据并转换为 (年份, 92天, lat_rcm, lon_rcm)
    ctl_file = f"{path_data}wrfout_2001-2023_{season}_daymean_nogravel.nc"
    ds_ctl = xr.open_dataset(ctl_file)
    s_time_dim = 'Times' if 'Times' in ds_ctl.dims else 'Time'
    var_ctl_reshaped = ds_ctl['PRAVG'].isel({s_time_dim: slice(0, n_years * 92)}).values.reshape(n_years, 92, lat2d.shape[0], lon2d.shape[1]) * 86400.0
    ds_ctl.close()

    exp_file = f"{path_data}wrfout_2001-2023_{season}_daymean_gravel.nc"
    ds_exp = xr.open_dataset(exp_file)
    var_exp_reshaped = ds_exp['PRAVG'].isel({s_time_dim: slice(0, n_years * 92)}).values.reshape(n_years, 92, lat2d.shape[0], lon2d.shape[1]) * 86400.0
    ds_exp.close()

    # 定义不同季节下包含的具体月份、观测索引以及 WRF 的自然天数切片边界
    if season == "MAM":
        month_configs = [
            (3, "March", 0, 0, 31),    # 3月有 31 天
            (4, "April", 1, 31, 61),   # 4月有 30 天
            (5, "May",   2, 61, 92)    # 5月有 31 天
        ]
    else:
        month_configs = [
            (6, "June",   0, 0, 30),   # 6月有 30 天
            (7, "July",   1, 30, 61),  # 7月有 31 天
            (8, "August", 2, 61, 92)   # 8月有 31 天
        ]

    for m_num, m_name, obs_idx, sim_s, sim_e in month_configs:
        print(f"正在抽取并插值 {m_name} ({m_num}月) 数据...")
        # 截取该月的观测
        obs_m = var_obs_reshaped[:, obs_idx, :, :]
        # 计算该月的日均值 (平均 30 或 31 天)
        ctl_m_raw = var_ctl_reshaped[:, sim_s:sim_e, :, :].mean(axis=1)
        exp_m_raw = var_exp_reshaped[:, sim_s:sim_e, :, :].mean(axis=1)

        ctl_m_regrid = np.zeros_like(obs_m)
        exp_m_regrid = np.zeros_like(obs_m)

        for y in range(n_years):
            ctl_m_regrid[y] = regrid_rcm2rgrid(ctl_m_raw[y], lat2d, lon2d, lon_mesh, lat_mesh)
            exp_m_regrid[y] = regrid_rcm2rgrid(exp_m_raw[y], lat2d, lon2d, lon_mesh, lat_mesh)

        data_all[m_num] = {
            'name': m_name,
            'obs': obs_m,
            'ctl': ctl_m_regrid,
            'exp': exp_m_regrid
        }

# ========================= 5. 离散色带环境配置 =========================
levels_obs = np.arange(0, 11, 1)
cmap_obs = plt.get_cmap('GnBu').copy()
norm_obs = mcolors.BoundaryNorm(levels_obs, ncolors=cmap_obs.N, extend='max')

levels_bias = np.arange(-100, 101, 20)  
cmap_bias = plt.get_cmap('RdBu').copy()
norm_bias = mcolors.BoundaryNorm(levels_bias, ncolors=cmap_bias.N, extend='both')

levels_diff = np.arange(-20, 21, 4)    
cmap_diff = plt.get_cmap('RdBu').copy()
norm_diff = mcolors.BoundaryNorm(levels_diff, ncolors=cmap_diff.N, extend='both')

# ========================= 6. 循环逐年逐月绘图输出 =========================

for y_idx in range(n_years):
    current_year = start_year + y_idx
    
    for m_num in months_list:
        m_name = data_all[m_num]['name']
        print(f"\n=========================================")
        print(f"正在处理并绘制 {current_year} 年 {m_name} ({m_num}月) 的数据...")
        
        # 按照每月一张图(1行4列)创建画板
        fig, axes = plt.subplots(1, 4, figsize=(16.5, 4.5), subplot_kw={'projection': proj})
        plt.subplots_adjust(wspace=0.01, hspace=0.01, left=0.05, right=0.98, bottom=0.25, top=0.88)

        # 提取当前年、月的数据
        obs_month = data_all[m_num]['obs'][y_idx]
        ctl_month = data_all[m_num]['ctl'][y_idx]
        exp_month = data_all[m_num]['exp'][y_idx]
        
        with np.errstate(divide='ignore', invalid='ignore'):
            ctl_bias_raw = (ctl_month - obs_month) / np.where(obs_month > 0.01, obs_month, np.nan) * 100.0
            exp_bias_raw = (exp_month - obs_month) / np.where(obs_month > 0.01, obs_month, np.nan) * 100.0
        
        diff_month = exp_month - ctl_month            
        diff_bias_raw = exp_bias_raw - ctl_bias_raw 
        
        valid_mask = ~np.isnan(obs_month)
        
        # 用有效掩膜过滤 raw 统计数据
        ctl_bias_raw = np.where(valid_mask, ctl_bias_raw, np.nan)
        exp_bias_raw = np.where(valid_mask, exp_bias_raw, np.nan)
        diff_bias_raw = np.where(valid_mask, diff_bias_raw, np.nan)

        # 转换为安全的掩膜数组，以便后续 pcolormesh 渲染
        obs_plot       = np.ma.masked_invalid(obs_month)
        ctl_bias_plot  = np.ma.masked_invalid(np.clip(ctl_bias_raw, -105.0, 105.0))
        exp_bias_plot  = np.ma.masked_invalid(np.clip(exp_bias_raw, -105.0, 105.0))
        diff_bias_plot = np.ma.masked_invalid(np.clip(diff_bias_raw, -22.0, 22.0))

        # ================= 统计打印逻辑 =================
        exclude_mask = (lat_mesh >= 25.0) & (lat_mesh <= 29.5) & (lon_mesh >= 91.0) & (lon_mesh <= 98.0)
        mask_all_calc = valid_mask & (~exclude_mask)
        mask_gravel_calc = (mean_gravel_025 >= 0.3) & (lon_mesh < 110.0) & valid_mask & (~exclude_mask)
        mask_sc_calc = (lon_mesh >= 105.0) & (lon_mesh <= 115.0) & (lat_mesh >= 20.0) & (lat_mesh <= 25.0) & valid_mask & (~exclude_mask)
        
        obs_min_all, obs_max_all = get_min_max(obs_month, mask_all_calc)
        obs_min_g, obs_max_g = get_min_max(obs_month, mask_gravel_calc)
        obs_mean_g = spatial_weighted_avg(obs_month, mask_gravel_calc, wgt_2d)
        
        diff_min_sc, diff_max_sc = get_min_max(diff_month, mask_sc_calc)
        diff_mean_sc = spatial_weighted_avg(diff_month, mask_sc_calc, wgt_2d)
        
        diff_b_min_sc, diff_b_max_sc = get_min_max(diff_bias_raw, mask_sc_calc)
        diff_b_mean_sc = spatial_weighted_avg(diff_bias_raw, mask_sc_calc, wgt_2d)
        
        print(f"  [{current_year} {m_name}] 降水特征统计:")
        print(f"    - OBS (全国) -> [最小值: {obs_min_all:.2f} | 最大值: {obs_max_all:.2f}] mm/day")
        print(f"    - OBS (砾石>=0.3 & 经度<110°E) -> [最值: {obs_min_g:.2f} ~ {obs_max_g:.2f} | 均值: {obs_mean_g:.2f}] mm/day")
        print(f"    - (EXP-CTL) (华南) -> 差值[均值: {diff_mean_sc:.2f}] mm/day | 偏差%差值[均值: {diff_b_mean_sc:.2f}%]")

        # ================= 制图布局 =================
        ax_obs, ax_ctl, ax_exp, ax_comp = axes

        # 1. OBS
        cf0 = ax_obs.pcolormesh(lon_mesh, lat_mesh, obs_plot, norm=norm_obs, cmap=cmap_obs, transform=ccrs.PlateCarree(), zorder=1, shading='auto')
        apply_reference_style(ax_obs, 'a', f"OBS ({current_year}-{m_num:02d})", is_left=True, is_bottom=True)
        ax_obs.text(-0.16, 0.5, m_name, va='center', ha='center', rotation=90, transform=ax_obs.transAxes, fontsize=12, fontweight='normal', color='black')

        # 2. CTL
        cf1 = ax_ctl.pcolormesh(lon_mesh, lat_mesh, ctl_bias_plot, norm=norm_bias, cmap=cmap_bias, transform=ccrs.PlateCarree(), zorder=1, shading='auto')
        apply_reference_style(ax_ctl, 'b', "CTL CPL", is_left=False, is_bottom=True)
        
        # 3. EXP
        cf2 = ax_exp.pcolormesh(lon_mesh, lat_mesh, exp_bias_plot, norm=norm_bias, cmap=cmap_bias, transform=ccrs.PlateCarree(), zorder=1, shading='auto')
        apply_reference_style(ax_exp, 'c', "EXP CPL", is_left=False, is_bottom=True)

        # 4. DIFF
        cf3 = ax_comp.pcolormesh(lon_mesh, lat_mesh, diff_bias_plot, norm=norm_diff, cmap=cmap_diff, transform=ccrs.PlateCarree(), zorder=1, shading='auto')
        apply_reference_style(ax_comp, 'd', "(EXP - CTL) CPL", is_left=False, is_bottom=True)

        # ================= 底部色标动态布局 =================
        fig.canvas.draw()

        pos_ax0 = ax_obs.get_position()
        pos_ax1 = ax_ctl.get_position()
        pos_ax2 = ax_exp.get_position()
        pos_ax3 = ax_comp.get_position()

        cb_height = 0.035
        cb_y = pos_ax0.y0 - 0.16

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

        # ================= 保存批量图表 (.png) =================
        save_path = f"./figs_mon/precip_bias_percent_raw_{current_year}_{m_num:02d}.png"
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"  --> 绘图完成，已保存为：{save_path}")

print("\n所有年份逐月批量拆分输出完毕！共生成 102 张图。")