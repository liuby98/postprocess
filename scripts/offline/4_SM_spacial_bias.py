import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, RegularGridInterpolator # === 引入 RegularGridInterpolator ===
from scipy.stats import ttest_rel
import warnings
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.ticker as mticker
import os

warnings.filterwarnings("ignore")

# ====================== 字体全局设置 ======================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ====================== 1. 路径与配置 ======================
path          = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
sim_ctl_file  = path + "colmoff_2001-2023_monmean_nogravel.nc"
sim_exp_file  = path + "colmoff_2001-2023_monmean_gravel.nc"
cn05_file     = "/share/home/dq117/CN05.1/CN05.1_Tm_2020_daily_025x025.nc"
shp_dir       = "../../shapefile_China/"
out_dir       = "../illustration/off/"

os.makedirs(out_dir, exist_ok=True)
os.environ['CARTOPY_OFFLINE'] = 'true'

# ====================== 2. 土壤水分特殊配置 ======================
seasons_config = {
    'MAM': {'name': 'Spring', 'suffix': '_MAM.nc', 'idx_start': 0, 'idx_end': 3},
    'JJA': {'name': 'Summer', 'suffix': '_JJA.nc', 'idx_start': 3, 'idx_end': 6}
}

target_layers = [
    {
        'name': '0-10cm', 
        'obs_var': 'SMs',
        'obs_prefix': '/share/home/dq135/reference/soil_moisture/SMs_1991-2023_GLEAM_v4.2a',
        'wgt': np.array([0.175, 0.276, 0.455, 0.094]),
        'sim_lev_end': 4
    }
]

# ====================== 3. 读取公共网格与 Mask ======================
print("加载公共网格与 Mask 模板...")
f_wrf = Dataset(wrfinput_file)
lat2d = f_wrf.variables['XLAT'][0, :, :]
lon2d = f_wrf.variables['XLONG'][0, :, :]
f_wrf.close()

f_obs_temp = Dataset(target_layers[0]['obs_prefix'] + '_MAM.nc')
lat1d = f_obs_temp.variables['lat'][:]
lon1d = f_obs_temp.variables['lon'][:]
f_obs_temp.close()

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

# ====================== 3.5 提取高分辨率砾石数据并生成条件掩膜 (新增) ======================
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

# ====================== 4. 核心计算函数 ======================
def regrid_rcm2rgrid(var2d, lat_in, lon_in, lat_out, lon_out):
    pts = np.column_stack((lon_in.ravel(), lat_in.ravel()))
    vals = var2d.ravel()
    gx, gy = np.meshgrid(lon_out, lat_out)
    interp = griddata(pts, vals, (gx, gy), method='linear')
    nan_mask = np.isnan(interp)
    if nan_mask.any():
        interp[nan_mask] = griddata(pts, vals, (gx, gy), method='nearest')[nan_mask]
    return interp

def get_obs_yearly_climatology(file_path, var_name, offset=30, nyear=17):
    ds = Dataset(file_path)
    data = ds.variables[var_name][offset : offset + nyear*3, :, :]
    if hasattr(ds.variables[var_name], '_FillValue'):
        data = np.where(data == ds.variables[var_name]._FillValue, np.nan, data)
    data = np.where((data < -999) | (np.abs(data) > 1e30), np.nan, data)
    ds.close()
    yearly_data = np.nanmean(data.reshape(nyear, 3, data.shape[1], data.shape[2]), axis=1)
    return yearly_data

def get_clean_model_weighted_yearly(filename, time_indices, sim_lev_end, wgt, nyear=17):
    ds = Dataset(filename)
    var = ds.variables['f_h2osoi']
    data = var[time_indices, :, :, 0:sim_lev_end]
    if hasattr(var, '_FillValue'):
        data = np.where(data == var._FillValue, np.nan, data)
    data = np.where((data < -999) | (np.abs(data) > 1e30), np.nan, data)
    weighted_data = np.nansum(data * wgt, axis=-1)
    all_nan = np.all(np.isnan(data), axis=-1)
    weighted_data[all_nan] = np.nan
    ds.close()
    yearly_data = np.nanmean(weighted_data.reshape(nyear, 3, weighted_data.shape[1], weighted_data.shape[2]), axis=1)
    return yearly_data

# === 新增：空间加权与极值提取函数 ===
def spatial_weighted_avg(data, mask, wgt):
    valid = mask & ~np.isnan(data) & ~np.isnan(wgt)
    if np.sum(valid) == 0: return np.nan
    return np.average(data[valid], weights=wgt[valid])

def get_min_max(data, mask):
    valid_data = data[mask & ~np.isnan(data)]
    if len(valid_data) == 0: 
        return np.nan, np.nan
    return np.min(valid_data), np.max(valid_data)

# ====================== 5. 准备绘图 (2 行 x 3 列) ======================
print("\n================ 开始计算与绘制 2x3 土壤水分面板 ================")
proj = ccrs.LambertConformal(central_longitude=110, central_latitude=40, standard_parallels=(30, 60))
lon_grid, lat_grid = np.meshgrid(lon1d, lat1d)

# 调整 figsize 比例，使其更符合 3 列 2 行的中国地图形态
fig, axes = plt.subplots(2, 3, figsize=(13, 8), subplot_kw={'projection': proj})
# 极度压缩 wspace 和 hspace，使视觉间距趋于一致
plt.subplots_adjust(wspace=0.01, hspace=0.03, left=0.05, right=0.95, bottom=0.12)

cf_obs = cf_bias = cf_diff = None
nyear = 17

for s_idx, (season_key, cfg) in enumerate(seasons_config.items()):
    
    time_indices = [6 * i + j for i in range(nyear) for j in range(3)] if season_key == 'MAM' else [6 * i + 3 + j for i in range(nyear) for j in range(3)]
    
    for l_idx, layer in enumerate(target_layers):
        row_idx = s_idx
        layer_name = layer['name']
        print(f"  -> 正在处理: {season_key} {layer_name} ...")
        
        obs_file = layer['obs_prefix'] + cfg['suffix']
        obs_yearly = get_obs_yearly_climatology(obs_file, layer['obs_var'], nyear=nyear)
        
        ctl_yearly_lambert = get_clean_model_weighted_yearly(sim_ctl_file, time_indices, layer['sim_lev_end'], layer['wgt'], nyear=nyear)
        exp_yearly_lambert = get_clean_model_weighted_yearly(sim_exp_file, time_indices, layer['sim_lev_end'], layer['wgt'], nyear=nyear)
        
        ctl_yearly_reg = np.zeros_like(obs_yearly)
        exp_yearly_reg = np.zeros_like(obs_yearly)
        
        for y in range(nyear):
            ctl_yearly_reg[y] = regrid_rcm2rgrid(ctl_yearly_lambert[y], lat2d, lon2d, lat1d, lon1d)
            exp_yearly_reg[y] = regrid_rcm2rgrid(exp_yearly_lambert[y], lat2d, lon2d, lat1d, lon1d)
            
        obs_mean = np.nanmean(obs_yearly, axis=0)
        ctl_reg = np.nanmean(ctl_yearly_reg, axis=0)
        exp_reg = np.nanmean(exp_yearly_reg, axis=0)
        
        _, pval_bias = ttest_rel(ctl_yearly_reg, obs_yearly, axis=0, nan_policy='omit')
        _, pval_diff = ttest_rel(exp_yearly_reg, ctl_yearly_reg, axis=0, nan_policy='omit')

        obs_masked = np.where(mask_target, obs_mean, np.nan)
        bias_reg = ctl_reg - obs_mean
        diff_reg = exp_reg - ctl_reg

        # ================== 新增：终端打印特定区域的均值与最值 ==================
        # 计算在砾石>=0.3且经度<=110°E区域的面积加权平均值
        sub_gravel_obs = spatial_weighted_avg(obs_mean, final_cond_mask, wgt_2d)
        sub_gravel_bias = spatial_weighted_avg(bias_reg, final_cond_mask, wgt_2d)
        sub_gravel_diff = spatial_weighted_avg(diff_reg, final_cond_mask, wgt_2d)
        
        # 计算区域内的最大、最小值
        obs_min, obs_max = get_min_max(obs_mean, final_cond_mask)
        bias_min, bias_max = get_min_max(bias_reg, final_cond_mask)
        diff_min, diff_max = get_min_max(diff_reg, final_cond_mask)
        
        print(f"    [Gravel>=0.3 & Lon<=110°E] {season_key} {layer_name} 区域特征:")
        print(f"        OBS     -> 加权平均: {sub_gravel_obs:.4f} m³/m³ | 最小值: {obs_min:.4f} m³/m³ | 最大值: {obs_max:.4f} m³/m³")
        print(f"        CTL-OBS -> 加权平均: {sub_gravel_bias:.4f} m³/m³ | 最小值: {bias_min:.4f} m³/m³ | 最大值: {bias_max:.4f} m³/m³")
        print(f"        EXP-CTL -> 加权平均: {sub_gravel_diff:.4f} m³/m³ | 最小值: {diff_min:.4f} m³/m³ | 最大值: {diff_max:.4f} m³/m³")
        # ======================================================================

        bias_masked = np.where(mask_target & (pval_bias < 0.05), ctl_reg - obs_masked, np.nan)
        diff_masked = np.where(mask_target & (pval_diff < 0.05), diff_reg, np.nan)
        
        # ================= 渲染当前行地图 =================
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
            
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False,
                              linewidth=0.6, color='gray', alpha=0.5, linestyle='--', zorder=2)
            gl.xlocator = mticker.FixedLocator(np.arange(60, 140, 10))
            gl.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
            gl.top_labels = False
            gl.right_labels = False
            gl.bottom_labels = True if row_idx == 1 else False 
            gl.left_labels = True if col_idx == 0 else False    
            gl.xlabel_style = {'size': 10, 'rotation': 0, 'ha': 'center', 'va': 'top'}
            gl.ylabel_style = {'size': 10, 'rotation': 0, 'ha': 'right', 'va': 'center'}

        if row_idx == 0:
            ax_row[0].set_title("OBS", loc='left', fontsize=11, pad=5)
            ax_row[1].set_title("CTL OFF - OBS", loc='left', fontsize=11, pad=5)
            ax_row[2].set_title("(EXP - CTL) OFF", loc='left', fontsize=11, pad=5)
            ax_row[2].set_title(r"m$^3$/m$^3$", loc='right', fontsize=11, pad=5)
            
        ax_row[0].text(-0.25, 0.5, cfg['name'], transform=ax_row[0].transAxes, fontsize=15, fontweight='normal', va='center', ha='center', rotation=90)

        # 第三列的 diff levels 调整为 -0.14 到 +0.14
        levels_obs  = np.arange(-0.525, 0.526, 0.05)
        levels_bias = np.arange(-0.21, 0.211, 0.02)
        levels_diff = np.arange(-0.14, 0.141, 0.01)
        
        cf_obs = ax_row[0].contourf(lon_grid, lat_grid, obs_masked, levels=levels_obs, cmap='BrBG', extend='both', transform=ccrs.PlateCarree(), zorder=1)
        cf_bias = ax_row[1].contourf(lon_grid, lat_grid, bias_masked, levels=levels_bias, cmap='BrBG', extend='both', transform=ccrs.PlateCarree(), zorder=1)
        cf_diff = ax_row[2].contourf(lon_grid, lat_grid, diff_masked, levels=levels_diff, cmap='BrBG', extend='both', transform=ccrs.PlateCarree(), zorder=1)

        base_char = 97 + row_idx * 3
        for i in range(3):
            ax_row[i].text(0.03, 0.96, f"({chr(base_char + i)}) {layer_name}", 
                           transform=ax_row[i].transAxes, fontsize=11, va='top', ha='left', zorder=5, 
                           bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

# ================= 底部统一独立 Colorbars 布局 =================
fig.canvas.draw()

# 调整了 y0 的偏移量，让 Colorbar 离子图更近
cbar_offset = 0.06

pos0 = axes[1, 0].get_position()
cax_obs = fig.add_axes([pos0.x0, pos0.y0 - cbar_offset, pos0.width, 0.018])
plt.colorbar(cf_obs, cax=cax_obs, orientation='horizontal', ticks=np.arange(-0.5, 0.51, 0.2))

pos1 = axes[1, 1].get_position()
cax_bias = fig.add_axes([pos1.x0, pos1.y0 - cbar_offset, pos1.width, 0.018])
plt.colorbar(cf_bias, cax=cax_bias, orientation='horizontal', ticks=np.arange(-0.2, 0.21, 0.05))

pos2 = axes[1, 2].get_position()
cax_diff = fig.add_axes([pos2.x0, pos2.y0 - cbar_offset, pos2.width, 0.018])
# 范围调整为 -0.14 到 +0.14，每隔 0.02 显示标签 (使用 np.round 避免浮点数精度带来的长尾数字)
diff_ticks = np.round(np.arange(-0.14, 0.141, 0.04), 2)
cbar_diff = plt.colorbar(cf_diff, cax=cax_diff, orientation='horizontal', ticks=diff_ticks)
cbar_diff.ax.tick_params(labelsize=9) # 刻度稍微缩小，防止数字密集重叠

save_name = os.path.join(out_dir, "FIG4.sm_diff_raw.pdf")
plt.savefig(save_name, dpi=300, bbox_inches='tight')
print(f"  -> 已保存至: {save_name}")
plt.close()