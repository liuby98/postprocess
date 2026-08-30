#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from scipy.interpolate import RegularGridInterpolator
import warnings
warnings.filterwarnings("ignore")

# ====================== 1. 基础环境设置 ======================
os.environ['CARTOPY_OFFLINE'] = 'true'
os.environ['PROJ_NETWORK'] = 'OFF'

input_file = "/share/home/dq013/zhwei/colm/data/CoLMrawdata/soil/vf_gravels_s.nc" 
cn05_file  = "/share/home/dq117/CN05.1/CN05.1_Tm_2020_daily_025x025.nc"
shp_dir    = "../../shapefile_China/"
out_dir    = "./figs_gravel_distribution"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# ====================== 2. 地图要素与投影 ======================
proj = ccrs.LambertConformal(
    central_longitude=110,
    central_latitude=40,
    standard_parallels=(30, 60)
)

COASTLINE_50M = cfeature.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='gray', facecolor='none')
LAKES_50M = cfeature.NaturalEarthFeature('physical', 'lakes', '50m', edgecolor='gray', facecolor='none')

def apply_style_and_shp(ax, title_str):
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

    add_shp(ax, "province.shp", lw=0.4, color='gray')
    add_shp(ax, "china.shp", lw=0.8, color='black')
    add_shp(ax, "south_china_sea.shp", lw=0.8, color='black')

    gl = ax.gridlines(
        crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False,
        linewidth=0.6, color='gray', alpha=0.5, linestyle='--', zorder=2
    )
    gl.xlocator = mticker.FixedLocator(np.arange(70, 150, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.xlabel_style = {'size': 10, 'rotation': 0, 'va': 'top', 'ha': 'center'}
    gl.ylabel_style = {'size': 10, 'rotation': 0, 'va': 'center', 'ha': 'right'}
    
    ax.set_title(title_str, loc='left', fontsize=12, fontweight='bold', pad=10)

# ====================== 3. 读取 CN05.1 提取精确掩膜 ======================
print("正在读取 CN05.1 掩膜数据...")
cn05 = Dataset(cn05_file)
obs_cn05 = cn05.variables['tm'][0, :, :] 
lat_cn = cn05.variables['lat'][:]
lon_cn = cn05.variables['lon'][:]
cn05.close()

# 【关键修复】：安全处理 MaskedArray 与极端缺省值
if hasattr(obs_cn05, 'mask'):
    obs_cn05 = np.ma.filled(obs_cn05, np.nan)
# 兜底：将异常大/小值也替换为 NaN
obs_cn05 = np.where(np.abs(obs_cn05) > 1000, np.nan, obs_cn05)
# 此时提取的 mask 才是完美只覆盖中国陆地的
mask_cn05 = ~np.isnan(obs_cn05)

# RegularGridInterpolator 要求坐标严格递增
if lat_cn[0] > lat_cn[-1]:
    lat_cn = lat_cn[::-1]
    mask_cn05 = mask_cn05[::-1, :]
if lon_cn[0] > lon_cn[-1]:
    lon_cn = lon_cn[::-1]
    mask_cn05 = mask_cn05[:, ::-1]

# 构建高速插值函数 (必须使用 nearest 防止布尔值在边界被模糊成小数)
mask_interp_func = RegularGridInterpolator(
    (lat_cn, lon_cn), mask_cn05.astype(float), 
    method='nearest', bounds_error=False, fill_value=0.0
)

# ====================== 4. 读取高分辨数据 & 匹配掩膜 ======================
print(f"正在读取目标数据: {input_file}")
ds = Dataset(input_file)

if 'longitude' in ds.variables:
    lon_all = ds.variables['longitude'][:]
    lat_all = ds.variables['latitude'][:]
else:
    lon_all = np.linspace(-180, 180, 86400)
    lat_all = np.linspace(90, -90, 43200)

lon_idx = np.where((lon_all >= 70) & (lon_all <= 135))[0]
lat_idx = np.where((lat_all >= 15) & (lat_all <= 55))[0]

lat_start, lat_end = np.min(lat_idx), np.max(lat_idx) + 1
lon_start, lon_end = np.min(lon_idx), np.max(lon_idx) + 1

stride = 1 
lon_subset = lon_all[lon_start:lon_end:stride]
lat_subset = lat_all[lat_start:lat_end:stride]
lon_grid, lat_grid = np.meshgrid(lon_subset, lat_subset)

print("正在将掩膜映射到高分辨网格...")
pts = np.stack((lat_grid, lon_grid), axis=-1)
mask_target_grid = mask_interp_func(pts) > 0.5

# ====================== 5. 循环填色并导出 ======================
# 高对比度色标，并将缺省/掩膜区域设为纯白
import matplotlib as mpl
cmap = mpl.colormaps.get_cmap('viridis').copy()
cmap.set_bad(color='white')

# 强制范围 0.0 - 0.5
vmin, vmax = 0.0, 0.55

for i in range(8):
    layer_num = i + 1
    var_name = f'vf_gravels_s_l{layer_num}'
    print(f"正在绘制 {var_name} ...")
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 7), subplot_kw={'projection': proj})
    plt.subplots_adjust(bottom=0.15)
    
    apply_style_and_shp(ax, f"Layer {layer_num}: Gravel Volumetric Fraction")
    
    # 读取这层数据
    data = ds.variables[var_name][lat_start:lat_end:stride, lon_start:lon_end:stride]
    
    if hasattr(data, 'mask'):
        data = np.ma.filled(data, np.nan)
        
    data = np.where(data > 1000, np.nan, data)
    data = np.where(data < 0, np.nan, data)
    
    # 【最后一步：应用掩膜】，非中国区全部变 NaN
    data = np.where(mask_target_grid, data, np.nan)
    
    cf = ax.pcolormesh(
        lon_grid, lat_grid, data,
        cmap=cmap, vmin=vmin, vmax=vmax,
        transform=ccrs.PlateCarree(), zorder=1,
        rasterized=True 
    )
    
    pos = ax.get_position()
    cbar_ax = fig.add_axes([pos.x0 + 0.05, 0.08, pos.width - 0.1, 0.03])
    cb = plt.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='max')
    cb.set_label('Gravel Volumetric Fraction', fontsize=11)
    
    save_name = f"{out_dir}/gravel_layer_{layer_num}_masked.pdf"
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.close(fig) 
    
    print(f"已保存: {save_name}")

ds.close()
print("所有图层处理完毕，完美掩膜生效！")