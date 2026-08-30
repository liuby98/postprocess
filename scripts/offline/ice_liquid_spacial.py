import os
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.spatial import cKDTree  # 引入更快的树形检索提速
import geopandas as gpd            # 引入预加载提速
import warnings
warnings.filterwarnings("ignore")

# ====================== 1. 环境变量 & 投影工具 ======================
os.environ['CARTOPY_OFFLINE'] = 'true'
os.environ['PROJ_NETWORK'] = 'OFF'

# Lambert Conformal 投影参数
proj = ccrs.LambertConformal(
    central_longitude=110,
    central_latitude=40,
    standard_parallels=(30, 60)
)

# ====================== 2. 控制开关 & 路径配置 ======================
season_switch = "DJF"  # 选项: "DJF", "MAM", "JJA"
var_switch    = "ice"  # 选项: "ice", "liquid"

path          = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
cn05_file     = "/share/home/dq117/CN05.1/CN05.1_Tm_2020_daily_025x025.nc"
shp_dir       = "../../shapefile_China/"

# 根据季节开关读取文件
if season_switch == "DJF":
    file_ctl = path + "colmoff_2001-2017_DJF_nogravel.nc"
    file_exp = path + "colmoff_2001-2017_DJF_gravel.nc"
    nyears = 17
else:
    file_ctl = path + "colmoff_2001-2023_monmean_nogravel.nc"
    file_exp = path + "colmoff_2001-2023_monmean_gravel.nc"
    nyears = 17 

# 变量参数设定
if var_switch == "ice":
    vname = "f_wice_soisno"
    lb_title_raw = "Ice Lens Content (kg/m²)"
    lb_title_dif = "Ice Lens Difference (kg/m²)"
    levels_raw = np.arange(0.0, 50.1, 0.1)
    levels_dif = np.arange(-80.0, 81.0, 1.0)
    cmap_raw_name = "Blues"
else:
    vname = "f_wliq_soisno"
    lb_title_raw = "Liquid Water Content (kg/m²)"
    lb_title_dif = "Liquid Water Difference (kg/m²)"
    levels_raw = np.arange(0.0, 155.0, 5.0)
    levels_dif = np.arange(-80.0, 81.0, 1.0)
    cmap_raw_name = "YlGnBu"

out_name = f"./figs/colmoff_{season_switch}_{var_switch}_spacial.pdf"

# ====================== 3. 读取模拟数据与坐标 ======================
print(f"正在读取 {season_switch} 数据并计算均值...")

f_wrf = Dataset(wrfinput_file)
lat2d = f_wrf.variables['XLAT'][0, :, :]
lon2d = f_wrf.variables['XLONG'][0, :, :]
f_wrf.close()

ds_ctl = Dataset(file_ctl)
ds_exp = Dataset(file_exp)
var_ctl = ds_ctl.variables[vname]
var_exp = ds_exp.variables[vname]

# 初始化 10 层的大小
ctl_season = np.zeros((nyears, lat2d.shape[0], lat2d.shape[1], 10))
exp_season = np.zeros((nyears, lat2d.shape[0], lat2d.shape[1], 10))

# 提取完整的 10 层 (索引 5:15)
for yr in range(nyears):
    if season_switch == "DJF":
        ctl_season[yr] = np.nanmean(var_ctl[3*yr : 3*yr+3, :, :, 5:15], axis=0)
        exp_season[yr] = np.nanmean(var_exp[3*yr : 3*yr+3, :, :, 5:15], axis=0)
    elif season_switch == "MAM":
        ctl_season[yr] = np.nanmean(var_ctl[6*yr : 6*yr+3, :, :, 5:15], axis=0)
        exp_season[yr] = np.nanmean(var_exp[6*yr : 6*yr+3, :, :, 5:15], axis=0)
    elif season_switch == "JJA":
        ctl_season[yr] = np.nanmean(var_ctl[6*yr+3 : 6*yr+6, :, :, 5:15], axis=0)
        exp_season[yr] = np.nanmean(var_exp[6*yr+3 : 6*yr+6, :, :, 5:15], axis=0)

ds_ctl.close()
ds_exp.close()

# 多年平均
raw_ctl = np.nanmean(ctl_season, axis=0).transpose(2, 0, 1)
raw_exp = np.nanmean(exp_season, axis=0).transpose(2, 0, 1)
raw_dif = raw_exp - raw_ctl

# ====================== 4. CN05.1 掩膜处理 (KDTree 加速) ======================
print("提取并映射观测掩膜...")
cn05 = Dataset(cn05_file)
obs_cn05 = np.nanmean(cn05.variables['tm'][:], axis=0)
lat_cn = cn05.variables['lat'][:]
lon_cn = cn05.variables['lon'][:]
cn05.close()

mask_cn05 = ~np.isnan(obs_cn05)
grid_x, grid_y = np.meshgrid(lon_cn, lat_cn)

points = np.column_stack((grid_x.ravel(), grid_y.ravel()))
values = mask_cn05.ravel().astype(float)

# 使用 cKDTree 加速掩膜映射
tree = cKDTree(points)
_, indices = tree.query(np.column_stack((lon2d.ravel(), lat2d.ravel())))
mask_interp = values[indices].reshape(lon2d.shape)
mask_target = mask_interp > 0.5

# 循环 10 层施加掩膜
for z in range(10):
    raw_ctl[z] = np.where(mask_target & (raw_ctl[z] > -1000.), raw_ctl[z], np.nan)
    raw_exp[z] = np.where(mask_target & (raw_exp[z] > -1000.), raw_exp[z], np.nan)
    raw_dif[z] = np.where(mask_target & (np.abs(raw_dif[z]) >= 0.001), raw_dif[z], np.nan)

# ====================== 5. 绘图样式配置 (预加载提速) ======================
print("预加载底图数据...")

def pre_load_shp(name):
    path = f"{shp_dir}{name}"
    if os.path.exists(path):
        try:
            gdf = gpd.read_file(path)
            if gdf.crs is None: 
                gdf = gdf.set_crs("EPSG:4326")
            return gdf.to_crs(proj)
        except: pass
    return None

china_gdf = pre_load_shp("china.shp")
scs_gdf   = pre_load_shp("south_china_sea.shp")
river_gdf = pre_load_shp("river.nc")

COASTLINE_50M = cfeature.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='gray', facecolor='none')
LAKES_50M = cfeature.NaturalEarthFeature('physical', 'lakes', '50m', edgecolor='gray', facecolor='none')

def apply_style_and_shp(ax, title_str, is_left, is_bottom):
    ax.set_extent([80, 130, 15, 55], crs=ccrs.PlateCarree())
    
    ax.add_feature(COASTLINE_50M, linewidth=0.6, zorder=2)
    ax.add_feature(LAKES_50M, linewidth=0.5, zorder=2)
    
    # 内存直读绘制，避免反复读取硬盘
    if china_gdf is not None: china_gdf.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=0.6, zorder=3)
    if scs_gdf is not None:   scs_gdf.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=0.6, zorder=3)
    if river_gdf is not None: river_gdf.plot(ax=ax, facecolor='none', edgecolor='blue', linewidth=0.4, zorder=3)

    gl = ax.gridlines(
        crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False,
        linewidth=0.8, color='grey', alpha=0.5, linestyle='--', zorder=2
    )
    gl.xlocator = mticker.FixedLocator(np.arange(60, 140, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
    
    gl.top_labels = False
    gl.right_labels = False
    # 【恢复动态轴标签，避免 10 层互相遮挡】
    gl.left_labels = True
    gl.bottom_labels = True
    
    gl.xlabel_style = {'size': 9, 'rotation': 0, 'va': 'top', 'ha': 'center'}
    gl.ylabel_style = {'size': 9, 'rotation': 0, 'va': 'center', 'ha': 'right'}
    
    ax.set_title(title_str, loc='left', fontsize=10, fontweight='bold', pad=4)

purples = plt.cm.Purples_r(np.linspace(0.1, 0.9, 42)) 
blues = plt.cm.Blues(np.linspace(0.1, 0.9, 42))
cmap_dif = mcolors.LinearSegmentedColormap.from_list('purple_blue', np.vstack((purples, blues)))
# cmap_dif = plt.get_cmap("BrBG")
cmap_raw = plt.get_cmap(cmap_raw_name)

# ====================== 6. 执行画图 ======================
print("正在绘制 10x3 面板空间图...")
# 画布拉长为 32 适应 10 行
fig, axes = plt.subplots(nrows=10, ncols=3, figsize=(15, 32), subplot_kw={'projection': proj})

# 保持紧凑的排列间距，为底部色标留出 8% 的空间
plt.subplots_adjust(wspace=0.05, hspace=0.18, left=0.12, right=0.90, bottom=0.08)

for z in range(10):
    ax_ctl = axes[z, 0]
    ax_exp = axes[z, 1]
    ax_dif = axes[z, 2]
    
    is_bottom_row = (z == 9)
    
    # CTL
    apply_style_and_shp(ax_ctl, f"CTL {season_switch}", is_left=True, is_bottom=is_bottom_row)
    cf_raw = ax_ctl.contourf(lon2d, lat2d, raw_ctl[z], levels=levels_raw, cmap=cmap_raw, extend='max', transform=ccrs.PlateCarree(), zorder=1)
    
    # 左侧层数标注恢复为 n=z+1
    ax_ctl.text(-0.18, 0.5, f"n={z+1}", va='center', ha='center', 
                rotation=90, transform=ax_ctl.transAxes, fontsize=10, fontweight='bold')

    # EXP
    apply_style_and_shp(ax_exp, f"EXP {season_switch}", is_left=False, is_bottom=is_bottom_row)
    ax_exp.contourf(lon2d, lat2d, raw_exp[z], levels=levels_raw, cmap=cmap_raw, extend='max', transform=ccrs.PlateCarree(), zorder=1)

    # Diff
    apply_style_and_shp(ax_dif, f"EXP-CTL {season_switch}", is_left=False, is_bottom=is_bottom_row)
    cf_dif = ax_dif.contourf(lon2d, lat2d, raw_dif[z], levels=levels_dif, cmap=cmap_dif, extend='both', transform=ccrs.PlateCarree(), zorder=1)

# ====================== 7. 统一色标 ======================
# 水平色标在 0.04 的高度，与图像保持适中距离
cbar_ax_raw = fig.add_axes([0.15, 0.04, 0.48, 0.01]) 
cb_raw = fig.colorbar(cf_raw, cax=cbar_ax_raw, orientation='horizontal')
cb_raw.set_label(lb_title_raw, fontsize=14, fontweight='bold')
cb_raw.ax.tick_params(labelsize=11)

# 垂直色标
cbar_ax_dif = fig.add_axes([0.92, 0.05, 0.015, 0.80]) 
cb_dif = fig.colorbar(cf_dif, cax=cbar_ax_dif, orientation='vertical')
cb_dif.set_label(lb_title_dif, fontsize=14, fontweight='bold')
cb_dif.ax.tick_params(labelsize=11)

# ====================== 8. 保存 ======================
os.makedirs("./figs", exist_ok=True)
plt.savefig(out_name, dpi=600, bbox_inches='tight')
plt.close()

print(f"绘图完成保存至：{out_name}")