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
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import os
import warnings

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
    central_longitude=105, 
    central_latitude=35,   
    standard_parallels=(30, 60)
)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ========================= 2. 数据处理函数 =========================

def get_precip_bias(season, case):
    print(f"正在处理 {case} - {season} 数据...")
    
    # 1. 读取观测数据 (OBS)
    obs_file = f"{path_ref}CN05.1_Pre_1991_2023_{season}_025x025.nc"
    ds_obs = xr.open_dataset(obs_file)
    obs_var = ds_obs['pre'].isel(time=slice(30, 30 + n_years * 3)) 
    
    obs_annual = obs_var.values.reshape(n_years, 3, obs_var.shape[1], obs_var.shape[2]).mean(axis=1)
    obs_mean = obs_annual.mean(axis=0)
    
    lat_obs, lon_obs = ds_obs.lat.values, ds_obs.lon.values
    lon_mesh, lat_mesh = np.meshgrid(lon_obs, lat_obs)

    # 2. 读取模拟数据 (SIM)
    sim_file = f"{path_data}wrfout_2001-2023_{season}_daymean_{case}.nc"
    ds_sim = xr.open_dataset(sim_file)
    s_time_dim = 'Times' if 'Times' in ds_sim.dims else 'Time'
    sim_raw = ds_sim['PRAVG'].isel({s_time_dim: slice(0, n_years * 92)}) 
    
    sim_annual = sim_raw.values.reshape(n_years, 92, sim_raw.shape[1], sim_raw.shape[2]).mean(axis=1)
    sim_mean_mm_s = sim_annual.mean(axis=0)
    sim_mean = sim_mean_mm_s * 86400.0

    # 3. 插值到观测网格 
    f_in_path = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
    if not os.path.exists(f_in_path):
         print(f"Error: wrfinput file not found at {f_in_path}. Cannot perform interpolation.")
         return None, None, None
    f_in = xr.open_dataset(f_in_path)
    
    i_time_dim = 'Time' if 'Time' in f_in['XLAT'].dims else 'Times'
    lat2d = f_in['XLAT'].isel({i_time_dim: 0}).values
    lon2d = f_in['XLONG'].isel({i_time_dim: 0}).values

    sim_ip = griddata((lon2d.flatten(), lat2d.flatten()), 
                      sim_mean.flatten(), (lon_mesh, lat_mesh), method='linear')

    # 4. 计算相对偏差百分比 
    with np.errstate(divide='ignore', invalid='ignore'):
        raw_bias = (sim_ip - obs_mean) / np.where(obs_mean > 0.01, obs_mean, np.nan) * 100.0
        bias = np.clip(raw_bias, -100, 100)
    
    # 5. 掩膜处理 
    bias = np.where(np.isnan(obs_mean), np.nan, bias)
    # bias[np.ix_((lat_obs >= 25) & (lat_obs <= 29.5), (lon_obs >= 91) & (lon_obs <= 98))] = np.nan
    
    return lon_mesh, lat_mesh, bias

# ========================= 3. 地图样式与聚合函数 =========================

def apply_reference_style(ax, letter, title_str=""):
    ax.set_extent([80, 130, 15, 55], crs=ccrs.PlateCarree())
    
    def add_local_shp(ax, name, lw=0.6, color='black', zorder=3):
        path = os.path.join(shp_dir, name)
        if os.path.exists(path):
            try:
                gdf = gpd.read_file(path)
                if gdf.crs is None:
                    gdf = gdf.set_crs("EPSG:4326")
                gdf.to_crs(ax.projection).plot(ax=ax, facecolor='none', edgecolor=color, linewidth=lw, zorder=zorder)
            except Exception as e:
                print(f"Warning: Could not plot {name}: {e}")

    add_local_shp(ax, "china.shp", lw=1.2, color='black', zorder=5) 
    add_local_shp(ax, "province.shp", lw=0.5, color='gray', zorder=4) 
    add_local_shp(ax, "south_china_sea.shp", lw=1.0, color='black', zorder=5) 

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False,
                      linewidth=0.8, color='grey', alpha=0.5, linestyle='--', zorder=2)
    gl.xlocator = mticker.FixedLocator(np.arange(70, 140, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
    
    # 强制每个子图的左侧和底部都显示经纬度
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = True
    gl.bottom_labels = True
    
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    gl.xlabel_style = {'size': 11, 'rotation': 0, 'va': 'top', 'ha': 'center'}
    gl.ylabel_style = {'size': 11, 'rotation': 0, 'va': 'center', 'ha': 'right'}

    if title_str:
       ax.set_title(title_str, loc='left', fontsize=12, fontweight='normal', pad=5)

    ax.text(0.02, 0.04, f"({letter})", transform=ax.transAxes, fontsize=12, 
            va='bottom', ha='left', zorder=10, 
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

def aggregate_scatter(mask, lon, lat, grid_size=0.75):
    valid_lon = lon[mask]
    valid_lat = lat[mask]
    if len(valid_lon) == 0:
        return np.array([]), np.array([]), np.array([])
        
    lon_bins = np.arange(70, 145, grid_size)
    lat_bins = np.arange(10, 65, grid_size)
    counts, xedges, yedges = np.histogram2d(valid_lon, valid_lat, bins=[lon_bins, lat_bins])
    
    xcenters = (xedges[:-1] + xedges[1:]) / 2
    ycenters = (yedges[:-1] + yedges[1:]) / 2
    X, Y = np.meshgrid(xcenters, ycenters, indexing='ij')
    
    valid = counts > 0
    return X[valid], Y[valid], counts[valid]

# ========================= 4. 图表核心组装 =========================

seasons_code = ["MAM", "JJA"]
seasons_name = ["Spring", "Summer"]
letters = [['a', 'b', 'c'], ['d', 'e', 'f']]

cmap_prgn = plt.get_cmap('PiYG_r') 
levels = [-100, -90, -80, -70, -60, -50, 50, 60, 70, 80, 90, 100] 
norm = mcolors.BoundaryNorm(levels, cmap_prgn.N)

# 按照 2x3 版式初始化图窗 
fig, axes = plt.subplots(2, 3, figsize=(15, 10), subplot_kw={'projection': proj})
# 【优化】因为每个子图都带有 Y 轴标签，因此加大 wspace(列距) 防止重叠
plt.subplots_adjust(wspace=0.12, hspace=0.12, left=0.05, right=0.98, bottom=0.15, top=0.92)

for r, season in enumerate(seasons_code):
    s_name = seasons_name[r]
    print(f"=========================================")
    print(f"正在分析 {season} 季节数据并生成综合图...")
    
    lon, lat, ctl_bias = get_precip_bias(season, "nogravel")
    _, _, exp_bias = get_precip_bias(season, "gravel")
    
    if lon is None or ctl_bias is None or exp_bias is None:
        print(f"Error: Skipping plotting for {season} due to missing data.")
        continue

    # ================= 第一列：CTL CPL =================
    ax_ctl = axes[r, 0]
    mask_ctl = (ctl_bias >= 50) | (ctl_bias <= -50)
    if np.any(mask_ctl):
        ax_ctl.scatter(lon[mask_ctl][::10], lat[mask_ctl][::10], c=ctl_bias[mask_ctl][::10], 
                       cmap=cmap_prgn, norm=norm, s=np.abs(ctl_bias[mask_ctl][::10])*0.7, 
                       alpha=0.9, transform=ccrs.PlateCarree(), zorder=6, 
                       edgecolors='black', linewidth=0.2)
    
    title_ctl = "CTL CPL" if r == 0 else ""
    apply_reference_style(ax_ctl, letters[r][0], title_ctl)
    
    # 在每一行第一列的最左侧(Y轴标签之外)标注 Spring 和 Summer
    ax_ctl.text(-0.15, 0.5, s_name, va='center', ha='center', rotation=90, transform=ax_ctl.transAxes, 
                fontsize=11, fontweight='normal', color='black')
    
    # ================= 第二列：EXP CPL =================
    ax_exp = axes[r, 1]
    mask_exp = (exp_bias >= 50) | (exp_bias <= -50)
    if np.any(mask_exp):
        ax_exp.scatter(lon[mask_exp][::10], lat[mask_exp][::10], c=exp_bias[mask_exp][::10], 
                       cmap=cmap_prgn, norm=norm, s=np.abs(exp_bias[mask_exp][::10])*0.7, 
                       alpha=0.9, transform=ccrs.PlateCarree(), zorder=6, 
                       edgecolors='black', linewidth=0.2)

    title_exp = "EXP CPL" if r == 0 else ""
    apply_reference_style(ax_exp, letters[r][1], title_exp)

    # ================= 第三列：Comparison =================
    ax_comp = axes[r, 2]
    
    mask_valid = ~np.isnan(ctl_bias) & ~np.isnan(exp_bias)
    mask_sig = (np.abs(ctl_bias) > 50) & (np.abs(exp_bias) > 50)
    mask_base = mask_valid & mask_sig
    
    mask_pink = mask_base & (ctl_bias > 0) & (exp_bias > 0) & (exp_bias - ctl_bias < 0)
    mask_yel = mask_base & (ctl_bias < 0) & (exp_bias < 0) & (exp_bias - ctl_bias > 0)
    mask_cyan = mask_base & (ctl_bias * exp_bias < 0) & (np.abs(exp_bias) - np.abs(ctl_bias) < 0)

    base_size = 15 
    size_scale = 4 
    comp_alpha = 0.85 # 稍微增加透明度让重叠区更透气
    comp_lw = 0.4     # 边框线稍微加粗一点点增加质感
    border_color = 'none'

    agg_lon_pink, agg_lat_pink, counts_pink = aggregate_scatter(mask_pink, lon, lat, grid_size=1.0)
    if len(agg_lon_pink) > 0:
        ax_comp.scatter(agg_lon_pink, agg_lat_pink, 
                        marker='v',  # 朝下三角形 (正偏差减小)
                        s=base_size + counts_pink * size_scale,
                        color='#b03a69', 
                        alpha=comp_alpha, edgecolors=border_color, linewidth=comp_lw,
                        transform=ccrs.PlateCarree(), zorder=7, 
                        label='Both >0 & Improved' if r==0 else "")

    agg_lon_yel, agg_lat_yel, counts_yel = aggregate_scatter(mask_yel, lon, lat, grid_size=1.0)
    if len(agg_lon_yel) > 0:
        ax_comp.scatter(agg_lon_yel, agg_lat_yel, 
                        marker='^',  # 朝上三角形 (负偏差回升)
                        s=base_size + counts_yel * size_scale,
                        color='darkolivegreen', 
                        alpha=comp_alpha, edgecolors=border_color, linewidth=comp_lw,
                        transform=ccrs.PlateCarree(), zorder=7, 
                        label='Both <0 & Improved' if r==0 else "")

    agg_lon_cyan, agg_lat_cyan, counts_cyan = aggregate_scatter(mask_cyan, lon, lat, grid_size=1.0)
    if len(agg_lon_cyan) > 0:
        ax_comp.scatter(agg_lon_cyan, agg_lat_cyan, 
                        marker='D',  # 饱满的菱形 (Diamond)
                        s=base_size + counts_cyan * size_scale,
                        color='#FDAE61', # 柔和的橘黄色
                        alpha=comp_alpha, edgecolors=border_color, linewidth=comp_lw,
                        transform=ccrs.PlateCarree(), zorder=7, 
                        label='Opposite sign & Improved' if r==0 else "")

    title_comp = "Bias Response(>50%)" if r == 0 else ""
    apply_reference_style(ax_comp, letters[r][2], title_comp)

# ========================= 5. 底部图例与色标布局 (SCI规范) =========================

# 1. 为第1列与第2列创建统一 Colorbar
cbar_ax = fig.add_axes([0.05, 0.08, 0.6, 0.015]) 
sm = cm.ScalarMappable(cmap=cmap_prgn, norm=norm)
sm.set_array([]) 
# 不再使用 cb.set_label(...) 以去除色标标题
cb = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal', spacing='uniform', extend='both')
cb.set_ticks(levels)
cb.ax.tick_params(labelsize=12)

# 2. 为第3列创建独立 Legend 图例
handles, labels = axes[0, 2].get_legend_handles_labels()
if handles:
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.84, 0.06), 
               fontsize=11, frameon=False, edgecolor=border_color)

save_path = "../illustration/cpl/precip_bias_percent.pdf"
plt.savefig(save_path, dpi=600, bbox_inches='tight')
plt.close()

print(f"绘图完成，合并图像已保存为：{save_path}")