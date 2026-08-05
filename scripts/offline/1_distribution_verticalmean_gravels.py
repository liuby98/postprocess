import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.ticker as mticker
from scipy.interpolate import RegularGridInterpolator
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader  # 【新增】参考脚本中的 shapefile 读取模块
import warnings
warnings.filterwarnings("ignore")

# ====================== 1. 基础环境与路径设置 ======================
os.environ['CARTOPY_OFFLINE'] = 'true'
os.environ['PROJ_NETWORK'] = 'OFF'

# 强制 Cartopy 读取你存放解压文件的个人离线目录
cartopy.config['data_dir'] = os.path.expanduser('~/.local/share/cartopy')

input_file = "/share/home/dq013/zhwei/colm/data/CoLMrawdata/soil/vf_gravels_s.nc" 
cn05_file  = "/share/home/dq117/CN05.1/CN05.1_Tm_2020_daily_025x025.nc"
out_dir    = "./figs_gravel_distribution"
shp_dir    = "../../shapefile_China/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# ====================== 2. 读取 CN05.1 提取中国陆地掩膜 ======================
print("正在读取 CN05.1 掩膜数据...")
cn05 = Dataset(cn05_file)
obs_cn05 = cn05.variables['tm'][0, :, :] 
lat_cn = cn05.variables['lat'][:]
lon_cn = cn05.variables['lon'][:]
cn05.close()

if hasattr(obs_cn05, 'mask'):
    obs_cn05 = np.ma.filled(obs_cn05, np.nan)
obs_cn05 = np.where(np.abs(obs_cn05) > 1000, np.nan, obs_cn05)
mask_cn05 = ~np.isnan(obs_cn05)

if lat_cn[0] > lat_cn[-1]:
    lat_cn = lat_cn[::-1]
    mask_cn05 = mask_cn05[::-1, :]
if lon_cn[0] > lon_cn[-1]:
    lon_cn = lon_cn[::-1]
    mask_cn05 = mask_cn05[:, ::-1]

mask_interp_func = RegularGridInterpolator(
    (lat_cn, lon_cn), mask_cn05.astype(float), 
    method='nearest', bounds_error=False, fill_value=0.0
)

# ====================== 3. 读取高分辨数据坐标与掩膜 ======================
print(f"正在读取目标数据坐标: {input_file}")
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

print("正在计算中国陆地高分辨掩膜...")
pts = np.stack((lat_grid, lon_grid), axis=-1)
mask_china_grid = mask_interp_func(pts) > 0.5

# ====================== 4. 逐层读取并计算模式10层垂直加权平均 ======================
# 1. 定义模式 10 层的厚度 (dz)
dz_values = np.array([0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058])

# 2. 计算 10 层的总厚度，并以此计算每一层的权重系数 (厚度/总厚度)
total_dz = np.sum(dz_values)
weight_coefs = dz_values / total_dz

# 3. 映射原始数据 8 层对应的权重系数
layer_weights = {
    1: weight_coefs[0] + weight_coefs[1],  # 对应模式层 1, 2
    2: weight_coefs[2],                    # 对应模式层 3
    3: weight_coefs[3],                    # 对应模式层 4
    4: weight_coefs[4],                    # 对应模式层 5
    5: weight_coefs[5],                    # 对应模式层 6
    6: weight_coefs[6],                    # 对应模式层 7
    7: weight_coefs[7],                    # 对应模式层 8
    8: weight_coefs[8] + weight_coefs[9]   # 对应模式层 9, 10
}

sum_gravel_weighted = np.zeros_like(lon_grid, dtype=np.float32)
sum_weights         = np.zeros_like(lon_grid, dtype=np.float32)

for i in range(1, 9):
    var_name = f'vf_gravels_s_l{i}'
    weight = layer_weights[i]
    print(f"正在读取第 {i} 层数据: {var_name} (权重系数: {weight:.4f}) ...")
    
    data = ds.variables[var_name][lat_start:lat_end:stride, lon_start:lon_end:stride]
    
    if hasattr(data, 'mask'):
        data = np.ma.filled(data, np.nan)
        
    data = np.where(data > 1000, np.nan, data)
    data = np.where(data < 0, np.nan, data)
    
    # 累加：有效值 * 权重系数
    valid_mask = ~np.isnan(data)
    sum_gravel_weighted[valid_mask] += data[valid_mask] * weight
    # 累加有效层的权重系数 
    sum_weights[valid_mask] += weight

ds.close()

print("正在计算加权平均砾石含量...")
with np.errstate(divide='ignore', invalid='ignore'):
    mean_gravel = np.where(sum_weights > 0, sum_gravel_weighted / sum_weights, np.nan)

# 【应用掩膜】非中国区变 NaN
mean_gravel = np.where(mask_china_grid, mean_gravel, np.nan)

# ====================== 5. 定义严格的离散色标 ======================
colors = [
    'lightsteelblue',  # 0.0-0.2 蓝色
    'mediumpurple',  # 0.2-0.3 紫色
    'gold',  # 0.3-0.4 黄色
    'palegreen',  # 0.4-0.5 绿色
    'lightcoral'   # > 0.5   红色
]
cmap = ListedColormap(colors)
cmap.set_bad(color='white') 
bounds = [0.0, 0.2, 0.3, 0.4, 0.5, 1.0]
norm = BoundaryNorm(bounds, cmap.N)

# ====================== 6. 绘制全国空间分布图 (兰伯特投影) ======================
print("正在绘制全国平均分布图(兰伯特投影)...")

proj = ccrs.LambertConformal(
    central_longitude=110,
    central_latitude=40,
    standard_parallels=(30, 60)
)

fig, ax = plt.subplots(1, 1, figsize=(10, 8), subplot_kw={'projection': proj})
ax.set_extent([80, 130, 10, 55], crs=ccrs.PlateCarree())

COASTLINE_50M = cfeature.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='gray', facecolor='none')
LAKES_50M = cfeature.NaturalEarthFeature('physical', 'lakes', '50m', edgecolor='gray', facecolor='none')
ax.add_feature(COASTLINE_50M, linewidth=0.8, zorder=2)
ax.add_feature(LAKES_50M, linewidth=0.5, zorder=2)

# 【修改点】参考 2_TS_spacial_bias.py 重写 Shapefile 加载逻辑，并包含南海九段线
shapefiles = [
    ("province.shp", 0.4, 'gray'),
    ("china.shp", 0.8, 'black'),
    ("south_china_sea.shp", 0.8, 'black')  # 加入九段线文件
]

for shp_file, lw, color in shapefiles:
    p_shp = os.path.join(shp_dir, shp_file)
    if os.path.exists(p_shp):
        try:
            reader = shpreader.Reader(p_shp)
            ax.add_geometries(reader.geometries(), crs=ccrs.PlateCarree(), 
                              facecolor='none', edgecolor=color, linewidth=lw, zorder=3)
        except Exception as e:
            print(f"无法加载 {shp_file}: {e}")
            pass

cf = ax.pcolormesh(
    lon_grid, lat_grid, mean_gravel,
    cmap=cmap, norm=norm,
    transform=ccrs.PlateCarree(), zorder=1,
    rasterized=True
)

gl = ax.gridlines(
        crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False,
        linewidth=0.6, color='gray', alpha=0.5, linestyle='--', zorder=2
    )
gl.xlocator = mticker.FixedLocator(np.arange(80, 140, 10))
gl.ylocator = mticker.FixedLocator(np.arange(10, 55, 10))
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = True
gl.left_labels = True
gl.xlabel_style = {'size': 10, 'rotation': 0, 'va': 'top', 'ha': 'center'}
gl.ylabel_style = {'size': 10, 'rotation': 0, 'va': 'center', 'ha': 'right'}

pos = ax.get_position()
cbar_ax = fig.add_axes([pos.x0 + 0.05, pos.y0 - 0.08, pos.width - 0.1, 0.03])
cb = plt.colorbar(cf, cax=cbar_ax, orientation='horizontal', spacing='proportional', ticks=[0.0, 0.2, 0.3, 0.4, 0.5])

cb.ax.set_xticklabels(['0.0', '0.2', '0.3', '0.4', '0.5'])
cb.ax.text(0.9, 0.5, '>0.5', transform=cb.ax.transAxes, ha='center', va='center', color='black', fontsize=10)

save_name = f"../illustration/off/FIG1.gravel_weighted_mean.pdf"
plt.savefig(save_name, dpi=300, bbox_inches='tight')
plt.close(fig)

print(f"已成功保存: {save_name}")