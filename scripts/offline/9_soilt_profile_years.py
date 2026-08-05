import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import geopandas as gpd
from matplotlib.path import Path

# ====================== 字体全局设置 ======================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ==========================================
# 1. 定义物理常量与路径
# ==========================================
base_path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"

soil_depths = np.array([0.0071, 0.0279, 0.0623, 0.1189, 0.2122, 0.3661, 0.6198, 1.0380, 1.7276, 2.8646])

var_name = "f_t_soisno" 
layer_slice = slice(5, 15)

print("读取 WRF 经纬度网格以生成青藏高原严格掩膜...")
ds_wrf = xr.open_dataset(wrfinput_file)
lat2d = ds_wrf['XLAT'].isel(Time=0).values
lon2d = ds_wrf['XLONG'].isel(Time=0).values
ds_wrf.close()

qtp_shp_path = "../../shapefile_China/TPBoundary_new_2021/TPBoundary_new(2021).shp" 
qtp_gdf = gpd.read_file(qtp_shp_path)

poly = qtp_gdf.geometry.union_all()
polygons = list(poly.geoms) if poly.geom_type == 'MultiPolygon' else [poly]

points = np.column_stack((lon2d.ravel(), lat2d.ravel()))
qtp_mask_np = np.zeros(lat2d.shape, dtype=bool).ravel()

for p in polygons:
    path = Path(np.asarray(p.exterior.coords))
    qtp_mask_np |= path.contains_points(points, radius=0)

qtp_mask_np = qtp_mask_np.reshape(lat2d.shape)
qtp_mask_da = xr.DataArray(qtp_mask_np, dims=['lat', 'lon'])

# ==========================================
# 2. 数据读取与预处理函数
# ==========================================
def get_qtp_regional_mean(filename, time_slice=None):
    file_path = os.path.join(base_path, filename)
    ds = xr.open_dataset(file_path)
    
    if time_slice:
        ds = ds.sel(time=time_slice)
        
    da = ds[var_name].isel(soilsnow=layer_slice)
    da = da.where(da > -1000.0)
    da = da - 273.15
    da_masked = da.where(qtp_mask_da)
    da_mean = da_masked.mean(dim=['lat', 'lon'])
    
    return da_mean

# ==========================================
# 3. 读取数据并拼接时间轴 (2001-2017)
# ==========================================
print("正在读取并计算区域平均，请稍候...")

ctl_djf = get_qtp_regional_mean(base_path+'colmoff_2001-2017_DJF_nogravel.nc')
exp_djf = get_qtp_regional_mean(base_path+'colmoff_2001-2017_DJF_gravel.nc')

ctl_mj = get_qtp_regional_mean(base_path+'colmoff_2001-2023_monmean_nogravel.nc', time_slice=slice('2001-01-01', '2017-12-31'))
exp_mj = get_qtp_regional_mean(base_path+'colmoff_2001-2023_monmean_gravel.nc', time_slice=slice('2001-01-01', '2017-12-31'))

ctl_12 = ctl_djf.isel(time=slice(0, None, 3)).mean(dim='time').expand_dims(month=[12])
ctl_01 = ctl_djf.isel(time=slice(1, None, 3)).mean(dim='time').expand_dims(month=[1])
ctl_02 = ctl_djf.isel(time=slice(2, None, 3)).mean(dim='time').expand_dims(month=[2])

exp_12 = exp_djf.isel(time=slice(0, None, 3)).mean(dim='time').expand_dims(month=[12])
exp_01 = exp_djf.isel(time=slice(1, None, 3)).mean(dim='time').expand_dims(month=[1])
exp_02 = exp_djf.isel(time=slice(2, None, 3)).mean(dim='time').expand_dims(month=[2])

ctl_03_08 = ctl_mj.groupby('time.month').mean(dim='time')
exp_03_08 = exp_mj.groupby('time.month').mean(dim='time')

ctl_profile = xr.concat([ctl_12, ctl_01, ctl_02, ctl_03_08], dim='month')
exp_profile = xr.concat([exp_12, exp_01, exp_02, exp_03_08], dim='month')

diff_profile = exp_profile - ctl_profile

months = ['Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug']
x_axis = np.arange(len(months))

# ==========================================
# 4. 绘制组图 (3 Panel Hovmöller)
# ==========================================
print("正在生成图像...")
fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)

# 去除子图之间的空白行间距
plt.subplots_adjust(hspace=0.0)

X, Y = np.meshgrid(x_axis, soil_depths)

temp_levels = np.arange(-12, 13, 1)       
diff_levels = np.arange(-1.0, 1.1, 0.05)     

text_bbox = dict(facecolor='white', alpha=0, edgecolor='none', pad=2)

# --- 图 a: CTL OFF ---
c1 = axes[0].contourf(X, Y, ctl_profile.T, levels=temp_levels, cmap='RdYlBu_r', extend='both')
axes[0].invert_yaxis()
axes[0].text(0.01, 0.05, '(a) CTL OFF', transform=axes[0].transAxes, ha='left', va='bottom', fontsize=11, bbox=text_bbox)

# --- 图 b: EXP OFF ---
c2 = axes[1].contourf(X, Y, exp_profile.T, levels=temp_levels, cmap='RdYlBu_r', extend='both')
axes[1].invert_yaxis()
axes[1].text(0.01, 0.05, '(b) EXP OFF', transform=axes[1].transAxes, ha='left', va='bottom', fontsize=11, bbox=text_bbox)

# ---> (a) 和 (b) 共享一个 colorbar <---
# 通过将 ax 设置为包含 axes[0] 和 axes[1] 的切片，colorbar 将跨越这两个子图的高度
cb1 = fig.colorbar(c1, ax=axes[:2], pad=0.02, shrink=0.95, aspect=35)

# --- 图 c: Difference (EXP - CTL) ---
c3 = axes[2].contourf(X, Y, diff_profile.T, levels=diff_levels, cmap='coolwarm', extend='both')
axes[2].invert_yaxis()
axes[2].text(0.01, 0.05, '(c) (EXP - CTL) OFF', transform=axes[2].transAxes, ha='left', va='bottom', fontsize=11, bbox=text_bbox)
fig.colorbar(c3, ax=axes[2], pad=0.02)

axes[1].text(1.1, 0.5, "℃", transform=axes[1].transAxes, ha='left', va='center', fontsize=12)
# ==========================================
# 格式化坐标轴 
# ==========================================

axes[1].set_ylabel('Soil Depth (m)', labelpad=10)

for ax in axes:
    ax.tick_params(axis='both', which='both', bottom=True, top=False, left=True, right=False, direction='in')
    
    # 强制在真实的物理深度打出所有的 10 根刻度短线
    ax.set_yticks(soil_depths)
    
    # 准备完整的标签列表
    labels = [f"{d:.4f}" for d in soil_depths]
    
    # 手动将第2层和第3层的文字
    labels[1] = ""  # 隐藏 0.0279
    labels[2] = ""  # 隐藏 0.0623
    
    # 应用修改后的标签
    ax.set_yticklabels(labels, fontsize=7)
    
    # 获取刚刚设置好的所有 Y 轴标签对象
    yticklabels = ax.get_yticklabels()
    
    # 1. 最顶部的 0.0071：设置 va='top'，让文字整体位于刻度线下方（向图内收）
    yticklabels[0].set_verticalalignment('top')
    
    # 2. 最底部的 2.8646：设置 va='bottom'，让文字整体位于刻度线上方（向图内收）
    yticklabels[-1].set_verticalalignment('bottom')

axes[2].set_xticks(x_axis)
axes[2].set_xticklabels(months, fontsize=9)
axes[2].set_xlabel('Month', fontsize=9)

# 将图片保存在你跑脚本的当前目录下
output_fig = '../illustration/off/FIG9.QTP_st_profile.pdf'
plt.savefig(output_fig, bbox_inches='tight', dpi=600)
print(f"绘图完成，图片已保存为: {output_fig}")