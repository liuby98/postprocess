import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
import warnings
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker

# 屏蔽空切片求均值的警告
warnings.filterwarnings('ignore', category=RuntimeWarning, message='Mean of empty slice')

# ========================= 1. 物理常数与设定 =========================
tfrz = 273.16  # freezing temperature [K]

# ========================= 2. 数据读取与处理 =========================
print("正在读取和处理数据...")
path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
f1 = xr.open_dataset(path + "colmoff_2001-2023_monmean_nogravel.nc")
f2 = xr.open_dataset(path + "colmoff_2001-2023_monmean_gravel.nc")

# 提取土壤温度 (soilsnow=5:14 对应 NCL 中的 5:14)
t_ctl_mon = f1['f_t_soisno'].isel(soilsnow=slice(5, 15)).values
t_exp_mon = f2['f_t_soisno'].isel(soilsnow=slice(5, 15)).values

n_years = 17
n_lat, n_lon, n_layers = t_ctl_mon.shape[1], t_ctl_mon.shape[2], 10

# 初始化季节平均数组
t_ctl_spr, t_exp_spr = np.zeros((n_layers, n_lat, n_lon)), np.zeros((n_layers, n_lat, n_lon))
t_ctl_sum, t_exp_sum = np.zeros((n_layers, n_lat, n_lon)), np.zeros((n_layers, n_lat, n_lon))

# 计算春季(MAM)和夏季(JJA)多年平均
for k in range(n_layers):
    for i in range(n_years):
        start_spr, end_spr = 6 * i, 6 * i + 3
        start_sum, end_sum = 6 * i + 3, 6 * i + 6
        
        t_ctl_spr[k] += np.nanmean(t_ctl_mon[start_spr:end_spr, :, :, k], axis=0) / n_years
        t_exp_spr[k] += np.nanmean(t_exp_mon[start_spr:end_spr, :, :, k], axis=0) / n_years
        t_ctl_sum[k] += np.nanmean(t_ctl_mon[start_sum:end_sum, :, :, k], axis=0) / n_years
        t_exp_sum[k] += np.nanmean(t_exp_mon[start_sum:end_sum, :, :, k], axis=0) / n_years

# 读取经纬度网格
f_input = xr.open_dataset("/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01")
lat2d = f_input['XLAT'].isel(Time=0).values
lon2d = f_input['XLONG'].isel(Time=0).values

# 读取土壤参数和已处理好的体积含水量
f3 = xr.open_dataset("interp_soil_params.nc")
theta_ctl_spr = f3['theta_ctl_spr_re'].values
theta_exp_spr = f3['theta_exp_spr_re'].values
theta_ctl_sum = f3['theta_ctl_sum_re'].values
theta_exp_sum = f3['theta_exp_sum_re'].values

ths_c_interp = f3['theta_s_ctl_interp'].values
ths_e_interp = f3['theta_s_exp_interp'].values
tkd_c_interp = f3['tkdry_ctl_interp'].values
tkd_e_interp = f3['tkdry_exp_interp'].values
tksu_c_interp = f3['tksatu_ctl_interp'].values
tksu_e_interp = f3['tksatu_exp_interp'].values
tksf_c_interp = f3['tksatf_ctl_interp'].values
tksf_e_interp = f3['tksatf_exp_interp'].values
alf_c_interp = f3['BA_alpha_ctl_interp'].values
alf_e_interp = f3['BA_alpha_exp_interp'].values
bet_c_interp = f3['BA_beta_ctl_interp'].values
bet_e_interp = f3['BA_beta_exp_interp'].values
vs_c_interp = f3['vf_sand_s_ctl_interp'].values
vs_e_interp = f3['vf_sand_s_exp_interp'].values
vom_c_interp = f3['vf_om_s_ctl_interp'].values
vom_e_interp = f3['vf_om_s_exp_interp'].values
vgr_c_interp = f3['vf_gravels_s_ctl_interp'].values
vgr_e_interp = f3['vf_gravels_s_exp_interp'].values

diff_spr = np.zeros((n_layers, n_lat, n_lon))
diff_sum = np.zeros((n_layers, n_lat, n_lon))

# ========================= 3. 计算热导率差 (Thermal Conductivity) =========================
print("正在计算土壤热导率...")

def calc_Sr(th, ths):
    """计算相对饱和度，并防除零与越界限制"""
    # 避免 ths 出现 0 或 nan 导致除零警告
    valid_mask = (ths != 0) & (~np.isnan(ths))
    Sr = np.zeros_like(th)
    Sr[valid_mask] = th[valid_mask] / ths[valid_mask]
    return np.clip(Sr, 1e-4, 1.0)

def calc_k(t, Sr, tkd, tksu, tksf, vom, vs, vgr, alf, bet):
    """根据物理状态计算 Kersten 数与热导率（修正版）"""
    ksat = np.where(t > tfrz, tksu, tksf)
    
    base = (1.0 / (1.0 + np.exp(-bet * Sr)))**3 - ((1.0 - Sr)/2.0)**3
    base = np.maximum(base, 1e-6) 
    
    # 修正：将 power 从乘积中分离
    power = 0.5 * (1.0 + vom - alf*vs - vgr)
    
    # 修正：独立计算两部分并相乘
    ke_unf = (Sr ** power) * (base ** (1.0 - vom))
    ke_frz = Sr ** (1.0 + vom)
    
    ke = np.where(t > tfrz, ke_unf, ke_frz)
    
    return (ksat - tkd) * ke + tkd

for k in range(10):
    # 建立 10层(k) 到 8层(p) 的索引映射关系
    p = 0 if k <= 1 else (7 if k >= 8 else k - 1)
    
    # 提取第 k 层状态
    tc_spr, tc_sum = t_ctl_spr[k], t_ctl_sum[k]
    te_spr, te_sum = t_exp_spr[k], t_exp_sum[k]
    th_c_spr, th_c_sum = theta_ctl_spr[k], theta_ctl_sum[k]
    th_e_spr, th_e_sum = theta_exp_spr[k], theta_exp_sum[k]
    
    # 提取第 p 层物理参数(CTL)
    ths_c, tkd_c = ths_c_interp[p], tkd_c_interp[p]
    tksu_c, tksf_c = tksu_c_interp[p], tksf_c_interp[p]
    vom_c, vs_c, vgr_c = vom_c_interp[p], vs_c_interp[p], vgr_c_interp[p]
    alf_c, bet_c = alf_c_interp[p], bet_c_interp[p]
    
    # 提取第 p 层物理参数(EXP)
    ths_e, tkd_e = ths_e_interp[p], tkd_e_interp[p]
    tksu_e, tksf_e = tksu_e_interp[p], tksf_e_interp[p]
    vom_e, vs_e, vgr_e = vom_e_interp[p], vs_e_interp[p], vgr_e_interp[p]
    alf_e, bet_e = alf_e_interp[p], bet_e_interp[p]

    # 计算相对饱和度
    Sr_c_spr, Sr_c_sum = calc_Sr(th_c_spr, ths_c), calc_Sr(th_c_sum, ths_c)
    Sr_e_spr, Sr_e_sum = calc_Sr(th_e_spr, ths_e), calc_Sr(th_e_sum, ths_e)
    
    # 计算热导率 k
    k_c_spr = calc_k(tc_spr, Sr_c_spr, tkd_c, tksu_c, tksf_c, vom_c, vs_c, vgr_c, alf_c, bet_c)
    k_c_sum = calc_k(tc_sum, Sr_c_sum, tkd_c, tksu_c, tksf_c, vom_c, vs_c, vgr_c, alf_c, bet_c)
    
    k_e_spr = calc_k(te_spr, Sr_e_spr, tkd_e, tksu_e, tksf_e, vom_e, vs_e, vgr_e, alf_e, bet_e)
    k_e_sum = calc_k(te_sum, Sr_e_sum, tkd_e, tksu_e, tksf_e, vom_e, vs_e, vgr_e, alf_e, bet_e)

    # 计算差值 (EXP - CTL)
    diff_spr[k, :, :] = k_e_spr - k_c_spr
    diff_sum[k, :, :] = k_e_sum - k_c_sum

# ========================= 4. 观测掩膜处理 =========================
print("正在进行掩膜处理...")
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
diff_spr_ip = np.where(np.isnan(t_mask_3d), np.nan, diff_spr)
diff_sum_ip = np.where(np.isnan(t_mask_3d), np.nan, diff_sum)

# ========================= 5. 绘图 =========================
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

fig, axes = plt.subplots(10, 2, figsize=(12, 40), subplot_kw={'projection': proj})
plt.subplots_adjust(wspace=0.05, hspace=0.15, left=0.05, right=0.95, bottom=0.05, top=0.98)

shp_dir = "../../shapefile_China/"

def apply_style_and_shp(ax, title_str, is_left, is_bottom):
    # 根据NCL代码调整了中国区域显示范围
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

    add_shp(ax, "china.shp", lw=0.6, color='black')
    add_shp(ax, "south_china_sea.shp", lw=0.6, color='black')
    add_shp(ax, "river.nc", lw=0.4, color='blue')

    gl = ax.gridlines(
        crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False,
        linewidth=0.8, color='grey', alpha=0.5, linestyle='--', zorder=2
    )
    gl.xlocator = mticker.FixedLocator(np.arange(60, 140, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
    
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = is_left
    gl.bottom_labels = is_bottom
    gl.xlabel_style = {'size': 10, 'rotation': 0, 'va': 'top', 'ha': 'center'}
    gl.ylabel_style = {'size': 10, 'rotation': 0, 'va': 'center', 'ha': 'right'}
    
    ax.set_title(title_str, loc='left', fontsize=12, fontweight='bold', pad=10)

# 对齐原NCL中 -0.6 到 0.6，间隔 0.05
levels_diff = np.arange(-0.6, 0.65, 0.05)
# 使用蓝-红双色调映射对应 NCL 中的 NCV_blu_red
cmap_diff = plt.get_cmap('RdBu_r') 

for j in range(10):
    for i, season_data in enumerate([diff_spr_ip[j], diff_sum_ip[j]]):
        ax = axes[j, i]
        
        season_str = "MAM" if i == 0 else "JJA"
        letter = chr(97 + j * 2 + i)
        title = f"({letter}) {season_str} (n={j+1})"
        
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
cb.set_label('Thermal Conductivity Difference (W m$^{-1}$ K$^{-1}$)', fontsize=14, fontweight='bold')
cb.set_ticks(np.arange(-0.6, 0.7, 0.1))
cb.ax.tick_params(labelsize=12)

# ====================== 7. 保存 ======================
save_name = "./figs/thermal_conductivity_diff.pdf"
os.makedirs("./figs", exist_ok=True)
plt.savefig(save_name, dpi=1200, bbox_inches='tight')
plt.close()

print(f"绘图完成保存至：{save_name}")