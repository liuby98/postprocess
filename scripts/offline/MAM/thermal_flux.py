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
import matplotlib.colors as mcolors
import matplotlib.cm as cm

# 屏蔽空切片求均值的警告
warnings.filterwarnings('ignore', category=RuntimeWarning, message='Mean of empty slice')

# ========================= 1. 物理常数与土壤深度设定 =========================
tfrz = 273.16  # freezing temperature [K]

# CoLM 10层土壤厚度 (delta z_i) [m]
dz = np.array([0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058])

# 节点深度 (z_i) [m]
z = np.array([0.0071, 0.0279, 0.0623, 0.1189, 0.2122, 0.3661, 0.6198, 1.0380, 1.7276, 2.8646])

# 界面深度 (z_{h,i}) [m]
zh = np.array([0.0175, 0.0451, 0.0906, 0.1655, 0.2891, 0.4929, 0.8289, 1.3828, 2.2961, 3.8019])

# ========================= 2. 数据读取与处理 =========================
print("正在读取和处理数据...")
path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
f1 = xr.open_dataset(path + "colmoff_2001-2023_monmean_nogravel.nc")
f2 = xr.open_dataset(path + "colmoff_2001-2023_monmean_gravel.nc")

# 提取土壤温度
t_ctl_mon = f1['f_t_soisno'].isel(soilsnow=slice(5, 15)).values
t_exp_mon = f2['f_t_soisno'].isel(soilsnow=slice(5, 15)).values

n_years = 17
n_lat, n_lon, n_layers = t_ctl_mon.shape[1], t_ctl_mon.shape[2], 10

# 初始化季节平均数组
t_ctl_spr, t_exp_spr = np.zeros((n_layers, n_lat, n_lon)), np.zeros((n_layers, n_lat, n_lon))
t_ctl_sum, t_exp_sum = np.zeros((n_layers, n_lat, n_lon)), np.zeros((n_layers, n_lat, n_lon))

# 计算春季(MAM)和夏季(JJA)多年平均
for k in range(n_layers):
    for yr in range(n_years):
        start_spr, end_spr = 6 * yr, 6 * yr + 3
        start_sum, end_sum = 6 * yr + 3, 6 * yr + 6
        
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

ths_c_interp, ths_e_interp = f3['theta_s_ctl_interp'].values, f3['theta_s_exp_interp'].values
tkd_c_interp, tkd_e_interp = f3['tkdry_ctl_interp'].values, f3['tkdry_exp_interp'].values
tksu_c_interp, tksu_e_interp = f3['tksatu_ctl_interp'].values, f3['tksatu_exp_interp'].values
tksf_c_interp, tksf_e_interp = f3['tksatf_ctl_interp'].values, f3['tksatf_exp_interp'].values
alf_c_interp, alf_e_interp = f3['BA_alpha_ctl_interp'].values, f3['BA_alpha_exp_interp'].values
bet_c_interp, bet_e_interp = f3['BA_beta_ctl_interp'].values, f3['BA_beta_exp_interp'].values
vs_c_interp, vs_e_interp = f3['vf_sand_s_ctl_interp'].values, f3['vf_sand_s_exp_interp'].values
vom_c_interp, vom_e_interp = f3['vf_om_s_ctl_interp'].values, f3['vf_om_s_exp_interp'].values
vgr_c_interp, vgr_e_interp = f3['vf_gravels_s_ctl_interp'].values, f3['vf_gravels_s_exp_interp'].values

# ========================= 3. 计算热传导通量 (Heat Flux) =========================
print("正在计算土壤热通量...")

def calc_Sr(th, ths):
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

# 先计算 10 层的热导率 k (对应公式中的 lambda_i)
k_ctl_spr, k_exp_spr = np.zeros((10, n_lat, n_lon)), np.zeros((10, n_lat, n_lon))
k_ctl_sum, k_exp_sum = np.zeros((10, n_lat, n_lon)), np.zeros((10, n_lat, n_lon))

for k_idx in range(10):
    p = 0 if k_idx <= 1 else (7 if k_idx >= 8 else k_idx - 1)
    
    # CTL
    Sr_c_spr = calc_Sr(theta_ctl_spr[k_idx], ths_c_interp[p])
    Sr_c_sum = calc_Sr(theta_ctl_sum[k_idx], ths_c_interp[p])
    k_ctl_spr[k_idx] = calc_k(t_ctl_spr[k_idx], Sr_c_spr, tkd_c_interp[p], tksu_c_interp[p], tksf_c_interp[p], vom_c_interp[p], vs_c_interp[p], vgr_c_interp[p], alf_c_interp[p], bet_c_interp[p])
    k_ctl_sum[k_idx] = calc_k(t_ctl_sum[k_idx], Sr_c_sum, tkd_c_interp[p], tksu_c_interp[p], tksf_c_interp[p], vom_c_interp[p], vs_c_interp[p], vgr_c_interp[p], alf_c_interp[p], bet_c_interp[p])
    
    # EXP
    Sr_e_spr = calc_Sr(theta_exp_spr[k_idx], ths_e_interp[p])
    Sr_e_sum = calc_Sr(theta_exp_sum[k_idx], ths_e_interp[p])
    k_exp_spr[k_idx] = calc_k(t_exp_spr[k_idx], Sr_e_spr, tkd_e_interp[p], tksu_e_interp[p], tksf_e_interp[p], vom_e_interp[p], vs_e_interp[p], vgr_e_interp[p], alf_e_interp[p], bet_e_interp[p])
    k_exp_sum[k_idx] = calc_k(t_exp_sum[k_idx], Sr_e_sum, tkd_e_interp[p], tksu_e_interp[p], tksf_e_interp[p], vom_e_interp[p], vs_e_interp[p], vgr_e_interp[p], alf_e_interp[p], bet_e_interp[p])

# 计算界面热导率和热通量 F_i
F_ctl_spr, F_exp_spr = np.zeros((10, n_lat, n_lon)), np.zeros((10, n_lat, n_lon))
F_ctl_sum, F_exp_sum = np.zeros((10, n_lat, n_lon)), np.zeros((10, n_lat, n_lon))

def calc_flux(lam_i, lam_i1, z_i, z_i1, zh_i, t_i, t_i1):
    # 计算公式 (10.10) 的 lambda[z_{h,i}]
    den = lam_i * (z_i1 - zh_i) + lam_i1 * (zh_i - z_i)
    den = np.where(den == 0, 1e-10, den) # 防止除以 0
    lam_int = (lam_i * lam_i1 * (z_i1 - z_i)) / den
    
    # 计算公式 (10.8) 的 F_i
    flux = lam_int * (t_i - t_i1) / (z_i1 - z_i)
    return flux

for i in range(9): # 对应第 1 层到第 9 层 (Python 索引 0-8)
    F_ctl_spr[i] = calc_flux(k_ctl_spr[i], k_ctl_spr[i+1], z[i], z[i+1], zh[i], t_ctl_spr[i], t_ctl_spr[i+1])
    F_ctl_sum[i] = calc_flux(k_ctl_sum[i], k_ctl_sum[i+1], z[i], z[i+1], zh[i], t_ctl_sum[i], t_ctl_sum[i+1])
    F_exp_spr[i] = calc_flux(k_exp_spr[i], k_exp_spr[i+1], z[i], z[i+1], zh[i], t_exp_spr[i], t_exp_spr[i+1])
    F_exp_sum[i] = calc_flux(k_exp_sum[i], k_exp_sum[i+1], z[i], z[i+1], zh[i], t_exp_sum[i], t_exp_sum[i+1])

# 第 10 层 (索引 9) 的界面热导率为 0，故通量为 0，数组初始化时已经是 0，无需额外计算

# 计算差值 (EXP - CTL)
diff_spr = F_exp_spr - F_ctl_spr
diff_sum = F_exp_sum - F_ctl_sum

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

diff_spr_ip = np.where(np.abs(diff_spr_ip) < 1e-6, np.nan, diff_spr_ip)
diff_sum_ip = np.where(np.abs(diff_sum_ip) < 1e-6, np.nan, diff_sum_ip)

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

# 注意：热通量的量级通常在几 W/m2 到十几 W/m2 之间。
# 此处预设色标范围为 -2 到 2，如果画出来的图过红或过蓝，请自行微调 levels_diff 的范围。
levels_diff = [-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]  # np.arange(-4.0, 4.5, 0.5) 
cmap_diff = plt.get_cmap('coolwarm') 
# colors_blue = cmap_base(np.linspace(0.0, 0.45, 6))
# colors_red  = cmap_base(np.linspace(0.55, 1.0, 6))

# cmap_diff = mcolors.ListedColormap(np.vstack((colors_blue, colors_red)))
# cmap_diff.set_under(cmap_base(0.0))
# cmap_diff.set_over(cmap_base(1.0))

norm_diff = mcolors.BoundaryNorm(levels_diff, cmap_diff.N)

for j in range(10):
    for i, season_data in enumerate([diff_spr_ip[j], diff_sum_ip[j]]):
        ax = axes[j, i]
        
        season_str = "MAM" if i == 0 else "JJA"
        letter = chr(97 + j * 2 + i)
        title = f"({letter}) {season_str} (n={j+1})"
        
        is_left_col = (i == 0)
        is_bottom_row = (j == 9)
        
        apply_style_and_shp(ax, title, is_left_col, is_bottom_row)
        
        # 为了美观，第 10 层全是 0 的情况可能会报错或画出纯色块，此处做一个简单的异常处理
        if not np.all(season_data == 0) and not np.all(np.isnan(season_data)):
            cf = ax.contourf(
                lon2d, lat2d, season_data,
                levels=levels_diff, cmap=cmap_diff, norm=norm_diff,
                extend='both',transform=ccrs.PlateCarree(), zorder=1
            )
        
        ax.text(0.98, 0.98, "EXP - CTL", transform=ax.transAxes, fontsize=9, 
                verticalalignment='top', horizontalalignment='right', 
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

# ====================== 6. 统一水平色标 ======================
cbar_ax = fig.add_axes([0.08, 0.02, 0.8, 0.008]) 
sm = cm.ScalarMappable(cmap=cmap_diff, norm=norm_diff)
sm.set_array([]) # Matplotlib 规定必须加这一行空数组

# 用新建的 sm 来画色标
cb = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal', 
                  spacing='uniform', extend='both')

cb.set_label('Downward Soil Heat Flux Difference (W m$^{-2}$)', fontsize=14, fontweight='bold')

# 强行打上正确的刻度标签，绝不允许它乱算
cb.set_ticks(levels_diff)
cb.ax.set_xticklabels([str(x) for x in levels_diff])
cb.ax.tick_params(labelsize=12)
# ====================== 7. 保存 ======================
save_name = "./figs/thermal_soil_heat_flux_diff.pdf"
os.makedirs("./figs", exist_ok=True)
plt.savefig(save_name, dpi=1200, bbox_inches='tight')
plt.close()

print(f"绘图完成保存至：{save_name}")