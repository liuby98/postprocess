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

warnings.filterwarnings('ignore', category=RuntimeWarning, message='Mean of empty slice')

# ========================= 1. 物理常数与设定 =========================
tfrz = 273.16
cpliq  = 4188.0
cpice  = 2117.27
dz = np.array([0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058])

# ========================= 2. 数据读取与处理 =========================
print("正在读取和处理数据...")
path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
f1 = xr.open_dataset(path + "colmoff_2001-2023_monmean_nogravel.nc")
f2 = xr.open_dataset(path + "colmoff_2001-2023_monmean_gravel.nc")

f1_djf = xr.open_dataset(path + "colmoff_2001-2017_DJF_nogravel.nc")
f2_djf = xr.open_dataset(path + "colmoff_2001-2017_DJF_gravel.nc")

f3 = xr.open_dataset("interp_soil_params.nc")
f_input = xr.open_dataset("/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01")

lat2d = f_input['XLAT'].isel(Time=0).values
lon2d = f_input['XLONG'].isel(Time=0).values
n_years = 17
n_lat, n_lon, n_layers = lat2d.shape[0], lat2d.shape[1], 10

def get_season_mean(da):
    """提取春夏季变量并计算多年平均"""
    var_mon = da.isel(soilsnow=slice(5, 15)).values
    spr_mean = np.zeros((n_layers, n_lat, n_lon))
    sum_mean = np.zeros((n_layers, n_lat, n_lon))
    for k in range(n_layers):
        for i in range(n_years):
            start_spr, end_spr = 6 * i, 6 * i + 3
            start_sum, end_sum = 6 * i + 3, 6 * i + 6
            spr_mean[k] += np.nanmean(var_mon[start_spr:end_spr, :, :, k], axis=0) / n_years
            sum_mean[k] += np.nanmean(var_mon[start_sum:end_sum, :, :, k], axis=0) / n_years
    return spr_mean, sum_mean

def get_djf_mean(ds, var):
    """提取冬季变量：每年截取3个月求均值，再求17年总平均"""
    var_mon = ds[var].isel(soilsnow=slice(5, 15)).values
    win_mean = np.zeros((n_layers, n_lat, n_lon))
    for k in range(n_layers):
        for i in range(n_years):
            start_win, end_win = 3 * i, 3 * i + 3
            win_mean[k] += np.nanmean(var_mon[start_win:end_win, :, :, k], axis=0) / n_years
    return win_mean

# 提取春夏季
t_c_spr, t_c_sum = get_season_mean(f1['f_t_soisno'])
t_e_spr, t_e_sum = get_season_mean(f2['f_t_soisno'])
w_c_spr, w_c_sum = get_season_mean(f1['f_wliq_soisno'])
w_e_spr, w_e_sum = get_season_mean(f2['f_wliq_soisno'])
i_c_spr, i_c_sum = get_season_mean(f1['f_wice_soisno'])
i_e_spr, i_e_sum = get_season_mean(f2['f_wice_soisno'])

# 提取冬季
t_c_win = get_djf_mean(f1_djf, 'f_t_soisno')
t_e_win = get_djf_mean(f2_djf, 'f_t_soisno')
w_c_win = get_djf_mean(f1_djf, 'f_wliq_soisno')
w_e_win = get_djf_mean(f2_djf, 'f_wliq_soisno')
i_c_win = get_djf_mean(f1_djf, 'f_wice_soisno')
i_e_win = get_djf_mean(f2_djf, 'f_wice_soisno')

def calc_Sr(th, ths):
    valid_mask = (ths != 0) & (~np.isnan(ths))
    Sr = np.zeros_like(th)
    Sr[valid_mask] = th[valid_mask] / ths[valid_mask]
    return np.clip(Sr, 1e-4, 1.0)

def calc_k(t, Sr, tkd, tksu, tksf, vom, vs, vgr, alf, bet):
    ksat = np.where(t > tfrz, tksu, tksf)
    base = (1.0 / (1.0 + np.exp(-bet * Sr)))**3 - ((1.0 - Sr)/2.0)**3
    base = np.maximum(base, 1e-6) 
    power = 0.5 * (1.0 + vom - alf*vs - vgr)
    ke_unf = (Sr ** power) * (base ** (1.0 - vom))
    ke_frz = Sr ** (1.0 + vom)
    ke = np.where(t > tfrz, ke_unf, ke_frz)
    return (ksat - tkd) * ke + tkd

diff_alpha_win = np.zeros((n_layers, n_lat, n_lon))
diff_alpha_spr = np.zeros((n_layers, n_lat, n_lon))
diff_alpha_sum = np.zeros((n_layers, n_lat, n_lon))

# ========================= 3. 计算热扩散率 =========================
print("正在计算热扩散率...")
for k in range(10):
    p = 0 if k <= 1 else (7 if k >= 8 else k - 1)
    
    th_c_win, th_c_spr, th_c_sum = f3['theta_ctl_win_re'].values[k], f3['theta_ctl_spr_re'].values[k], f3['theta_ctl_sum_re'].values[k]
    th_e_win, th_e_spr, th_e_sum = f3['theta_exp_win_re'].values[k], f3['theta_exp_spr_re'].values[k], f3['theta_exp_sum_re'].values[k]
    
    ths_c, ths_e = f3['theta_s_ctl_interp'].values[p], f3['theta_s_exp_interp'].values[p]
    tkd_c, tkd_e = f3['tkdry_ctl_interp'].values[p], f3['tkdry_exp_interp'].values[p]
    tksu_c, tksu_e = f3['tksatu_ctl_interp'].values[p], f3['tksatu_exp_interp'].values[p]
    tksf_c, tksf_e = f3['tksatf_ctl_interp'].values[p], f3['tksatf_exp_interp'].values[p]
    vom_c, vom_e = f3['vf_om_s_ctl_interp'].values[p], f3['vf_om_s_exp_interp'].values[p]
    vs_c, vs_e = f3['vf_sand_s_ctl_interp'].values[p], f3['vf_sand_s_exp_interp'].values[p]
    vgr_c, vgr_e = f3['vf_gravels_s_ctl_interp'].values[p], f3['vf_gravels_s_exp_interp'].values[p]
    alf_c, alf_e = f3['BA_alpha_ctl_interp'].values[p], f3['BA_alpha_exp_interp'].values[p]
    bet_c, bet_e = f3['BA_beta_ctl_interp'].values[p], f3['BA_beta_exp_interp'].values[p]

    Sr_c_win, Sr_c_spr, Sr_c_sum = calc_Sr(th_c_win, ths_c), calc_Sr(th_c_spr, ths_c), calc_Sr(th_c_sum, ths_c)
    Sr_e_win, Sr_e_spr, Sr_e_sum = calc_Sr(th_e_win, ths_e), calc_Sr(th_e_spr, ths_e), calc_Sr(th_e_sum, ths_e)
    
    k_c_win = calc_k(t_c_win[k], Sr_c_win, tkd_c, tksu_c, tksf_c, vom_c, vs_c, vgr_c, alf_c, bet_c)
    k_e_win = calc_k(t_e_win[k], Sr_e_win, tkd_e, tksu_e, tksf_e, vom_e, vs_e, vgr_e, alf_e, bet_e)
    k_c_spr = calc_k(t_c_spr[k], Sr_c_spr, tkd_c, tksu_c, tksf_c, vom_c, vs_c, vgr_c, alf_c, bet_c)
    k_e_spr = calc_k(t_e_spr[k], Sr_e_spr, tkd_e, tksu_e, tksf_e, vom_e, vs_e, vgr_e, alf_e, bet_e)
    k_c_sum = calc_k(t_c_sum[k], Sr_c_sum, tkd_c, tksu_c, tksf_c, vom_c, vs_c, vgr_c, alf_c, bet_c)
    k_e_sum = calc_k(t_e_sum[k], Sr_e_sum, tkd_e, tksu_e, tksf_e, vom_e, vs_e, vgr_e, alf_e, bet_e)

    c_sol_c, c_sol_e = f3['csol_ctl_interp'].values[p], f3['csol_exp_interp'].values[p]
    
    wc_win, ic_win = np.nan_to_num(w_c_win[k]), np.nan_to_num(i_c_win[k])
    we_win, ie_win = np.nan_to_num(w_e_win[k]), np.nan_to_num(i_e_win[k])
    wc_spr, ic_spr = np.nan_to_num(w_c_spr[k]), np.nan_to_num(i_c_spr[k])
    we_spr, ie_spr = np.nan_to_num(w_e_spr[k]), np.nan_to_num(i_e_spr[k])
    wc_sum, ic_sum = np.nan_to_num(w_c_sum[k]), np.nan_to_num(i_c_sum[k])
    we_sum, ie_sum = np.nan_to_num(w_e_sum[k]), np.nan_to_num(i_e_sum[k])

    c_c_win = c_sol_c + (wc_win / dz[k]) * cpliq + (ic_win / dz[k]) * cpice
    c_e_win = c_sol_e + (we_win / dz[k]) * cpliq + (ie_win / dz[k]) * cpice
    c_c_spr = c_sol_c + (wc_spr / dz[k]) * cpliq + (ic_spr / dz[k]) * cpice
    c_e_spr = c_sol_e + (we_spr / dz[k]) * cpliq + (ie_spr / dz[k]) * cpice
    c_c_sum = c_sol_c + (wc_sum / dz[k]) * cpliq + (ic_sum / dz[k]) * cpice
    c_e_sum = c_sol_e + (we_sum / dz[k]) * cpliq + (ie_sum / dz[k]) * cpice

    alpha_c_win = (k_c_win / c_c_win) * 1e6
    alpha_e_win = (k_e_win / c_e_win) * 1e6
    alpha_c_spr = (k_c_spr / c_c_spr) * 1e6
    alpha_e_spr = (k_e_spr / c_e_spr) * 1e6
    alpha_c_sum = (k_c_sum / c_c_sum) * 1e6
    alpha_e_sum = (k_e_sum / c_e_sum) * 1e6

    diff_alpha_win[k] = alpha_e_win - alpha_c_win
    diff_alpha_spr[k] = alpha_e_spr - alpha_c_spr
    diff_alpha_sum[k] = alpha_e_sum - alpha_c_sum

# ========================= 4. 观测掩膜处理 =========================
print("正在进行掩膜处理...")
f_obs = xr.open_dataset("/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc")
lon1d_mesh, lat1d_mesh = np.meshgrid(f_obs['lon'].values, f_obs['lat'].values)
t_mask_2d = griddata(
    (lon1d_mesh.flatten(), lat1d_mesh.flatten()), 
    f_obs['tm'].isel(time=0).values.flatten(), 
    (lon2d, lat2d), method='nearest'
)

t_mask_3d = np.broadcast_to(t_mask_2d, diff_alpha_spr.shape)
diff_win_ip = np.where(np.isnan(t_mask_3d), np.nan, diff_alpha_win)
diff_spr_ip = np.where(np.isnan(t_mask_3d), np.nan, diff_alpha_spr)
diff_sum_ip = np.where(np.isnan(t_mask_3d), np.nan, diff_alpha_sum)

# ========================= 5. 绘图 =========================
print("正在绘图...")
os.environ['CARTOPY_OFFLINE'] = 'true'
os.environ['PROJ_NETWORK'] = 'OFF'

COASTLINE_50M = cfeature.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='gray', facecolor='none')
LAKES_50M = cfeature.NaturalEarthFeature('physical', 'lakes', '50m', edgecolor='gray', facecolor='none')

proj = ccrs.LambertConformal(central_longitude=110, central_latitude=40, standard_parallels=(30, 60))
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

levels_diff = np.arange(-0.3, 0.31, 0.01)
cmap_diff = plt.get_cmap('RdBu_r') 

for j in range(10):
    for i, season_data in enumerate([diff_win_ip[j], diff_spr_ip[j], diff_sum_ip[j]]):
        ax = axes[j, i]
        
        season_str = ["DJF", "MAM", "JJA"][i]
        title = f"{season_str} (n={j+1})"
        
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
cb.set_label('Thermal Diffusivity Difference ($10^{-6}$ m$^{2}$ s$^{-1}$)', fontsize=14, fontweight='bold')
cb.ax.tick_params(labelsize=12)

# ====================== 7. 保存 ======================
save_name = "./figs/thermal_diffusivity_diff.pdf"
os.makedirs("./figs", exist_ok=True)
plt.savefig(save_name, dpi=1200, bbox_inches='tight')
plt.close()

print(f"绘图完成保存至：{save_name}")