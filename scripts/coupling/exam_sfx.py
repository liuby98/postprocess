import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import warnings
warnings.filterwarnings("ignore")
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.ticker as mticker
import os

# ====================== 1. 全局配置 ======================
PLOT_EXTENT = [80, 130, 10, 55]
PROJ_TYPE = 'LambertConformal' 
PROJ_CENTRAL_LON = 110.0        
PROJ_CENTRAL_LAT = 40.0         
PROJ_STD_PARALLELS = (30.0, 60.0) 

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 基础数据路径设置
path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
obs_mask_file = "/share/home/dq135/reference/CN05.1_Tm_1991_2023_JJA_025x025.nc"
shp_dir = "../../shapefile_China/"

# 绘制变量
variables_to_plot = ['AHFX']

# 时间轴相关配置
nyears = 17
ndays_per_season = 92  
ntimes_per_day = 4
hours = ['00 UTC', '06 UTC', '12 UTC', '18 UTC']

# 季节列表
seasons = ['MAM', 'JJA']

# ========== 新增控制开关 ==========
OUTPUT_MONTHLY = False   # 设为 True 则额外输出各月平均图

# 定义各季节内月份的天数切片（基于公历）
season_month_slices = {
    'MAM': [('Mar', slice(0, 31)), ('Apr', slice(31, 61)), ('May', slice(61, 92))],
    'JJA': [('Jun', slice(0, 30)), ('Jul', slice(30, 61)), ('Aug', slice(61, 92))]
}

# ====================== 2. 核心函数 ======================
def regrid_rcm2rgrid(var2d, lat2d, lon2d, lat1d, lon1d):
    pts = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    vals = var2d.ravel()
    gx, gy = np.meshgrid(lon1d, lat1d)
    interp = griddata(pts, vals, (gx, gy), method='linear')
    nan_mask = np.isnan(interp)
    if nan_mask.any():
        interp[nan_mask] = griddata(pts, vals, (gx, gy), method='nearest')[nan_mask]
    return interp

# 辅助函数：绘制单个子图
def draw_single_panel(ax, data, title, is_diff, lon_grid, lat_grid, mean_levels, diff_levels, show_left_label, show_bottom_label):
    ax.set_extent(PLOT_EXTENT, crs=ccrs.PlateCarree())
    
    # 添加 Shapefile 边界
    for shp_file, lw in [("province.shp", 0.4), ("china.shp", 0.6), ("south_china_sea.shp", 0.8)]:
        p_shp = os.path.join(shp_dir, shp_file)
        if os.path.exists(p_shp):
            color = 'blue' if 'river' in shp_file else 'black'
            try:
                reader = shpreader.Reader(p_shp)
                ax.add_geometries(reader.geometries(), crs=ccrs.PlateCarree(),
                                  facecolor='none', edgecolor=color, linewidth=lw, zorder=3)
            except Exception:
                pass
    
    # 添加网格线
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                      x_inline=False, y_inline=False, 
                      linewidth=0.6, color='gray', alpha=0.5, linestyle='--', zorder=2)
    gl.xlocator = mticker.FixedLocator(np.arange(70, 135, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(10, 55, 10))
    gl.top_labels = gl.right_labels = False
    gl.left_labels = show_left_label
    gl.bottom_labels = show_bottom_label
    gl.xlabel_style = {'size': 10, 'rotation': 0, 'va': 'top', 'ha': 'center'}
    gl.ylabel_style = {'size': 10, 'rotation': 0, 'va': 'center', 'ha': 'right'}
    ax.set_title(title, loc='center', fontsize=14, pad=10)

    # 统一使用 'RdBu_r' 色板
    levels = diff_levels if is_diff else mean_levels
    cf = ax.contourf(lon_grid, lat_grid, data, levels=levels, cmap='RdBu_r', 
                     extend='both', transform=ccrs.PlateCarree(), zorder=1)
    return cf

# 辅助函数：统一添加和格式化 Colorbar
def add_formatted_colorbar(fig, cf, rect, label):
    cbar_ax = fig.add_axes(rect)
    cb = plt.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='both', format='%.1f')
    cb.locator = mticker.MaxNLocator(nbins=6)
    cb.update_ticks()
    cb.ax.tick_params(labelsize=11, rotation=15)
    cb.set_label(label, fontsize=12)
    return cb

if __name__ == '__main__':
    # ====================== 3. 读取网格与掩码数据 ======================
    with Dataset(wrfinput_file) as f_wrf:
        lat2d = f_wrf.variables['XLAT'][0, :, :]
        lon2d = f_wrf.variables['XLONG'][0, :, :]
        nlat, nlon = lat2d.shape

    with Dataset(obs_mask_file) as f_obs:
        lat1d = f_obs.variables['lat'][:]
        lon1d = f_obs.variables['lon'][:]
        tm_obs = f_obs.variables['tm'][0, :, :]
        mask_obs = ~np.ma.getmaskarray(tm_obs) & ~np.isnan(tm_obs)

    lon_grid, lat_grid = np.meshgrid(lon1d, lat1d)

    os.makedirs("./figs", exist_ok=True)
    os.environ['CARTOPY_OFFLINE'] = 'true'

    proj = ccrs.LambertConformal(
        central_longitude=PROJ_CENTRAL_LON, 
        central_latitude=PROJ_CENTRAL_LAT, 
        standard_parallels=PROJ_STD_PARALLELS
    )

    # ====================== 4. 跨季节数据处理与绘图 ======================
    # 固定组图的色板范围
    mean_levels = np.arange(-200, 205, 5)   # 原始分布：±200，间隔5
    diff_levels = np.arange(-15, 15.5, 0.5) # 差异分布：±15，间隔0.5

    for var_key in variables_to_plot:
        print(f"========== 正在处理变量: {var_key} ==========")
        
        for season in seasons:
            print(f"\n>> 开始处理 {season} 季节数据 (各时刻多年平均)...")
            ctl_file = os.path.join(path, f"wrfout_2001-2017_{season}_4monents_nogravel.nc")
            exp_file = os.path.join(path, f"wrfout_2001-2017_{season}_4monents_gravel.nc")

            with Dataset(ctl_file) as nc_ctl, Dataset(exp_file) as nc_exp:
                c_raw = nc_ctl.variables[var_key][:]
                e_raw = nc_exp.variables[var_key][:]
            
            if np.ma.isMaskedArray(c_raw): c_raw = c_raw.filled(np.nan)
            if np.ma.isMaskedArray(e_raw): e_raw = e_raw.filled(np.nan)
            c_raw[c_raw < -1e30] = np.nan
            e_raw[e_raw < -1e30] = np.nan

            c_data = c_raw.reshape(nyears, ndays_per_season, ntimes_per_day, nlat, nlon)
            e_data = e_raw.reshape(nyears, ndays_per_season, ntimes_per_day, nlat, nlon)

            # ----- 计算当前季节各个时刻的多年均值（季节平均） -----
            c_season_moments = np.nanmean(np.nanmean(c_data, axis=1), axis=0)
            e_season_moments = np.nanmean(np.nanmean(e_data, axis=1), axis=0)
            d_season_moments = e_season_moments - c_season_moments

            # ===== 季节平均绘图（原代码） =====
            c_plot_data, e_plot_data, d_plot_data = [], [], []
            for t in range(4):
                c_plot_data.append(np.where(mask_obs, regrid_rcm2rgrid(c_season_moments[t], lat2d, lon2d, lat1d, lon1d), np.nan))
                e_plot_data.append(np.where(mask_obs, regrid_rcm2rgrid(e_season_moments[t], lat2d, lon2d, lat1d, lon1d), np.nan))
                d_plot_data.append(np.where(mask_obs, regrid_rcm2rgrid(d_season_moments[t], lat2d, lon2d, lat1d, lon1d), np.nan))
            
            print(f"  正在生成 {season} 的 4个时刻组图（季节平均）...")
            fig_m, axes_m = plt.subplots(4, 3, figsize=(18, 20), subplot_kw={'projection': proj})
            plt.subplots_adjust(wspace=0.1, hspace=0.15)
            
            cf_mean_m, cf_diff_m = None, None
            for i in range(4):     
                for j in range(3): 
                    ax = axes_m[i, j]
                    is_diff = (j == 2)
                    
                    if j == 0:
                        data, title = c_plot_data[i], f'{season} {hours[i]} CTL ({var_key})'
                    elif j == 1:
                        data, title = e_plot_data[i], f'{season} {hours[i]} EXP ({var_key})'
                    else:
                        data, title = d_plot_data[i], f'{season} {hours[i]} (EXP - CTL)'

                    cf = draw_single_panel(ax, data, title, is_diff, lon_grid, lat_grid, 
                                           mean_levels, diff_levels, show_left_label=(j==0), show_bottom_label=(i==3))
                    if j == 1: cf_mean_m = cf
                    if j == 2: cf_diff_m = cf

            add_formatted_colorbar(fig_m, cf_mean_m, [0.15, 0.08, 0.4, 0.015], f'{var_key} Mean Flux (W/m²)')
            add_formatted_colorbar(fig_m, cf_diff_m, [0.68, 0.08, 0.2, 0.015], f'{var_key} Difference (EXP - CTL)')

            plt.savefig(f"./figs/{var_key}_{season}_4moments_Mean.png", dpi=300, bbox_inches='tight')
            plt.close()

            # ===== 若开关开启，额外输出每月平均图 =====
            if OUTPUT_MONTHLY:
                month_list = season_month_slices.get(season, [])
                for month_name, month_slice in month_list:
                    print(f"    > 正在处理 {season} 内 {month_name} 月平均...")
                    # 从完整季节数据中切出该月份的天数
                    c_month = c_data[:, month_slice, :, :, :]   # (nyears, ndays_month, ntimes, nlat, nlon)
                    e_month = e_data[:, month_slice, :, :, :]
                    # 多年平均（先平均年份，再平均天数）
                    c_month_moments = np.nanmean(np.nanmean(c_month, axis=1), axis=0)
                    e_month_moments = np.nanmean(np.nanmean(e_month, axis=1), axis=0)
                    d_month_moments = e_month_moments - c_month_moments

                    # 插值并应用掩码
                    c_plot_m, e_plot_m, d_plot_m = [], [], []
                    for t in range(4):
                        c_plot_m.append(np.where(mask_obs, regrid_rcm2rgrid(c_month_moments[t], lat2d, lon2d, lat1d, lon1d), np.nan))
                        e_plot_m.append(np.where(mask_obs, regrid_rcm2rgrid(e_month_moments[t], lat2d, lon2d, lat1d, lon1d), np.nan))
                        d_plot_m.append(np.where(mask_obs, regrid_rcm2rgrid(d_month_moments[t], lat2d, lon2d, lat1d, lon1d), np.nan))

                    # 绘图（结构与季节图一致）
                    fig_mo, axes_mo = plt.subplots(4, 3, figsize=(18, 20), subplot_kw={'projection': proj})
                    plt.subplots_adjust(wspace=0.1, hspace=0.15)
                    cf_mean_mo, cf_diff_mo = None, None
                    for i in range(4):
                        for j in range(3):
                            ax = axes_mo[i, j]
                            is_diff = (j == 2)
                            if j == 0:
                                data, title = c_plot_m[i], f'{month_name} {hours[i]} CTL ({var_key})'
                            elif j == 1:
                                data, title = e_plot_m[i], f'{month_name} {hours[i]} EXP ({var_key})'
                            else:
                                data, title = d_plot_m[i], f'{month_name} {hours[i]} (EXP - CTL)'
                            cf = draw_single_panel(ax, data, title, is_diff, lon_grid, lat_grid,
                                                   mean_levels, diff_levels, show_left_label=(j==0), show_bottom_label=(i==3))
                            if j == 1: cf_mean_mo = cf
                            if j == 2: cf_diff_mo = cf

                    add_formatted_colorbar(fig_mo, cf_mean_mo, [0.15, 0.08, 0.4, 0.015], f'{var_key} Mean Flux (W/m²)')
                    add_formatted_colorbar(fig_mo, cf_diff_mo, [0.68, 0.08, 0.2, 0.015], f'{var_key} Difference (EXP - CTL)')

                    plt.savefig(f"./figs/{var_key}_{season}_{month_name}_4moments_Mean.png", dpi=300, bbox_inches='tight')
                    plt.close()

    print("\n所有绘图任务顺利执行完毕！")