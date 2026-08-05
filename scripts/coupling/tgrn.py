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

# 路径配置
path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
ctl_file = path + "colmrun_2001-2023_monmean_nogravel.nc"
exp_file = path + "colmrun_2001-2023_monmean_gravel.nc"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
obs_mask_file = "/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc"
shp_dir = "../../shapefile_China/"

# 变量配置 (仅包含地表温度)
var_name = 'f_t_grnd'

# 差异场等值线级别 (-1.5℃ 到 1.5℃)
diff_levels = np.linspace(-1.5, 1.5, 61)
diff_cmap = 'RdBu_r' 

# 绝对温度(平均态)的等值线范围配置 (转换为摄氏度: -15℃ 到 40℃)
mean_levels = np.linspace(-15, 40, 56)
mean_cmap = 'Spectral_r'

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

if __name__ == '__main__':
    # ====================== 3. 读取网格与掩码数据 ======================
    with Dataset(wrfinput_file) as f_wrf:
        lat2d = f_wrf.variables['XLAT'][0, :, :]
        lon2d = f_wrf.variables['XLONG'][0, :, :]

    with Dataset(obs_mask_file) as f_obs:
        lat1d = f_obs.variables['lat'][:]
        lon1d = f_obs.variables['lon'][:]
        tm_obs = f_obs.variables['tm'][0, :, :]
        mask_obs = ~np.ma.getmaskarray(tm_obs) & ~np.isnan(tm_obs)

    lon_grid, lat_grid = np.meshgrid(lon1d, lat1d)

    # 每年提取6个月，分离出 2001-2017年的春季(索引0,1,2)和夏季(索引3,4,5)
    nyears = 17
    nmonths = 6
    spring_idx = [nmonths * i + j for i in range(nyears) for j in range(3)]
    summer_idx = [nmonths * i + j for i in range(nyears) for j in range(3, 6)]

    # ====================== 4. 数据处理 ======================
    print(f"处理变量: {var_name} [Spring & Summer Panels - Converted to ℃ - No Significance Test]...")
    
    nc_ctl = Dataset(ctl_file)
    nc_exp = Dataset(exp_file)

    def process_season(indices):
        """处理单一季节，返回转换为℃后的插值CTL, EXP, 原始Diff网格数据"""
        c_raw = nc_ctl.variables[var_name][indices, :, :]
        e_raw = nc_exp.variables[var_name][indices, :, :]
        
        if np.ma.isMaskedArray(c_raw): c_raw = c_raw.filled(np.nan)
        if np.ma.isMaskedArray(e_raw): e_raw = e_raw.filled(np.nan)
        c_raw[c_raw < -1e30] = np.nan
        e_raw[e_raw < -1e30] = np.nan

        # 将单位从 Kelvin 转换为摄氏度 ℃
        c_raw = c_raw - 273.15
        e_raw = e_raw - 273.15

        nlat, nlon = c_raw.shape[1], c_raw.shape[2]
        
        # 计算多年气候态均值 (直接在所有选定时间步上求平均)
        c_mean = np.nanmean(c_raw, axis=0)
        e_mean = np.nanmean(e_raw, axis=0)
        
        # 计算差异
        d_mean = e_mean - c_mean 

        # 插值到高分辨率观测网格
        c_ip = regrid_rcm2rgrid(c_mean, lat2d, lon2d, lat1d, lon1d)
        e_ip = regrid_rcm2rgrid(e_mean, lat2d, lon2d, lat1d, lon1d)
        d_ip = regrid_rcm2rgrid(d_mean, lat2d, lon2d, lat1d, lon1d)

        # 应用观测数据陆地掩码
        c_masked = np.where(mask_obs, c_ip, np.nan)
        e_masked = np.where(mask_obs, e_ip, np.nan)
        d_masked = np.where(mask_obs, d_ip, np.nan)
        
        # 不做检验，直接返回原始差异场
        return c_masked, e_masked, d_masked

    # 分别获取春秋结果
    c_spring, e_spring, d_spring = process_season(spring_idx)
    c_summer, e_summer, d_summer = process_season(summer_idx)
    
    nc_ctl.close()
    nc_exp.close()

    # ====================== 5. 绘图环境准备 ======================
    os.makedirs("../illustration/cpl", exist_ok=True)
    os.environ['CARTOPY_OFFLINE'] = 'true'

    proj = ccrs.LambertConformal(
        central_longitude=PROJ_CENTRAL_LON, 
        central_latitude=PROJ_CENTRAL_LAT, 
        standard_parallels=PROJ_STD_PARALLELS
    )

    # 调整 2x3 组图的画布大小
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), subplot_kw={'projection': proj})
    plt.subplots_adjust(wspace=0.1, hspace=0.15)
    
    plot_data = [
        [c_spring, e_spring, d_spring],
        [c_summer, e_summer, d_summer]
    ]
    
    titles = [
        ['Spring CTL', 'Spring EXP', 'Spring (EXP - CTL)'],
        ['Summer CTL', 'Summer EXP', 'Summer (EXP - CTL)']
    ]

    # ====================== 6. 循环绘制面板 ======================
    for i in range(2):      # 行 (0: Spring, 1: Summer)
        for j in range(3):  # 列 (0: CTL, 1: EXP, 2: DIFF)
            ax = axes[i, j]
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
            
            # 添加经纬网格线
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                              x_inline=False, y_inline=False, 
                              linewidth=0.6, color='gray', alpha=0.5, linestyle='--', zorder=2)
            gl.xlocator = mticker.FixedLocator(np.arange(70, 135, 10))
            gl.ylocator = mticker.FixedLocator(np.arange(10, 55, 10))
            gl.top_labels = gl.right_labels = False
            
            # 仅在最底层/最左侧显示标签，避免重叠
            gl.left_labels = (j == 0)
            gl.bottom_labels = (i == 1) 
            gl.xlabel_style = {'size': 10, 'rotation': 0, 'va': 'top', 'ha': 'center'}
            gl.ylabel_style = {'size': 10, 'rotation': 0, 'va': 'center', 'ha': 'right'}

            data = plot_data[i][j]
            ax.set_title(titles[i][j], loc='center', fontsize=14, pad=10)

            # 区分均值场和差异场使用不同的 colormap 和 levels
            if j < 2:  # CTL 或 EXP
                cf_mean = ax.contourf(lon_grid, lat_grid, data, 
                                      levels=mean_levels, cmap=mean_cmap, 
                                      extend='both', transform=ccrs.PlateCarree(), zorder=1)
            else:      # DIFF (EXP-CTL) (原始填色，无检验屏蔽)
                cf_diff = ax.contourf(lon_grid, lat_grid, data, 
                                      levels=diff_levels, cmap=diff_cmap, 
                                      vmin=-1.5, vmax=1.5,
                                      extend='both', transform=ccrs.PlateCarree(), zorder=1)

    # ====================== 7. Colorbar 配置 ======================
    # 为前两列共享绝对温度均值色条 (放置于下方，单位 ℃)
    cbar_ax_mean = fig.add_axes([0.15, 0.05, 0.4, 0.02]) 
    cb_mean = plt.colorbar(cf_mean, cax=cbar_ax_mean, orientation='horizontal', extend='both')
    cb_mean.ax.tick_params(labelsize=11)
    cb_mean.set_label('Ground Surface Temperature (℃)', fontsize=12)

    # 为第三列配置差异色条 (放置于下方，单位 ℃)
    cbar_ax_diff = fig.add_axes([0.68, 0.05, 0.2, 0.02])
    cb_ticks_diff = np.linspace(-1.5, 1.5, 7)
    cb_diff = plt.colorbar(cf_diff, cax=cbar_ax_diff, orientation='horizontal', extend='both', ticks=cb_ticks_diff)
    cb_diff.ax.tick_params(labelsize=11)
    cb_diff.set_label('Difference (℃)', fontsize=12)

    # 保存图片
    save_name = "./ground_temp_diff.png"
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"2x3 地表温度多面板对比图(原始未检验差异)已成功导出至: {save_name}")