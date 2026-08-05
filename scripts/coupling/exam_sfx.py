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
ctl_file = path + "add_colmrun_2001-2017_monmean_nogravel.nc"
exp_file = path + "add_colmrun_2001-2017_monmean_gravel.nc"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
obs_mask_file = "/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc"
shp_dir = "../../shapefile_China/"

# 需要绘制的感热通量变量列表
variables_to_plot = ['f_fsena', 'f_fsenl', 'f_fseng']

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

    # 2001-2017年共17年，每年6个月（3至8月），总计102个时间步
    # 按照每年的月份顺序：3, 4, 5, 6, 7, 8月
    nyears = 17
    nmonths = 6
    spring_idx = [nmonths * i + j for i in range(nyears) for j in range(3)]     # 3, 4, 5月
    summer_idx = [nmonths * i + j for i in range(nyears) for j in range(3, 6)]  # 6, 7, 8月

    # ====================== 4. 数据处理与循环绘图 ======================
    nc_ctl = Dataset(ctl_file)
    nc_exp = Dataset(exp_file)

    os.makedirs("../illustration/cpl", exist_ok=True)
    os.environ['CARTOPY_OFFLINE'] = 'true'

    proj = ccrs.LambertConformal(
        central_longitude=PROJ_CENTRAL_LON, 
        central_latitude=PROJ_CENTRAL_LAT, 
        standard_parallels=PROJ_STD_PARALLELS
    )

    def process_and_plot(var_key):
        print(f"正在处理并绘制变量: {var_key} ...")
        
        def extract_and_calc(indices):
            """在指定时间步索引上提取数据，并在原始时间步上完成处理"""
            # 直接读取对应的通量变量
            c_raw = nc_ctl.variables[var_key][indices, :, :]
            e_raw = nc_exp.variables[var_key][indices, :, :]

            # 处理缺测值
            if np.ma.isMaskedArray(c_raw): c_raw = c_raw.filled(np.nan)
            if np.ma.isMaskedArray(e_raw): e_raw = e_raw.filled(np.nan)
            c_raw[c_raw < -1e30] = np.nan
            e_raw[e_raw < -1e30] = np.nan
            
            # 对原始时间步计算多年季节平均态
            c_mean = np.nanmean(c_raw, axis=0)
            e_mean = np.nanmean(e_raw, axis=0)
            d_mean = e_mean - c_mean 

            # 插值到高分辨率观测网格
            c_ip = regrid_rcm2rgrid(c_mean, lat2d, lon2d, lat1d, lon1d)
            e_ip = regrid_rcm2rgrid(e_mean, lat2d, lon2d, lat1d, lon1d)
            d_ip = regrid_rcm2rgrid(d_mean, lat2d, lon2d, lat1d, lon1d)

            # 应用观测数据陆地掩码
            return (
                np.where(mask_obs, c_ip, np.nan),
                np.where(mask_obs, e_ip, np.nan),
                np.where(mask_obs, d_ip, np.nan)
            )

        # 获取春秋季结果
        c_spring, e_spring, d_spring = extract_and_calc(spring_idx)
        c_summer, e_summer, d_summer = extract_and_calc(summer_idx)

        # ------------------ 开始绘图 (2x3 面板) ------------------
        fig, axes = plt.subplots(2, 3, figsize=(18, 12), subplot_kw={'projection': proj})
        plt.subplots_adjust(wspace=0.1, hspace=0.15)
        
        plot_data = [
            [c_spring, e_spring, d_spring],
            [c_summer, e_summer, d_summer]
        ]
        
        titles = [
            [f'Spring CTL ({var_key})', f'Spring EXP ({var_key})', f'Spring (EXP - CTL)'],
            [f'Summer CTL ({var_key})', f'Summer EXP ({var_key})', f'Summer (EXP - CTL)']
        ]

        # 动态设置色标范围（利用百分位数防止极值拉平颜色）
        all_vals = np.concatenate([c_spring.ravel(), e_spring.ravel(), c_summer.ravel(), e_summer.ravel()])
        v_min, v_max = np.nanpercentile(all_vals, 2), np.nanpercentile(all_vals, 98)
        mean_levels = np.linspace(v_min, v_max, 56)
        
        all_diffs = np.concatenate([d_spring.ravel(), d_summer.ravel()])
        d_absmax = np.nanpercentile(np.abs(all_diffs), 98)
        if d_absmax == 0 or np.isnan(d_absmax): d_absmax = 5.0  # 感热通量差异可能稍大，给个保底值5.0
        diff_levels = np.linspace(-d_absmax, d_absmax, 61)

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
                gl.left_labels = (j == 0)
                gl.bottom_labels = (i == 1) 
                gl.xlabel_style = {'size': 10, 'rotation': 0, 'va': 'top', 'ha': 'center'}
                gl.ylabel_style = {'size': 10, 'rotation': 0, 'va': 'center', 'ha': 'right'}

                data = plot_data[i][j]
                ax.set_title(titles[i][j], loc='center', fontsize=14, pad=10)

                # 填色图绘制
                if j < 2:  # CTL 或 EXP
                    cf_mean = ax.contourf(lon_grid, lat_grid, data, 
                                          levels=mean_levels, cmap='Spectral_r', 
                                          extend='both', transform=ccrs.PlateCarree(), zorder=1)
                else:      # DIFF (EXP-CTL)
                    cf_diff = ax.contourf(lon_grid, lat_grid, data, 
                                          levels=diff_levels, cmap='RdBu_r', 
                                          extend='both', transform=ccrs.PlateCarree(), zorder=1)

        # Colorbar 配置
        cbar_ax_mean = fig.add_axes([0.15, 0.05, 0.4, 0.02]) 
        cb_mean = plt.colorbar(cf_mean, cax=cbar_ax_mean, orientation='horizontal', extend='both')
        cb_mean.ax.tick_params(labelsize=11)
        cb_mean.set_label(f'{var_key} Mean Flux (W/m²)', fontsize=12)

        cbar_ax_diff = fig.add_axes([0.68, 0.05, 0.2, 0.02])
        cb_diff = plt.colorbar(cf_diff, cax=cbar_ax_diff, orientation='horizontal', extend='both')
        cb_diff.ax.tick_params(labelsize=11)
        cb_diff.set_label(f'{var_key} Difference (EXP - CTL)', fontsize=12)

        # 保存图片
        save_name = f"./{var_key}_spring_summer_diff.png"
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"{var_key} 对比图已成功导出至: {save_name}")

    # 循环生成 fsena, fsenl, fseng 的对比图
    for var in variables_to_plot:
        process_and_plot(var)

    nc_ctl.close()
    nc_exp.close()
    print("全部变量绘制完成！")