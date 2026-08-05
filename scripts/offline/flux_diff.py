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
ctl_file = path + "colmoff_2001-2023_monmean_nogravel.nc"
exp_file = path + "colmoff_2001-2023_monmean_gravel.nc"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
obs_mask_file = "/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc"
shp_dir = "../../shapefile_China/"

# 变量配置 (2行5列)
seasons = ['MAM', 'JJA']
season_names = ['Spring', 'Summer']
variables = ['f_fgrnd', 'f_olrg', 'f_rnet', 'f_fsena', 'f_lfevpa']
var_titles = ['Ground Heat Flux', 'Outgoing Longwave Radiation', 
              'Net Radiation', 'Sensible Heat Flux', 'Latent Heat Flux']

# 设置colorbar层级及colormap (参考NCL的 -3 到 3 和 MPL_BrBG)
diff_levels = np.linspace(-3, 3, 21)
diff_cmap = 'PRGn' # NCL中的MPL_BrBG中间是白色，两端棕蓝/绿

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
    # WRF 经纬度
    with Dataset(wrfinput_file) as f_wrf:
        lat2d = f_wrf.variables['XLAT'][0, :, :]
        lon2d = f_wrf.variables['XLONG'][0, :, :]

    # CN05.1 (观测) 经纬度及陆地掩码
    with Dataset(obs_mask_file) as f_obs:
        lat1d = f_obs.variables['lat'][:]
        lon1d = f_obs.variables['lon'][:]
        tm_obs = f_obs.variables['tm'][0, :, :]
        mask_obs = ~np.ma.getmaskarray(tm_obs) & ~np.isnan(tm_obs)

    lon_grid, lat_grid = np.meshgrid(lon1d, lat1d)

    # 构造时间掩码 (假设17年，每年6个月数据，0-2为春季，3-5为夏季)
    nyears = 17
    mam_idx = [6 * i + j for i in range(nyears) for j in range(3)]
    jja_idx = [6 * i + j + 3 for i in range(nyears) for j in range(3)]

    # ====================== 4. 绘图环境准备 ======================
    os.makedirs("../illustration/off", exist_ok=True)
    os.environ['CARTOPY_OFFLINE'] = 'true'

    proj = ccrs.LambertConformal(
        central_longitude=PROJ_CENTRAL_LON, 
        central_latitude=PROJ_CENTRAL_LAT, 
        standard_parallels=PROJ_STD_PARALLELS
    )

    # 创建 2行 x 5列 画布
    fig = plt.figure(figsize=(24, 8.5)) 
    plt.subplots_adjust(wspace=0.1, hspace=0.15, left=0.05, right=0.98, bottom=0.15, top=0.92)

    all_axes = []
    cf_last = None

    # ====================== 5. 循环变量与季节绘图 ======================
    nc_ctl = Dataset(ctl_file)
    nc_exp = Dataset(exp_file)

    plot_idx = 1
    for row_idx, season in enumerate(seasons):
        time_idx = mam_idx if season == 'MAM' else jja_idx
        
        for col_idx, var_name in enumerate(variables):
            print(f"处理变量: {var_name} [{season}]...")
            
            # 提取并计算时间平均
            ctl_data = np.nanmean(nc_ctl.variables[var_name][time_idx, :, :], axis=0)
            exp_data = np.nanmean(nc_exp.variables[var_name][time_idx, :, :], axis=0)
            
            # 计算差值并插值到 0.25x0.25 网格
            diff_raw = exp_data - ctl_data
            diff_ip = regrid_rcm2rgrid(diff_raw, lat2d, lon2d, lat1d, lon1d)
            
            # 施加 CN05.1 的陆地掩码
            diff_masked = np.where(mask_obs, diff_ip, np.nan)

            # --- 绘图 ---
            ax = fig.add_subplot(2, 5, plot_idx, projection=proj)
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
            # 仅第一列显示Y轴标签，仅最后一行显示X轴标签，以保持紧凑
            gl.left_labels = True if col_idx == 0 else False
            gl.bottom_labels = True if row_idx == 1 else False
            gl.xlabel_style = {'size': 9, 'rotation': 0, 'va': 'top', 'ha': 'center'}
            gl.ylabel_style = {'size': 9, 'rotation': 0, 'va': 'center', 'ha': 'right'}

            # 绘制填色图
            cf_last = ax.contourf(lon_grid, lat_grid, diff_masked, 
                                  levels=diff_levels, cmap=diff_cmap, 
                                  extend='both', transform=ccrs.PlateCarree(), zorder=1)

            # 编号标注 (a), (b)...
            char = chr(97 + plot_idx - 1)
            ax.text(0.02, 0.04, f"({char})", transform=ax.transAxes, fontsize=11, 
                    va='bottom', ha='left', zorder=5, 
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1.5))
            
            # 标题标注
            if row_idx == 0: # 只有第一行标出变量名
                ax.set_title(var_titles[col_idx], loc='center', fontsize=12, pad=6)
            
            # 左侧季节标注
            if col_idx == 0:
                ax.text(-0.15, 0.5, season_names[row_idx], transform=ax.transAxes, 
                        fontsize=14, fontweight='normal', va='center', ha='center', rotation=90)

            plot_idx += 1

    nc_ctl.close()
    nc_exp.close()

    # ====================== 6. 共享 Colorbar 与保存 ======================
    # 调整布局，为底部的 Colorbar 留出空间
    pos_left = fig.axes[-5].get_position()
    pos_right = fig.axes[-1].get_position()
    
    # 在所有子图正下方创建一个共享的 colorbar 轴
    cbar_ax = fig.add_axes([pos_left.x0 + 0.1, 0.06, (pos_right.x1 - pos_left.x0) - 0.2, 0.02])
    cb = plt.colorbar(cf_last, cax=cbar_ax, orientation='horizontal', extend='both')
    # cb.set_label('EXP - CTL Difference (W/m²)', fontsize=12, fontweight='normal')
    cb.ax.tick_params(labelsize=10)
    
    # 保存图片
    save_name = "../illustration/off/heat_fluxes_diff.pdf"
    plt.savefig(save_name, dpi=600, bbox_inches='tight')
    plt.close()
    
    print(f"已成功导出至: {save_name}")