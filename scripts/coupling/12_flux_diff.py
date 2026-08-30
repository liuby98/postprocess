import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata
from scipy.stats import ttest_rel  
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

# 变量配置
variables = ['f_fgrnd', 'f_olrg', 'f_rnet', 'f_fsena', 'f_lfevpa']
var_titles = ['Ground Heat Flux', 'Outgoing Longwave Radiation', 
              'Net Radiation', 'Sensible Heat Flux', 'Latent Heat Flux']

diff_levels = np.linspace(-3, 3, 61)
diff_cmap = 'PRGn' 

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

    # 构造包含 17年 x 6个月 = 102个时间步 的索引
    nyears = 17
    nmonths = 6
    ss_idx = [nmonths * i + j for i in range(nyears) for j in range(nmonths)]

    # ====================== 4. 绘图环境准备 ======================
    os.makedirs("../illustration/cpl", exist_ok=True)
    os.environ['CARTOPY_OFFLINE'] = 'true'

    proj = ccrs.LambertConformal(
        central_longitude=PROJ_CENTRAL_LON, 
        central_latitude=PROJ_CENTRAL_LAT, 
        standard_parallels=PROJ_STD_PARALLELS
    )

    fig = plt.figure(figsize=(18, 10)) 
    gs = gridspec.GridSpec(2, 6, figure=fig)
    plt.subplots_adjust(wspace=0.4, hspace=0.25, left=0.05, right=0.95, bottom=0.15, top=0.95)
    
    cf_last = None
    ax_list = []  

    # ====================== 5. 循环变量绘图 ======================
    nc_ctl = Dataset(ctl_file)
    nc_exp = Dataset(exp_file)

    for col_idx, var_name in enumerate(variables):
        print(f"处理变量: {var_name} [Spring & Summer - Multi-step Average]...")
        
        ctl_raw = nc_ctl.variables[var_name][ss_idx, :, :]
        exp_raw = nc_exp.variables[var_name][ss_idx, :, :]
        nlat, nlon = ctl_raw.shape[1], ctl_raw.shape[2]
        
        # 计算年际平均值以进行显著性检验
        ctl_reshaped = ctl_raw.reshape(nyears, nmonths, nlat, nlon)
        exp_reshaped = exp_raw.reshape(nyears, nmonths, nlat, nlon)
        ctl_yearly = np.nanmean(ctl_reshaped, axis=1) 
        exp_yearly = np.nanmean(exp_reshaped, axis=1) 

        # 配对样本 T 检验
        t_stat, p_value = ttest_rel(exp_yearly, ctl_yearly, axis=0, nan_policy='omit')
        
        # 计算多年气候态差异均值
        diff_monthly = exp_raw - ctl_raw 
        diff_reshaped = diff_monthly.reshape(nyears, nmonths, nlat, nlon)
        diff_yearly_diff = np.nanmean(diff_reshaped, axis=1)  
        diff_raw = np.nanmean(diff_yearly_diff, axis=0)       
        
        # 差异场和 p-value 场插值
        diff_ip = regrid_rcm2rgrid(diff_raw, lat2d, lon2d, lat1d, lon1d)
        pval_ip = regrid_rcm2rgrid(p_value, lat2d, lon2d, lat1d, lon1d)

        # 叠加陆地掩码
        diff_masked = np.where(mask_obs, diff_ip, np.nan)
        pval_masked = np.where(mask_obs, pval_ip, np.nan)

        # --- 绘图 ---
        if col_idx == 0:
            ax = fig.add_subplot(gs[0, 0:2], projection=proj)  
        elif col_idx == 1:
            ax = fig.add_subplot(gs[0, 2:4], projection=proj)  
        elif col_idx == 2:
            ax = fig.add_subplot(gs[0, 4:6], projection=proj)  
        elif col_idx == 3:
            ax = fig.add_subplot(gs[1, 1:3], projection=proj)  
        elif col_idx == 4:
            ax = fig.add_subplot(gs[1, 3:5], projection=proj)  
            
        ax_list.append(ax)
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
        
        gl.left_labels = True if col_idx in [0, 3] else False
        gl.bottom_labels = True 
        
        gl.xlabel_style = {'size': 9, 'rotation': 0, 'va': 'top', 'ha': 'center'}
        gl.ylabel_style = {'size': 9, 'rotation': 0, 'va': 'center', 'ha': 'right'}

        # 【修改 1】：改回使用全场原始分布 diff_masked 进行底图填色
        cf_last = ax.contourf(lon_grid, lat_grid, diff_masked, 
                              levels=diff_levels, cmap=diff_cmap, 
                              vmin=-3, vmax=3,
                              extend='both', transform=ccrs.PlateCarree(), zorder=1)

        # 【修改 2】：打点标注显著性区域 (p < 0.05)
        # 获取显著性通过的布尔掩码
        sig_mask = pval_masked < 0.05
        # 设置降采样步长（step数字越大，点越稀疏。当前 0.25 度网格下，取 3 较合适）
        step = 5 
        sig_lons = lon_grid[::step, ::step][sig_mask[::step, ::step]]
        sig_lats = lat_grid[::step, ::step][sig_mask[::step, ::step]]
        
        # 叠加白点散点图（s控制大小，可微调）
        ax.scatter(sig_lons, sig_lats, color='pink', s=8, marker='o', 
                   edgecolors='none', transform=ccrs.PlateCarree(), zorder=2)

        # 编号标注 (a), (b)...
        char = chr(97 + col_idx)
        ax.text(0.03, 0.96, f"({char})", transform=ax.transAxes, fontsize=11, 
                va='top', ha='left', zorder=5, 
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1.5))
        
        # 标题标注
        ax.set_title(var_titles[col_idx], loc='center', fontsize=12, pad=6)

    nc_ctl.close()
    nc_exp.close()

    # ====================== 6. 共享 Colorbar 与保存 ======================
    fig.canvas.draw()  
    
    pos_left = ax_list[0].get_position()  
    pos_right = ax_list[2].get_position() 
    
    cbar_ax = fig.add_axes([pos_left.x0, 0.06, pos_right.x1 - pos_left.x0, 0.02])

    cb_ticks = np.linspace(-3, 3, 11)
    cb = plt.colorbar(cf_last, cax=cbar_ax, orientation='horizontal', extend='both', ticks=cb_ticks)
    cb.ax.tick_params(labelsize=10)
    cb.ax.set_title(r'W m$^{-2}$', fontsize=11, pad=5)
    
    save_name = "../illustration/cpl/FIG12.heat_fluxes_diff.pdf"
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"已成功导出至: {save_name}")