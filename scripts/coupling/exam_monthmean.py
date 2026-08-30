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
PROJ_CENTRAL_LON = 110.0        
PROJ_CENTRAL_LAT = 40.0         
PROJ_STD_PARALLELS = (30.0, 60.0) 

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ========== 控制开关 ==========
# True  → 绘制 T_grnd - T_xy 的 CTL/EXP 及差异（默认）
# False → 绘制感热通量 f_fsena 的 CTL/EXP 及差异
PLOT_TDIFF =  False  # 根据需要改为 True

# 路径配置
path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
ctl_file = path + "colmrun_2001-2017_nogravel.nc"
exp_file = path + "colmrun_2001-2017_gravel.nc"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
obs_mask_file = "/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc"
shp_dir = "../../shapefile_China/"

# ====================== 2. 变量配置（根据开关动态选择）======================
if PLOT_TDIFF:
    # 配置：温度差 (T_grnd - T_xy)
    var_cfg = {
        'mean_levels': np.linspace(-10, 10, 41),      # CTL/EXP 用 ±10℃，间隔0.5
        'diff_levels': np.linspace(-1, 1, 21),        # 差异 ±1℃，间隔0.1
        'mean_cmap': 'RdBu_r',
        'diff_cmap': 'RdBu_r',
        'mean_label': '$T_{grnd} - T_{xy}$  (°C)',
        'diff_label': 'EXP - CTL (°C)',
        'title_var': '$T_{grnd}-T_{xy}$',
        'fname': 'tdiff.png',
        'var1': 'f_t_grnd',
        'var2': 'f_xy_t',
        'is_diff_var': True,
        # 手动指定色条刻度，使显示更清晰
        'mean_ticks': np.linspace(-10, 10, 9),        # -10,-7.5,...,10
        'diff_ticks': np.linspace(-1, 1, 5)           # -1,-0.5,0,0.5,1
    }
else:
    # 配置：感热通量 (f_fsena)
    var_cfg = {
        'mean_levels': np.linspace(-100, 100, 41),    # CTL/EXP 用 ±100，间隔5
        'diff_levels': np.linspace(-15, 15, 31),      # 差异 ±15，间隔1
        'mean_cmap': 'RdBu_r',
        'diff_cmap': 'RdBu_r',
        'mean_label': 'Sensible Heat Flux (W/m²)',
        'diff_label': 'EXP - CTL (W/m²)',
        'title_var': 'Sensible Heat Flux',
        'fname': 'fsena.png',
        'var1': 'f_fsena',
        'var2': None,
        'is_diff_var': False,
        'mean_ticks': np.linspace(-100, 100, 9),      # -100,-75,...,100
        'diff_ticks': np.linspace(-15, 15, 7)         # -15,-10,-5,0,5,10,15
    }

# ====================== 3. 核心函数 ======================
def regrid_rcm2rgrid(var2d, lat2d, lon2d, lat1d, lon1d):
    pts = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    vals = var2d.ravel()
    gx, gy = np.meshgrid(lon1d, lat1d)
    interp = griddata(pts, vals, (gx, gy), method='linear')
    nan_mask = np.isnan(interp)
    if nan_mask.any():
        interp[nan_mask] = griddata(pts, vals, (gx, gy), method='nearest')[nan_mask]
    return interp

def get_clean_var(nc_obj, var_name, indices):
    raw = nc_obj.variables[var_name][indices, :, :]
    if np.ma.isMaskedArray(raw): 
        raw = raw.filled(np.nan)
    raw[raw <= -0.99e36] = np.nan
    return raw

def process_season(nc_ctl, nc_exp, indices, cfg, mask_obs, lat2d, lon2d, lat1d, lon1d):
    """
    根据配置 cfg 处理单个季节的数据，返回 CTL, EXP, Diff 三个插值掩码场
    """
    if cfg['is_diff_var']:
        c_v1 = get_clean_var(nc_ctl, cfg['var1'], indices)
        c_v2 = get_clean_var(nc_ctl, cfg['var2'], indices)
        e_v1 = get_clean_var(nc_exp, cfg['var1'], indices)
        e_v2 = get_clean_var(nc_exp, cfg['var2'], indices)
        c_data = c_v1 - c_v2
        e_data = e_v1 - e_v2
    else:
        c_data = get_clean_var(nc_ctl, cfg['var1'], indices)
        e_data = get_clean_var(nc_exp, cfg['var1'], indices)

    c_mean = np.nanmean(c_data, axis=0)
    e_mean = np.nanmean(e_data, axis=0)
    d_mean = e_mean - c_mean

    c_ip = regrid_rcm2rgrid(c_mean, lat2d, lon2d, lat1d, lon1d)
    e_ip = regrid_rcm2rgrid(e_mean, lat2d, lon2d, lat1d, lon1d)
    d_ip = regrid_rcm2rgrid(d_mean, lat2d, lon2d, lat1d, lon1d)

    c_masked = np.where(mask_obs, c_ip, np.nan)
    e_masked = np.where(mask_obs, e_ip, np.nan)
    d_masked = np.where(mask_obs, d_ip, np.nan)
    return c_masked, e_masked, d_masked

if __name__ == '__main__':
    # ====================== 4. 读取网格与掩码数据 ======================
    with Dataset(wrfinput_file) as f_wrf:
        lat2d = f_wrf.variables['XLAT'][0, :, :]
        lon2d = f_wrf.variables['XLONG'][0, :, :]

    with Dataset(obs_mask_file) as f_obs:
        lat1d = f_obs.variables['lat'][:]
        lon1d = f_obs.variables['lon'][:]
        tm_obs = f_obs.variables['tm'][0, :, :]
        mask_obs = ~np.ma.getmaskarray(tm_obs) & ~np.isnan(tm_obs)

    lon_grid, lat_grid = np.meshgrid(lon1d, lat1d)

    # 每年 6 个月，春季（0,1,2）和夏季（3,4,5）
    nyears = 17
    nmonths = 6
    spring_idx = [nmonths * i + j for i in range(nyears) for j in range(3)]
    summer_idx = [nmonths * i + j for i in range(nyears) for j in range(3, 6)]

    # ====================== 5. 数据处理 ======================
    print(f"当前绘制模式: {'温度差 (T_grnd - T_xy)' if PLOT_TDIFF else '感热通量 (f_fsena)'}")
    nc_ctl = Dataset(ctl_file)
    nc_exp = Dataset(exp_file)

    c_sp, e_sp, d_sp = process_season(nc_ctl, nc_exp, spring_idx, var_cfg, mask_obs, lat2d, lon2d, lat1d, lon1d)
    c_su, e_su, d_su = process_season(nc_ctl, nc_exp, summer_idx, var_cfg, mask_obs, lat2d, lon2d, lat1d, lon1d)

    nc_ctl.close()
    nc_exp.close()

    # ====================== 6. 绘图 ======================
    out_dir = "./figs"
    os.makedirs(out_dir, exist_ok=True)
    os.environ['CARTOPY_OFFLINE'] = 'true'

    proj = ccrs.LambertConformal(
        central_longitude=PROJ_CENTRAL_LON, 
        central_latitude=PROJ_CENTRAL_LAT, 
        standard_parallels=PROJ_STD_PARALLELS
    )

    fig, axes = plt.subplots(2, 3, figsize=(18, 12), subplot_kw={'projection': proj})
    plt.subplots_adjust(wspace=0.1, hspace=0.15)

    plot_data = [
        [c_sp, e_sp, d_sp],
        [c_su, e_su, d_su]
    ]
    row_names = ['Spring', 'Summer']
    col_names = [
        f'CTL ({var_cfg["title_var"]})',
        f'EXP ({var_cfg["title_var"]})',
        'EXP - CTL'
    ]

    cf_mean = None
    cf_diff = None

    for i in range(2):
        for j in range(3):
            ax = axes[i, j]
            ax.set_extent(PLOT_EXTENT, crs=ccrs.PlateCarree())

            # 边界
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

            # 网格线
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

            ax.set_title(f'{row_names[i]} {col_names[j]}', fontsize=14, pad=10)

            data = plot_data[i][j]
            if j < 2:
                cf = ax.contourf(lon_grid, lat_grid, data,
                                 levels=var_cfg['mean_levels'], cmap=var_cfg['mean_cmap'],
                                 extend='both', transform=ccrs.PlateCarree(), zorder=1)
                cf_mean = cf
            else:
                cf = ax.contourf(lon_grid, lat_grid, data,
                                 levels=var_cfg['diff_levels'], cmap=var_cfg['diff_cmap'],
                                 extend='both', transform=ccrs.PlateCarree(), zorder=1)
                cf_diff = cf

    # ====================== 7. 动态对齐 Colorbar ======================
    ax00 = axes[0, 0]
    ax01 = axes[0, 1]
    ax02 = axes[0, 2]

    pos00 = ax00.get_position()
    pos01 = ax01.get_position()
    pos02 = ax02.get_position()

    left_mean = pos00.x0
    right_mean = pos01.x1
    width_mean = right_mean - left_mean
    left_diff = pos02.x0
    right_diff = pos02.x1
    width_diff = right_diff - left_diff

    bottom = 0.05
    height = 0.02

    # 第一个色条（CTL/EXP 共享）
    cbar_ax_mean = fig.add_axes([left_mean, bottom, width_mean, height])
    cb_mean = plt.colorbar(cf_mean, cax=cbar_ax_mean, orientation='horizontal', extend='both')
    cb_mean.set_ticks(var_cfg['mean_ticks'])
    cb_mean.ax.tick_params(labelsize=11)
    cb_mean.set_label(var_cfg['mean_label'], fontsize=12)

    # 第二个色条（差异）
    cbar_ax_diff = fig.add_axes([left_diff, bottom, width_diff, height])
    cb_diff = plt.colorbar(cf_diff, cax=cbar_ax_diff, orientation='horizontal', extend='both')
    cb_diff.set_ticks(var_cfg['diff_ticks'])
    cb_diff.ax.tick_params(labelsize=11)
    cb_diff.set_label(var_cfg['diff_label'], fontsize=12)

    save_name = os.path.join(out_dir, var_cfg['fname'])
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"组图已保存至: {save_name}")