import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from wrf import interplevel
import warnings
import concurrent.futures

# 屏蔽常规的无效值警告
warnings.filterwarnings("ignore")

# ==========================================
# 0. 环境变量与静态设置
# ==========================================
os.environ['CARTOPY_OFFLINE'] = 'true'
os.environ['PROJ_NETWORK'] = 'OFF'

COASTLINE_50M = cfeature.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='gray', facecolor='none')
LAKES_50M = cfeature.NaturalEarthFeature('physical', 'lakes', '50m', edgecolor='gray', facecolor='none')

# Shapefile 所在目录 (请确保该路径正确)
shp_dir = "../../shapefile_China/"

# Lambert 投影参数
proj = ccrs.LambertConformal(
    central_longitude=110,
    central_latitude=40,
    standard_parallels=(30, 60)
)

# ==========================================
# 1. 核心底图绘制函数
# ==========================================
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

# ==========================================
# 2. 逐年循环数据处理 (850 hPa 风场)
# ==========================================
def get_850hPa_wind_manual(wrfout_file, pb_static):
    """
    切片计算每年的春季(MAM)或夏季(JJA)均值，并插值到 850 hPa
    """
    fname = os.path.basename(wrfout_file)
    print(f"  [进程启动] 正在读取和计算: {fname} ...")
    
    nyear = 17
    days_per_year = 92  # MAM 或 JJA 均为 92 天
    
    nc = Dataset(wrfout_file, 'r')
    var_U = nc.variables['U']
    var_V = nc.variables['V']
    var_P = nc.variables['P']
    
    shape_u = var_U.shape[1:]
    shape_v = var_V.shape[1:]
    shape_p = var_P.shape[1:]
    
    u_sum = np.zeros(shape_u, dtype=np.float32)
    v_sum = np.zeros(shape_v, dtype=np.float32)
    p_sum = np.zeros(shape_p, dtype=np.float32)
    
    for yr in range(nyear):
        t_start = yr * days_per_year
        t_end   = t_start + days_per_year
        
        u_sum += np.nanmean(var_U[t_start:t_end, :, :, :], axis=0)
        v_sum += np.nanmean(var_V[t_start:t_end, :, :, :], axis=0)
        p_sum += np.nanmean(var_P[t_start:t_end, :, :, :], axis=0)
        
    nc.close()
    
    u_clim = u_sum / nyear
    v_clim = v_sum / nyear
    p_clim = p_sum / nyear
    
    # 空间退交错 (Destagger)
    u_unstag = 0.5 * (u_clim[..., :-1] + u_clim[..., 1:])
    v_unstag = 0.5 * (v_clim[..., :-1, :] + v_clim[..., 1:, :])
    
    wspd_clim = np.sqrt(u_unstag**2 + v_unstag**2)
    p_tot_clim = (p_clim + pb_static) / 100.0

    print(f"  [进程节点] {fname} 正在进行 850 hPa 垂直插值...")
    u_850 = interplevel(u_unstag, p_tot_clim, 850.0)
    v_850 = interplevel(v_unstag, p_tot_clim, 850.0)
    wspd_850 = interplevel(wspd_clim, p_tot_clim, 850.0)
    
    return u_850, v_850, wspd_850

def process_task(args):
    key, fpath, pb_static = args
    u_dat, v_dat, wspd_dat = get_850hPa_wind_manual(fpath, pb_static)
    return key, u_dat, v_dat, wspd_dat

# ==========================================
# 3. 主程序
# ==========================================
def main():
    wrfout_dir = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
    wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/ICBC_2002/wrfinput_d01"
    
    files = {
        'CTL_MAM': "wrfout_2001-2023_MAM_daymean_nogravel.nc",
        'CTL_JJA': "wrfout_2001-2023_JJA_daymean_nogravel.nc",
        'EXP_MAM': "wrfout_2001-2023_MAM_daymean_gravel.nc",
        'EXP_JJA': "wrfout_2001-2023_JJA_daymean_gravel.nc"
    }

    print("读取 wrfinput 静态坐标与基础气压...")
    nc_static = Dataset(wrfinput_file, 'r')
    lons = nc_static.variables['XLONG'][0, :, :]
    lats = nc_static.variables['XLAT'][0, :, :]
    pb_static = nc_static.variables['PB'][0, :, :, :]
    nc_static.close()

    print(" >>> 正在启动多进程并发处理 WRF 大文件...")
    tasks = []
    for key, fname in files.items():
        fpath = os.path.join(wrfout_dir, fname)
        tasks.append((key, fpath, pb_static))
        
    data = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        results = executor.map(process_task, tasks)
        for key, u_dat, v_dat, wspd_dat in results:
            data[key] = (u_dat, v_dat, wspd_dat)
            print(f"  ✅ [完成] {key} 数据全部就绪！")

    # --- 计算 Diff 差值 (EXP - CTL) ---
    print("  > 所有文件读取完毕，正在计算 EXP - CTL 差异场...")
    data['DIFF_MAM'] = (
        data['EXP_MAM'][0] - data['CTL_MAM'][0],  # U diff
        data['EXP_MAM'][1] - data['CTL_MAM'][1],  # V diff
        data['EXP_MAM'][2] - data['CTL_MAM'][2]   # Wspd diff
    )
    data['DIFF_JJA'] = (
        data['EXP_JJA'][0] - data['CTL_JJA'][0],
        data['EXP_JJA'][1] - data['CTL_JJA'][1],
        data['EXP_JJA'][2] - data['CTL_JJA'][2]
    )

    print("初始化画布...")
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(16, 18), subplot_kw={'projection': proj})
    axes = axes.flatten() 

    panels = [
        {"key": "CTL_MAM", "label": "(a) CTL MAM", "type": "abs"},
        {"key": "CTL_JJA", "label": "(b) CTL JJA", "type": "abs"},
        {"key": "EXP_MAM", "label": "(c) EXP MAM", "type": "abs"},
        {"key": "EXP_JJA", "label": "(d) EXP JJA", "type": "abs"},
        {"key": "DIFF_MAM", "label": "(e) EXP MAM - CTL MAM", "type": "diff"},
        {"key": "DIFF_JJA", "label": "(f) EXP JJA - CTL JJA", "type": "diff"}
    ]

    # 850 hPa 风速层级 (通常在 2~16 m/s 之间)
    wspd_levels = np.arange(2, 17, 2) 
    # 850 hPa 差异风速层级
    diff_wspd_levels = np.arange(-2.0, 2.2, 0.2) 

    # 矢量抽稀比例 (X和Y方向每隔12个格点画一个箭头，可根据网格分辨率微调)
    skip = (slice(None, None, 12), slice(None, None, 12))

    cf_abs, cf_diff = None, None

    for i, panel in enumerate(panels):
        ax = axes[i]
        
        is_left_col = (i % 2 == 0)
        is_bottom_row = (i >= 4)
        
        apply_style_and_shp(ax, panel["label"], is_left_col, is_bottom_row)
        u_dat, v_dat, wspd_dat = data[panel["key"]]

        if panel["type"] == "abs":
            # 风速阴影
            cf = ax.contourf(lons, lats, wspd_dat, levels=wspd_levels, cmap='YlOrRd', 
                             transform=ccrs.PlateCarree(), extend='both', zorder=1)
            cf_abs = cf
            
            # 风向矢量箭头 (使用 quvier)
            q = ax.quiver(lons[skip], lats[skip], u_dat[skip], v_dat[skip], 
                          transform=ccrs.PlateCarree(), color='black', 
                          scale=150, width=0.003, headwidth=4, zorder=4)
            # 在子图左下角添加图例参考箭头 (10 m/s)
            ax.quiverkey(q, X=0.08, Y=0.05, U=10, label='10 m/s', labelpos='E', coordinates='axes')

        else:
            # 差异场风速阴影
            cf = ax.contourf(lons, lats, wspd_dat, levels=diff_wspd_levels, cmap='RdBu_r', 
                             transform=ccrs.PlateCarree(), extend='both', zorder=1)
            cf_diff = cf
            
            # 差异场风矢量箭头
            q = ax.quiver(lons[skip], lats[skip], u_dat[skip], v_dat[skip], 
                          transform=ccrs.PlateCarree(), color='black', 
                          scale=50, width=0.003, headwidth=4, zorder=4)
            # 在子图左下角添加差异图例参考箭头 (2 m/s)
            ax.quiverkey(q, X=0.08, Y=0.05, U=2, label='2 m/s', labelpos='E', coordinates='axes')

    # Colorbar 排版
    fig.subplots_adjust(bottom=0.08, top=0.95, hspace=0.35, wspace=0.05)
    
    cbar_ax_abs = fig.add_axes([0.25, 0.35, 0.5, 0.012]) 
    cbar_abs = fig.colorbar(cf_abs, cax=cbar_ax_abs, orientation='horizontal')
    cbar_abs.set_label('850 hPa Wind Speed (m s$^{-1}$)', fontsize=13, fontweight='bold')
    
    cbar_ax_diff = fig.add_axes([0.25, 0.03, 0.5, 0.012])
    cbar_diff = fig.colorbar(cf_diff, cax=cbar_ax_diff, orientation='horizontal')
    cbar_diff.set_label('Wind Speed Difference (m s$^{-1}$)', fontsize=13, fontweight='bold')

    fig.suptitle("850 hPa Mean Wind Field (2001-2017)", fontsize=18, y=0.98, fontweight='bold')

    os.makedirs("./figs", exist_ok=True)
    save_name = "./figs/Wind_850hPa_Loop.png" 
    plt.savefig(save_name, dpi=300, bbox_inches='tight') 
    plt.close()
    
    print(f"已保存为: {save_name}")

if __name__ == "__main__":
    main()