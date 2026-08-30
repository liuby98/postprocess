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
# 2. 逐年循环数据处理 (500 hPa 位势高度)
# ==========================================
def get_500hPa_data_manual(wrfout_file, pb_static, phb_static):
    """
    切片计算每年的春季(MAM)或夏季(JJA)均值，求17年气候态后插值到500hPa位势高度
    """
    fname = os.path.basename(wrfout_file)
    print(f"  [进程启动] 正在读取和计算: {fname} ...")
    
    nyear = 17
    days_per_year = 92  # MAM (31+30+31=92) 或 JJA (30+31+31=92)
    
    nc = Dataset(wrfout_file, 'r')
    var_P = nc.variables['P']
    var_PH = nc.variables['PH']  # 扰动位势
    
    # 获取变量形状以初始化累加数组
    shape_p = var_P.shape[1:]
    shape_ph = var_PH.shape[1:]
    
    p_sum = np.zeros(shape_p, dtype=np.float32)
    ph_sum = np.zeros(shape_ph, dtype=np.float32)
    
    # --- 逐年循环累加 ---
    for yr in range(nyear):
        t_start = yr * days_per_year
        t_end   = t_start + days_per_year
        
        p_sum += np.nanmean(var_P[t_start:t_end, :, :, :], axis=0)
        ph_sum += np.nanmean(var_PH[t_start:t_end, :, :, :], axis=0)
        
    nc.close()
    
    # 计算 17 年的多年气候态平均
    p_clim = p_sum / nyear
    ph_clim = ph_sum / nyear
    
    # 计算全气压 (hPa)
    p_tot_clim = (p_clim + pb_static) / 100.0
    
    # 计算总位势高度 (m) -> 重力加速度取 WRF 标准的 9.81
    z_tot_clim = (ph_clim + phb_static) / 9.81
    
    # --- 垂直空间退交错 (Destagger Z) ---
    # PH 和 PHB 在 bottom_top_stag 层上，退交错至质量点(bottom_top)与 P 层级对齐
    z_unstag = 0.5 * (z_tot_clim[:-1, :, :] + z_tot_clim[1:, :, :])
    
    print(f"  [进程节点] {fname} 正在进行 500 hPa 垂直插值...")
    z_500 = interplevel(z_unstag, p_tot_clim, 500.0)
    
    return z_500

def process_task(args):
    key, fpath, pb_static, phb_static = args
    z_dat = get_500hPa_data_manual(fpath, pb_static, phb_static)
    return key, z_dat

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

    print("读取 wrfinput 静态坐标、基础气压(PB)与基础位势(PHB)...")
    nc_static = Dataset(wrfinput_file, 'r')
    lons = nc_static.variables['XLONG'][0, :, :]
    lats = nc_static.variables['XLAT'][0, :, :]
    pb_static = nc_static.variables['PB'][0, :, :, :]
    phb_static = nc_static.variables['PHB'][0, :, :, :]
    nc_static.close()

    print(" >>> 正在启动多进程 (4个核心) 并发处理 WRF 大文件...")
    tasks = []
    for key, fname in files.items():
        fpath = os.path.join(wrfout_dir, fname)
        tasks.append((key, fpath, pb_static, phb_static))
        
    data = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        results = executor.map(process_task, tasks)
        for key, z_dat in results:
            data[key] = z_dat
            print(f"  ✅ [完成] {key} 数据全部就绪！")

    # --- 计算 Diff 差值 (EXP - CTL) ---
    print("  > 所有文件读取完毕，正在计算 EXP - CTL 差异场...")
    data['DIFF_MAM'] = data['EXP_MAM'] - data['CTL_MAM']
    data['DIFF_JJA'] = data['EXP_JJA'] - data['CTL_JJA']

    print("初始化画布...")
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(16, 18), subplot_kw={'projection': proj})
    axes = axes.flatten() 

    panels = [
        {"key": "CTL_MAM", "label": "(a) CTL MAM", "type": "abs"},
        {"key": "CTL_JJA", "label": "(b) CTL JJA", "type": "abs"},
        {"key": "EXP_MAM", "label": "(c) EXP MAM", "type": "abs"},
        {"key": "EXP_JJA", "label": "(d) EXP JJA", "type": "abs"},
        {"key": "DIFF_MAM", "label": "(e) EXP MAM - CTL MAM", "type": "diff", "ctl_key": "CTL_MAM", "exp_key": "EXP_MAM"},
        {"key": "DIFF_JJA", "label": "(f) EXP JJA - CTL JJA", "type": "diff", "ctl_key": "CTL_JJA", "exp_key": "EXP_JJA"}
    ]

    # 等值线层级设置 (根据东亚 500hPa 高度场的典型范围调整)
    z_levels = np.arange(5400, 5940, 20)   # 绝对场阴影：5400 到 5920 gpm
    diff_z_levels = np.arange(-10, 11, 1)  # 差异场阴影：-30 到 30 gpm 左右

    cf_abs, cf_diff = None, None

    for i, panel in enumerate(panels):
        ax = axes[i]
        is_left_col = (i % 2 == 0)
        is_bottom_row = (i >= 4)
        
        apply_style_and_shp(ax, panel["label"], is_left_col, is_bottom_row)
        z_dat = data[panel["key"]]

        if panel["type"] == "abs":
            # 绝对高度场颜色阴影
            cf = ax.contourf(lons, lats, z_dat, levels=z_levels, cmap='RdYlBu_r', 
                             transform=ccrs.PlateCarree(), extend='both', zorder=1)
            cf_abs = cf
            
            # 黄绿色 5880 gpm 等值线(副高脊线)
            cs = ax.contour(lons, lats, z_dat, levels=[5880], colors='yellowgreen', 
                            linewidths=2.5, transform=ccrs.PlateCarree(), zorder=4)
            ax.clabel(cs, inline=True, fontsize=12, fmt='%1.0f')
            
        else:
            # 高度场差异颜色阴影
            cf = ax.contourf(lons, lats, z_dat, levels=diff_z_levels, cmap='RdBu_r', 
                             transform=ccrs.PlateCarree(), extend='both', zorder=1)
            cf_diff = cf
            
            # (扩展功能) 在差异图中画出 CTL 和 EXP 的 5880 脊线对比，直观展示位移
            ctl_z = data[panel["ctl_key"]]
            exp_z = data[panel["exp_key"]]
            
            # CTL的588使用蓝色虚线
            ax.contour(lons, lats, ctl_z, levels=[5880], colors='blue', linestyles='--', 
                       linewidths=2.0, transform=ccrs.PlateCarree(), zorder=4)
            # EXP的588使用红色实线
            ax.contour(lons, lats, exp_z, levels=[5880], colors='red', linestyles='-', 
                       linewidths=2.0, transform=ccrs.PlateCarree(), zorder=4)

    # Colorbar 排版
    fig.subplots_adjust(bottom=0.08, top=0.95, hspace=0.35, wspace=0.05)
    
    cbar_ax_abs = fig.add_axes([0.25, 0.35, 0.5, 0.012]) 
    cbar_abs = fig.colorbar(cf_abs, cax=cbar_ax_abs, orientation='horizontal')
    cbar_abs.set_label('500 hPa Geopotential Height (gpm)', fontsize=13, fontweight='bold')
    
    cbar_ax_diff = fig.add_axes([0.25, 0.03, 0.5, 0.012])
    cbar_diff = fig.colorbar(cf_diff, cax=cbar_ax_diff, orientation='horizontal')
    cbar_diff.set_label('Geopotential Height Difference (gpm)', fontsize=13, fontweight='bold')

    fig.suptitle("500 hPa Geopotential Height and 5880 gpm Subtropical High (2001-2017)", fontsize=18, y=0.98, fontweight='bold')

    os.makedirs("./figs", exist_ok=True)
    save_name = "./figs/Geopotential_500hPa_Loop.pdf"
    plt.savefig(save_name, dpi=300, bbox_inches='tight') 
    plt.close()
    
    print(f"已保存为: {save_name}")

if __name__ == "__main__":
    main()