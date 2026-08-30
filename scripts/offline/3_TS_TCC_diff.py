import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata
import warnings
warnings.filterwarnings("ignore")
from scipy.stats import pearsonr, t as t_dist
from scipy.signal import detrend
from scipy.ndimage import gaussian_filter1d
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import csv

# ====================== 字体全局设置 ======================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ====================== 1. 投影工具与路径配置 ======================
try:
    from pyproj import Proj
except ImportError:
    pass
p = Proj(proj='lcc', lat_1=30, lat_2=60, lat_0=40, lon_0=110, datum='WGS84')

path          = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
sim1_file     = path + "colmoff_2001-2023_monmean_nogravel.nc"
sim2_file     = path + "colmoff_2001-2023_monmean_gravel.nc"
cn05_file     = "/share/home/dq117/CN05.1/CN05.1_Tm_2020_daily_025x025.nc"
shp_dir       = "../../shapefile_China/"
nyear         = 17

seasons_config = {
    'MAM': {
        'obs_file': "/share/home/dq135/reference/HOM_TS0_40cm_025x025_1991-2017_MAM.nc",
        'model_idx_start': 0, 'model_idx_end': 3, 'name': 'Spring'
    },
    'JJA': {
        'obs_file': "/share/home/dq135/reference/HOM_TS0_40cm_025x025_1991-2017_JJA.nc",
        'model_idx_start': 3, 'model_idx_end': 6, 'name': 'Summer'
    }
}

target_layers = [
    {'name': '0cm', 'obs_idx': 0, 'colm_depth': 0.00},
    {'name': '15cm', 'obs_idx': 3, 'colm_depth': 0.15},
    {'name': '40cm', 'obs_idx': 5, 'colm_depth': 0.40}
]
obs_lev = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.40])
colm_lev = np.array([0.0071, 0.0279, 0.0623, 0.1189, 0.2122, 0.3661, 0.6198, 1.0380, 1.7276, 2.8646])

table_metrics_data = []

# ====================== 2. 计算 17 年样本的单尾 p=0.05 临界值 ======================
# 对于 17 年序列，自由度 df = 15，计算单尾 0.05 的阈值
df = nyear - 2
t_crit = t_dist.ppf(0.95, df) # 单尾 0.05 对应的 95% 分位数
r_crit = np.sqrt(t_crit**2 / (df + t_crit**2)) # 计算结果约为 0.412
print(f"[*] 设定的单尾显著性检验 (p=0.05) 最小 TCC 临界值为: r = {r_crit:.3f}")

# ====================== 3. 读取公共网格与核心函数 ======================
f_wrf = Dataset(wrfinput_file)
lat2d = f_wrf.variables['XLAT'][0, :, :]
lon2d = f_wrf.variables['XLONG'][0, :, :]
f_wrf.close()

ff = Dataset(seasons_config['MAM']['obs_file'])
lat1d = ff.variables['latitude'][:]
lon1d = ff.variables['longitude'][:]
ff.close()

cn05 = Dataset(cn05_file)
obs_cn05 = np.nanmean(cn05.variables['tm'][:], axis=0)
lat_cn = cn05.variables['lat'][:]
lon_cn = cn05.variables['lon'][:]
cn05.close()

mask_cn05 = ~np.isnan(obs_cn05)
grid_x, grid_y = np.meshgrid(lon1d, lat1d)
points = np.column_stack((np.tile(lon_cn, len(lat_cn)), np.repeat(lat_cn, len(lon_cn))))
values = mask_cn05.ravel().astype(float)
mask_target = griddata(points, values, (grid_x, grid_y), method='nearest') > 0.5

rad = np.pi / 180
cos_lat = np.cos(lat1d * rad)
dlon = lon1d[1] - lon1d[0]
wgt_2d = np.outer(cos_lat, np.ones_like(lon1d)) * dlon
wgt_2d = np.where(mask_target, wgt_2d, np.nan)

def regrid_rcm2rgrid(var2d, lat2d, lon2d, lat1d, lon1d):
    pts = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    vals = var2d.ravel()
    gx, gy = np.meshgrid(lon1d, lat1d)
    interp = griddata(pts, vals, (gx, gy), method='linear')
    nan_mask = np.isnan(interp)
    if nan_mask.any():
        interp[nan_mask] = griddata(pts, vals, (gx, gy), method='nearest')[nan_mask]
    return interp

def spatial_weighted_avg(data, weight):
    mask = ~np.isnan(data) & ~np.isnan(weight)
    if np.sum(mask) == 0: return np.nan
    return np.average(data[mask], weights=weight[mask])

def weighted_spatial_correlation(obs, sim, wgt):
    o, s, w = obs.flatten(), sim.flatten(), wgt.flatten()
    valid = ~np.isnan(o) & ~np.isnan(s) & ~np.isnan(w)
    if np.sum(valid) < 3: return np.nan, np.nan
    o_v, s_v, w_v = o[valid], s[valid], w[valid]
    
    mean_o = np.average(o_v, weights=w_v)
    mean_s = np.average(s_v, weights=w_v)
    cov = np.average((o_v - mean_o) * (s_v - mean_s), weights=w_v)
    var_o = np.average((o_v - mean_o)**2, weights=w_v)
    var_s = np.average((s_v - mean_s)**2, weights=w_v)
    
    scc = cov / np.sqrt(var_o * var_s)
    return scc

def apply_detrend_3d(data_3d, mask_2d):
    out = np.copy(data_3d)
    vy, vx = np.where(mask_2d)
    for y, x in zip(vy, vx):
        ts = data_3d[:, y, x]
        valid = ~np.isnan(ts)
        if valid.sum() > 2:
            out[valid, y, x] = detrend(ts[valid]) + np.mean(ts[valid])
    return out

# ====================== 4. 大图综合排版与绘制 ======================
os.makedirs("../illustration/off", exist_ok=True)
os.environ['CARTOPY_OFFLINE'] = 'true'
os.environ['PROJ_NETWORK'] = 'OFF'

proj = ccrs.LambertConformal(central_longitude=110, central_latitude=40, standard_parallels=(30, 60))
lon_grid, lat_grid = np.meshgrid(lon1d, lat1d)

fig = plt.figure(figsize=(20, 10))
gs = gridspec.GridSpec(2, 4, width_ratios=[1.2, 1, 1, 1])

# 间距与外沿设置
plt.subplots_adjust(wspace=0.14, hspace=0.06, left=0.04, right=0.96, bottom=0.15, top=0.95)

colors = ['forestgreen', 'darkorange', 'rebeccapurple']
cf_diff_last = None
ax_dist_list = []
ax_map_list = [[], []]

def apply_style_and_shp(ax):
    ax.set_extent([80, 130, 10, 55], crs=ccrs.PlateCarree())
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
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False,
        linewidth=0.6, color='gray', alpha=0.5, linestyle='--', zorder=2
    )
    gl.xlocator = mticker.FixedLocator(np.arange(70, 140, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
    gl.top_labels = False; gl.right_labels = False
    gl.bottom_labels = True; gl.left_labels = True
    gl.xlabel_style = {'size': 9, 'rotation': 0, 'va': 'top', 'ha': 'center'}
    gl.ylabel_style = {'size': 9, 'rotation': 0, 'va': 'center', 'ha': 'right'}
    try:
        gl.padding = 6
    except Exception:
        pass

seasons_list = ['MAM', 'JJA']
for row_idx, season in enumerate(seasons_list):
    cfg = seasons_config[season]
    season_name = cfg['name']
    print(f"\n================ 开始处理季节: {season} ================")
    
    ax_dist = fig.add_subplot(gs[row_idx, 0])
    ax_dist_list.append(ax_dist)
    
    for col_idx, layer in enumerate(target_layers):
        layer_name = layer['name']
        obs_idx = layer['obs_idx']
        colm_depth = layer['colm_depth']
        print(f"  -> [{layer_name}] 开始计算...")

        # 读取观测
        ff = Dataset(cfg['obs_file'])
        var_obs = ff.variables['TS']
        data_obs = var_obs[30:81, obs_idx, :, :] 
        if hasattr(var_obs, '_FillValue'): data_obs[data_obs == var_obs._FillValue] = np.nan
        data_obs[np.abs(data_obs) > 1e30] = np.nan
        ff.close()

        temp_17 = np.zeros((nyear, lat1d.size, lon1d.size))
        for yr in range(nyear): temp_17[yr] = np.nanmean(data_obs[3*yr:3*yr+3], axis=0)
        obs_mean = np.nanmean(temp_17, axis=0) 
        mask_obs = ~np.isnan(obs_mean)

        sim1_17, sim2_17 = np.zeros((nyear, lat1d.size, lon1d.size)), np.zeros((nyear, lat1d.size, lon1d.size))
        idx_s, idx_e = cfg['model_idx_start'], cfg['model_idx_end']
        def interp_to_depth(temps): return np.interp(colm_depth, colm_lev, temps)

        for yr in range(nyear):
            ds1 = Dataset(sim1_file)
            d1 = ds1.variables['f_t_soisno'][6*yr+idx_s:6*yr+idx_e, :, :, 5:15]
            if hasattr(ds1.variables['f_t_soisno'], '_FillValue'): d1[d1 == ds1.variables['f_t_soisno']._FillValue] = np.nan
            d1[np.abs(d1) > 1e30] = np.nan
            sim1_17[yr] = regrid_rcm2rgrid(np.apply_along_axis(interp_to_depth, axis=2, arr=np.nanmean(d1, axis=0)-273.15), lat2d, lon2d, lat1d, lon1d)
            ds1.close()

            ds2 = Dataset(sim2_file)
            d2 = ds2.variables['f_t_soisno'][6*yr+idx_s:6*yr+idx_e, :, :, 5:15]
            if hasattr(ds2.variables['f_t_soisno'], '_FillValue'): d2[d2 == ds2.variables['f_t_soisno']._FillValue] = np.nan
            d2[np.abs(d2) > 1e30] = np.nan
            sim2_17[yr] = regrid_rcm2rgrid(np.apply_along_axis(interp_to_depth, axis=2, arr=np.nanmean(d2, axis=0)-273.15), lat2d, lon2d, lat1d, lon1d)
            ds2.close()

        for yr in range(nyear):
            sim1_17[yr] = np.where(~np.isnan(temp_17[yr]), sim1_17[yr], np.nan)
            sim2_17[yr] = np.where(~np.isnan(temp_17[yr]), sim2_17[yr], np.nan)

        # 线性去势
        temp_17 = apply_detrend_3d(temp_17, mask_obs)
        sim1_17 = apply_detrend_3d(sim1_17, mask_obs)
        sim2_17 = apply_detrend_3d(sim2_17, mask_obs)
        obs_mean = np.nanmean(temp_17, axis=0)

        # 【核心修正】：精简了统计循环，去掉了不必要的差值检验
        TCC1, TCC2 = np.full_like(obs_mean, np.nan), np.full_like(obs_mean, np.nan)
        for i in range(lat1d.size):
            for j in range(lon1d.size):
                if mask_obs[i,j]:
                    v = ~np.isnan(sim1_17[:, i, j]) & ~np.isnan(temp_17[:, i, j]) & ~np.isnan(sim2_17[:, i, j])
                    n = v.sum()
                    if n >= 3: 
                        o = temp_17[:, i, j][v]
                        s1 = sim1_17[:, i, j][v]
                        s2 = sim2_17[:, i, j][v]
                        TCC1[i,j], _ = pearsonr(s1, o)
                        TCC2[i,j], _ = pearsonr(s2, o)

        TCC1 = np.where(mask_target, TCC1, np.nan)
        TCC2 = np.where(mask_target, TCC2, np.nan)
        
        # 【核心修正】：计算单尾显著性并集掩码 (任意一个模式通过单尾正相关检验即填色)
        sig_ctl = mask_target & ~np.isnan(TCC1) & (TCC1 > r_crit)
        sig_exp = mask_target & ~np.isnan(TCC2) & (TCC2 > r_crit)
        sig_union = sig_ctl | sig_exp  # 取并集
        
        # 只在显著并集的区域上计算差值并填色
        diff = np.where(sig_union, TCC2 - TCC1, np.nan)

        # 扣除气候态求距平并计算 ACC
        obs_anom = temp_17 - obs_mean
        sim1_anom = sim1_17 - np.nanmean(sim1_17, axis=0)
        sim2_anom = sim2_17 - np.nanmean(sim2_17, axis=0)

        acc1_list, acc2_list = [], []
        for yr in range(nyear):
            acc1_list.append(weighted_spatial_correlation(obs_anom[yr], sim1_anom[yr], wgt_2d))
            acc2_list.append(weighted_spatial_correlation(obs_anom[yr], sim2_anom[yr], wgt_2d))

        tcc1_domain_avg = spatial_weighted_avg(TCC1, wgt_2d)
        tcc2_domain_avg = spatial_weighted_avg(TCC2, wgt_2d)
        acc1_time_avg = np.nanmean(acc1_list)
        acc2_time_avg = np.nanmean(acc2_list)
        rmse1_avg = spatial_weighted_avg(np.sqrt(np.nanmean((sim1_17 - temp_17)**2, axis=0)), wgt_2d)
        rmse2_avg = spatial_weighted_avg(np.sqrt(np.nanmean((sim2_17 - temp_17)**2, axis=0)), wgt_2d)

        table_metrics_data.append([season, layer_name, "CTL OFF", f"{tcc1_domain_avg:.3f}", f"{acc1_time_avg:.3f}", f"{rmse1_avg:.3f}"])
        table_metrics_data.append([season, layer_name, "EXP OFF", f"{tcc2_domain_avg:.3f}", f"{acc2_time_avg:.3f}", f"{rmse2_avg:.3f}"])

        # ------------------ (1) 绘制平滑概率分布图 ------------------
        tcc1_raw = TCC1[mask_target & ~np.isnan(TCC1)]
        tcc2_raw = TCC2[mask_target & ~np.isnan(TCC2)]
        
        valid1_sig = mask_target & ~np.isnan(TCC1) & (TCC1 > r_crit)
        valid2_sig = mask_target & ~np.isnan(TCC2) & (TCC2 > r_crit)
        
        tcc1_sig_mean = np.nanmean(TCC1[valid1_sig]) if np.sum(valid1_sig) > 0 else np.nan
        tcc2_sig_mean = np.nanmean(TCC2[valid2_sig]) if np.sum(valid2_sig) > 0 else np.nan
        
        c = colors[col_idx]
        if len(tcc1_raw) > 1 and len(tcc2_raw) > 1:
            bins_edges = np.linspace(-1, 1.1, 50) 
            hist1_counts, _ = np.histogram(tcc1_raw, bins=bins_edges, density=False)
            hist2_counts, _ = np.histogram(tcc2_raw, bins=bins_edges, density=False)
            
            prob1 = hist1_counts / len(tcc1_raw)
            prob2 = hist2_counts / len(tcc2_raw)
            
            # 高斯平滑
            prob1_smooth = gaussian_filter1d(prob1, sigma=1.2)
            prob2_smooth = gaussian_filter1d(prob2, sigma=1.2)
            bin_centers = 0.5 * (bins_edges[1:] + bins_edges[:-1])
            
            # 【核心修改点】：在图例中显示括号里的均值 μ
            label_ctl = f'{layer_name} CTL ($\mu$={tcc1_sig_mean:.2f})'
            label_exp = f'{layer_name} EXP ($\mu$={tcc2_sig_mean:.2f})'
            
            ax_dist.plot(bin_centers, prob1_smooth, color=c, ls='-', lw=1.8, alpha=0.9, label=label_ctl)
            ax_dist.plot(bin_centers, prob2_smooth, color=c, ls='--', lw=1.8, alpha=0.9, label=label_exp)

        # ------------------ (2) 绘制差异图 (取消显著性打点，只保留显著并集填色) ------------------
        ax_map = fig.add_subplot(gs[row_idx, col_idx + 1], projection=proj)
        ax_map_list[row_idx].append(ax_map)
        apply_style_and_shp(ax_map)
        
        cf_diff_last = ax_map.contourf(lon_grid, lat_grid, diff, levels=np.arange(-0.2, 0.21, 0.05), cmap='RdBu_r', extend='both', transform=ccrs.PlateCarree(), zorder=1)
        
        char_map = chr(97 + row_idx*4 + col_idx + 1)
        ax_map.text(0.03, 0.96, f"({char_map}) {layer_name}", transform=ax_map.transAxes, fontsize=11, va='top', ha='left', zorder=5)
        if row_idx == 0 and col_idx == 0:
            ax_map.set_title("TCC Difference", loc='left', fontsize=11, pad=5)

        # ------------------ (3) 嵌入比例条形图 ------------------
        # 因为 diff 只在 sig_union 的网格上有数据，所以面积统计天然就是只针对“所有填色的显著并集网格”
        valid_diff_mask = ~np.isnan(diff) & mask_target 
        total_wgt = np.nansum(wgt_2d[valid_diff_mask])
        area = [0, 0]
        if total_wgt > 0:
            area[0] = np.nansum(wgt_2d[(diff < 0) & valid_diff_mask]) / total_wgt * 100.0
            area[1] = np.nansum(wgt_2d[(diff > 0) & valid_diff_mask]) / total_wgt * 100.0

        ax_inset = inset_axes(ax_map, width="70%", height="4.5%", loc='lower left', bbox_to_anchor=(0.05, 0.08, 0.9, 0.9), bbox_transform=ax_map.transAxes, borderpad=0)
        ax_inset.patch.set_alpha(0.0)
        left = 0
        for val, col in zip(area, ['lightblue', 'lightcoral']):
            ax_inset.barh(0, val, left=left, color=col, height=0.7, edgecolor='black', linewidth=0.6)
            if val > 5: ax_inset.text(left + val/2, 0, f"{val:.1f}%", va='center', ha='center', color='black', fontsize=7)
            left += val
        for spine in ax_inset.spines.values(): spine.set_visible(False)
        ax_inset.set(xlim=(0, 100), ylim=(-0.5, 0.5), xticks=[], yticks=[])
        ax_inset.set_title('Areal Fraction (%)', fontsize=7, pad=2)
        ax_inset.legend(handles=[Rectangle((0,0),1,1, facecolor=c, edgecolor='black', lw=0.6) for c in ['lightblue', 'lightcoral']],
                        labels=['Negative', 'Positive'], loc='upper center', bbox_to_anchor=(0.5, -0.03), ncol=2, fontsize=7, frameon=False)


    # ------------------ (4) 概率图外挂单尾阈值线与统一格式 ------------------
    # 绘制灰色虚线，代表单尾 p=0.05 的及格线
    ax_dist.axvline(x=r_crit, color='red', ls=':', lw=2.0, zorder=1)
    trans = ax_dist.get_xaxis_transform()
    # ax_dist.text(r_crit + 0.02, 0.95, f'$p=0.05$\n($r={r_crit:.3f}$)', transform=trans, color='dimgray', fontsize=9, va='top', ha='left')

    ax_dist.set_xlim(0.0, 1.02)
    ax_dist.set_xlabel('')
    ax_dist.set_ylabel('Probability', fontsize=11, labelpad=8)
    
    char_dist = chr(97 + row_idx*4)
    if row_idx == 0:
        ax_dist.set_title(f"TCC", loc='left', fontsize=11, pad=5)
        
    ax_dist.text(0.03, 0.96, f"({char_dist}) {season_name}", transform=ax_dist.transAxes, 
                fontsize=11, va='top', ha='left', zorder=5)
    
    handles, labels = ax_dist.get_legend_handles_labels()
    order = [0, 1, 2, 3, 4, 5]
    if len(handles) == 6:
        ax_dist.legend([handles[i] for i in order], [labels[i] for i in order], 
                      loc='upper left', ncol=1, frameon=False, fontsize=8, bbox_to_anchor=(0.02, 0.92))
    else:
        ax_dist.legend(loc='upper left', ncol=1, frameon=False, fontsize=8, bbox_to_anchor=(0.02, 0.92))

    ax_dist.grid(True, ls='--', alpha=0.4)


# ====================== 5. 绝对坐标系动态对齐 ======================
plt.draw() 

pos_a = ax_dist_list[0].get_position()
pos_b = ax_map_list[0][0].get_position()
pos_e = ax_dist_list[1].get_position()
pos_f = ax_map_list[1][0].get_position() 
pos_h = ax_map_list[1][2].get_position() 

cbar_y = pos_f.y0 - 0.065
cbar_height = 0.012
cbar_ax = fig.add_axes([pos_f.x0, cbar_y, pos_h.x1 - pos_f.x0, cbar_height]) 
plt.colorbar(cf_diff_last, cax=cbar_ax, orientation='horizontal')

ax_dist_list[0].set_position([pos_a.x0, pos_b.y0, pos_a.width, pos_b.height])
ax_dist_list[1].set_position([pos_e.x0, cbar_y, pos_e.width, pos_f.y1 - cbar_y])

# ====================== 6. 保存与导出 ======================
save_name = "../illustration/off/FIG3.st_TCC_combined.pdf"
plt.savefig(save_name, dpi=300, bbox_inches='tight')
plt.close()

csv_file = "TS_metrics_table.csv"
with open(csv_file, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Season", "Layer", "Experiment", "TCC", "ACC", "RMSE"])
    writer.writerows(table_metrics_data)
print(f"\n所有数据处理完成！已使用单尾 TCC 并集遮罩导出至 {save_name}，表格至 {csv_file}")