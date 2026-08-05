import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import warnings
warnings.filterwarnings("ignore")
from scipy.stats import pearsonr, gaussian_kde, ttest_rel, wilcoxon
from scipy.signal import detrend
import os
import concurrent.futures

# ====================== 1. 全局配置与客制化控制 ======================
# -------------------------------------------------------------------
SC_EXTENT = [104, 118, 20, 27] 
KDE_BANDWIDTH = 0.1
NUM_WORKERS = 10 # 并行进程数
# -------------------------------------------------------------------

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"

VAR_CONFIG = {
    't2m': {
        'obs_var': 'tm', 'mod_var': 'AT2M', 'convert': lambda x: x - 273.15,
        'obs_file_mam': "/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc",
        'obs_file_jja': "/share/home/dq135/reference/CN05.1_Tm_1991_2023_JJA_025x025.nc",
        'mod_file_ctl_mam': path + "wrfout_2001-2023_MAM_daymean_nogravel.nc",
        'mod_file_exp_mam': path + "wrfout_2001-2023_MAM_daymean_gravel.nc",
        'mod_file_ctl_jja': path + "wrfout_2001-2023_JJA_daymean_nogravel.nc",
        'mod_file_exp_jja': path + "wrfout_2001-2023_JJA_daymean_gravel.nc",
        'title': 'T2m'
    },
    'pre': {
        'obs_var': 'pre', 'mod_var': 'PRAVG', 'convert': lambda x: x * 86400.0,
        'obs_file_mam': "/share/home/dq135/reference/CN05.1_Pre_1991_2023_MAM_025x025.nc",
        'obs_file_jja': "/share/home/dq135/reference/CN05.1_Pre_1991_2023_JJA_025x025.nc",
        'mod_file_ctl_mam': path + "wrfout_2001-2023_MAM_daymean_nogravel.nc",
        'mod_file_exp_mam': path + "wrfout_2001-2023_MAM_daymean_gravel.nc",
        'mod_file_ctl_jja': path + "wrfout_2001-2023_JJA_daymean_nogravel.nc",
        'mod_file_exp_jja': path + "wrfout_2001-2023_JJA_daymean_gravel.nc",
        'title': 'Pre'
    }
}
variables = ['t2m', 'pre']
seasons = ['MAM', 'JJA']
nyears = 17 

# ====================== 2. 核心计算函数 ======================
def regrid_rcm2rgrid(var2d, lat2d, lon2d, lat1d, lon1d):
    pts = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    vals = var2d.ravel()
    gx, gy = np.meshgrid(lon1d, lat1d)
    interp = griddata(pts, vals, (gx, gy), method='linear')
    nan_mask = np.isnan(interp)
    if nan_mask.any():
        interp[nan_mask] = griddata(pts, vals, (gx, gy), method='nearest')[nan_mask]
    return interp

def get_yearly_mean(file_path, var_name, is_obs, convert_func=None):
    nc = Dataset(file_path)
    data = nc.variables[var_name][:]
    nc.close()
    if is_obs:
        subset = data[30:30+nyears*3]
        yearly = np.nanmean(subset.reshape(nyears, 3, subset.shape[1], subset.shape[2]), axis=1)
    else:
        subset = data[:nyears*92]
        yearly = np.nanmean(subset.reshape(nyears, 92, subset.shape[1], subset.shape[2]), axis=1)
        if convert_func: yearly = convert_func(yearly)
    return yearly

def spatial_acc(anom1, anom2, wgt):
    v1, v2, w = anom1.ravel(), anom2.ravel(), wgt.ravel()
    valid = ~np.isnan(v1) & ~np.isnan(v2) & ~np.isnan(w)
    if np.sum(valid) < 3: return np.nan
    v1, v2, w = v1[valid], v2[valid], w[valid]
    mean1 = np.average(v1, weights=w); mean2 = np.average(v2, weights=w)
    cov = np.average((v1 - mean1) * (v2 - mean2), weights=w)
    var1 = np.average((v1 - mean1)**2, weights=w); var2 = np.average((v2 - mean2)**2, weights=w)
    return cov / np.sqrt(var1 * var2)

def calc_row_tcc(args):
    # 接收条件变量 is_t2m
    i, c_mat, e_mat, o_mat, m_row, is_t2m = args
    t1 = np.full(m_row.shape[0], np.nan); t2 = np.full(m_row.shape[0], np.nan)
    for j in range(m_row.shape[0]):
        if m_row[j]:
            c_ts, e_ts, o_ts = c_mat[:, j], e_mat[:, j], o_mat[:, j]
            v1 = ~np.isnan(c_ts) & ~np.isnan(o_ts); v2 = ~np.isnan(e_ts) & ~np.isnan(o_ts)
            
            # 气温数据先去除线性趋势，降水数据直接计算
            if v1.sum() >= 3: 
                if is_t2m:
                    t1[j], _ = pearsonr(detrend(c_ts[v1]), detrend(o_ts[v1]))
                else:
                    t1[j], _ = pearsonr(c_ts[v1], o_ts[v1])
            if v2.sum() >= 3: 
                if is_t2m:
                    t2[j], _ = pearsonr(detrend(e_ts[v2]), detrend(o_ts[v2]))
                else:
                    t2[j], _ = pearsonr(e_ts[v2], o_ts[v2])
                
    return i, t1, t2

def get_sig_asterisks(p_val):
    if p_val < 0.01: return '**'
    elif p_val < 0.05: return '*'
    else: return 'ns'

def format_pval(p):
    return "<0.001" if p < 0.001 else f"{p:.3f}"

def add_significance_line(ax, x1, x2, y_max, p_val):
    y, h = y_max + 0.03, 0.02
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.2, c='k')
    p_str = format_pval(p_val)
    sig_text = f"p={p_str}"
    ax.text((x1+x2)*.5, y+h, sig_text, ha='center', va='bottom', fontsize=10)

# ====================== 3. 主绘图流程 ======================
if __name__ == '__main__':
    f_wrf = Dataset(wrfinput_file)
    lat2d = f_wrf.variables['XLAT'][0, :, :]; lon2d = f_wrf.variables['XLONG'][0, :, :]
    f_wrf.close()

    f_obs_ref = Dataset(VAR_CONFIG['t2m']['obs_file_mam'])
    lat1d = f_obs_ref.variables['lat'][:]; lon1d = f_obs_ref.variables['lon'][:]
    f_obs_ref.close()

    lon_grid, lat_grid = np.meshgrid(lon1d, lat1d)
    mask_sc = (lon_grid >= SC_EXTENT[0]) & (lon_grid <= SC_EXTENT[1]) & (lat_grid >= SC_EXTENT[2]) & (lat_grid <= SC_EXTENT[3])
    wgt_2d = np.outer(np.cos(lat1d * np.pi / 180), np.ones_like(lon1d))

    os.makedirs("../illustration/cpl", exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    plt.subplots_adjust(wspace=0.25, hspace=0.25)

    # 定义所有子图标签内容
    panel_texts = [
        "(a) T2m (Spring)", "(b) T2m (Summer)", "(c)",
        "(d) Precip (Spring)", "(e) Precip (Summer)", "(f)"
    ]

    for row_idx, var_name in enumerate(variables):
        cfg = VAR_CONFIG[var_name]
        acc_data = {'MAM': {}, 'JJA': {}}
        is_t2m = (var_name == 't2m') # 判断当前是否处理温度数据
        
        for col_idx, season in enumerate(seasons):
            print(f"Processing {var_name} {season}...")
            obs_17 = get_yearly_mean(cfg[f'obs_file_{season.lower()}'], cfg['obs_var'], True)
            ctl_17_raw = get_yearly_mean(cfg[f'mod_file_ctl_{season.lower()}'], cfg['mod_var'], False, cfg['convert'])
            exp_17_raw = get_yearly_mean(cfg[f'mod_file_exp_{season.lower()}'], cfg['mod_var'], False, cfg['convert'])
            
            ctl_17 = np.zeros((nyears, lat1d.size, lon1d.size))
            exp_17 = np.zeros((nyears, lat1d.size, lon1d.size))
            for yr in range(nyears):
                ctl_17[yr] = regrid_rcm2rgrid(ctl_17_raw[yr], lat2d, lon2d, lat1d, lon1d)
                exp_17[yr] = regrid_rcm2rgrid(exp_17_raw[yr], lat2d, lon2d, lat1d, lon1d)
                ctl_17[yr] = np.where(~np.isnan(obs_17[yr]), ctl_17[yr], np.nan)
                exp_17[yr] = np.where(~np.isnan(obs_17[yr]), exp_17[yr], np.nan)
                
            valid_mask = ~np.isnan(np.nanmean(obs_17, axis=0)) & mask_sc
            
            # MAE 计算
            mae_ctl = np.nanmean(np.abs(ctl_17[:, valid_mask] - obs_17[:, valid_mask]))
            mae_exp = np.nanmean(np.abs(exp_17[:, valid_mask] - obs_17[:, valid_mask]))
            
            # TCC 并行计算 (传入 is_t2m 参数)
            TCC1 = np.full_like(valid_mask, np.nan, dtype=float); TCC2 = np.full_like(valid_mask, np.nan, dtype=float)
            row_args = [(i, ctl_17[:, i, :], exp_17[:, i, :], obs_17[:, i, :], valid_mask[i, :], is_t2m) for i in range(lat1d.size)]
            with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_WORKERS) as ex:
                for i, t1, t2 in ex.map(calc_row_tcc, row_args):
                    TCC1[i, :], TCC2[i, :] = t1, t2
            
            common_valid = valid_mask & ~np.isnan(TCC1) & ~np.isnan(TCC2)
            tcc1_v = TCC1[common_valid] 
            tcc2_v = TCC2[common_valid]
            
            # ACC 计算 (基于原时间序列去除长期平均态得到Anomaly)
            if is_t2m:
                # 气温：在时间轴 (17年) 上逐像元去除线性趋势，detrend 处理后均值自动为 0，即为异常场
                obs_anom = np.full_like(obs_17, np.nan)
                ctl_anom = np.full_like(ctl_17, np.nan)
                exp_anom = np.full_like(exp_17, np.nan)
                for i in range(lat1d.size):
                    for j in range(lon1d.size):
                        if valid_mask[i, j]:
                            obs_anom[:, i, j] = detrend(obs_17[:, i, j])
                            ctl_anom[:, i, j] = detrend(ctl_17[:, i, j])
                            exp_anom[:, i, j] = detrend(exp_17[:, i, j])
            else:
                # 降水：仅减去多年平均态得到异常场，不进行去趋势处理
                obs_anom = obs_17 - np.nanmean(obs_17, axis=0)
                ctl_anom = ctl_17 - np.nanmean(ctl_17, axis=0)
                exp_anom = exp_17 - np.nanmean(exp_17, axis=0)
            wgt_acc = np.where(valid_mask, wgt_2d, np.nan)
            acc_data[season]['CTL'] = [spatial_acc(obs_anom[y], ctl_anom[y], wgt_acc) for y in range(nyears)]
            acc_data[season]['EXP'] = [spatial_acc(obs_anom[y], exp_anom[y], wgt_acc) for y in range(nyears)]

            # KDE 绘图
            ax = axes[row_idx, col_idx]
            if len(tcc1_v) > 1:
                kde1 = gaussian_kde(tcc1_v, bw_method=KDE_BANDWIDTH); kde2 = gaussian_kde(tcc2_v, bw_method=KDE_BANDWIDTH)
                xr = np.linspace(-1, 1, 200)
                ax.plot(xr, kde1(xr), color='royalblue', label=f'CTL CPL (MAE={mae_ctl:.2f})', lw=2)
                ax.plot(xr, kde2(xr), color='indianred', label=f'EXP CPL (MAE={mae_exp:.2f})', lw=2, ls='--')
                ax.fill_between(xr, kde1(xr), alpha=0.15, color='royalblue'); ax.fill_between(xr, kde2(xr), alpha=0.15, color='indianred')
                
                # TCC 统一使用双尾配对 t 检验
                t_pval = ttest_rel(tcc1_v, tcc2_v)[1]
                p_str = format_pval(t_pval)
                # ax.text(0.04, 0.78, f'p={p_str} ({get_sig_asterisks(t_pval)})', transform=ax.transAxes, fontsize=10)
                ax.text(0.04, 0.78, f'p={p_str}', transform=ax.transAxes, fontsize=10)
            
            ax.set_xlim(-0.8, 1); ax.set_xlabel('TCC', fontsize=11)
            if col_idx == 0: ax.set_ylabel('Probability Density(KDE)', fontsize=11)
            ax.legend(loc='upper left', frameon=False, fontsize=9)
            ax.grid(True, ls='--', alpha=0.4)
            # 恢复框线
            for s in ax.spines.values(): s.set_visible(True)

        # 箱线图
        ax_box = axes[row_idx, 2]
        d_box = [acc_data['MAM']['CTL'], acc_data['MAM']['EXP'], acc_data['JJA']['CTL'], acc_data['JJA']['EXP']]
        pos = [1, 2, 4, 5]
        bp = ax_box.boxplot(d_box, positions=pos, widths=0.6, patch_artist=True, medianprops=dict(color="black", linewidth=1.5))
        colors = ['royalblue', 'indianred', 'royalblue', 'indianred']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color); patch.set_alpha(0.6)
        for i, d in enumerate(d_box):
            ax_box.scatter(np.random.normal(pos[i], 0.05, len(d)), d, color=colors[i], s=15, alpha=0.7)
        
        # ACC 检验：按变量类型决定使用 Wilcoxon 或 ttest
        pm = ttest_rel(acc_data['MAM']['CTL'], acc_data['MAM']['EXP'])[1]
        pj = ttest_rel(acc_data['JJA']['CTL'], acc_data['JJA']['EXP'])[1]
            
        add_significance_line(ax_box, 1, 2, max(max(d_box[0]), max(d_box[1])), pm)
        add_significance_line(ax_box, 4, 5, max(max(d_box[2]), max(d_box[3])), pj)
        
        ax_box.set_xticks([1.5, 4.5]); ax_box.set_xticklabels(['Spring', 'Summer'], fontsize=11)
        ax_box.set_ylabel('ACC', fontsize=11)
        
        from matplotlib.lines import Line2D
        # 图例仅在子图 c (row_idx == 0) 显示
        if row_idx == 0:
            # bbox_to_anchor=(0.2, 0.08), 
            ax_box.legend([Line2D([0],[0], color='royalblue', lw=4, alpha=0.6), Line2D([0],[0], color='indianred', lw=4, alpha=0.6)], 
                          ['CTL CPL', 'EXP CPL'], loc='best', frameon=False)
        ax_box.grid(True, axis='y', ls='--', alpha=0.4)
        for s in ax_box.spines.values(): s.set_visible(True)

        # 设置外部左上角复合标签 (a)-(f)
        for i, ax_sub in enumerate(axes[row_idx, :]):
            idx = row_idx * 3 + i
            ax_sub.text(0.01, 0.1, panel_texts[idx], transform=ax_sub.transAxes, fontsize=10, fontweight='normal', va='bottom', ha='left')

    plt.savefig("../illustration/cpl/SouthChina_TCC_KDE_ACC_Boxplot_Stats.pdf", dpi=600, bbox_inches='tight')
    print("Done! Plot saved.")