import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import scipy.stats
import warnings
import os
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore")

# ====================== 1. 全局配置 (Nature / JGR Style) ======================
TARGET_VAR = 'pre'
PLOT_SEASONS = ['MAM', 'JJA']
NYEARS = 17

# 学术期刊级绘图全局设置 (无3D透视，清晰的无衬线字体，极简线框)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['font.size'] = 10

# 路径配置
path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"

VAR_CONFIG = {
    'obs_var': 'pre', 'mod_var': 'PRAVG', 'convert': lambda x: x * 86400.0,
    'obs_file_mam': "/share/home/dq135/reference/CN05.1_Pre_1991_2023_MAM_025x025.nc",
    'obs_file_jja': "/share/home/dq135/reference/CN05.1_Pre_1991_2023_JJA_025x025.nc",
    'mod_file_ctl_mam': path + "wrfout_2001-2023_MAM_daymean_nogravel.nc",
    'mod_file_exp_mam': path + "wrfout_2001-2023_MAM_daymean_gravel.nc",
    'mod_file_ctl_jja': path + "wrfout_2001-2023_JJA_daymean_nogravel.nc",
    'mod_file_exp_jja': path + "wrfout_2001-2023_JJA_daymean_gravel.nc",
}

# 饼图类别与 NPG (Nature Publishing Group) 高级配色方案
CATEGORIES = [
    r'Positive TCC, $p < 0.01$',
    r'Positive TCC, $0.01 \leq p < 0.05$',
    r'Positive TCC, $0.05 \leq p < 0.1$',
    r'Positive TCC, $p \geq 0.1$',
    r'Negative TCC'
]
# 色卡：深红 -> 珊瑚红 -> 浅桃红 -> 浅灰 (不显著) -> 湖蓝 (负相关)
COLORS = ['#E64B35', '#F39B7F', '#FDDBC7', '#E0E0E0', '#4DBBD5']

# ====================== 2. 数据处理函数 ======================
def regrid_rcm2rgrid(var2d, lat2d, lon2d, lat1d, lon1d):
    """双线性插值结合最近邻插值填补 NaN"""
    pts = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    vals = var2d.ravel()
    gx, gy = np.meshgrid(lon1d, lat1d)
    interp = griddata(pts, vals, (gx, gy), method='linear')
    nan_mask = np.isnan(interp)
    if nan_mask.any():
        interp[nan_mask] = griddata(pts, vals, (gx, gy), method='nearest')[nan_mask]
    return interp

def get_yearly_mean(file_path, var_name, is_obs, convert_func=None):
    """读取网格数据并计算逐年均值"""
    nc = Dataset(file_path)
    data = nc.variables[var_name][:]
    nc.close()
    if is_obs:
        subset = data[30:30+NYEARS*3]
        yearly = np.nanmean(subset.reshape(NYEARS, 3, subset.shape[1], subset.shape[2]), axis=1)
    else:
        subset = data[:NYEARS*92]
        yearly = np.nanmean(subset.reshape(NYEARS, 92, subset.shape[1], subset.shape[2]), axis=1)
        if convert_func:
            yearly = convert_func(yearly)
    return yearly

def calc_tcc_and_p(obs_3d, mod_3d, mask):
    """
    计算格点的时间相关系数 (TCC) 及其双尾 p-value
    注意：此处严格使用双尾检验，方向判断在后续 categorize_tcc 中完成
    """
    tcc = np.full(mask.shape, np.nan)
    pval = np.full(mask.shape, np.nan)
    
    vy, vx = np.where(mask)
    for y, x in zip(vy, vx):
        o_ts = obs_3d[:, y, x]
        m_ts = mod_3d[:, y, x]
        
        valid = ~np.isnan(o_ts) & ~np.isnan(m_ts)
        if np.sum(valid) > 3:
            r, p = scipy.stats.pearsonr(m_ts[valid], o_ts[valid])
            tcc[y, x] = r
            pval[y, x] = p
            
    return tcc, pval

def categorize_tcc(tcc, pval, mask):
    """根据 TCC 的正负及双尾 P 值进行等级划分，统计网格点占比"""
    valid_grids = mask & ~np.isnan(tcc) & ~np.isnan(pval)
    total_valid = np.sum(valid_grids)
    
    if total_valid == 0:
        return [0]*5
    
    tcc_v = tcc[valid_grids]
    pval_v = pval[valid_grids]
    
    # 核心逻辑：区分正向模拟技巧与负相关
    pos_mask = tcc_v > 0
    neg_mask = tcc_v <= 0
    
    counts = [
        np.sum(pos_mask & (pval_v < 0.01)),                    # 正相关且极显著
        np.sum(pos_mask & (pval_v >= 0.01) & (pval_v < 0.05)), # 正相关且显著
        np.sum(pos_mask & (pval_v >= 0.05) & (pval_v < 0.1)),  # 正相关且边缘显著
        np.sum(pos_mask & (pval_v >= 0.1)),                    # 正相关但不显著
        np.sum(neg_mask)                                       # 负相关 (无论是否显著)
    ]
    
    percentages = [c / total_valid * 100 for c in counts]
    return percentages

# ====================== 3. 主程序 ======================
if __name__ == '__main__':
    # 读取基础网格信息
    f_wrf = Dataset(wrfinput_file)
    lat2d = f_wrf.variables['XLAT'][0, :, :]
    lon2d = f_wrf.variables['XLONG'][0, :, :]
    f_wrf.close()

    f_obs_ref = Dataset(VAR_CONFIG['obs_file_mam'])
    lat1d = f_obs_ref.variables['lat'][:]
    lon1d = f_obs_ref.variables['lon'][:]
    f_obs_ref.close()

    os.makedirs("../illustration/cpl", exist_ok=True)
    
    # 初始化画布 2x2
    fig = plt.figure(figsize=(10, 9))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.15, wspace=0.05)
    axes = [[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(2)]
    
    row_labels = ['Spring', 'Summer']
    col_labels = ['CTL', 'EXP']
    letters = [['a', 'b'], ['c', 'd']]

    for i, season in enumerate(PLOT_SEASONS):
        print(f"---> 开始处理降水 TCC 统计，季节: {season}")
        
        # 1. 加载数据
        obs_17 = get_yearly_mean(VAR_CONFIG[f'obs_file_{season.lower()}'], VAR_CONFIG['obs_var'], True)
        ctl_17_raw = get_yearly_mean(VAR_CONFIG[f'mod_file_ctl_{season.lower()}'], VAR_CONFIG['mod_var'], False, VAR_CONFIG['convert'])
        exp_17_raw = get_yearly_mean(VAR_CONFIG[f'mod_file_exp_{season.lower()}'], VAR_CONFIG['mod_var'], False, VAR_CONFIG['convert'])
        
        # 2. 空间插值对齐
        ctl_17 = np.zeros_like(obs_17)
        exp_17 = np.zeros_like(obs_17)
        
        print("     正在将模式数据插值到观测网格...")
        for yr in range(NYEARS):
            ctl_17[yr] = regrid_rcm2rgrid(ctl_17_raw[yr], lat2d, lon2d, lat1d, lon1d)
            exp_17[yr] = regrid_rcm2rgrid(exp_17_raw[yr], lat2d, lon2d, lat1d, lon1d)
        
        # 提取陆地观测掩膜 (针对均值场存在有效值的区域)
        mask_obs = ~np.isnan(np.nanmean(obs_17, axis=0))
        
        # 3. 计算 TCC 与 P-value
        print("     正在计算 TCC 与双尾显著性...")
        tcc_ctl, pval_ctl = calc_tcc_and_p(obs_17, ctl_17, mask_obs)
        tcc_exp, pval_exp = calc_tcc_and_p(obs_17, exp_17, mask_obs)
        
        # 4. 统计全国网格点各类占比
        pct_ctl = categorize_tcc(tcc_ctl, pval_ctl, mask_obs)
        pct_exp = categorize_tcc(tcc_exp, pval_exp, mask_obs)
        
        # 5. 组图绘制逻辑
        for j, (pct_data, exp_name) in enumerate(zip([pct_ctl, pct_exp], col_labels)):
            ax = axes[i][j]
            
            # 绘制 Donut Chart (环形图)
            wedges, texts, autotexts = ax.pie(
                pct_data, 
                colors=COLORS, 
                autopct='%1.1f%%', 
                startangle=90,          # 从正上方开始顺时针绘制
                pctdistance=0.75,       # 百分比数值的位置半径
                counterclock=False,
                wedgeprops=dict(width=0.4, edgecolor='white', linewidth=1.5) # width控制环宽，产生甜甜圈效果
            )
            
            # 优化标签文本字体与可见度
            for autotext in autotexts:
                autotext.set_fontsize(9)
                autotext.set_weight('normal')
                autotext.set_color('black')
                
                # 若某类占比太小(<3%)，隐藏文本避免数字挤压重叠
                text_val = float(autotext.get_text().strip('%'))
                if text_val < 3.0:
                    autotext.set_visible(False)

            # 顶部标题 (仅第一行显示 CTL/EXP)
            if i == 0:
                ax.set_title(exp_name, fontsize=13, fontweight='normal', pad=15)
            
            # 侧边季节标签 (仅第一列显示 MAM/JJA)
            if j == 0:
                ax.text(-1.4, 0, row_labels[i], fontsize=13, fontweight='normal', 
                        va='center', ha='center', rotation=90)
                
            # 子图序号 a/b/c/d
            ax.text(-1.2, 1.1, f"({letters[i][j]})", fontsize=12, fontweight='normal')
            
            # 保证环形为正圆
            ax.axis('equal') 

    # ====================== 4. 统一图例与保存 ======================
    # 在底部统一添加无边框图例
    fig.legend(
        wedges, CATEGORIES, 
        loc='lower center', 
        ncol=3, 
        fontsize=11, 
        frameon=False, 
        bbox_to_anchor=(0.5, 0.02)
    )
    
    # 底部留白以容纳图例
    plt.subplots_adjust(bottom=0.15)
    
    save_path = "../illustration/cpl/FIG8.pre_TCC_Donut.pdf"
    plt.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"\n✅ 降水 TCC 统计环形图已成功保存至：{save_path}")