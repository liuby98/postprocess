import warnings
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# 忽略全 NaN 切片求均值时产生的无害警告
warnings.filterwarnings('ignore', category=RuntimeWarning)

# ================= 1. 设置路径和读取数据 =================
path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
f1_path = path + "colmoff_2001-2023_monmean_nogravel.nc"
f2_path = path + "colmoff_2001-2023_monmean_gravel.nc"
f_input_path = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"

print("正在读取数据...")
ds_ctl = xr.open_dataset(f1_path)
ds_exp = xr.open_dataset(f2_path)
ds_in = xr.open_dataset(f_input_path)

lat2d = ds_in['XLAT'].values[0, :, :]
lon2d = ds_in['XLONG'].values[0, :, :]

# NCL 脚本中使用的冰透镜变量名
ICE_VAR = 'f_wice_soisno'  

try:
    # 提取 5:15 索引 (即第 1-10 层)
    i_ctl = ds_ctl[ICE_VAR].values[:, :, :, 5:15]
    i_exp = ds_exp[ICE_VAR].values[:, :, :, 5:15]
except KeyError:
    print(f"错误: 找不到变量 {ICE_VAR}，请检查 nc 文件中的真实变量名！")
    exit()

# 将异常值 / 缺测值剔除 (< -1000)，直接使用原始 kg/m2，不减 273.15
i_ctl = np.where(i_ctl < -1000., np.nan, i_ctl)
i_exp = np.where(i_exp < -1000., np.nan, i_exp)

# ================= 2. 获取 3, 4, 5 月份索引及 TP 掩膜 =================
mar_idx = [6 * i + 0 for i in range(17)]
apr_idx = [6 * i + 1 for i in range(17)]
may_idx = [6 * i + 2 for i in range(17)]

# 设定青藏高原 (Tibetan Plateau) 区域的经纬度范围
lat_min, lat_max = 26, 39
lon_min, lon_max = 73, 105
tp_mask = (lat2d >= lat_min) & (lat2d <= lat_max) & (lon2d >= lon_min) & (lon2d <= lon_max)

# 2. 设定需要抠除的特定区域 (80-90E, 30-33N)
ex_lat_min, ex_lat_max = 30, 33
ex_lon_min, ex_lon_max = 80, 90
exclude_mask = (lat2d >= ex_lat_min) & (lat2d <= ex_lat_max) & (lon2d >= ex_lon_min) & (lon2d <= ex_lon_max)

# 3. 最终掩膜：在青藏高原内，且（&）不在（~）抠除区域内
final_tp_mask = tp_mask & (~exclude_mask)

# ================= 3. 定义区域平均计算函数 =================
def calc_tp_profile_ice(data, time_indices, mask):
    """提取特定月份17年均值，计算掩膜内空间平均"""
    data_time_mean = np.nanmean(data[time_indices, :, :, :], axis=0)
    # 扩展 2D mask 到 3D
    mask_3d = np.repeat(mask[:, :, np.newaxis], data_time_mean.shape[-1], axis=2)
    # 只针对掩膜区域取值并求空间均值
    tp_spatial_mean = np.nanmean(np.where(mask_3d, data_time_mean, np.nan), axis=(0, 1))
    return tp_spatial_mean

print("正在计算青藏高原区域平均...")
# 计算各月 CTL 冰透镜剖面
ctl_mar = calc_tp_profile_ice(i_ctl, mar_idx, final_tp_mask)
ctl_apr = calc_tp_profile_ice(i_ctl, apr_idx, final_tp_mask)
ctl_may = calc_tp_profile_ice(i_ctl, may_idx, final_tp_mask)

# 计算各月 EXP 冰透镜剖面
exp_mar = calc_tp_profile_ice(i_exp, mar_idx, final_tp_mask)
exp_apr = calc_tp_profile_ice(i_exp, apr_idx, final_tp_mask)
exp_may = calc_tp_profile_ice(i_exp, may_idx, final_tp_mask)

# ================= 4. 绘制剖面曲线图 =================
print("正在绘制曲线图...")
layers = np.arange(1, 11)  # 对应土壤 1-10 层

fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharey=True)
plt.subplots_adjust(hspace=0.3)
# 'June', 'July', 'August'
months = ['March', 'April', 'May']
ctl_data = [ctl_mar, ctl_apr, ctl_may]
exp_data = [exp_mar, exp_apr, exp_may]

# --- 核心修改部分：计算统一的 X 轴范围 ---
diff_data = [exp_data[i] - ctl_data[i] for i in range(3)]
diff_min = min([np.nanmin(d) for d in diff_data])
diff_max = max([np.nanmax(d) for d in diff_data])

# 稍微加一点 padding 让曲线不要紧贴边框
padding = (diff_max - diff_min) * 0.1
xlim_min = diff_min - padding 
xlim_max = diff_max + padding

# 如果需要强制包含0刻度线，可以取消下面这两行的注释
# xlim_min = min(xlim_min, -0.01)
# xlim_max = max(xlim_max, 0.01)
# ----------------------------------------

# 第一排: CTL 和 EXP 的原始土壤冰含量
for i in range(3):
    ax = axes[0, i]
    ax.plot(ctl_data[i], layers, marker='o', label='CTL', color='#1f77b4', linewidth=2)
    ax.plot(exp_data[i], layers, marker='s', label='EXP', color='#d62728', linewidth=2)
    
    if i == 0:
        ax.invert_yaxis()  # 反转 Y 轴，让表层在上面
        ax.set_ylabel('Soil Layer', fontsize=12)
    
    ax.set_title(f'{months[i]} Soil Ice Content', fontsize=13, fontweight='bold')
    ax.set_xlabel('Ice Content ($kg/m^2$)', fontsize=12)
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.6)

# 第二排: EXP - CTL 的差值变化
for i in range(3):
    ax = axes[1, i]
    diff = diff_data[i]
    
    ax.plot(diff, layers, marker='^', label='EXP - CTL', color='green', linewidth=2)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1) 
    
    # --- 统一设置第二排的 X 轴范围 ---
    ax.set_xlim(xlim_min, xlim_max)
    
    if i == 0:
        ax.set_ylabel('Soil Layer', fontsize=12)
        
    ax.set_title(f'{months[i]} Ice Diff (EXP-CTL)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Ice Diff ($kg/m^2$)', fontsize=12)
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.6)

fig.suptitle('Month Mean Soil Ice Lens Profiles over Qinghai-Tibet Plateau (26°N-39°N, 73°E-105°E)', 
             fontsize=16, fontweight='bold', y=0.96)

out_name = "TP_soil_ice_profile.pdf"
plt.savefig(out_name, dpi=300, bbox_inches='tight')
print(f"绘图完成，结果已保存为: {out_name}")

ds_ctl.close()
ds_exp.close()
ds_in.close()