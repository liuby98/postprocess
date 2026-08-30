import warnings
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
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

# 获取经纬度二维网格数据
lat2d = ds_in['XLAT'].values[0, :, :]
lon2d = ds_in['XLONG'].values[0, :, :]

# 【重要提醒】: 根据液态水变量的真实名称进行修改
# NCL脚本中误用了 f_t_soisno，这里假设液态水变量名为 f_wliq_soisno
VAR_NAME = 'f_wliq_soisno' 
# 若确定你的液态水存的就是 f_t_soisno，请将其改回 'f_t_soisno'

try:
    water_ctl = ds_ctl[VAR_NAME].values[:, :, :, 5:15]
    water_exp = ds_exp[VAR_NAME].values[:, :, :, 5:15]
except KeyError:
    print(f"警告：找不到变量 {VAR_NAME}，请检查 nc 文件中的真实变量名！")
    # 作为退回选项，暂时用 f_t_soisno 代替演示
    water_ctl = ds_ctl['f_wliq_soisno'].values[:, :, :, 5:15]
    water_exp = ds_exp['f_wliq_soisno'].values[:, :, :, 5:15]

# 将异常缺测值剔除 (< -1000)
water_ctl = np.where(water_ctl < -1000., np.nan, water_ctl)
water_exp = np.where(water_exp < -1000., np.nan, water_exp)

# ================= 2. 获取 3, 4, 5 月份索引及 TP 掩膜 =================
mar_idx = [6 * i + 0 for i in range(17)]
apr_idx = [6 * i + 1 for i in range(17)]
may_idx = [6 * i + 2 for i in range(17)]

# 设定青藏高原 (Tibetan Plateau) 区域的经纬度范围
lat_min, lat_max = 26, 39
lon_min, lon_max = 73, 105

# 生成二维布尔掩膜
tp_mask = (lat2d >= lat_min) & (lat2d <= lat_max) & (lon2d >= lon_min) & (lon2d <= lon_max)

# 2. 设定需要抠除的特定区域 (80-90E, 30-33N)
ex_lat_min, ex_lat_max = 30, 33
ex_lon_min, ex_lon_max = 80, 90
exclude_mask = (lat2d >= ex_lat_min) & (lat2d <= ex_lat_max) & (lon2d >= ex_lon_min) & (lon2d <= ex_lon_max)

# 3. 最终掩膜：在青藏高原内，且（&）不在（~）抠除区域内
final_tp_mask = tp_mask & (~exclude_mask)

# ================= 3. 定义区域平均计算函数 =================
def calc_tp_profile(data, time_indices, mask):
    """
    提取特定月份的17年数据，计算时间均值后，再求掩膜内的区域空间平均
    返回: shape 为 (10,) 的1维数组
    """
    # 提取时间切片并求时间平均 (结果形状: lat, lon, level)
    data_time_mean = np.nanmean(data[time_indices, :, :, :], axis=0)
    
    # 将二维 mask (lat, lon) 扩展到三维 (lat, lon, level) 以便与数据对齐
    mask_3d = np.repeat(mask[:, :, np.newaxis], data_time_mean.shape[-1], axis=2)
    
    # 仅保留掩膜(TP)区域的格点，求空间均值 (剔除 NaN)
    tp_spatial_mean = np.nanmean(np.where(mask_3d, data_time_mean, np.nan), axis=(0, 1))
    
    # 注意：此处不再减去 273.15，保持原本的水含量 kg/m2 的量级
    return tp_spatial_mean

print("正在计算区域平均...")
# 计算各月 CTL 的剖面
ctl_mar = calc_tp_profile(water_ctl, mar_idx, final_tp_mask)
ctl_apr = calc_tp_profile(water_ctl, apr_idx, final_tp_mask)
ctl_may = calc_tp_profile(water_ctl, may_idx, final_tp_mask)

# 计算各月 EXP 的剖面
exp_mar = calc_tp_profile(water_exp, mar_idx, final_tp_mask)
exp_apr = calc_tp_profile(water_exp, apr_idx, final_tp_mask)
exp_may = calc_tp_profile(water_exp, may_idx, final_tp_mask)

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

# 绘制上排 (CTL 和 EXP 的原始液态水)
for i in range(3):
    ax = axes[0, i]
    ax.plot(ctl_data[i], layers, marker='o', label='CTL', color='#1f77b4', linewidth=2)
    ax.plot(exp_data[i], layers, marker='s', label='EXP', color='#d62728', linewidth=2)
    
    if i == 0:
        ax.invert_yaxis()  # 反转 Y 轴，让第一层在最上面
        ax.set_ylabel('Soil Layer', fontsize=12)
    
    ax.set_title(f'{months[i]} Liquid Water', fontsize=13, fontweight='bold')
    ax.set_xlabel('Liquid Water (kg/m²)', fontsize=12)
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.6)

# 绘制下排 (EXP - CTL 差值图)
for i in range(3):
    ax = axes[1, i]
    diff = diff_data[i]
    
    ax.plot(diff, layers, marker='^', label='EXP - CTL', color='green', linewidth=2)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1) # 添加 0 变化参考线
    
    # --- 统一设置第二排的 X 轴范围 ---
    ax.set_xlim(xlim_min, xlim_max)

    if i == 0:
        ax.set_ylabel('Soil Layer', fontsize=12)
        
    ax.set_title(f'{months[i]} Water Diff (EXP-CTL)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Liquid Water Diff (kg/m²)', fontsize=12)
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.6)

# 添加总标题
fig.suptitle('Month Mean Soil Liquid Water Profiles over Qinghai-Tibet Plateau (26°N-39°N, 73°E-105°E)', 
             fontsize=16, fontweight='bold', y=0.96)

out_name = "TP_liquid_water_profile.pdf"
plt.savefig(out_name, dpi=300, bbox_inches='tight')
print(f"绘图完成，结果已保存为: {out_name}")

# 关闭 Dataset 释放资源
ds_ctl.close()
ds_exp.close()
ds_in.close()