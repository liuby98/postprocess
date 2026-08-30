import os
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import RegularGridInterpolator
import warnings
warnings.filterwarnings("ignore")

# ====================== 1. 基础环境与路径设置 ======================
input_file = "/share/home/dq013/zhwei/colm/data/CoLMrawdata/soil/vf_gravels_s.nc" 
cn05_file  = "/share/home/dq117/CN05.1/CN05.1_Tm_2020_daily_025x025.nc"
out_dir    = "./figs_gravel_statistics"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# ====================== 2. 读取 CN05.1 提取中国陆地掩膜 ======================
print("正在读取 CN05.1 掩膜数据...")
cn05 = Dataset(cn05_file)
obs_cn05 = cn05.variables['tm'][0, :, :] 
lat_cn = cn05.variables['lat'][:]
lon_cn = cn05.variables['lon'][:]
cn05.close()

if hasattr(obs_cn05, 'mask'):
    obs_cn05 = np.ma.filled(obs_cn05, np.nan)
obs_cn05 = np.where(np.abs(obs_cn05) > 1000, np.nan, obs_cn05)
mask_cn05 = ~np.isnan(obs_cn05)

if lat_cn[0] > lat_cn[-1]:
    lat_cn = lat_cn[::-1]
    mask_cn05 = mask_cn05[::-1, :]
if lon_cn[0] > lon_cn[-1]:
    lon_cn = lon_cn[::-1]
    mask_cn05 = mask_cn05[:, ::-1]

mask_interp_func = RegularGridInterpolator(
    (lat_cn, lon_cn), mask_cn05.astype(float), 
    method='nearest', bounds_error=False, fill_value=0.0
)

# ====================== 3. 读取高分辨数据坐标与生成掩膜 ======================
print(f"正在读取目标数据坐标: {input_file}")
ds = Dataset(input_file)

if 'longitude' in ds.variables:
    lon_all = ds.variables['longitude'][:]
    lat_all = ds.variables['latitude'][:]
else:
    lon_all = np.linspace(-180, 180, 86400)
    lat_all = np.linspace(90, -90, 43200)

lon_idx = np.where((lon_all >= 70) & (lon_all <= 135))[0]
lat_idx = np.where((lat_all >= 15) & (lat_all <= 55))[0]

lat_start, lat_end = np.min(lat_idx), np.max(lat_idx) + 1
lon_start, lon_end = np.min(lon_idx), np.max(lon_idx) + 1

# 采样提取
stride = 3 
lon_subset = lon_all[lon_start:lon_end:stride]
lat_subset = lat_all[lat_start:lat_end:stride]
lon_grid, lat_grid = np.meshgrid(lon_subset, lat_subset)

print("正在计算中国陆地高分辨掩膜...")
pts = np.stack((lat_grid, lon_grid), axis=-1)
mask_china_grid = mask_interp_func(pts) > 0.5

# ====================== 4. 第一遍循环：计算 1-8 层垂直平均值 ======================
sum_gravel = np.zeros_like(lon_grid, dtype=np.float32)
count_gravel = np.zeros_like(lon_grid, dtype=np.float32)

for i in range(1, 9):
    var_name = f'vf_gravels_s_l{i}'
    print(f"第一阶段：正在读取第 {i} 层数据计算均值...")
    data = ds.variables[var_name][lat_start:lat_end:stride, lon_start:lon_end:stride]
    
    if hasattr(data, 'mask'):
        data = np.ma.filled(data, np.nan)
    data = np.where(data > 1000, np.nan, data)
    data = np.where(data < 0, np.nan, data)
    
    valid_mask = ~np.isnan(data)
    sum_gravel[valid_mask] += data[valid_mask]
    count_gravel[valid_mask] += 1

with np.errstate(divide='ignore', invalid='ignore'):
    mean_gravel = np.where(count_gravel > 0, sum_gravel / count_gravel, np.nan)

# ====================== 5. 动态生成全图 5 大区域掩膜 (核心逻辑) ======================
print("正在根据 0.3 阈值和气候经纬度，动态划分全国 5 大区域...")

# 数据驱动的阈值界线
mask_low = mean_gravel < 0.3
mask_high = mean_gravel >= 0.3

# 宏观气候地理界线 (以 105°E 和 40°N 划分东西部与西北)
mask_west = lon_grid < 105
mask_east = lon_grid >= 105
mask_north = lat_grid >= 40
mask_south = lat_grid < 40

# 严格交叉组合，确保所有格点 100% 被分配且不重叠
region_masks = {
    '1. Desert Areas (Gravel content < 0.3)': mask_china_grid & mask_low & mask_west & mask_north,
    '2. Northern Tibetan Plateau lake region (Gravel content < 0.3)': mask_china_grid & mask_low & mask_west & mask_south,
    '3. Eastern Plains & Coast (Gravel content < 0.3)': mask_china_grid & mask_low & mask_east,
    '4. Qinghai-Tibet Plateau (Gravel content $\geq$ 0.3)': mask_china_grid & mask_high & mask_west,
    '5. Greater Khingan Mts - Songnen Plain transition zone (Gravel content $\geq$ 0.3)': mask_china_grid & mask_high & mask_east
}

# ====================== 6. 第二遍循环：提取逐层分区统计数据 ======================
df_list = []
mean_data_list = []

for i in range(1, 9):
    var_name = f'vf_gravels_s_l{i}'
    print(f"第二阶段：正在提取第 {i} 层的分区统计数据...")
    data = ds.variables[var_name][lat_start:lat_end:stride, lon_start:lon_end:stride]
    
    if hasattr(data, 'mask'):
        data = np.ma.filled(data, np.nan)
    data = np.where(data > 1000, np.nan, data)
    data = np.where(data < 0, np.nan, data)
    
    for r_name, r_mask in region_masks.items():
        r_data = data[r_mask]
        r_data = r_data[~np.isnan(r_data)]
        
        # 计算该区域所有有效网格的绝对均值
        exact_mean = np.mean(r_data) if len(r_data) > 0 else np.nan
        mean_data_list.append({'Layer': f"Layer {i}", 'Region': r_name, 'Mean_Fraction': exact_mean})
        
        # 降采样作图以防内存爆炸 (保持分布特征不变)
        # if len(r_data) > 5000:
        #     r_data = np.random.choice(r_data, 5000, replace=False)
            
        df_list.append(pd.DataFrame({'Layer': i, 'Region': r_name, 'Gravel Fraction': r_data}))

ds.close()
print("正在合并全量数据...")
df_all = pd.concat(df_list, ignore_index=True)
df_means = pd.DataFrame(mean_data_list).pivot(index='Region', columns='Layer', values='Mean_Fraction')

# 配色：前3个低值区用冷色系(蓝)，后2个高值区用暖色系(红/橙)
sns.set_theme(style="whitegrid")
custom_palette = {
    '1. Desert Areas (Gravel content < 0.3)': '#85C1E9',          
    '2. Northern Tibetan Plateau lake region (Gravel content < 0.3)': '#3498DB',      
    '3. Eastern Plains & Coast (Gravel content < 0.3)': '#21618C',
    '4. Qinghai-Tibet Plateau (Gravel content $\geq$ 0.3)': '#E74C3C',
    '5. Greater Khingan Mts - Songnen Plain transition zone (Gravel content $\geq$ 0.3)': '#E67E22'
}

# ====================== 7. 绘制 分组箱线图 ======================
print("正在绘制 5 区域箱线图...")
fig, ax = plt.subplots(figsize=(14, 7))
sns.boxplot(
    data=df_all, x='Layer', y='Gravel Fraction', hue='Region', 
    palette=custom_palette, showfliers=False, ax=ax, width=0.7, linewidth=1.2
)

# 添加 0.3 基准线
ax.axhline(0.3, color='gray', linestyle='--', linewidth=1.5, zorder=0, label='Threshold (0.3)')

ax.set_title('Gravel Volumetric Fraction Distribution across Layers', fontsize=15, fontweight='bold', pad=15)
ax.set_xlabel('Soil Layer', fontsize=12)
ax.set_ylabel('Gravel Volumetric Fraction', fontsize=12)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels, title='Regions', fontsize=10, loc='upper left', bbox_to_anchor=(1.02, 1))

plt.tight_layout()
box_save = f"{out_dir}/gravel_5regions_dynamic_boxplot.pdf"
plt.savefig(box_save, dpi=300, bbox_inches='tight')
plt.close()

# ====================== 8. 绘制 折线图 ======================
print("正在绘制 5 区域折线图...")
fig, ax = plt.subplots(figsize=(12, 7))
sns.lineplot(
    data=df_all, x='Layer', y='Gravel Fraction', hue='Region', 
    palette=custom_palette, marker='o', markersize=8, linewidth=2.5, 
    errorbar='sd', ax=ax
)

ax.axhline(0.3, color='gray', linestyle='--', linewidth=1.5, zorder=0, label='Threshold (0.3)')

ax.set_title('Mean Gravel Volumetric Fraction Trend Across Layers', fontsize=15, fontweight='bold', pad=15)
ax.set_xlabel('Soil Layer', fontsize=12)
ax.set_ylabel('Mean Gravel Volumetric Fraction', fontsize=12)
ax.set_xticks(range(1, 9))

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels, title='Regions', fontsize=10, loc='upper left', bbox_to_anchor=(1.02, 1))

plt.tight_layout()
line_save = f"{out_dir}/gravel_5regions_dynamic_lineplot.pdf"
plt.savefig(line_save, dpi=300, bbox_inches='tight')
plt.close()

# ====================== 9. 绘制 均值统计表格图 ======================
print("正在绘制 均值统计三线表...")
df_means_rounded = df_means.round(4).astype(str)

# 重置索引，将原先的行名变成真实的第一列，并赋予英文表头
df_means_rounded = df_means_rounded.reset_index()
df_means_rounded.rename(columns={'Region': 'Region Classification (Threshold: 0.3)'}, inplace=True)

# 调整画布比例，稍微缩减高度以拉近标题和表格的距离
fig, ax = plt.subplots(figsize=(16, 2.5 + len(df_means)*0.3))
ax.axis('tight')
ax.axis('off')

# 渲染表格内容
table = ax.table(
    cellText=df_means_rounded.values,
    colLabels=df_means_rounded.columns,
    cellLoc='center', loc='center'
)

# 【核心修改 1】：自动调整列宽，彻底解决第一列长文本溢出重叠的问题
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1.0, 1.8) # X轴设为1.0，让 auto_set_column_width 接管宽度
table.auto_set_column_width(col=list(range(len(df_means_rounded.columns))))

# 核心逻辑：移除所有默认边框，只画上中下三条线
for (row, col), cell in table.get_celld().items():
    # 1. 初始化：清除所有边框，把背景设为纯白
    cell.set_edgecolor('black')
    cell.set_facecolor('white')
    cell.visible_edges = ''
    
    # 2. 绘制三线表的“第一条”和“第二条”线 (即表头的顶部和底部)
    if row == 0:
        cell.visible_edges = 'BT'  
        cell.set_linewidth(1.5)
        
    # 3. 绘制三线表的“第三条”线 (即表格最后一行数据的底部)
    elif row == len(df_means_rounded):
        cell.visible_edges = 'B'   
        cell.set_linewidth(1.5)
        
    # 4. 表头和最左侧行名加粗
    if row == 0 or col == 0:
        cell.set_text_props(weight='bold')
        
    # 5. 第一列靠左对齐更美观，数字列保持居中
    if col == 0:
        # 必须同时修改单元格锚点和底层文本对齐属性
        cell._loc = 'left'
        cell.set_text_props(ha='left')
        
        # 增加左侧空格缩进（防止文字死死贴着左边边缘）
        old_text = cell.get_text().get_text()
        # 避免多次循环重复加空格，判断一下
        if not old_text.startswith('   '):
            cell.get_text().set_text(f'   {old_text}')

# 【核心修改 2 & 3】：更新翻译后的英文标题，并将 pad 设为 5 拉近距离
title_str = "Regional Mean of Gravel Volumetric Fraction Across Different Regions and Soil Layers"
plt.title(title_str, fontsize=14, fontweight='bold')

table_save = f"{out_dir}/gravel_5regions_dynamic_mean_table.pdf"
plt.savefig(table_save, dpi=300, bbox_inches='tight')
plt.close()

print("统计图全部处理完毕！")