import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# -------------------------- 核心参数配置 --------------------------
csv_file = '/share/home/dq135/draw/scripts/gravel_SMrz_TCC.csv'  # CSV数据文件路径
x_col = 'Gravel content'  # CSV中自变量X的列名
y_col = 'Soil Moisture Root Zone TCC'  # CSV中因变量Y的列名
group_col = 'Group'  # CSV中分组列的列名
x_label = 'Gravel content'  # 图中X轴的标签
y_label = 'Soil Moisture Root Zone TCC'  # 图中Y轴的标签
output_filename = 'curve_SMrz_TCC_gravel.png'  # 输出图片的文件名
figsize = (10, 7)  # 图片尺寸（宽, 高），单位英寸（加宽以适应更多区间）
colors = ['#B2A3DD', "#FFD900", "#96CDEAFF", "#F87079"]  # 定义分组颜色
smoothing_factor = 0.5  # 平滑因子（UnivariateSpline的s参数，越大越平滑；0为无平滑）

# 1. 读取数据
df = pd.read_csv(csv_file)
# 若未指定轴标签，用CSV列名作为标签
x_label = x_label or x_col
y_label = y_label or y_col
# 获取所有分组（假设4个）
groups = df[group_col].unique()
if len(groups) != 4:
    print(f"警告：CSV文件中分组数为 {len(groups)}，预期为4。")

# 2. 定义分桶区间和中心点（缩小到0.01步长）
bins = np.arange(0, 1.01, 0.01)  # 区间边界：0, 0.01, 0.02, ..., 1.00
bin_centers = np.arange(0.005, 1.005, 0.01)  # 每个区间的中心点：0.005, 0.015, ..., 0.995

# 3. 创建图形
fig, ax = plt.subplots(figsize=figsize)

# 4. 分组处理：计算平均值、平滑并绘制曲线
for i, group in enumerate(groups):
    # 提取当前组的数据（添加.copy()以避免警告）
    group_data = df[df[group_col] == group].copy()
    # 添加分桶列（使用pd.cut，labels为bin_centers，便于groupby）
    group_data['bin'] = pd.cut(group_data[x_col], bins=bins, labels=bin_centers, include_lowest=True)
    # 计算每个区间的TCC平均值（设置observed=True以避免警告）
    mean_tcc = group_data.groupby('bin', observed=True)[y_col].mean().reindex(bin_centers)  # 确保所有bin_centers都有值（无数据为NaN）
    # 当前组的颜色
    color = colors[i % len(colors)]
    
    # 绘制原折线（可选，半透明以突出平滑曲线）
    # ax.plot(bin_centers, mean_tcc, color=color, linewidth=1.0, linestyle='--', alpha=0.8, label=f"{group} (original)") #
    
    # 平滑处理：去除NaN，只对有数据的点进行样条插值
    valid_mask = ~mean_tcc.isna()
    if valid_mask.sum() > 1:  # 至少2个点才能插值
        x_valid = bin_centers[valid_mask]
        y_valid = mean_tcc[valid_mask]
        # 使用UnivariateSpline进行平滑
        spl = UnivariateSpline(x_valid, y_valid, s=smoothing_factor)
        # 生成平滑曲线（使用更多点以使曲线光滑）
        x_smooth = np.linspace(x_valid.min(), x_valid.max(), 500)
        y_smooth = spl(x_smooth)
        # 绘制平滑曲线
        ax.plot(x_smooth, y_smooth, color=color, linewidth=2.0, label=f"{group}")
    else:
        print(f"警告：组 {group} 有效数据点不足，无法平滑。")
    
    # 添加数据点（可选，避免曲线太光滑）
    # ax.scatter(bin_centers[valid_mask], mean_tcc[valid_mask], color=color, s=10)

# 5. 图的美化和标注
# 添加图例（显示每个分组）
ax.legend(title='Groups', loc='upper right', fontsize=10)
# 坐标轴标签
ax.set_xlabel(x_label, fontsize=16, fontweight='bold')
ax.set_ylabel(y_label, fontsize=16, fontweight='bold')
# 设置x轴范围和刻度（覆盖0到1，每0.1一个主刻度）
ax.set_xlim(0, 0.9)
ax.set_xticks(np.arange(0, 1.0, 0.1))
# 添加网格（可选，便于阅读）
ax.grid(True, linestyle='--', alpha=0.5)

# 6. 保存为高分辨率图片并显示
plt.savefig(output_filename, dpi=600)
plt.show()