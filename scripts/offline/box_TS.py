import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy import stats
from scipy.stats import gaussian_kde

# -------------------------- 核心参数配置 --------------------------
csv_file = 'gravel_TS_deep_TCC.csv' # CSV数据文件路径
x_col = 'Gravel content' # CSV中自变量X的列名
y_col = 'Soil Temperature TCC (Deep)' # CSV中因变量Y的列名
group_col = 'Group' # CSV中分组列的列名
x_label = 'Gravel content' # 图中X轴的标签
y_label = 'Soil Temperature TCC (Deep)' # 图中Y轴的标签
output_filename = 'regression_TS_deep_TCC_kde.png' # 输出图片的文件名
figsize = (7, 7)  # 图片尺寸（宽, 高），单位英寸
confidence_level = 0.95 # 回归置信区间（默认95%）
version = 2 # 1=箱线图边缘版；2=密度图边缘版
colors = [ '#B2A3DD', "#FFD900", "#96CDEAB7", '#ED949A']  # 定义分组颜色(3种颜色,超过3组会循环使用,可自行添加)
fontcols = ["#9985D1","#FFD900","#50BBF5B7",  "#EC4853EE"]

# 1. 读取数据
df = pd.read_csv(csv_file)
# 若未指定轴标签，用CSV列名作为标签
x_label = x_label or x_col
y_label = y_label or y_col
# 获取所有分组
groups = df[group_col].unique()
# 存储回归结果
results = []

# 2. 创建图形和子图（3个子图：主图+顶部边缘图+右侧边缘图）
fig = plt.figure(figsize=figsize)
# 子图位置参数（left=左距, bottom=下距, width=宽, height=高, spacing=间距）
left, width, bottom, height, spacing = 0.15, 0.60, 0.15, 0.60, 0.005
# 主图：放散点+回归线+置信区间
ax_scatter = plt.axes([left, bottom, width, height])
# 顶部边缘图：放X轴分布（箱线图/密度图）
ax_histx = plt.axes([left, bottom + height + spacing, width, 0.15])
# 右侧边缘图：放Y轴分布（箱线图/密度图）
ax_histy = plt.axes([left + width + spacing, bottom, 0.15, height])
# 预先设置主图坐标轴范围（加8% padding，避免数据贴边，边缘图也能对齐）
x_padding = (df[x_col].max() - df[x_col].min()) * 0.08
y_padding = (df[y_col].max() - df[y_col].min()) * 0.08
ax_scatter.set_xlim(df[x_col].min() - x_padding, df[x_col].max() + x_padding)
ax_scatter.set_ylim(df[y_col].min() - y_padding, df[y_col].max() + y_padding)
# 存储分组数据（后续画边缘图用，避免重复筛选）
group_data_list = []

# 3. 分组处理：散点图+回归分析+边缘图
for i, group in enumerate(groups):   
    # 提取当前组的数据    
    group_data = df[df[group_col] == group]   
    # 当前组的颜色    
    color = colors[i % len(colors)]   
    fontcol = fontcols[i % len(colors)]
    # 存储分组数据    
    group_data_list.append((group, group_data, color))
    # 3.1 主图：绘制分组散点    
    ax_scatter.scatter(
            group_data[x_col], group_data[y_col], 
            color=color, s=10, alpha=0.4,  # 大小40，透明度0.4        
            edgecolor='black', linewidth=0,  # 黑色描边，突出散点    
            )
    # 3.2 线性回归分析(样本量>1才做回归)
    if len(group_data) > 1:       
       # 准备自变量X（需为2D数组）和因变量Y        
       X = group_data[[x_col]].values        
       y = group_data[y_col].values       
       # 拟合线性回归模型        
       model = LinearRegression().fit(X, y)        
       y_pred = model.predict(X)
       # 计算回归统计量        
       r2 = r2_score(y, y_pred)  # R²值（拟合优度）        
       slope = model.coef_[0]  # 回归斜率        
       intercept = model.intercept_  # 回归截距        
       n = len(X)  # 样本量       
       # 存储回归结果        
       results.append({           
           'group': group, 'r2': r2, 'slope': slope,           
           'intercept': intercept, 'n': n, 'color': color        
           })
       # 绘制回归线和95%置信区间       
       # 生成回归线的X取值（从X最小值到最大值，100个点，使线平滑）        
       x_vals = np.linspace(X.min(), X.max(), 100)       
       # 回归线的Y值        
       y_vals = model.predict(x_vals.reshape(-1, 1))       
       # 计算95%置信区间（基于t分布）        
       mse = np.sum((y - y_pred) ** 2) / (n - 2)  # 均方误差       
       # 标准误差（置信区间宽度的核心）        
       se = np.sqrt(mse * (1 / n + (x_vals - X.mean()) ** 2 / np.sum((X - X.mean()) ** 2))) # t分布临界值（自由度n-2，95%置信水平）        
       t_val = stats.t.ppf((1 + confidence_level) / 2, n - 2)       
       # 绘制置信区间（半透明填充）        
       ax_scatter.fill_between(
            x_vals, y_vals - t_val * se, y_vals + t_val * se,
            color=fontcol, alpha=0.2, linewidth=0 # 透明度0.2，无轮廓线
            )
       # 绘制回归线（虚线，宽度2.5）      
       ax_scatter.plot(            
            x_vals, y_vals, color=fontcol, linestyle='--',             
            linewidth=2.0, alpha=1.0,    
            )
    else:
       print(f"警告：组 {group} 样本数为 {len(group_data)}，跳过回归分析。")
    # 3.3 版本1：边缘图-箱线图   
    if version == 1:       
       # 顶部边缘图（X轴分布，水平箱线图）        
       bp_x = ax_histx.boxplot(            
            [group_data[x_col].values], positions=[i + 1],  # 位置对应分组            
             vert=False, patch_artist=True,  # 水平箱线图，填充颜色            
             widths=0.6, showfliers=False,  # 宽度0.6，隐藏异常点           
             # 箱线图样式（颜色与主图一致）            
             boxprops=dict(linewidth=2, color=color),            
             medianprops=dict(linewidth=2, color=color),            
             whiskerprops=dict(linewidth=1.5, color=color),            
             capprops=dict(linewidth=1.5, color=color)        
            )
       # 填充箱线图颜色        
       bp_x['boxes'][0].set_facecolor(color)        
       bp_x['boxes'][0].set_alpha(0.7)       
       # 右侧边缘图（Y轴分布，垂直箱线图）        
       bp_y = ax_histy.boxplot(            
            [group_data[y_col].values], positions=[i + 1],            
             vert=True, patch_artist=True,            
             widths=0.6, showfliers=False,            
             boxprops=dict(linewidth=2, color=color),            
             medianprops=dict(linewidth=2, color=color),            
             whiskerprops=dict(linewidth=1.5, color=color),            
             capprops=dict(linewidth=1.5, color=color)        
            )       
       # 填充箱线图颜色        
       bp_y['boxes'][0].set_facecolor(color)        
       bp_y['boxes'][0].set_alpha(0.7)       
       # 调整箱线图的坐标轴范围（避免分组重叠）        
       ax_histx.set_ylim(0.5, len(groups) + 0.5)        
       ax_histy.set_xlim(0.5, len(groups) + 0.5)
    # 3.4 版本2：边缘图-密度图（在循环外绘制，确保使用完整的坐标轴范围）
    if version == 2:   
       for group, group_data, color in group_data_list:       
           # 计算X和Y的核密度估计        
           x_kde = gaussian_kde(group_data[x_col])        
           y_kde = gaussian_kde(group_data[y_col])       
           # 使用主图的完整范围来生成密度图        
           x_plot = np.linspace(ax_scatter.get_xlim()[0], ax_scatter.get_xlim()[1], 200)        
           y_plot = np.linspace(ax_scatter.get_ylim()[0], ax_scatter.get_ylim()[1], 200)       
           # 顶部边缘图（X轴密度图）        
           # ax_histx.plot(x_plot, x_kde(x_plot), color=color, linewidth=2)        
           # ax_histx.fill_between(x_plot, x_kde(x_plot), alpha=0.3, color=color)       
           # 右侧边缘图（Y轴密度图）        
           ax_histy.plot(y_kde(y_plot), y_plot, color=color, linewidth=2)        
           ax_histy.fill_betweenx(y_plot, y_kde(y_plot), alpha=0.3, color=color)
# 4. 图的美化和标注
# 4.1 主图：添加R²统计信息（左上角，分组排列）
x_pos = 0.6
y_pos = 0.2 # 起始位置
for result in results: 
    ax_scatter.text(
            x_pos, y_pos,
            f"{result['group']}: $R^2$ ={result['r2']:.3f}",
            color=result['color'], transform=ax_scatter.transAxes, 
            fontsize=10, fontweight='bold'
            )
    y_pos -= 0.06 # 每组文本向下移
# 4.2 主图：坐标轴标签
ax_scatter.set_xlabel(x_label, fontsize=16, fontweight='bold')  # X轴标签（16号字）
ax_scatter.set_ylabel(y_label, fontsize=16, fontweight='bold')  # Y轴标签（16号字）
# 4.3 边缘图：隐藏边框
for ax in [ax_histx, ax_histy]:
    for spine in ax.spines.values():
        spine.set_visible(False)    
    ax.set_xticks([])    
    ax.set_yticks([])
# 4.4 同步主图和边缘图的坐标轴范围
ax_histx.set_xlim(ax_scatter.get_xlim())
ax_histy.set_ylim(ax_scatter.get_ylim())
# 4.5保存为高分辨率图片
plt.savefig(output_filename, dpi=600)
plt.show()
