import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import geopandas as gpd
from matplotlib.path import Path
import warnings
warnings.filterwarnings('ignore')

# ====================== 字体全局设置 ======================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ==========================================
# 1. 物理常数、深度设定与层级映射
# ==========================================
tfrz = 273.16  
cpliq = 4188.0     
cpice = 2117.27    
dz = np.array([0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058])
z  = np.array([0.0071, 0.0279, 0.0623, 0.1189, 0.2122, 0.3661, 0.6198, 1.0380, 1.7276, 2.8646])
zh = np.array([0.0175, 0.0451, 0.0906, 0.1655, 0.2891, 0.4929, 0.8289, 1.3828, 2.2961, 3.8019])
layer_map = [0, 0, 1, 2, 3, 4, 5, 6, 7, 7] 

# ==========================================
# 2. 生成青藏高原 1D 掩膜索引
# ==========================================
print("初始化青藏高原掩膜...")
base_path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
wrfinput_file = "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
shp_path = "../../shapefile_China/TPBoundary_new_2021/TPBoundary_new(2021).shp"

ds_wrf = xr.open_dataset(wrfinput_file)
lat2d = ds_wrf['XLAT'].isel(Time=0).values
lon2d = ds_wrf['XLONG'].isel(Time=0).values
ds_wrf.close()

qtp_gdf = gpd.read_file(shp_path)
poly = qtp_gdf.geometry.union_all()
polygons = list(poly.geoms) if poly.geom_type == 'MultiPolygon' else [poly]
points = np.column_stack((lon2d.ravel(), lat2d.ravel()))
qtp_mask_1d = np.zeros(lat2d.size, dtype=bool)
for p in polygons:
    path = Path(np.asarray(p.exterior.coords))
    qtp_mask_1d |= path.contains_points(points, radius=0)
mask_2d = qtp_mask_1d.reshape(lat2d.shape)

# ==========================================
# 3. 数据提取与物理张量运算
# ==========================================
def get_fast_seasonal_means(ds_djf, ds_mj, var_name):
    da_djf = ds_djf[var_name].isel(soilsnow=slice(5, 15))
    da_djf_np = da_djf.transpose('time', 'soilsnow', 'lat', 'lon').values
    da_mj = ds_mj[var_name].isel(time=slice(0, 102), soilsnow=slice(5, 15))
    da_mj_np = da_mj.transpose('time', 'soilsnow', 'lat', 'lon').values
    
    qtp_djf = np.where(da_djf_np[:, :, mask_2d] < -1000, np.nan, da_djf_np[:, :, mask_2d])
    qtp_mj  = np.where(da_mj_np[:, :, mask_2d] < -1000, np.nan, da_mj_np[:, :, mask_2d])

    win = np.nanmean(qtp_djf, axis=0) 
    qtp_mj_reshaped = qtp_mj.reshape(17, 6, 10, -1)
    spr = np.nanmean(qtp_mj_reshaped[:, 0:3, :, :], axis=(0, 1))
    sum_ = np.nanmean(qtp_mj_reshaped[:, 3:6, :, :], axis=(0, 1))
    return win, spr, sum_

def extract_param(f3, var_name, is_8_layer=True):
    da_qtp = f3[var_name].values[:, mask_2d]
    return da_qtp[layer_map, :] if is_8_layer else da_qtp

print("提取数据并运算...")
f1 = xr.open_dataset(os.path.join(base_path, "colmoff_2001-2023_monmean_nogravel.nc"))
f2 = xr.open_dataset(os.path.join(base_path, "colmoff_2001-2023_monmean_gravel.nc"))
f1_djf = xr.open_dataset(os.path.join(base_path, "colmoff_2001-2017_DJF_nogravel.nc"))
f2_djf = xr.open_dataset(os.path.join(base_path, "colmoff_2001-2017_DJF_gravel.nc"))
f3 = xr.open_dataset("interp_soil_params.nc")

tc_w, tc_sp, tc_su = get_fast_seasonal_means(f1_djf, f1, 'f_t_soisno')
te_w, te_sp, te_su = get_fast_seasonal_means(f2_djf, f2, 'f_t_soisno')
wc_w, wc_sp, wc_su = get_fast_seasonal_means(f1_djf, f1, 'f_wliq_soisno')
we_w, we_sp, we_su = get_fast_seasonal_means(f2_djf, f2, 'f_wliq_soisno')
ic_w, ic_sp, ic_su = get_fast_seasonal_means(f1_djf, f1, 'f_wice_soisno')
ie_w, ie_sp, ie_su = get_fast_seasonal_means(f2_djf, f2, 'f_wice_soisno')

def calc_Sr_vec(th, ths):
    valid_mask = (ths != 0) & (~np.isnan(ths))
    Sr = np.zeros_like(th)
    Sr[valid_mask] = th[valid_mask] / ths[valid_mask]
    return np.clip(Sr, 1e-4, 1.0)

def calc_k_vec(t, Sr, tkd, tksu, tksf, vom, vs, vgr, alf, bet):
    ksat = np.where(t > tfrz, tksu, tksf)
    base = np.maximum((1.0 / (1.0 + np.exp(-bet * Sr)))**3 - ((1.0 - Sr)/2.0)**3, 1e-6)
    power = 0.5 * (1.0 + vom - alf*vs - vgr)
    ke_unf = (Sr ** power) * (base ** (1.0 - vom))
    ke_frz = Sr ** (1.0 + vom)
    ke = np.where(t > tfrz, ke_unf, ke_frz)
    return (ksat - tkd) * ke + tkd

tkd_c = extract_param(f3, 'tkdry_ctl_interp'); tksu_c = extract_param(f3, 'tksatu_ctl_interp'); tksf_c = extract_param(f3, 'tksatf_ctl_interp'); vom_c = extract_param(f3, 'vf_om_s_ctl_interp'); vs_c = extract_param(f3, 'vf_sand_s_ctl_interp'); vgr_c = extract_param(f3, 'vf_gravels_s_ctl_interp'); alf_c = extract_param(f3, 'BA_alpha_ctl_interp'); bet_c = extract_param(f3, 'BA_beta_ctl_interp'); csol_c = extract_param(f3, 'csol_ctl_interp')
tkd_e = extract_param(f3, 'tkdry_exp_interp'); tksu_e = extract_param(f3, 'tksatu_exp_interp'); tksf_e = extract_param(f3, 'tksatf_exp_interp'); vom_e = extract_param(f3, 'vf_om_s_exp_interp'); vs_e = extract_param(f3, 'vf_sand_s_exp_interp'); vgr_e = extract_param(f3, 'vf_gravels_s_exp_interp'); alf_e = extract_param(f3, 'BA_alpha_exp_interp'); bet_e = extract_param(f3, 'BA_beta_exp_interp'); csol_e = extract_param(f3, 'csol_exp_interp')

K_c_w = calc_k_vec(tc_w, calc_Sr_vec(extract_param(f3, 'theta_ctl_win_re', False), extract_param(f3, 'theta_s_ctl_interp')), tkd_c, tksu_c, tksf_c, vom_c, vs_c, vgr_c, alf_c, bet_c)
K_c_sp = calc_k_vec(tc_sp, calc_Sr_vec(extract_param(f3, 'theta_ctl_spr_re', False), extract_param(f3, 'theta_s_ctl_interp')), tkd_c, tksu_c, tksf_c, vom_c, vs_c, vgr_c, alf_c, bet_c)
K_c_su = calc_k_vec(tc_su, calc_Sr_vec(extract_param(f3, 'theta_ctl_sum_re', False), extract_param(f3, 'theta_s_ctl_interp')), tkd_c, tksu_c, tksf_c, vom_c, vs_c, vgr_c, alf_c, bet_c)

K_e_w = calc_k_vec(te_w, calc_Sr_vec(extract_param(f3, 'theta_exp_win_re', False), extract_param(f3, 'theta_s_exp_interp')), tkd_e, tksu_e, tksf_e, vom_e, vs_e, vgr_e, alf_e, bet_e)
K_e_sp = calc_k_vec(te_sp, calc_Sr_vec(extract_param(f3, 'theta_exp_spr_re', False), extract_param(f3, 'theta_s_exp_interp')), tkd_e, tksu_e, tksf_e, vom_e, vs_e, vgr_e, alf_e, bet_e)
K_e_su = calc_k_vec(te_su, calc_Sr_vec(extract_param(f3, 'theta_exp_sum_re', False), extract_param(f3, 'theta_s_exp_interp')), tkd_e, tksu_e, tksf_e, vom_e, vs_e, vgr_e, alf_e, bet_e)

dz_mat = dz[:, np.newaxis]
C_c_w = csol_c + (np.nan_to_num(wc_w)/dz_mat)*cpliq + (np.nan_to_num(ic_w)/dz_mat)*cpice
C_c_sp = csol_c + (np.nan_to_num(wc_sp)/dz_mat)*cpliq + (np.nan_to_num(ic_sp)/dz_mat)*cpice
C_c_su = csol_c + (np.nan_to_num(wc_su)/dz_mat)*cpliq + (np.nan_to_num(ic_su)/dz_mat)*cpice

C_e_w = csol_e + (np.nan_to_num(we_w)/dz_mat)*cpliq + (np.nan_to_num(ie_w)/dz_mat)*cpice
C_e_sp = csol_e + (np.nan_to_num(we_sp)/dz_mat)*cpliq + (np.nan_to_num(ie_sp)/dz_mat)*cpice
C_e_su = csol_e + (np.nan_to_num(we_su)/dz_mat)*cpliq + (np.nan_to_num(ie_su)/dz_mat)*cpice

Alpha_c_w, Alpha_c_sp, Alpha_c_su = K_c_w/C_c_w, K_c_sp/C_c_sp, K_c_su/C_c_su
Alpha_e_w, Alpha_e_sp, Alpha_e_su = K_e_w/C_e_w, K_e_sp/C_e_sp, K_e_su/C_e_su

def mean1D(arr): 
    return np.nanmean(arr, axis=1)

# ==============================================================================
# 4. 客制化控制台
# ==============================================================================

# 【配置1】: 选择你想展示的土壤层数，范围是 1 到 10
CUSTOM_LAYERS = [2, 5, 8]  

# 【配置2】: 选择你想绘制的图表组合 ("water_ice", "thermal", 或 "both")
PLOT_CHOICE = "both"  

# ==============================================================================

# 真实层数转换为代码索引
target_idx = [layer - 1 for layer in CUSTOM_LAYERS]

def get_target_data(var_w, var_sp, var_su):
    full_data = np.vstack([mean1D(var_w), mean1D(var_sp), mean1D(var_su)]).T
    return full_data[target_idx, :]

data_dict = {
    "Liquid Water Content": (get_target_data(wc_w, wc_sp, wc_su), get_target_data(we_w, we_sp, we_su)),
    "Ice Content":          (get_target_data(ic_w, ic_sp, ic_su), get_target_data(ie_w, ie_sp, ie_su)),
    "Heat Capacity":        (get_target_data(C_c_w, C_c_sp, C_c_su) * 1e-6, get_target_data(C_e_w, C_e_sp, C_e_su) * 1e-6),
    "Thermal Conductivity": (get_target_data(K_c_w, K_c_sp, K_c_su), get_target_data(K_e_w, K_e_sp, K_e_su)),
    "Thermal Diffusivity":  (get_target_data(Alpha_c_w, Alpha_c_sp, Alpha_c_su) * 1e6, get_target_data(Alpha_e_w, Alpha_e_sp, Alpha_e_su) * 1e6)
}

# ==========================================
# 5. 可视化绘图主函数
# ==========================================
def plot_custom_grid(vars_info, filename_suffix, figsize):
    print(f"正在生成 {filename_suffix}...")
    n_cols = len(vars_info)
    fig, axes = plt.subplots(2, n_cols, figsize=figsize)
    # 增加 hspace 以容纳第二行子图的外置标题
    plt.subplots_adjust(hspace=0.35, wspace=0.25, top=0.92) 
    
    if n_cols == 1:
        axes = np.atleast_2d(axes).T
        
    x_vals = np.array([0, 1, 2])
    x_labels = ['DJF', 'MAM', 'JJA']
    
    # 动态调色板和标记
    color_palette = plt.cm.tab10(np.linspace(0, 1, 10))
    marker_palette = ['o', '^', 's', 'D', 'v', 'p', '*', 'h', 'X', '<']
    
    for col_idx, v_info in enumerate(vars_info):
        var_name = v_info["name"]
        unit_str = v_info["unit"]
        
        ax_raw = axes[0, col_idx]
        ax_diff = axes[1, col_idx]
        
        ctl_data, exp_data = data_dict[var_name]
        diff_data = exp_data - ctl_data
        
        for i, original_layer_idx in enumerate(target_idx):
            color = color_palette[i % len(color_palette)]
            marker = marker_palette[i % len(marker_palette)]
            
            # 画线
            ax_raw.plot(x_vals, ctl_data[i], color=color, ls='-', marker=marker, ms=7, lw=1.5)
            ax_raw.plot(x_vals, exp_data[i], color=color, ls='--', marker=marker, ms=7, lw=1.5, markerfacecolor='white')
            
            ax_diff.plot(x_vals, diff_data[i], color=color, ls='-.', marker=marker, ms=7, lw=1.5)
            
        # ==========================================
        # 设置外部标题与单位 (替代原先的 ylabel 和中央 title)
        # ==========================================
        label_raw = chr(97 + col_idx)            # (a), (b), (c)...
        label_diff = chr(97 + col_idx + n_cols)  # (d), (e), (f)...
        
        # 第一行子图: 左上角放 "(序号) 变量名", 右上角放 "单位"
        ax_raw.set_title(f"({label_raw}) {var_name}", loc='left', fontsize=12, fontweight='normal')
        ax_raw.set_title(unit_str, loc='right', fontsize=11, fontweight='normal')
        
        # 第二行子图(差值): 左上角放 "(序号) Δ 变量名", 右上角不放单位 (按要求只在第一行显示)
        ax_diff.set_title(f"({label_diff}) $\Delta$ {var_name}", loc='left', fontsize=12, fontweight='normal')
        
        
        # ==========================================
        # 仅在第一列 (col_idx == 0) 生成原格式图例
        # ==========================================
        if col_idx == 0:
            # 1. 状态图 (CTL/EXP) 分类说明
            style_raw_handles = [
                Line2D([0], [0], color='dimgray', linestyle='-', lw=1.5, label='CTL OFF'),
                Line2D([0], [0], color='dimgray', linestyle='--', lw=1.5, label='EXP OFF')
            ]
            
            # 2. 差值图说明
            style_diff_handles = [
                Line2D([0], [0], color='dimgray', linestyle='-.', lw=1.5, label='EXP - CTL')
            ]
            
            # 3. 共用的土壤层数说明
            layer_handles = [
                Line2D([0], [0], color=color_palette[i % len(color_palette)], 
                       marker=marker_palette[i % len(marker_palette)], lw=1.5, ms=6, 
                       label=f'L{CUSTOM_LAYERS[i]} ({z[target_idx[i]]:.3f}m)')
                for i in range(len(target_idx))
            ]
            
            # 强制图例为单列排版 (ncol=1)
            ax_raw.legend(handles=style_raw_handles + layer_handles, loc='best', ncol=1, 
                          frameon=True, facecolor='white', framealpha=0.85, edgecolor='lightgray', fontsize=8.5)
                          
            ax_diff.legend(handles=style_diff_handles + layer_handles, loc='best', ncol=1, 
                           frameon=True, facecolor='white', framealpha=0.85, edgecolor='lightgray', fontsize=8.5)
        # ==========================================

        for ax in [ax_raw, ax_diff]:
            ax.set_xticks(x_vals)
            ax.set_xticklabels(x_labels, fontsize=11)
            ax.grid(axis='y', linestyle='--', alpha=0.4)
            ax.tick_params(axis='both', labelsize=10, direction='in')
            
            if ax == ax_diff:
                ax.axhline(0, color='black', linestyle=':', lw=1.5, alpha=0.5, zorder=0)
                y_min, y_max = ax.get_ylim()
                bound = max(abs(y_min), abs(y_max))
                ax.set_ylim(-bound, bound)

    output_dir = '../illustration/off/'
    os.makedirs(output_dir, exist_ok=True)
    output_fig = os.path.join(output_dir, f'{filename_suffix}.pdf')
    plt.savefig(output_fig, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"成功保存至: {output_fig}")

# ==========================================
# 6. 执行绘图任务
# ==========================================
# 注意这里重新分离了 name 和 unit 以适配新排版
vars_water_ice = [
    {"name": "Liquid Water Content", "unit": "kg m$^{-2}$"},
    {"name": "Ice Content",          "unit": "kg m$^{-2}$"}
]

vars_thermal = [
    {"name": "Heat Capacity",        "unit": "10$^6$ J m$^{-3}$K$^{-1}$"},
    {"name": "Thermal Conductivity", "unit": "W m$^{-1}$K$^{-1}$"},
    {"name": "Thermal Diffusivity",  "unit": "10$^{-6}$ m$^{2}$s$^{-1}$"}
]

if PLOT_CHOICE in ["water_ice", "both"]:
    plot_custom_grid(vars_water_ice, "QTP_Water_Ice", figsize=(10, 8))

if PLOT_CHOICE in ["thermal", "both"]:
    plot_custom_grid(vars_thermal, "QTP_Thermal", figsize=(15, 8))

print("任务全部执行完毕！")