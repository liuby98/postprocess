#!/bin/bash
monday=(31 28 31 30 31 30 31 31 30 31 30 31)
base_input_path="/share/home/dq135/CRESM_0307S92/SCRIPTS/ICBCScripts/PrepCoLM/CN_15km/CoLMrun_gravel/unstructured_cwrf_CN_15km_gravel/history/unstructured_cwrf_CN_15km_gravel_hist_"
base_output_path="/share/home/dq135/draw/process/temporary/"

out_list=''
# 年份修改为 2001 到 2017
for year in {2001..2017}; do
    # 计算上一年的年份
    prev_year=$((year - 1))
    
    # 直接遍历前年12月、今年01月、今年02月
    for ym in "${year}-03" "${year}-04" "${year}-05" "${year}-06" "${year}-07" "${year}-08"; do
        input_name="${base_input_path}${ym}_remap.nc"
        output_name="${base_output_path}colmoff_${ym}_gravel.nc"
        
        #cdo selvar,time,f_xy_t,f_xy_q,f_xy_prc,f_xy_prl,f_xy_pbot,f_xy_frl,f_xy_solarin,f_sr,f_fsena,f_lfevpa,f_fevpa,f_fseng,f_fevpg,f_fgrnd,f_olrg,f_rnet,f_t_grnd,f_scv,f_fsno,f_sigf,f_tref,f_qref,f_rss,f_t_soisno,f_wliq_soisno,f_wice_soisno,f_h2osoi,f_us10m,f_vs10m "$input_name" "$output_name"
        cdo selvar,time,f_xy_us,f_xy_vs,f_ustar,f_ustar2 "$input_name" "$output_name"
        out_list="$out_list $output_name"
    done
done

# 将合并后的输出文件名更改为符合当前时间范围和季节(DJF)的名称
cdo mergetime $out_list colmoff_2001-2017_uv_gravel.nc

# 删除临时文件
rm $out_list
