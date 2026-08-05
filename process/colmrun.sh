#!/bin/bash

# 定义输入和输出基础路径
base_input_path="/share/home/dq135/CRESM_0307S92/CASES/RUN_"
base_output_path="/share/home/dq135/draw/process/temporary/"

out_list = ''
# 开始年份到结束年份的循环
for year in {2001..2017}; do  # 这里你可以更改年份范围
    # 处理6月到8月
    for mon in {3..8}; do
        # 格式化月份和天数为两位数
        mon_str=$(printf "%02d" $mon)
        # 拼接输入和输出文件路径
        input_name="${base_input_path}${year}_nogravel/CoLMrun/unstructured_cwrf_CN_15km_nogravel/history/unstructured_cwrf_CN_15km_nogravel_hist_${year}-${mon_str}_remap.nc"
        output_name="${base_output_path}colmrun_${year}-${mon_str}_nogravel.nc"

        # cdo命令提取变量
        # cdo selname,time,f_xy_t,f_xy_q,f_xy_prc,f_xy_prl,f_xy_pbot,f_xy_frl,f_xy_solarin,f_sr,f_fsena,f_lfevpa,f_fevpa,f_fseng,f_fevpg,f_fgrnd,f_olrg,f_rnet,f_t_grnd,f_scv,f_fsno,f_sigf,f_tref,f_qref,f_rss,f_t_soisno,f_wliq_soisno,f_wice_soisno,f_h2osoi,f_us10m,f_vs10m "$input_name" "$output_name"
        
	cdo selname,time,f_ustar,f_zol,f_fh,f_z0m,f_fsena,f_fseng,f_fsenl "$input_name" "$output_name"

        # 将输出文件名添加到列表
        out_list="$out_list $output_name"
    done
done
# 打印所有生成的输出文件
# echo "Files to merge: $out_list"
cdo mergetime $out_list add_colmrun_2001-2017_monmean_nogravel.nc
rm $out_list
