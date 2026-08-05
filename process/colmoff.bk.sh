#!/bin/bash
monday=(31 28 31 30 31 30 31 31 30 31 30 31)
base_input_path="/share/home/dq135/CRESM_0307S92/SCRIPTS/ICBCScripts/PrepCoLM/CN_15km/CoLMrun_gravel/unstructured_cwrf_CN_15km_gravel/history/unstructured_cwrf_CN_15km_gravel_hist_"
base_output_path="/share/home/dq135/draw/process/temporary/"

out_list=''
for year in {2001..2023}; do
    for mon in {3..8}; do
	mon_str=$(printf "%02d" $mon)
        input_name="${base_input_path}${year}-${mon_str}_remap.nc"
        output_name="${base_output_path}colmoff_${year}-${mon_str}_gravel.nc"
	cdo selvar,time,f_xy_t,f_xy_q,f_xy_prc,f_xy_prl,f_xy_pbot,f_xy_frl,f_xy_solarin,f_sr,f_fsena,f_lfevpa,f_fevpa,f_fseng,f_fevpg,f_fgrnd,f_olrg,f_rnet,f_t_grnd,f_scv,f_fsno,f_sigf,f_tref,f_qref,f_rss,f_t_soisno,f_wliq_soisno,f_wice_soisno,f_h2osoi,f_us10m,f_vs10m "$input_name" "$output_name"
        out_list="$out_list $output_name"
    done
done

cdo mergetime $out_list colmoff_2001-2023_monmean_gravel.nc
rm $out_list


