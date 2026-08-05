#!/bin/bash

# 定义每个月的天数数组
monday=(31 28 31 30 31 30 31 31 30 31 30 31)

# 定义输入和输出基础路径
base_input_path="/share/home/dq135/CRESM_0307S92/CASES/RUN_"
base_output_path="/share/home/dq135/draw/scripts/temporary/"

out_list = ''
# 开始年份到结束年份的循环
for year in {2001..2023}; do  # 这里你可以更改年份范围
    # 判断是否为闰年，如果是闰年则二月有29天
    if (( year % 4 == 0 )); then
        monday[1]=29
    else
        monday[1]=28
    fi
    
    # 处理6月到8月
    for mon in {3..5}; do
        days_in_month=${monday[$((mon-1))]}  # 获取该月的天数

        # 遍历每一天
        for day in $(seq 1 $days_in_month); do
            # 格式化月份和天数为两位数
            mon_str=$(printf "%02d" $mon)
            day_str=$(printf "%02d" $day)

            # 拼接输入和输出文件路径
            input_name="${base_input_path}${year}_gravel/wrfout_d01_${year}-${mon_str}-${day_str}_00:00:00"
            output_name="${base_output_path}wrfout_${year}-${mon_str}-${day_str}_gravel.nc"
           
            # 使用cdo命令提取T2M变量（执行实际操作时取消注释）
            cdo daymean "$input_name" "$output_name"

	    # 将输出文件名添加到列表
            out_list="$out_list $output_name"
        done
    done
done

# 打印所有生成的输出文件
# echo "Files to merge: $out_list"
cdo mergetime $out_list wrfout_2001-2023_MAM_daymean_gravel.nc
#rm $out_list

