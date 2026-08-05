#!/bin/bash

# 设置路径
temp_dir="/share/home/dq135/draw/scripts/temp/"
target_dir="/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
final_output="${target_dir}wrfout_2001-2023_MAM_daymean_nogravel.nc"

# 清理可能遗留的旧年度文件和最终文件（可选）
#rm -f ${temp_dir}wrfout_*_MAM_daymean_gravel.nc
rm -f "$final_output"

# 用于收集每年合并后的文件列表
yearly_files=""

echo "开始分年合并 MAM 日平均文件..."

for year in {2001..2023}; do
    # 本年 MAM 的文件列表
    year_mam_files=""

    # 查找该年 3月、4月、5月的所有 _gravel.nc 文件
    for month in 03 04 05; do
        year_mam_files="$year_mam_files ${temp_dir}wrfout_${year}-${month}-??_nogravel.nc"
    done

    # 年度输出文件
    yearly_output="${temp_dir}wrfout_${year}_MAM_daymean_nogravel.nc"

    echo "正在合并 $year 年 MAM ..."
    cdo mergetime $year_mam_files "$yearly_output"

    # 检查是否合并成功
    if [ $? -eq 0 ] && [ -f "$yearly_output" ]; then
        echo "$year 年合并成功 -> $yearly_output"
        yearly_files="$yearly_files $yearly_output"

        # 可选：删除已合并的每日文件，节省空间
         rm -f $year_mam_files
        # echo "已删除 $year 年每日文件"
    else
        echo "警告：$year 年合并失败，请检查文件是否存在！"
    fi
done

# 最终合并所有年的 MAM 文件
echo ""
echo "正在进行最终合并：2001-2023 MAM ..."
cdo mergetime $yearly_files "$final_output"

if [ $? -eq 0 ]; then
    echo "全部完成！"
    echo "最终文件生成：$final_output"
    echo "包含 2001-2023 年所有 MAM 季节的日平均数据"
else
    echo "最终合并失败，请检查年度文件是否完整。"
fi

# 可选：清理年度中间文件（如果你只想要最终文件）
 echo "正在清理年度中间文件..."
 rm -f $yearly_files
 echo "清理完成，仅保留最终文件。"
