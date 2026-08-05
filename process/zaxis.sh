#!/bin/bash

# 定义气压层（按降序排列）
levels=(1000 925 850 700 600 500 400 300)

# 获取所有 wrfout_1991_2023_u*_MAM.nc 文件并按气压层排序
files=()
for level in "${levels[@]}"; do
    file="q${level}_nog_1991_2023_MAM_seasmean.nc"
    if [ -f "$file" ]; then
        files+=("$file")
    else
        echo "Warning: File $file not found, skipping..."
    fi
done

# 检查是否找到文件
if [ ${#files[@]} -eq 0 ]; then
    echo "Error: No u*_1991_2023_MAM_seasmean.nc!"
    exit 1
fi

# 为每个文件添加 level 维度
for i in ${!files[@]}; do
    file=${files[$i]}
    level=${levels[$i]}
    
    echo "Processing $file with level $level hPa"

    # 修改 zaxis 文件，设置 level
    sed -i "s/^levels[ \t]*=[ \t]*[0-9]*/levels    = $level/" zaxis

    # 使用 cdo setzaxis 添加 level 维度
    cdo setzaxis,zaxis "$file" "temp_u${level}.nc"
    
    # 检查 cdo 命令是否成功
    if [ $? -ne 0 ]; then
        echo "Error: Failed to process $file with cdo"
        exit 1
    fi
done

# 合并所有临时文件（按 levels 顺序）
echo "Merging temporary files..."
temp_files=()
for level in "${levels[@]}"; do
    if [ -f "temp_u${level}.nc" ]; then
        temp_files+=("temp_u${level}.nc")
    fi
done

# 检查临时文件是否完整
if [ ${#temp_files[@]} -eq 0 ]; then
    echo "Error: No temporary files found for merging!"
    exit 1
fi

# 合并文件
cdo merge "${temp_files[@]}" wrfpost_q_nog_1991_2023_MAM_seasmean.nc

# 检查合并是否成功
if [ $? -ne 0 ]; then
    echo "Error: Failed to merge files"
    exit 1
fi

# 确保 level 顺序正确（从高到低）
#echo "Setting level order..."
#rm -f wrfpost_U_1991_2023_MAM.nc
#cdo invertlev wrfpost_U_1991_2023_MAM_test.nc wrfpost_U_1991_MAM_seasmon.nc

# 检查倒序是否成功
if [ $? -ne 0 ]; then
    echo "Error: Failed to invert levels"
    exit 1
fi

# 清理临时文件
echo "Cleaning up temporary files..."
rm -f temp_*.nc wrfpost_*_test.nc

#pls=(1000 700 600 400 300)
#for pl in "${pls[@]}"; do
#    rm -f "wrfout_1991_u${pl}_MAM.nc"
#done

echo "Processing completed successfully. Output: wrfpost_q_nog_1991_2023_MAM_seasmean.nc"
