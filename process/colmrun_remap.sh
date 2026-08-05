#!/bin/bash
# 使用相对路径 python ../../../../../SCRIPTS/PostProcScripts/CoLM_Remap.py

START_YEAR=2001
END_YEAR=2023

# 脚本所在目录
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
LOGFILE="${BASE_DIR}/Remap_2001_2023_$(date +%Y%m%d_%H%M%S).log"

echo "================= CoLM_Remap 批量处理 ==================" > $LOGFILE
echo "开始时间: $(date)" >> $LOGFILE
echo "年份范围: $START_YEAR - $END_YEAR" >> $LOGFILE
echo "日志文件: $LOGFILE" >> $LOGFILE
echo "========================================================" >> $LOGFILE

for year in $(seq $START_YEAR $END_YEAR); do
    for type in nogravel ; do
        # 目标目录
        TARGET_DIR="/share/home/dq135/CRESM_0307S92/CASES/RUN_${year}_${type}/CoLMrun/unstructured_cwrf_CN_15km_${type}/history"
        
        if [ -d "$TARGET_DIR" ]; then
            echo "[$year $type] 正在处理..." | tee -a $LOGFILE
            cd "$TARGET_DIR" || continue

            # 相对路径调用 python 脚本
            python ../../../../../SCRIPTS/PostProcScripts/CoLM_Remap.py \
                   -lats 342 -lons 462 -cpus 1 \
                   &>> $LOGFILE

            if [ $? -eq 0 ]; then
                echo "[$year $type] 成功完成" | tee -a $LOGFILE
            else
                echo "[$year $type] 执行失败！请查看日志" | tee -a $LOGFILE
            fi

            # 回到当前目录，继续下一个
            cd "$BASE_DIR"
            echo "--------------------------------------------------" >> $LOGFILE
        else
            echo "[$year $type] 目录不存在，跳过: $TARGET_DIR" | tee -a $LOGFILE
        fi
    done
done

echo "全部年份处理完毕！结束时间: $(date)" | tee -a $LOGFILE
echo "完整日志见: $LOGFILE"

