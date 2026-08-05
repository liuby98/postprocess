#!/bin/bash
#BSUB -J wind             # 作业名称，改为描述当前任务的名称
#BSUB -q normal          # 队列名
#BSUB -o logout          # 标准输出文件
#BSUB -e logout          # 错误输出文件，合并到标准输出
#BSUB -n 48               # 核数
#BSUB -R span[ptile=48]   # 每个节点分配 48 个核即可

# 1. 加载你的环境变量（确保里面包含环境变量和路径）
conda activate plot

# 2. 记录一下开始时间（可选）
echo "python script started at: $(date)"

# 3. 运行脚本，并将屏幕上的打印信息输出到log文件中
python u200_spacial.py &>log

# 4. 记录结束时间（可选）
echo "python script finished at: $(date)"