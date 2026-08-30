#!/bin/bash
#BSUB -J ncl             # 作业名称，改为描述当前任务的名称
#BSUB -q normal          # 队列名
#BSUB -o %J.out          # 标准输出文件
#BSUB -e %J.out          # 错误输出文件，合并到标准输出
#BSUB -n 48               # 核数
#BSUB -R span[ptile=48]   # 每个节点分配 48 个核即可

# 1. 加载你的环境变量（确保里面包含了 NCL 的环境变量和路径）
source ~/.bashrc_zhwei

# 2. 记录一下开始时间（可选，方便看处理了多久）
echo "NCL script started at: $(date)"

# 3. 运行 NCL 脚本，并将屏幕上的打印信息输出到 ncl_runlog 文件中
ncl parameter_interp.ncl &> ncllog

# 4. 记录结束时间（可选）
echo "NCL script finished at: $(date)"