from netCDF4 import Dataset
import numpy as np

# 输入和输出文件路径
input_file = 'wrfout_nog_1991_2023_MAM_seasmean_years.nc'
output_file = 'wrfout_nog_1991_2023_MAM_seasmean_years_mod.nc'

# 打开原始 NetCDF 文件
src = Dataset(input_file, 'r')

# 创建新的 NetCDF 文件，覆盖已有文件
dst = Dataset(output_file, 'w', format='NETCDF4', clobber=True)

# 跟踪已创建的维度，避免重复
created_dims = set()

# 复制维度，替换 lev 和 lev_2
for dimname, dim in src.dimensions.items():
    new_dimname = dimname
    if dimname == 'lev':
        new_dimname = 'bottom_top'
    elif dimname == 'lev_2':
        new_dimname = 'bottom_top_stag'

    # 仅在维度未创建时创建
    if new_dimname not in created_dims:
        dst.createDimension(new_dimname, len(dim) if not dim.isunlimited() else None)
        created_dims.add(new_dimname)

# 复制变量并调整维度
for varname, var in src.variables.items():
    # 替换变量的维度名称
    new_dims = tuple(
        'bottom_top' if d == 'lev' else 'bottom_top_stag' if d == 'lev_2' else d
        for d in var.dimensions
    )

    # 创建变量
    var_out = dst.createVariable(varname, var.dtype, new_dims)

    # 复制变量数据
    var_out[:] = var[:]

    # 复制变量属性，并检查是否需要添加 cell_methods
    for attr in var.ncattrs():
        var_out.setncattr(attr, var.getncattr(attr))
    
    # 为 HGT, PB, TB 添加 cell_methods 属性（如果目标格式需要）
    if varname in ['HGT', 'PB', 'TB'] and 'cell_methods' not in var.ncattrs():
        var_out.setncattr('cell_methods', 'Times: mean')

# 关闭文件
src.close()
dst.close()

print(f"成功将文件转换为新格式，输出文件为：{output_file}")
