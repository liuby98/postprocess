load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
begin
  ; 打开 WRF 输出文件
  ; a = addfile("wrfout_1991_2023_MAM_seasmean.nc", "r")
  ; ua = a->T(:,:,:,:)
  ; printVarSummary(ua)
  ; printMinMax(ua, 0)
  ; exit
  ; 定义需要处理的气压层（hPa）
  ; 17 pressure levels (hPa): 1000.0, 925.0, 850.0, 700.0, 600.0, 500.0, 400.0, \
  ;                           300.0, 250.0, 200.0, 150.0, 100.0, 70.0, 50.0, 30.0, \
  ;                           20.0, 10.0
  pressure_levels = (/1000, 925, 850, 700, 600, 500, 400, 300, 200/) ; 

  f_input = addfile("/share/home/dq135/wrfpost/wrfinput_d01", "r") ;wrfinput_1991_2023.nc
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; year=ispan(1991,2023,1)
  ; do k=1,2
  ; f = addfile("/share/home/dq135/wrfpost/wrfout_"+year(k)+"_MAM_seasmean.nc", "w")
  ; f = addfile("wrfout_1991_2023_MAM_seasmean_years.nc", "w")
  ; hgt = f_input->HGT(:,:,:)
  ; pb = f_input->PB(:,:,:,:) 
  ; tb = f_input->TB(:,:,:,:)
  ; printVarSummary(hgt)
  ; printVarSummary(pb)
  ; printVarSummary(tb)
  ; f->HGT = hgt
  ; f->PB = pb
  ; f->TB = tb
  ; end do
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  f = addfile("/share/home/dq135/wrfpost/wrfout_nog_1991_2023_MAM_seasmean.nc", "r")
  times = f->Times
  ntimes = dimsizes(times)

  ; 获取网格的纬度和经度
  lat = f_input->XLAT(0,:,0)
  lon = f_input->XLONG(0,0,:)
  lon@standard_name = "longitude"
  lon@long_name = "longitude"
  lon@units = "degrees_east"
  lon@axis = "X"

  lat@standard_name = "latitude"
  lat@long_name = "latitude"
  lat@units = "degrees_north"
  lat@axis = "Y"

  ; 获取纬度、经度的维度长度
  dim_lat = dimsizes(lat)    
  dim_lon = dimsizes(lon)    
; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;          U V W风场后处理              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; u0 = f->V(0,0,:,:)
  ; ; 循环处理每个气压层
  ;   do ip = 0, dimsizes(pressure_levels)-1
  ;     plevel = pressure_levels(ip)
  ;     print("Processing pressure level: " + plevel + " hPa")

  ;     ; 获取当前时间步的 P 和 PB（单位：Pa）
  ;     P = wrf_user_getvar(f,"P",-1)
  ;     PB = wrf_user_getvar(f,"PB",-1)
  ;     printVarSummary(P)
  ;     printVarSummary(PB)
  ;     ; 计算总气压（单位：hPa）
  ;     pressure = (P + PB)/100.
  ;     copy_VarCoords(P, pressure)
  ;     pressure@units = "hPa"
  ;     pressure@cell_methods = "Times: mean"
  ;     pressure@FieldType = 104
  ;     pressure@MemoryOrder = "XYZ"
  ;     pressure@description = "pressure"
  ;     ; printVarSummary(pressure)

  ;     ; 获取当前时间步的 U 风速分量
  ;     u_stag = wrf_user_getvar(f,"V",-1)  
  ;     printVarSummary(u_stag)
  ;     ; 将 U 从交错网格转换为非交错网格
  ;     u = wrf_user_unstagger(u_stag, u_stag@stagger)  
  ;     printVarSummary(u)

  ;     ; 检查地表压力，判断插值可行性
  ;     psfc = wrf_user_getvar(f, "PSFC", -1) / 100.0  ; 转换为 hPa
  ;     if (min(psfc) .lt. plevel) then
  ;         print("Warning: Some grid points have surface pressure below " + plevel + " hPa")
  ;     end if
      
  ;     u_interp = new((/342,462/),float)
  ;     ; 将 U 插值到指定的气压层
  ;     ; u_interp = wrf_user_interp_level(u, pressure, plevel, False)
  ;     ; 或使用带外推法函数将 U 插值到指定的气压层:1000-600hPa(珠峰)
  ;     opts             = True
  ;     opts@extrapolate = True 
  ;     vert_coord       = "pressure"
  ;     uu = wrf_user_vert_interp(f, u, vert_coord, plevel, opts)
  ;     u_interp(:,:) = uu(0,0,:,:)
  ;     copy_VarCoords(u0, u_interp)
  ;     delete(u_interp@_FillValue)     
  ;     printVarSummary(u_interp)
  ;     printMinMax(u_interp, 0)

  ;     ; 生成输出文件名
  ;     output_file = "v" + toint(plevel) + "_1991_2023_MAM_seasmean.nc" ;
  ;     system("rm -f " + output_file)

  ;     ; 将插值结果保存为新的 NetCDF 文件
  ;     ofile = addfile(output_file, "c")
  ;     ofile->V = u_interp
  ;     delete(u_interp)
  ;     delete(uu)
  ;   end do

  ; print("All pressure levels processed.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;              位势高度场            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; p0 = f->PH(0,0,:,:)
  ; ; 循环处理每个气压层
  ; do ip = 0, dimsizes(pressure_levels)-1
  ;   plevel = pressure_levels(ip)
  ;   print("Processing pressure level: " + plevel + " hPa")

  ;     ; 获取当前时间步的 P 和 PB（单位：Pa）
  ;     P = wrf_user_getvar(f,"P",-1)
  ;     PB = wrf_user_getvar(f,"PB",-1)
  ;     printVarSummary(P)
  ;     printVarSummary(PB)
  ;     ; 计算总气压（单位：hPa）
  ;     pressure = (P + PB)/100.
  ;     copy_VarCoords(P, pressure)
  ;     pressure@units = "hPa"
  ;     pressure@cell_methods = "Times: mean"
  ;     pressure@FieldType = 104
  ;     pressure@MemoryOrder = "XYZ"
  ;     pressure@description = "pressure"
  ;     ; printVarSummary(pressure)
      
  ;     ; 获取当前时间步的总位势
  ;     ph = f->PH
  ;     phb = f->PHB
  ;     ; printVarSummary(ph)

  ;     geopt_w = ph
  ;     geopt_w = ph+phb
  ;     ; 从交错网格转换为非交错网格
  ;     geopt = wrf_user_unstagger(geopt_w, geopt_w@stagger)
  ;     ght = geopt/9.8

  ;     ; 插值到指定的气压层
  ;     ; ght_interp = wrf_interp_3d_z(ght_stag, pressure, plevel)
  ;     ; h_interp = wrf_user_interp_level(ght, pressure, plevel, False)

  ;     opts             = True
  ;     opts@extrapolate = True 
  ;     vert_coord       = "pressure"
  ;     hh = wrf_user_vert_interp(f, ght, vert_coord, plevel, opts)
  ;     h_interp = new((/342,462/),float)
  ;     h_interp(:,:) = hh(0,0,:,:)
  ;     copy_VarCoords(p0, h_interp)
  ;     delete(h_interp@_FillValue)
  ;     h_interp@units = "gpm" ;
  ;     h_interp@FieldType = 104 ;
  ;     h_interp@description = "geopotential height" ;
  ;     h_interp@stagger = "" ;
  ;     h_interp@coordinates = "XLONG XLAT" ;

  ;   printVarSummary(h_interp)
  ;   printMinMax(h_interp, 0)

  ;   ; 生成输出文件名
  ;   output_file = "h" + toint(plevel) + "_1991_2023_JJA_seasmean.nc" ;
  ;   system("rm -f " + output_file)

  ;   ; 将插值结果保存为新的 NetCDF 文件
  ;   ofile = addfile(output_file, "c")
  ;   ofile->GH = h_interp
  ;   delete(h_interp)
  ;   ; delete(hh)
  ; end do

  ; print("All pressure levels processed.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;              温度场           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; t0 = f->T(0,0,:,:)
  ; ; 循环处理每个气压层
  ; do ip = 0, dimsizes(pressure_levels)-1
  ;   plevel = pressure_levels(ip)
  ;   print("Processing pressure level: " + plevel + " hPa")

  ;     ; 获取当前时间步的 P 和 PB（单位：Pa）
  ;     P = wrf_user_getvar(f,"P",-1)
  ;     PB = wrf_user_getvar(f,"PB",-1)
  ;     printVarSummary(P)
  ;     printVarSummary(PB)
  ;     ; 计算总气压（单位：hPa）
  ;     pressure = (P + PB)/100.
  ;     copy_VarCoords(P, pressure)
  ;     pressure@units = "hPa"
  ;     pressure@cell_methods = "Times: mean"
  ;     pressure@FieldType = 104
  ;     pressure@MemoryOrder = "XYZ"
  ;     pressure@description = "pressure"
  ;     ; printVarSummary(pressure)
      
  ;     ; 获取当前时间步的总温度
  ;     t  = wrf_user_getvar(f,"T",-1)
  ;     tb = wrf_user_getvar(f,"TB",-1)
  ;     printVarSummary(t)
  ;     printVarSummary(tb)
  ;     air = t
  ;     air = t+tb

  ;     t_interp = new((/342,462/),float)
  ;     ; 插值到指定的气压层
  ;     ; t_interp = wrf_user_interp_level(air, pressure, plevel, False)
  ;     ; 外推法
  ;     opts             = True
  ;     opts@extrapolate = True 
  ;     vert_coord       = "pressure"
  ;     tt = wrf_user_vert_interp(f, air, vert_coord, plevel, opts)
  ;     t_interp(:,:) = tt(0,0,:,:)
  ;     copy_VarCoords(t0, t_interp)
  ;     delete(t_interp@_FillValue)
  ;     t_interp@units = "K" ;
  ;     t_interp@FieldType = 104 ;
  ;     t_interp@description = "temperature" 
  ;     t_interp@stagger = "" ;
  ;     t_interp@coordinates = "XLONG XLAT" 

  ;   printVarSummary(t_interp)
  ;   printMinMax(t_interp, 0)

  ;   ; 生成输出文件名
  ;   output_file = "t" + toint(plevel) + "_1991_2023_MAM_seasmean.nc" ;
  ;   system("rm -f " + output_file)

  ;   ; 将插值结果保存为新的 NetCDF 文件
  ;   ofile = addfile(output_file, "c")
  ;   ofile->air = t_interp
  ;   delete(t_interp)
  ; end do

  ; print("All pressure levels processed.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;              比湿场           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  qv0 = f->QVAPOR(0,0,:,:)
  ; 循环处理每个气压层
  do ip = 0, dimsizes(pressure_levels)-1
    plevel = pressure_levels(ip)
    print("Processing pressure level: " + plevel + " hPa")

      ; 获取当前时间步的 P 和 PB（单位：Pa）
      P = wrf_user_getvar(f,"P",-1)
      PB = wrf_user_getvar(f,"PB",-1)
      printVarSummary(P)
      printVarSummary(PB)
      ; 计算总气压（单位：hPa）
      pressure = (P + PB)/100.
      copy_VarCoords(P, pressure)
      pressure@units = "hPa"
      pressure@cell_methods = "Times: mean"
      pressure@FieldType = 104
      pressure@MemoryOrder = "XYZ"
      pressure@description = "pressure"
      ; printVarSummary(pressure)
      
      ; 获取当前时间步的比湿
      qv  = wrf_user_getvar(f,"QVAPOR",-1)

      q_interp = new((/342,462/),float)
      ; 插值到指定的气压层
      ; q_interp = wrf_user_interp_level(qv, pressure, plevel, False)
      ; 外推法
      opts             = True
      opts@extrapolate = True 
      vert_coord       = "pressure"
      q = wrf_user_vert_interp(f, qv, vert_coord, plevel, opts)
      q_interp(:,:) = q(0,0,:,:)
      copy_VarCoords(qv0, q_interp)
      delete(q_interp@_FillValue)
      q_interp@units = "kg kg-1"
      q_interp@cell_methods = "Times: mean" 
      q_interp@FieldType = 104 
      q_interp@MemoryOrder = "XYZ"
      q_interp@description = "water vapor mixing ratio" 
      q_interp@coordinates = "XLONG XLAT" 
   
    printVarSummary(q_interp)
    printMinMax(q_interp, 0)

    ; 生成输出文件名
    output_file = "q" + toint(plevel) + "_nog_1991_2023_MAM_seasmean.nc" ;
    system("rm -f " + output_file)

    ; 将插值结果保存为新的 NetCDF 文件
    ofile = addfile(output_file, "c")
    ofile->QVAPOR = q_interp
    delete(q_interp)
  end do

  print("All pressure levels processed.")

end
