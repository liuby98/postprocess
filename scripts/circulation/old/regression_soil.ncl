begin
  f1 = addfile("/share/home/dq135/draw/data/unstructured_cwrf_CN_15km_1991-2023_JJA_soivar.nc","r")
  f2 = addfile("/share/home/dq135/draw/data/unstructured_cwrf_CN_15km_nogravel_1991-2023_JJA_soivar.nc","r")
  f3 = addfile("/share/home/dq135/wrfpost/wrfout_1991_2023_JJA_seasmean_years.nc","r")
  f4 = addfile("/share/home/dq135/wrfpost/wrfout_nog_1991_2023_JJA_seasmean_years.nc","r")
  f_input = addfile("/share/home/dq135/wrfpost/wrfinput_d01", "r")
  lat2d = f_input->XLAT(0,:,:)
  lon2d = f_input->XLONG(0,:,:)
  lon = f_input->XLONG(0,0,:)
  lat = f_input->XLAT(0,:,0)
  lon@standard_name = "longitude"
  lon@long_name = "longitude"
  lon@units = "degrees_east"
  lon@axis = "X"

  lat@standard_name = "latitude"
  lat@long_name = "latitude"
  lat@units = "degrees_north"
  lat@axis = "Y"

  f = addfile("/share/home/dq135/reference/era5_1991_2023_monthly_025x025.nc","r")
  lat1d = f->latitude({15:55})
  lon1d = f->longitude({70:140})
  z_obs_mon = f->z(:,{200},{15:55},{70:140})
  z = new((/33,161,281/),float)
  do i=0,32
     z(i,:,:) = dim_avg_n_Wrap(z_obs_mon(6*i+3:6*i+5,:,:),0)
  end do
  hgt = z/9.81
  copy_VarCoords(z,hgt)
  hgt!0 = "time"
  hgt@_FillValue = 9.96921e+36
  printVarSummary(hgt)
  printMinMax(hgt, 0)
  ; exit
; 求区域平均、距平、标准化得到时间序列
  t_ctl_mon = f1->f_t_soisno(:,:,:,12)
  t_ctl_mon = (/t_ctl_mon - 273.112/)
  t_ctl_mon@_FillValue = 9.96921e+36 
  t_ctl_mon@missing_value = 9.96921e+36
  t_ctl_mon = where(abs(t_ctl_mon) .gt. 1.e+10, t_ctl_mon@_FillValue, t_ctl_mon)
  t_ctl = new((/33,342,462/),double)

  t_exp_mon = f2->f_t_soisno(:,:,:,12)
  t_exp_mon = (/t_exp_mon - 273.15/)
  t_exp_mon@_FillValue = 9.96921e+36
  t_exp_mon@missing_value = 9.96921e+36
  t_exp_mon = where(abs(t_exp_mon) .gt. 1.e+10, t_exp_mon@_FillValue, t_exp_mon)
  t_exp = new((/33,342,462/),double)

  do i=0,32
     t_ctl(i,:,:) = dim_avg_n_Wrap(t_ctl_mon(3*i:3*i+2,:,:),0)
     t_exp(i,:,:) = dim_avg_n_Wrap(t_exp_mon(3*i:3*i+2,:,:),0)
  end do
  t_ctl@lat2d = lat2d
  t_ctl@lon2d = lon2d
  t_ctl_ip = rcm2rgrid_Wrap(lat2d,lon2d,t_ctl,lat1d,lon1d,1)

  t_exp@lat2d = lat2d
  t_exp@lon2d = lon2d
  t_exp_ip = rcm2rgrid_Wrap(lat2d,lon2d,t_exp,lat1d,lon1d,1)

  sim_ctl = dim_avg_n_Wrap(t_ctl_ip(:,{25:40},{75:105}),(/1,2/))
  ctl_rmv = dim_rmvmean_n_Wrap(sim_ctl,0)
  ti1 = dim_standardize_Wrap(ctl_rmv, 1)

  sim_exp = dim_avg_n_Wrap(t_exp_ip(:,{25:40},{75:105}),(/1,2/))
  exp_rmv = dim_rmvmean_n_Wrap(sim_exp,0)
  ti2 = dim_standardize_Wrap(exp_rmv, 1)

; 位势高度
  p0 = f3->PH(:,0,:,:)
      ; 获取当前时间步的总位势
      ph = f3->PH
      phb = f3->PHB
      ph_exp = f4->PH
      phb_exp = f4->PHB
      ; printVarSummary(ph)

      ; 从交错网格转换为非交错网格
      geopt_w = ph
      geopt_w = ph+phb
      geopt = wrf_user_unstagger(geopt_w, geopt_w@stagger)
      ght = geopt/9.8

      geopt_w_exp = ph_exp
      geopt_w_exp = ph_exp+phb_exp
      geopt_exp = wrf_user_unstagger(geopt_w_exp, geopt_w_exp@stagger)
      ght_exp = geopt_exp/9.8

  ; 循环处理每个气压层
  pressure_levels = (/200/) ; 
  do ip = 0, dimsizes(pressure_levels)-1
      plevel = pressure_levels(ip)
      print("Processing pressure level: " + plevel + " hPa")
      ; 插值到指定的气压层
      opts             = True
      opts@extrapolate = True 
      vert_coord       = "pressure"
      hh = wrf_user_vert_interp(f3, ght, vert_coord, plevel, opts)
      h_interp = new((/33,342,462/),float)
      h_interp(:,:,:) = hh(:,0,:,:)
      copy_VarCoords(p0, h_interp)
      delete(h_interp@_FillValue)
      h_interp!0 = "time" ;
      h_interp@units = "gpm" 

      hh_exp = wrf_user_vert_interp(f4, ght_exp, vert_coord, plevel, opts)
      h_interp_exp = new((/33,342,462/),float)
      h_interp_exp(:,:,:) = hh_exp(:,0,:,:)
      copy_VarCoords(p0, h_interp_exp)
      delete(h_interp_exp@_FillValue)
      h_interp_exp!0 = "time" ;
      h_interp_exp@units = "gpm" 
    ; delete(hh)
  end do
  h_interp_anom = dim_rmvmean_n_Wrap(h_interp,0)
  h_interp_anom@lat2d = lat2d
  h_interp_anom@lon2d = lon2d
  h_interp_ip = rcm2rgrid_Wrap(lat2d,lon2d,h_interp_anom,lat1d,lon1d,1)

  h_interp_anom_exp = dim_rmvmean_n_Wrap(h_interp_exp,0)
  h_interp_anom_exp@lat2d = lat2d
  h_interp_anom_exp@lon2d = lon2d
  h_interp_ip_exp = rcm2rgrid_Wrap(lat2d,lon2d,h_interp_anom_exp,lat1d,lon1d,1)
  printVarSummary(h_interp_anom)
  printMinMax(h_interp_anom, 0)

  hgt_anom = dim_rmvmean_n_Wrap(hgt,0)
; exit
; 回归   
  rc = regCoef_n(ti1, hgt_anom, 0, 0)
  rc@long_name = "regression coefficient"
  copy_VarCoords(hgt(0,:,:),rc)
  rchgt = regCoef_n(ti1, h_interp_ip, 0, 0)
  rchgt@long_name = "regression coefficient"
  rchgt_exp = regCoef_n(ti2, h_interp_ip_exp, 0, 0)
  rchgt_exp@long_name = "regression coefficient"
  diff = rchgt - rchgt_exp
  copy_VarCoords(rc,diff)
  printVarSummary(rc)
  printMinMax(rc,0)

  diff!0 = "lat"
  diff!1 = "lon"
  diff&lat = lat1d
  diff&lon = lon1d
  printVarSummary(diff)
  printMinMax(diff, 0)
;t-test
  confi_hgt = rc
  p = student_t(rc@tval, rc@nptxy - 2)
  dims = dimsizes(hgt)
  confi_hgt = onedtond(p, (/dims(1),dims(2)/))
  copy_VarCoords(hgt(0,:,:),confi_hgt)
  confi_hgt!0 = "lat"
  confi_hgt!1 = "lon"
  confi_hgt&lat = lat1d
  confi_hgt&lon = lon1d

  printVarSummary(confi_hgt)
  printMinMax(confi_hgt, 0)
; exit
  wks  = gsn_open_wks("pdf","sub_reg_dep_hgt200_JJA_diff")
  gsn_define_colormap(wks,"BlueWhiteOrangeRed")          

res = True
  res@gsnDraw  = False         
  res@gsnFrame = False 
  res@gsnAddCyclic   = False       
  res@tiMainString   = " "
  res@gsnLeftString  = " "
  res@gsnRightString = " "
   
  res@mpDataBaseVersion      = "MediumRes"        
  res@pmTickMarkDisplayMode  = "Always"
  res@mpOutlineOn = True   ;开启海岸线绘制
  res@mpGeophysicalLineColor  = "grey60"   ;设置海岸线颜色
  res@mpGeophysicalLineThicknessF = 0.3  
  res@mpOutlineBoundarySets = "GeophysicalAndUSStates"  ; 改善中国边界（如果适用）     
  res@mpFillOn       = False       
  res@mpMinLatF      = 10
  res@mpMaxLatF      = 60          
  res@mpCenterLonF   = 105
  res@mpMinLonF      = 70
  res@mpMaxLonF      = 140
  res@mpGeophysicalLineThicknessF = 0.5
  res@pmTickMarkDisplayMode = "Always"

  res@cnFillOn              = True   
  res@cnInfoLabelOn         = False
  res@cnLinesOn             = False
  res@cnLineLabelsOn        = False
  res@cnSmoothingOn         = True
  res@lbLabelBarOn          = True
  ; res@cnSmoothingDistanceF  = 0.001
  ; res@cnSmoothingTensionF   = -2
  res@pmLabelBarOrthogonalPosF = 0.2
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnMinLevelValF        = -3.0
  res@cnMaxLevelValF        =  3.0
  res@cnLevelSpacingF       =  0.1
  plot = gsn_csm_contour_map(wks, diff, res)

res0 = True   
  res0@gsnDraw  = False         
  res0@gsnFrame = False 
  res0@gsnAddCyclic   = False       
  res0@tiMainString   = " "
  res0@gsnLeftString  = " "
  res0@gsnRightString = " "       
  res0@cnFillOn             = True         
  res0@cnLinesOn            = False        
  res0@cnLevelSelectionMode = "ExplicitLevels"         
  res0@cnMonoFillPattern    = False       
  res0@cnLevels             = (/0.05, 0.1/)
  res0@cnFillPatterns       = (/17, 14, -1/)
  res0@cnFillScaleF         = 1.5
  res0@cnFillDotSizeF       = 0.003
  res0@cnFillColors         = (/"gray60","gray80","transparent"/)
  res0@cnFillDrawOrder      = "postdraw"        
  res0@cnLineLabelsOn       = False      
  res0@cnInfoLabelOn        = False 
  res0@lbLabelBarOn         = False
  base0 = gsn_csm_contour(wks, confi_hgt, res0)
  ; overlay(plot,base0)

dir         = "../shapefile_China/"
  fShape1Name = dir + "china.shp"
  fShape2Name = dir + "south_china_sea.shp"
  fShape3Name = dir + "province.shp"
  resLine                  = True
  resLine@gsLineThicknessF = 0.8
  resLine@gsLineColor      = "black"
  line = gsn_add_shapefile_polylines(wks, plot, fShape1Name,resLine)
  scs  = gsn_add_shapefile_polylines(wks, plot, fShape2Name,resLine)
  province = gsn_add_shapefile_polylines(wks, plot, fShape3Name,resLine)

  draw(plot)
  frame(wks)
end