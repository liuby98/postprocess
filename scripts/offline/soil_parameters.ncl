begin
  path="/share/home/dq135/CRESM_0307S92/SCRIPTS/ICBCScripts/PrepCoLM/CN_15km/"
  f1=addfile(path+"CoLMrun_nogravel/unstructured_cwrf_CN_15km_nogravel/landdata/diag/soil_parameters_2005.nc","r")
  f2=addfile(path+"CoLMrun_gravel/unstructured_cwrf_CN_15km_gravel/landdata/diag/soil_parameters_2005.nc","r")
  lat1d = f1->lat({15:55})  ;8:60
  lon1d = f1->lon({70:140}) ;58:162
  lat1d@units = "degrees_north"
  lon1d@units = "degrees_east"
  printVarSummary(lat1d)
  printVarSummary(lon1d)
  diff = new((/8,400,700/),float)
  diff!0 = "layer"
  diff!1 = "lat"
  diff!2 = "lon"
  diff&lat = lat1d
  diff&lon = lon1d
  param = (/"theta_s_l1","theta_s_l2","theta_s_l3","theta_s_l4","theta_s_l5","theta_s_l6","theta_s_l7","theta_s_l8"/)

  ; f3=addfile("/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/colmoff_2001-2023_monmean_nogravel.nc","r")
  ; f4=addfile("/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/colmoff_2001-2023_monmean_gravel.nc","r")
  ; t_ctl_mon = tofloat(f3->f_h2osoi(:,:,:,0:7))
  ; t_exp_mon = tofloat(f4->f_h2osoi(:,:,:,0:7))

;   f_input = addfile("/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01","r")
;   nlat = f_input->XLAT(0,:,0)
;   nlon = f_input->XLONG(0,0,:)
;   lat2d = f_input->XLAT(0,:,:)
;   lon2d = f_input->XLONG(0,:,:)
;   printVarSummary(lat2d)
;   printMinMax(lat2d,0)
;   printVarSummary(lon2d)
;   printMinMax(lon2d,0)

;   ctl_spr = new((/17,342,462,8/),float)
;   exp_spr = new((/17,342,462,8/),float)
;   ctl_sum = new((/17,342,462,8/),float)
;   exp_sum = new((/17,342,462,8/),float)
;  do i=0,16
;    ctl_spr(i,:,:,:)=dim_avg_n_Wrap(t_ctl_mon(6*i:6*i+2,:,:,:),0)
;    exp_spr(i,:,:,:)=dim_avg_n_Wrap(t_exp_mon(6*i:6*i+2,:,:,:),0)
;    ctl_sum(i,:,:,:)=dim_avg_n_Wrap(t_ctl_mon(6*i+3:6*i+5,:,:,:),0)
;    exp_sum(i,:,:,:)=dim_avg_n_Wrap(t_exp_mon(6*i+3:6*i+5,:,:,:),0)
;  end do

;   theta_ctl_spr = dim_avg_n_Wrap(ctl_spr(:,:,:,:),0)  
;   theta_exp_spr = dim_avg_n_Wrap(exp_spr(:,:,:,:),0)  
;   theta_ctl_sum = dim_avg_n_Wrap(ctl_sum(:,:,:,:),0)
;   theta_exp_sum = dim_avg_n_Wrap(exp_sum(:,:,:,:),0)

;   theta_ctl_spr = where(theta_ctl_spr .le. -1e+35, theta_ctl_spr@_FillValue, theta_ctl_spr)
;   theta_exp_spr = where(theta_exp_spr .le. -1e+35, theta_exp_spr@_FillValue, theta_exp_spr)
;   theta_ctl_sum = where(theta_ctl_sum .le. -1e+35, theta_ctl_sum@_FillValue, theta_ctl_sum)
;   theta_exp_sum = where(theta_exp_sum .le. -1e+35, theta_exp_sum@_FillValue, theta_exp_sum)

; theta_ctl_spr!0 = "south_north"
;   theta_ctl_spr!1 = "west_east"
;   theta_ctl_spr!2 = "layer"
;   theta_ctl_spr&south_north@units = "degrees_north"
;   theta_ctl_spr&west_east@units = "degrees_east"

;   theta_exp_spr!0 = "south_north"
;   theta_exp_spr!1 = "west_east"
;   theta_exp_spr!2 = "layer"
;   theta_exp_spr&south_north@units = "degrees_north"
;   theta_exp_spr&west_east@units = "degrees_east"

;   theta_ctl_sum!0 = "south_north"
;   theta_ctl_sum!1 = "west_east"
;   theta_ctl_sum!2 = "layer"
;   theta_ctl_sum&south_north@units = "degrees_north"
;   theta_ctl_sum&west_east@units = "degrees_east"

;   theta_exp_sum!0 = "south_north"
;   theta_exp_sum!1 = "west_east"
;   theta_exp_sum!2 = "layer"
;   theta_exp_sum&south_north@units = "degrees_north"
;   theta_exp_sum&west_east@units = "degrees_east"

; theta_ctl_spr_re = theta_ctl_spr(layer|:, south_north|:, west_east|:)
;   theta_exp_spr_re = theta_exp_spr(layer|:, south_north|:, west_east|:)
;   theta_ctl_sum_re = theta_ctl_sum(layer|:, south_north|:, west_east|:)
;   theta_exp_sum_re = theta_exp_sum(layer|:, south_north|:, west_east|:)
;   printVarSummary(theta_ctl_spr_re_re)
;   printMinMax(theta_ctl_spr_re,0)
 
; theta_s_ctl_interp = new((/8,342,462/), float)
;   theta_s_exp_interp = new((/8,342,462/), float)
;   theta_s_ctl_interp!0 = "layer"
;   theta_s_ctl_interp!1 = "south_north"
;   theta_s_ctl_interp!2 = "west_east"
;   theta_s_exp_interp!0 = "layer"
;   theta_s_exp_interp!1 = "south_north"
;   theta_s_exp_interp!2 = "west_east"

;   tkdry_ctl_interp = new((/8,342,462/), float)
;   tkdry_exp_interp = new((/8,342,462/), float)
;   tkdry_ctl_interp!0 = "layer"
;   tkdry_ctl_interp!1 = "south_north"
;   tkdry_ctl_interp!2 = "west_east"
;   tkdry_exp_interp!0 = "layer"
;   tkdry_exp_interp!1 = "south_north"
;   tkdry_exp_interp!2 = "west_east"

;   tksatu_ctl_interp = new((/8,342,462/), float)
;   tksatu_exp_interp = new((/8,342,462/), float)
;   tksatu_ctl_interp!0 = "layer"
;   tksatu_ctl_interp!1 = "south_north"
;   tksatu_ctl_interp!2 = "west_east"
;   tksatu_exp_interp!0 = "layer"
;   tksatu_exp_interp!1 = "south_north"
;   tksatu_exp_interp!2 = "west_east"

;   tksatf_ctl_interp = new((/8,342,462/), float)
;   tksatf_exp_interp = new((/8,342,462/), float)
;   tksatf_ctl_interp!0 = "layer"
;   tksatf_ctl_interp!1 = "south_north"
;   tksatf_ctl_interp!2 = "west_east"
;   tksatf_exp_interp!0 = "layer"
;   tksatf_exp_interp!1 = "south_north"
;   tksatf_exp_interp!2 = "west_east"

do k=0,7
  param_ctl = tofloat(f1->$param(k)$({15:55},{70:140},:))
  param_exp = tofloat(f2->$param(k)$({15:55},{70:140},:))
  
  param_ctl_avg = dim_avg_n_Wrap(param_ctl, 2) ; *(0.039*(dim_avg_n_Wrap(param_ctl, 2)^(-2.2)))
  param_exp_avg = dim_avg_n_Wrap(param_exp, 2) ; *(0.039*(dim_avg_n_Wrap(param_exp, 2)^(-2.2)))

  diff(k,:,:) = param_exp_avg - param_ctl_avg
  delete([/param_ctl,param_exp/])
  delete([/param_ctl_avg, param_exp_avg/])

end do

;   printVarSummary(theta_s_ctl_interp)
;   printMinMax(theta_s_ctl_interp,0)
; exit
  f_mask=addfile("/share/home/dq117/CN05.1/CN05.1_Tm_2020_daily_025x025.nc","r")
  nlat = f_mask->lat
  nlon = f_mask->lon
  obs_mon = f_mask->tm
  obs = dim_avg_n_Wrap(obs_mon,0)
  mask_cn05 = where(ismissing(obs), 0, 1)
  mask_target = linint2_Wrap(nlon, nlat, mask_cn05, False, lon1d, lat1d, 0)
  mask_target = where(mask_target.gt.0.5, 1, 0)

  diff_masked = mask(diff, mask_target, 1)
  diff_masked@_FillValue = default_fillvalue(typeof(diff))
  copy_VarCoords(diff, diff_masked)
  diff_masked@units = "cm³/cm³" ;cm/day W/(m·K)
  printVarSummary(diff_masked)
  printMinMax(diff_masked,0)
  ; exit

  plot = new(8,graphic)
  wks = gsn_open_wks("pdf","./figs/theta_s_diff") ;heat_capacity ksat_hydraulic
do k=0,7
  res = True                            ; res probability plots
  res@gsnDraw             = False       ; Do not draw plot
  res@gsnFrame            = False       ; Do not advance frome
  res@tiMainString        = ""
  res@gsnRightString      = "cm~S~3~N~/cm~S~3~N~" ;W/(m K) KJ/(m3 K) cm/day
  res@gsnLeftString       = "l="+(k+1)
  res@gsnLeftStringOrthogonalPosF = -0.003
  res@gsnRightStringOrthogonalPosF = -0.003

  res@cnLinesOn           = False
  res@mpFillOn            = False
  res@cnFillOn            = True
  res@lbLabelBarOn        = True
  res@mpOutlineOn         = True
  res@gsnAddCyclic        = False
;   res@vpWidthF            = 0.6
;   res@vpHeightF           = 0.7

  res@tmXTOn                 = False  
  res@tmYROn                 = False

  res@mpMinLonF               = 73
  res@mpMaxLonF               = 136
  res@mpCenterLonF            = 105  
  res@mpMinLatF               = 10
  res@mpMaxLatF               = 60

  res@mpDataBaseVersion      = "MediumRes"        
  res@pmTickMarkDisplayMode  = "Always"
  res@mpOutlineOn            = True       ;开启海岸线绘制
  res@mpGeophysicalLineColor = "grey60"   ;设置海岸线颜色
  res@mpGeophysicalLineThicknessF = 0.3  
  res@mpOutlineBoundarySets = "GeophysicalAndUSStates"  ; 改善中国边界（如果适用）

  res@cnMissingValFillColor = -1         ; 缺失值透明
  res@cnFillDrawOrder       = "PreDraw"  ; 先绘制填充，再叠加 shapefile
  res@cnFillPalette         = "MPL_BrBG" ;cmocean_delta cmp_b2r BlueWhiteOrangeRed WhiteYellowOrangeRed 
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnMinLevelValF        = -0.13         ; set min contour level
  res@cnMaxLevelValF        =  0.13        ; set max contour level
  res@cnLevelSpacingF       =  0.01          ; set contour spacing
  res@cnInfoLabelOn         = False
  plot(k) = gsn_csm_contour_map(wks,diff_masked(k,:,:),res)

end do

  dir                  = "../../shapefile_China/"
  fShape1Name          = dir + "china.shp"
  fShape2Name          = dir + "south_china_sea.shp"
  fShape3Name          = dir + "province.shp"

  resLine                     = True
  resLine@gsLineThicknessF    = 0.8
  resLine@gsLineColor         = "black"
  line  = gsn_add_shapefile_polylines(wks, plot, fShape1Name,resLine)
  scs   = gsn_add_shapefile_polylines(wks, plot, fShape2Name,resLine)
  province  = gsn_add_shapefile_polylines(wks, plot, fShape3Name,resLine)

  resP = True
  resP@gsnPanelMainString    = "Saturated Water Content" ;Saturated Hydraulic Conductivity Thermal conductivity of Saturated unfrozen soil
  resP@gsnPanelFigureStrings = (/"a)","b)","c)","d)","e)","f)","g)","h)"/) ; 设置具体标签内容
  resP@gsnPanelFigureStringsFontHeightF = 0.015       ; 设置标签字体大小
  resP@gsnPanelFigureStringsPerimOn     = False       ; 去掉标签周围的边框 (可选)
  resP@gsnPanelFigureStringsBackgroundFillColor = -1  ; 背景透明 (可选, -1为透明)


  ; 设置标签的位置 (默认在左上角，常用 "TopLeft")
  resP@amJust = "TopLeft"
  ; 微调标签位置 (可选)
  ; resP@amOrthogonalPosF = -0.5   ; 垂直方向移动 (-0.5 为向上移)
  ; resP@amParallelPosF   = -0.5   ; 水平方向移动 (-0.5 为向左移)
  gsn_panel(wks, plot, (/2,4/), resP)

draw(plot)
frame(wks)
end
