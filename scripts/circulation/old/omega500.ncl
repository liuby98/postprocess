  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
  begin

  G = 9.81
  Rd = 287.04
  Cp = 1004.
  RCP = Rd/Cp
  p0 = 100000.
  EPS = 0.622
  f = addfile("/share/home/dq135/wrfpost/wrfout_1991_2023_MAM_seasmean.nc", "r")
  f1= addfile("/share/home/dq135/wrfpost/wrfpost_w_1991_2023_MAM_seasmean.nc","r")
  f2= addfile("/share/home/dq135/wrfpost/wrfpost_t_1991_2023_MAM_seasmean.nc","r")
  f3= addfile("/share/home/dq135/wrfpost/wrfpost_q_1991_2023_MAM_seasmean.nc","r")
  ff= addfile("/share/home/dq135/wrfpost/wrfout_nog_1991_2023_MAM_seasmean.nc", "r")
  f4= addfile("/share/home/dq135/wrfpost/wrfpost_w_nog_1991_2023_MAM_seasmean.nc","r")
  f5= addfile("/share/home/dq135/wrfpost/wrfpost_t_nog_1991_2023_MAM_seasmean.nc","r")
  f6= addfile("/share/home/dq135/wrfpost/wrfpost_q_nog_1991_2023_MAM_seasmean.nc","r")
  f_input = addfile("/share/home/dq135/wrfpost/wrfinput_d01", "r")
  ; 获取网格的纬度和经度
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
  
  t = f2->air({500},:,:)
  qv = f3->QVAPOR({500},:,:)
  t_exp = f5->air({500},:,:)
  qv_exp = f6->QVAPOR({500},:,:)
  pres = 50000.
  tk = (t+300.) * (( pres / p0 )^RCP)
  tk_exp = (t_exp+300.) * (( pres / p0 )^RCP)

  w = f1->W({500},:,:)
  w_exp = f4->W({500},:,:)
  omega = -G*pres/ (Rd*((tk*(EPS + qv))/ (EPS*(1. + qv))))*w
  omega_exp = -G*pres/ (Rd*((tk_exp*(EPS + qv_exp))/ (EPS*(1. + qv_exp))))*w_exp
  ; omega=omega*1.0e+2
  omega_diff = omega_exp - omega
  omega_diff!0 = "lat"
  omega_diff!1 = "lon"
  omega_diff&lat = lat
  omega_diff&lon = lon
  printVarSummary(omega_diff)
  printMinMax(omega_diff,0)
  ; exit
  wks  = gsn_open_wks("png","omega500_diff_MAM")        
  res0                  = True
  res0@tiMainString     = " "
  res0@gsnDraw          = False
  res0@gsnFrame         = False  
  res0@mpFillOn         = False
  res0@mpOutlineOn      = True
  res0@gsnAddCyclic     = False

  res0@tfDoNDCOverlay          = True    ; 原生投影：不变换数据 

  res0@tmXTOn                 = False  
  res0@tmYROn                 = False

  res0@mpLimitMode            = "Corners" 
  res0@mpLeftCornerLatF       = 48
  res0@mpLeftCornerLonF       = 60
  res0@mpRightCornerLatF      = 10
  res0@mpRightCornerLonF      = 135

  res0@mpProjection           = "LambertConformal"
  res0@mpLambertParallel1F    = 30.
  res0@mpLambertParallel2F    = 60.
  res0@mpLambertMeridianF     = 110.
  res0@mpCenterLatF           = 35.1778 
  res0@mpCenterLonF           = 110.   
  res0@mpDataBaseVersion      = "MediumRes"        
  res0@pmTickMarkDisplayMode  = "Always"
  res0@mpOutlineOn = True   ;开启海岸线绘制
  res0@mpGeophysicalLineColor  = "grey60"   ;设置海岸线颜色
  res0@mpGeophysicalLineThicknessF = 0.3  
  res0@mpOutlineBoundarySets = "GeophysicalAndUSStates"  ; 改善中国边界（如果适用）
 
  res0@cnLinesOn         = False  
  res0@cnFillOn          = True
  res0@cnFillPalette     = "BlWhRe"  
  res0@cnLineLabelsOn    = True 
  res0@cnLineLabelsOn    = False 
  ; res0@lbOrientation    = "Vertical"
  res0@cnLevelSelectionMode = "ManualLevels"  
  res0@cnLineLabelInterval  =  1
  res0@cnMinLevelValF       = -0.08
  res0@cnMaxLevelValF       =  0.08 
  res0@cnLevelSpacingF      =  0.01
  res0@gsnLeftString     = " 500hPa omega "
  res0@gsnRightString    ="Pa/s" ;*10~S~-2~N~
  plot=gsn_csm_contour_map(wks,omega_diff,res0)

  dir                    = "../shapefile_China/"
  fShape1Name            = dir + "china.shp"
  fShape2Name            = dir + "south_china_sea.shp"
  ; fShape3Name            = dir + "province.shp"

  resLine                     = True
  resLine@gsLineThicknessF    = 1.2
  resLine@gsLineColor         = "black"
  line  = gsn_add_shapefile_polylines(wks, plot, fShape1Name,resLine)
  scs   = gsn_add_shapefile_polylines(wks, plot, fShape2Name,resLine)
  ; province  = gsn_add_shapefile_polylines(wks, plot, fShape3Name,resLine)

  draw(plot)
  frame(wks)          
end