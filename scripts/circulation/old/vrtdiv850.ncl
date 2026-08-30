begin
  f1= addfile("/share/home/dq135/wrfpost/wrfpost_u_1991_2023_MAM_seasmean.nc","r")
  f2= addfile("/share/home/dq135/wrfpost/wrfpost_v_1991_2023_MAM_seasmean.nc","r")
  f3= addfile("/share/home/dq135/wrfpost/wrfpost_u_nog_1991_2023_MAM_seasmean.nc","r")
  f4= addfile("/share/home/dq135/wrfpost/wrfpost_v_nog_1991_2023_MAM_seasmean.nc","r")
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

  u=f1->U({850},:,:)
  v=f2->V({850},:,:)
  u_exp=f3->U({850},:,:)
  v_exp=f4->V({850},:,:)

  div=uv2dv_cfd(u, v, lat, lon, 2)*10^5
  vrt=uv2vr_cfd(u, v, lat, lon, 2)*10^5
  div_exp=uv2dv_cfd(u_exp, v_exp, lat, lon, 2)*10^5
  vrt_exp=uv2vr_cfd(u_exp, v_exp, lat, lon, 2)*10^5

  vrt!0 = "lat"
  vrt!1 = "lon"
  vrt&lat = lat
  vrt&lon = lon
  div!0 = "lat"
  div!1 = "lon"
  div&lat = lat
  div&lon = lon
  vrt_exp!0 = "lat"
  vrt_exp!1 = "lon"
  vrt_exp&lat = lat
  vrt_exp&lon = lon
  div_exp!0 = "lat"
  div_exp!1 = "lon"
  div_exp&lat = lat
  div_exp&lon = lon
  vrt_diff = vrt_exp - vrt
  div_diff = div_exp - div
  copy_VarCoords(vrt, vrt_diff)
  copy_VarCoords(div, div_diff)
printVarSummary(div_diff)
printMinMax(div_diff,0)
; exit
  ; plot=new(2,graphic)
  wks  = gsn_open_wks("png","div850_diff_MAM")

  res0              = True
  res0@gsnDraw      = False
  res0@gsnFrame     = False  
  res0@mpFillOn     = False
  res0@mpOutlineOn  = True
  res0@gsnAddCyclic = False  
  ; res0@mpMinLonF         = 85
  ; res0@mpMaxLonF         = 135
  ; res0@mpMinLatF         = 10
  ; res0@mpMaxLatF         = 55
  ; res0@mpCenterLonF      = 105

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

  res0@cnFillOn                = True
  res0@cnLinesOn               = False
  res0@mpFillOn                = False
  res0@mpOutlineOn             = True
  ; res0@lbOrientation         = "Vertical"
  res0@cnLinesOn         = False  
  res0@cnFillOn 	       = True
  res0@cnFillPalette     = "CBR_coldhot"  ;cmocean_balance
  ;res0@cnLineLabelsOn         =True 
  res0@cnLevelSelectionMode = "ManualLevels"

  ; res0@cnMinLevelValF       = -10.
  ; res0@cnMaxLevelValF       =  10. 
  ; res0@cnLevelSpacingF      =  1.  
  ; res0@cnLineLabelInterval  = 1
  ; vrt@long_name = "(a) ctl_vorticity"
  ; vrt@units     = "*1e-05 /s" 
  ; plot = gsn_csm_contour_map(wks,vrt,res0)

  res0@cnMinLevelValF       = -4.   
  res0@cnMaxLevelValF       =  4.
  res0@cnLevelSpacingF      =  1.
  div_diff@long_name = "(c) divergence difference"
  div_diff@units     = "*1e-05 /s"  
  plot=gsn_csm_contour_map(wks,div_diff,res0)
  
  ; plres=True
  ; plres@txString     = "MAM"
  ; gsn_panel(wks, plot, (/2,1/), plres)

  dir                    = "../shapefile_China/"
  fShape1Name            = dir + "china.shp"
  fShape2Name            = dir + "south_china_sea.shp"
;   fShape3Name            = dir + "province.shp"

resLine                     = True
  resLine@gsLineThicknessF    = 1.2
  resLine@gsLineColor         = "black"
  line  = gsn_add_shapefile_polylines(wks, plot, fShape1Name,resLine)
  scs   = gsn_add_shapefile_polylines(wks, plot, fShape2Name,resLine)
;   province  = gsn_add_shapefile_polylines(wks, plot, fShape3Name,resLine)
draw(plot)
frame(wks)               
end

