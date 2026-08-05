begin
  path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
  f1=addfile(path+"wrfout_2001-2023_MAM_daymean_nogravel.nc","r")
  f2=addfile(path+"wrfout_2001-2023_MAM_daymean_gravel.nc","r")
  f3=addfile(path+"wrfout_2001-2023_JJA_daymean_nogravel.nc","r")
  f4=addfile(path+"wrfout_2001-2023_JJA_daymean_gravel.nc","r")
  ; AHFX:"average or accumulative sensible  heat flux" ;
  ; ALFX:"average or accumulative latent    heat flux" ;
  ; AGFX:"average or accumulative ground    heat flux" ;
  ; ALWUPS:"average upwelling   lw flux at surface"    ;
  t_ctl_spr = tofloat(f1->AHFX(:,:,:))
  t_exp_spr = tofloat(f2->AHFX(:,:,:))
  t_ctl_sum = tofloat(f3->AHFX(:,:,:))
  t_exp_sum = tofloat(f4->AHFX(:,:,:))

  f_input = addfile("/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01","r")
  nlat = f_input->XLAT(0,:,0)
  nlon = f_input->XLONG(0,0,:)
  lat2d = f_input->XLAT(0,:,:)
  lon2d = f_input->XLONG(0,:,:)

  ctl_spr = new((/17,342,462/),float)
  exp_spr = new((/17,342,462/),float)
  ctl_sum = new((/17,342,462/),float)
  exp_sum = new((/17,342,462/),float)
 do i=0,16
   ctl_spr(i,:,:)=dim_avg_n_Wrap(t_ctl_spr(92*i:92*i+91,:,:),0)
   exp_spr(i,:,:)=dim_avg_n_Wrap(t_exp_spr(92*i:92*i+91,:,:),0)
   ctl_sum(i,:,:)=dim_avg_n_Wrap(t_ctl_sum(92*i:92*i+91,:,:),0)
   exp_sum(i,:,:)=dim_avg_n_Wrap(t_exp_sum(92*i:92*i+91,:,:),0)
 end do

  diff = new((/2, 342, 462/), float)
  diff(0,:,:) = dim_avg_n_Wrap(exp_spr,0) - dim_avg_n_Wrap(ctl_spr,0)
  diff(1,:,:) = dim_avg_n_Wrap(exp_sum,0) - dim_avg_n_Wrap(ctl_sum,0)

  diff!1 = "lat"
  diff!2 = "lon"
  diff&lat = nlat
  diff&lon = nlon
  diff&lat@units = "degrees_north"
  diff&lon@units = "degrees_east"

  f=addfile("/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc","r")
  t_obs_mon = f->tm
  lat1d = f->lat
  lon1d = f->lon
  diff@lat2d = lat2d
  diff@lon2d = lon2d
  diff_ip = rcm2rgrid_Wrap(lat2d,lon2d,diff,lat1d,lon1d,1)
  mask_2d = t_obs_mon(0,:,:)
  mask_3d = conform(diff_ip, mask_2d, (/1,2/))
  diff_ip = where(ismissing(mask_3d), diff_ip@_FillValue, diff_ip)
  printVarSummary(diff_ip)
  printMinMax(diff_ip, 0)

; f_wrf=addfile("/share/home/dq117/CN05.1/CN05.1_Tm_2020_daily_025x025.nc","r")
;   lat1d = f_wrf->lat
;   lon1d = f_wrf->lon
;   t_obs_mon = f_wrf->tm
;   t_obs = dim_avg_n_Wrap(t_obs_mon,0)
;   lat = t_obs_mon(0,:,:)
;   lon = t_obs_mon(0,:,:)
;   mask_cn05 = where(ismissing(t_obs), 0, 1)
;   mask_target = linint2_Wrap(lon1d, lat1d, mask_cn05, False, nlon, nlat, 0)
;   mask_target = where(mask_target.gt.0.5, 1, 0)

;   diff_masked = mask(diff, mask_target, 1)
;   diff_masked@_FillValue = default_fillvalue(typeof(diff))
;   copy_VarCoords(diff, diff_masked)
;   printVarSummary(diff_masked)
;   printMinMax(diff_masked, 0)
  ; exit
  wks = gsn_open_wks("pdf","./wrfout_latent_heat_diff")
  plot= new(2,graphic)
  res                         = True
  res@gsnDraw                 = False
  res@gsnFrame                = False
  res@gsnAddCyclic            = False
  ; res@gsnMaximize             = True               ; enlarge image  

  res@mpLimitMode            = "Corners"          ; choose region of map
  res@mpLeftCornerLatF       = 48
  res@mpLeftCornerLonF       = 60
  res@mpRightCornerLatF      = 10
  res@mpRightCornerLonF      = 135
  res@tmXTOn                 = False  
  res@tmYROn                 = False

  res@mpProjection           = "LambertConformal"
  res@mpLambertParallel1F    = 30.
  res@mpLambertParallel2F    = 60.
  res@mpLambertMeridianF     = 110.

  res@cnFillOn              = True
  res@cnLinesOn             = False
  res@mpFillOn              = False
  res@mpOutlineOn           = True
  res@mpGridAndLimbOn       = True
  res@mpGridSpacingF        = 10          ; 可选 5、10、15
  res@mpGridLineDashPattern = 2           ; 2=短虚线（最好看）
  res@mpGridLineColor       = "grey80"    ; 灰色
  res@mpGridLineThicknessF  = 1.0

  res@mpDataBaseVersion      = "MediumRes"        
  res@pmTickMarkDisplayMode  = "Always"
  res@mpOutlineOn = True                   ;开启海岸线绘制
  res@mpGeophysicalLineColor  = "grey60"   ;设置海岸线颜色
  res@mpGeophysicalLineThicknessF = 0.3
  res@mpOutlineDrawOrder    = "PreDraw"

  res@lbLabelBarOn            = False
  res@lbBoxLinesOn            = False
  res@lbOrientation           = "Horizontal"
  ; res@lbBoxEndCapStyle        = "TriangleBothEnds"
  res@lbTitlePosition         = "Bottom"
  res@lbTitleDirection        = "Across"
  res@lbTitleAngleF           =  0
  res@lbTitleString           = " "
  ; res@cnFillPalette           = "MPL_BrBG" ; 
  gsn_define_colormap(wks, "NCV_blu_red")   
  res@gsnSpreadColorStart = 2
  res@gsnSpreadColorEnd   = 256
  res@gsnRightString          = "W/m~S~2~N~" ;~S~o~N~C
  ; res@gsnRightStringFontHeightF= 0.025
  ; res@gsnLeftStringFontHeightF= 0.025
  ;res@gsnLeftStringOrthogonalPosF = -0.001
  res@tiYAxisFontHeightF      = 0.02
  res@tiYAxisAngleF           = 90
  res@cnLineLabelsOn        = False
  res@cnLevelSelectionMode    = "ManualLevels"
  res@cnMaxLevelValF          =  3.0  
  res@cnMinLevelValF          = -3.0  
  res@cnLevelSpacingF         =  0.1
  ; res@cnLevelSelectionMode   = "ExplicitLevels"
  ; res@cnLevels               = (/-0.5,-0.4,-0.3,-0.2,-0.1,0,0.5,1.0,1.5,2.0,2.5/)
  res@gsnLeftString           = "MAM"
  plot(0) = gsn_csm_contour_map(wks,diff_ip(0,:,:),res)
  res@gsnLeftString           = "JJA"
  plot(1) = gsn_csm_contour_map(wks,diff_ip(1,:,:),res)

  txres = True
  txres@txFontHeightF = 0.02
  txres@txFontColor = "black"
  c1 =  gsn_add_text(wks, plot(0), "EXP-CTL",lon1d({127}),lat1d({17}), txres)
  c2 =  gsn_add_text(wks, plot(1), "EXP-CTL",lon1d({127}),lat1d({17}), txres)

  dir                    = "../../shapefile_China/"
  fShape1Name            = dir + "china.shp"
  fShape2Name            = dir + "south_china_sea.shp"
  fShape3Name            = dir + "river.nc"

  resLine                     = True
  resLine@gsLineThicknessF    = 0.4
  resLine@gsLineColor         = "black"
  line  = gsn_add_shapefile_polylines(wks, plot, fShape1Name,resLine)
  scs   = gsn_add_shapefile_polylines(wks, plot, fShape2Name,resLine)
  resLine@gsLineColor         = "blue"
  river  = gsn_add_shapefile_polylines(wks, plot, fShape3Name,resLine)

  resP = True
  ; 设置标签的位置 (默认在左上角，常用 "TopLeft")
  ; resP@amJust = "TopLeft"   
  ; resP@gsnPanelFigureStrings = (/"a)","b)"/) ; 设置具体标签内容
  ; resP@gsnPanelFigureStringsFontHeightF = 0.015       ; 设置标签字体大小
  ; resP@gsnPanelFigureStringsPerimOn     = False       ; 去掉标签周围的边框 (可选)
  ; resP@gsnPanelFigureStringsBackgroundFillColor = -1  ; 背景透明 (可选, -1为透明)
  ; 微调标签位置 (可选)
  ; resP@amOrthogonalPosF = -0.5   ; 垂直方向移动 (-0.5 为向上移)
  ; resP@amParallelPosF   = -0.5   ; 水平方向移动 (-0.5 为向左移)

  resP@gsnPanelMainString    = " " ;Upwelling Longwave flux at surface(W/m2) /Sensible/ /Ground 
  resP@gsnPanelLabelBar      = True
  resP@lbBoxLinesOn          = False      ; 去掉色块边框
  resP@pmLabelBarWidthF      = 0.8        ; 色标条宽度 (0.0-1.0)
  resP@pmLabelBarHeightF     = 0.06       ; 色标条高度
  resP@pmLabelBarOrthogonalPosF = -0.01   ; 垂直位置微调 (负数向上，正数向下)
  resP@lbBoxEndCapStyle         = "TriangleBothEnds"
  ; --- 设置单位标题 ---
  resP@lbTitleString         = "Latent Heat Flux Difference" ; 标题
  resP@lbTitlePosition       = "Bottom"      ; 标题放在右侧
  resP@lbTitleFontHeightF    = 0.012        ; 标题字体大小
  resP@lbTitleDirection      = "Across"     ; 标题方向
  resP@lbTitleOffsetF        = 0.025         ; 标题与色标的距离
  gsn_panel(wks, plot, (/1,2/), resP)

end
