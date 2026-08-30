begin
  path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
  f1=addfile(path+"wrfpost_gh200_2001_2017_MAM_nogravel.nc","r")
  f2=addfile(path+"wrfpost_gh200_2001_2017_JJA_nogravel.nc","r")
  f3=addfile(path+"wrfpost_gh200_2001_2017_MAM_gravel.nc","r")
  f4=addfile(path+"wrfpost_gh200_2001_2017_JJA_gravel.nc","r")

  f_input = addfile("/share/home/dq135/wrfpost/wrfinput_d01", "r")
  lon2d = f_input->XLONG(0,:,:)
  lat2d = f_input->XLAT(0,:,:)

  h_ctl_mam=f1->geopt(:,0,:,:) / 9.8
  h_ctl_jja=f2->geopt(:,0,:,:) / 9.8
  h_exp_mam=f3->geopt(:,0,:,:) / 9.8
  h_exp_jja=f4->geopt(:,0,:,:) / 9.8

  h_ctl_spr = new((/17,342,462/),float)
  h_exp_spr = new((/17,342,462/),float)
  h_ctl_sum = new((/17,342,462/),float)
  h_exp_sum = new((/17,342,462/),float)
 do i=0,16
   h_ctl_spr(i,:,:)=dim_avg_n_Wrap(h_ctl_mam(92*i:92*i+91,:,:),0)
   h_exp_spr(i,:,:)=dim_avg_n_Wrap(h_exp_mam(92*i:92*i+91,:,:),0)
   h_ctl_sum(i,:,:)=dim_avg_n_Wrap(h_ctl_jja(92*i:92*i+91,:,:),0)
   h_exp_sum(i,:,:)=dim_avg_n_Wrap(h_exp_jja(92*i:92*i+91,:,:),0)
 end do

  h_c_spr = dim_avg_n_Wrap(h_ctl_spr,0)
  h_e_spr = dim_avg_n_Wrap(h_exp_spr,0)
  h_c_sum = dim_avg_n_Wrap(h_ctl_sum,0)
  h_e_sum = dim_avg_n_Wrap(h_exp_sum,0)

  h_c_spr@lat2d = lat2d
  h_e_spr@lat2d = lat2d
  h_c_sum@lat2d = lat2d
  h_e_sum@lat2d = lat2d

  h_c_spr@lon2d = lon2d
  h_e_spr@lon2d = lon2d
  h_c_sum@lon2d = lon2d
  h_e_sum@lon2d = lon2d

  diff_spr  = h_e_spr - h_c_spr
  copy_VarCoords(h_c_spr, diff_spr)

  diff_sum  = h_e_sum - h_c_sum
  copy_VarCoords(h_c_sum, diff_sum)

  printVarSummary(h_c_spr) 
  printMinMax(h_c_spr,0)
  printVarSummary(diff_sum) 
  printMinMax(diff_sum,0)
; exit

plot = new(6,graphic)
wks=gsn_open_wks("pdf","./figs/hgt200_diff")
  res                   = True
  res@tiMainString      = ""
  res@gsnDraw                 = False
  res@gsnFrame                = False
  res@gsnAddCyclic            = False
  res@tfDoNDCOverlay          = True    ; 原生投影：不变换数据 

  res@tmXTOn                 = False  
  res@tmYROn                 = False

  res@mpLimitMode            = "Corners"          ; choose region of map
  res@mpLeftCornerLatF       = 48
  res@mpLeftCornerLonF       = 60
  res@mpRightCornerLatF      = 10
  res@mpRightCornerLonF      = 135

  res@mpProjection           = "LambertConformal"
  res@mpLambertParallel1F    = 30.
  res@mpLambertParallel2F    = 60.
  res@mpLambertMeridianF     = 110.
  res@mpCenterLatF           = 35.1778 
  res@mpCenterLonF           = 110.   

  res@mpFillOn                = False
  res@mpOutlineOn             = True

  res@mpDataBaseVersion      = "MediumRes"        
  res@pmTickMarkDisplayMode  = "Always"
  res@mpOutlineOn = True   ;开启海岸线绘制
  res@mpGeophysicalLineColor  = "grey60"   ;设置海岸线颜色
  res@mpGeophysicalLineThicknessF = 0.3  
  res@mpOutlineBoundarySets = "GeophysicalAndUSStates"  ; 改善中国边界（如果适用）
  
  res@cnFillOn             = True
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnInfoLabelOn        = False
  res@lbLabelBarOn         = True      ; 为差异图打开色标
  res@lbOrientation        = "Horizontal"
  res@pmLabelBarOrthogonalPosF = 0.1   ; 调整色标的上下位置
  
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       =  11350. 
  res@cnMaxLevelValF       =  12600.
  res@cnLevelSpacingF      =  25.
  res@cnFillPalette        = "BlAqGrYeOrReVi200" ;BlueWhiteOrangeRed 适合差异图的红蓝发散色表
  res@gsnRightString       = "gpm"
  res@gsnLeftString           = "(a)CTL MAM"
  plot(0) = gsn_csm_contour_map(wks,h_c_spr,res)
  res@gsnLeftString           = "(b)CTL JJA"
  plot(1) = gsn_csm_contour_map(wks,h_c_sum,res)
  res@gsnLeftString           = "(c)EXP MAM"
  plot(2) = gsn_csm_contour_map(wks,h_e_spr,res)
  res@gsnLeftString           = "(d)EXP JJA"
  plot(3) = gsn_csm_contour_map(wks,h_e_sum,res)

  ; res@lbLabelBarOn         = False
  res@cnMinLevelValF       = -3.0  
  res@cnMaxLevelValF       =  3.0
  res@cnLevelSpacingF      =  0.1
  res@cnFillPalette        = "NCV_blu_red"
  res@gsnLeftString           = "(e)EXP MAM - CTL MAM"
  plot(4) = gsn_csm_contour_map(wks, diff_spr, res)
  res@gsnLeftString           = "(f)EXP JJA - CTL JJA"
  plot(5) = gsn_csm_contour_map(wks, diff_sum, res)

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

  resP@gsnPanelMainString    = " "
  resP@gsnPanelLabelBar      = False
  resP@lbBoxLinesOn          = False      ; 去掉色块边框
  resP@pmLabelBarWidthF      = 0.8        ; 色标条宽度 (0.0-1.0)
  resP@pmLabelBarHeightF     = 0.06       ; 色标条高度
  resP@pmLabelBarOrthogonalPosF = -0.01   ; 垂直位置微调 (负数向上，正数向下)
  resP@lbBoxEndCapStyle         = "TriangleBothEnds"
  ; --- 设置单位标题 ---
  resP@lbTitleString         = "200 hPa Mean Geopotential Height Field Difference(2001-2017)" ; 标题
  resP@lbTitlePosition       = "Bottom"      ; 标题放在右侧
  resP@lbTitleFontHeightF    = 0.012        ; 标题字体大小
  resP@lbTitleDirection      = "Across"     ; 标题方向
  resP@lbTitleOffsetF        = 0.025         ; 标题与色标的距离
  gsn_panel(wks, plot, (/3,2/), resP)

end
