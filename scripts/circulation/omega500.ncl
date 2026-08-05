begin
  path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
  f1=addfile(path+"wrfpost_omega500_2001_2017_MAM_nogravel.nc","r")
  f2=addfile(path+"wrfpost_omega500_2001_2017_JJA_nogravel.nc","r")
  f3=addfile(path+"wrfpost_omega500_2001_2017_MAM_gravel.nc","r")
  f4=addfile(path+"wrfpost_omega500_2001_2017_JJA_gravel.nc","r")

  f_input = addfile("/share/home/dq135/wrfpost/wrfinput_d01", "r")
  lon2d = f_input->XLONG(0,:,:)
  lat2d = f_input->XLAT(0,:,:)

  ; --- 读取垂直速度 Omega (Pa/s) ---

  w_ctl_mam = f1->omega(:,0,:,:)
  w_ctl_jja = f2->omega(:,0,:,:)
  w_exp_mam = f3->omega(:,0,:,:)
  w_exp_jja = f4->omega(:,0,:,:)

  w_ctl_spr = new((/17,342,462/),float)
  w_exp_spr = new((/17,342,462/),float)
  w_ctl_sum = new((/17,342,462/),float)
  w_exp_sum = new((/17,342,462/),float)
  
  do i=0,16
    w_ctl_spr(i,:,:)=dim_avg_n_Wrap(w_ctl_mam(92*i:92*i+91,:,:),0)
    w_exp_spr(i,:,:)=dim_avg_n_Wrap(w_exp_mam(92*i:92*i+91,:,:),0)
    w_ctl_sum(i,:,:)=dim_avg_n_Wrap(w_ctl_jja(92*i:92*i+91,:,:),0)
    w_exp_sum(i,:,:)=dim_avg_n_Wrap(w_exp_jja(92*i:92*i+91,:,:),0)
  end do

  w_c_spr = dim_avg_n_Wrap(w_ctl_spr,0)
  w_e_spr = dim_avg_n_Wrap(w_exp_spr,0)
  w_c_sum = dim_avg_n_Wrap(w_ctl_sum,0)
  w_e_sum = dim_avg_n_Wrap(w_exp_sum,0)

  w_c_spr@lat2d = lat2d
  w_e_spr@lat2d = lat2d
  w_c_sum@lat2d = lat2d
  w_e_sum@lat2d = lat2d

  w_c_spr@lon2d = lon2d
  w_e_spr@lon2d = lon2d
  w_c_sum@lon2d = lon2d
  w_e_sum@lon2d = lon2d

  diff_spr  = w_e_spr - w_c_spr
  copy_VarCoords(w_c_spr, diff_spr)

  diff_sum  = w_e_sum - w_c_sum
  copy_VarCoords(w_c_sum, diff_sum)

  ; ======================================================================
  ; 添加掩膜
  ; ======================================================================
  f_mask = addfile("/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc","r")
  t_obs_mon = f_mask->tm(0,:,:)
  lat1d = f_mask->lat
  lon1d = f_mask->lon
  
  ; 将一维网格的掩膜插值到 WRF 的二维网格上
  t_mask_2d = rgrid2rcm_Wrap(lat1d, lon1d, t_obs_mon, lat2d, lon2d, 0)
  
  ; 对 4 个基础场应用掩膜
  w_c_spr = where(ismissing(t_mask_2d), w_c_spr@_FillValue, w_c_spr)
  w_c_sum = where(ismissing(t_mask_2d), w_c_sum@_FillValue, w_c_sum)
  w_e_spr = where(ismissing(t_mask_2d), w_e_spr@_FillValue, w_e_spr)
  w_e_sum = where(ismissing(t_mask_2d), w_e_sum@_FillValue, w_e_sum)
  
  ; 对 2 个差异场应用掩膜
  diff_spr = where(ismissing(t_mask_2d), diff_spr@_FillValue, diff_spr)
  diff_sum = where(ismissing(t_mask_2d), diff_sum@_FillValue, diff_sum)
  ; ======================================================================
  w_c_spr@lat2d = lat2d
  w_c_spr@lon2d = lon2d
  w_c_sum@lat2d = lat2d
  w_c_sum@lon2d = lon2d
  
  w_e_spr@lat2d = lat2d
  w_e_spr@lon2d = lon2d
  w_e_sum@lat2d = lat2d
  w_e_sum@lon2d = lon2d

  diff_spr@lat2d = lat2d
  diff_spr@lon2d = lon2d
  diff_sum@lat2d = lat2d
  diff_sum@lon2d = lon2d

  printVarSummary(w_c_spr) 
  printMinMax(w_c_spr,0)

  printVarSummary(w_c_spr) 
  printMinMax(w_c_spr,0)
  printVarSummary(diff_sum) 
  printMinMax(diff_sum,0)

; exit
  ; ========== 开始绘图 ==========
  plot = new(6,graphic)
  
  wks=gsn_open_wks("pdf","./figs/omega500_diff")
  
  res                         = True
  res@tiMainString            = ""
  res@gsnDraw                 = False
  res@gsnFrame                = False
  res@gsnAddCyclic            = False
  ; res@tfDoNDCOverlay          = True    

  res@tmXTOn                  = False  
  res@tmYROn                  = False

  res@mpLimitMode             = "Corners"
  res@mpLeftCornerLatF        = 48
  res@mpLeftCornerLonF        = 60
  res@mpRightCornerLatF       = 10
  res@mpRightCornerLonF       = 135

  res@mpProjection            = "LambertConformal"
  res@mpLambertParallel1F     = 30.
  res@mpLambertParallel2F     = 60.
  res@mpLambertMeridianF      = 110.
  ; res@mpCenterLatF            = 35.1778 
  ; res@mpCenterLonF            = 110.   

  res@mpFillOn                = False
  res@mpOutlineOn             = True
  res@mpDataBaseVersion       = "MediumRes"        
  res@pmTickMarkDisplayMode   = "Always"
  res@mpGeophysicalLineColor  = "grey60"
  res@mpGeophysicalLineThicknessF = 0.3  
  res@mpOutlineBoundarySets   = "GeophysicalAndUSStates"
  
  res@cnFillOn             = True
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnInfoLabelOn        = False
  res@lbLabelBarOn         = True      
  res@lbOrientation        = "Horizontal"
  res@pmLabelBarOrthogonalPosF = 0.1   
  
  ; --- 500hPa 垂直速度 基础场阴影设置 ---
  ; 提示：数值需根据 printMinMax 的结果进行调整，此处提供常用默认值
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -0.1
  res@cnMaxLevelValF       =  0.1
  res@cnLevelSpacingF      =  0.05
  res@cnFillPalette        = "MPL_RdBu"  ;BlueWhiteOrangeRed 蓝红渐变色表（负值上升蓝，正值下沉红）
  res@gsnRightString       = "Pa/s"

  ; --- 绘制 a-d 基础场 ---
  res@gsnLeftString = "(a)CTL MAM"
  plot(0) = gsn_csm_contour_map(wks,w_c_spr,res)

  res@gsnLeftString = "(b)CTL JJA"
  plot(1) = gsn_csm_contour_map(wks,w_c_sum,res)

  res@gsnLeftString = "(c)EXP MAM"
  plot(2) = gsn_csm_contour_map(wks,w_e_spr,res)

  res@gsnLeftString = "(d)EXP JJA"
  plot(3) = gsn_csm_contour_map(wks,w_e_sum,res)

  ; --- 500hPa 垂直速度 差异场阴影设置 ---
  res@lbLabelBarOn         =  True
  res@cnMinLevelValF       = -0.05  
  res@cnMaxLevelValF       =  0.05
  res@cnLevelSpacingF      =  0.01
  res@cnFillPalette        = "NCV_blu_red"
  
  ; (e) 差异场图叠加
  res@gsnLeftString = "(e)EXP MAM - CTL MAM"
  plot(4) = gsn_csm_contour_map(wks, diff_spr, res)

  ; (f) 差异场图叠加
  res@gsnLeftString = "(f)EXP JJA - CTL JJA"
  plot(5) = gsn_csm_contour_map(wks, diff_sum, res)

  ; --- 添加 Shapefile ---
  dir                    = "../../shapefile_China/"
  fShape1Name            = dir + "china.shp"
  fShape2Name            = dir + "south_china_sea.shp"
  fShape3Name            = dir + "river.nc"

  resLine                = True
  resLine@gsLineThicknessF    = 0.4
  
  line  = new(6, graphic)
  scs   = new(6, graphic)
  river = new(6, graphic)
  
  do i=0,5
    resLine@gsLineColor = "black"
    line(i) = gsn_add_shapefile_polylines(wks, plot(i), fShape1Name, resLine)
    scs(i)  = gsn_add_shapefile_polylines(wks, plot(i), fShape2Name, resLine)
    resLine@gsLineColor = "blue"
    river(i) = gsn_add_shapefile_polylines(wks, plot(i), fShape3Name, resLine)
  end do

  ; --- Panel 设置 ---
  resP = True
  resP@gsnPanelMainString    = " "
  resP@gsnPanelLabelBar      = False
  resP@lbBoxLinesOn          = False      
  resP@pmLabelBarWidthF      = 0.8        
  resP@pmLabelBarHeightF     = 0.06       
  resP@pmLabelBarOrthogonalPosF = -0.01   
  resP@lbBoxEndCapStyle      = "TriangleBothEnds"
  resP@lbTitleString         = "500 hPa Mean Vertical Velocity (Omega) Difference(2001-2017)" 
  resP@lbTitlePosition       = "Bottom"      
  resP@lbTitleFontHeightF    = 0.012        
  resP@lbTitleDirection      = "Across"     
  resP@lbTitleOffsetF        = 0.025         
  gsn_panel(wks, plot, (/3,2/), resP)

end
