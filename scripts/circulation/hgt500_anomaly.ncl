begin
  path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
  ; 【修改点1】: 文件名由 gh200 改为 gh500
  f1=addfile(path+"wrfpost_gh500_2001_2017_MAM_nogravel.nc","r")
  f2=addfile(path+"wrfpost_gh500_2001_2017_JJA_nogravel.nc","r")
  f3=addfile(path+"wrfpost_gh500_2001_2017_MAM_gravel.nc","r")
  f4=addfile(path+"wrfpost_gh500_2001_2017_JJA_gravel.nc","r")

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

  ; 1. 计算气候态平均
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

  ; ==========================================================
  ; 2. 计算纬向距平异常 (去掉纬向平均)
  ; ==========================================================
  h_c_spr_ano = dim_rmvmean_n_Wrap(h_c_spr, 1)
  h_e_spr_ano = dim_rmvmean_n_Wrap(h_e_spr, 1)
  h_c_sum_ano = dim_rmvmean_n_Wrap(h_c_sum, 1)
  h_e_sum_ano = dim_rmvmean_n_Wrap(h_e_sum, 1)

  ; ==========================================================
  ; 3. 计算距平异常的差值 (EXP - CTL)
  ; ==========================================================
  diff_spr_ano  = h_e_spr_ano - h_c_spr_ano
  copy_VarCoords(h_c_spr_ano, diff_spr_ano)

  diff_sum_ano  = h_e_sum_ano - h_c_sum_ano
  copy_VarCoords(h_c_sum_ano, diff_sum_ano)

  printVarSummary(h_c_spr_ano) 
  printMinMax(h_c_spr_ano,0)
  printVarSummary(diff_sum_ano) 
  printMinMax(diff_sum_ano,0)

  ; 【修改点2】: 修改输出文件名
  plot = new(6,graphic)
  wks=gsn_open_wks("pdf","./figs/hgt500_anomaly_diff")
  
  res                   = True
  res@tiMainString      = ""
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@gsnAddCyclic      = False
  res@tfDoNDCOverlay    = True    

  res@tmXTOn                 = False  
  res@tmYROn                 = False  

  res@mpLimitMode            = "Corners"          
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
  res@mpOutlineOn = True   
  res@mpGeophysicalLineColor  = "grey60"   
  res@mpGeophysicalLineThicknessF = 0.3  
  res@mpOutlineBoundarySets = "GeophysicalAndUSStates"  
  
  res@cnFillOn             = True
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnInfoLabelOn        = False
  res@lbLabelBarOn         = True      
  res@lbOrientation        = "Horizontal"
  res@pmLabelBarOrthogonalPosF = 0.1   
  
  ; ==========================================================
  ; 【修改点3】: 调整500 hPa距平场的等值线范围
  ; 500hPa的高度场振幅较小，此处设为 -100 到 100
  ; ==========================================================
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -100. 
  res@cnMaxLevelValF       =  100.
  res@cnLevelSpacingF      =  10.
  res@cnFillPalette        = "NCV_blu_red" 
  res@gsnRightString       = "gpm"

  res@gsnLeftString           = "(a)CTL MAM Anomaly"
  plot(0) = gsn_csm_contour_map(wks, h_c_spr_ano, res)
  res@gsnLeftString           = "(b)CTL JJA Anomaly"
  plot(1) = gsn_csm_contour_map(wks, h_c_sum_ano, res)
  res@gsnLeftString           = "(c)EXP MAM Anomaly"
  plot(2) = gsn_csm_contour_map(wks, h_e_spr_ano, res)
  res@gsnLeftString           = "(d)EXP JJA Anomaly"
  plot(3) = gsn_csm_contour_map(wks, h_e_sum_ano, res)

  ; ==========================================================
  ; 5. 调整差值图的等值线范围 (EXP - CTL)
  ; 暂时保留原先的设置，若500hPa差异更小/大，可调整这里
  ; ==========================================================
;   res@lbLabelBarOn         = False
  res@cnMinLevelValF       = -1.2  
  res@cnMaxLevelValF       =  1.2
  res@cnLevelSpacingF      =  0.1
  res@cnFillPalette        = "NCV_blu_red"
  
  res@gsnLeftString           = "(e)EXP - CTL MAM Anomaly"
  plot(4) = gsn_csm_contour_map(wks, diff_spr_ano, res)
  res@gsnLeftString           = "(f)EXP - CTL JJA Anomaly"
  plot(5) = gsn_csm_contour_map(wks, diff_sum_ano, res)

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
  resP@gsnPanelMainString    = " "
  resP@gsnPanelLabelBar      = False
  resP@lbBoxLinesOn          = False      
  resP@pmLabelBarWidthF      = 0.8        
  resP@pmLabelBarHeightF     = 0.06       
  resP@pmLabelBarOrthogonalPosF = -0.01   
  resP@lbBoxEndCapStyle         = "TriangleBothEnds"
  
  ; 【修改点4】: 修改单位标题为 500 hPa
  resP@lbTitleString         = "500 hPa Zonal Anomaly Difference (2001-2017)" 
  resP@lbTitlePosition       = "Bottom"      
  resP@lbTitleFontHeightF    = 0.012        
  resP@lbTitleDirection      = "Across"     
  resP@lbTitleOffsetF        = 0.025         
  
  gsn_panel(wks, plot, (/3,2/), resP)

end
