begin
  path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
  f1=addfile(path+"wrfpost_gh500_2001_2017_MAM_nogravel.nc","r")
  f2=addfile(path+"wrfpost_gh500_2001_2017_JJA_nogravel.nc","r")
  f3=addfile(path+"wrfpost_gh500_2001_2017_MAM_gravel.nc","r")
  f4=addfile(path+"wrfpost_gh500_2001_2017_JJA_gravel.nc","r")

  f_input = addfile("/share/home/dq135/wrfpost/wrfinput_d01", "r")
  lon2d = f_input->XLONG(0,:,:)
  lat2d = f_input->XLAT(0,:,:)

  ; 读取 geopt 并直接换算为位势米 (gpm)
  h_ctl_mam = f1->geopt(:,0,:,:) / 9.8
  h_ctl_jja = f2->geopt(:,0,:,:) / 9.8
  h_exp_mam = f3->geopt(:,0,:,:) / 9.8
  h_exp_jja = f4->geopt(:,0,:,:) / 9.8

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
  printVarSummary(h_c_sum) 
  printMinMax(h_c_sum,0)
  printVarSummary(diff_sum) 
  printMinMax(diff_sum,0)
; exit
  ; ========== 开始绘图 ==========
  plot = new(6,graphic)
  plot_wpsh = new(6,graphic)  ; 用于存放脊线图层
  
  wks=gsn_open_wks("pdf","./figs/hgt500_diff")
  
  res                         = True
  res@tiMainString            = ""
  res@gsnDraw                 = False
  res@gsnFrame                = False
  res@gsnAddCyclic            = False
  res@tfDoNDCOverlay          = True    

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
  res@mpCenterLatF            = 35.1778 
  res@mpCenterLonF            = 110.   

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
  
  ; --- 500hPa 阴影场设置 ---
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       =  5300.  ;1330.  
  res@cnMaxLevelValF       =  5900.  ;1580. 
  res@cnLevelSpacingF      =  10.     ;5. 
  res@cnFillPalette        = "BlAqGrYeOrReVi200" 
  res@gsnRightString       = "gpm"

  ; --- 副高脊线 (5880 gpm) 属性设置 ---
  res_line = True
  res_line@gsnDraw = False
  res_line@gsnFrame = False
  res_line@cnFillOn = False
  res_line@cnLinesOn = True
  res_line@cnLineLabelsOn = True  ; 开启等值线标签显示 "5880"
  res_line@cnLineLabelBackgroundColor = -1 ; 标签背景透明
  res_line@cnInfoLabelOn = False
  res_line@cnLevelSelectionMode = "ExplicitLevels"
  res_line@cnLevels = (/ 5880. /) ; 选取 5880 脊线
  res_line@cnLineColor = "greenyellow"  ; 黄绿色
  res_line@cnLineThicknessF = 0.6
  res_line@gsnLeftString = ""
  res_line@gsnRightString = ""

  ; --- 绘制 a-d 基础场及叠加其脊线 ---
  res@gsnLeftString = "(a)CTL MAM"
  plot(0) = gsn_csm_contour_map(wks,h_c_spr,res)
  ; plot_wpsh(0) = gsn_csm_contour(wks,h_c_spr,res_line)
  ; overlay(plot(0), plot_wpsh(0))

  res@gsnLeftString = "(b)CTL JJA"
  plot(1) = gsn_csm_contour_map(wks,h_c_sum,res)
  ; plot_wpsh(1) = gsn_csm_contour(wks,h_c_sum,res_line)
  ; overlay(plot(1), plot_wpsh(1))

  res@gsnLeftString = "(c)EXP MAM"
  plot(2) = gsn_csm_contour_map(wks,h_e_spr,res)
  ; plot_wpsh(2) = gsn_csm_contour(wks,h_e_spr,res_line)
  ; overlay(plot(2), plot_wpsh(2))

  res@gsnLeftString = "(d)EXP JJA"
  plot(3) = gsn_csm_contour_map(wks,h_e_sum,res)
  ; plot_wpsh(3) = gsn_csm_contour(wks,h_e_sum,res_line)
  ; overlay(plot(3), plot_wpsh(3))

  ; --- 500hPa 差异场图设置 ---
  res@lbLabelBarOn         =  True
  res@cnMinLevelValF       = -1.8  
  res@cnMaxLevelValF       =  1.8
  res@cnLevelSpacingF      =  0.1
  res@cnFillPalette        = "NCV_blu_red"
  
  ; (e) 差异场图叠加：CTL 为黑色对比线，EXP 为黄绿色目标线
  res@gsnLeftString = "(e)EXP MAM - CTL MAM"
  plot(4) = gsn_csm_contour_map(wks, diff_spr, res)
  res_line@cnLineColor = "black"
  res_line@cnLineThicknessF = 0.5
  plot_wpsh_ctl_spr = gsn_csm_contour(wks, h_c_spr, res_line)
  res_line@cnLineColor = "greenyellow"
  res_line@cnLineThicknessF = 0.5
  plot_wpsh_exp_spr = gsn_csm_contour(wks, h_e_spr, res_line)
  ; overlay(plot(4), plot_wpsh_ctl_spr)
  ; overlay(plot(4), plot_wpsh_exp_spr)

  ; (f) 差异场图叠加
  res@gsnLeftString = "(f)EXP JJA - CTL JJA"
  plot(5) = gsn_csm_contour_map(wks, diff_sum, res)
  res_line@cnLineColor = "black"
  res_line@cnLineThicknessF = 0.5
  plot_wpsh_ctl_sum = gsn_csm_contour(wks, h_c_sum, res_line)
  res_line@cnLineColor = "greenyellow"
  res_line@cnLineThicknessF = 0.5
  plot_wpsh_exp_sum = gsn_csm_contour(wks, h_e_sum, res_line)
  ; overlay(plot(5), plot_wpsh_ctl_sum)
  ; overlay(plot(5), plot_wpsh_exp_sum)

  ; --- 添加 Shapefile ---
  dir                    = "../../shapefile_China/"
  fShape1Name            = dir + "china.shp"
  fShape2Name            = dir + "south_china_sea.shp"
  fShape3Name            = dir + "river.nc"

  resLine                     = True
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
  resP@lbBoxEndCapStyle         = "TriangleBothEnds"
  resP@lbTitleString         = "500 hPa Mean Geopotential Height Field Difference(2001-2017)" 
  resP@lbTitlePosition       = "Bottom"      
  resP@lbTitleFontHeightF    = 0.012        
  resP@lbTitleDirection      = "Across"     
  resP@lbTitleOffsetF        = 0.025         
  gsn_panel(wks, plot, (/3,2/), resP)

end
