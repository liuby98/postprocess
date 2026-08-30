begin
  path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
  
  f1=addfile(path+"wrfpost_qx850_2001_2017_MAM_nogravel.nc","r")
  f2=addfile(path+"wrfpost_qx850_2001_2017_JJA_nogravel.nc","r")
  f3=addfile(path+"wrfpost_qx850_2001_2017_MAM_gravel.nc","r")
  f4=addfile(path+"wrfpost_qx850_2001_2017_JJA_gravel.nc","r")

  lon2d = f1->lon
  lat2d = f1->lat

  uq_ctl_mam=f1->uq(:,0,:,:)
  vq_ctl_mam=f1->vq(:,0,:,:)
  uq_ctl_jja=f2->uq(:,0,:,:)
  vq_ctl_jja=f2->vq(:,0,:,:)

  uq_exp_mam=f3->uq(:,0,:,:)
  vq_exp_mam=f3->vq(:,0,:,:)
  uq_exp_jja=f4->uq(:,0,:,:)
  vq_exp_jja=f4->vq(:,0,:,:)

  uq_ctl_spr = new((/17,342,462/),float)
  uq_exp_spr = new((/17,342,462/),float)
  uq_ctl_sum = new((/17,342,462/),float)
  uq_exp_sum = new((/17,342,462/),float)

  vq_ctl_spr = new((/17,342,462/),float)
  vq_exp_spr = new((/17,342,462/),float)
  vq_ctl_sum = new((/17,342,462/),float)
  vq_exp_sum = new((/17,342,462/),float)

 do i=0,16
   uq_ctl_spr(i,:,:)=dim_avg_n_Wrap(uq_ctl_mam(92*i:92*i+91,:,:),0)
   uq_exp_spr(i,:,:)=dim_avg_n_Wrap(uq_exp_mam(92*i:92*i+91,:,:),0)
   uq_ctl_sum(i,:,:)=dim_avg_n_Wrap(uq_ctl_jja(92*i:92*i+91,:,:),0)
   uq_exp_sum(i,:,:)=dim_avg_n_Wrap(uq_exp_jja(92*i:92*i+91,:,:),0)

   vq_ctl_spr(i,:,:)=dim_avg_n_Wrap(vq_ctl_mam(92*i:92*i+91,:,:),0)
   vq_exp_spr(i,:,:)=dim_avg_n_Wrap(vq_exp_mam(92*i:92*i+91,:,:),0)
   vq_ctl_sum(i,:,:)=dim_avg_n_Wrap(vq_ctl_jja(92*i:92*i+91,:,:),0)
   vq_exp_sum(i,:,:)=dim_avg_n_Wrap(vq_exp_jja(92*i:92*i+91,:,:),0)
 end do

  scale = 1000.0

  uq_c_spr = dim_avg_n(uq_ctl_spr,0) * scale
  uq_e_spr = dim_avg_n(uq_exp_spr,0) * scale
  uq_c_sum = dim_avg_n(uq_ctl_sum,0) * scale
  uq_e_sum = dim_avg_n(uq_exp_sum,0) * scale

  vq_c_spr = dim_avg_n(vq_ctl_spr,0) * scale
  vq_e_spr = dim_avg_n(vq_exp_spr,0) * scale
  vq_c_sum = dim_avg_n(vq_ctl_sum,0) * scale
  vq_e_sum = dim_avg_n(vq_exp_sum,0) * scale

  uq_c_spr@lat2d = lat2d
  uq_e_spr@lat2d = lat2d
  uq_c_sum@lat2d = lat2d
  uq_e_sum@lat2d = lat2d
  vq_c_spr@lat2d = lat2d
  vq_e_spr@lat2d = lat2d
  vq_c_sum@lat2d = lat2d
  vq_e_sum@lat2d = lat2d
  uq_c_spr@lon2d = lon2d
  uq_e_spr@lon2d = lon2d
  uq_c_sum@lon2d = lon2d
  uq_e_sum@lon2d = lon2d
  vq_c_spr@lon2d = lon2d
  vq_e_spr@lon2d = lon2d
  vq_c_sum@lon2d = lon2d
  vq_e_sum@lon2d = lon2d


  spd_c_spr = sqrt(uq_c_spr^2 + vq_c_spr^2)
  spd_e_spr = sqrt(uq_e_spr^2 + vq_e_spr^2)
  diff_spr  = spd_e_spr - spd_c_spr
  copy_VarCoords(uq_c_spr, diff_spr)

  spd_c_sum = sqrt(uq_c_sum^2 + vq_c_sum^2)
  spd_e_sum = sqrt(uq_e_sum^2 + vq_e_sum^2)
  diff_sum  = spd_e_sum - spd_c_sum
  copy_VarCoords(uq_c_sum, diff_sum)

; ========================= 绘图部分 =========================
  plot = new(6,graphic)
  wks=gsn_open_wks("pdf","./figs/qflux850_diff")
  
  res                   = True
  res@tiMainString      = ""
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@gsnAddCyclic      = False
  
  res@tfDoNDCOverlay    = True     

  res@tmXTOn            = False  
  res@tmYROn            = False

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
  res@mpDataBaseVersion       = "MediumRes"        
  res@pmTickMarkDisplayMode   = "Always"  
  res@mpGeophysicalLineColor  = "grey60"   
  res@mpGeophysicalLineThicknessF = 0.3  
  res@mpOutlineBoundarySets   = "GeophysicalAndUSStates"  

  res_diff = res

  ; === a-d面板：矢量箭头参数 ===
  res@vcGlyphStyle              = "CurlyVector"   
  res@vcMonoLineArrowColor      = False           ; 让矢量自动根据自身大小(模长)赋予渐变色
  res@vcLevelPalette            = "GMT_drywet"    
  res@lbLabelBarOn              = False           
  
  res@vcLevelSelectionMode      = "ManualLevels"
  res@vcMinLevelValF            = 20.0            
  res@vcMaxLevelValF            = 200.0           
  res@vcLevelSpacingF           = 20.0            

  res@vcRefMagnitudeF           = 50.0            
  res@vcRefLengthF              = 0.05            
  res@vcRefAnnoOrthogonalPosF   = -1.085           
  res@vcRefAnnoParallelPosF     =  1.0            
  ; res@vcRefAnnoPerimOn          = False           
  ; res@vcRefAnnoBackgroundColor  = "Transparent"   
  res@vcRefAnnoArrowLineColor   = "red"           
  res@vcRefAnnoArrowUseVecColor = False           
  res@vcMinDistanceF            = 0.02            
  res@vcLineArrowThicknessF     = 2.0             
  res@vcLineArrowHeadMaxSizeF   = 0.01
  res@vcRefAnnoString2On        = False
  
  res@gsnRightString            = "g kg~S~-1~N~ m s~S~-1~N~"
  res@gsnLeftString           = "(a)CTL MAM"
  plot(0) = gsn_csm_vector_map(wks,uq_c_spr,vq_c_spr,res)
  res@gsnLeftString           = "(b)CTL JJA"
  plot(1) = gsn_csm_vector_map(wks,uq_c_sum,vq_c_sum,res)
  res@gsnLeftString           = "(c)EXP MAM"
  plot(2) = gsn_csm_vector_map(wks,uq_e_spr,vq_e_spr,res)
  res@gsnLeftString           = "(d)EXP JJA"
  plot(3) = gsn_csm_vector_map(wks,uq_e_sum,vq_e_sum,res)

  ; === e-f面板：差异填色图 ===
  res_diff@cnFillOn             = True
  res_diff@cnLinesOn            = False
  res_diff@cnLineLabelsOn       = False
  res_diff@cnInfoLabelOn        = False
  res_diff@lbLabelBarOn         = False           
  
  res_diff@cnMissingValFillColor= -1              
  res_diff@cnFillDrawOrder      = "PreDraw"       
  res_diff@cnFillPalette        = "MPL_BrBG"      
  
  res_diff@cnLevelSelectionMode = "ManualLevels"
  res_diff@cnMinLevelValF       = -5.0            
  res_diff@cnMaxLevelValF       =  5.0            
  res_diff@cnLevelSpacingF      =  0.5            

  res_diff@gsnRightString       = "g kg~S~-1~N~ m s~S~-1~N~"
  
  res_diff@gsnLeftString        = "(e)EXP - CTL (MAM)"
  plot(4) = gsn_csm_contour_map(wks, diff_spr, res_diff)

  res_diff@gsnLeftString        = "(f)EXP - CTL (JJA)"
  plot(5) = gsn_csm_contour_map(wks, diff_sum, res_diff)

  ; === 添加边界及河流 ===
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


  ; === 分割面板绘制 ===
  resP1 = True
  resP1@gsnFrame             = False        
  resP1@gsnPanelMainString   = "850 hPa Mean Moisture Flux Difference (2001-2017)"
  
  resP1@gsnPanelLabelBar     = True         
  resP1@pmLabelBarWidthF     = 0.6        
  resP1@pmLabelBarHeightF    = 0.05
  resP1@gsnPanelTop          = 0.95         
  resP1@gsnPanelBottom       = 0.36         
  
  gsn_panel(wks, plot(0:3), (/2,2/), resP1)

  resP2 = True
  resP2@gsnFrame             = True         
  
  resP2@gsnPanelLabelBar     = True         
  resP2@pmLabelBarWidthF     = 0.6        
  resP2@pmLabelBarHeightF    = 0.05
  resP2@gsnPanelTop          = 0.34         
  resP2@gsnPanelBottom       = 0.02         
  
  gsn_panel(wks, plot(4:5), (/1,2/), resP2)

end
