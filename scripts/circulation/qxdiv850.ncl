begin
  path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
  
  ; 1. 读取水汽通量文件
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

  ; 先将量级统一放大到 g/kg
  scale = 1000.0

  uq_c_spr = dim_avg_n(uq_ctl_spr,0) * scale
  uq_e_spr = dim_avg_n(uq_exp_spr,0) * scale
  uq_c_sum = dim_avg_n(uq_ctl_sum,0) * scale
  uq_e_sum = dim_avg_n(uq_exp_sum,0) * scale

  vq_c_spr = dim_avg_n(vq_ctl_spr,0) * scale
  vq_e_spr = dim_avg_n(vq_exp_spr,0) * scale
  vq_c_sum = dim_avg_n(vq_ctl_sum,0) * scale
  vq_e_sum = dim_avg_n(vq_exp_sum,0) * scale

  ; ==========================================================
  ; === 核心新增：计算水汽通量散度 ===
  ; ==========================================================
  ; 提取区域中心的一维经纬度用于 CFD 中心差分计算，确保坐标单调递增
  lat1d = lat2d(:, 230)  
  lon1d = lon2d(170, :)
  
  div_scale = 100000.0   ; 放大 10^5 倍以便于画图展示
  
  div_c_spr = uv2dv_cfd(uq_c_spr, vq_c_spr, lat1d, lon1d, 2) * div_scale
  div_c_sum = uv2dv_cfd(uq_c_sum, vq_c_sum, lat1d, lon1d, 2) * div_scale
  div_e_spr = uv2dv_cfd(uq_e_spr, vq_e_spr, lat1d, lon1d, 2) * div_scale
  div_e_sum = uv2dv_cfd(uq_e_sum, vq_e_sum, lat1d, lon1d, 2) * div_scale

  ; 补齐绘图所需的二维经纬度坐标属性
  div_c_spr@lat2d = lat2d
  div_c_spr@lon2d = lon2d
  div_e_spr@lat2d = lat2d
  div_e_spr@lon2d = lon2d
  div_c_sum@lat2d = lat2d
  div_c_sum@lon2d = lon2d
  div_e_sum@lat2d = lat2d
  div_e_sum@lon2d = lon2d

  ; 计算散度差异 (EXP - CTL)
  diff_spr = div_e_spr - div_c_spr
  diff_sum = div_e_sum - div_c_sum
  
  diff_spr@lat2d = lat2d
  diff_spr@lon2d = lon2d
  diff_sum@lat2d = lat2d
  diff_sum@lon2d = lon2d

; ========================= 绘图部分 =========================
  plot = new(6,graphic)
  wks=gsn_open_wks("pdf","./figs/qdiv850_diff")  ; 保存为散度图的pdf
  
  res                   = True
  res@tiMainString      = ""
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@gsnAddCyclic      = False
  
  ; 继续使用原生贴图模式确保与 WRF 弯曲网格完全贴合
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

  ; === a-d面板：气候态散度填色图参数 ===
  res_div                       = res
  res_div@cnFillOn              = True
  res_div@cnLinesOn             = False
  res_div@cnLineLabelsOn        = False
  res_div@cnInfoLabelOn         = False
  res_div@lbLabelBarOn          = False           
  res_div@cnMissingValFillColor = -1              
  res_div@cnFillDrawOrder       = "PreDraw"       
  
  ; 读取并翻转色表：负值(辐合)为绿，正值(辐散)为棕
  cmap = read_colormap_file("MPL_BrBG")
  res_div@cnFillPalette         = cmap(::-1,:)
  
  res_div@cnLevelSelectionMode  = "ManualLevels"
  ; 建议根据实际散度极大极小值微调以下范围
  res_div@cnMinLevelValF        = -15.0            
  res_div@cnMaxLevelValF        =  15.0            
  res_div@cnLevelSpacingF       =   3.0            

  res_div@gsnRightString        = "10~S~-5~N~ g kg~S~-1~N~ s~S~-1~N~ m~S~-1~N~"
  
  res_div@gsnLeftString           = "(a)CTL MAM"
  plot(0) = gsn_csm_contour_map(wks, div_c_spr, res_div)
  res_div@gsnLeftString           = "(b)CTL JJA"
  plot(1) = gsn_csm_contour_map(wks, div_c_sum, res_div)
  res_div@gsnLeftString           = "(c)EXP MAM"
  plot(2) = gsn_csm_contour_map(wks, div_e_spr, res_div)
  res_div@gsnLeftString           = "(d)EXP JJA"
  plot(3) = gsn_csm_contour_map(wks, div_e_sum, res_div)

  ; === e-f面板：散度差异图参数 ===
  res_diff                      = res_div
  res_diff@cnMinLevelValF       = -5.0            
  res_diff@cnMaxLevelValF       =  5.0            
  res_diff@cnLevelSpacingF      =  0.5            

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
  resP1@gsnPanelMainString   = "850 hPa Mean Moisture Flux Divergence Difference(2001-2017)"
  
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
