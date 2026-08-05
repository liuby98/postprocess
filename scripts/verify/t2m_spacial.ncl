begin
  path = "/share/home/dq135/draw/data/CN15_graveltests_old_bk/"
  CASE=(/"nogravel","gravel"/)
  legend = (/"CTL","EXP"/)

  diff = new((/2,163,283/),float)
  data = new((/2,3/),float)

  f=addfile("/share/home/dq135/reference/CN05.1_Tm_1991_2023_JJA_025x025.nc","r")
  t_obs_mon = f->tm(30::,:,:)
  temp = new((/21,163,283/),float)
  do i=0,20
     temp(i,:,:) = dim_avg_n_Wrap(t_obs_mon(3*i:3*i+2,:,:),0)
  end do
  t2m_obs = dim_avg_n_Wrap(temp,0)

  f1 = addfile("/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01","r")
  lat2d = f1->XLAT(0,:,:)
  lon2d = f1->XLONG(0,:,:)

  lat1d = f->lat
  lon1d = f->lon

 do k=0,1

  f=addfile(path+"wrfout_1991-2023_JJA_daily_"+CASE(k)+"_t2m.nc","r")
  t_sim_mon = f->AT2M(920::,:,:)

     sim = new((/21,342,462/),float)
     do i=0,20
        sim(i,:,:) = dim_avg_n_Wrap(t_sim_mon(92*i:92*i+91,:,:),0)
     end do
     t2m = dim_avg_n_Wrap(sim,0)
     t2m = (/t2m-273.15/)


  t2m@lat2d = lat2d
  t2m@lon2d = lon2d
  t2m_ip    = rcm2rgrid_Wrap(lat2d,lon2d,t2m,lat1d,lon1d,1)
  ; t2m_ip({25:29.5},{91:98}) = t2m_ip@_FillValue
  t2m_ip    = where(ismissing(t2m_obs),t2m_ip@_FillValue,t2m_ip)
   
  ;printVarSummary(t2m_ip)
  printMinMax(t2m_ip,0)

  ; t2m_obs({25:29.5},{91:98}) = t2m_obs@_FillValue

  wgt = tofloat(cos(0.01745329*lat1d))
  tay_stat = taylor_stats(t2m_ip, t2m_obs, wgt, 1)

  data(k,0) = tay_stat(0)
  data(k,1) = tay_stat(1)
  data(k,2) = tay_stat(7)

  diff(k,:,:) = t2m_ip
  copy_VarCoords(temp(0,:,:), diff(k,:,:))

  delete(t_sim_mon)
  delete(t2m_ip)
end do
  taylor = taylor_stats(t2m_obs, t2m_obs, wgt, 1)
  m = taylor(1)
  n = taylor(7)
; exit
  wks = gsn_open_wks("pdf","wrfout_t2m_JJA_spacial")
  plot = new(3,graphic)
do k=0,1
  res                         = True
  res@gsnDraw                 = False
  res@gsnFrame                = False
  res@vpWidthF                = 0.7

  res@gsnAddCyclic            = False
  res@cnFillOn                = True
  res@cnLinesOn               = False
  res@mpFillOn                = False

  res@lbBoxLinesOn            = False
  res@lbOrientation           = "vertical"
  res@lbLabelAutoStride       = False
  ;res@lbLabelStride           = 1
  res@lbBoxEndCapStyle        = "TriangleBothEnds"
  res@lbTitleString           =  "~S~o~N~C"
  res@lbTitlePosition         = "Right"
  res@lbTitleDirection        = "Across"
  res@lbTitleAngleF           =   90

  res@mpLimitMode            = "Corners"
  res@mpLeftCornerLatF       = 47
  res@mpLeftCornerLonF       = 67
  res@mpRightCornerLatF      = 13
  res@mpRightCornerLonF      = 130
  res@tmXTOn                 = False
  res@tmYROn                 = False

  res@mpProjection           = "LambertConformal"
  res@mpLambertParallel1F    = 30.
  res@mpLambertParallel2F    = 60.
  res@mpLambertMeridianF     = 110.
  res@mpGridAndLimbOn        = True
  res@mpGridSpacingF        = 10          ; 可选 5、10、15
  res@mpGridLineDashPattern = 2           ; 2=短虚线（最好看）
  res@mpGridLineColor       = "grey80"    ; 灰色
  res@mpGridLineThicknessF  = 1.0

  res@mpDataBaseVersion      = "MediumRes"        
  res@pmTickMarkDisplayMode  = "Always"      
  res@mpGeophysicalLineColor  = "grey30"  
  res@mpGeophysicalLineThicknessF = 0.05 
  res@mpOutlineDrawOrder    = "PreDraw"
  ; res@mpDataSetName              = "../../CASES1991-2021/ncl_cnmap/Earth..4"
  ;res@gsnMaximize            = True

  res@gsnLeftString           = ""  
  res@gsnRightString          = ""

  res@cnFillPalette          = "WhiteBlueGreenYellowRed"
  res@cnLevelSelectionMode    = "ManualLevels"
  res@cnMaxLevelValF          =  25.0
  res@cnMinLevelValF          = -10.0
  res@cnLevelSpacingF         =  1.0
  plot(k+1) = gsn_csm_contour_map(wks,diff(k,:,:),res)
  txres = True
  txres@txFontHeightF = 0.02  
  txres@txFontColor = "red"   
  corplot    =  gsn_add_text(wks, plot(k+1), "CC:"+sprintf("%5.2f", data(k,0)), lon1d({92}),lat1d({21}), txres)
  stdplot    =  gsn_add_text(wks, plot(k+1), "STD ratio:"+sprintf("%5.2f", data(k,1)), lon1d({90}),lat1d({17.8}), txres)
  rmseplot   =  gsn_add_text(wks, plot(k+1), "RMSE:"+sprintf("%5.2f", data(k,2)),lon1d({91.5}),lat1d({15}), txres)
  txres@txFontColor = "Black"
  ctl        =  gsn_add_text(wks, plot(k+1), legend(k),lon1d({127}),lat1d({17}), txres)
end do

res2                         = res
  res2@gsnLeftString           = ""  
  res2@gsnRightString          = ""
  plot(0)     = gsn_csm_contour_map(wks,t2m_obs,res2)

  txres@txFontColor = "Black" 
  sim_diff1    = gsn_add_text(wks, plot(0), "OBS",lon1d({127}),lat1d({17}), txres)

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


  resP                       = True
  gsn_panel(wks,plot, (/2,3/), resP)
  draw(plot)

end

