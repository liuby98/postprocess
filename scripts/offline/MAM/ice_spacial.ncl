begin
  ; ========================================================
  ; --- 核心修改：控制开关 ---
  ; 选项: "raw" (输出3,4,5月的CTL和EXP原始分布) 
  ;       "diff" (输出3,4,5月的 EXP-CTL 差值分布)
  plot_type = "raw" 
  ; ========================================================

  path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
  f1=addfile(path+"colmoff_2001-2023_monmean_nogravel.nc","r")
  f2=addfile(path+"colmoff_2001-2023_monmean_gravel.nc","r")
  
  ; 读取 f_wice_soisno:"ice lens in soil layers" kg/m2
  t_ctl_mon = tofloat(f1->f_wice_soisno(:,:,:,5:14))
  t_exp_mon = tofloat(f2->f_wice_soisno(:,:,:,5:14))

  f_input = addfile("/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01","r")
  nlat = f_input->XLAT(0,:,0)
  nlon = f_input->XLONG(0,0,:)
  lat2d = f_input->XLAT(0,:,:)
  lon2d = f_input->XLONG(0,:,:)

  ; 提取 3,4,5 月数据
  ctl_mar = new((/17,342,462,10/),float)
  exp_mar = new((/17,342,462,10/),float)
  ctl_apr = new((/17,342,462,10/),float)
  exp_apr = new((/17,342,462,10/),float)
  ctl_may = new((/17,342,462,10/),float)
  exp_may = new((/17,342,462,10/),float)
  
  do i=0,16
     ctl_mar(i,:,:,:)=t_ctl_mon(6*i+3,:,:,:)
     exp_mar(i,:,:,:)=t_exp_mon(6*i+3,:,:,:)
     
     ctl_apr(i,:,:,:)=t_ctl_mon(6*i+4,:,:,:)
     exp_apr(i,:,:,:)=t_exp_mon(6*i+4,:,:,:)
     
     ctl_may(i,:,:,:)=t_ctl_mon(6*i+5,:,:,:)
     exp_may(i,:,:,:)=t_exp_mon(6*i+5,:,:,:)
  end do

  ; 创建 raw 数组保存 6 种原始状态
  raw = new((/6, 342, 462, 10/), float)
  raw(0,:,:,:) = dim_avg_n_Wrap(ctl_mar,0)
  raw(1,:,:,:) = dim_avg_n_Wrap(exp_mar,0)
  raw(2,:,:,:) = dim_avg_n_Wrap(ctl_apr,0)
  raw(3,:,:,:) = dim_avg_n_Wrap(exp_apr,0)
  raw(4,:,:,:) = dim_avg_n_Wrap(ctl_may,0)
  raw(5,:,:,:) = dim_avg_n_Wrap(exp_may,0)

  raw!0 = "config"
  raw!1 = "lat"
  raw!2 = "lon"
  raw!3 = "level"
  raw&lat = nlat
  raw&lon = nlon
  raw&lat@units = "degrees_north"
  raw&lon@units = "degrees_east"
  raw_re = raw(config|:, level|:, lat|:, lon|:)

  f=addfile("/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc","r")
  t_obs_mon = f->tm
  lat1d = f->lat
  lon1d = f->lon
  raw_re@lat2d = lat2d
  raw_re@lon2d = lon2d
  raw_ip = rcm2rgrid_Wrap(lat2d,lon2d,raw_re,lat1d,lon1d,1)
  
  mask_2d = t_obs_mon(0,:,:)
  mask_3d = conform(raw_ip, mask_2d, (/2,3/))
  raw_ip = where(ismissing(mask_3d), raw_ip@_FillValue, raw_ip)
  
  ; 过滤极小异常值 (注意这里不需要减 273.15，因为是 kg/m2 不是温度)
  raw_ip = where(raw_ip .lt. -1000., raw_ip@_FillValue, raw_ip)

  ; ========================================================
  ; --- 核心修改：根据开关处理数据、复制元数据和设置维度 ---
  dims = dimsizes(raw_ip)
  if (plot_type .eq. "diff") then
      plot_data = new((/3, dims(1), dims(2), dims(3)/), float)
      copy_VarMeta(raw_ip(0:2,:,:,:), plot_data) 
      
      plot_data(0,:,:,:) = raw_ip(1,:,:,:) - raw_ip(0,:,:,:)  ; 3月 EXP-CTL
      plot_data(1,:,:,:) = raw_ip(3,:,:,:) - raw_ip(2,:,:,:)  ; 4月 EXP-CTL
      plot_data(2,:,:,:) = raw_ip(5,:,:,:) - raw_ip(4,:,:,:)  ; 5月 EXP-CTL
      
      ; 原脚本中的极小值过滤 (仅对差值过滤)
      plot_data = where(abs(plot_data) .lt. 0.001, plot_data@_FillValue, plot_data)
      
      ncols = 3
      ntot  = 30
      out_name = "colmoff_JJA_ice_spacial_diff"
      lb_title = "Ice Lens Difference"
  else
      plot_data = raw_ip
      ncols = 6
      ntot  = 60
      out_name = "colmoff_JJA_ice_spacial_raw"
      lb_title = "Ice Lens Content"
  end if
  ; ========================================================

  printVarSummary(plot_data)
  printMinMax(plot_data, 0)

  wks = gsn_open_wks("pdf", out_name)
  plot= new(ntot,graphic)
  res                         = True
  res@gsnDraw                 = False
  res@gsnFrame                = False
  res@gsnAddCyclic            = False

  res@mpLimitMode            = "Corners"          
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
  res@mpGridAndLimbOn       = True
  res@mpGridSpacingF        = 10          
  res@mpGridLineDashPattern = 2           
  res@mpGridLineColor       = "grey80"    
  res@mpGridLineThicknessF  = 0.8

  res@mpDataBaseVersion      = "MediumRes"        
  res@pmTickMarkDisplayMode  = "Always"
  res@mpOutlineOn           = True                   
  res@mpGeophysicalLineColor  = "grey60"   
  res@mpGeophysicalLineThicknessF = 0.3
  res@mpOutlineDrawOrder    = "PreDraw"

  res@lbLabelBarOn            = False
  res@lbBoxLinesOn            = False
  res@lbOrientation           = "Horizontal"
  res@lbTitlePosition         = "Bottom"
  res@lbTitleDirection        = "Across"
  res@lbTitleAngleF           =  0
  res@lbTitleString           = " "

  ; ========================================================
  ; --- 核心修改：为原始图和差值图分别配置色板和范围 ---
  res@cnLevelSelectionMode  = "ManualLevels"
  
  if (plot_type .eq. "diff") then
      ; --- 保持你原来的紫蓝拼接色板 ---
      purple = read_colormap_file("MPL_Purples")     
      blue   = read_colormap_file("MPL_Blues")       
      purple_r = purple(::-1, :)
      neg_raw = purple_r(10::2, :)     
      pos_raw = blue(5::2, :)        

      n_neg = 42                     
      n_pos = 42                     
      n_neg := min((/n_neg, dimsizes(neg_raw(:,0))/))
      n_pos := min((/n_pos, dimsizes(pos_raw(:,0))/))
      neg_part = neg_raw(0:n_neg-1, :)
      pos_part = pos_raw(0:n_pos-1, :)
      ncols_pal = dimsizes(neg_part(0,:))          

      total_colors = n_neg + n_pos        
      spliced = new( (/total_colors, ncols_pal/), float )
      spliced(0:n_neg-1, :)          = neg_part
      spliced(n_neg:total_colors-1, :) = pos_part

      res@cnFillPalette           = spliced 
      res@cnMaxLevelValF          =  80.0  
      res@cnMinLevelValF          = -80.0  
      res@cnLevelSpacingF         =   1.0
  else
      ; --- 原始含量图：使用单向蓝色渐变 ---
      gsn_define_colormap(wks, "MPL_Blues")
      res@gsnSpreadColorStart = 2
      res@gsnSpreadColorEnd   = 128
      ; 注意：这里的原始值范围可能需要你根据 printMinMax 的结果进行调整
      res@cnMaxLevelValF          =  300.0  
      res@cnMinLevelValF          =  100.0  
      res@cnLevelSpacingF         =    2.0
  end if
  ; ========================================================

  res@gsnRightString          = "kg/m~S~2~N~" 
  res@tiYAxisFontHeightF      = 0.02
  res@tiYAxisAngleF           = 90
  res@cnLineLabelsOn        = False
  
  do j=0,9
     k = j * ncols 
     res@tiYAxisString = "n=" + (j+1)
     
     if (plot_type .eq. "diff") then
         res@gsnLeftString = "EXP-CTL Jun"
         plot(k)   = gsn_csm_contour_map(wks,plot_data(0,j,:,:),res)
         res@tiYAxisString = ""
         res@gsnLeftString = "EXP-CTL Jul"
         plot(k+1) = gsn_csm_contour_map(wks,plot_data(1,j,:,:),res)
         res@gsnLeftString = "EXP-CTL Aug"
         plot(k+2) = gsn_csm_contour_map(wks,plot_data(2,j,:,:),res)
     else
         res@gsnLeftString = "CTL Jun"
         plot(k)   = gsn_csm_contour_map(wks,plot_data(0,j,:,:),res)
         res@tiYAxisString = ""
         res@gsnLeftString = "EXP Jun"
         plot(k+1) = gsn_csm_contour_map(wks,plot_data(1,j,:,:),res)
         res@gsnLeftString = "CTL Jul"
         plot(k+2) = gsn_csm_contour_map(wks,plot_data(2,j,:,:),res)
         res@gsnLeftString = "EXP Jul"
         plot(k+3) = gsn_csm_contour_map(wks,plot_data(3,j,:,:),res)
         res@gsnLeftString = "CTL Aug"
         plot(k+4) = gsn_csm_contour_map(wks,plot_data(4,j,:,:),res)
         res@gsnLeftString = "EXP Aug"
         plot(k+5) = gsn_csm_contour_map(wks,plot_data(5,j,:,:),res)
     end if
  end do

  ; 添加文本 (可选，这里仅在 diff 模式下添加，且位置可能需要根据 3 列排版微调)
  if (plot_type .eq. "diff") then
      txres = True
      txres@txFontHeightF = 0.02
      txres@txFontColor = "black"
      txt = new(30, graphic)
      do p=0,29
          txt(p) = gsn_add_text(wks, plot(p), "EXP-CTL", lon1d({127}), lat1d({17}), txres)
      end do
  end if

  dir                    = "../../../shapefile_China/"
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
  resP@gsnPanelBottom        = 0.08 
  resP@gsnPanelLabelBar      = True
  resP@lbBoxLinesOn          = False      
  resP@pmLabelBarWidthF      = 0.8        
  resP@pmLabelBarHeightF     = 0.06       
  resP@pmLabelBarOrthogonalPosF = -0.01   
  resP@lbBoxEndCapStyle         = "TriangleBothEnds"
  resP@lbTitleString         = lb_title
  resP@lbTitlePosition       = "Bottom"      
  resP@lbTitleFontHeightF    = 0.012        
  resP@lbTitleDirection      = "Across"     
  resP@lbTitleOffsetF        = 0.025         
  resP@gsnPanelYWhiteSpace   = 5.0  
  resP@gsnPanelXWhiteSpace   = 5.0  

  ; ========================================================
  ; --- 核心修改：根据模式自适应 panel 排版和页数 ---
  if (plot_type .eq. "diff") then
      gsn_panel(wks, plot(0:14), (/5,3/), resP)
      gsn_panel(wks, plot(15:29), (/5,3/), resP)
  else
      gsn_panel(wks, plot(0:29), (/5,6/), resP)
      gsn_panel(wks, plot(30:59), (/5,6/), resP)
  end if
  ; ========================================================
end