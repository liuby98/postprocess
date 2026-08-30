begin
  path = "/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/"
  f1=addfile(path+"colmoff_2001-2023_monmean_nogravel.nc","r")
  f2=addfile(path+"colmoff_2001-2023_monmean_gravel.nc","r")
  ; f_t_soisno(:,:,:,5:14):"soil temperature"
  ; f_h2osoi:"volumetric water in soil layers"
  ; f_wliq_soisno(:,:,:,5:14):"liquid water in soil layers" kg/m2
  ; f_wice_soisno(:,:,:,5:14):"ice lens in soil layers" kg/m2

  t_ctl_mon = tofloat(f1->f_wice_soisno(:,:,:,5:14))
  t_exp_mon = tofloat(f2->f_wice_soisno(:,:,:,5:14))

  f_input = addfile("/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01","r")
  nlat = f_input->XLAT(0,:,0)
  nlon = f_input->XLONG(0,0,:)
  lat2d = f_input->XLAT(0,:,:)
  lon2d = f_input->XLONG(0,:,:)

  ctl_spr = new((/17,342,462,10/),float)
  exp_spr = new((/17,342,462,10/),float)
  ctl_sum = new((/17,342,462,10/),float)
  exp_sum = new((/17,342,462,10/),float)
 do i=0,16
   ctl_spr(i,:,:,:)=dim_avg_n_Wrap(t_ctl_mon(6*i:6*i+2,:,:,:),0)
   exp_spr(i,:,:,:)=dim_avg_n_Wrap(t_exp_mon(6*i:6*i+2,:,:,:),0)
   ctl_sum(i,:,:,:)=dim_avg_n_Wrap(t_ctl_mon(6*i+3:6*i+5,:,:,:),0)
   exp_sum(i,:,:,:)=dim_avg_n_Wrap(t_exp_mon(6*i+3:6*i+5,:,:,:),0)
 end do

  diff = new((/2, 342, 462, 10/), float)
  diff(0,:,:,:) = dim_avg_n_Wrap(exp_spr,0) - dim_avg_n_Wrap(ctl_spr,0)
  diff(1,:,:,:) = dim_avg_n_Wrap(exp_sum,0) - dim_avg_n_Wrap(ctl_sum,0)

  diff!0 = "season"
  diff!1 = "lat"
  diff!2 = "lon"
  diff!3 = "level"
  diff&lat = nlat
  diff&lon = nlon
  diff&lat@units = "degrees_north"
  diff&lon@units = "degrees_east"
  diff_re = diff(season|:, level|:, lat|:, lon|:)

  f=addfile("/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc","r")
  t_obs_mon = f->tm
  lat1d = f->lat
  lon1d = f->lon
  diff_re@lat2d = lat2d
  diff_re@lon2d = lon2d
  diff_ip = rcm2rgrid_Wrap(lat2d,lon2d,diff_re,lat1d,lon1d,1)
  mask_2d = t_obs_mon(0,:,:)
  mask_3d = conform(diff_ip, mask_2d, (/2,3/))
  diff_ip = where(ismissing(mask_3d), diff_ip@_FillValue, diff_ip)
  diff_ip = where(abs(diff_ip) .lt. 0.001, diff_ip@_FillValue, diff_ip)
  printVarSummary(diff_ip)
  printMinMax(diff_ip, 0)
  ; exit
  wks = gsn_open_wks("pdf","./figs/colmoff_ice_diff")
  plot= new(20,graphic)
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
  res@mpGridLineThicknessF  = 0.8

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

  ; ===== 1. 读取色板 =====
  purple = read_colormap_file("MPL_Purples")     ; 负值紫色（可换 GMT_purple / cmocean_matter）
  blue   = read_colormap_file("MPL_Blues")       ; 正值蓝色（可换 cmocean_dense / BlueDarkOrange18）

  ; ===== 2. 安全切片 + 预分配拼接 =====
  purple_r = purple(::-1, :)
  neg_raw = purple_r(10::2, :)     ; 深紫 → 浅紫
  pos_raw = blue(5::2, :)        ; 浅蓝 → 深蓝

  n_neg = 42                     ; ← 你想要的负颜色数量（30~60 随便调）
  n_pos = 42                     ; ← 正颜色数量（建议和 n_neg 一样或差 1~2）

  ; 防止下标越界（自动截取）
  n_neg := min((/n_neg, dimsizes(neg_raw(:,0))/))
  n_pos := min((/n_pos, dimsizes(pos_raw(:,0))/))

  neg_part = neg_raw(0:n_neg-1, :)
  pos_part = pos_raw(0:n_pos-1, :)

  ncols = dimsizes(neg_part(0,:))          ; 自动适配 3列(RGB) 或 4列(RGBA)

  ; === 核心：预分配数组（彻底告别 append 赋值不匹配）===
  ; total_colors = n_neg + 1 + n_pos         ; 负 + 白 + 正
  ; spliced = new( (/total_colors, ncols/), float )

  ; spliced(0:n_neg-1, :)           = neg_part
  ; spliced(n_neg, :)               = 1.0    ; 第 n_neg 行是纯白色（对应 0 值）
  ; spliced(n_neg+1:total_colors-1, :) = pos_part

  total_colors = n_neg + n_pos         ; 负 + 正
  spliced = new( (/total_colors, ncols/), float )

  spliced(0:n_neg-1, :)           = neg_part
  spliced(n_neg:total_colors-1, :) = pos_part

  res@cnFillPalette           = spliced ;""

  ; gsn_define_colormap(wks, "MPL_coolwarm")   ;MPL_BrBG  NCV_blu_red
  ; res@gsnSpreadColorStart = 2
  ; res@gsnSpreadColorEnd   = 256
  res@gsnRightString          = "kg/m~S~2~N~" ;~S~o~N~C m~S~3~N~/m~S~3~N~ 
  ; res@gsnRightStringFontHeightF= 0.025
  ; res@gsnLeftStringFontHeightF= 0.025
  ;res@gsnLeftStringOrthogonalPosF = -0.001
  res@tiYAxisFontHeightF      = 0.02
  res@tiYAxisAngleF           = 90
  res@cnLineLabelsOn        = False
  res@cnLevelSelectionMode    = "ManualLevels"
  res@cnMaxLevelValF          =  50.0  
  res@cnMinLevelValF          = -50.0  
  res@cnLevelSpacingF         =  0.5
  ; res@cnLevelSelectionMode   = "ExplicitLevels"
  ; res@cnLevels               = (/-0.5,-0.4,-0.3,-0.2,-0.1,0,0.5,1.0,1.5,2.0,2.5/)
  res@tiYAxisFontHeightF      = 0.02
  res@tiYAxisAngleF           = 90
  do j=0,9
     k = j * 2 
     res@tiYAxisString = "n=" + (j+1)
     res@gsnLeftString = "MAM"
     plot(k)   = gsn_csm_contour_map(wks,diff_ip(0,j,:,:),res)

     res@tiYAxisString = ""
     res@gsnLeftString = "JJA"
     plot(k+1) = gsn_csm_contour_map(wks,diff_ip(1,j,:,:),res)
  end do
  txres = True
  txres@txFontHeightF = 0.02
  txres@txFontColor = "black"

  txt = new(20, graphic)
  do p=0,19
     txt(p) = gsn_add_text(wks, plot(p), "EXP-CTL", lon1d({127}), lat1d({17}), txres)
  end do

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
  resP@gsnPanelBottom        = 0.08 
  ;/Surface Soil Temperature/Volumetric Water in Soil (m3/m3)/Soil Moisture/Liquid Water
  resP@gsnPanelLabelBar      = True
  resP@lbBoxLinesOn          = False      ; 去掉色块边框
  resP@pmLabelBarWidthF      = 0.8        ; 色标条宽度 (0.0-1.0)
  resP@pmLabelBarHeightF     = 0.06       ; 色标条高度
  resP@pmLabelBarOrthogonalPosF = -0.01   ; 垂直位置微调 (负数向上，正数向下)
  resP@lbBoxEndCapStyle         = "TriangleBothEnds"
  ; --- 设置单位标题 ---
  resP@lbTitleString         = "Ice Lens Difference" ; 标题
  resP@lbTitlePosition       = "Bottom"      ; 标题放在右侧
  resP@lbTitleFontHeightF    = 0.012        ; 标题字体大小
  resP@lbTitleDirection      = "Across"     ; 标题方向
  resP@lbTitleOffsetF        = 0.025         ; 标题与色标的距离
  gsn_panel(wks, plot, (/10,2/), resP)

end
