  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
  ;================read data 
begin
f = addfile("/share/home/dq135/reference/era5_1991_2023_monthly_025x025.nc","r")
u_obs_mon = f->u(:,{200},:,:)
u200 = new((/33,221,281/),float)
do i=0,32
   u200(i,:,:) = dim_avg_n_Wrap(u_obs_mon(6*i:6*i+2,:,:),0)
end do
umam = dim_avg_n_Wrap(u200,0)
; copy_VarCoords(u200,umam)
printVarSummary(umam)
printMinMax(umam,0)
; exit

; lat2d = f->u(0,0,:,:)
; lon2d = f->u(0,0,:,:)
; f_input = addfile("/share/home/dq135/wrfpost/wrfinput_d01", "r")
;   ; 获取网格的纬度和经度
;   lon = f_input->XLONG(0,0,:)
;   lat = f_input->XLAT(0,:,0)
;   lon@standard_name = "longitude"
;   lon@long_name = "longitude"
;   lon@units = "degrees_east"
;   lon@axis = "X"

;   lat@standard_name = "latitude"
;   lat@long_name = "latitude"
;   lat@units = "degrees_north"
;   lat@axis = "Y"

; umam_ip = rcm2rgrid_Wrap(lat2d,lon2d,umam,lat,lon,1)
; printVarSummary(umam_ip)
; printMinMax(umam_ip,0)
; exit
wks=gsn_open_wks("png","era5_MAM_u200")
res = True                            ; res probability plots
  res@gsnDraw             = False       ; Do not draw plot
  res@gsnFrame            = False       ; Do not advance frome
  res@tiMainString        = "MAM"
  res@gsnRightString      = "m/s"
  res@gsnLeftString       = "u200"
  res@mpFillOn            = False
  res@cnFillOn            = True
  res@lbLabelBarOn        = True
  res@mpOutlineOn         = True
  res@gsnAddCyclic        = False

  res@tfDoNDCOverlay         = True    ; 原生投影：不变换数据 
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

  res@mpDataBaseVersion      = "MediumRes"        
  res@pmTickMarkDisplayMode  = "Always"
  res@mpOutlineOn = True   ;开启海岸线绘制
  res@mpGeophysicalLineColor  = "grey60"   ;设置海岸线颜色
  res@mpGeophysicalLineThicknessF = 0.3  
  res@mpOutlineBoundarySets = "GeophysicalAndUSStates"  ; 改善中国边界（如果适用）
  res@cnMissingValFillColor = -1  ; 缺失值透明
  res@cnFillDrawOrder       = "PreDraw"  ; 先绘制 填充，再叠加 shapefile
  res@cnFillPalette         = "WhiteBlueGreenYellowRed"
  res@cnLinesOn             = False
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnMinLevelValF        = -20.        ; set min contour level
  res@cnMaxLevelValF        =  40.        ; set max contour level
  res@cnLevelSpacingF       =  1.       ; set contour spacing
  res@cnInfoLabelOn         = False
  plot = gsn_csm_contour_map(wks,umam(::-1,:),res)

  
;   res=True
;   res@gsnDraw             = False  
;   res@gsnFrame            =False
;   res@gsnRightString         = ""
;   res@gsnLeftString        = ""
;   res@cnLineLabelsOn = False 
;   res@cnInfoLabelOn       = False 
;   res@lbLabelBarOn   = False
;   res@cnFillOn=True
;   res@cnLinesOn          = False
;   res@cnLevelSelectionMode = "ExplicitLevels"
;   res@cnMonoFillPattern    = False
;   res@cnLevels          = (/90.,95./)
;   res@cnFillPatterns      = (/-1,16,17/)  
;   plot=gsn_csm_contour(wks,prob1,res)
;   overlay(plot2,plot)
  

; res1=True
; res1@gsnDraw             = False       ; Do not draw plot
; res1@gsnFrame            = False
; res1@gsnRightString      = ""
; res1@gsnLeftString       = ""
; ; res1@cnFillOn           = True
; res1@cnLevelSelectionMode = "ExplicitLevels"
; res1@cnLevels            = (/ 20.,30./)
; res1@cnInfoLabelOn       = False
; res1@cnLineLabelInterval = 1
; res1@cnLineLabelDensityF = 0.5 
; res1@cnLineColor ="magenta3"
; res1@cnLineThicknessF = 3
; res1@cnLineLabelPlacementMode = "Computed"  
; plot1=gsn_csm_contour(wks,umam,res1)
; overlay(plot,plot1)

dir                    = "../shapefile_China/"
  fShape1Name          = dir + "china.shp"
  fShape2Name          = dir + "south_china_sea.shp"
  fShape3Name          = dir + "province.shp"

  resLine                     = True
  resLine@gsLineThicknessF    = 0.8
  resLine@gsLineColor         = "black"
  line  = gsn_add_shapefile_polylines(wks, plot, fShape1Name,resLine)
  scs   = gsn_add_shapefile_polylines(wks, plot, fShape2Name,resLine)
  province  = gsn_add_shapefile_polylines(wks, plot, fShape3Name,resLine)


draw(plot)
frame(wks)
end