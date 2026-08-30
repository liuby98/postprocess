  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
  ;================read data 
begin
fs = addfile("hgt.sfc.nc","r")
sur=short2flt(fs->hgt)
sur1=sur(0,{5:60},{70:140})
f = addfile("/share/home/dq135/reference/era5_1991_2023_monthly_025x025.nc","r")
u_obs_mon = f->u(:,{850},{60:5},{70:140})
v_obs_mon = f->v(:,{850},{60:5},{70:140})
z_obs_mon = f->z(:,{850},{60:5},{70:140})
u = new((/33,221,281/),float)
v = new((/33,221,281/),float)
z = new((/33,221,281/),float)
do i=0,32
   u(i,:,:) = dim_avg_n_Wrap(u_obs_mon(6*i:6*i+2,:,:),0)
   v(i,:,:) = dim_avg_n_Wrap(v_obs_mon(6*i:6*i+2,:,:),0)
   z(i,:,:) = dim_avg_n_Wrap(z_obs_mon(6*i:6*i+2,:,:),0)
end do
umam = dim_avg_n_Wrap(u,0)
vmam = dim_avg_n_Wrap(v,0)
zmam = dim_avg_n_Wrap(z,0)
hgt = zmam/9.81
printVarSummary(umam)
printMinMax(umam,0)
printVarSummary(vmam)
printMinMax(vmam,0)
printVarSummary(hgt)
printMinMax(hgt,0)
; exit
umam@long_name = ""
umam@units     = ""
vmam@long_name = ""
vmam@units     = ""
wks=gsn_open_wks("png","era5_uvhgt850_JJA")
res0                   = True
res0@tiMainString      = ""
res0@gsnLeftString     = "OBS"
res0@gsnRightString    = "JJA"
res0@gsnDraw                 = False
res0@gsnFrame                = False
res0@gsnAddCyclic            = False
res0@tfDoNDCOverlay          = True    ; 原生投影：不变换数据 

res0@tmXTOn                 = False  
res0@tmYROn                 = False

res0@mpLimitMode            = "Corners"          ; choose region of map
res0@mpLeftCornerLatF       = 48
res0@mpLeftCornerLonF       = 60
res0@mpRightCornerLatF      = 10
res0@mpRightCornerLonF      = 135

res0@mpProjection           = "LambertConformal"
res0@mpLambertParallel1F    = 30.
res0@mpLambertParallel2F    = 60.
res0@mpLambertMeridianF     = 110.
res0@mpCenterLatF           = 35.1778 
res0@mpCenterLonF           = 110.   

res0@cnFillOn                = True
res0@cnLinesOn               = False
res0@mpFillOn                = False
res0@mpOutlineOn             = True

res0@lbBoxLinesOn            = False
res0@lbOrientation           = "Horizontal"
; res0@lbBoxEndCapStyle        = "TriangleBothEnds"
res0@lbTitlePosition         = "Bottom"
res0@lbTitleDirection        = "Across"
res0@lbTitleAngleF           =  0

res0@mpDataBaseVersion      = "MediumRes"        
res0@pmTickMarkDisplayMode  = "Always"
res0@mpOutlineOn = True   ;开启海岸线绘制
res0@mpGeophysicalLineColor  = "grey60"   ;设置海岸线颜色
res0@mpGeophysicalLineThicknessF = 0.3  
res0@mpOutlineBoundarySets = "GeophysicalAndUSStates"  ; 改善中国边界（如果适用）

res0@cnMissingValFillColor = -1  ; 缺失值透明
res0@cnFillDrawOrder    = "PreDraw"  ; 先绘制填充，再叠加 shapefile
res0@cnLineThicknessF  = 3
res0@cnLevelSelectionMode = "ManualLevels"
res0@cnMaxLevelValF       =  1550
res0@cnMinLevelValF       =  1350
res0@cnLevelSpacingF      =  5
; res0@cnLineLabelDensityF  =  0.5                       ;increase the number of line labels/line
; res0@cnLineLabelInterval  =  1
res0@cnFillPalette        = "WhiteBlueGreenYellowRed" ;
;=====================================================draw contour plot 
plot0=gsn_csm_contour_map(wks,hgt(::-1,:),res0)

;===============================================================res vector
res                           = True               ; plot mods desired
res@gsnDraw                   = False
res@gsnFrame                  = False
res@vcRefMagnitudeF           = 5.0             ; define vector ref mag
res@vcRefLengthF              = 0.05            ; define length of vec 
res@vcRefAnnoOrthogonalPosF   = -1.39
res@vcRefAnnoParallelPosF     = 1
res@vcRefAnnoArrowLineColor   = "red"           ; change ref vector color
res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
res@vcMinDistanceF            = 0.03           ; thin vectors
res@vcLineArrowColor          = "black"           ; change vector color
res@vcLineArrowThicknessF     = 2.0             ; change vector thickness
res@vcLineArrowHeadMaxSizeF   = 0.01
res@vcRefAnnoString2On        = True
res@vcRefAnnoString2          = "m/s"

;=====================================================draw vector plot 
plot=gsn_csm_vector(wks,umam,vmam,res)
overlay(plot0,plot)


; res1=True
; res1@gsnDraw             = False       ; Do not draw plot
; res1@gsnFrame            = False
; res1@gsnRightString      = ""
; res1@gsnLeftString       = ""
; ; res1@cnFillOn           = True
; res1@trGridType = "TriangularMesh"
; res1@cnLevelSelectionMode = "ExplicitLevels"
; res1@cnLevels            = (/5800/)
; res1@cnInfoLabelOn       = False
; res1@cnLineLabelInterval = 1
; res1@cnLineLabelDensityF = 0.5 
; res1@cnLineColor ="magenta3"
; res1@cnLineThicknessF = 3
; res1@cnLineLabelPlacementMode = "Computed"  
; plot1=gsn_csm_contour(wks,hgt,res1)
; overlay(plot0,plot1)

;=============================================================
res1                    = True
res1@gsnDraw            = False
res1@gsnFrame           = False
res1@cnFillOn           = True
res1@cnLinesOn          = False
;res1@gsnAddCyclic      = False
;res1@cnFillDrawOrder   = "PreDraw"
res1@cnLevelSelectionMode = "ExplicitLevels";set explicit contour levels
res1@cnLevels             = (/2400/)
res1@cnLineLabelsOn       = False
res1@cnFillColors         = (/"Transparent","Gray"/)
res1@cnInfoLabelOn        = False
res1@lbLabelAutoStride    = False
res1@lbLabelBarOn         = False
sur1@long_name=""
sur1@units=""
plot11=gsn_csm_contour(wks,sur1,res1)
overlay(plot0,plot11)

dir                    = "../shapefile_China/"
  fShape1Name          = dir + "china.shp"
  fShape2Name          = dir + "south_china_sea.shp"
  fShape3Name          = dir + "province.shp"

  resLine                     = True
  resLine@gsLineThicknessF    = 0.8
  resLine@gsLineColor         = "black"
  line  = gsn_add_shapefile_polylines(wks, plot0, fShape1Name,resLine)
  scs   = gsn_add_shapefile_polylines(wks, plot0, fShape2Name,resLine)
  province  = gsn_add_shapefile_polylines(wks, plot0, fShape3Name,resLine)

draw(plot0)
frame(wks)
end