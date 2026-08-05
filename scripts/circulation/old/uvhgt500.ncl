  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
  begin
  ;========================================================= read data
f1= addfile("/share/home/dq135/wrfpost/wrfpost_h_1991_2023_JJA_seasmean.nc","r")
f2= addfile("/share/home/dq135/wrfpost/wrfpost_u_1991_2023_JJA_seasmean.nc","r")
f3= addfile("/share/home/dq135/wrfpost/wrfpost_v_1991_2023_JJA_seasmean.nc","r")
; f4= addfile("/share/home/dq135/wrfpost/hgt.sfc.nc","r")
f4= addfile("/share/home/dq135/wrfpost/wrfpost_h_nog_1991_2023_JJA_seasmean.nc","r")
f5= addfile("/share/home/dq135/wrfpost/wrfpost_u_nog_1991_2023_JJA_seasmean.nc","r")
f6= addfile("/share/home/dq135/wrfpost/wrfpost_v_nog_1991_2023_JJA_seasmean.nc","r")
fs = addfile("hgt.sfc.nc","r")
sur=short2flt(fs->hgt)
sur1=sur(0,{55:10},{70:150})
f_input = addfile("/share/home/dq135/wrfpost/wrfinput_d01", "r")
  ; 获取网格的纬度和经度
  lon = f_input->XLONG(0,0,:)
  lat = f_input->XLAT(0,:,0)
  lon@standard_name = "longitude"
  lon@long_name = "longitude"
  lon@units = "degrees_east"
  lon@axis = "X"

  lat@standard_name = "latitude"
  lat@long_name = "latitude"
  lat@units = "degrees_north"
  lat@axis = "Y"

hgt=short2flt(f1->GH({500},:,:))
u=short2flt(f2->U({500},:,:))
v=short2flt(f3->V({500},:,:))
hgt_exp=short2flt(f4->GH({500},:,:))
u_exp=short2flt(f5->U({500},:,:))
v_exp=short2flt(f6->V({500},:,:))

hgt!0 = "lat"
hgt!1 = "lon"
hgt&lat = lat
hgt&lon = lon
hgt_exp!0 = "lat"
hgt_exp!1 = "lon"
hgt_exp&lat = lat
hgt_exp&lon = lon
hgt_diff = hgt - hgt_exp
u_diff = u - u_exp
v_diff = v - v_exp
copy_VarCoords(hgt, hgt_diff)
copy_VarCoords(u, u_diff)
copy_VarCoords(v, v_diff)
printVarSummary(hgt_diff) 
printMinMax(hgt_diff,0)
; exit
hgt@long_name = ""
hgt@units = ""
u@long_name = ""
u@units = ""
v@long_name = ""
v@units = ""
hgt_exp@long_name = ""
hgt_exp@units = ""
u_exp@long_name = ""
u_exp@units = ""
v_exp@long_name = ""
v_exp@units = ""
;================================================================= plot
wks=gsn_open_wks("png","uvhgt500_diff_JJA")
res0                   = True
res0@tiMainString      = ""
res0@gsnLeftString     = "JJA"
res0@gsnRightString    = "CTL-EXP"
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
res0@cnMaxLevelValF       =  10;6020
res0@cnMinLevelValF       = -10;5320
res0@cnLevelSpacingF      =   1;15
res0@cnLineLabelDensityF  = 0.5                       ;increase the number of line labels/line
res0@cnLineLabelInterval  = 1
res0@cnFillPalette        = "BlWhRe" ;WhiteBlueGreenYellowRed
;=====================================================draw contour plot 
plot0=gsn_csm_contour_map(wks,hgt_diff,res0)

;===============================================================res vector
res                           = True               ; plot mods desired
res@gsnDraw                   = False
res@gsnFrame                  = False
res@vcRefMagnitudeF           = 1.0             ; define vector ref mag
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
plot=gsn_csm_vector(wks,u_diff,v_diff,res)
overlay(plot0,plot)

;=============================================================
; res1                    = True
; res1@gsnDraw            = False
; res1@gsnFrame           = False
; res1@cnFillOn           = True
; res1@cnLinesOn          = False
; ;res1@gsnAddCyclic      = False
; ;res1@cnFillDrawOrder   = "PreDraw"
; res1@cnLevelSelectionMode = "ExplicitLevels";set explicit contour levels
; res1@cnLevels             = (/2400/)
; res1@cnLineLabelsOn       = False
; res1@cnFillColors         = (/"Transparent","Gray"/)
; res1@cnInfoLabelOn        = False
; res1@lbLabelAutoStride    = False
; res1@lbLabelBarOn         = False
; sur1@long_name=""
; sur1@units=""
; plot11=gsn_csm_contour(wks,sur1,res1)
; overlay(plot0,plot11)

  dir                    = "../shapefile_China/"
  fShape1Name            = dir + "china.shp"
  fShape2Name            = dir + "south_china_sea.shp"
  fShape3Name            = dir + "province.shp"

resLine                     = True
  resLine@gsLineThicknessF    = 0.8
  resLine@gsLineColor         = "black"
  line  = gsn_add_shapefile_polylines(wks, plot0, fShape1Name,resLine)
  scs   = gsn_add_shapefile_polylines(wks, plot0, fShape2Name,resLine)
  province  = gsn_add_shapefile_polylines(wks, plot0, fShape3Name,resLine)

draw(plot0)
frame(wks)
end