load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
begin
f1 = addfile("/share/home/dq135/wrfpost/wrfpost_u_1991_2023_JJA_seasmean.nc","r")
f2 = addfile("/share/home/dq135/wrfpost/wrfpost_v_1991_2023_JJA_seasmean.nc","r")
f3 = addfile("/share/home/dq135/wrfpost/wrfpost_q_1991_2023_JJA_seasmean.nc","r")
f4 = addfile("/share/home/dq135/wrfpost/wrfout_1991_2023_JJA_seasmean.nc","r")

f5 = addfile("/share/home/dq135/wrfpost/wrfpost_u_nog_1991_2023_JJA_seasmean.nc","r")
f6 = addfile("/share/home/dq135/wrfpost/wrfpost_v_nog_1991_2023_JJA_seasmean.nc","r")
f7 = addfile("/share/home/dq135/wrfpost/wrfpost_q_nog_1991_2023_JJA_seasmean.nc","r")
f8 = addfile("/share/home/dq135/wrfpost/wrfout_nog_1991_2023_JJA_seasmean.nc","r")

f_wrf=addfile("/share/home/dq117/CN05.1/CN05.1_Tm_2020_daily_025x025.nc","r")
f_input=addfile("/share/home/dq135/CRESM_0307S92/CASES/RUN_1991/wrfinput_d01", "r")
t_obs_mon = f_wrf->tm
t_obs = dim_avg_n_Wrap(t_obs_mon,0)
lat1d = f_wrf->lat
lon1d = f_wrf->lon
lat2d = f_input->XLAT(0,:,:)
lon2d = f_input->XLONG(0,:,:)

u  = f1->U(2,:,:) 
v  = f2->V(2,:,:) 
q  = f3->QVAPOR (2,:,:) 
ps = f4->PSFC(0,:,:)

u_exp  = f5->U(2,:,:) 
v_exp  = f6->V(2,:,:) 
q_exp  = f7->QVAPOR (2,:,:) 
ps_exp = f8->PSFC(0,:,:)

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

;=====850hPa的水平水汽通量和水汽通量散度及水汽通量散度的垂直积分=========
q       = q*1000
q_exp   = q_exp*1000
q@units = "g/kg"
q_exp@units = "g/kg"
uq      = q*u/9.8
vq      = q*v/9.8
uq_exp  = q_exp*u_exp/9.8
vq_exp  = q_exp*v_exp/9.8

; uq!0 = "level"
uq!0 = "lat"
uq!1 = "lon"
uq&lat = lat
uq&lon = lon
; vq!0 = "level"
vq!0 = "lat"
vq!1 = "lon"
vq&lat = lat
vq&lon = lon

; uq_exp!0 = "level"
uq_exp!0 = "lat"
uq_exp!1 = "lon"
uq_exp&lat = lat
uq_exp&lon = lon
; vq_exp!0 = "level"
vq_exp!0 = "lat"
vq_exp!1 = "lon"
vq_exp&lat = lat
vq_exp&lon = lon

qv_magnitude = sqrt(uq^2 + vq^2)
qv_magnitude!0 = "lat"
qv_magnitude!1 = "lon"
qv_magnitude&lat = lat
qv_magnitude&lon = lon
qv_magnitude@units = "g/(cm*s)" ;g/(cm s hPa)
qv_magnitude_exp = sqrt(uq_exp^2 + vq_exp^2)
qv_magnitude_exp!0 = "lat"
qv_magnitude_exp!1 = "lon"
qv_magnitude_exp&lat = lat
qv_magnitude_exp&lon = lon
qv_magnitude_exp@units = "g/(cm*s)" ;g/(cm s hPa)

uq@lat2d = lat2d
uq@lon2d = lon2d
uq_ip    = rcm2rgrid_Wrap(lat2d,lon2d,uq,lat1d,lon1d,1)
uq_ip    = where(ismissing(t_obs),uq_ip@_FillValue,uq_ip)
vq@lat2d = lat2d
vq@lon2d = lon2d
vq_ip    = rcm2rgrid_Wrap(lat2d,lon2d,vq,lat1d,lon1d,1)
vq_ip    = where(ismissing(t_obs),vq_ip@_FillValue,vq_ip)
uq_exp@lat2d = lat2d
uq_exp@lon2d = lon2d
uq_exp_ip    = rcm2rgrid_Wrap(lat2d,lon2d,uq_exp,lat1d,lon1d,1)
uq_exp_ip    = where(ismissing(t_obs),uq_exp_ip@_FillValue,uq_exp_ip)
vq_exp@lat2d = lat2d
vq_exp@lon2d = lon2d
vq_exp_ip    = rcm2rgrid_Wrap(lat2d,lon2d,vq_exp,lat1d,lon1d,1)
vq_exp_ip    = where(ismissing(t_obs),vq_exp_ip@_FillValue,vq_exp_ip)

qv_magnitude@lat2d = lat2d
qv_magnitude@lon2d = lon2d
qv_magnitude_ip    = rcm2rgrid_Wrap(lat2d,lon2d,qv_magnitude,lat1d,lon1d,1)
qv_magnitude_ip    = where(ismissing(t_obs),qv_magnitude_ip@_FillValue,qv_magnitude_ip)
qv_magnitude_exp@lat2d = lat2d
qv_magnitude_exp@lon2d = lon2d
qv_magnitude_exp_ip    = rcm2rgrid_Wrap(lat2d,lon2d,qv_magnitude_exp,lat1d,lon1d,1)
qv_magnitude_exp_ip    = where(ismissing(t_obs),qv_magnitude_exp_ip@_FillValue,qv_magnitude_exp_ip)
diff= qv_magnitude_ip
diff= qv_magnitude_ip - qv_magnitude_exp_ip

printVarSummary(diff)
printMinMax(diff,0)
; exit
 
wks  = gsn_open_wks("pdf","vaporflux_diff_JJA")   
res0 = True                            ; res0 probability plots
res0@gsnDraw             = False       ; Do not draw plot
res0@gsnFrame            = False       ; Do not advance frome
res0@tiMainString        = " "
res0@gsnRightString      = "g/(cm*s)"
res0@gsnLeftString       = "CTL-EXP"
res0@gsnLeftStringOrthogonalPosF = -0.008
res0@gsnRightStringOrthogonalPosF = -0.008
res0@cnLinesOn           = False
res0@mpFillOn            = False
res0@cnFillOn            = True
res0@lbLabelBarOn        = True
res0@mpOutlineOn         = True
res0@gsnAddCyclic        = False
; res0@vpWidthF            = 0.6
; res0@vpHeightF           = 0.8

res0@tmXTOn                 = False  
res0@tmYROn                 = False

res0@mpMinLonF               = 73 
res0@mpMaxLonF               = 136 
res0@mpCenterLonF            = 105  
res0@mpMinLatF               = 10 
res0@mpMaxLatF               = 60 

res0@lbBoxLinesOn            = False
res0@lbOrientation           = "Vertical"
res0@lbLabelAutoStride       = False
res0@lbLabelStride           = 1
res0@lbBoxEndCapStyle        = "TriangleBothEnds"
res0@lbTitleString           = " "
res0@lbTitlePosition         = "Right"
res0@lbTitleDirection        = "Across"
res0@lbTitleAngleF           =  90

res0@mpDataBaseVersion      = "MediumRes"        
res0@pmTickMarkDisplayMode  = "Always"
res0@mpOutlineOn            = True       ;开启海岸线绘制
res0@mpGeophysicalLineColor = "grey60"   ;设置海岸线颜色
res0@mpGeophysicalLineThicknessF = 0.3  
res0@mpOutlineBoundarySets = "GeophysicalAndUSStates"  ; 改善中国边界（如果适用）

res0@cnMissingValFillColor = -1         ; 缺失值透明
res0@cnFillDrawOrder       = "PreDraw"  ; 先绘制填充，再叠加 shapefile
res0@cnFillPalette         = "MPL_BrBG"  ; WhiteBlue
res0@cnLevelSelectionMode  = "ManualLevels"      ; set manual contour levels
res0@cnMinLevelValF        =-1.0                  ; set min contour level
res0@cnMaxLevelValF        = 1.0                   ; set max contour level
res0@cnLevelSpacingF       = 0.2
res0@cnLineLabelInterval   = 1

plot = gsn_csm_contour_map(wks, diff, res0)

; res                           = True             ; plot mods desired
; res@gsnDraw                   = False
; res@gsnFrame                  = False
; res@lbLabelBarOn              = True
; res@vcRefMagnitudeF           = 2.0             ; define vector ref mag
; res@vcRefLengthF              = 0.05             ; define length of vec 
; res@vcRefAnnoOrthogonalPosF   = -1.             ; move ref vector
; ;res@vcLineArrowHeadMaxSizeF   = 0.01
; res@vcRefAnnoArrowLineColor   = "red"            ; change ref vector color
; res@vcRefAnnoArrowUseVecColor = False            ; don't use vec color for ref
; res@vcMinDistanceF            = 0.02             ; thin vectors
; res@vcGlyphStyle              = "CurlyVector"    ; turn on curley vectors
; res@vcLevelPalette             = "GMT_drywet"  ; 向量颜色hotcolr_19lev
; ;res@vcLineArrowColor          = "red"           ; change vector color
; ;res@vcLineArrowThicknessF     = 2.0             ; change vector thickness
; res@vcVectorDrawOrder         = "PostDraw"      ; draw vectors last
; res@vcRefAnnoString2On        = True
; res@vcRefAnnoString2          = " " 
; res@vcRefAnnoBackgroundColor  = "Transparent"
; res@vcRefAnnoFontHeightF     = 0.012
; res@gsnLeftString             = " "
; res@gsnRightString            = " "
; res@vcRefAnnoPerimOn          = False
; ; vector = gsn_csm_vector(wks, uq(:,:), vq(:,:), res)
; vector = gsn_csm_vector_scalar(wks, uq_exp_ip, vq_exp_ip, qv_magnitude_exp_ip, res)
; ;---Overlay the vectors on the contour/map plot
; overlay(plot,vector)

  dir                    = "../shapefile_China/"
  fShape1Name            = dir + "china.shp"
  fShape2Name            = dir + "south_china_sea.shp"
  fShape3Name            = dir + "province.shp"

  resLine                     = True
  resLine@gsLineThicknessF    = 0.8
  resLine@gsLineColor         = "black"
  line  = gsn_add_shapefile_polylines(wks, plot, fShape1Name,resLine)
  scs   = gsn_add_shapefile_polylines(wks, plot, fShape2Name,resLine)

draw(plot)
frame(wks)
end