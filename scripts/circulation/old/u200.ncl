  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
  ;================read data 
begin
f1=addfile("/share/home/dq135/wrfpost/u200_1991_2023_JJA_seasmean.nc","r")
f2=addfile("/share/home/dq135/wrfpost/u200_nog_1991_2023_JJA_seasmean.nc","r")
f3=addfile("uwnd.mon.ltm.nc","r")
time=f3->time
tim= cd_calendar(time,-1)
u3=f3->uwnd(:,{200},{10:55},{65:155})
ujja=month_to_season (u3, "JJA")
u200cli=dim_avg_n_Wrap(ujja,0)
printVarSummary(u200cli)
printMinMax(u200cli,0)

u1=short2flt(f1->U(:,:))
u2=short2flt(f2->U(:,:))

udiff=u1-u2
copy_VarCoords(u1,udiff)

f = addfile("/share/home/dq135/reference/era5_1991_2023_monthly_025x025.nc","r")
u_obs_mon = f->u(:,{200},:,:)
u200 = new((/33,221,281/),float)
do i=0,32
   u200(i,:,:) = dim_avg_n_Wrap(u_obs_mon(6*i+3:6*i+5,:,:),0)
end do
umam = dim_avg_n_Wrap(u200,0)

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
  
  u1!0 = "lat"
  u1!1 = "lon"
  u1&lat = lat
  u1&lon = lon
  udiff!0 = "lat"
  udiff!1 = "lon"
  udiff&lat = lat
  udiff&lon = lon
;   printMinMax(udiff,0)
; exit
wks=gsn_open_wks("png","u200_JJA")
res2 = True                            ; res2 probability plots
  res2@gsnDraw             = False       ; Do not draw plot
  res2@gsnFrame            = False       ; Do not advance frome
  res2@tiMainString        = "JJA"
  res2@gsnRightString      = "m/s"
  res2@gsnLeftString       = "u200:CTL minus EXP"
  res2@mpFillOn            = False
  res2@cnFillOn            = True
  res2@lbLabelBarOn        = True
  res2@mpOutlineOn         = True
  res2@gsnAddCyclic        = False

  res2@tfDoNDCOverlay         = True    ; 原生投影：不变换数据 
  res2@tmXTOn                 = False  
  res2@tmYROn                 = False

  res2@mpLimitMode            = "Corners"          ; choose region of map
  res2@mpLeftCornerLatF       = 48
  res2@mpLeftCornerLonF       = 60
  res2@mpRightCornerLatF      = 10
  res2@mpRightCornerLonF      = 135

  res2@mpProjection           = "LambertConformal"
  res2@mpLambertParallel1F    = 30.
  res2@mpLambertParallel2F    = 60.
  res2@mpLambertMeridianF     = 110.
  res2@mpCenterLatF           = 35.1778 
  res2@mpCenterLonF           = 110.  

  res2@mpDataBaseVersion      = "MediumRes"        
  res2@pmTickMarkDisplayMode  = "Always"
  res2@mpOutlineOn = True   ;开启海岸线绘制
  res2@mpGeophysicalLineColor  = "grey60"   ;设置海岸线颜色
  res2@mpGeophysicalLineThicknessF = 0.3  
  res2@mpOutlineBoundarySets = "GeophysicalAndUSStates"  ; 改善中国边界（如果适用）

  res2@cnMissingValFillColor = -1  ; 缺失值透明
  res2@cnFillDrawOrder       = "PreDraw"  ; 先绘制填充，再叠加 shapefile
  res2@cnFillPalette         = "temp_19lev"
  res2@cnLevelSelectionMode  = "ManualLevels"
  res2@cnMinLevelValF        = -1.8        ; set min contour level
  res2@cnMaxLevelValF        =  1.8        ; set max contour level
  res2@cnLevelSpacingF       = 0.2       ; set contour spacing
  res2@cnInfoLabelOn         = False
;   res2@cnLineLabelsOn        = True
;   res2@cnLineLabelInterval   = 1
;   res2@cnLineLabelPlacementMode  = "Computed" 
  plot2=gsn_csm_contour_map(wks,udiff,res2)
  

res3=True
res3@gsnDraw             = False       ; Do not draw plot
res3@gsnFrame            = False
res3@gsnRightString         = ""
res3@gsnLeftString        = ""
;res@cnFillOn=True
res3@cnLineLabelsOn        = True
res3@cnLevelSelectionMode = "ExplicitLevels"
res3@cnLevels    = (/ 20.,30./)
res3@cnInfoLabelOn         = False
res3@cnLineLabelInterval   = 1
res3@cnLineLabelDensityF   = 0.5 
res3@cnLineColor           = "magenta3"
res3@cnLineThicknessF       = 3
res3@cnLineLabelPlacementMode   = "Computed"  
plot3=gsn_csm_contour(wks,u200cli,res3)
overlay(plot2,plot3)


dir                    = "../shapefile_China/"
  fShape1Name          = dir + "china.shp"
  fShape2Name          = dir + "south_china_sea.shp"
  fShape3Name          = dir + "province.shp"

  resLine                     = True
  resLine@gsLineThicknessF    = 0.8
  resLine@gsLineColor         = "black"
  line  = gsn_add_shapefile_polylines(wks, plot2, fShape1Name,resLine)
  scs   = gsn_add_shapefile_polylines(wks, plot2, fShape2Name,resLine)
  province  = gsn_add_shapefile_polylines(wks, plot2, fShape3Name,resLine)


draw(plot2)
frame(wks)
end