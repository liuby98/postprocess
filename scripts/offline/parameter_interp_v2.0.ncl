begin
  path="/share/home/dq135/CRESM_0307S92/SCRIPTS/ICBCScripts/PrepCoLM/CN_15km/"
  f1=addfile(path+"CoLMrun_nogravel/unstructured_cwrf_CN_15km_nogravel/landdata/diag/soil_parameters_2005.nc","r")
  f2=addfile(path+"CoLMrun_gravel/unstructured_cwrf_CN_15km_gravel/landdata/diag/soil_parameters_2005.nc","r")
  lat1d = f1->lat({8:60})  ;15:55
  lon1d = f1->lon({58:162}) ;70:140
  lat1d@units = "degrees_north"
  lon1d@units = "degrees_east"

  ;===================================================================
  ; 1. 变量名列表定义 
  ;===================================================================
  tkdry   = (/"tkdry_l1","tkdry_l2","tkdry_l3","tkdry_l4","tkdry_l5","tkdry_l6","tkdry_l7","tkdry_l8"/)
  tksatu  = (/"tksatu_l1","tksatu_l2","tksatu_l3","tksatu_l4","tksatu_l5","tksatu_l6","tksatu_l7","tksatu_l8"/)
  tksatf  = (/"tksatf_l1","tksatf_l2","tksatf_l3","tksatf_l4","tksatf_l5","tksatf_l6","tksatf_l7","tksatf_l8"/)
  theta_s = (/"theta_s_l1","theta_s_l2","theta_s_l3","theta_s_l4","theta_s_l5","theta_s_l6","theta_s_l7","theta_s_l8"/)
  
  BA_alpha     = (/"BA_alpha_l1","BA_alpha_l2","BA_alpha_l3","BA_alpha_l4","BA_alpha_l5","BA_alpha_l6","BA_alpha_l7","BA_alpha_l8"/)
  BA_beta      = (/"BA_beta_l1","BA_beta_l2","BA_beta_l3","BA_beta_l4","BA_beta_l5","BA_beta_l6","BA_beta_l7","BA_beta_l8"/)
  vf_sand_s    = (/"vf_sand_s_l1","vf_sand_s_l2","vf_sand_s_l3","vf_sand_s_l4","vf_sand_s_l5","vf_sand_s_l6","vf_sand_s_l7","vf_sand_s_l8"/)
  vf_om_s      = (/"vf_om_s_l1","vf_om_s_l2","vf_om_s_l3","vf_om_s_l4","vf_om_s_l5","vf_om_s_l6","vf_om_s_l7","vf_om_s_l8"/)
  vf_gravels_s = (/"vf_gravels_s_l1","vf_gravels_s_l2","vf_gravels_s_l3","vf_gravels_s_l4","vf_gravels_s_l5","vf_gravels_s_l6","vf_gravels_s_l7","vf_gravels_s_l8"/)
  vf_quartz    = (/"vf_quartz_mineral_s_l1","vf_quartz_mineral_s_l2","vf_quartz_mineral_s_l3","vf_quartz_mineral_s_l4","vf_quartz_mineral_s_l5","vf_quartz_mineral_s_l6","vf_quartz_mineral_s_l7","vf_quartz_mineral_s_l8"/)
  vf_clay_s    = (/"vf_clay_s_l1","vf_clay_s_l2","vf_clay_s_l3","vf_clay_s_l4","vf_clay_s_l5","vf_clay_s_l6","vf_clay_s_l7","vf_clay_s_l8"/)

  ; heat capacity of soil solids [J/(m3 K)]
  csol = (/"csol_l1","csol_l2","csol_l3","csol_l4","csol_l5","csol_l6","csol_l7","csol_l8"/)
  
  ;===================================================================
  ; 2. 读取参照变量及缺测处理
  ;===================================================================
  f3=addfile("/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/colmoff_2001-2023_monmean_nogravel.nc","r")
  f4=addfile("/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/colmoff_2001-2023_monmean_gravel.nc","r")
  t_ctl_mon = tofloat(f3->f_h2osoi(:,:,:,:))
  t_exp_mon = tofloat(f4->f_h2osoi(:,:,:,:))

  ; 读取冬季(DJF)数据
  f3_djf = addfile("/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/colmoff_2001-2017_DJF_nogravel.nc", "r")
  f4_djf = addfile("/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL/colmoff_2001-2017_DJF_gravel.nc", "r")
  t_ctl_win = tofloat(f3_djf->f_h2osoi(:,:,:,:))
  t_exp_win = tofloat(f4_djf->f_h2osoi(:,:,:,:))

  f_input = addfile("/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01","r")
  lat2d = f_input->XLAT(0,:,:)
  lon2d = f_input->XLONG(0,:,:)

  ctl_win = new((/17,342,462,10/),float)
  exp_win = new((/17,342,462,10/),float)
  ctl_spr = new((/17,342,462,10/),float)
  exp_spr = new((/17,342,462,10/),float)
  ctl_sum = new((/17,342,462,10/),float)
  exp_sum = new((/17,342,462,10/),float)
  do i=0,16
    ctl_win(i,:,:,:)=dim_avg_n_Wrap(t_ctl_win(3*i:3*i+2,:,:,:),0)
    exp_win(i,:,:,:)=dim_avg_n_Wrap(t_exp_win(3*i:3*i+2,:,:,:),0)
    ctl_spr(i,:,:,:)=dim_avg_n_Wrap(t_ctl_mon(6*i:6*i+2,:,:,:),0)
    exp_spr(i,:,:,:)=dim_avg_n_Wrap(t_exp_mon(6*i:6*i+2,:,:,:),0)
    ctl_sum(i,:,:,:)=dim_avg_n_Wrap(t_ctl_mon(6*i+3:6*i+5,:,:,:),0)
    exp_sum(i,:,:,:)=dim_avg_n_Wrap(t_exp_mon(6*i+3:6*i+5,:,:,:),0)
  end do

  ; 计算各季节多年平均
  theta_ctl_win = dim_avg_n_Wrap(ctl_win, 0)
  theta_exp_win = dim_avg_n_Wrap(exp_win, 0)
  theta_ctl_spr = dim_avg_n_Wrap(ctl_spr(:,:,:,:),0)  
  theta_exp_spr = dim_avg_n_Wrap(exp_spr(:,:,:,:),0)  
  theta_ctl_sum = dim_avg_n_Wrap(ctl_sum(:,:,:,:),0)
  theta_exp_sum = dim_avg_n_Wrap(exp_sum(:,:,:,:),0)

  ; 缺测处理
  theta_ctl_win = where(theta_ctl_win .le. -1e+35, theta_ctl_win@_FillValue, theta_ctl_win)
  theta_exp_win = where(theta_exp_win .le. -1e+35, theta_exp_win@_FillValue, theta_exp_win)
  theta_ctl_spr = where(theta_ctl_spr .le. -1e+35, theta_ctl_spr@_FillValue, theta_ctl_spr)
  theta_exp_spr = where(theta_exp_spr .le. -1e+35, theta_exp_spr@_FillValue, theta_exp_spr)
  theta_ctl_sum = where(theta_ctl_sum .le. -1e+35, theta_ctl_sum@_FillValue, theta_ctl_sum)
  theta_exp_sum = where(theta_exp_sum .le. -1e+35, theta_exp_sum@_FillValue, theta_exp_sum)

  ; 赋予维度和坐标属性
  vars_3d = [/theta_ctl_win, theta_exp_win, theta_ctl_spr, theta_exp_spr, theta_ctl_sum, theta_exp_sum/]
  do i = 0, ListCount(vars_3d)-1
    vars_3d[i]!0 = "south_north"
    vars_3d[i]!1 = "west_east"
    vars_3d[i]!2 = "layer"
  end do

  ; 维度重排
  theta_ctl_win_re = theta_ctl_win(layer|:, south_north|:, west_east|:)
  theta_exp_win_re = theta_exp_win(layer|:, south_north|:, west_east|:)
  theta_ctl_spr_re = theta_ctl_spr(layer|:, south_north|:, west_east|:)
  theta_exp_spr_re = theta_exp_spr(layer|:, south_north|:, west_east|:)
  theta_ctl_sum_re = theta_ctl_sum(layer|:, south_north|:, west_east|:)
  theta_exp_sum_re = theta_exp_sum(layer|:, south_north|:, west_east|:)
  
  printVarSummary(theta_ctl_spr_re)
  printMinMax(theta_ctl_spr_re,0)

  ;===================================================================
  ; 3. [优化点] 初始化 4D 源数据数组 (用于批量插值)
  ;===================================================================
  nlat_src = dimsizes(lat1d)
  nlon_src = dimsizes(lon1d)
  ; 创建一个包含 24 个变量、8 层、加上经纬度的 4D 大数组
  all_vars_src = new((/24, 8, nlat_src, nlon_src/), float)
  all_vars_src@_FillValue = -1.e+36

  ;===================================================================
  ; 4. 循环读取、按主导 PFT 选取、并存入 4D 数组
  ;===================================================================
  ; CTL (nogravel)
  f_pft_ctl = addfile("/share/home/dq135/CRESM_0307S92/SCRIPTS/ICBCScripts/PrepCoLM/CN_15km/CoLMrun_nogravel/unstructured_cwrf_CN_15km_nogravel/landdata/diag/patchfrac_elm.nc", "r")
  pctpfts_ctl = tofloat(f_pft_ctl->patchfrac_elm({8:60},{58:162},:))   

  ; EXP (gravel)
  f_pft_exp = addfile("/share/home/dq135/CRESM_0307S92/SCRIPTS/ICBCScripts/PrepCoLM/CN_15km/CoLMrun_gravel/unstructured_cwrf_CN_15km_gravel/landdata/diag/patchfrac_elm.nc", "r")
  pctpfts_exp = tofloat(f_pft_exp->patchfrac_elm({8:60},{58:162},:))

  ; 基础缺测处理：将缺省值或极小的异常值设为 0
  pctpfts_ctl = where(ismissing(pctpfts_ctl), 0.0, pctpfts_ctl)
  pctpfts_exp = where(ismissing(pctpfts_exp), 0.0, pctpfts_exp)

  ; 剔除特定地表覆盖类型：将 15(永久冰雪) 和 17(水体) 的权重强制设为 0
  pctpfts_ctl(:,:,15) = 0.0
  pctpfts_ctl(:,:,17) = 0.0
  
  pctpfts_exp(:,:,15) = 0.0
  pctpfts_exp(:,:,17) = 0.0

  print("pctpfts_ctl 维度: " + dimsizes(pctpfts_ctl))
  print("pctpfts_exp 维度: " + dimsizes(pctpfts_exp))
  delete([/f_pft_ctl, f_pft_exp/])

  do k=0,7
     ; ====================== 读取参数 3D 数据 ======================
     theta_s_ctl_3d = tofloat(f1->$theta_s(k)$({8:60},{58:162},:))
     theta_s_exp_3d = tofloat(f2->$theta_s(k)$({8:60},{58:162},:))
     tkdry_ctl_3d   = tofloat(f1->$tkdry(k)$({8:60},{58:162},:))
     tkdry_exp_3d   = tofloat(f2->$tkdry(k)$({8:60},{58:162},:))
     tksatu_ctl_3d  = tofloat(f1->$tksatu(k)$({8:60},{58:162},:))
     tksatu_exp_3d  = tofloat(f2->$tksatu(k)$({8:60},{58:162},:))
     tksatf_ctl_3d  = tofloat(f1->$tksatf(k)$({8:60},{58:162},:))
     tksatf_exp_3d  = tofloat(f2->$tksatf(k)$({8:60},{58:162},:))

     BA_alpha_ctl_3d  = tofloat(f1->$BA_alpha(k)$({8:60},{58:162},:))
     BA_alpha_exp_3d  = tofloat(f2->$BA_alpha(k)$({8:60},{58:162},:))
     BA_beta_ctl_3d   = tofloat(f1->$BA_beta(k)$({8:60},{58:162},:))
     BA_beta_exp_3d   = tofloat(f2->$BA_beta(k)$({8:60},{58:162},:))
     vf_sand_s_ctl_3d = tofloat(f1->$vf_sand_s(k)$({8:60},{58:162},:))
     vf_sand_s_exp_3d = tofloat(f2->$vf_sand_s(k)$({8:60},{58:162},:))
     vf_om_s_ctl_3d   = tofloat(f1->$vf_om_s(k)$({8:60},{58:162},:))
     vf_om_s_exp_3d   = tofloat(f2->$vf_om_s(k)$({8:60},{58:162},:))
     vf_gravels_s_ctl_3d = tofloat(f1->$vf_gravels_s(k)$({8:60},{58:162},:))
     vf_gravels_s_exp_3d = tofloat(f2->$vf_gravels_s(k)$({8:60},{58:162},:))
     vf_quartz_ctl_3d = tofloat(f1->$vf_quartz(k)$({8:60},{58:162},:))
     vf_quartz_exp_3d = tofloat(f2->$vf_quartz(k)$({8:60},{58:162},:))
     vf_clay_s_ctl_3d = tofloat(f1->$vf_clay_s(k)$({8:60},{58:162},:))
     vf_clay_s_exp_3d = tofloat(f2->$vf_clay_s(k)$({8:60},{58:162},:))

     csol_ctl_3d = tofloat(f1->$csol(k)$({8:60},{58:162},:))
     csol_exp_3d = tofloat(f2->$csol(k)$({8:60},{58:162},:))

     ; ====================== 加权平均 ======================
     fillv = -1.e+36
     dims2d = (/dimsizes(lat1d), dimsizes(lon1d)/)

     ; ---------------- CTL 加权 ----------------
     sum_frac_ctl = dim_sum_n_Wrap(pctpfts_ctl, 2)
     sum_frac_ctl@_FillValue = fillv  
     sum_frac_ctl = where(sum_frac_ctl .le. 0.0, fillv, sum_frac_ctl)

     theta_s_ctl_avg  = dim_sum_n_Wrap(theta_s_ctl_3d   * pctpfts_ctl, 2) / sum_frac_ctl
     tkdry_ctl_avg    = dim_sum_n_Wrap(tkdry_ctl_3d     * pctpfts_ctl, 2) / sum_frac_ctl
     tksatu_ctl_avg   = dim_sum_n_Wrap(tksatu_ctl_3d    * pctpfts_ctl, 2) / sum_frac_ctl
     tksatf_ctl_avg   = dim_sum_n_Wrap(tksatf_ctl_3d    * pctpfts_ctl, 2) / sum_frac_ctl
     BA_alpha_ctl_avg = dim_sum_n_Wrap(BA_alpha_ctl_3d  * pctpfts_ctl, 2) / sum_frac_ctl
     BA_beta_ctl_avg  = dim_sum_n_Wrap(BA_beta_ctl_3d   * pctpfts_ctl, 2) / sum_frac_ctl
     vf_sand_s_ctl_avg= dim_sum_n_Wrap(vf_sand_s_ctl_3d * pctpfts_ctl, 2) / sum_frac_ctl
     vf_om_s_ctl_avg  = dim_sum_n_Wrap(vf_om_s_ctl_3d   * pctpfts_ctl, 2) / sum_frac_ctl
     vf_gravels_s_ctl_avg= dim_sum_n_Wrap(vf_gravels_s_ctl_3d* pctpfts_ctl, 2) / sum_frac_ctl
     vf_quartz_ctl_avg= dim_sum_n_Wrap(vf_quartz_ctl_3d * pctpfts_ctl, 2) / sum_frac_ctl
     vf_clay_s_ctl_avg= dim_sum_n_Wrap(vf_clay_s_ctl_3d * pctpfts_ctl, 2) / sum_frac_ctl
     csol_ctl_avg     = dim_sum_n_Wrap(csol_ctl_3d      * pctpfts_ctl, 2) / sum_frac_ctl

     ; ---------------- EXP 加权 ----------------
     sum_frac_exp = dim_sum_n_Wrap(pctpfts_exp, 2)
     sum_frac_exp@_FillValue = fillv
     sum_frac_exp = where(sum_frac_exp .le. 0.0, fillv, sum_frac_exp)

     theta_s_exp_avg  = dim_sum_n_Wrap(theta_s_exp_3d   * pctpfts_exp, 2) / sum_frac_exp
     tkdry_exp_avg    = dim_sum_n_Wrap(tkdry_exp_3d     * pctpfts_exp, 2) / sum_frac_exp
     tksatu_exp_avg   = dim_sum_n_Wrap(tksatu_exp_3d    * pctpfts_exp, 2) / sum_frac_exp
     tksatf_exp_avg   = dim_sum_n_Wrap(tksatf_exp_3d    * pctpfts_exp, 2) / sum_frac_exp
     BA_alpha_exp_avg = dim_sum_n_Wrap(BA_alpha_exp_3d  * pctpfts_exp, 2) / sum_frac_exp
     BA_beta_exp_avg  = dim_sum_n_Wrap(BA_beta_exp_3d   * pctpfts_exp, 2) / sum_frac_exp
     vf_sand_s_exp_avg= dim_sum_n_Wrap(vf_sand_s_exp_3d * pctpfts_exp, 2) / sum_frac_exp
     vf_om_s_exp_avg  = dim_sum_n_Wrap(vf_om_s_exp_3d   * pctpfts_exp, 2) / sum_frac_exp
     vf_gravels_s_exp_avg= dim_sum_n_Wrap(vf_gravels_s_exp_3d* pctpfts_exp, 2) / sum_frac_exp
     vf_quartz_exp_avg= dim_sum_n_Wrap(vf_quartz_exp_3d * pctpfts_exp, 2) / sum_frac_exp
     vf_clay_s_exp_avg= dim_sum_n_Wrap(vf_clay_s_exp_3d * pctpfts_exp, 2) / sum_frac_exp
     csol_exp_avg     = dim_sum_n_Wrap(csol_exp_3d      * pctpfts_exp, 2) / sum_frac_exp

     ; ====================== [优化点] 存入 4D 数组 (替代原有的逐个插值) ======================
     all_vars_src(0, k, :, :) = theta_s_ctl_avg
     all_vars_src(1, k, :, :) = theta_s_exp_avg
     all_vars_src(2, k, :, :) = tkdry_ctl_avg
     all_vars_src(3, k, :, :) = tkdry_exp_avg
     all_vars_src(4, k, :, :) = tksatu_ctl_avg
     all_vars_src(5, k, :, :) = tksatu_exp_avg
     all_vars_src(6, k, :, :) = tksatf_ctl_avg
     all_vars_src(7, k, :, :) = tksatf_exp_avg

     all_vars_src(8, k, :, :) = BA_alpha_ctl_avg
     all_vars_src(9, k, :, :) = BA_alpha_exp_avg
     all_vars_src(10, k, :, :) = BA_beta_ctl_avg
     all_vars_src(11, k, :, :) = BA_beta_exp_avg
     all_vars_src(12, k, :, :) = vf_sand_s_ctl_avg
     all_vars_src(13, k, :, :) = vf_sand_s_exp_avg
     all_vars_src(14, k, :, :) = vf_om_s_ctl_avg
     all_vars_src(15, k, :, :) = vf_om_s_exp_avg
     all_vars_src(16, k, :, :) = vf_gravels_s_ctl_avg
     all_vars_src(17, k, :, :) = vf_gravels_s_exp_avg
     all_vars_src(18, k, :, :) = vf_quartz_ctl_avg
     all_vars_src(19, k, :, :) = vf_quartz_exp_avg
     all_vars_src(20, k, :, :) = vf_clay_s_ctl_avg
     all_vars_src(21, k, :, :) = vf_clay_s_exp_avg

     all_vars_src(22, k, :, :) = csol_ctl_avg
     all_vars_src(23, k, :, :) = csol_exp_avg

     ; ====================== 清理内存 ======================
     delete([/theta_s_ctl_3d,theta_s_exp_3d,tkdry_ctl_3d,tkdry_exp_3d,tksatu_ctl_3d,tksatu_exp_3d, \
              tksatf_ctl_3d,tksatf_exp_3d,BA_alpha_ctl_3d,BA_alpha_exp_3d,BA_beta_ctl_3d,BA_beta_exp_3d, \
              vf_sand_s_ctl_3d,vf_sand_s_exp_3d,vf_om_s_ctl_3d,vf_om_s_exp_3d,vf_gravels_s_ctl_3d, \
              vf_gravels_s_exp_3d,vf_quartz_ctl_3d,vf_quartz_exp_3d,vf_clay_s_ctl_3d,vf_clay_s_exp_3d, \
              csol_ctl_3d,csol_exp_3d/])

     delete([/theta_s_ctl_avg,theta_s_exp_avg,tkdry_ctl_avg,tkdry_exp_avg,tksatu_ctl_avg,tksatu_exp_avg, \
              tksatf_ctl_avg,tksatf_exp_avg,BA_alpha_ctl_avg,BA_alpha_exp_avg,BA_beta_ctl_avg,BA_beta_exp_avg, \
              vf_sand_s_ctl_avg,vf_sand_s_exp_avg,vf_om_s_ctl_avg,vf_om_s_exp_avg,vf_gravels_s_ctl_avg, \
              vf_gravels_s_exp_avg,vf_quartz_ctl_avg,vf_quartz_exp_avg,vf_clay_s_ctl_avg,vf_clay_s_exp_avg, \
              csol_ctl_avg,csol_exp_avg/])
  end do

  ;===================================================================
  ; 4.5 [优化核心] 矢量化批量插值与数组解包
  ;===================================================================
  print("开始进行 4D 数组的批量插值...")
  all_vars_interp = rgrid2rcm(lat1d, lon1d, all_vars_src, lat2d, lon2d, 0)
  print("批量插值完成！")

  ; 将插值后的 4D 数组解包回各个独立的 3D 变量 (层数, target_lat, target_lon)
  theta_s_ctl_interp = all_vars_interp(0, :, :, :)
  theta_s_exp_interp = all_vars_interp(1, :, :, :)
  tkdry_ctl_interp   = all_vars_interp(2, :, :, :)
  tkdry_exp_interp   = all_vars_interp(3, :, :, :)
  tksatu_ctl_interp  = all_vars_interp(4, :, :, :)
  tksatu_exp_interp  = all_vars_interp(5, :, :, :)
  tksatf_ctl_interp  = all_vars_interp(6, :, :, :)
  tksatf_exp_interp  = all_vars_interp(7, :, :, :)

  BA_alpha_ctl_interp = all_vars_interp(8, :, :, :)
  BA_alpha_exp_interp = all_vars_interp(9, :, :, :)
  BA_beta_ctl_interp  = all_vars_interp(10, :, :, :)
  BA_beta_exp_interp  = all_vars_interp(11, :, :, :)
  vf_sand_s_ctl_interp= all_vars_interp(12, :, :, :)
  vf_sand_s_exp_interp= all_vars_interp(13, :, :, :)
  vf_om_s_ctl_interp  = all_vars_interp(14, :, :, :)
  vf_om_s_exp_interp  = all_vars_interp(15, :, :, :)
  vf_gravels_s_ctl_interp= all_vars_interp(16, :, :, :)
  vf_gravels_s_exp_interp= all_vars_interp(17, :, :, :)
  vf_quartz_ctl_interp= all_vars_interp(18, :, :, :)
  vf_quartz_exp_interp= all_vars_interp(19, :, :, :)
  vf_clay_s_ctl_interp= all_vars_interp(20, :, :, :)
  vf_clay_s_exp_interp= all_vars_interp(21, :, :, :)

  csol_ctl_interp = all_vars_interp(22, :, :, :)
  csol_exp_interp = all_vars_interp(23, :, :, :)

  ;===================================================================
  ; 5. 复制坐标属性并准备输出
  ;===================================================================
  interp_vars = [/theta_s_ctl_interp, theta_s_exp_interp, tkdry_ctl_interp, tkdry_exp_interp, \
                 tksatu_ctl_interp, tksatu_exp_interp, tksatf_ctl_interp, tksatf_exp_interp, \
                 BA_alpha_ctl_interp, BA_alpha_exp_interp, BA_beta_ctl_interp, BA_beta_exp_interp, \
                 vf_sand_s_ctl_interp, vf_sand_s_exp_interp, vf_om_s_ctl_interp, vf_om_s_exp_interp, \
                 vf_gravels_s_ctl_interp, vf_gravels_s_exp_interp, vf_quartz_ctl_interp, vf_quartz_exp_interp, \
                 vf_clay_s_ctl_interp, vf_clay_s_exp_interp, csol_ctl_interp, csol_exp_interp/]
                 
  do i = 0, ListCount(interp_vars)-1
    interp_vars[i]!0 = "level"
    interp_vars[i]!1 = "south_north"
    interp_vars[i]!2 = "west_east"
    interp_vars[i]@lat2d = lat2d
    ; interp_vars[i]@lon2d = lon2d
  end do

  ;===================================================================
  ; 6. 输出到新的 NetCDF 数据文件
  ;===================================================================
  out_filename = "/share/home/dq135/draw/scripts/offline/interp_soil_params.nc"
  system("rm -f " + out_filename)  
  fout = addfile(out_filename, "c")
  
  fout@description = "Interpolated soil parameters to CWRF Grid"
  
  ; 写入二维经纬度坐标 
  fout->XLAT  = lat2d
  fout->XLONG = lon2d

  ; 写入计算过的重排变量 (新增冬季数据)
  fout->theta_ctl_win_re = theta_ctl_win_re
  fout->theta_exp_win_re = theta_exp_win_re
  fout->theta_ctl_spr_re = theta_ctl_spr_re
  fout->theta_exp_spr_re = theta_exp_spr_re
  fout->theta_ctl_sum_re = theta_ctl_sum_re
  fout->theta_exp_sum_re = theta_exp_sum_re

  ; 写入插值后的变量
  fout->theta_s_ctl_interp = theta_s_ctl_interp
  fout->theta_s_exp_interp = theta_s_exp_interp
  fout->tkdry_ctl_interp   = tkdry_ctl_interp
  fout->tkdry_exp_interp   = tkdry_exp_interp
  fout->tksatu_ctl_interp  = tksatu_ctl_interp
  fout->tksatu_exp_interp  = tksatu_exp_interp
  fout->tksatf_ctl_interp  = tksatf_ctl_interp
  fout->tksatf_exp_interp  = tksatf_exp_interp
  
  fout->BA_alpha_ctl_interp     = BA_alpha_ctl_interp
  fout->BA_alpha_exp_interp     = BA_alpha_exp_interp
  fout->BA_beta_ctl_interp      = BA_beta_ctl_interp
  fout->BA_beta_exp_interp      = BA_beta_exp_interp
  fout->vf_sand_s_ctl_interp    = vf_sand_s_ctl_interp
  fout->vf_sand_s_exp_interp    = vf_sand_s_exp_interp
  fout->vf_om_s_ctl_interp      = vf_om_s_ctl_interp
  fout->vf_om_s_exp_interp      = vf_om_s_exp_interp
  fout->vf_gravels_s_ctl_interp = vf_gravels_s_ctl_interp
  fout->vf_gravels_s_exp_interp = vf_gravels_s_exp_interp
  fout->vf_quartz_ctl_interp    = vf_quartz_ctl_interp
  fout->vf_quartz_exp_interp    = vf_quartz_exp_interp
  fout->vf_clay_s_ctl_interp    = vf_clay_s_ctl_interp
  fout->vf_clay_s_exp_interp    = vf_clay_s_exp_interp

  fout->csol_ctl_interp         = csol_ctl_interp
  fout->csol_exp_interp         = csol_exp_interp

  print("====== Data successfully written to: " + out_filename + " ======")

end