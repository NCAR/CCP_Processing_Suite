;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "attributes_ccsm.ncl"
;******************************************************************************

begin
  n1 = integertochar(10)
;--------------------------------------------------------------------
; Get input file env vars
;--------------------------------------------------------------------
  ncin       = getenv("NCL_FILE")
  case       = getenv("CASE")
  type       = getenv("TYPE")
  tval       = getenv("NCL_TVAL")
  nclscript  = getenv("NCL_SCRIPT")
  nclversion = getenv("NCL_VERSION")
  prochost   = getenv("PROCHOST") 
;--------------------------------------------------------------------
; Define global attributes
;--------------------------------------------------------------------
  glbAtt = True
;--------------------------------------------------------------------
; Set various attributes
;--------------------------------------------------------------------
  forcing = ""
  expid   = ""
  attributes_ccsm(case,glbAtt,expid)
;
  table   = "O1"
  glbAtt@table_id = "Table "+table
;--------------------------------------------------------------------
; Define table and ipcc_vars
;--------------------------------------------------------------------
  table     = "O1"
  ipcc_vars = (/"hfogo","stfmmc"/)
;--------------------------------------------------------------------
; Read input file
;--------------------------------------------------------------------
  in_file = addfile (ncin,"r")
  nbnd = 2
;--------------------------------------------------------------------
; Read input file X and T coordinates
;--------------------------------------------------------------------
  time  = in_file->time
  if (isfilevar(in_file,"lat_aux_grid")) then
    lat  = (/in_file->lat_aux_grid/)
  else
    lat  = (/in_file->lat/)
  end if
  ntim  = dimsizes(time)
  nlat  = dimsizes(lat)
;--------------------------------------------------------------------
; Calculate time bounds
; Note that bound(X,1) is the same as bound(X-1,2)
;--------------------------------------------------------------------
; Month beg, month end, in terms of day-of-year
;-----------------------------------------------------------------
  mbeg = (/  1, 32, 60, 91,121,152,182,213,244,274,305,335/)
  mend = (/ 31, 59, 90,120,151,181,212,243,273,304,334,365/)
  time_bnds = new((/ntim,nbnd/),double)
  tnew      = new((/ntim/),double)
  fyr       = ((time-15)/365)
  iyr       = floor(fyr)
  rmn       = fyr - iyr
  imn       = floor(rmn*365)
  jmn       = doubletointeger(floor((imn/(365/12))))
  time_bnds(0:ntim-1,0) = ((iyr*365)+mbeg(jmn)-1)
  time_bnds(0:ntim-1,1) = ( iyr*365)+mend(jmn)
  tnew = (time_bnds(:,0)+time_bnds(:,1))/2
  delete(tnew@_FillValue)
  delete(time_bnds@_FillValue)
;--------------------------------------------------------------------
; Read input field(s), using IPCC required fields as parser
;--------------------------------------------------------------------
  nvars = dimsizes(ipcc_vars)
  do v = 0,nvars - 1
    if (ipcc_vars(v).eq."hfogo") then
      if (isfilevar(in_file,"N_HEAT")) then
;
; In-POP calculation
;
; assumes structure of:
; N_HEAT(time,transport_reg,transport_comp,lat_aux_grid)
; 
; transport_reg(0) = GLOBAL OCEANS
; transport_reg(1) = ATLANTIC
; 
; transport_comp(0)      = Total
; transport_comp(1)      = Eulerian-Mean Advection
; transport_comp(2)      = Eddy-Induced Advection (bolus) + Diffusion
        nht_total  = (/in_file->N_HEAT(:,0,0,:)/)
        bolus_diff = (/in_file->N_HEAT(:,0,2,:)/)

;  fix the sign of the bolus + diffusion term 
        bolus_diff = (-1) * bolus_diff 

;  remove erroneous bolus+diffusion contribution to total
        nht_total  = nht_total + (2*bolus_diff)

        nregion           = 4
        raw                = new((/ntim,nregion,nlat/),float)
        raw(:,0,:)         = missing_value
        raw(:,1,:)         = missing_value
        raw(:,2,:)         = missing_value
        raw(:,3,:)         = (/nht_total/)*1.e15
        raw@standard_name  = "northward_ocean_heat_transport"
        raw@original_name  = "N_HEAT"
        raw@original_units = "PW"
        raw@units          = "W"
        raw@comment        = "Only global_ocean"
        raw@history        = "Scaled by 1.e^15"
        raw@cell_methods   = "time: mean (interval: 1 month)"
      end if

      if (isfilevar(in_file,"gmht")) then
        global   = (/in_file->gmht/)
        atlantic = (/in_file->amht/)
;
; Offline calculation
;
        nregion           = 4
        raw                = new((/ntim,nregion,nlat/),float)
        raw(:,0,:)         = (/atlantic/)*1.e15
        raw(:,1,:)         = missing_value
        raw(:,2,:)         = missing_value
        raw(:,3,:)         = (/global/)*1.e15
        raw@standard_name  = "northward_ocean_heat_transport"
        raw@original_name  = "amht, gmht"
        raw@original_units = "PW"
        raw@units          = "W"
        raw@comment        = "Offline calculation; eddy component not available"
        raw@history        = "Scaled by 1.e^15"
        raw@cell_methods   = "time: annual (interval: 1 year)"
      end if
    end if

    if (ipcc_vars(v).eq."stfmmc") then
      if (isfilevar(in_file,"MOC")) then
        depth = (/in_file->moc_z/)/100
        ndep  = dimsizes(depth)
        global_moc   = in_file->MOC(:,0,:,:,:)
        atlantic_moc = global_moc      ; initialization is important
        indopac_moc  = global_moc      ; for this case
        atlantic_moc = in_file->MOC(:,1,:,:,:)
;======================================================================
; filter out erroneous values and topographically mask
;======================================================================
; reduce to one dimension for ind function)
        global_moc1d   = ndtooned(global_moc)
        atlantic_moc1d = ndtooned(atlantic_moc)

; remove near zero numbers (erroneous)
        global_moc1d(ind(global_moc1d.gt.-10e-5 .and. global_moc1d.lt.10e-5)) \
                    = global_moc@_FillValue
        atlantic_moc1d(ind(atlantic_moc1d.gt.-10e-6 .and. atlantic_moc1d.lt.10e-6)) \
                    = atlantic_moc@_FillValue

; expand to multiple dimensions
        global_moc   = onedtond(global_moc1d,dimsizes(global_moc))
        atlantic_moc = onedtond(atlantic_moc1d,dimsizes(atlantic_moc))

; cut off all values from South Pole to index 83 (missing vals)
        atlantic_moc(:,:,:,0:83) = atlantic_moc@_FillValue

;======================================================================
; Calculate Indo-Pacific Basin MOC
;======================================================================
;        indopac_moc = global_moc - atlantic_moc
;        indopac_moc(:,:,:,349::) = indopac_moc@_FillValue  ; no MOC contribution 
        					           ; from points north
						           ; of Bering Strait or
        nregion           = 4
        raw               = new((/ntim,nregion,ndep,nlat/),float)
        raw(:,0,:,:)      = (/atlantic_moc/)*1.e6
        raw(:,1,:,:)      = missing_value
        raw(:,2,:,:)      = missing_value
        raw(:,3,:,:)      = (/global_moc/)*1.e6
        raw@standard_name  = "ocean_meridional_overturning_streamfunction"
        raw@original_name  = "MOC"
        raw@original_units = "Sv"
        raw@units          = "m 3 s-1"
        raw@history        = "Scaled by 1.e^6"
        raw@comment        = "Only global_ocean and atlantic_ocean"
        raw@cell_methods   = "time: mean (interval: 1 month)"
      end if
    end if
;--------------------------------------------------------------------
; Field-invariant metadata
;--------------------------------------------------------------------
    if (isvar("raw")) then
     raw@long_name      = raw@standard_name
     raw@_FillValue     = 1.e+20
     raw@missing_value  = 1.e+20
     raw@comment        = "Created using NCL version "+nclversion+n1+\
                          " NCL script: "+nclscript+" on "+n1+\
                          "processing host "+prochost+" by "+systemfunc("whoami")+" at "+n1+\
                          systemfunc("date")
;--------------------------------------------------------------------
; Create output file
; [file name as per requirements document]
;--------------------------------------------------------------------
     nccf = "./IPCC/"+ipcc_vars(v)+"_"+table+"."+expid+".CCSM."+type+"."+tval+".nc"
     system("csh -c 'rm -f "+nccf+"'")
     cf_file = addfile(nccf, "c")
     fileattdef(cf_file, glbAtt)

     if (ipcc_vars(v).eq."hfogo") then
;--------------------------------------------------------------------
; Calculate Y bounds
; Note that bound(X,1) is the same as bound(X-1,2)
;-----------------------------------------------------------------
       lat_bnds  = new((/nlat,nbnd/),double)
       lat_bnds(0,0)      = lat(0)-((lat(1)-lat(0))/2.)
       lat_bnds(nlat-1,1) =  90.
       do j = 0,nlat-2
         lat_bnds(j,1) = ((lat(j)+lat(j+1))/2.)
       end do
       do j = 1,nlat-1
         lat_bnds(j,0) =  lat_bnds(j-1,1)
       end do
       delete(lat_bnds@_FillValue)
       lat@bounds    = "lat_bnds"

       var_y           = (/ "lat" /)
       var_y_type      = (/ "double" /)
       var_y_bnds      = (/ "lat_bnds" /)
       var_y_bnds_type = (/ "double" /)
       var_r           = (/ "region" /)
       var_r_type      = (/ "character" /)
       var_t_bnds      = (/ "time_bnds" /)
       var_t_bnds_type = (/ "double" /)
       var_t           = (/ "time" /)
       var_t_type      = (/ "double" /)

       dimNames = (/ "region", "strlen", "lat", "bnds", "time" /)
       dimSizes = (/ nregion,    stglen,  nlat,   nbnd,     -1 /)
       dimUnlim = (/    False,    False, False,  False,   True /)
       filedimdef(cf_file, dimNames, dimSizes, dimUnlim )
       filevardef(cf_file, var_r, var_r_type, (/"region","strlen"/))
       filevardef(cf_file, var_y, var_y_type, (/"lat"/))
       filevardef(cf_file, var_t, var_t_type, (/"time"/))
       filevardef(cf_file, var_y_bnds, var_y_bnds_type, (/"lat","bnds"/) )
       filevardef(cf_file, var_t_bnds, var_t_bnds_type, (/"time","bnds"/) )
;--------------------------------------------------------------------
; Write coordinates and their bounds to output file
;--------------------------------------------------------------------
       regions = (/"atlantic_ocean","indian_ocean  ","pacific_ocean ","global_ocean"/)
       region = stringtocharacter(regions)
       region!0 = "region"
       region!1 = "strlen"
       region@standard_name = "region"
       lat@axis             = "Y"
       lat@units            = "degrees_north"
       lat@standard_name    = "latitude"
       
       raw!0              = "time"
       raw!1              = "region"
       raw!2              = "lat"
       raw&lat            = lat
       raw&time           = tnew
       tnew@bounds        = "time_bnds"
       tnew@units         = "days since 0000-1-1"
       tnew@axis          = "T"
       tnew@standard_name = "time"
       tnew@calendar      = "noleap"

       cf_file->region     = region
       cf_file->lat        = lat
       cf_file->lat_bnds   = lat_bnds
       cf_file->time       = tnew
       cf_file->time_bnds  = time_bnds
       cf_file->$ipcc_vars(v)$ = raw
     end if

     if (ipcc_vars(v).eq."stfmmc") then
;--------------------------------------------------------------------
; Calculate Y bounds
; Note that bound(X,1) is the same as bound(X-1,2)
;-----------------------------------------------------------------
       lat_bnds  = new((/nlat,nbnd/),double)
       lat_bnds(0,0)      = lat(0)-((lat(1)-lat(0))/2.)
       lat_bnds(nlat-1,1) =  90.
       do j = 0,nlat-2
         lat_bnds(j,1) = ((lat(j)+lat(j+1))/2.)
       end do
       do j = 1,nlat-1
         lat_bnds(j,0) =  lat_bnds(j-1,1)
       end do
       delete(lat_bnds@_FillValue)
       lat@bounds    = "lat_bnds"

       depth = (/in_file->moc_z/)/100
       ndep  = dimsizes(depth)
       depth_bnds  = new((/ndep,nbnd/),double)
       depth_bnds(:,0) = (/0 ,5.00622,  15.06873,  25.28343,  35.75849,  46.6127, \
                             57.98098,  70.02135,  82.92405,  96.92411, 112.3189, \
                            129.4936 , 148.9582 , 171.4044 , 197.7919 , 229.4842, \
                            268.4617 , 317.6501 , 381.3864 , 465.9132 , 579.3073, \
                            729.3514 , 918.3725 ,1139.154  ,1378.574  ,1625.701 , \
                           1875.107  ,2125.012  ,2375      ,2624.999  ,2874.999 , \
                           3124.999  ,3374.999  ,3624.999  ,3874.999  ,4125     , \
                           4375      ,4625      ,4874.999  ,5124.999  ,5374.999 /)
       depth_bnds(:,1) = (/   5.00622,  15.06873,  25.28343,  35.75849,  46.6127, \
                             57.98098,  70.02135,  82.92405,  96.92411, 112.3189, \
                            129.4936 , 148.9582 , 171.4044 , 197.7919 , 229.4842, \
                            268.4617 , 317.6501 , 381.3864 , 465.9132 , 579.3073, \
                            729.3514 , 918.3725 ,1139.154  ,1378.574  ,1625.701 , \
                           1875.107  ,2125.012  ,2375      ,2624.999  ,2874.999 , \
                           3124.999  ,3374.999  ,3624.999  ,3874.999  ,4125     , \
                           4375      ,4625      ,4874.999  ,5124.999  ,5374.999 , \
                           5750      /)
       delete(depth_bnds@_FillValue)

       var_y           = (/ "lat" /)
       var_y_type      = (/ "double" /)
       var_y_bnds      = (/ "lat_bnds" /)
       var_y_bnds_type = (/ "double" /)
       var_z           = (/ "depth" /)
       var_z_type      = (/ "double" /)
       var_z_bnds      = (/ "depth_bnds" /)
       var_z_bnds_type = (/ "double" /)
       var_r           = (/ "region" /)
       var_r_type      = (/ "character" /)
       var_t_bnds      = (/ "time_bnds" /)
       var_t_bnds_type = (/ "double" /)
       var_t           = (/ "time" /)
       var_t_type      = (/ "double" /)

       dimNames = (/ "region", "strlen", "lat", "depth", "bnds", "time" /)
       dimSizes = (/ nregion,    stglen,  nlat,    ndep,   nbnd,     -1 /)
       dimUnlim = (/    False,    False, False,   False,  False,   True /)
       filedimdef(cf_file, dimNames, dimSizes, dimUnlim )
       filevardef(cf_file, var_r, var_r_type, (/"region","strlen"/))
       filevardef(cf_file, var_y, var_y_type, (/"lat"/))
       filevardef(cf_file, var_z, var_z_type, (/"depth"/))
       filevardef(cf_file, var_t, var_t_type, (/"time"/))
       filevardef(cf_file, var_y_bnds, var_y_bnds_type, (/"lat","bnds"/) )
       filevardef(cf_file, var_z_bnds, var_z_bnds_type, (/"depth","bnds"/) )
       filevardef(cf_file, var_t_bnds, var_t_bnds_type, (/"time","bnds"/) )
;--------------------------------------------------------------------
; Write coordinates and their bounds to output file
;--------------------------------------------------------------------
       regions = (/"atlantic_ocean","indian_ocean  ","pacific_ocean ","global_ocean"/)
       region = stringtocharacter(regions)
       region!0 = "region"
       region!1 = "strlen"
       region@standard_name = "region"
       lat@units            = "degrees_north"
       lat@standard_name    = "latitude"
       depth@bounds         = "depth_bnds"
       depth@units          = "m"
       depth@positive       = "down"
       depth@standard_name  = "depth"
       depth@units          = "m"
       depth@axis           = "Z"
       
       raw!0        = "time"
       raw!1        = "region"
       raw!2        = "depth"
       raw!3        = "lat"
       raw&lat      = lat
       raw&depth    = depth
       raw&time     = tnew
       tnew@bounds        = "time_bnds"
       tnew@units         = "days since 0000-1-1"
       tnew@axis          = "T"
       tnew@standard_name = "time"
       tnew@calendar      = "noleap"

       cf_file->region = region
       cf_file->lat        = lat
       cf_file->lat_bnds   = lat_bnds
       cf_file->depth_bnds = depth_bnds
       cf_file->time       = tnew
       cf_file->time_bnds  = time_bnds
       cf_file->$ipcc_vars(v)$ = raw
     end if
     delete(raw)
     delete(dimNames)
     delete(dimSizes)
     delete(dimUnlim)
    end if
  end do
end
