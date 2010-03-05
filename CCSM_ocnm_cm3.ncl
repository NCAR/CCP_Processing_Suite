;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "attributes_ccsm.ncl"
;******************************************************************************

begin
;--------------------------------------------------------------------
; Constants
;--------------------------------------------------------------------
; Density of water, 1000 kg m-3
  density_h2o = 1000
; grams to kilograms
  g_to_kg     = 1000.
; Stefan-Boltzmann constant, W m-2 K-4
  sb          = 5.67E-8 

  n1 = integertochar(10)
;--------------------------------------------------------------------
; Get input file env vars
;--------------------------------------------------------------------
  ncin       = getenv("NCL_FILE")
  case       = getenv("CASE")
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
;--------------------------------------------------------------------
; Define table and ipcc_vars
;--------------------------------------------------------------------
  table           = "O1"
  glbAtt@table_id = "Table "+table
;
  ipcc_vars = (/"so","thetao","rhopoto","uo","vo","wo","tos","zos"/)
;--------------------------------------------------------------------
; Read input file
;--------------------------------------------------------------------
  in_file = addfile (ncin,"r")
  nbnd = 2
;--------------------------------------------------------------------
; Read input file X and T coordinates
;--------------------------------------------------------------------
  time = in_file->time
  ntim = dimsizes(time)
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
  fyr       = ((time-1)/365)
  iyr       = floor(fyr)
  rmn       = fyr - iyr
  imn       = floor(rmn*365)-15
  jmn       = floattoint(floor((imn/(365/12))))
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
    if (ipcc_vars(v).eq."zos") then
      if (isfilevar(in_file,"SSH")) then
        raw                = (/in_file->SSH/)/100
        raw@standard_name  = "sea_surface_height_above_geoid"
        raw@original_name  = "SSH"
        raw@original_units = "cm"
        raw@units          = "m"
        raw@history        = "Interpolated to regular grid from dipole grid"+n1+\
                             "SSH/100"
        raw@cell_methods   = "time: mean (interval: 1 month)"
      end if
    end if
    if (ipcc_vars(v).eq."tos") then
      if (isfilevar(in_file,"TEMP")) then
        raw                = (/in_file->TEMP(:,0,:,:)/)+273.15
        raw@standard_name  = "sea_surface_temperature"
        raw@original_name  = "TEMP"
        raw@original_units = "C"
        raw@units          = "K"
        raw@history        = "Interpolated to regular grid from dipole grid,"+n1+\
                             "TEMP+273.15"
        raw@cell_methods   = "time: mean (interval: 1 month)"
      end if
    end if
    if (ipcc_vars(v).eq."so") then
      if (isfilevar(in_file,"SALT")) then
        raw                = (/in_file->SALT/)*1000.
        raw@standard_name  = "sea_water_salinity"
        raw@original_name  = "SALT"
        raw@original_units = raw@units
        raw@units          = "1e-3"
        raw@history        = "Interpolated to regular grid from dipole grid"
        raw@cell_methods   = "time: mean (interval: 1 month)"
      end if
    end if
    if (ipcc_vars(v).eq."thetao") then
      if (isfilevar(in_file,"TEMP")) then
        raw                = (/in_file->TEMP/)+273.15
        raw@standard_name  = "sea_water_potential_temperature"
        raw@original_name  = "TEMP"
        raw@original_units = "C"
        raw@units          = "K"
        raw@history        = "Interpolated to regular grid from dipole grid,"+n1+\
                             "TEMP+273.15"
        raw@cell_methods   = "time: mean (interval: 1 month)"
      end if
    end if
    if (ipcc_vars(v).eq."rhopoto") then
      if (isfilevar(in_file,"PD")) then
        raw                =  (/in_file->PD/)*1.e3
        raw@standard_name  = "sea_water_potential_density"
        raw@original_name  = "PD"
        raw@original_units = "g cm-3"
        raw@units          = "kg m-3"
        raw@history        = "Interpolated to regular grid from dipole grid,"+n1+\
                             "PD*1.e3"
        raw@cell_methods   = "time: mean (interval: 1 month)"
      end if
    end if
    if (ipcc_vars(v).eq."uo") then
      if (isfilevar(in_file,"UVEL")) then
        raw                =  (/in_file->UVEL/)/100.
        raw@standard_name  = "eastward_sea_water_velocity"
        raw@original_name  = "UVEL"
        raw@original_units = "cm s-1"
        raw@units          = "m s-1"
        raw@history        = "Interpolated and rotated to regular grid from dipole grid,"+n1+\
                             "UVEL/100"
        raw@cell_methods   = "time: mean (interval: 1 month)"
      end if
    end if
    if (ipcc_vars(v).eq."vo") then
      if (isfilevar(in_file,"VVEL")) then
        raw                =  (/in_file->VVEL/)/100.
        raw@standard_name  = "northward_sea_water_velocity"
        raw@original_name  = "VVEL"
        raw@original_units = "cm s-1"
        raw@units          = "m s-1"
        raw@history        = "Interpolated and rotated to regular grid from dipole grid,"+n1+\
                             "VVEL/100"
        raw@cell_methods   = "time: mean (interval: 1 month)"
      end if
    end if
    if (ipcc_vars(v).eq."wo") then
      if (isfilevar(in_file,"WVEL")) then
        raw                =  (/in_file->WVEL/)/100.
        raw@standard_name  = "upward_sea_water_velocity"
        raw@original_name  = "WVEL"
        raw@original_units = "cm s-1"
        raw@units          = "m s-1"
        raw@history        = "Interpolated and rotated to regular grid from dipole grid,"+n1+\
                             "WVEL/100"
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
     nccf = "./CMIP/"+ipcc_vars(v)+"_"+table+"."+expid+".CCSM.ocnm."+tval+".nc"
     system("csh -c 'rm -f "+nccf+"'")
     cf_file = addfile(nccf, "c")
     fileattdef(cf_file, glbAtt)
;--------------------------------------------------------------------
; Setups for output file
;--------------------------------------------------------------------
; Coordinates
;--------------------------------------------------------------------
      if ((ipcc_vars(v).eq."so").or.\
          (ipcc_vars(v).eq."thetao").or.\
          (ipcc_vars(v).eq."rhopoto").or.\
          (ipcc_vars(v).eq."uo").or.\
          (ipcc_vars(v).eq."vo").or.\
          (ipcc_vars(v).eq."wo")) then
;--------------------------------------------------------------------
; Read input file coordinates
;--------------------------------------------------------------------
        nlon  = dimsizes(in_file->lon)
        nlat  = dimsizes(in_file->lat)
        lat   = new((/nlat/),double)
        lat   = in_file->lat
        lon   = new((/nlon/),double)
        lon   = in_file->lon
        delete(lat@_FillValue)
        delete(lon@_FillValue)
        if (ipcc_vars(v) .eq. "wo") then
          depth = (/in_file->z_w/)
        else
          depth = (/in_file->z_t/)
        end if
        nlon  = dimsizes(lon)
        nlat  = dimsizes(lat)
        ndep  = dimsizes(depth)
;--------------------------------------------------------------------
; Calculate X/Y/Z bounds
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
        lon_bnds  = new((/nlon,nbnd/),double)
        lon_bnds(0,     0) = (lon(0) - ((lon(1)-lon(0))/2.))
        lon_bnds(nlon-1,1) = (lon(nlon-1) + ((lon(nlon-1)-lon(nlon-2))/2.))
        do i = 0,nlon-2
          lon_bnds(i,1) = ((lon(i)+lon(i+1))/2.)
        end do
        do i = 1,nlon-1
          lon_bnds(i,0) =  lon_bnds(i-1,1)
        end do
        depth_bnds = new((/ndep,nbnd/),double)
        if (ipcc_vars(v) .eq. "wo") then
          depth_bnds(:,0) = (/0 ,5.00622,  15.06873,  25.28343,  35.75849,  46.6127, \
                                57.98098,  70.02135,  82.92405,  96.92411, 112.3189, \
                               129.4936 , 148.9582 , 171.4044 , 197.7919 , 229.4842, \
                               268.4617 , 317.6501 , 381.3864 , 465.9132 , 579.3073, \
                               729.3514 , 918.3725 ,1139.154  ,1378.574  ,1625.701 , \
                              1875.107  ,2125.012  ,2375      ,2624.999  ,2874.999 , \
                              3124.999  ,3374.999  ,3624.999  ,3874.999  ,4125     , \
                              4375      ,4625      ,4874.999  ,5124.999  /)
          depth_bnds(:,1) = (/   5.00622,  15.06873,  25.28343,  35.75849,  46.6127, \
                                57.98098,  70.02135,  82.92405,  96.92411, 112.3189, \
                               129.4936 , 148.9582 , 171.4044 , 197.7919 , 229.4842, \
                               268.4617 , 317.6501 , 381.3864 , 465.9132 , 579.3073, \
                               729.3514 , 918.3725 ,1139.154  ,1378.574  ,1625.701 , \
                              1875.107  ,2125.012  ,2375      ,2624.999  ,2874.999 , \
                              3124.999  ,3374.999  ,3624.999  ,3874.999  ,4125     , \
                              4375      ,4625      ,4874.999  ,5124.999  ,5374.999 /)
        else
          depth_bnds(:,0) = (/        0 ,  10.0124 ,  20.125  ,  30.4418 ,  41.0751, \
                                52.1503 ,  63.8117 ,  76.231  ,  89.6171 , 104.231 , \
                               120.407  , 138.58   , 159.336  , 183.473  , 212.111 , \
                               246.857  , 290.066  , 345.234  , 417.539  , 514.288 , \
                               644.327  , 814.376  ,1022.37   ,1255.94   ,1501.21  , \
                              1750.19   ,2000.02   ,2250.00   ,2500.00   ,2750.00  , \
                              3000.00   ,3250.00   ,3500.00   ,3750.00   ,4000.00  , \
                              4250.00   ,4500.00   ,4750.00   ,5000.00   ,5250.00  /)
          depth_bnds(:,1) = (/             10.0124 ,  20.125  ,  30.4418 ,  41.0751, \
                                52.1503 ,  63.8117 ,  76.231  ,  89.6171 , 104.231 , \
                               120.407  , 138.58   , 159.336  , 183.473  , 212.111 , \
                               246.857  , 290.066  , 345.234  , 417.539  , 514.288 , \
                               644.327  , 814.376  ,1022.37   ,1255.94   ,1501.21  , \
                              1750.19   ,2000.02   ,2250.00   ,2500.00   ,2750.00  , \
                              3000.00   ,3250.00   ,3500.00   ,3750.00   ,4000.00  , \
                              4250.00   ,4500.00   ,4750.00   ,5000.00   ,5250.00  , \
                              5374.999  /)
        end if

        delete(lat_bnds@_FillValue)
        delete(lon_bnds@_FillValue)
        delete(depth_bnds@_FillValue)
        lat@bounds    = "lat_bnds"
        lon@bounds    = "lon_bnds"
        depth@bounds  = "depth_bnds"

        var_x           = (/ "lon" /)
        var_x_type      = (/ "double" /)
        var_x_bnds      = (/ "lon_bnds" /)
        var_x_bnds_type = (/ "double" /)
        var_y           = (/ "lat" /)
        var_y_type      = (/ "double" /)
        var_y_bnds      = (/ "lat_bnds" /)
        var_y_bnds_type = (/ "double" /)
        var_z           = (/ "depth" /)
        var_z_type      = (/ "double" /)
        var_z_bnds_type = (/ "double" /)
        var_z_bnds      = (/ "depth_bnds" /)
        var_t_bnds      = (/ "time_bnds" /)
        var_t_bnds_type = (/ "double" /)
        var_t           = (/ "time" /)
        var_t_type      = (/ "double" /)

        dimNames = (/ "lon", "lat", "depth", "bnds", "time" /)
        dimSizes = (/  nlon,  nlat, ndep ,  nbnd, -1      /)
        dimUnlim = (/ False, False,  False, False, True    /)
        filedimdef(cf_file, dimNames, dimSizes, dimUnlim )
        filevardef(cf_file, var_x_bnds, var_x_bnds_type, (/"lon","bnds"/) )
        filevardef(cf_file, var_y_bnds, var_y_bnds_type, (/"lat","bnds"/) )
        filevardef(cf_file, var_z_bnds, var_z_bnds_type, (/"depth","bnds"/) )
        filevardef(cf_file, var_t_bnds, var_t_bnds_type, (/"time","bnds"/) )
        filevardef(cf_file, var_x, var_x_type, (/"lon" /))
        filevardef(cf_file, var_y, var_y_type, (/"lat" /))
        filevardef(cf_file, var_z, var_z_type, (/"depth"/))
        filevardef(cf_file, var_t, var_t_type, (/"time"/))
;--------------------------------------------------------------------
; Write coordinates and their bounds to output file
;--------------------------------------------------------------------
        lat@standard_name    = "latitude"
        lat@axis             = "Y"
        lon@standard_name    = "longitude"
        lon@axis             = "X"
        depth@bounds         = "depth_bnds"
        depth@units          = "m"
        depth@positive       = "down"
        depth@standard_name  = "depth"
        depth@units          = "m"
        depth@axis           = "Z"
;
        raw!0        = "time"
        raw&time     = tnew
        tnew@bounds        = "time_bnds"
        tnew@units         = "days since 0000-1-1"
        tnew@axis          = "T"
        tnew@standard_name = "time"
        tnew@calendar      = "noleap"
        raw!1        = "depth"
        raw&depth    = depth
        raw!2        = "lat"
        raw&lat      = lat
        raw!3        = "lon"
        raw&lon      = lon

        cf_file->lat = lat
        cf_file->lat_bnds = lat_bnds
        cf_file->lon = lon
        cf_file->lon_bnds = lon_bnds
        cf_file->depth      = depth
        cf_file->depth_bnds = depth_bnds
        cf_file->$ipcc_vars(v)$ = raw
      end if
      if ((ipcc_vars(v).eq."zos").or.\
          (ipcc_vars(v).eq."tos")) then
;--------------------------------------------------------------------
; Read input file coordinates
;--------------------------------------------------------------------
        nlon  = dimsizes(in_file->lon)
        nlat  = dimsizes(in_file->lat)
        lat   = new((/nlat/),double)
        lat   = in_file->lat
        lon   = new((/nlon/),double)
        lon   = in_file->lon
;--------------------------------------------------------------------
; Calculate X/Y bounds
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
        lon_bnds  = new((/nlon,nbnd/),double)
        lon_bnds(0,     0) = (lon(0) - ((lon(1)-lon(0))/2.))
        lon_bnds(nlon-1,1) = (lon(nlon-1) + ((lon(nlon-1)-lon(nlon-2))/2.))
        do i = 0,nlon-2
          lon_bnds(i,1) = ((lon(i)+lon(i+1))/2.)
        end do
        do i = 1,nlon-1
          lon_bnds(i,0) =  lon_bnds(i-1,1)
        end do
        delete(lat_bnds@_FillValue)
        delete(lon_bnds@_FillValue)
        lat@bounds    = "lat_bnds"
        lon@bounds    = "lon_bnds"

        var_x           = (/ "lon" /)
        var_x_type      = (/ "double" /)
        var_x_bnds      = (/ "lon_bnds" /)
        var_x_bnds_type = (/ "double" /)
        var_y           = (/ "lat" /)
        var_y_type      = (/ "double" /)
        var_y_bnds      = (/ "lat_bnds" /)
        var_y_bnds_type = (/ "double" /)
        var_t_bnds      = (/ "time_bnds" /)
        var_t_bnds_type = (/ "double" /)
        var_t           = (/ "time" /)
        var_t_type      = (/ "double" /)

        dimNames = (/ "lon", "lat", "bnds", "time" /)
        dimSizes = (/  nlon,  nlat,   nbnd, -1     /)
        dimUnlim = (/ False, False,  False, True    /)
        filedimdef(cf_file, dimNames, dimSizes, dimUnlim )
        filevardef(cf_file, var_x_bnds, var_x_bnds_type, (/"lon","bnds"/) )
        filevardef(cf_file, var_y_bnds, var_y_bnds_type, (/"lat","bnds"/) )
        filevardef(cf_file, var_t_bnds, var_t_bnds_type, (/"time","bnds"/) )
        filevardef(cf_file, var_t, var_t_type, (/"time"/))
;--------------------------------------------------------------------
; Write coordinates and their bounds to output file
;--------------------------------------------------------------------
        lat@standard_name    = "latitude"
        lat@axis             = "Y"
        lon@standard_name    = "longitude"
        lon@axis             = "X"
;
        raw!0        = "time"
        raw&time     = tnew
        tnew@bounds        = "time_bnds"
        tnew@units         = "days since 0000-1-1"
        tnew@axis          = "T"
        tnew@standard_name = "time"
        tnew@calendar      = "noleap"
        raw!1        = "lat"
        raw&lat      = lat
        raw!2        = "lon"
        raw&lon      = lon

        cf_file->lat = lat
        cf_file->lat_bnds = lat_bnds
        cf_file->lon = lon
        cf_file->lon_bnds = lon_bnds
        cf_file->$ipcc_vars(v)$ = raw
      end if

      cf_file->time      = tnew
      cf_file->time_bnds = time_bnds
;
      delete(raw)
      delete(dimNames)
      delete(dimSizes)
      delete(dimUnlim)
    end if
  end do
end
