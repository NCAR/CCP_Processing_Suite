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
; Latent heat of evaporation, J kg-1
  latent_ht_v = 2.501e6
; Latent heat of fusion, J kg-1
  latent_ht_f = 3.337e5
; grams to kilograms
  g_to_kg     = 1000.
;
  n1 = integertochar(10)
;--------------------------------------------------------------------
; Get input file env vars
;--------------------------------------------------------------------
  ncin = getenv("CCSM_FILE")
  case = getenv("CCSM_CASE")
  comp = getenv("CCSM_COMP")
  tval = getenv("CCSM_TIME")
  scpt = getenv("NCL_SCPT")
;--------------------------------------------------------------------
; Define global attributes
;--------------------------------------------------------------------
  glbAtt = True
;--------------------------------------------------------------------
; Set various attributes
;--------------------------------------------------------------------
  parent  = ""
  year1   = ""
  machine = ""
  forcing = ""
  expid   = ""
  attributes_ccsm(case,glbAtt,parent,year1,machine,forcing,expid)
;--------------------------------------------------------------------
; Define table and ipcc_vars
;--------------------------------------------------------------------
  table     = "O1"
  ipcc_vars = (/"sic","sit","usi","vsi"/)
;--------------------------------------------------------------------
; Read input file
;--------------------------------------------------------------------
  in_file = addfile (ncin,"r")
  nbnd = 2
;--------------------------------------------------------------------
; Read input file T coordinates
;--------------------------------------------------------------------
  time = (/in_file->time/)
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
  iyr       = (floor(fyr))+1
  rmn       = fyr - (iyr - 1)
  imn       = floor(rmn*365)-15
  jmn       = floattoint(floor((imn/(365/12))))
  time_bnds(0:ntim-1,0) = ((iyr*365)+mbeg(jmn)-1)
  time_bnds(0:ntim-1,1) = ( iyr*365)+mend(jmn)
  tnew = (time_bnds(:,0)+time_bnds(:,1))/2
  do n = 0,ntim-1
    print(time_bnds(n,0)+" "+tnew(n)+" "+time_bnds(n,1))
  end do
  delete(tnew@_FillValue)
  delete(time_bnds@_FillValue)
;--------------------------------------------------------------------
; Read input field(s), using IPCC required fields as parser
;--------------------------------------------------------------------
  nvars = dimsizes(ipcc_vars)
  do v = 0,nvars - 1
    if (ipcc_vars(v).eq."sic") then
      if (isfilevar(in_file,"aice")) then
        raw                = (/in_file->aice/)
        raw                = raw < 1.e20
        raw@standard_name  = "sea_ice_area_fraction"
        raw@original_name  = "aice"
        raw@original_units = "%"
        raw@units          = "%"
        raw@history        = "Interpolated to regular grid from dipole grid"
        raw@cell_methods   = "time: mean (interval: 1 month)"
        glbAtt@table_id    = "Table "+table
        nccf               = "./IPCC/"+ipcc_vars(v)+"_"+table+"."+expid+".CCSM."+comp+"."+tval+".nc"
      end if
    end if
    if (ipcc_vars(v).eq."sit") then
      if (isfilevar(in_file,"hi")) then
        raw                = (/in_file->hi/)
        raw                = raw < 1.e20
        raw@standard_name  = "sea_ice_thickness"
        raw@original_name  = "hi"
        raw@original_units = "m"
        raw@units          = "m"
        raw@history        = "Interpolated to regular grid from dipole grid"
        raw@cell_methods   = "time: mean (interval: 1 month)"
        glbAtt@table_id    = "Table "+table
        nccf               = "./IPCC/"+ipcc_vars(v)+"_"+table+"."+expid+".CCSM."+comp+"."+tval+".nc"
      end if
    end if
    if (ipcc_vars(v).eq."usi") then
      if (isfilevar(in_file,"u")) then
        raw                = (/in_file->u/)/100
        raw                = raw < 1.e20
        raw@standard_name  = "eastward_sea_ice_velocity"
        raw@original_name  = "u"
        raw@original_units = "cm s-1"
        raw@units          = "m s-1"
        raw@history        = "Interpolated and rotated to regular grid from dipole grid"+n1+"u/100"
        raw@cell_methods   = "time: mean (interval: 1 month)"
        glbAtt@table_id    = "Table "+table
        nccf               = "./IPCC/"+ipcc_vars(v)+"_"+table+"."+expid+".CCSM."+comp+"."+tval+".nc"
      end if
    end if
    if (ipcc_vars(v).eq."vsi") then
      if (isfilevar(in_file,"v")) then
        raw                = (/in_file->v/)/100
        raw                = raw < 1.e20
        raw@standard_name  = "northward_sea_ice_velocity"
        raw@original_name  = "v"
        raw@original_units = "cm s-1"
        raw@units          = "m s-1"
        raw@history        = "Interpolated and rotated to regular grid from dipole grid"+n1+"v/100"
        raw@cell_methods   = "time: mean (interval: 1 month)"
        glbAtt@table_id    = "Table "+table
        nccf               = "./IPCC/"+ipcc_vars(v)+"_"+table+"."+expid+".CCSM."+comp+"."+tval+".nc"
      end if
    end if
;--------------------------------------------------------------------
; Field-invariant metadata
;--------------------------------------------------------------------
    if (isvar("raw")) then
     raw@long_name      = raw@standard_name
     raw@_FillValue     = 1.e+20
     raw@missing_value  = 1.e+20
     raw@comment        = "Created using NCL code "+scpt+" on"+n1+\
                         " machine "+systemfunc("hostname")
;--------------------------------------------------------------------
; Create output file
; [file name as per requirements document]
;--------------------------------------------------------------------
     system("csh -c 'rm -f "+nccf+"'")
     cf_file = addfile(nccf, "c")
     fileattdef(cf_file, glbAtt)
;--------------------------------------------------------------------
; Setups for output file
;--------------------------------------------------------------------
; Coordinates
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
     delete(lat@_FillValue)
     delete(lon@_FillValue)
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
     raw!0                = "time"
     raw&time             = tnew
     tnew@bounds          = "time_bnds"
     tnew@units           = "days since 0000-1-1"
     tnew@axis            = "T"
     tnew@standard_name   = "time"
     tnew@calendar        = "noleap"
     raw!1                = "lat"
     raw&lat              = lat
     raw!2                = "lon"
     raw&lon              = lon

     cf_file->lat            = lat
     cf_file->lat_bnds       = lat_bnds
     cf_file->lon            = lon
     cf_file->lon_bnds       = lon_bnds
     cf_file->time           = tnew
     cf_file->time_bnds      = time_bnds
     cf_file->$ipcc_vars(v)$ = raw
;
     delete(raw)
     delete(dimNames)
     delete(dimSizes)
     delete(dimUnlim)
    end if
  end do
end
