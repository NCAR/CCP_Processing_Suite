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
;--------------------------------------------------------------------
; Read input file
;--------------------------------------------------------------------
  in_file = addfile (ncin,"r")
;--------------------------------------------------------------------
; Read input file coordinates
;--------------------------------------------------------------------
  time = in_file->time
  date = in_file->mcdate
  lat  = in_file->lat
  lon  = in_file->lon
  nlon = dimsizes(lon)
  nlat = dimsizes(lat)
  ntim = dimsizes(time)
  nbnd = 2
;--------------------------------------------------------------------
; Calculate lon/lat bounds
; Note that bound(X,1) is the same as bound(X-1,2)
;-----------------------------------------------------------------
  lat_bnds  = new((/nlat,nbnd/),double)
  lat_bnds(0,0)      = -90.
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
  do n = 0,ntim-1
    iyr = date(n)/10000
    imn = (date(n)-(iyr*10000))/100
    if (imn.eq.1) then
      iyr = iyr - 1
      imn = 12
    else
      imn = imn - 1
    end if
    time_bnds(n,0) = ((iyr*365)+mbeg(imn-1))-1
    time_bnds(n,1) = ((iyr*365)+mend(imn-1))
    tnew(n) = (time_bnds(n,0)+time_bnds(n,1))/2
  end do
  delete(tnew@_FillValue)
  delete(time_bnds@_FillValue)
  delete(lat_bnds@_FillValue)
  delete(lon_bnds@_FillValue)
  lat@bounds  = "lat_bnds"
  lon@bounds  = "lon_bnds"
  tnew@bounds = "time_bnds"
;--------------------------------------------------------------------
; Depth and depth bounds
;--------------------------------------------------------------------
  depth = new(1,double)
  depth = 0.05d0
  depth!0 = "ncl_scalar"
  depth@standard_name = "depth"
  depth@long_name     = "depth"
  depth@units         = "m"
  depth@axis          = "Z"
  depth@positive      = "up"

  depth_bnds = new((/nbnd/),double)
  depth_bnds = (/0.d0,0.1d0/)
  depth_bnds@standard_name = "depth_bnds"
  depth_bnds@long_name     = "depth_bnds"
  depth_bnds@units         = "m"
  depth_bnds@axis          = "Z"
  depth_bnds@positive      = "up"
  depth@bounds             = "depth_bnds"

  delete(depth@_FillValue)
  delete(depth_bnds@_FillValue)
;--------------------------------------------------------------------
; Read input field(s), using IPCC required fields as parser
;--------------------------------------------------------------------
  table           = "A1"
  glbAtt@table_id = "Table "+table

  ipcc_vars = (/"mrso","mrsos","mrfso","mrros","mrro","snw","snc","snm","snd"/);  

  nvars = dimsizes(ipcc_vars)
  do v = 0,nvars - 1
    if (ipcc_vars(v).eq."mrsos") then
      if ((isfilevar(in_file,"SOILLIQ")).and.(isfilevar(in_file,"SOILICE"))) then
        tmp                = in_file->SOILLIQ(:,0:2,:,:)
        tmp                = in_file->SOILLIQ(:,0:2,:,:)+in_file->SOILICE(:,0:2,:,:)
        raw                = dim_sum_Wrap(tmp(time|:,lat|:,lon|:,levsoi|:))
        delete(tmp)
        raw@standard_name  = "moisture_content_of_soil_layer"
        raw@original_name  = "SOILLIQ,SOILICE"
        raw@original_units = "kg m-2"
        raw@units          = "kg m-2"
        raw@history        = "SOILLIQ+SOILICE, 0-0.1m, added depth coordinate"
        raw@coordinates    = "depth"
      end if
    end if
    if (ipcc_vars(v).eq."mrso") then
      if ((isfilevar(in_file,"SOILLIQ")).and.(isfilevar(in_file,"SOILICE"))) then
        tmp                = in_file->SOILLIQ(:,:,:,:)
        tmp                = in_file->SOILLIQ(:,:,:,:)+in_file->SOILICE(:,:,:,:)
        raw                = dim_sum_Wrap(tmp(time|:,lat|:,lon|:,levsoi|:))
        delete(tmp)
        raw@standard_name  = "soil_moisture_content"
        raw@original_name  = "SOILLIQ,SOILICE"
        raw@original_units = "kg m-2"
        raw@units          = "kg m-2"
        raw@history        = "SOILLIQ+SOILICE, added depth coordinate"
      end if
    end if
    if (ipcc_vars(v).eq."mrfso") then
      if (isfilevar(in_file,"SOILICE")) then
        tmp                = in_file->SOILICE(:,:,:,:)
        raw                = dim_sum_Wrap(tmp(time|:,lat|:,lon|:,levsoi|:))
        delete(tmp)
        raw@standard_name  = "soil_frozen_water_content"
        raw@original_name  = "SOILICE"
        raw@original_units = "kg m-2"
        raw@units          = "kg m-2"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."mrros") then
      if (isfilevar(in_file,"QOVER")) then
        raw                = in_file->QOVER
        raw@standard_name  = "surface_runoff_flux"
        raw@original_name  = "QOVER"
        raw@original_units = "mm s-1"
        raw@units          = "kg m-2 s-1"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."snd") then
      if (isfilevar(in_file,"SNOWDP")) then
        raw                = in_file->SNOWDP
        raw@standard_name  = "surface_snow_thickness"
        raw@original_name  = "SNOWDP"
        raw@original_units = raw@units
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."mrro") then
      if ((isfilevar(in_file,"QOVER")).and.\
          (isfilevar(in_file,"QRGWL")).and.\
          (isfilevar(in_file,"QDRAI"))) then
        raw                = in_file->QOVER
        raw                = in_file->QOVER+in_file->QRGWL+in_file->QDRAI
        raw@standard_name  = "runoff_flux"
        raw@original_name  = "QOVER+QRGWL+QDRAI"
        raw@original_units = "mm s-1"
        raw@units          = "kg m-2 s-1"
        raw@history        = "QOVER+QRGWL+QDRAI"
      end if
    end if
    if (ipcc_vars(v).eq."snw") then
      if ((isfilevar(in_file,"SNOWLIQ")).and.(isfilevar(in_file,"SNOWICE"))) then
        raw                = in_file->SNOWLIQ
        raw                = in_file->SNOWLIQ+in_file->SNOWICE
        raw@standard_name  = "surface_snow_amount_where_land"
        raw@original_name  = "SNOWICE,SNOWLIQ"
        raw@original_units = "kg m-2"
        raw@units          = "kg m-2"
        raw@history        = "SNOWICE+SNOWLIQ"
      end if
    end if
    if (ipcc_vars(v).eq."snc") then
      if (isfilevar(in_file,"FSNO")) then
        raw                = in_file->FSNO
        raw                = in_file->FSNO*100.
        raw@standard_name  = "surface_snow_area_fraction_where_land"
        raw@original_name  = "FSNO"
        raw@original_units = "unitless"
        raw@units          = "percent"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."snm") then
      if (isfilevar(in_file,"QMELT")) then
        raw                = in_file->QMELT
        raw@standard_name  = "surface_snow_melt_flux_where_land"
        raw@original_name  = "QMELT"
        raw@original_units = "mm s-1"
        raw@units          = "kg m-2 s-1"
        raw@history        = "No change"
      end if
    end if

    if (isvar("raw")) then

;--------------------------------------------------------------------
; Field-invariant metadata
;--------------------------------------------------------------------
      raw@long_name      = raw@standard_name
      raw@cell_methods   = "time: mean (interval: 1 month)"
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
      now = systemfunc("date '+%y%m%d'")
      nccf = "./IPCC/"+ipcc_vars(v)+"_"+table+"."+expid+".CCSM."+type+"."+tval+"_c"+now+".nc"
      system("csh -c 'rm -f "+nccf+"'")
      cf_file = addfile(nccf, "c")
      fileattdef(cf_file, glbAtt)
;--------------------------------------------------------------------
; Setups for output file
;--------------------------------------------------------------------
; Coordinates
;--------------------------------------------------------------------
      var_x_bnds      = (/ "lon_bnds" /)
      var_x_bnds_type = (/ "double" /)
      var_y_bnds      = (/ "lat_bnds" /)
      var_y_bnds_type = (/ "double" /)
      var_t_bnds      = (/ "time_bnds" /)
      var_t_bnds_type = (/ "double" /)
      var_x           = (/ "lon" /)
      var_x_type      = (/ "double" /)
      var_y           = (/ "lat" /)
      var_y_type      = (/ "double" /)
      var_t           = (/ "time" /)
      var_t_type      = (/ "double" /)
;
      dimNames = (/ "lon", "lat", "bnds", "time" /)
      dimSizes = (/  nlon,  nlat,   nbnd, -1     /)
      dimUnlim = (/ False, False,  False, True    /)
      filedimdef(cf_file, dimNames, dimSizes, dimUnlim )
      filevardef(cf_file, var_x_bnds, var_x_bnds_type, (/"lon","bnds"/) )
      filevardef(cf_file, var_y_bnds, var_y_bnds_type, (/"lat","bnds"/) )
      filevardef(cf_file, var_t_bnds, var_t_bnds_type, (/"time","bnds"/) )
      filevardef(cf_file, var_x, var_x_type, (/"lon" /))
      filevardef(cf_file, var_y, var_y_type, (/"lat" /))
      filevardef(cf_file, var_t, var_t_type, (/"time"/))
;--------------------------------------------------------------------
; Write coordinates and their bounds to output file
;--------------------------------------------------------------------
      lat@standard_name  = "latitude"
      lat@axis           = "Y"
      lon@standard_name  = "longitude"
      lon@axis           = "X"
      tnew@units         = "days since 0000-1-1"
      tnew@axis          = "T"
      tnew@standard_name = "time"
      tnew@calendar      = "noleap"
;
      delete(raw&time)
      raw&time           = tnew
;
      cf_file->time = tnew
      cf_file->time_bnds = time_bnds
      cf_file->lat = lat
      cf_file->lat_bnds = lat_bnds
      cf_file->lon = lon
      cf_file->lon_bnds = lon_bnds
;
      if (ipcc_vars(v).eq."mrsos") then
        depth_bnds!0             = "bnds"
        cf_file->depth      = depth
        cf_file->depth_bnds = depth_bnds
      end if
;--------------------------------------------------------------------
; The big deal...
;--------------------------------------------------------------------
      cf_file->$ipcc_vars(v)$ = raw
      delete(raw)
      delete(dimNames)
      delete(dimSizes)
      delete(dimUnlim)
    end if
  end do
end
