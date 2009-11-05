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
;
  table   = "A2"
  glbAtt@table_id = "Table "+table
;--------------------------------------------------------------------
; Read input file coordinates
;--------------------------------------------------------------------
  time = in_file->time
  date = in_file->date
  lat  = in_file->lat
  lev  = (/in_file->lev/)
  ilev = (/in_file->ilev/)
  lon  = in_file->lon
  nlon = dimsizes(lon)
  nlat = dimsizes(lat)
  nlev = dimsizes(lev)
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
  mbeg = (/ 0, 31, 59, 90,120,151,181,212,243,273,304,334/)
  time_bnds = new((/ntim,nbnd/),double)
  tnew      = new((/ntim/),double)
  do n = 0,ntim-1
    iyr = date(n)/10000
    imn = (date(n)-(iyr*10000))/100
    idy = date(n)-(iyr*10000)-(imn*100)
    time_bnds(n,0) = ((iyr*365)+mbeg(imn-1))+idy-1
    time_bnds(n,1) = ((iyr*365)+mbeg(imn-1))+idy
    tnew(n) = (time_bnds(n,0)+time_bnds(n,1))/2
  end do
  del_FillValue(tnew)
  del_FillValue(time_bnds)
  del_FillValue(lat_bnds)
  del_FillValue(lon_bnds)
  lat@bounds  = "lat_bnds"
  lon@bounds  = "lon_bnds"
  tnew@bounds = "time_bnds"
;--------------------------------------------------------------------
; Fields needed for vertical interpolation to IPCC pressure levels
;--------------------------------------------------------------------
  if (isfilevar(in_file, "PS")) then
    hyam = in_file->hyam
    hybm = in_file->hybm
    hyai = in_file->hyai
    hybi = in_file->hybi
    ps   = in_file->PS
    P0mb = new(1,float)
    del_FillValue(P0mb)
    P0mb = 0.01*(/dble2flt(in_file->P0)/)
    phis = in_file->PHIS
    if (isfilevar(in_file, "TEMP")) then
      tbot = in_file->TEMP(:,nlev-1,:,:)
    end if
    if (isfilevar(in_file, "T")) then
      tbot = in_file->T(:,nlev-1,:,:)
    end if
  end if
;--------------------------------------------------------------------
; IPCC-required pressure levels
;--------------------------------------------------------------------
  plev               = (/1000.,925.,850.,700.,600.,500.,400.,300.,250.,200.,\
                          150.,100., 70., 50., 30.,20., 10./)
  plev@long_name     = "pressure"
; units are hPa at this point but will be written as Pa
  plev@units         = "Pa"
  plev@standard_name = "air_pressure"
  plev@positive      = "down"
  plev@axis          = "Z"
  plev!0             = "plev"
  plev&plev          = plev
  nplev              = dimsizes(plev)
;--------------------------------------------------------------------
; Read input field(s), using IPCC required fields as parser
  height = 2.d0
  height!0 = "ncl_scalar"
  height@standard_name = "height"
  height@long_name = "height"
  height@units = "m"
  height@axis = "Z"
  height@positive = "up"
  del_FillValue(height)

  ipcc_vars = (/"psl","pr","tas","tasmin","tasmax",\
                "hfls","hfss","rlds","rlus","rsds","rsus","rlut"/)
  nvars = dimsizes(ipcc_vars)
  do v = 0,nvars - 1
    if (ipcc_vars(v).eq."psl") then
      if (isfilevar(in_file,"PSL")) then
        raw                = in_file->PSL
        raw@standard_name  = "air_pressure_at_sea_level"
        raw@original_name  = "PSL"
        raw@original_units = raw@units
        raw@cell_methods   = "time: mean (interval: 1 day)"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."pr") then
      if ((isfilevar(in_file,"PRECC")).and.(isfilevar(in_file,"PRECL"))) then
        raw                =  in_file->PRECC
        raw                = (in_file->PRECC+in_file->PRECL)*density_h2o
        raw@standard_name  = "precipitation_flux"
        raw@original_name  = "PRECC,PRECL"
        raw@original_units = "m-1 s-1"
        raw@units          = "kg m-2 s-1"
        raw@cell_methods   = "time: mean (interval: 1 day)"
        raw@history        = "(PRECC+PRECL)*density_h2o"
      end if
    end if
    if (ipcc_vars(v).eq."tas") then
      if (isfilevar(in_file,"TREFHT")) then
        raw                = in_file->TREFHT
        raw@standard_name  = "air_temperature"
        raw@original_name  = "TREFHT"
        raw@original_units = raw@units
        raw@coordinates    = "height"
        raw@cell_methods   = "time: mean (interval: 1 day)"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."tasmin") then
      if (isfilevar(in_file,"TREFHTMN")) then
        raw                = in_file->TREFHTMN
        raw@standard_name  = "air_temperature"
        raw@original_name  = "TREFHTMN"
        raw@original_units = raw@units
        raw@coordinates    = "height"
        raw@cell_methods   = "time: minimum within days"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."tasmax") then
      if (isfilevar(in_file,"TREFHTMX")) then
        raw                = in_file->TREFHTMX
        raw@standard_name  = "air_temperature"
        raw@original_name  = "TREFHTMX"
        raw@original_units = raw@units
        raw@coordinates    = "height"
        raw@cell_methods   = "time: maximum within days"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."hfls") then
      if (isfilevar(in_file,"LHFLX")) then
        raw                = in_file->LHFLX
        raw@standard_name  = "surface_upward_latent_heat_flux"
        raw@original_name  = "LHFLX"
        raw@original_units = "W m-2"
        raw@units          = "W m-2"
        raw@cell_methods   = "time: mean (interval: 1 day)"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."hfss") then
      if (isfilevar(in_file,"SHFLX")) then
        raw                = in_file->SHFLX
        raw@standard_name  = "surface_upward_sensible_heat_flux"
        raw@original_name  = "SHFLX"
        raw@original_units = "W m-2"
        raw@units          = "W m-2"
        raw@cell_methods   = "time: mean (interval: 1 day)"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."rlds") then
      if ((isfilevar(in_file,"TS")).and.(isfilevar(in_file,"FLNS"))) then
        raw                = in_file->TS
        raw                = sb*((in_file->TS)^4)-(in_file->FLNS)
        raw@standard_name  = "surface_downwelling_longwave_flux_in_air"
        raw@original_name  = "TS,FLNS"
        raw@original_units = "K,W m-2"
        raw@units          = "W m-2"
        raw@history        = "SB*TS^4-FLNS"
        raw@cell_methods   = "time: mean (interval: 1 day)"
      end if
    end if
    if (ipcc_vars(v).eq."rlus") then
      if (isfilevar(in_file,"TS")) then
        raw                = in_file->TS
        raw                = sb*((in_file->TS)^4)
        raw@standard_name  = "surface_upwelling_longwave_flux_in_air"
        raw@original_name  = "TS"
        raw@original_units = "K"
        raw@units          = "W m-2"
        raw@history        = "SB*TS^4"
        raw@cell_methods   = "time: mean (interval: 1 day)"
      end if
    end if
    if (ipcc_vars(v).eq."rsds") then
      if (isfilevar(in_file,"FSDS")) then
        raw                = in_file->FSDS
        raw@standard_name  = "surface_downwelling_shortwave_flux_in_air"
        raw@original_name  = "FSDS"
        raw@original_units = raw@units
        raw@units          = "W m-2"
        raw@cell_methods   = "time: mean (interval: 1 day)"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."rsus") then
      if ((isfilevar(in_file,"FSDS")).and.(isfilevar(in_file,"FSNS"))) then
        raw                = in_file->FSDS
        raw                = in_file->FSDS-in_file->FSNS
        raw@standard_name  = "surface_upwelling_shortwave_flux_in_air"
        raw@original_name  = "FSDS,FSNS"
        raw@original_units = "W m-2"
        raw@units          = "W m-2"
        raw@cell_methods   = "time: mean (interval: 1 day)"
        raw@original_name  = "FSDS-FSNS"
      end if
    end if
    if (ipcc_vars(v).eq."rlut") then
      if (isfilevar(in_file,"FLUT")) then
        raw                =  in_file->FLUT
        raw@standard_name  = "toa_outgoing_longwave_flux"
        raw@original_name  = "FLUT"
        raw@original_units = raw@units
        raw@units          = "W m-2"
        raw@cell_methods   = "time: mean (interval: 1 day)"
        raw@history        = "No change"
      end if
    end if

    if (isvar("raw")) then

;--------------------------------------------------------------------
; Field-invariant metadata
;--------------------------------------------------------------------
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
      nccf = "./CMIP/"+ipcc_vars(v)+"_"+table+"."+expid+".CCSM."+type+"."+tval+".nc"
      system("csh -c 'rm -f "+nccf+"'")
      cf_file = addfile(nccf, "c")
      fileattdef(cf_file, glbAtt)
;--------------------------------------------------------------------
; Setups for output file
;--------------------------------------------------------------------
; Coordinates
;--------------------------------------------------------------------
      if ( num(dimsizes(raw)) .eq. 4) then
        var_x_bnds      = (/ "lon_bnds" /)
        var_x_bnds_type = (/ "double" /)
        var_y_bnds      = (/ "lat_bnds" /)
        var_y_bnds_type = (/ "double" /)
        var_t_bnds      = (/ "time_bnds" /)
        var_t_bnds_type = (/ "double" /)
        var_t           = (/ "time" /)
        var_t_type      = (/ "double" /)
        raw@comment        = "Interpolated to standard pressure levels, "+n1+\
                             " added PCMDI-required metadata and renamed "+n1+\
                             " field to CF 1.0 compliance"
        var_z           = (/ "plev" /)
        var_z_type      = (/ "double" /)
        dimNames = (/ "lon", "lat", "plev", "bnds", "time" /)
        dimSizes = (/  nlon,  nlat, nplev ,  nbnd, -1      /)
        dimUnlim = (/ False, False,  False, False, True    /)
        filedimdef(cf_file, dimNames, dimSizes, dimUnlim )
        filevardef(cf_file, var_x_bnds, var_x_bnds_type, (/"lon","bnds"/) )
        filevardef(cf_file, var_y_bnds, var_y_bnds_type, (/"lat","bnds"/) )
        filevardef(cf_file, var_t_bnds, var_t_bnds_type, (/"time","bnds"/) )
        filevardef(cf_file, var_t, var_t_type, (/"time"/))
        filevardef(cf_file, var_z, var_z_type, (/"plev"/))
      else
        var_x_bnds      = (/ "lon_bnds" /)
        var_x_bnds_type = (/ "double" /)
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
      end if
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
      raw&time     = tnew
;
      cf_file->time = tnew
      cf_file->time_bnds = time_bnds
      cf_file->lat = lat
      cf_file->lat_bnds = lat_bnds
      cf_file->lon = lon
      cf_file->lon_bnds = lon_bnds
;
      if ((ipcc_vars(v).eq."tas").or.\
          (ipcc_vars(v).eq."tasmin").or.\
          (ipcc_vars(v).eq."tasmax")) then
        cf_file->height = height
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
;********************************************************************************
; del_time
; get rid of a bogus time attribute from x
;
undef("del_time")
procedure del_time (x)
; usage:     del_time(T)
;            del_time(ps)
begin
    if (isdim(x,"time")) then
        delete(x&time)
    end if
end   
