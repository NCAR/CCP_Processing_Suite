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
;
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
;
  table   = "A1"
  glbAtt@table_id = "Table "+table
;--------------------------------------------------------------------
; Read input file
;--------------------------------------------------------------------
  in_file = addfile (ncin,"r")
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
  hyam = in_file->hyam
  hybm = in_file->hybm
  hyai = in_file->hyai
  hybi = in_file->hybi
  ps   = in_file->PS
;  P0mb = dble2flt(0.01*in_file->P0)
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
;--------------------------------------------------------------------
; IPCC-required pressure levels
;--------------------------------------------------------------------
  plev               = (/1000.,925.,850.,700.,600.,500.,400.,300.,250.,200.,\
                          150.,100., 70., 50., 30.,20., 10./)
  plev@long_name     = "pressure"
; units are hPa at this point but will be written as Pa
  plev@units	     = "Pa"
  plev@standard_name = "air_pressure"
  plev@positive      = "down"
  plev@axis          = "Z"
  plev!0             = "plev"
  plev&plev          = plev
  nplev              = dimsizes(plev)
;
  height = 2.d0
  height!0 = "ncl_scalar"
  height@standard_name = "height"
  height@long_name = "height"
  height@units = "m"
  height@axis = "Z"
  height@positive = "up"
  del_FillValue(height)
;--------------------------------------------------------------------
; Read input field(s), using IPCC required fields as parser
;--------------------------------------------------------------------
  ipcc_vars = (/"psl","pr","tas","tasmin","tasmax","tauu","tauv",\
                "rlds","rlus","rsds","rsus","ts","ps","prsn","prc",\
                "prw","huss","rsdt","snd","hfls","hfss",\
                "rsut","rlut","rtmt","rsdscs","rsuscs","rldscs","rlutcs",\
                "rsutcs","clt","clwvi","clivi",\
                "cl",\
                "ta","ua","va","hus","wap","zg","hur","tro3",\
                "taE","uaE","vaE","husE","wapE","zgE","hurE","tro3E"/)
  nvars = dimsizes(ipcc_vars)
  do v = 0,nvars - 1
    if (ipcc_vars(v).eq."psl") then
      if (isfilevar(in_file,"PSL")) then
        raw                = in_file->PSL
        raw@standard_name  = "air_pressure_at_sea_level"
        raw@original_name  = "PSL"
        raw@original_units = raw@units
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."pr") then
      if ((isfilevar(in_file,"PRECC")).and.(isfilevar(in_file,"PRECL"))) then
        raw                =  in_file->PRECC
        raw                = (in_file->PRECC+in_file->PRECL)*density_h2o
        raw@standard_name  = "precipitation_flux"
        raw@original_name  = "PRECC, PRECL"
        raw@original_units = "m-1 s-1"
        raw@units          = "kg m-2 s-1"
        raw@history        = "(PRECC+PRECL)*r[h2o]"
      end if
    end if
    if (ipcc_vars(v).eq."tas") then
      if (isfilevar(in_file,"TREFHT")) then
        raw                = in_file->TREFHT
        raw@standard_name  = "air_temperature"
        raw@original_name  = "TREFHT"
        raw@original_units = raw@units
        raw@coordinates	 = "height"
        raw@history        = "Added height coordinate"
      end if
    end if
    if (ipcc_vars(v).eq."tasmin") then
      if (isfilevar(in_file,"TREFMNAV")) then
        raw                = in_file->TREFMNAV
        raw@standard_name  = "air_temperature"
        raw@original_name  = "TREFMNAV"
        raw@original_units = raw@units
        raw@coordinates	 = "height"
        raw@history        = "Added height coordinate"
        raw@cell_methods   = "time: minimum within days time: mean over days"
      end if
    end if
    if (ipcc_vars(v).eq."tasmax") then
      if (isfilevar(in_file,"TREFMXAV")) then
        raw                = in_file->TREFMXAV
        raw@standard_name  = "air_temperature"
        raw@original_name  = "TREFMXAV"
        raw@original_units = raw@units
        raw@coordinates	 = "height"
        raw@history        = "Added height coordinate"
        raw@cell_methods   = "time: maximum within days time: mean over days"
      end if
    end if
    if (ipcc_vars(v).eq."tauu") then
      if (isfilevar(in_file,"TAUX")) then
        raw                = in_file->TAUX
        raw@standard_name  = "surface_downward_eastward_stress"
        raw@original_name  = "TAUX"
        raw@original_units = raw@units
        raw@units          = "N m-2"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."tauv") then
      if (isfilevar(in_file,"TAUY")) then
        raw                = in_file->TAUY
        raw@standard_name  = "surface_downward_northward_stress"
        raw@original_name  = "TAUY"
        raw@original_units = raw@units
        raw@units          = "N m-2"
        raw@history        = "No change"
      end if
    end if
;    if (ipcc_vars(v).eq."snd") then
;      if (isfilevar(in_file,"SNOWHLND")) then
;        raw                = in_file->SNOWHLND
;        raw@standard_name  = "surface_snow_thickness"
;        raw@original_name  = "SNOWHLND"
;        raw@original_units = raw@units
;        raw@history        = "No change"
;      end if
;    end if
; --------------------------------------------------------------------------------
; Correct
; --------------------------------------------------------------------------------
    if (ipcc_vars(v).eq."hfls") then
      if (isfilevar(in_file,"LHFLX")) then
        raw                = in_file->LHFLX
        raw@standard_name  = "surface_upward_latent_heat_flux"
        raw@original_name  = "LHFLX"
        raw@original_units = "W m-2"
        raw@units          = "W m-2"
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
      end if
    end if
    if (ipcc_vars(v).eq."rsds") then
      if (isfilevar(in_file,"FSDS")) then
        raw                = in_file->FSDS
        raw@standard_name  = "surface_downwelling_shortwave_flux_in_air"
        raw@original_name  = "FSDS"
        raw@original_units = raw@units
        raw@units          = "W m-2"
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
        raw@history        = "FSDS-FSNS"
      end if
    end if
    if (ipcc_vars(v).eq."ts") then
      if (isfilevar(in_file,"TS")) then
        raw                = in_file->TS
        raw@standard_name  = "surface_temperature"
        raw@original_name  = "TS"
        raw@original_units = raw@units
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."ps") then
      if (isfilevar(in_file,"PS")) then
        raw                = in_file->PS
        raw@standard_name  = "surface_air_pressure"
        raw@original_name  = "PS"
        raw@original_units = raw@units
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."prsn") then
      if ((isfilevar(in_file,"PRECSC")).and.(isfilevar(in_file,"PRECSL"))) then
        raw                =  in_file->PRECSC
        raw                = (in_file->PRECSC+in_file->PRECSL)*density_h2o
        raw@standard_name  = "snowfall_flux"
        raw@original_name  = "PRECSC,PRECSL"
        raw@original_units = "m-1 s-1"
        raw@units          = "kg m-2 s-1"
        raw@history        = "(PRECSC+PRECSL)*r[h2o]"
      end if
    end if
    if (ipcc_vars(v).eq."prc") then
      if (isfilevar(in_file,"PRECC")) then
        raw                =  in_file->PRECC
        raw                =  raw*density_h2o
        raw@standard_name  = "convective_precipitation_flux"
        raw@original_name  = "PRECC"
        raw@original_units = "m-1 s-1"
        raw@units          = "kg m-2 s-1"
        raw@history        = "PRECC*r[h2o]"
      end if
    end if
    if (ipcc_vars(v).eq."prw") then
      if (isfilevar(in_file,"TMQ")) then
        raw                = in_file->TMQ
        raw@standard_name  = "atmosphere_water_vapor_content"
        raw@original_name  = "TMQ"
        raw@original_units = raw@units
        raw@units          = "kg m-2"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."huss") then
      if (isfilevar(in_file,"QREFHT")) then
        raw                =  in_file->QREFHT
        raw@standard_name  = "specific_humidity"
        raw@original_name  = "QREFHT"
        raw@original_units = "kg/kg"
        raw@units          = "kg kg-1"
        raw@coordinates	 = "height"
        raw@history	 = "Added height coordinate"
        raw@note   	 = "Use 'huss' with caution. There are incorrect values over"+n1+\
                           " Antarctica in January. This data is from the CCSM3 model"+n1+\
                           " coupler code, not the atmospheric model code."
      end if
    end if
    if (ipcc_vars(v).eq."rsdt") then
      if (isfilevar(in_file,"SOLIN")) then
        raw                =  in_file->SOLIN
        raw@standard_name  = "toa_incoming_shortwave_flux"
        raw@original_name  = "SOLIN"
        raw@original_units = raw@units
        raw@units          = "W m-2"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."rsut") then
      if ((isfilevar(in_file,"SOLIN")).and.(isfilevar(in_file,"FSNTOA"))) then
        raw                =  in_file->SOLIN
        raw                =  in_file->SOLIN-in_file->FSNTOA
        raw@standard_name  = "toa_outgoing_shortwave_flux"
        raw@original_name  = "SOLIN,FSNTOA"
        raw@original_units = "W m-2"
        raw@units          = "W m-2"
        raw@history        = "SOLIN-FSNTOA"
      end if
    end if
    if (ipcc_vars(v).eq."rlut") then
      if (isfilevar(in_file,"FLUT")) then
        raw                =  in_file->FLUT
        raw@standard_name  = "toa_outgoing_longwave_flux"
        raw@original_name  = "FLUT"
        raw@original_units = raw@units
        raw@units          = "W m-2"
        raw@history        = "No change"
      end if
    end if
; --------------------------------------------------------------------------------
; Correct
; --------------------------------------------------------------------------------
    if (ipcc_vars(v).eq."rtmt") then
      if ((isfilevar(in_file,"FSNT")).and.(isfilevar(in_file,"FLNT"))) then
        raw                = in_file->FSNT
        raw                = in_file->FSNT-in_file->FLNT
        raw@standard_name  = "net_downward_radiative_flux_at_top_of_atmosphere_model"
        raw@original_name  = "FSNT,FLNT"
        raw@original_units = "W m-2"
        raw@units          = "W m-2"
        raw@history        = "FSNT-FLNT"
      end if
    end if
    if (ipcc_vars(v).eq."rsdscs") then
      if (isfilevar(in_file,"FSDSC")) then
        raw                = in_file->FSDSC
        raw@standard_name  = "surface_downwelling_shortwave_flux_in_air_assuming_clear_sky"
        raw@original_name  = "FSDSC"
        raw@original_units = raw@units
        raw@units          = "W m-2"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."rsuscs") then
      if ((isfilevar(in_file,"FSDSC")).and.(isfilevar(in_file,"FSNSC"))) then
        raw                = in_file->FSDSC
        raw                = in_file->FSDSC-in_file->FSNSC
        raw@standard_name  = "surface_upwelling_shortwave_flux_in_air_assuming_clear_sky"
        raw@original_name  = "FSDSC,FSNSC"
        raw@original_units = "W m-2"
        raw@units          = "W m-2"
        raw@history        = "FSDSC-FSNSC"
      end if
    end if
    if (ipcc_vars(v).eq."rldscs") then
      if ((isfilevar(in_file,"TS")).and.(isfilevar(in_file,"FLNSC"))) then
        raw                = in_file->TS
        raw                = sb*((in_file->TS)^4)-(in_file->FLNSC)
        raw@standard_name  = "surface_downwelling_longwave_flux_in_air_assuming_clear_sky"
        raw@original_name  = "TS,FLNSC"
        raw@original_units = "K,W m-2"
        raw@units          = "W m-2"
        raw@history        = "SB*TS^4-FLNSC"
      end if
    end if
    if (ipcc_vars(v).eq."rlutcs") then
      if (isfilevar(in_file,"FLUTC")) then
        raw                = in_file->FLUTC
        raw@standard_name  = "toa_outgoing_longwave_flux_assuming_clear_sky"
        raw@original_name  = "FLUTC"
        raw@original_units = raw@units
        raw@units          = "W m-2"
        raw@history        = "No change"
      end if
    end if
    if (ipcc_vars(v).eq."rsutcs") then
      if ((isfilevar(in_file,"SOLIN")).and.(isfilevar(in_file,"FSNTOAC"))) then
        raw                = in_file->SOLIN
        raw                = in_file->SOLIN-in_file->FSNTOAC
        raw@standard_name  = "toa_outgoing_shortwave_flux_assuming_clear_sky"
        raw@original_name  = "SOLIN,FSNTOAC"
        raw@original_units = "W m-2"
        raw@units          = "W m-2"
        raw@history        = "SOLIN-FSNTOAC"
      end if
    end if
    if (ipcc_vars(v).eq."clt") then
      if (isfilevar(in_file,"CLDTOT")) then
        raw                = in_file->CLDTOT
        raw                = raw*100
        raw@standard_name  = "cloud_area_fraction"
        raw@original_name  = "CLDTOT"
        raw@original_units = "fraction"
        raw@units          = "percent"
        raw@history        = "CLDTOT*100"
      end if
    end if
    if (ipcc_vars(v).eq."clwvi") then
      if (isfilevar(in_file,"TGCLDLWP")) then
        raw                = in_file->TGCLDLWP
        raw                = in_file->TGCLDLWP/g_to_kg
        raw@standard_name  = "atmosphere_cloud_condensed_water_content"
        raw@original_name  = "TGCLDLWP"
        raw@original_units = "g m-2"
        raw@units          = "kg m-2"
        raw@history        = "TGCLDLWP/g2kg"
      end if
    end if
    if (ipcc_vars(v).eq."clivi") then
      if (isfilevar(in_file,"TGCLDIWP")) then
        raw                = in_file->TGCLDIWP
        raw                = in_file->TGCLDIWP/g_to_kg
        raw@standard_name  = "atmosphere_cloud_ice_content"
        raw@original_name  = "TGCLDIWP"
        raw@original_units = "g m-2"
        raw@units          = "kg m-2"
        raw@history        = "TGCLDIWP/g2kg"
     end if
   end if
   if (ipcc_vars(v).eq."ta") then
      if (isfilevar(in_file, "TEMP")) then
        raw = in_file->TEMP
      end if
      if (isfilevar(in_file, "T")) then
        raw = in_file->T
      end if
      rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,False,1,tbot,phis)
      delete(raw)
      raw                = rawp
      delete(rawp)
      delete(raw&plev)
      raw&plev           = plev
      raw&plev		 = plev
      raw&plev           = plev*100.
      raw@standard_name  = "air_temperature"
      raw@original_name  = "T"
      raw@original_units = "K"
      raw@units          = "K"
      raw@history        = "Interpolated T with NCL 'vinth2p_ecmwf'"
    end if
   if (ipcc_vars(v).eq."taE") then
      if (isfilevar(in_file, "TEMP")) then
        raw = in_file->TEMP
      end if
      if (isfilevar(in_file, "T")) then
        raw = in_file->T
      end if
      rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,True,1,tbot,phis)
      delete(raw)
      raw                = rawp
      delete(rawp)
      delete(raw&plev)
      raw&plev           = plev
      raw&plev		 = plev
      raw&plev           = plev*100.
      raw@standard_name  = "air_temperature"
      raw@original_name  = "T"
      raw@original_units = "K"
      raw@units          = "K"
      raw@history        = "Interpolated T with NCL 'vinth2p_ecmwf', extrapolation below ground"
    end if
    if (ipcc_vars(v).eq."cl") then
      if (isfilevar(in_file, "CLOUD")) then
        raw                = in_file->CLOUD(:,::-1,:,:)
        raw                = raw*100
        delete(raw&lev)
        lev 		 = lev(::-1)*0.001
        ilev 		 = ilev(::-1)*0.001
        raw&lev		 = lev
        raw@standard_name  = "cloud_area_fraction_in_atmosphere_layer"
        raw@original_name  = "CLOUD"
        raw@original_units = "fraction"
        raw@history        = "CLOUD*100"
        raw@units          = "percent"
        ps@standard_name   = "surface_air_pressure"
        ps@original_name   = "PS"
        ps@original_units  = ps@units
      end if
    end if
    if (ipcc_vars(v).eq."ua") then
      if (isfilevar(in_file,"U")) then
        raw = in_file->U
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,False,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "eastward_wind"
        raw@original_name  = "U"
        raw@original_units = "m s-1"
        raw@units          = "m s-1"
        raw@history        = "Interpolated U with NCL 'vinth2p_ecmwf'"
      end if
    end if
    if (ipcc_vars(v).eq."uaE") then
      if (isfilevar(in_file,"U")) then
        raw = in_file->U
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,True,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "eastward_wind"
        raw@original_name  = "U"
        raw@original_units = "m s-1"
        raw@units          = "m s-1"
        raw@history        = "Interpolated U with NCL 'vinth2p_ecmwf', extrapolation below ground"
      end if
    end if
    if (ipcc_vars(v).eq."va") then
      if (isfilevar(in_file,"V")) then
        raw = in_file->V
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,False,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "northward_wind"
        raw@original_name  = "V"
        raw@original_units = "m s-1"
        raw@units          = "m s-1"
        raw@history        = "Interpolated V with NCL 'vinth2p_ecmwf'"
      end if
    end if
    if (ipcc_vars(v).eq."vaE") then
      if (isfilevar(in_file,"V")) then
        raw = in_file->V
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,True,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "northward_wind"
        raw@original_name  = "V"
        raw@original_units = "m s-1"
        raw@units          = "m s-1"
        raw@history        = "Interpolated V with NCL 'vinth2p_ecmwf', extrapolation below ground"
      end if
    end if
    if (ipcc_vars(v).eq."hus") then
      if (isfilevar(in_file,"Q")) then
        raw = in_file->Q
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,False,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "specific_humidity"
        raw@original_name  = "Q"
        raw@original_units = "kg kg-1"
        raw@units          = "1"
        raw@history        = "Interpolated Q with NCL 'vinth2p_ecmwf'"
      end if
    end if
    if (ipcc_vars(v).eq."husE") then
      if (isfilevar(in_file,"Q")) then
        raw = in_file->Q
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,True,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "specific_humidity"
        raw@original_name  = "Q"
        raw@original_units = "kg kg-1"
        raw@units          = "1"
        raw@history        = "Interpolated Q with NCL 'vinth2p_ecmwf', extrapolation below ground"
      end if
    end if
    if (ipcc_vars(v).eq."wap") then
      if (isfilevar(in_file,"OMEGA")) then
        raw = in_file->OMEGA
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,False,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "lagrangian_tendency_of_air_pressure"
        raw@original_name  = "OMEGA"
        raw@original_units = "Pa s-1"
        raw@units          = "Pa s-1"
        raw@history        = "Interpolated OMEGA with NCL 'vinth2p_ecmwf'"
      end if
    end if
    if (ipcc_vars(v).eq."wapE") then
      if (isfilevar(in_file,"OMEGA")) then
        raw = in_file->OMEGA
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,True,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "lagrangian_tendency_of_air_pressure"
        raw@original_name  = "OMEGA"
        raw@original_units = "Pa s-1"
        raw@units          = "Pa s-1"
        raw@history        = "Interpolated OMEGA with NCL 'vinth2p_ecmwf', extrapolation below ground"
      end if
    end if
    if (ipcc_vars(v).eq."zg") then
      if (isfilevar(in_file,"Z3")) then
        raw = in_file->Z3
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,False,-1,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "geopotential_height"
        raw@original_name  = "Z3"
        raw@original_units = "m"
        raw@units          = "m"
        raw@history        = "Interpolated Z3 with NCL 'vinth2p_ecmwf'"
      end if
    end if
    if (ipcc_vars(v).eq."zgE") then
      if (isfilevar(in_file,"Z3")) then
        raw = in_file->Z3
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,True,-1,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "geopotential_height"
        raw@original_name  = "Z3"
        raw@original_units = "m"
        raw@units          = "m"
        raw@history        = "Interpolated Z3 with NCL 'vinth2p_ecmwf', extrapolation below ground"
      end if
    end if
    if (ipcc_vars(v).eq."hur") then
      if (isfilevar(in_file,"RELHUM")) then
        raw = in_file->RELHUM
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,False,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "relative_humidity"
        raw@original_name  = "RELHUM"
        raw@original_units = "percent"
        raw@units          = "percent"
        raw@history        = "Interpolated RELHUM with NCL 'vinth2p_ecmwf'"
      end if
    end if
    if (ipcc_vars(v).eq."hurE") then
      if (isfilevar(in_file,"RELHUM")) then
        raw = in_file->RELHUM
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,True,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "relative_humidity"
        raw@original_name  = "RELHUM"
        raw@original_units = "percent"
        raw@units          = "percent"
        raw@history        = "Interpolated RELHUM with NCL 'vinth2p_ecmwf', extrapolation below ground"
      end if
    end if
    if (ipcc_vars(v).eq."tro3") then
      if (isfilevar(in_file,"O3")) then
        raw = in_file->O3
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,False,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "mole_fraction_of_o3_in_air"
        raw@original_name  = "O3VMR"
        raw@original_units = "fraction"
        raw@units          = "1e-9"
        raw@history        = "Interpolated O3 with NCL 'vinth2p_ecmwf'"
      end if
    end if
    if (ipcc_vars(v).eq."tro3E") then
      if (isfilevar(in_file,"O3")) then
        raw = in_file->O3
        rawp = vinth2p_ecmwf(raw,hyam,hybm,plev,ps,1,P0mb,1,True,0,tbot,phis)
        delete(raw)
        raw                = rawp
        delete(rawp)
        delete(raw&plev)
        raw&plev		 = plev
        raw&plev		 = plev*100.
        raw@standard_name  = "mole_fraction_of_o3_in_air"
        raw@original_name  = "O3VMR"
        raw@original_units = "fraction"
        raw@units          = "1e-9"
        raw@history        = "Interpolated O3 with NCL 'vinth2p_ecmwf', extrapolation below ground"
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
;--------------------------------------------------------------------
; Create output file
; [file name as per requirements document]
;--------------------------------------------------------------------
      nccf = "./CMIP/"+ipcc_vars(v)+"_"+table+"."+expid+".CCSM.atmm."+tval+".nc"
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
       if (ipcc_vars(v).eq."cl") then
         var_z           = (/ "lev" /)
         var_z_type      = (/ "double" /)
         var_z_bnds      = (/ "lev_bnds" /)
         var_z_bnds_type = (/ "double" /)
         dimNames = (/ "lon", "lat",  "lev", "bnds", "time" /)
         dimSizes = (/  nlon,  nlat,  nlev ,  nbnd, -1      /)
         dimUnlim = (/ False, False,  False, False, True    /)
         filedimdef(cf_file, dimNames, dimSizes, dimUnlim )
         filevardef(cf_file, var_x_bnds, var_x_bnds_type, (/"lon","bnds"/) )
         filevardef(cf_file, var_y_bnds, var_y_bnds_type, (/"lat","bnds"/) )
         filevardef(cf_file, var_z_bnds, var_z_bnds_type, (/"lev","bnds"/) )
         filevardef(cf_file, var_t_bnds, var_t_bnds_type, (/"time","bnds"/) )
         filevardef(cf_file, var_t, var_t_type, (/"time"/))
         filevardef(cf_file, var_z, var_z_type, (/"lev"/))
         lev@standard_name = "atmosphere_hybrid_sigma_pressure_coordinate"
         lev@axis          = "Z"
         lev@units	   = "1"
         lev@bounds        = "lev_bnds"
         lev@positive      = "down"
         lev@formula	   = "p(n,k,j,i) = a(k)*p0 + b(k)*ps(n,j,i)"
         lev@formula_terms = "p0: p0 a: a b: b ps: ps"
         p0                = P0mb*100.
         p0!0 	           = "ncl_scalar"
         p0@units          = "Pa"
         p0@long_name      = "reference pressure for hybrid sigma coordinate"
         a                 = new((/nlev/),double)
         a                 = (/hyam/)
         a!0               = "lev"
         a&lev             = lev
         a@long_name       = "hybrid sigma coordinate A coefficient for layer"
         del_FillValue(a)
         a_bnds            = new((/nlev,nbnd/),double)
         a_bnds!0          = "lev"
         a_bnds!1          = "bnds"
         a_bnds&lev        = lev
         a_bnds(:,0)       = (/hyai(nlev-1:0:-1)/)
         a_bnds(:,1)       = (/hyai(nlev:1:-1)/)
         a_bnds@long_name  = "hybrid sigma coordinate A coefficient for layer bounds"
         del_FillValue(a_bnds)
         b                 = new((/nlev/),double)
         b                 = (/hybm/)
         b!0               = "lev"
         b&lev             = lev
         b@long_name       = "hybrid sigma coordinate B coefficient for layer"
         del_FillValue(b)
         b_bnds            = new((/nlev,nbnd/),double)
         b_bnds!0               = "lev"
         b_bnds!1          = "bnds"
         b_bnds&lev             = lev
         b_bnds(:,0)       = (/hybi(nlev-1:0:-1)/)
         b_bnds(:,1)       = (/hybi(nlev:1:-1)/)
         b_bnds@long_name  = "hybrid sigma coordinate B coefficient for layer bounds"
         del_FillValue(b_bnds)
         lev_bnds               =  new((/nlev,nbnd/),double)
         lev_bnds(:,0)          =  (/ilev(nlev-1:0:-1)/)
         lev_bnds(:,1)          =  (/ilev(nlev:1:-1)/)
         lev_bnds@standard_name =  lev@standard_name
         lev_bnds@formula       =  lev@formula
         lev_bnds@formula_terms = "p0: p0 a: a_bnds b: b_bnds ps: ps"
         del_FillValue(lev_bnds)
         del_FillValue(ps)
         cf_file->p0 = p0
         cf_file->a  = a
         cf_file->b  = b
         cf_file->a_bnds   = a_bnds
         cf_file->b_bnds   = b_bnds
         cf_file->lev_bnds = lev_bnds
         cf_file->lev = lev
         cf_file->ps = ps
        else
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
        end if
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
          (ipcc_vars(v).eq."tasmax").or.\
          (ipcc_vars(v).eq."huss")) then
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
