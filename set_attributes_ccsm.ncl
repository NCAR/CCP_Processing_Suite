;********************************************************************************
; set_attributes
; set a whole bunch of mainly text attributes
;
undef("set_attributes_ccsm")
procedure set_attributes_ccsm(case,glbAtt,parent,year1,machine,forcing,expid)
;
begin
  n1 = integertochar(10)
;--------------------------------------------------------------------
; initialize arguments
;--------------------------------------------------------------------
  parent  = ""
  year1   = ""
  machine = ""
  forcing = ""
  expid   = ""
;--------------------------------------------------------------------
; experiment-invariant global attributes
;--------------------------------------------------------------------
  glbAtt@title          = "model output prepared for IPCC AR4"
  glbAtt@institution    = "NCAR (National Center for Atmospheric "+n1+\
                          "Research, Boulder, CO, USA)"
  glbAtt@source         = "CCSM3.0, version beta19 (2004): "+n1+\
                          "atmosphere: CAM3.0, T85L26;"+n1+\
                          "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                          "sea ice   : CSIM5.0, gx1v3;"+n1+\
                          "land      : CLM3.0, T85"
  glbAtt@contact        = "ccsm@ucar.edu"
  glbAtt@project_id     = "IPCC Fourth Assessment"
  glbAtt@Conventions    = "CF-1.0"
  glbAtt@references     = "Collins, W.D., et al., 2005:"+n1+\
                          " The Community Climate System Model, Version 3"+n1+\
                          " Journal of Climate"+n1+\
                          " "+n1+\
                          " Main website: http://www.ccsm.ucar.edu"
  glbAtt@acknowledgment = " Any use of CCSM data should acknowledge the contribution"+n1+\
                          " of the CCSM project and CCSM sponsor agencies with the "+n1+\
                          " following citation:"+n1+\
                          " 'This research uses data provided by the Community Climate"+n1+\
                          " System Model project (www.ccsm.ucar.edu), supported by the"+n1+\
                          " Directorate for Geosciences of the National Science Foundation"+n1+\
                          " and the Office of Biological and Environmental Research of"+n1+\
                          " the U.S. Department of Energy.'"+n1+\
                          "In addition, the words 'Community Climate System Model' and"+n1+\
                          " 'CCSM' should be included as metadata for webpages referencing"+n1+\
                          " work using CCSM data or as keywords provided to journal or book"+n1+\
                          "publishers of your manuscripts."+n1+\
                          "Users of CCSM data accept the responsibility of emailing"+n1+\
                          " citations of publications of research using CCSM data to"+n1+\
                          " ccsm@ucar.edu."+n1+\
                          "Any redistribution of CCSM data must include this data"+n1+\
                          " acknowledgement statement."
;--------------------------------------------------------------------
; IPCC experiment ID - changes with case name
; Allowed values                          Abbreviation
;--------------------------------------------------------------------
; pre-industrial control experiment       PIcntrl
; present-day control experiment          PDcntrl
; climate of the 20th Century experiment  20C3M
; committed climate change experiment     Commit
; SRES A2 experiment                      SRESA2
; 720 ppm stabilization experiment        SRESA1B
; 550 ppm stabilization experiment        SRESB1
; 1%/year CO2 increase experiment, to 2x  1%_to2x
; 1%/year CO2 increase experiment, to 4x  1%_to4x
; slab ocean control experiment           Slabcntl
; 2xCO2 equilibrium experiment            2xCO2
; AMIP experiment                         AMIP
;--------------------------------------------------------------------
; CCSP level 2 (CO2 550 ppm @ 2300) stabilization experiment (CCSPL2)
;--------------------------------------------------------------------
  if (case.eq."b30.240a") then
    parent  = "b30.030a"
    year1   = "2000"
    machine = "jaguar.ccs.ornl.gov"
    forcing = "ozone forcing    : ozone_WLEV2-53_128x64_L18_1990-2100_c080701.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_WLEV2-53_1765-2300_c080630.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc (2005 value ONLY)"+n1+\
              "GHGs             : ghg_WLEV2-53_1765-2300_c080708.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_WLEV2-53_256x128L2_1990-2100_c080701.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "CCSP level 2 (CO2 550 ppm @ 2300) stabilization experiment (CCSPL2)"
    expid   = "CCSPL2_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; black carbon (bc)
;--------------------------------------------------------------------
  if (case.eq."b30.053a")then
    parent  = "b30.020"
    year1   = "360"
    machine = "seaborg.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 only"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_bc_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
  if (case.eq."b30.053b")then
    parent  = "b30.020"
    year1   = "380"
    machine = "seaborg.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 ONLY"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_bc_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
  if (case.eq."b30.053c")then
    parent  = "b30.020"
    year1   = "400"
    machine = "seaborg.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 ONLY"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 3
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_bc_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
  if (case.eq."b30.053d")then
    parent  = "b30.020"
    year1   = "420"
    machine = "seaborg.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 ONLY"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 4
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_bc_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
;--------------------------------------------------------------------
; A1FI experiment
;--------------------------------------------------------------------
  if (case.eq."b30.099a")then
    parent  = "b30.020"
    year1   = "360"
    machine = "jaguar.ornl.ccs.gov"
    forcing = "ozone forcing    : A1FI.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1FI_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1FI_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1FI_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "fossil-fuel intensive experiment (A1FI)"
    expid   = "A1FI_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; Anthropogenic-only 20th century experiments
;--------------------------------------------------------------------
  if (case.eq."b30.090a")then
    parent  = "b30.020"
    year1   = "360"
    machine = "phoenix.ornl.ccs.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at at 1366.5 W m-2 (1870 Lean et.al. value)"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "anthropogenic forcings only for the 20th Century (20C3M_anthro)"
    expid   = "20C3M_anthro_"+glbAtt@realization
  end if
  if (case.eq."b30.090b") then
    parent  = "b30.020"
    year1   = "380"
    machine = "phoenix.ornl.ccs.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at at 1366.5 W m-2 (1870 Lean et.al. value)"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "anthropogenic forcings only for the 20th Century (20C3M_anthro)"
    expid   = "20C3M_anthro_"+glbAtt@realization
  end if
  if (case.eq."b30.090c") then
    parent  = "b30.020"
    year1   = "400"
    machine = "phoenix.ornl.ccs.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at at 1366.5 W m-2 (1870 Lean et.al. value)"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 3
    glbAtt@experiment_id  = "anthropogenic forcings only for the 20th Century (20C3M_anthro)"
    expid   = "20C3M_anthro_"+glbAtt@realization
  end if
  if (case.eq."b30.090d") then
    parent  = "b30.020"
    year1   = "420"
    machine = "phoenix.ornl.ccs.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at at 1366.5 W m-2 (1870 Lean et.al. value)"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 4
    glbAtt@experiment_id  = "anthropogenic forcings only for the 20th Century (20C3M_anthro)"
    expid   = "20C3M_anthro_"+glbAtt@realization
  end if
  if (case.eq."b30.090e") then
    parent  = "b30.020"
    year1   = "440"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at at 1366.5 W m-2 (1870 Lean et.al. value)"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 5
    glbAtt@experiment_id  = "anthropogenic forcings only for the 20th Century (20C3M_anthro)"
    expid   = "20C3M_anthro_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; Natural-forcings-only 20th century experiments
;--------------------------------------------------------------------
  if (case.eq."b30.072a")then
    parent  = "b30.020"
    year1   = "360"
    machine = "phoenix.ornl.ccs.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c060324.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "natural forcings only for the 20th Century (NatOnly)"
    expid   = "20C3M_nat_"+glbAtt@realization
  end if
  if (case.eq."b30.072b") then
    parent  = "b30.020"
    year1   = "380"
    machine = "phoenix.ornl.ccs.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c060324.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "natural forcings only for the 20th Century (20C3M_nat)"
    expid   = "20C3M_nat_"+glbAtt@realization
  end if
  if (case.eq."b30.072c") then
    parent  = "b30.020"
    year1   = "400"
    machine = "phoenix.ornl.ccs.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c060324.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"
    glbAtt@realization    = 3
    glbAtt@experiment_id  = "natural forcings only for the 20th Century (20C3M_nat)"
    expid   = "20C3M_nat_"+glbAtt@realization
  end if
  if (case.eq."b30.072d") then
    parent  = "b30.020"
    year1   = "420"
    machine = "phoenix.ornl.ccs.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c060324.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"
    glbAtt@realization    = 4
    glbAtt@experiment_id  = "natural forcings only for the 20th Century (20C3M_nat)"
    expid   = "20C3M_nat_"+glbAtt@realization
  end if
  if (case.eq."b30.072e") then
    parent  = "b30.020"
    year1   = "440"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c060324.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"
    glbAtt@realization    = 5
    glbAtt@experiment_id  = "natural forcings only for the 20th Century (20C3M_nat)"
    expid   = "20C3M_nat_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; slab ocean control experiment           Slabcntl
;--------------------------------------------------------------------
  if (case.eq."eul128x256_d50som")then
    parent  = "none"
    year1   = "0"
    machine = "bluesky.ucar.edu"
    forcing = "See "+n1+\
              "http://www.ccsm.ucar.edu/models/atm-cam/docs/description/node28.html"+n1+\
              "for details."
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "slab ocean control experiment"
    glbAtt@source         = "CCSM3.0, version beta08 (2003): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : N/A"+n1+\
                            "sea ice   : N/A"+n1+\
                            "land      : CLM3.0, T85"
    expid   = "Slabcntl_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; 2xCO2 equilibrium experiment            2xCO2
;--------------------------------------------------------------------
  if (case.eq."eul128x256_d50somd")then
    parent  = "none"
    year1   = "0"
    machine = "bluesky.ucar.edu"
    forcing = "See "+n1+\
              "http://www.ccsm.ucar.edu/models/atm-cam/docs/description/node28.html"+n1+\
              "for details."
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "2xCO2 equilibrium experiment"
    glbAtt@source         = "CCSM3.0, version beta08 (2003): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : N/A"+n1+\
                            "sea ice   : N/A"+n1+\
                            "land      : CLM3.0, T85"
    expid   = "2XCO2_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; present-day control experiment          PDcntrl
;--------------------------------------------------------------------
  if (case.eq."b30.009")then
    parent  = "none"
    year1   = "0"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : pcmdio3.r8.64x1_L60_clim_c970515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c030529.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : none"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 2000 values only"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 2000 values only"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040130.nc, 2000 values only"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "present-day control experiment"
    glbAtt@source         = "CCSM3.0, version beta08 (2003): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    expid   = "PDcntrl_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; 1% CO2 to 4x (and 2x - MESSY!)
;--------------------------------------------------------------------
  if (case.eq."b30.026b")then
    parent  = "b30.026.ES01"
    year1   = "0549"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : pcmdio3.r8.64x1_L60_clim_c970515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c030529.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : none"+n1+\
              "solar forcing    : constant"+n1+\
              "GHGs             : none"+n1+\
              "GHG loss rates   : none"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : none"+n1+\
              "oxidants         : none"+n1+\
              "SOx emissions    : none"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "1%/year CO2 increase experiment (to quadrupling)"
    glbAtt@source         = "CCSM3.0, version beta08 (2003): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    expid   = "1%_to4x_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; 1% CO2 to 4x (and 2x - MESSY!)
;--------------------------------------------------------------------
  if (case.eq."b30.026b.ES01")then
    parent  = "b30.009"
    year1   = "0400"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : pcmdio3.r8.64x1_L60_clim_c970515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c030529.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : none"+n1+\
              "solar forcing    : constant"+n1+\
              "GHGs             : none"+n1+\
              "GHG loss rates   : none"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : none"+n1+\
              "oxidants         : none"+n1+\
              "SOx emissions    : none"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "1%/year CO2 increase experiment (to quadrupling)"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    expid   = "1%_to2x_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; 1% CO2 to 4x (and 2x - MESSY!)
;--------------------------------------------------------------------
  if (case.eq."b30.026a.ES01")then
    parent  = "b30.026"
    year1   = "0479"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : pcmdio3.r8.64x1_L60_clim_c970515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c030529.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : none"+n1+\
              "solar forcing    : constant"+n1+\
              "GHGs             : none"+n1+\
              "GHG loss rates   : none"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : none"+n1+\
              "oxidants         : none"+n1+\
              "SOx emissions    : none"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "1%/year CO2 increase experiment (to doubling)"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    expid   = "1%_to2x_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; pre-industrial control experiment       PIcntrl
;--------------------------------------------------------------------
  if (case.eq."b30.020")then
    parent  = "b30.017"
    year1   = "280"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870cyc_c040123.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-1999_c040114.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040130.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "pre-industrial control experiment"
    expid   = "PIcntrl_"+glbAtt@realization
  end if
  if (case.eq."b30.020.ES01")then
    parent  = "b30.020"
    year1   = "300"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870cyc_c040123.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-1999_c040114.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040130.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "pre-industrial control experiment"
    expid   = "PIcntrl_"+glbAtt@realization
  end if
  if (case.eq."b30.020.ES02")then
    parent  = "b30.020.ES01"
    year1   = "430"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870cyc_c040123.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-1999_c040114.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040130.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "pre-industrial control experiment"
    expid   = "PIcntrl_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; single-forcing 20th Century experiments (20C3M_XX)
;--------------------------------------------------------------------
; sulfate aerosol direct effects (sd)
;--------------------------------------------------------------------
  if (case.eq."b30.102a")then
    parent  = "b30.020"
    year1   = "360"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_sd_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
  if (case.eq."b30.102b")then
    parent  = "b30.020"
    year1   = "380"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    glbAtt@source         = "CCSM3.0.1, version beta24"
    expid   = "20C3M_sd_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; well-mixed greenhouse gasses (gg)
;--------------------------------------------------------------------
  if (case.eq."b30.216a")then
    parent  = "b30.020"
    year1   = "360"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 ONLY"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_gg_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
  if (case.eq."b30.216b")then
    parent  = "b30.020"
    year1   = "380"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 ONLY"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_gg_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
;--------------------------------------------------------------------
; solar irradiance (sl)
;--------------------------------------------------------------------
  if (case.eq."b30.217a")then
    parent  = "b30.020"
    year1   = "360"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_sl_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
  if (case.eq."b30.217b")then
    parent  = "b30.020"
    year1   = "380"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_sl_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
;--------------------------------------------------------------------
; volcanic aerosols (vl)
;--------------------------------------------------------------------
  if (case.eq."b30.218a")then
    parent  = "b30.020"
    year1   = "360"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 ONLY"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_vl_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
  if (case.eq."b30.218b")then
    parent  = "b30.020"
    year1   = "360"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 ONLY"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_vl_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
;--------------------------------------------------------------------
; tropospheric and stratospheric ozone (oz; to + so)
;--------------------------------------------------------------------
  if (case.eq."b30.219a")then
    parent  = "b30.020"
    year1   = "360"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 ONLY"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_oz_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
  if (case.eq."b30.219b")then
    parent  = "b30.020"
    year1   = "380"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 ONLY"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "volcanic forcing : none"+n1+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_oz_"+glbAtt@realization
    glbAtt@source         = "CCSM3.0.1, version beta24"
  end if
;--------------------------------------------------------------------
; sulfate aerosol indirect effects (si)
;--------------------------------------------------------------------
  if (case.eq."b30.220a")then
    parent  = "b30.020"
    year1   = "360"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_si_"+glbAtt@realization
  end if
  if (case.eq."b30.220b")then
    parent  = "b30.020"
    year1   = "360"
    machine = "bassi.nersc.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_si_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; climate of the 20th Century experiment (20C3M)
;--------------------------------------------------------------------
  if (case.eq."b30.030a")then
    parent  = "b30.020"
    year1   = "360"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_"+glbAtt@realization
  end if
  if (case.eq."b30.030b") then
    parent  = "b30.020"
    year1   = "380"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_"+glbAtt@realization
  end if
  if (case.eq."b30.030c") then
    parent  = "b30.020"
    year1   = "400"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 3
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_"+glbAtt@realization
  end if
  if (case.eq."b30.030d") then
    parent  = "b30.020"
    year1   = "420"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 4
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_"+glbAtt@realization
  end if
  if (case.eq."b30.030e") then
    parent  = "b30.020"
    year1   = "440"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 5
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_"+glbAtt@realization
  end if
  if (case.eq."b30.030b.ES01")then
    parent  = "b30.020.ES01"
    year1   = "380"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 6
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_"+glbAtt@realization
  end if
  if (case.eq."b30.030f.ES01")then
    parent  = "b30.020.ES01"
    year1   = "410"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 7
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_"+glbAtt@realization
  end if
  if (case.eq."b30.030g.ES01")then
    parent  = "b30.020.ES01"
    year1   = "460"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 8
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_"+glbAtt@realization
  end if
  if (case.eq."b30.030h.ES01")then
    parent  = "b30.020b.ES01"
    year1   = "540"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : VolcanicMass_1870-1999_64x1_L18_c040123.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 9
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    expid   = "20C3M_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; 720 ppm stabilization experiment (SRES A1B)
;--------------------------------------------------------------------
  if (case.eq."b30.040a") then
    parent  = "b30.030a"
    year1   = "2000"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.040b") then
    parent  = "b30.030b"
    year1   = "2000"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.040c") then
    parent  = "b30.030c"
    year1   = "2000"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 3
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.040d") then
    parent  = "b30.030d"
    year1   = "2000"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 4
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.040e") then
    parent  = "b30.030e"
    year1   = "2000"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 5
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.040g") then
    parent  = "b30.030g.ES01"
    year1   = "2000"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.044a") then
    parent  = "b30.040a"
    year1   = "2100"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.044b") then
    parent  = "b30.040b"
    year1   = "2100"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.044c") then
    parent  = "b30.040c"
    year1   = "2100"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 3
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.044d") then
    parent  = "b30.040d"
    year1   = "2100"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 4
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.044e") then
    parent  = "b30.040e"
    year1   = "2100"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 5
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.040b.ES01") then
    parent  = "b30.030b.ES01"
    year1   = "2000"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 6
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.040f.ES01")then
    parent  = "b30.030f.ES01"
    year1   = "2000"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 7
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.040g.ES01")then
    parent  = "b30.030g.ES01"
    year1   = "2000"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 8
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.044b.ES01") then
    parent  = "b30.040b.ES01"
    year1   = "2100"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 6
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.044f.ES01")then
    parent  = "b30.040f.ES01"
    year1   = "2100"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 7
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
  if (case.eq."b30.044g.ES01")then
    parent  = "b30.040g.ES01"
    year1   = "2100"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 8
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    expid   = "SRESA1B_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; 550 ppm stabilization experiment (SRES B1)
;--------------------------------------------------------------------
  if (case.eq."b30.041a") then
    parent  = "b30.030a"
    year1   = "2000"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.041b") then
    parent  = "b30.030b"
    year1   = "2000"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.041c") then
    parent  = "b30.030c"
    year1   = "2000"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 3
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.041d") then
    parent  = "b30.030d"
    year1   = "2000"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 4
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.041e") then
    parent  = "b30.030e"
    year1   = "2000"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 5
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.041g") then
    parent  = "b30.030g.ES01"
    year1   = "2000"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 9
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.045a") then
    parent  = "b30.041a"
    year1   = "2100"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.045b") then
    parent  = "b30.041b"
    year1   = "2100"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.045c") then
    parent  = "b30.041c"
    year1   = "2100"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 3
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.045d") then
    parent  = "b30.041d"
    year1   = "2100"
    machine = "cheetah.ccs.ornl.gov"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 4
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.045e") then
    parent  = "b30.041e"
    year1   = "2100"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 5
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.041b.ES01")then
    parent  = "b30.030b.ES01"
    year1   = "2000"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 6
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.041f.ES01")then
    parent  = "b30.030f.ES01"
    year1   = "2000"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 7
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.041g.ES01")then
    parent  = "b30.030g.ES01"
    year1   = "2000"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 8
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.045b.ES01")then
    parent  = "b30.041b.ES01"
    year1   = "2100"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 6
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.045f.ES01")then
    parent  = "b30.041f.ES01"
    year1   = "2100"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 7
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
  if (case.eq."b30.045g.ES01")then
    parent  = "b30.041g.ES01"
    year1   = "2100"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : B1.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_B1_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_B1_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_B1_128x256_L2_1990-2100_c040608.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 8
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
    expid   = "SRESB1_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; SRESA2 experiment
;--------------------------------------------------------------------
  if (case.eq."b30.042a") then
    parent  = "b30.030a"
    year1   = "2000"
    machine = "seaborg.nersc.gov"
    forcing = "ozone forcing    : A2.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A2_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A2_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "SRES A2 experiment"
    expid   = "SRESA2_"+glbAtt@realization
  end if
  if (case.eq."b30.042b") then
    parent  = "b30.030b"
    year1   = "2000"
    machine = "seaborg.nersc.gov"
    forcing = "ozone forcing    : A2.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A2_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A2_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "SRES A2 experiment"
    expid   = "SRESA2_"+glbAtt@realization
  end if
  if (case.eq."b30.042c") then
    parent  = "b30.030c"
    year1   = "2000"
    machine = "seaborg.nersc.gov"
    forcing = "ozone forcing    : A2.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A2_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A2_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 3
    glbAtt@experiment_id  = "SRES A2 experiment"
    expid   = "SRESA2_"+glbAtt@realization
  end if
  if (case.eq."b30.042d") then
    parent  = "b30.030d"
    year1   = "2000"
    machine = "seaborg.nersc.gov"
    forcing = "ozone forcing    : A2.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A2_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A2_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 4
    glbAtt@experiment_id  = "SRES A2 experiment"
    expid   = "SRESA2_"+glbAtt@realization
  end if
  if (case.eq."b30.042e") then
    parent  = "b30.030e"
    year1   = "2000"
    machine = "seaborg.nersc.gov"
    forcing = "ozone forcing    : A2.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A2_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A2_128x256_L2_1990-2100_c040608.nc"
    glbAtt@realization    = 5
    glbAtt@experiment_id  = "SRES A2 experiment"
    expid   = "SRESA2_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; committed climate change experiment
;--------------------------------------------------------------------
  if (case.eq."b30.036a") then
    parent  = "b30.030a"
    year1   = "2000"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_2000_c040526.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 1
    glbAtt@experiment_id  = "committed climate change experiment"
    expid   = "Commit_"+glbAtt@realization
  end if
  if (case.eq."b30.036b") then
    parent  = "b30.030b"
    year1   = "2000"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_2000_c040526.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 2
    glbAtt@experiment_id  = "committed climate change experiment"
    expid   = "Commit_"+glbAtt@realization
  end if
  if (case.eq."b30.036c") then
    parent  = "b30.030c"
    year1   = "2000"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_2000_c040526.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 3
    glbAtt@experiment_id  = "committed climate change experiment"
    expid   = "Commit_"+glbAtt@realization
  end if
  if (case.eq."b30.036d") then
    parent  = "b30.030d"
    year1   = "2000"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_2000_c040526.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 4
    glbAtt@experiment_id  = "committed climate change experiment"
    expid   = "Commit_"+glbAtt@realization
  end if
  if (case.eq."b30.036e") then
    parent  = "b30.030e"
    year1   = "2000"
    machine = "bluesky.ucar.edu"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_2000_c040526.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@realization    = 5
    glbAtt@experiment_id  = "committed climate change experiment"
    expid   = "Commit_"+glbAtt@realization
  end if
  if (case.eq."b30.036b.ES01")then
    parent  = "b30.030b.ES01"
    year1   = "2000"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_2000_c040526.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 6
    glbAtt@experiment_id  = "committed climate change experiment"
    expid   = "Commit_"+glbAtt@realization
  end if
  if (case.eq."b30.036f.ES01")then
    parent  = "b30.030f.ES01"
    year1   = "2000"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_2000_c040526.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 7
    glbAtt@experiment_id  = "committed climate change experiment"
    expid   = "Commit_"+glbAtt@realization
  end if
  if (case.eq."b30.036g.ES01")then
    parent  = "b30.030g.ES01"
    year1   = "2000"
    machine = "Earth Simulator Center, JAMSTEC"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_2000_c040526.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "volcanic forcing : none"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
    glbAtt@source         = "CCSM3.0, version vector05 (2004): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    glbAtt@realization    = 8
    glbAtt@experiment_id  = "committed climate change experiment"
    expid   = "Commit_"+glbAtt@realization
  end if
;--------------------------------------------------------------------
; Check to make sure we got values!
;--------------------------------------------------------------------
  if (expid.eq."") then
    print(case+" NOT DEFINED IN set_attributes_ccsm.ncl *PLEASE CORRECT*")
    status_exit(1)
  end if
  if (parent.eq."") then
    print(case+" NOT DEFINED IN set_attributes_ccsm.ncl *PLEASE CORRECT*")
    status_exit(1)
  end if
  if (machine.eq."") then
    print(case+" NOT DEFINED IN set_attributes_ccsm.ncl *PLEASE CORRECT*")
    status_exit(1)
  end if
  if (forcing.eq."") then
    print(case+" NOT DEFINED IN set_attributes_ccsm.ncl *PLEASE CORRECT*")
    status_exit(1)
  end if
  if (year1.eq."") then
    print(case+" NOT DEFINED IN set_attributes_ccsm.ncl *PLEASE CORRECT*")
    status_exit(1)
  end if
;--------------------------------------------------------------------
; experiment-invariant global attributes
;--------------------------------------------------------------------
  glbAtt@history        = "Created from CCSM3 case "+case+n1+\
                          " by "+systemfunc("whoami")+"@ucar.edu"+n1+\
                          " on "+systemfunc("date")+n1+\
                          " "+n1+\
                          " For all data, added IPCC requested metadata"
  glbAtt@comment        = "This simulation was initiated from year "+year1+" of "+n1+\
                          " CCSM3 model run "+parent+" and executed on "+n1+\
                          " hardware "+machine+". The input external forcings are"+n1+\
                          forcing+n1+\
                           " Physical constants used for derived data:"+n1+\
                           " Lv (latent heat of evaporation): 2.501e6 J kg-1"+n1+\
                           " Lf (latent heat of fusion     ): 3.337e5 J kg-1"+n1+\
                           " r[h2o] (density of water      ): 1000 kg m-3"+n1+\
                           " g2kg   (grams to kilograms    ): 1000 g kg-1"+n1+\
                          " "+n1+\
                          " Integrations were performed by NCAR and CRIEPI with support"+n1+\
                          " and facilities provided by NSF, DOE, MEXT and ESC/JAMSTEC."
end
;********************************************************************************
; del_FillValue
; get rid of a bogus FillValue attribute from x
;
undef("del_FillValue")
procedure del_FillValue (x)
; usage:     del_FillValue(T)
;            del_FillValue(ps)
begin
    if (isatt(x,"_FillValue")) then
        delete(x@_FillValue)
    end if
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
