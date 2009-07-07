;********************************************************************************
; attributes
; set a whole bunch of mainly text attributes
;
undef("attributes_ccsm")
procedure attributes_ccsm(case,glbAtt,expid)
;
begin
  n1 = integertochar(10)
;--------------------------------------------------------------------
; initialize arguments
;--------------------------------------------------------------------
  expid   = ""
;--------------------------------------------------------------------
; environment variables
;--------------------------------------------------------------------
  hardware    = getenv("HARDWARE")
  location    = getenv("LOCATION")
  realization = getenv("REALIZN")
  resolution  = getenv("RES")
  parent      = getenv("PARENT")
  branch      = getenv("BRANCH")
  description = getenv("DESCRIP")
;--------------------------------------------------------------------
; experiment-invariant global attributes
;--------------------------------------------------------------------
  glbAtt@title          = "model output prepared for IPCC AR4"
  glbAtt@institution    = "NCAR (National Center for Atmospheric "+n1+\
                          "Research, Boulder, CO, USA)"
  if (hardware.eq."JAMSTEC") then
    machine = "Earth Simulator Center, JAMSTEC"
    glbAtt@source = "CCSM3.0, version vector05 (2004): "+n1+\
                    "atmosphere: CAM3.0, T85L26;"+n1+\
                    "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                    "sea ice   : CSIM5.0, gx1v3;"+n1+\
                    "land      : CLM3.0, T85"
  else
    machine = hardware
    glbAtt@source = "CCSM3.0, version beta19 (2004): "+n1+\
                    "atmosphere: CAM3.0, T85L26;"+n1+\
                    "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                    "sea ice   : CSIM5.0, gx1v3;"+n1+\
                    "land      : CLM3.0, T85"
  end if
  glbAtt@realization    = realization
  expid                 = description+"_"+glbAtt@realization
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
; Define forcings - common to description:
;--------------------------------------------------------------------
  if (description.eq."1%_to2x") then
    glbAtt@experiment_id  = "1%/year CO2 increase experiment (to doubling)"
    glbAtt@source         = "CCSM3.0, version beta08 (2003): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    forcing = "ozone forcing    : pcmdio3.r8.64x1_L60_clim_c970515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c030529.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"
  end if
  if (description.eq."1%_to4x") then
    glbAtt@experiment_id  = "1%/year CO2 increase experiment (to quadrupling)"
    glbAtt@source         = "CCSM3.0, version beta08 (2003): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    forcing = "ozone forcing    : pcmdio3.r8.64x1_L60_clim_c970515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c030529.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"
  end if
  if (description.eq."20C3M") then
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
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
  end if
  if (description.eq."20C3M_anthro") then
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    glbAtt@source         = "CCSM3.0.1, version beta24"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : Fixed at at 1366.5 W m-2 (1870 Lean et.al. value)"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"
  end if
  if (description.eq."20C3M_nat") then
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    glbAtt@source         = "CCSM3.0.1, version beta24"
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
  end if
  if (description.eq."20C3M_bc") then
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    glbAtt@source         = "CCSM3.0.1, version beta24"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : carbonscaling_1870-2000_c040225.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 only"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
  end if
  if (description.eq."20C3M_gg") then
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    glbAtt@source         = "CCSM3.0.1, version beta24"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 ONLY"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
  end if
  if (description.eq."20C3M_oz") then
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    glbAtt@source         = "CCSM3.0.1, version beta24"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 1870 ONLY"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
  end if
  if (description.eq."20C3M_sd") then
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    glbAtt@source         = "CCSM3.0.1, version beta24"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
  end if
  if (description.eq."20C3M_sl") then
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    glbAtt@source         = "CCSM3.0.1, version beta24"
    forcing = "ozone forcing    : mozart.o3.128x64_L18_1870-2000_c040515.nc, 1870 ONLY"+n1+\
              "carbon scaling   : Fixed at 0.3"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 1870 ONLY"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc, 1870 ONLY"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"
  end if
  if (description.eq."20C3M_vl") then
    glbAtt@experiment_id  = "climate of the 20th Century experiment (20C3M)"
    glbAtt@source         = "CCSM3.0.1, version beta24"
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
  end if
  if (description.eq."2XCO2") then
    glbAtt@experiment_id  = "2xCO2 equilibrium experiment"
    glbAtt@source         = "CCSM3.0, version beta08 (2003): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : N/A"+n1+\
                            "sea ice   : N/A"+n1+\
                            "land      : CLM3.0, T85"
    forcing = "See "+n1+\
              "http://www.ccsm.ucar.edu/models/atm-cam/docs/descriptiontion/node28.html"+n1+\
              "for details."
  end if
  if (description.eq."AMIP") then
    glbAtt@experiment_id  = "AMIP experiment"
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
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040321.nc"+n1+\
              "SSTs             : Hurrell, J. W., J. J. Hack, D. Shea, J. Caron, "+n1+\
              " and J. Rosinski, 2005: A new sea surface  temperature and sea ice"+n1+\
              " boundary data set for the Community Atmospheric Model.  J. Climate, submitted."
  end if
  if (description.eq."CCSP_Level2") then
    glbAtt@experiment_id  = "CCSP level 2 (CO2 550 ppm @ 2300) stabilization experiment (CCSP_Level2)"
    forcing = "ozone forcing    : ozone_WLEV2-53_128x64_L18_1990-2100_c080701.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_WLEV2-53_1765-2300_c080630.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc (2005 value ONLY)"+n1+\
              "GHGs             : ghg_WLEV2-53_1765-2300_c080708.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_WLEV2-53_256x128L2_1990-2100_c080701.nc"
  end if  
  if (description.eq."CCSP_Level3") then
    glbAtt@experiment_id  = "CCSP level 3 (CO2 650 ppm @ 2300) stabilization experiment (CCSP_Level3)"
    forcing = "ozone forcing    : ozone_WLEV3-53_128x64_L18_1990-2100_c080701.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_WLEV3-53_1765-2300_c080630.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc (2005 value ONLY)"+n1+\
              "GHGs             : ghg_WLEV3-53_1765-2300_c080708.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_WLEV3-53_256x128L2_1990-2100_c080701.nc"
  end if  
  if (description.eq."CCSP_Level4") then
    glbAtt@experiment_id  = "CCSP level 3 (CO2 750 ppm @ 2300) stabilization experiment (CCSP_Level4)"
    forcing = "ozone forcing    : ozone_WLEV4-53_128x64_L18_1990-2100_c080701.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_WLEV4-53_1765-2300_c080630.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc (2005 value ONLY)"+n1+\
              "GHGs             : ghg_WLEV4-53_1765-2300_c080708.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_WLEV4-53_256x128L2_1990-2100_c080701.nc"
  end if
  if (description.eq."CCSP_Overshoot") then
    glbAtt@experiment_id  = "CCSP overshoot experiment (CCSP_Overshoot)"
    glbAtt@experiment_id  = "CCSP overshoot experiment"
    forcing = "ozone forcing    : ozone_Overshoot_128x64_L18_1990-2100_c070425.nc"+n1+\
              "                 : ozone_Overshoot_128x64_L18_2100-2300_c081210.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_Overshoot_1765-2300_c070515.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc (2005 value ONLY)"+n1+\
              "GHGs             : ghg_Overshoot_1765-2300_c070425.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_Overshoot_256x128L2_1990-2100_c070608.nc"+n1+\
              "                 : SOx_emissions_Overshoot_256x128L2_2100-2300_c081208.nc"
  end if
  if (description.eq."CCSP_Reference") then
    glbAtt@experiment_id  = "CCSP reference (CO2 1000 ppm @ 2300) experiment (CCSP_Reference)"
    forcing = "ozone forcing    : ozone_Reference_128x64_L18_1990-2100_c070425.nc"+n1+\
              "                 : ozone_Reference_128x64_L18_2100-2300_c081210.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_Reference_1765-2300_c070515.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc (2005 value ONLY)"+n1+\
              "GHGs             : ghg_Reference_1765-2300_c070425.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_Reference_256x128L2_1990-2100_c070608.nc"+n1+\
              "                 : SOx_emissions_Reference_256x128L2_2100-2300_c081208.nc"
  end if
  if (description.eq."CCSP_Stab450") then
    glbAtt@experiment_id  = "CCSP stabilization (CO2 max 450ppm) experiment (CCSP_Stab450)"
    forcing = "ozone forcing    : ozone_Stab450_128x64_L18_1990-2100_c070425.nc"+n1+\
              "                 : ozone_Stab450_128x64_L18_2100-2300_c081210.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_Stab450_1765-2300_c070515.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc (2005 value ONLY)"+n1+\
              "GHGs             : ghg_Stab450_1765-2300_c070425.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_Stab450_256x128L2_1990-2100_c070608.nc"+n1+\
              "                 : SOx_emissions_Stab450_256x128L2_2100-2300_c081208.nc"
  end if
  if (description.eq."Commit") then
    glbAtt@experiment_id  = "committed climate change experiment"
    if (resolution.eq."T85_gx1") then
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
    end if
    if (resolution.eq."T42_gx1") then
      forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040521.nc, 2000 value only"+n1+\
                "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
                "aerosol MMR      : AerosolMass_V_64x128_clim_c031022.nc"+n1+\
                "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc, 2000 value only"+n1+\
                "solar forcing    : scon_lean_1870-2100_c040123.nc, 2000 value only"+n1+\
                "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc, 2000 values only"+n1+\
                "GHG loss rates   : noaamisc.r8.nc"+n1+\
                "DMS emissions    : DMS_emissions_64x128_c030722.nc"+n1+\
                "oxidants         : oxid_64x128_L26_clim_c050816.nc"+n1+\
                "SOx emissions    : SOx_emissions_A1_64x128_L2_1990-2100_c050912.nc, 2000 values only"
    end if
  end if
  if (description.eq."Geoengineering") then
  end if
  if (description.eq."PDcntrl") then
    glbAtt@experiment_id  = "present-day control experiment"
    glbAtt@source         = "CCSM3.0, version beta08 (2003): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : POP1.4.3 (modified), gx1v3"+n1+\
                            "sea ice   : CSIM5.0, gx1v3;"+n1+\
                            "land      : CLM3.0, T85"
    forcing = "ozone forcing    : pcmdio3.r8.64x1_L60_clim_c970515.nc"+n1+\
              "aerosol optics   : AerosolOptics_c030529.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "solar forcing    : scon_lean_1870-2100_c040123.nc, 2000 values only"+n1+\
              "GHGs             : ghg_1870_2100_c040122.nc, 2000 values only"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_128x256_L2_1850-2000_c040130.nc, 2000 values only"
  end if
  if (description.eq."PIcntrl") then
    glbAtt@experiment_id  = "pre-industrial control experiment"
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
  end if
  if (description.eq."Slabcntl") then
    glbAtt@experiment_id  = "slab ocean control experiment"
    glbAtt@source         = "CCSM3.0, version beta08 (2003): "+n1+\
                            "atmosphere: CAM3.0, T85L26;"+n1+\
                            "ocean     : N/A"+n1+\
                            "sea ice   : N/A"+n1+\
                            "land      : CLM3.0, T85"
  end if
  if (description.eq."SRESA1B") then
    glbAtt@experiment_id  = "720 ppm stabilization experiment (SRES A1B)"
    if (resolution.eq."T85_gx1") then
      forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
                "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
                "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
                "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
                "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
                "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
                "GHG loss rates   : noaamisc.r8.nc"+n1+\
                "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
                "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
                "SOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc"
    end if
    if (resolution.eq."T42_gx1") then
      forcing = "ozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
                "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
                "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
                "carbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc"+n1+\
                "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
                "GHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc"+n1+\
                "GHG loss rates   : noaamisc.r8.nc"+n1+\
                "DMS emissions    : DMS_emissions_64x128_c030722.nc"+n1+\
                "oxidants         : oxid_64x128_L26_clim_c050816.nc"+n1+\
                "SOx emissions    : SOx_emissions_A1_64x128_L2_1990-2100_c050912.nc"
    end if
  end if
  if (description.eq."SRESA1FI") then
    glbAtt@experiment_id  = "fossil-fuel intensive experiment (A1FI)"
    forcing = "ozone forcing    : A1FI.ozone.128x64_L18_1991-2100_c040528.nc"+n1+\
              "aerosol optics   : AerosolOptics_c040105.nc"+n1+\
              "aerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc"+n1+\
              "carbon scaling   : carbonscaling_A1FI_1990-2100_c040609.nc"+n1+\
              "solar forcing    : Fixed at 1366.5 W m-2"+n1+\
              "GHGs             : ghg_ipcc_A1FI_1870-2100_c040521.nc"+n1+\
              "GHG loss rates   : noaamisc.r8.nc"+n1+\
              "DMS emissions    : DMS_emissions_128x256_clim_c040122.nc"+n1+\
              "oxidants         : oxid_128x256_L26_clim_c040112.nc"+n1+\
              "SOx emissions    : SOx_emissions_A1FI_128x256_L2_1990-2100_c040608.nc"
  end if
  if (description.eq."SRESA2") then
    glbAtt@experiment_id  = "SRES A2 experiment"
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
  end if
  if (description.eq."SRESB1") then
    glbAtt@experiment_id  = "550 ppm stabilization experiment (SRES B1)"
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
  end if
;--------------------------------------------------------------------
; Check to make sure we got values!
;--------------------------------------------------------------------
  if (expid.eq."") then
    print("expid NOT DEFINED IN attributes_ccsm.ncl *PLEASE CORRECT*")
    status_exit(1)
  end if
  if (forcing.eq."") then
    print("forcing NOT DEFINED IN attributes_ccsm.ncl *PLEASE CORRECT*")
    status_exit(1)
  end if
;--------------------------------------------------------------------
; additional experiment-invariant global attributes
;--------------------------------------------------------------------
  glbAtt@history        = "Created from CCSM3 case "+case+n1+\
                          " by "+systemfunc("whoami")+n1+\
                          " on "+systemfunc("date")+n1+\
                          " "+n1+\
                          " For all data, added IPCC requested metadata"
  glbAtt@comment        = "This simulation was initiated from date "+branch+" of "+n1+\
                          " CCSM3 model run "+parent+" and executed on "+n1+\
                          " hardware "+machine+". The input external forcings are"+n1+\
                          forcing+n1+\
                           " Physical constants used for derived data:"+n1+\
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
