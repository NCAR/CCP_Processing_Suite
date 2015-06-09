subroutine get_cmor_info
  ! Routine to initialize needed CMOR information
  use mycmor_info
  use exp_info
  !
  implicit none
  !
!!$  mycmor%table_file(1:) = ' '
!!$  mycmor%outpath(1:) = ' '
!!$  mycmor%experiment_id(1:) = ' '
!!$  mycmor%institution(1:) = ' '
!!$  mycmor%source(1:) = ' '
!!$  mycmor%calendar(1:) = ' '
!!$  mycmor%contact(1:) = ' '
!!$  mycmor%history(1:) = ' '
!!$  mycmor%comment(1:) = ' '
!!$  mycmor%references(1:) = ' '
!!$  mycmor%model_id(1:) = ' '
!!$  mycmor%forcing(1:) = ' '
!!$  mycmor%institute_id(1:) = ' '
!!$  mycmor%parent_experiment_id(1:) = ' '
!!$  mycmor%parent_experiment_rip(1:) = ' '
!!$  mycmor%positive(1:) = ' '
!!$  mycmor%ack_NC(1:) = ' '
!!$  mycmor%ack_OR(1:) = ' '
!!$  mycmor%ack_NE(1:) = ' '
!!$  mycmor%forcing_note(1:) = ' '
!!$  mycmor%realization           = -99
!!$  mycmor%initialization_method = -99
!!$  mycmor%physics_version       = -99
!!$  mycmor%branch_time           = -99.d0
  !
  mycmor%ack_NC = 'The CESM project is supported by the National Science Foundation and the Office of Science (BER) of the U.S. Department of Energy. '//&
       'NCAR is sponsored by the National Science Foundation. '//&
       'Computing resources were provided by the Climate Simulation Laboratory at the NCAR Computational and Information Systems Laboratory (CISL), '//&
       'sponsored by the National Science Foundation and other agencies.'
  mycmor%ack_NE = 'The CESM project is supported by the National Science Foundation and the Office of Science (BER) of the U.S. Department of Energy. '//&
       'NCAR is sponsored by the National Science Foundation. '//&
       'This research used resources of the National Energy Research Scientific Computing Center, which is supported by the Office of Science (BER) '//&
       'of the U.S. Department of Energy under Contract No. DE-AC02-05CH11231'
  mycmor%ack_OR = 'The CESM project is supported by the National Science Foundation and the Office of Science (BER) of the U.S. Department of Energy. '//&
       'NCAR is sponsored by the National Science Foundation. '//&
       'This research used resources of the Oak Ridge Leadership Computing Facility, located in the National Center for Computational Sciences '//&
       'at Oak Ridge National Laboratory, which is supported by the Office of Science (BER) of the Department of Energy under Contract DE-AC05-00OR22725.'
  !
  mycmor%forcing_note = 'Additional information on the external forcings used in this experiment can be found at '//&
       'http://www.cesm.ucar.edu/CMIP5/forcing_information'
  !
  ! Define arguments to 'cmor_dataset' - set by load_exp and init routines 
  !
  mycmor%model_id      = trim(exp(exp_found)%model_id)
  mycmor%outpath       = 'CMOR'
  mycmor%experiment_id = exp(exp_found)%expt_id(1:)
  mycmor%source        = trim(exp(exp_found)%model_id)
  mycmor%calendar      = 'noleap'
  select case (exp(exp_found)%expt_id)
  case ('AEROCOM-A2-CTRL')
     mycmor%calendar   = 'gregorian'
  case default
     mycmor%calendar   = 'noleap'
  end select
  select case (exp(exp_found)%case)
  case ('f.e11.TSREFC1SD.f19.f19.ccmi23.001','cesm111ccmi23_geos5_htap_base','cesm111ccmi23_geos5_htap_eural')
     mycmor%calendar   = 'gregorian'
  case default
     mycmor%calendar   = 'noleap'
  end select

  mycmor%contact       = 'cesm_data@ucar.edu'
  mycmor%history       = ' '
  mycmor%comment       = ' '
  mycmor%positive      = ' '
  !
  ! References
  !
  select case (exp(exp_found)%model_id)
  case ('CCSM4')
     mycmor%references    = 'Gent P. R., et.al. 2011: The Community Climate System Model version 4. J. Climate, doi: 10.1175/2011JCLI4083.1'
     mycmor%institute_id  = 'NCAR'
     mycmor%institution   = 'NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-CAM5')
     mycmor%references    = 'Neale, R., et.al. 2012: Coupled simulations from CESM1 using the Community Atmosphere Model version 5: (CAM5).'//&
          ' See also http://www.cesm.ucar.edu/publications'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-BGC')
     mycmor%references    = 'Lindsay K., et al.: Preindustrial Control and 20th Century Experiments with the Earth System Model CESM1-(BGC) (in preparation for Journal of Climate).\n'//&
          ' See also http://www.cesm.ucar.edu/publications'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-FASTCHEM')
     mycmor%references    = 'Cameron-Smith, P., et.al. 2012: CESM1 with superfast chemistry.'//&
          ' See also http://www.cesm.ucar.edu/publications'
     mycmor%references    = 'TBD'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-WACCM')
     mycmor%references    = 'Marsh, D., et.al. 2012: WACCM4 simulations of atmospheric trends from 1850 to present.'//&
          ' See also http://www.cesm.ucar.edu/publications'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-CAM4Chem')
     mycmor%references    = 'Tilmes S. in preparation'//&
          ' See also http://www.cesm.ucar.edu/publications'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-CAM4ChemSD')
     mycmor%references    = 'Tilmes S. in preparation'//&
          ' See also http://www.cesm.ucar.edu/publications'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-BASE')
     mycmor%references    = 'Tilmes S. in preparation'//&
          ' See also http://www.cesm.ucar.edu/publications'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-EASALL')
     mycmor%references    = 'Tilmes S. in preparation'//&
          ' See also http://www.cesm.ucar.edu/publications'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
  case ('CESM1-EURALL')
     mycmor%references    = 'Tilmes S. in preparation'//&
          ' See also http://www.cesm.ucar.edu/publications'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-CAM5.1-FV2')
     mycmor%references    = 'Neale, R. B. et al. (2011a) Description of the NCAR Community Atmosphere Model (CAM5), Technical Report NCAR/TN-486+STR.\n'//&
          'Gent P. R., et.al. 2011: The Community Climate System Model version 4. J. Climate, doi: 10.1175/2011JCLI4083.1'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'PNNL (Pacific Northwest National Laboratory) Richland, WA, USA/NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case default
     write(*,*) 'Unknown model_id: ',trim(adjustl(exp(exp_found)%model_id)),' Stopping.'
     stop
  end select
  !
  ! Global comment
  !
  select case (exp(exp_found)%expt_id)
  case ('tamip200810','tamip200901','tamip200904','tamip200907')
     mycmor%comment = 'Land, aerosols, and the atmospheric fast processes were spun up as follows: \n'//&
          'A 5-year (2003-01-01 to 2008-05-01) free-running (no nudging) AMIP run was made to bring the \n'//&
          'land/atmosphere state to reasonable equilibrium (using monthly mean AMIP SSTs); \n'//&
          'Beginning from 2008-05-01 and running through the entire Transpose-AMIP period, the simulation \n'//&
          'was then nudged using 6hr YOTC analyses in the CAPT protocol (described in Phillips, et al, 2004: "Evaluating \n'//&
          'Parameterizations in General Circulation Models: Climate Simulation Meets Weather Prediction", BAMS, 85, 1903-1915). \n'//&
          'Daily YOTC SSTs and sea-ice cover data were used.'
  case default
     mycmor%comment = 'CESM home page: http://www.cesm.ucar.edu'
  end select
  !leap_year    =
  !leap_month   =
  !month_lengths=
  mycmor%forcing       = exp(exp_found)%forcing(1:)
  !
end subroutine get_cmor_info
