subroutine cmor_info
  ! Routine to initialize needed CMOR information
  use mycmor_info
  implicit none
  !
  mycmor%realization           = -99
  mycmor%initialization_method = -99
  mycmor%physics_version       = -99
  !
  mycmor%ack_NC = 'The CESM project is supported by the National Science Foundation and the Office of Science (BER) of the U.S. Department of Energy.\n'//&
       'NCAR is sponsored by the National Science Foundation.\n'//&
       'Computing resources were provided by the Climate Simulation Laboratory at the NCAR Computational and Information Systems Laboratory (CISL),\n'//&
       'sponsored by the National Science Foundation and other agencies.'
  mycmor%ack_NE = 'The CESM project is supported by the National Science Foundation and the Office of Science (BER) of the U.S. Department of Energy.\n'//&
       'NCAR is sponsored by the National Science Foundation.\n'//&
       'This research used resources of the National Energy Research Scientific Computing Center, which is supported by the Office of Science (BER)\n'//&
       'of the U.S. Department of Energy under Contract No. DE-AC02-05CH11231'
  mycmor%ack_OR = 'The CESM project is supported by the National Science Foundation and the Office of Science (BER) of the U.S. Department of Energy.\n'//&
       'NCAR is sponsored by the National Science Foundation.\n'//&
       'This research used resources of the Oak Ridge Leadership Computing Facility, located in the National Center for Computational Sciences\n'//&
       'at Oak Ridge National Laboratory, which is supported by the Office of Science (BER) of the Department of Energy under Contract DE-AC05-00OR22725.'
  !
  mycmor%forcing_note = 'Additional information on the external forcings used in this experiment can be found at\n'//&
       'http://www.cesm.ucar.edu/CMIP5/forcing_information'
  !
  ! Define arguments to 'cmor_dataset' - set by load_exp and init routines 
  !
  mycmor%model_id      = trim(exp(exp_found)%model_id)
  mycmor%outpath       = 'CMOR'
  mycmor%experiment_id = exp(exp_found)%expt_id(1:)
  mycmor%institution   = 'NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  mycmor%source        = trim(exp(exp_found)%model_id)
  mycmor%calendar      = 'noleap'
  mycmor%contact       = 'cesm_data@ucar.edu'
  mycmor%history       = ' '
  mycmor%comment       = ' '
  !
  ! References
  !
  select case (trim(adjustl(exp(exp_found)%model_id)))
  case ('CCSM4')
     mycmor%references = 'Gent P. R., et.al. 2011: The Community Climate System Model version 4. J. Climate, doi: 10.1175/2011JCLI4083.1'
  case ('CESM1')
     mycmor%references = 'TBD'
  case ('CCSM4-BGC')
     mycmor%references = 'TBD'
  case ('CCSM4-FSCHEM')
     mycmor%references = 'TBD'
  case ('CCSM4-WACCM')
     mycmor%references = 'TBD'
  case default
     write(*,*) 'Unknown model_id: ',trim(adjustl(exp(exp_found)%model_id)),' Stopping.'
     stop
  end select
  !leap_year    =
  !leap_month   =
  !month_lengths=
  mycmor%forcing       = exp(exp_found)%forcing(1:)
  mycmor%institute_id  = 'NCAR'
  !
end subroutine cmor_info
