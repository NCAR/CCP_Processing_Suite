PROGRAM Amon_CMOR
  ! Convert CCSM4 atmospheric monthly (cam2.h0) data from single-field format
  ! to CMOR-compliant format
  !
  ! NOTE: 'model_id' and first part of 'source' MUST MATCH or CMOR will throw error
  !
  use cmor_users_functions
  use counters_netcdf_jfl
  use interfaces_netcdf_jfl
  use definitions_netcdf_jfl
  use exp_info
  use table_info
  use xwalk_info
  use atm_grid_info
  !
  implicit none
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  INTEGER::error_flag
  INTEGER::znondim_id,zfactor_id
  INTEGER,DIMENSION(1000)::var2d_ids
  REAL,DIMENSION(:,:),ALLOCATABLE::data2d
  REAL,DIMENSION(:,:,:),ALLOCATABLE::data3d
  real::allmax,allmin
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::bnds_time
  DOUBLE PRECISION,DIMENSION(1)  ::tval
  DOUBLE PRECISION,DIMENSION(2,1)::tbnd
  !
  ! CMOR-needed variables
  !
  CHARACTER(len=256)::table_file,outpath,experiment_id,institution,source,calendar,contact,history
  CHARACTER(len=256)::comment,references,model_id,forcing,institute_id
  CHARACTER(len=256)::parent_experiment_id,parent_experiment_rip,positive
  INTEGER::realization,initialization_method,physics_version
  DOUBLE PRECISION::branch_time
  !
  ! other variables:
  !
  character(len=256)::time_units,exp_file,xwalk_file,ncfile
  character(len=256)::case,comp,svar,tstr
  character(len=512)::ack_NC,ack_OR,ack_NE,forcing_note
  integer::i,j,n,tcount,it,m,length,iexp,jexp,exp_found,parent_found,xw_found,var_found,tab_found
  integer::ncid,ntimes
  integer::ilon,ilat,ipres,ilev,itim,itim2,ilon2,ilat2
  !
  allmax = -1.e36
  allmin =  1.e36
  realization           = -99
  initialization_method = -99
  physics_version       = -99
  !
  ack_NC = 'The CESM project is supported by the National Science Foundation and the Office of Science (BER) of the U.S. Department of Energy.\n'//&
       'NCAR is sponsored by the National Science Foundation.\n'//&
       'Computing resources were provided by the Climate Simulation Laboratory at the NCAR Computational and Information Systems Laboratory (CISL),\n'//&
       'sponsored by the National Science Foundation and other agencies.'
  ack_NE = 'The CESM project is supported by the National Science Foundation and the Office of Science (BER) of the U.S. Department of Energy.\n'//&
       'NCAR is sponsored by the National Science Foundation.\n'//&
       'This research used resources of the National Energy Research Scientific Computing Center, which is supported by the Office of Science (BER)\n'//&
       'of the U.S. Department of Energy under Contract No. DE-AC02-05CH11231'
  ack_OR = 'The CESM project is supported by the National Science Foundation and the Office of Science (BER) of the U.S. Department of Energy.\n'//&
       'NCAR is sponsored by the National Science Foundation.\n'//&
       'This research used resources of the Oak Ridge Leadership Computing Facility, located in the National Center for Computational Sciences\n'//&
       'at Oak Ridge National Laboratory, which is supported by the Office of Science (BER) of the Department of Energy under Contract DE-AC05-00OR22725.'
  !
  forcing_note = 'Additional information on the external forcings used in this experiment can be found at\n'//&
       'http://www.cesm.ucar.edu/CMIP5/forcing_information'
  !
  write(*,*) len_trim(ack_NC),' ',len_trim(ack_NE),' ',len_trim(ack_OR)
  !
  ! ================================
  !  Execution begins here:
  ! ================================
  !
  read(*,*) ncfile
  call parse_ncfile(ncfile,case,comp,svar,tstr)
  !
  table_file = 'Tables/CMIP5_Amon'
  call load_table(table_file)
  !
  ! Get experiment information
  !
  exp_file = 'experiments.txt'
  call load_exp(exp_file)
  write(*,'(''Number of experiments loaded: '',i5)') num_exp
  !
  exp_found = 0
  do iexp = 1,num_exp
     if (trim(exp(iexp)%case) == trim(case)) then
        exp_found = iexp
     endif
  enddo
  if (exp_found /= 0) then
     write(*,*) 'MATCH experiment: ',trim(exp(exp_found)%case),' ',trim(case),' ',trim(exp(exp_found)%model_id)
     write(*,*) trim(exp(exp_found)%run_refcase),' ',trim(exp(exp_found)%repotag),' ',trim(exp(exp_found)%compset),' ',trim(exp(exp_found)%rip_code)
  else
     write(*,*) 'case ',case(1:len_trim(case)),' NOT FOUND. Dying.'
     stop
  endif
  !
  parent_found = 0
  do iexp = 1,num_exp
     if (trim(exp(exp_found)%run_refcase) == trim(exp(iexp)%case)) then
        parent_found = iexp
     endif
  enddo
  if (parent_found /= 0) then
     write(*,*) 'MATCH parent: ',trim(case),' ',trim(exp(parent_found)%case),' ',trim(exp(parent_found)%rip_code)
     parent_experiment_id  = exp(parent_found)%expt_id(1:)
     parent_experiment_rip = exp(parent_found)%rip_code(1:)
     read(exp(exp_found)%run_refdate(1:4),'(f4.0)') branch_time
  else
     if (trim(exp(exp_found)%run_refcase) == "N/A") then
        parent_experiment_id  = "N/A"
        parent_experiment_rip = "N/A"
        branch_time = 0.
     else
        write(*,*) 'parent ',case(1:len_trim(case)),' NOT FOUND. Dying.'
        stop
     endif
  endif
  !
  ! Get grid information
  !
  call get_atm_grid
  allocate(data2d(nlons,nlats))
  !
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
  xwalk_file = 'xwalk_Amon.txt'
  call load_xwalk(xwalk_file)
  !
  ! Parse RIP code into components
  !
  write(*,'(''RIP 0: '',a,'' '',3i3)') trim(exp(exp_found)%rip_code),realization,initialization_method,physics_version
  call parse_rip(exp(exp_found)%rip_code,realization,initialization_method,physics_version)
  write(*,'(''RIP 1: '',a,'' '',3i3)') trim(exp(exp_found)%rip_code),realization,initialization_method,physics_version
  !
  ! Open actual data file and get information
  !
  call open_cdf(ncid,'data/'//ncfile,.true.)
  call get_dims(ncid)
  call get_vars(ncid)
  !
  do n=1,dim_counter
     length = len_trim(dim_info(n)%name)
     if(dim_info(n)%name(:length).eq.'time') then
        ntimes = dim_info(n)%length
     endif
  enddo
  time_units(1:) = ' '
  call read_att_text(ncid,'time','units',time_units)
  !
  ALLOCATE(time(ntimes),bnds_time(2,ntimes))
  !
  do n=1,ntimes
     time_counter = n
     call read_var(ncid,'time_bnds',bnds_time(:,n))
     time(n) = (bnds_time(1,n)+bnds_time(2,n))/2.
  enddo
  !
  var_found    = 0
   xw_found    = 0
  positive(1:) = ' '
  do n = 1,var_counter
     select case (trim(adjustl(var_info(n)%name)))
     case ('FSDTOA','FSN200','FSN200C','FSNIRTOA','FSNRTOAC','FSNRTOAS','FLUT','FLUTC','FSNS','FSNSC','FSNT','FSNTC','FSNTOA','FSNTOAC','FSUTOA')
        positive = 'up'
        write(*,*) trim(adjustl(var_info(n)%name)),' up'
     case ('FSDS','FSDSC','FDL','FDLC','FLDS','FLDSC','FLN200','FLN200C','FLNS','FLNSC','FLNT','FLNTC','FUL','FULC')
        positive = 'down'
        write(*,*) trim(adjustl(var_info(n)%name)),' down'
     case default
     end select
     !
     do j = 1,num_xw
        if (trim(var_info(n)%name) == trim(xw(j)%varin2d)) then
           var_found = n
            xw_found = j
            if ((trim(var_info(n)%name) == 'RHREFHT') .and. (trim(var_info(n)%units) == 'fraction')) var_info(n)%units = '%'
            if ((trim(var_info(n)%name) == 'FREQZM')  .and. (trim(var_info(n)%units) == 'fraction')) var_info(n)%units = '1'
        endif
     enddo
  enddo
  if ((xw_found /= 0).and.(var_found /= 0)) then
     write(*,*) 'MATCH field: ',trim(var_info(var_found)%name),' ',trim(var_info(var_found)%units),' ',trim(xw(xw_found)%varin2d),&
          ' ',var_info(var_found)%missing_value,' ',var_info(var_found)%FillValue,' ',positive
  endif
  !
  tab_found = 0
  do n = 1,num_tab
     if (trim(adjustl(table(n)%out_name)) == trim(adjustl(xw(xw_found)%entry2d))) then
        tab_found = n
     endif
  enddo
  if (tab_found /= 0) then
     write(*,*) 'MATCH table # ',trim(table(tab_found)%out_name),' ',trim(table(tab_found)%units),' ',trim(xw(xw_found)%entry2d)
  endif
  !
  ! Specify path where tables can be found and indicate that existing netCDF files should be overwritten.
  !
  error_flag = cmor_setup(inpath='CMOR', netcdf_file_action='replace')
  !
  ! Define arguments to 'cmor_dataset'
  ! Set by load_exp and init routines above
  !
  model_id      = trim(exp(exp_found)%model_id)
  outpath       = 'CMOR'
  experiment_id = exp(exp_found)%expt_id(1:)
  institution   = 'NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  source        = trim(exp(exp_found)%model_id)//' (repository tag: '//trim(exp(exp_found)%repotag)//' compset: '//trim(exp(exp_found)%compset)//')'
  calendar      = 'noleap'
  contact       = 'cesm_data@ucar.edu'
  history       = ' '
  comment       = ' '
  !
  ! References
  !
  select case (trim(adjustl(exp(exp_found)%model_id)))
  case ('CCSM4')
     references = 'Gent P. R., et.al. 2011: The Community Climate System Model version 4. J. Climate, doi: 10.1175/2011JCLI4083.1'
  case ('CESM1')
     references = 'TBD'
  case ('CCSM4-BGC')
     references = 'TBD'
  case ('CCSM4-FSCHEM')
     references = 'TBD'
  case ('CCSM4-WACCM')
     references = 'TBD'
  case default
     write(*,*) 'Unknown model_id: ',trim(adjustl(exp(exp_found)%model_id)),' Stopping.'
     stop
  end select
  !leap_year    =
  !leap_month   =
  !month_lengths=
  forcing       = exp(exp_found)%forcing(1:)
  institute_id  = 'NCAR'
  !
  !
  error_flag = cmor_dataset(                       &
       outpath=outpath,                            &
       experiment_id=experiment_id,                &
       institution=institution,                    &
       source=source,                              &
       calendar=calendar,                          &
       realization=realization,                    &
       contact=contact,                            &
       history=history,                            &
       comment=comment,                            &
       references=references,                      &
       model_id=model_id,                          &
       forcing=forcing,                            &
       initialization_method=initialization_method,&
       physics_version=physics_version,            &
       institute_id=institute_id,                  &
       parent_experiment_id=parent_experiment_id,  &
       parent_experiment_rip=parent_experiment_rip,&
       branch_time=branch_time)
  !
  ! Add acknowledgements
  !
  write(*,*) exp(exp_found)%loc(1:2)
  if (exp(exp_found)%loc(1:2) == 'NC') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(ack_NC))
  if (exp(exp_found)%loc(1:2) == 'NE') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(ack_NE))
  if (exp(exp_found)%loc(1:2) == 'OR') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(ack_OR))
  !
  ! Add grid information
  !
  if (exp(exp_found)%grid(1:1) /= ' ') error_flag = cmor_set_cur_dataset_attribute("resolution",trim(adjustl(exp(exp_found)%grid)))
  !
  ! Add additional forcing information
  !
  if (forcing_note(1:1) /= ' ') error_flag = cmor_set_cur_dataset_attribute("forcing_note",trim(adjustl(forcing_note)))
  !
  !  Define all axes that will be needed
  !
  ilat = cmor_axis(                  &
       table=table_file,    &
       table_entry='latitude',       &
       units='degrees_north',        &
       length=SIZE(alats),           &
       coord_vals=alats,             &
       cell_bounds=bnds_lat)

  ilon = cmor_axis(  &
       table=table_file,    &
       table_entry='longitude',      &
       length=SIZE(alons),           &
       units='degrees_east',         &
       coord_vals=alons,             &
       cell_bounds=bnds_lon)

!  ipres = cmor_axis(  &
!       table=table_file,    &
!       table_entry='plevs',          &
!       units='Pa',                   &
!       length=SIZE(plevs),           &
!       coord_vals=plevs)

  !   note that the time axis is defined next, but the time coordinate
  !   values and bounds will be passed to cmor through function
  !   cmor_write (later, below).

  itim = cmor_axis(  &
       table=table_file,    &
       table_entry='time',           &
       units=time_units,             &
       length=ntimes,                &
       interval='30 days')

  !  define model eta levels (although these must be provided, they will
  !    actually be replaced by a+b before writing the netCDF file)

!!$  ilev = cmor_axis(  &
!!$       table=table_file,    &
!!$       table_entry='standard_hybrid_sigma',       &
!!$       units='1', &
!!$       length=SIZE(zlevs),           &
!!$       coord_vals=zlevs,             &
!!$       cell_bounds=zlev_bnds)
!!$
!!$  !   define z-factors needed to transform from model level to pressure
!!$  call read_var(ncid,'P0',p0)
!!$  p0 = p0 * 100
!!$  call read_var(ncid,'hyam',a_coeff)
!!$  call read_var(ncid,'hybm',b_coeff)
!!$  call read_var(ncid,'hyai',a_coeff_bnds)
!!$  call read_var(ncid,'hybi',b_coeff_bnds)
!!$
!!$  error_flag = cmor_zfactor(  &
!!$       zaxis_id=ilev,                      &
!!$       zfactor_name='p0',                  &
!!$       units='Pa',                         &
!!$       zfactor_values = p0)
!!$
!!$  error_flag = cmor_zfactor(  &
!!$       zaxis_id=ilev,                       &
!!$       zfactor_name='b',                    &
!!$       axis_ids= (/ ilev /),                &
!!$       zfactor_values = b_coeff,            &
!!$       zfactor_bounds = b_coeff_bnds  )
!!$
!!$  error_flag = cmor_zfactor(  &
!!$       zaxis_id=ilev,                       &
!!$       zfactor_name='a',                    &
!!$       axis_ids= (/ ilev /),                &
!!$       zfactor_values = a_coeff,            &
!!$       zfactor_bounds = a_coeff_bnds )
!!$
!!$  zfactor_id = cmor_zfactor(  &
!!$       zaxis_id=ilev,                         &
!!$       zfactor_name='ps',                     &
!!$       axis_ids=(/ ilon, ilat, itim /),       &
!!$       units='Pa' )
!!$  !  Define the only field to be written that is a function of model level
!!$  !    (appearing in IPCC table A1c)
!!$
!!$  var3d_ids(1) = cmor_variable(    &
!!$       table=table_file,  &
!!$       table_entry=entry3d(1),     &
!!$       units=units3d(1),           &
!!$       axis_ids=(/ ilon, ilat, ilev, itim /),  &
!!$       missing_value=1.0e28, &
!!$       original_name=varin3d(1))
!!$

!!$  !  Define variables appearing in IPCC table A1c that are a function of pressure
!!$  !         (3-d variables)
!!$
!!$  DO m=2,n3d
!!$     var3d_ids(m) = cmor_variable(    &
!!$          table=table_file,  &
!!$          table_entry=entry3d(m),     &
!!$          units=units3d(m),           &
!!$          axis_ids=(/ ilon, ilat, ipres, itim /), &
!!$          missing_value=1.0e28,       &
!!$          original_name=varin3d(m))
!!$  ENDDO
!!$
!!$
!!$  !  Define variables appearing in IPCC table A1a (2-d variables)
  DO m=1,1
     write(*,*) &
          'table=',trim(table_file),'     ',                              &
          'table_entry=',trim(xw(xw_found)%entry2d),'     ',              &
          'units=',trim(var_info(var_found)%units),'     ',               &
          'missing_value=',var_info(var_found)%missing_value,'     ',     &
          'positive=',trim(positive),'     ',                             &
          'original_name=',xw(xw_found)%varin2d

     var2d_ids(m) = cmor_variable(                        &
          table=table_file,                               &
          table_entry=xw(xw_found)%entry2d,               &
          units=var_info(var_found)%units,                &
          axis_ids=(/ ilon, ilat, itim /),                &
          missing_value=var_info(var_found)%missing_value,&
          positive=positive,                              &
          original_name=xw(xw_found)%varin2d)
  ENDDO
  PRINT*, ' '
  PRINT*, 'completed everything up to writing output fields '
  PRINT*, ' '

  !  Loop through history files (each containing several different fields,
  !       but only a single month of data, averaged over the month).  Then
  !       extract fields of interest and write these to netCDF files (with
  !       one field per file, but all months included in the loop).

  time_loop: DO it=1, ntimes

     ! Cycle through the 2-d fields, retrieve the requested variable and
     ! append each to the appropriate netCDF file.

     ! The user must write the code that fills the arrays of data
     ! that will be passed to CMOR.  The following line is simply a
     ! a place-holder for the user's code, which should replace it.
     !
     ! append a single time sample of data for a single field to
     ! the appropriate netCDF file.

     time_counter = it
     call read_var(ncid,var_info(var_found)%name,data2d)
     tval(1)   = time(it)
     tbnd(1,1) = bnds_time(1,it)
     tbnd(2,1) = bnds_time(2,it)
     m = 1
     error_flag = cmor_write(                                  &
          var_id        = var2d_ids(m),                        &
          data          = data2d,                              &
          ntimes_passed = 1,                                   &
          time_vals     = tval,                                &
          time_bnds     = tbnd)
!     write(*,'(i6,2g12.4)') it,minval(data2d),maxval(data2d)
     allmax = max(allmax,maxval(data2d))
     allmin = min(allmin,minval(data2d))

     IF (error_flag < 0) THEN
        ! write diagnostic messages to standard output device
        write(*,*) ' Error encountered writing IPCC Table A1a field ',xw(xw_found)%entry2d, ', which I call ', xw(xw_found)%varin2d
        write(*,*) ' Was processing time sample: ', time
     END IF

  END DO time_loop
  !
  ! Close all files opened by CMOR.
  !
  error_flag = cmor_close()
  write(*,*) '********************************************************************************'
  write(*,*) 'code executed to completion '
  write(*,*) allmin,allmax
  write(*,*) '********************************************************************************'

END PROGRAM Amon_CMOR
