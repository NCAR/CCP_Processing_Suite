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
  use cmor_info
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
  ! Other variables
  !
  character(len=256)::time_units,exp_file,xwalk_file,ncfile
  character(len=256)::case_read,comp_read,svar,tstr
  integer::i,j,n,tcount,it,m,length,iexp,jexp,var_found
  integer::ncid,ntimes,nsamps_read
  integer::ilon,ilat,ipres,ilev,itim,itim2,ilon2,ilat2
  !
  allmax = -1.e36
  allmin =  1.e36
  !
  ! GO!
  !
  cmor%table_file = 'Tables/CMIP5_Amon'
  call load_table(cmor%table_file)
  !
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
  xwalk_file = 'xwalk_Amon.txt'
  call load_xwalk(xwalk_file)
  !
  ! Get experiment information
  !
  exp_file = 'experiments.txt'
  call load_exp(exp_file)
  write(*,'(''Number of experiments loaded: '',i5)') num_exp
  !
  read(*,*) case_read
  read(*,*) comp_read
  read(*,*) nsamps_read
  do i = 1,nsamps_read
     read(*,*) tstr
  enddo
  !  call parse_ncfile(ncfile,case,comp,svar,tstr)
  !
  ! Get experiment metadata from exp table and input case information
  !
  call get_exp_metadata(case_read,comp_read,nsamps_read)
  !
  ! Get grid information
  !
  call get_atm_grid
  !
  ! Parse RIP code into components
  !
  call parse_rip
  !
  ! Set up CMOR subroutine arguments
  !
  call get_cmor_info
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
  xw_found     = 0
  do j = 1,num_xw
     do n = 1,var_counter
        if ((trim(var_info(n)%name) == 'RHREFHT') .and. (trim(var_info(n)%units) == 'fraction')) var_info(n)%units = '%'
        if ((trim(var_info(n)%name) == 'FREQZM')  .and. (trim(var_info(n)%units) == 'fraction')) var_info(n)%units = '1'
        if (trim(var_info(n)%name) == trim(xw(j)%varin2d)) then
           var_found = n
           xw_found = j
        endif
     enddo
  enddo
  if ((xw_found /= 0).and.(var_found /= 0)) then
     write(*,*) 'MATCH field: ',trim(var_info(var_found)%name),' ',trim(var_info(var_found)%units),' ',trim(xw(xw_found)%varin2d),&
          ' ',var_info(var_found)%missing_value,' ',var_info(var_found)%FillValue,' ',cmor%positive
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
  error_flag = cmor_dataset(                            &
       outpath=cmor%outpath,                            &
       experiment_id=cmor%experiment_id,                &
       institution=cmor%institution,                    &
       source=cmor%source,                              &
       calendar=cmor%calendar,                          &
       realization=cmor%realization,                    &
       contact=cmor%contact,                            &
       history=cmor%history,                            &
       comment=cmor%comment,                            &
       references=cmor%references,                      &
       model_id=cmor%model_id,                          &
       forcing=cmor%forcing,                            &
       initialization_method=cmor%initialization_method,&
       physics_version=cmor%physics_version,            &
       institute_id=cmor%institute_id,                  &
       parent_experiment_id=cmor%parent_experiment_id,  &
       parent_experiment_rip=cmor%parent_experiment_rip,&
       branch_time=cmor%branch_time)
  !
  ! Add acknowledgements
  !
  write(*,*) exp(exp_found)%loc(1:2)
  if (exp(exp_found)%loc(1:2) == 'NC') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(cmor%ack_NC))
  if (exp(exp_found)%loc(1:2) == 'NE') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(cmor%ack_NE))
  if (exp(exp_found)%loc(1:2) == 'OR') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(cmor%ack_OR))
  !
  ! Add grid information
  !
  if (exp(exp_found)%grid(1:1) /= ' ') error_flag = cmor_set_cur_dataset_attribute("resolution",trim(adjustl(exp(exp_found)%grid)))
  !
  ! Add additional forcing information
  !
  if (cmor%forcing_note(1:1) /= ' ') error_flag = cmor_set_cur_dataset_attribute("forcing_note",trim(adjustl(cmor%forcing_note)))
  !
  !  Define all axes that will be needed
  !
  ilat = cmor_axis(                  &
       table=cmor%table_file,    &
       table_entry='latitude',       &
       units='degrees_north',        &
       length=SIZE(alats),           &
       coord_vals=alats,             &
       cell_bounds=bnds_lat)

  ilon = cmor_axis(  &
       table=cmor%table_file,    &
       table_entry='longitude',      &
       length=SIZE(alons),           &
       units='degrees_east',         &
       coord_vals=alons,             &
       cell_bounds=bnds_lon)

!  ipres = cmor_axis(  &
!       table=cmor%table_file,    &
!       table_entry='plevs',          &
!       units='Pa',                   &
!       length=SIZE(plevs),           &
!       coord_vals=plevs)

  !   note that the time axis is defined next, but the time coordinate
  !   values and bounds will be passed to cmor through function
  !   cmor_write (later, below).

  itim = cmor_axis(  &
       table=cmor%table_file,    &
       table_entry='time',           &
       units=time_units,             &
       length=ntimes,                &
       interval='30 days')

  !  define model eta levels (although these must be provided, they will
  !    actually be replaced by a+b before writing the netCDF file)

!!$  ilev = cmor_axis(  &
!!$       table=cmor%table_file,    &
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
!!$       table=cmor%table_file,  &
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
!!$          table=cmor%table_file,  &
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
          'table=',trim(cmor%table_file),'     ',                         &
          'table_entry=',trim(xw(xw_found)%entry2d),'     ',              &
          'units=',trim(var_info(var_found)%units),'     ',               &
          'missing_value=',var_info(var_found)%missing_value,'     ',     &
          'positive=',trim(cmor%positive),'     ',                        &
          'original_name=',xw(xw_found)%varin2d

     var2d_ids(m) = cmor_variable(                        &
          table=cmor%table_file,                          &
          table_entry=xw(xw_found)%entry2d,               &
          units=var_info(var_found)%units,                &
          axis_ids=(/ ilon, ilat, itim /),                &
          missing_value=var_info(var_found)%missing_value,&
          positive=cmor%positive,                         &
          original_name=xw(xw_found)%varin2d)
  ENDDO
  PRINT*, ' '
  PRINT*, 'completed everything up to writing output fields '
  PRINT*, ' '

  !  Loop through history files (each containing several different fields,
  !       but only a single month of data, averaged over the month).  Then
  !       extract fields of interest and write these to netCDF files (with
  !       one field per file, but all months included in the loop).
  !
  allocate(data2d(nlons,nlats))
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
