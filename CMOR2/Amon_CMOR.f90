program Amon_CMOR
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
  use mycmor_info
  !
  implicit none
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  INTEGER::error_flag,var_ids
  REAL,DIMENSION(:,:),ALLOCATABLE::indat2a,indat2b,indat2c,cmordat
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::bnds_time
  DOUBLE PRECISION,DIMENSION(1)  ::tval
  DOUBLE PRECISION,DIMENSION(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,time_units,original_name
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,itab,ixw,ntimes
  integer::ilon,ilat,ipres,ilev,itim,itim2,ilon2,ilat2
  logical::all_continue
  !
  character(len=256),dimension(10)::ncfile
  real,dimension(10)::allmax,allmin,scale_factor
  integer,dimension(10)::ncid,var_found
  logical,dimension(10)::continue
  !
  ! GO!
  !
  mycmor%table_file = 'Tables/CMIP5_Amon'
  call load_table
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
  !
  read(*,*) case_read
  read(*,*) comp_read
  !
  ! Get experiment metadata from exp table and input case information
  !
  call get_exp_metadata
  !
  ! Get grid information
  !
  call get_atm_grid
  !
  ! Set up CMOR subroutine arguments
  !
  call get_cmor_info
  !
  ! Parse RIP code into components
  !
  call parse_rip
  !
  ! Step through CMOR table entries to see what CESM fields we can read and in process, and if so, do it!
  !
  table_loop: do itab = 1,num_tab
     xwalk_loop: do ixw = 1,num_xw
        mycmor%positive = ' '
        time_counter = 0
        var_counter  = 0
        error_flag   = 0
        var_found    = 0
        scale_factor = 1.
        allmax       = -1.e36
        allmin       =  1.e36
        all_continue = .true.
        continue(:)  = .false.
        time_units   = ' '
        original_name= ' '
!
! The meaty part
!
        if (xw(ixw)%entry == table(itab)%variable_entry) then
           do ivar = 1,xw(ixw)%ncesm_vars
              if ((trim(xw(ixw)%cesm_vars(ivar)) == 'UNKNOWN').or.(trim(xw(ixw)%cesm_vars(ivar)) == 'UNAVAILABLE')) then
                 write(*,'('' UNAVAILABLE/UNKNOWN: '',a,'' == '',a)') trim(xw(ixw)%entry),trim(table(itab)%variable_entry)
              else
                 write(ncfile(ivar),'(''data/'',a,''.'',a,''.'',a,''.'',a,''01-'',a,''12.nc'')') &
                      trim(case_read),&
                      trim(comp_read),&
                      trim(xw(ixw)%cesm_vars(ivar)),&
                      exp(exp_found)%begin_end(1:4),&
                      exp(exp_found)%begin_end(6:9)
                 inquire(file=trim(ncfile(ivar)),exist=continue(ivar))
                 if (.not.(continue(ivar))) then
                    write(ncfile(ivar),'(''data/'',a,''.'',a,''.'',a,''.'',a,''-01_cat_'',a,''-12.nc'')') &
                         trim(case_read),&
                         trim(comp_read),&
                         trim(xw(ixw)%cesm_vars(ivar)),&
                         exp(exp_found)%begin_end(1:4),&
                         exp(exp_found)%begin_end(6:9)
                    inquire(file=trim(ncfile(ivar)),exist=continue(ivar))
                 endif
                 write(*,'('' GOOD TO GO : '',a,'' == '',a,'' from CESM file: '',a)') &
                      trim(xw(ixw)%entry),&
                      trim(table(itab)%variable_entry),&
                      trim(ncfile(ivar))
                 if (.not.(continue(ivar))) write(*,*) trim(ncfile(ivar)),' NOT FOUND.'
              endif
              !
              ! Check and make sure all files available
              !
              all_continue = all_continue.and.(continue(ivar))
           enddo
           !
           ! Open CESM file(s) and get information(s)
           !
           if (all_continue) then
              do ivar = 1,xw(ixw)%ncesm_vars
                 call open_cdf(ncid(ivar),trim(ncfile(ivar)),.true.)
                 write(*,'(''OPENING: '',a80,'' ncid: '',i10)') trim(ncfile(ivar)),ncid(ivar)
                 call get_dims(ncid(ivar))
                 call get_vars(ncid(ivar))
                 !
                 do n=1,dim_counter
                    length = len_trim(dim_info(n)%name)
                    if(dim_info(n)%name(:length).eq.'time') then
                       ntimes = dim_info(n)%length
                    endif
                 enddo
                 call read_att_text(ncid(1),'time','units',time_units)
                 !
                 do n=1,var_counter
                    if (trim(var_info(n)%name) == trim(xw(ixw)%cesm_vars(ivar))) then
                       var_found(ivar) = n
                    endif
                 enddo
                 if (var_found(ivar) == 0) then
                    write(*,*) trim(xw(ixw)%cesm_vars(ivar)),' NEVER FOUND. STOP.'
                    stop
                 endif
                 !
                 allocate(time(ntimes),bnds_time(2,ntimes))
                 !
                 do n=1,ntimes
                    time_counter = n
                    call read_var(ncid(ivar),'time_bnds',bnds_time(:,n))
                    time(n) = (bnds_time(1,n)+bnds_time(2,n))/2.
                 enddo
              enddo
           endif
           if (all_continue) then
              !
              ! Specify path where tables can be found and indicate that existing netCDF files should be overwritten.
              !
              error_flag = cmor_setup(inpath='CMOR',netcdf_file_action=CMOR_REPLACE,logfile='LOG_CMOR.'//trim(xw(ixw)%entry))
              !
              error_flag = cmor_dataset(                              &
                   outpath=mycmor%outpath,                            &
                   experiment_id=mycmor%experiment_id,                &
                   institution=mycmor%institution,                    &
                   source=mycmor%source,                              &
                   calendar=mycmor%calendar,                          &
                   realization=mycmor%realization,                    &
                   contact=mycmor%contact,                            &
                   history=mycmor%history,                            &
                   comment=mycmor%comment,                            &
                   references=mycmor%references,                      &
                   model_id=mycmor%model_id,                          &
                   forcing=mycmor%forcing,                            &
                   initialization_method=mycmor%initialization_method,&
                   physics_version=mycmor%physics_version,            &
                   institute_id=mycmor%institute_id,                  &
                   parent_experiment_id=mycmor%parent_experiment_id,  &
                   parent_experiment_rip=mycmor%parent_experiment_rip,&
                   branch_time=mycmor%branch_time)
              if (error_flag < 0) then
                 write(*,*) 'Error on cmor_dataset!'
                 write(*,*) 'outpath=',mycmor%outpath
                 write(*,*) 'experiment_id=',mycmor%experiment_id
                 write(*,*) 'institution=',mycmor%institution
                 write(*,*) 'source=',mycmor%source
                 write(*,*) 'calendar=',mycmor%calendar
                 write(*,*) 'realization=',mycmor%realization
                 write(*,*) 'contact=',mycmor%contact
                 write(*,*) 'history=',mycmor%history
                 write(*,*) 'comment=',mycmor%comment
                 write(*,*) 'references=',mycmor%references
                 write(*,*) 'model_id=',mycmor%model_id
                 write(*,*) 'forcing=',mycmor%forcing
                 write(*,*) 'initialization_method=',mycmor%initialization_method
                 write(*,*) 'physics_version=',mycmor%physics_version
                 write(*,*) 'institute_id=',mycmor%institute_id
                 write(*,*) 'parent_experiment_id=',mycmor%parent_experiment_id
                 write(*,*) 'parent_experiment_rip=',mycmor%parent_experiment_rip
                 write(*,*) 'branch_time=',mycmor%branch_time
              endif
              !
              ! Add global metadata
              !
              call add_global_metadata
              !
              ! Define all axes that will be needed
              !
              ilat = cmor_axis(                  &
                   table=mycmor%table_file,  &
                   table_entry='latitude',       &
                   units='degrees_north',        &
                   length=SIZE(alats),           &
                   coord_vals=alats,             &
                   cell_bounds=bnds_lat)

              ilon = cmor_axis(  &
                   table=mycmor%table_file,  &
                   table_entry='longitude',      &
                   length=SIZE(alons),           &
                   units='degrees_east',         &
                   coord_vals=alons,             &
                   cell_bounds=bnds_lon)
              !
              ! Note that the time axis is defined next, but the time coordinate
              ! values and bounds will be passed to cmor through function
              ! cmor_write (later, below).
              !
              itim = cmor_axis(  &
                   table=mycmor%table_file,      &
                   table_entry='time',           &
                   units=time_units,             &
                   length=ntimes,                &
                   interval='30 days')

              write(*,*) 'CMOR axes defined'
              ! 
              ! Make manual alterations so that CMOR works. Silly code!
              !
              !
              if (xw(ixw)%ncesm_vars == 1) then
                 allocate(indat2a(nlons,nlats),cmordat(nlons,nlats))
                 write(original_name,'(a)') xw(ixw)%cesm_vars(1)
              endif
              if (xw(ixw)%ncesm_vars == 2) then
                 allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),cmordat(nlons,nlats))
                 write(original_name,'(a,'','',a)') (trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
              endif
              if (xw(ixw)%ncesm_vars == 3) then
                 allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats),cmordat(nlons,nlats))
                 write(original_name,'(a,'','',a,'','',a)') (trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
              endif
              !
              select case (xw(ixw)%entry)
              case ('tauu','tauv','hfss','rlut','rlutcs','hfls','rlus','rsus','rsuscs','rsut','rsutcs')
                 mycmor%positive = 'up'
              case ('rlds','rldscs','rsds','rsdscs','rsdt','rtmt')
                 mycmor%positive = 'down'
              case ('hurs','clt','ci')
                 var_info(var_found(1))%units = '1'
              case ('prc','pr','prsn')
                 var_info(var_found(1))%units = 'kg m-2 s-1'
              end select
                 !
              var_ids = cmor_variable(                     &
                   table=mycmor%table_file,                  &
                   table_entry=xw(ixw)%entry,               &
                   units=var_info(var_found(1))%units,              &
                   axis_ids=(/ ilon, ilat, itim /),              &
                   missing_value=var_info(var_found(1))%missing_value,&
                   positive=mycmor%positive,                       &
                   original_name=original_name)
              !
              write(*,*) 'cmor_variable:'
              write(*,*) 'varids=',var_ids
              write(*,*) 'table=',trim(mycmor%table_file)
              write(*,*) 'table_entry=',trim(xw(ixw)%entry)
              write(*,*) 'units=',trim(var_info(var_found(ivar))%units)
              write(*,*) 'missing_value=',var_info(var_found(ivar))%missing_value
              write(*,*) 'positive=',trim(mycmor%positive)
              write(*,*) 'original_name=',trim(original_name)
              !
              ! Cycle through time
              !
              time_loop: DO it=1, ntimes
                 time_counter = it
                 if (xw(ixw)%ncesm_vars == 1) then
                    call read_var(ncid(1),var_info(var_found(1))%name,indat2a)
!                    write(*,'(''Reading '',a20,'' T= '',i10)') trim(var_info(var_found(1))%name),it
                    allmax(1) = max(allmax(1),maxval(indat2a)) ; allmin(1) = min(allmin(1),minval(indat2a))
                    cmordat = indat2a
                 endif
                 if (xw(ixw)%ncesm_vars == 2) then
                    call read_var(ncid(1),var_info(var_found(1))%name,indat2a)
!                    write(*,'(''Reading '',a20,'' T= '',i10)') trim(var_info(var_found(1))%name),it
                    call read_var(ncid(2),var_info(var_found(2))%name,indat2b)
!                    write(*,'(''Reading '',a20,'' T= '',i10)') trim(var_info(var_found(2))%name),it
                    allmax(1) = max(allmax(1),maxval(indat2a)) ; allmin(1) = min(allmin(1),minval(indat2a))
                    allmax(2) = max(allmax(2),maxval(indat2b)) ; allmin(2) = min(allmin(2),minval(indat2b))
                    select case (xw(ixw)%entry)
                    case ('pr','prsn')
                       var_info(var_found(ivar))%units = 'kg m-2 s-1'
                       cmordat = 1000.*(indat2a + indat2b)
                    case ('rlus')
                       cmordat = indat2a + indat2b
                    case ('rsus','rsuscs','rsut','rsutcs','rtmt')
                       cmordat = indat2a - indat2b
                    case default
                       cmordat = indat2a
                    end select
                 endif
                 if (xw(ixw)%ncesm_vars == 3) then
                    call read_var(ncid(1),var_info(var_found(1))%name,indat2a)
                    write(*,'(''Reading '',a20,'' T= '',i10)') trim(var_info(var_found(1))%name),it
                    call read_var(ncid(2),var_info(var_found(2))%name,indat2b)
                    write(*,'(''Reading '',a20,'' T= '',i10)') trim(var_info(var_found(2))%name),it
                    call read_var(ncid(3),var_info(var_found(3))%name,indat2b)
                    write(*,'(''Reading '',a20,'' T= '',i10)') trim(var_info(var_found(3))%name),it
                    allmax(1) = max(allmax(1),maxval(indat2a)) ; allmin(1) = min(allmin(1),minval(indat2a))
                    allmax(2) = max(allmax(2),maxval(indat2b)) ; allmin(2) = min(allmin(2),minval(indat2b))
                    allmax(3) = max(allmax(3),maxval(indat2c)) ; allmin(3) = min(allmin(3),minval(indat2c))
                    cmordat = indat2a
                 endif
                 !
                 tval(1)   = time(it)
                 tbnd(1,1) = bnds_time(1,it)
                 tbnd(2,1) = bnds_time(2,it)
                 error_flag = cmor_write(      &
                      var_id        = var_ids, &
                      data          = cmordat, &
                      ntimes_passed = 1,       &
                      time_vals     = tval,    &
                      time_bnds     = tbnd)
                 if (error_flag < 0) then
                    write(*,*) 'Error writing ',xw(ixw)%entry, ', which I call ', xw(ixw)%cesm_vars
                    write(*,*) 'Processing time sample: ', time
                    stop
                 endif
              end do time_loop
              do ivar = 1,xw(ixw)%ncesm_vars
                 call close_cdf(ncid(ivar))
              enddo
              !
              ! Close all files opened by CMOR.
              !
              error_flag = cmor_close()
              write(*,'(''********************************************************************************'')')
              write(*,'(''********************************************************************************'')')
              write(*,'(''CMOR executed to completion; T#: '',i5,'' X#: '',i5,'' EXT: '',3(2g10.4))') &
                   itab,ixw,allmin(1:xw(ixw)%ncesm_vars),allmax(1:xw(ixw)%ncesm_vars)
              write(*,'(''********************************************************************************'')')
              write(*,'(''********************************************************************************'')')
              time_counter = 0
              var_counter  = 0
              error_flag   = 0
              var_found    = 0
              scale_factor = 1.
              allmax       = -1.e36
              allmin       =  1.e36
              continue(:)  = .false.
              mycmor%positive = ' '
              original_name= ' '
              !
              if (xw(ixw)%ncesm_vars == 1) deallocate(indat2a,cmordat)
              if (xw(ixw)%ncesm_vars == 2) deallocate(indat2a,indat2b,cmordat)
              if (xw(ixw)%ncesm_vars == 3) deallocate(indat2a,indat2b,indat2c,cmordat)
           endif
        endif
     enddo xwalk_loop
  enddo table_loop
  error_flag = cmor_close()
end program Amon_CMOR
