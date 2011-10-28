program OImon_CMOR
  ! Convert CCSM4 sea ice monthly (cice.h) data from single-field format
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
  use grid_info
  use mycmor_info
  !
  implicit none
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  integer::error_flag,var_ids
  real,dimension(:,:),allocatable::indat2anh,indat2ash,indat2bnh,indat2bsh,cmordat
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::time_bnds
  double precision,dimension(1)  ::tval
  double precision,dimension(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,itab,ixw
  integer,dimension(:),allocatable::i_indices,j_indices
  logical::all_continue
  !
  character(len=256),dimension(10)::ncfilenh,ncfilesh
  real,dimension(10)::allmax,allmin,scale_factor
  integer,dimension(10)::ncidnh,ncidsh,var_found
  logical,dimension(10)::continue
  !
  ! GO!
  !
  mycmor%table_file = 'Tables/CMIP5_OImon'
  call load_table_info
  !
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
  xwalk_file = 'xwalk_OImon.txt'
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
  call get_ice_grid
  allocate(i_indices(nlons),j_indices(nlats))
  do i = 1,nlons
     i_indices(i) = i
  enddo
  do j = 1,nlats
     j_indices(j) = j
  enddo
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
                 write(ncfilenh(ivar),'(''data/'',a,''.'',a,''.'',a,''_nh.'',a,''01-'',a,''12.nc'')') &
                      trim(case_read),&
                      trim(comp_read),&
                      trim(xw(ixw)%cesm_vars(ivar)),&
                      exp(exp_found)%begin_end(1:4),&
                      exp(exp_found)%begin_end(6:9)
                 inquire(file=trim(ncfilenh(ivar)),exist=continue(ivar))
                 write(ncfilesh(ivar),'(''data/'',a,''.'',a,''.'',a,''_sh.'',a,''01-'',a,''12.nc'')') &
                      trim(case_read),&
                      trim(comp_read),&
                      trim(xw(ixw)%cesm_vars(ivar)),&
                      exp(exp_found)%begin_end(1:4),&
                      exp(exp_found)%begin_end(6:9)
                 inquire(file=trim(ncfilenh(ivar)),exist=continue(ivar))
                 if (.not.(continue(ivar))) then
                    write(ncfilenh(ivar),'(''data/'',a,''.'',a,''.'',a,''_nh.'',a,''-01_cat_'',a,''-12.nc'')') &
                         trim(case_read),&
                         trim(comp_read),&
                         trim(xw(ixw)%cesm_vars(ivar)),&
                         exp(exp_found)%begin_end(1:4),&
                         exp(exp_found)%begin_end(6:9)
                    inquire(file=trim(ncfilenh(ivar)),exist=continue(ivar))
                    write(ncfilesh(ivar),'(''data/'',a,''.'',a,''.'',a,''_sh.'',a,''-01_cat_'',a,''-12.nc'')') &
                         trim(case_read),&
                         trim(comp_read),&
                         trim(xw(ixw)%cesm_vars(ivar)),&
                         exp(exp_found)%begin_end(1:4),&
                         exp(exp_found)%begin_end(6:9)
                    inquire(file=trim(ncfilesh(ivar)),exist=continue(ivar))
                 endif
                 if (.not.(continue(ivar))) then
                    write(*,'(a,'' and/or '',a,'' NOT FOUND.'')') trim(ncfilenh(ivar)),trim(ncfilesh(ivar))
                 else
                    write(*,'('' GOOD TO GO : '',a,'' == '',a,'' from CESM files: '',a,'' and '',a)') &
                         trim(xw(ixw)%entry),&
                         trim(table(itab)%variable_entry),&
                         trim(ncfilenh(ivar)),&
                         trim(ncfilesh(ivar))
                 endif
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
                 !
                 ! Get NH data and time values
                 !
                 call open_cdf(ncidnh(ivar),trim(ncfilenh(ivar)),.true.)
                 write(*,'(''OPENING: '',a80,'' ncid: '',i10)') trim(ncfilenh(ivar)),ncidnh(ivar)
                 call get_dims(ncidnh(ivar))
                 call get_vars(ncidnh(ivar))
                 !
                 do n=1,dim_counter
                    length = len_trim(dim_info(n)%name)
                    if(dim_info(n)%name(:length).eq.'time') then
                       ntimes = dim_info(n)%length
                    endif
                 enddo
                 call read_att_text(ncidnh(1),'time','units',time_units)
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
                 if (.not.(allocated(time)))      then
                    allocate(time(ntimes))
                    write(*,*) 'allocate(time(ntimes))'
                 endif
                 if (.not.(allocated(time_bnds))) then
                    allocate(time_bnds(2,ntimes))
                    write(*,*) 'allocate(time_bnds(2,ntimes))'
                 endif
                 !
                 do n=1,ntimes
                    time_counter = n
                    call read_var(ncidnh(ivar),'time_bounds',time_bnds(:,n))
                    time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
                 enddo
                 !
                 ! Get SH data info
                 !
                 call open_cdf(ncidsh(ivar),trim(ncfilesh(ivar)),.true.)
                 write(*,'(''OPENING: '',a80,'' ncid: '',i10)') trim(ncfilesh(ivar)),ncidsh(ivar)
                 call get_dims(ncidsh(ivar))
                 call get_vars(ncidsh(ivar))
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
              enddo
           endif
           if (all_continue) then
              !
              ! Specify path where tables can be found and indicate that existing netCDF files should be overwritten.
              !
              write(logfile,'(''log_cmor.'',a,''.'',a,''_'',a)') &
                   trim(mycmor%experiment_id),&
                   trim(exp(exp_found)%rip_code),&
                   trim(xw(ixw)%entry)
              error_flag = cmor_setup(inpath='CMOR',&
                   netcdf_file_action=CMOR_REPLACE,&
                   logfile=logfile)
              table_ids(1) = cmor_load_table(mycmor%table_file)
              call cmor_set_table(table_ids(1))
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
              ! Define axes via 'cmor_axis'
              !
              table_ids(2) = cmor_load_table('Tables/CMIP5_grids')
              call cmor_set_table(table_ids(2))
              call define_ice_axes(table(itab)%dimensions)
              call cmor_set_table(table_ids(1))
              !
              ! Make manual alterations so that CMOR works. Silly code!
              !
              allocate(indat2anh(nlons,104),indat2ash(nlons,76),cmordat(nlons,384))
              if (xw(ixw)%ncesm_vars == 1) then
                 write(original_name,'(a)') xw(ixw)%cesm_vars(1)
              endif
              if (xw(ixw)%ncesm_vars .ge. 2) then
                 allocate(indat2bnh(nlons,104),indat2bsh(nlons,76))
                 write(original_name,'(a,'','',a)') (trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
              endif
              !
              ! Define unit renamings and/or "positive" as needed
              !
              select case (xw(ixw)%entry)
              case ('evap')
                 mycmor%positive = 'up'
                 var_info(var_found(1))%units = 'kg m-2 s-1'
              case ('bmelt','grCongel','grFrazil','pr','prsn','snomelt','snoToIce','tmelt')
                 var_info(var_found(1))%units = 'kg m-2 s-1'
              end select
              !
              write(*,*) 'cmor_variable:'
              write(*,*) 'varids=',var_ids
              write(*,*) 'table=',trim(mycmor%table_file)
              write(*,*) 'table_entry=',trim(xw(ixw)%entry)
              write(*,*) 'dimensions=',trim(table(itab)%dimensions)
              write(*,*) 'axis_ids=',axis_ids
              write(*,*) 'grid_id=',grid_id
              write(*,*) 'units=',trim(var_info(var_found(1))%units)
              write(*,*) 'missing_value=',var_info(var_found(1))%missing_value
              write(*,*) 'positive=',trim(mycmor%positive)
              write(*,*) 'original_name=',trim(original_name)
              !
              var_ids = cmor_variable(                                &
                   table=mycmor%table_file,                           &
                   table_entry=xw(ixw)%entry,                         &
                   units=var_info(var_found(1))%units,                &
                   axis_ids=(/grid_id(1),axis_ids(3)/),               &
                   missing_value=var_info(var_found(1))%missing_value,&
                   positive=mycmor%positive,                          &
                   original_name=original_name,                       &
                   comment=xw(ixw)%comment)
              !
              ! Cycle through time
              !
              time_loop: DO it=1, ntimes
                 time_counter = it
                 if (xw(ixw)%ncesm_vars == 1) then
                    call read_var(ncidnh(1),var_info(var_found(1))%name,indat2anh)
                    call read_var(ncidsh(1),var_info(var_found(1))%name,indat2ash)
                    write(*,'(''Reading NH and SH '',a20,'' T= '',i10)') trim(var_info(var_found(1))%name),it
                    cmordat(:,  1: 76) = indat2ash(:,1: 76)
                    cmordat(:,281:384) = indat2anh(:,1:104)
                 endif
                 !
                 ! Perform necessary derivations
                 !
                 select case (xw(ixw)%entry)
                 case ('bmelt','evap','grCongel','grFrazil','pr','prsn','snomelt','snoToIce','tmelt')
                    ! Convert cm day-1 to kg m-2 s-1
                    cmordat = cmordat * (1./86400.) * (1./100.) * 1000.
                 end select
                 !
                 ! Write 'em, Dano!
                 !
                 tval(1)   = time(it)
                 tbnd(1,1) = time_bnds(1,it)
                 tbnd(2,1) = time_bnds(2,it)
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
                 call close_cdf(ncidnh(ivar))
                 call close_cdf(ncidsh(ivar))
              enddo
              !
              ! Close all files opened by CMOR.
              !
              error_flag = cmor_close()
              write(*,'(''********************************************************************************'')')
              write(*,'(''********************************************************************************'')')
              write(*,'(''CMOR executed to completion; T#: '',i5,'' X#: '',i5)') itab,ixw
              write(*,'(''********************************************************************************'')')
              write(*,'(''********************************************************************************'')')
              time_counter = 0
              var_counter  = 0
              error_flag   = 0
              var_found    = 0
              scale_factor = 1.
              continue(:)  = .false.
              mycmor%positive = ' '
              original_name= ' '
              !
              if (allocated(time)) then
                 deallocate(time)
                 write(*,*) 'deallocate(time)'
              endif
              if (allocated(time_bnds)) then
                 deallocate(time_bnds)
                 write(*,*) 'deallocate(time_bnds)'
              endif
              if (allocated(indat2anh)) then
                 deallocate(indat2anh)
                 write(*,*) 'deallocate(indat2anh)'
              endif
              if (allocated(indat2ash)) then
                 deallocate(indat2ash)
                 write(*,*) 'deallocate(indat2ash)'
              endif
              if (allocated(indat2bnh)) then
                 deallocate(indat2bnh)
                 write(*,*) 'deallocate(indat2bnh)'
              endif
              if (allocated(indat2bsh)) then
                 deallocate(indat2bsh)
                 write(*,*) 'deallocate(indat2bsh)'
              endif
           endif
        endif
     enddo xwalk_loop
  enddo table_loop
end program OImon_CMOR
