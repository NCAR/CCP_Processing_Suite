program fx_CMOR
  ! Convert CCSM4 data from single-field format to CMOR-compliant format
  !
  ! NOTE: 'model_id' and first part of 'source' MUST MATCH or CMOR will throw error
  !
  use cmor_users_functions
  use counters_netcdf_jfl
  use interfaces_netcdf_jfl
  use definitions_netcdf_jfl
  use exp_info
  use files_info
  use table_info
  use xwalk_info
  use grid_info
  use mycmor_info
  use output_times_info
  !
  implicit none
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  integer::error_flag,cmor_var_id
  real,dimension(:,:),  allocatable::indat2a,cmordat2d
  real,dimension(:,:,:),allocatable::cmordat3d
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,itab,ixw
  real::spval
  logical,dimension(10)::continue
  !
  ! GO!
  !
  mycmor%table_file = 'Tables/CMIP5_fx'
  call load_table_info
  !
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
  xwalk_file = 'xwalk_fx.txt'
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
  !  call get_atm_grid
  call get_ocn_grid
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
        all_continue = .true.
        continue(:)  = .false.
        original_name= ' '
!
! The meaty part
!
        if (xw(ixw)%entry == table(itab)%variable_entry) then
           do ivar = 1,xw(ixw)%ncesm_vars
              if ((trim(xw(ixw)%cesm_vars(ivar)) == 'UNKNOWN').or.(trim(xw(ixw)%cesm_vars(ivar)) == 'UNAVAILABLE')) then
                 write(*,'('' UNAVAILABLE/UNKNOWN: '',a,'' == '',a)') trim(xw(ixw)%entry),trim(table(itab)%variable_entry)
              else
                 write(ncfile(1,ivar),'(''data/'',a,''.'',a,''.'',a,''.'',a,''01-'',a,''12.nc'')') &
                      trim(case_read),&
                      trim(comp_read),&
                      trim(xw(ixw)%cesm_vars(ivar)),&
                      exp(exp_found)%begin_end(1:4),&
                      exp(exp_found)%begin_end(6:9)
                 inquire(file=trim(ncfile(1,ivar)),exist=continue(ivar))
                 if (.not.(continue(ivar))) then
                    write(ncfile(1,ivar),'(''data/'',a,''.'',a,''.'',a,''.'',a,''-01_cat_'',a,''-12.nc'')') &
                         trim(case_read),&
                         trim(comp_read),&
                         trim(xw(ixw)%cesm_vars(ivar)),&
                         exp(exp_found)%begin_end(1:4),&
                         exp(exp_found)%begin_end(6:9)
                    inquire(file=trim(ncfile(1,ivar)),exist=continue(ivar))
                 endif
                 if (.not.(continue(ivar))) then
                    write(*,*) trim(ncfile(1,ivar)),' NOT FOUND.'
                 else
                    write(*,'('' GOOD TO GO : '',a,'' == '',a,'' from CESM file: '',a)') &
                         trim(xw(ixw)%entry),&
                         trim(table(itab)%variable_entry),&
                         trim(ncfile(1,ivar))
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
                 call open_cdf(myncid(1,ivar),trim(ncfile(1,ivar)),.true.)
                 write(*,'(''OPENING: '',a80,'' myncid: '',i10)') trim(ncfile(1,ivar)),myncid(1,ivar)
                 call get_dims(myncid(1,ivar))
                 call get_vars(myncid(1,ivar))
                 !
                 do n=1,var_counter
                    if (trim(var_info(n)%name) == trim(xw(ixw)%cesm_vars(ivar))) then
                       var_found(1,ivar) = n
                    endif
                 enddo
                 if (var_found(1,ivar) == 0) then
                    write(*,'(''NEVER FOUND: '',a,'' STOP. '')') trim(xw(ixw)%cesm_vars(ivar))
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
              !
              error_flag = cmor_dataset(                              &
                   outpath=mycmor%outpath,                            &
                   experiment_id=mycmor%experiment_id,                &
                   institution=mycmor%institution,                    &
                   source=mycmor%source,                              &
                   calendar=mycmor%calendar,                          &
                   realization=0,                                     &
                   contact=mycmor%contact,                            &
                   history=mycmor%history,                            &
                   comment=mycmor%comment,                            &
                   references=mycmor%references,                      &
                   model_id=mycmor%model_id,                          &
                   forcing=mycmor%forcing,                            &
                   initialization_method=0,                           &
                   physics_version=0,                                 &
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
              call define_ocn_axes(table(itab)%dimensions)
              call cmor_set_table(table_ids(1))
              !
              ! Define axes via 'cmor_axis'
              !
              !              call define_atm_axes(table(itab)%dimensions)
              if (xw(ixw)%ncesm_vars == 1) write(original_name,'(a)') xw(ixw)%cesm_vars(1)
              if (xw(ixw)%ncesm_vars == 2) write(original_name,'(a,'','',a)') (trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
              if (xw(ixw)%ncesm_vars == 3) &
                   write(original_name,'(a,'','',a,'','',a)') (trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
              ! 
              ! Make manual alterations so that CMOR works. Silly code!
              !
              select case (xw(ixw)%entry)
              case ('sftlf')
                 var_info(var_found(1,1))%units = '1'
              case ('volcello')
                 var_info(var_found(1,1))%units = 'm3'
              end select
              !
              spval=var_info(var_found(1,1))%missing_value
              !
              write(*,*) 'calling cmor_variable:'
              write(*,*) 'table         = ',trim(mycmor%table_file)
              write(*,*) 'table_entry   = ',trim(xw(ixw)%entry)
              write(*,*) 'dimensions    = ',trim(table(itab)%dimensions)
              write(*,*) 'units         = ',trim(var_info(var_found(1,1))%units)
              write(*,*) 'axis_ids      = ',axis_ids(1:4)
              write(*,*) 'missing_value = ',var_info(var_found(1,1))%missing_value
              write(*,*) 'positive      = ',trim(mycmor%positive)
              write(*,*) 'original_name = ',trim(original_name)
              !
              select case (xw(ixw)%entry)
              case ('deptho','areacello')
                 cmor_var_id = cmor_variable(                &
                      table=mycmor%table_file,               &
                      table_entry=xw(ixw)%entry,             &
                      units=var_info(var_found(1,1))%units,  &
                      axis_ids=(/grid_id(1)/),               &
                      missing_value=spval,                   &
                      positive=mycmor%positive,              &
                      original_name=original_name,           &
                      comment=xw(ixw)%comment)
              case ('volcello')
                 cmor_var_id = cmor_variable(                &
                      table=mycmor%table_file,               &
                      table_entry=xw(ixw)%entry,             &
                      units=var_info(var_found(1,1))%units,  &
                      axis_ids=(/grid_id(1),axis_ids(3)/),   &
                      missing_value=spval,                   &
                      positive=mycmor%positive,              &
                      original_name=original_name,           &
                      comment=xw(ixw)%comment)
              end select
!!$              !
!!$              cmor_var_id = cmor_variable(                                &
!!$                   table=mycmor%table_file,                           &
!!$                   table_entry=xw(ixw)%entry,                         &
!!$                   units=var_info(var_found(1,ivar))%units,                &
!!$                   axis_ids=(/ axis_ids(1), axis_ids(2) /),    &
!!$                   missing_value=var_info(var_found(1,ivar))%missing_value,&
!!$                   positive=mycmor%positive,                          &
!!$                   original_name=original_name,                       &
!!$                   comment=xw(ixw)%comment)
              select case (xw(ixw)%entry)
              case ('deptho','areacello')
                 allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
                 write(*,*) 'TO READ: ',trim(var_info(var_found(1,1))%name)
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 write(*,*) 'READ: ',trim(var_info(var_found(1,1))%name),maxval(indat2a),minval(indat2a)
                 where (kmt == 0)
                    cmordat2d = spval
                 elsewhere
                    cmordat2d = indat2a
                 endwhere
                 error_flag = cmor_write(      &
                      var_id        = cmor_var_id, &
                      data          = cmordat2d)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              case ('volcello')
                 allocate(indat2a(nlons,nlats),cmordat3d(nlons,nlats,nlevs))
                 write(*,*) 'TO READ: ',trim(var_info(var_found(1,1))%name)
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 cmordat3d = spval
                 do k = 1,nlevs
                    do j = 1,nlats
                       do i = 1,nlons
                          if (kmt(i,j).ge.k) then
                             cmordat3d(i,j,k) = (indat2a(i,j)*ocn_levs(k))/(100. * 100. * 100.)
                          else
                             cmordat3d(i,j,k) = spval
                          endif
                       enddo
                    enddo
                 enddo
                 !
                 error_flag = cmor_write(      &
                      var_id        = cmor_var_id, &
                      data          = cmordat3d)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              end select
              do ivar = 1,xw(ixw)%ncesm_vars
                 call close_cdf(myncid(1,ivar))
              enddo
              !
              ! Close all files opened by CMOR.
              !
              !
              error_flag = cmor_close()
              if (error_flag < 0) then
                 write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
              else
                 write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
              endif
              time_counter = 0
              var_counter  = 0
              error_flag   = 0
              var_found    = 0
              continue(:)  = .false.
              mycmor%positive = ' '
              original_name= ' '
              !
              if (allocated(indat2a))   deallocate(indat2a)
              if (allocated(cmordat2d)) deallocate(cmordat2d)
              if (allocated(cmordat3d)) deallocate(cmordat3d)
           endif
        endif
     enddo xwalk_loop
  enddo table_loop
end program fx_CMOR
