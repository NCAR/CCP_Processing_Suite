program Lmon_CMOR
  ! Convert CCSM4 land monthly (clm2.h0) data from single-field format
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
  real,dimension(:,:)  ,allocatable::indat2a,indat2b,indat2c,cmordat2d
  real,dimension(:,:,:),allocatable::indat3a,indat3b,indat3c,cmordat3d,work3da,work3db
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::time_bnds
  double precision,dimension(1)  ::tval
  double precision,dimension(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,itab,ixw,ilev
  real::spval
  logical::all_continue
  !
  character(len=256),dimension(10)::ncfile
  real,dimension(10)::allmax,allmin,scale_factor
  integer,dimension(10)::ncid,var_found
  logical,dimension(10)::continue
  !
  ! GO!
  !
  mycmor%table_file = 'Tables/CMIP5_Lmon'
  call load_table_info
  !
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
  xwalk_file = 'xwalk_Lmon.txt'
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
  call get_lnd_grid
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
                 if (.not.(continue(ivar))) then
                    write(*,*) trim(ncfile(ivar)),' NOT FOUND.'
                 else
                    write(*,'('' GOOD TO GO : '',a,'' == '',a,'' from CESM file: '',a)') &
                         trim(xw(ixw)%entry),&
                         trim(table(itab)%variable_entry),&
                         trim(ncfile(ivar))
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
                    call read_var(ncid(ivar),'time_bounds',time_bnds(:,n))
                    if (n == 1) time_bnds(1,n) = 0.
                    time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
!                    write(*,'(''TIMES: '',3f12.4)') time_bnds(1,n),time(n),time_bnds(2,n)
                 enddo
                 write(*,*) 'time_bnds calculated'
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
              call define_lnd_axes(table(itab)%dimensions)
              ! 
              ! Make manual alterations so that CMOR works. Silly code!
              !
              if (xw(ixw)%ncesm_vars == 1) then
                 write(original_name,'(a)') xw(ixw)%cesm_vars(1)
              endif
              if (xw(ixw)%ncesm_vars == 2) then
                 write(original_name,'(a,'','',a)') (trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
              endif
              if (xw(ixw)%ncesm_vars == 3) then
                 write(original_name,'(a,'','',a,'','',a)') (trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
              endif
              !
              ! Modify units as necessary to accomodate udunits' inability to convert 
              !
              select case (xw(ixw)%entry)
              case ('cVeg','cLitter','cSoil','cProduct')
                 var_info(var_found(1))%units = 'kg m-2'
              case ('gpp')
                 mycmor%positive = 'down'
                 var_info(var_found(1))%units = 'kg m-2 s-1'
              case ('evspsblveg','evspsblsoi','tran')
                 ! mm s-1 is the same as kg m-2 s-1
                 var_info(var_found(1))%units = 'kg m-2 s-1'
                 mycmor%positive = 'up'
              case ('mrro','mrros','prveg')
                 ! mm s-1 is the same as kg m-2 s-1
                 var_info(var_found(1))%units = 'kg m-2 s-1'
              case ('burntArea')
                 ! Units 'proportion' replaced by 'something'
                 var_info(var_found(1))%units = '%'
              case ('lai')
                 ! Units 'none' replaced by '1'
                 var_info(var_found(1))%units = '1'
              end select
              !
              spval=var_info(var_found(1))%missing_value
              select case (xw(ixw)%entry)
                 case ('tsl')
                    var_ids = cmor_variable(                                &
                         table=mycmor%table_file,                           &
                         table_entry=xw(ixw)%entry,                         &
                         units=var_info(var_found(1))%units,                &
                         axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3),axis_ids(4)/),  &
                         missing_value=var_info(var_found(1))%missing_value,&
                         positive=mycmor%positive,                          &
                         original_name=original_name,                       &
                         comment=xw(ixw)%comment)
                 case default
                    var_ids = cmor_variable(                                &
                         table=mycmor%table_file,                           &
                         table_entry=xw(ixw)%entry,                         &
                         units=var_info(var_found(1))%units,                &
                         axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3)/),  &
                         missing_value=var_info(var_found(1))%missing_value,&
                         positive=mycmor%positive,                          &
                         original_name=original_name,                       &
                         comment=xw(ixw)%comment)
                 end select
              !
              write(*,*) 'cmor_variable:'
              write(*,*) 'varids=',var_ids
              write(*,*) 'table=',trim(mycmor%table_file)
              write(*,*) 'table_entry=',trim(xw(ixw)%entry)
              write(*,*) 'dimensions=',trim(table(itab)%dimensions)
              write(*,*) 'units=',trim(var_info(var_found(1))%units)
              write(*,*) 'missing_value=',var_info(var_found(1))%missing_value
              write(*,*) 'positive=',trim(mycmor%positive)
              write(*,*) 'original_name=',trim(original_name)
              !
              ! Cycle through time
              !
              time_loop: DO it=1, ntimes
                 time_counter = it
                 !
                 ! Perform derivations
                 !
                 select case (xw(ixw)%entry)
                 case ('evspsblveg','evspsblsoi','tran','mrros','prveg','lai')
                    !
                    ! No change
                    !
                    allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
                    call read_var(ncid(1),var_info(var_found(1))%name,indat2a)
                    allmax(1) = max(allmax(1),maxval(indat2a)) ; allmin(1) = min(allmin(1),minval(indat2a))
                    ! 
                    cmordat2d = indat2a
                 case ('cVeg','cLitter','cSoil','cProduct','gpp')
                    !
                    ! Unit change - grams to kg
                    !
                    allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
                    call read_var(ncid(1),var_info(var_found(1))%name,indat2a)
                    allmax(1) = max(allmax(1),maxval(indat2a)) ; allmin(1) = min(allmin(1),minval(indat2a))
                    ! 
                    cmordat2d = indat2a/1000.
                 case ('burntArea')
                    !
                    ! Unit change - 'proportion' to percentage
                    !
                    allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
                    call read_var(ncid(1),var_info(var_found(1))%name,indat2a)
                    allmax(1) = max(allmax(1),maxval(indat2a)) ; allmin(1) = min(allmin(1),minval(indat2a))
                    !
                    cmordat2d = indat2a*100.
                 case ('mrro')
                    !
                    ! Add QOVER + QDRAI + QRGWL, units result in no numerical change (mm s-1 to kg m-2 s-1)
                    !
                    allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats))
                    allocate(cmordat2d(nlons,nlats))
                    call read_var(ncid(1),var_info(var_found(1))%name,indat2a)
                    call read_var(ncid(2),var_info(var_found(2))%name,indat2b)
                    call read_var(ncid(3),var_info(var_found(3))%name,indat2c)
                    allmax(1) = max(allmax(1),maxval(indat2a,mask=indat2a/=spval)) ; allmin(1) = min(allmin(1),minval(indat2a,mask=indat2a/=spval))
                    allmax(2) = max(allmax(2),maxval(indat2b,mask=indat2b/=spval)) ; allmin(2) = min(allmin(2),minval(indat2b,mask=indat2b/=spval))
                    allmax(3) = max(allmax(3),maxval(indat2c,mask=indat2c/=spval)) ; allmin(3) = min(allmin(3),minval(indat2c,mask=indat2c/=spval))
                    ! 
                    where ((indat2a /= spval).and.(indat2b /= spval).and.(indat2c /= spval))
                       cmordat2d = indat2a + indat2b + indat2c
                    elsewhere
                       cmordat2d = spval
                    endwhere
                 case ('mrsos')
                    !
                    ! Integrate SOILICE and SOILIQ over top 10 cm
                    !
                    allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs))
                    allocate(work3da(nlons,nlats,nlevs),work3db(nlons,nlats,nlevs))
                    allocate(cmordat2d(nlons,nlats))
                    call read_var(ncid(1),var_info(var_found(1))%name,indat3a)
                    call read_var(ncid(2),var_info(var_found(2))%name,indat3b)
                    allmax(1) = max(allmax(1),maxval(indat3a)) ; allmin(1) = min(allmin(1),minval(indat3a))
                    allmax(2) = max(allmax(2),maxval(indat3b)) ; allmin(2) = min(allmin(2),minval(indat3b))
                    work3da = 0. ; work3db = 0.
                    do k = 1,4
                       do j = 1,nlats
                          do i = 1,nlons
                             if (indat3a(i,j,k) /= spval) work3da(i,j,k) = indat3a(i,j,k)*lnd_dzsoi(i,j,k)
                             if (indat3b(i,j,k) /= spval) work3db(i,j,k) = indat3b(i,j,k)*lnd_dzsoi(i,j,k)
                          enddo
                       enddo
                    enddo
                    cmordat2d = (sum(work3da,dim=3) + sum(work3db,dim=3))/sum(lnd_levs(1:4))
                    write(*,*) 'mrsos at ',it,' X, N: ',maxval(cmordat2d,mask=cmordat2d/=spval),minval(cmordat2d,mask=cmordat2d/=spval)
                 case ('mrso')
                    !
                    ! Integrate SOILICE and SOILIQ over all layers
                    !
                    allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs))
                    allocate(work3da(nlons,nlats,nlevs),work3db(nlons,nlats,nlevs))
                    allocate(cmordat2d(nlons,nlats))
                    call read_var(ncid(1),var_info(var_found(1))%name,indat3a)
                    call read_var(ncid(2),var_info(var_found(2))%name,indat3b)
                    allmax(1) = max(allmax(1),maxval(indat3a)) ; allmin(1) = min(allmin(1),minval(indat3a))
                    allmax(2) = max(allmax(2),maxval(indat3b)) ; allmin(2) = min(allmin(2),minval(indat3b))
                    work3da = 0. ; work3db = 0.
                    do k = 1,nlevs
                       do j = 1,nlats
                          do i = 1,nlons
                             if (indat3a(i,j,k) /= spval) work3da(i,j,k) = indat3a(i,j,k)*lnd_dzsoi(i,j,k)
                             if (indat3b(i,j,k) /= spval) work3db(i,j,k) = indat3b(i,j,k)*lnd_dzsoi(i,j,k)
                          enddo
                       enddo
                    enddo
                    cmordat2d = (sum(work3da,dim=3) + sum(work3db,dim=3))/sum(lnd_dzsoi,dim=3)
                    write(*,*) 'mrso at ',it,' X, N: ',maxval(cmordat2d,mask=cmordat2d/=spval),minval(cmordat2d,mask=cmordat2d/=spval)
                 case ('mrfso')
                    !
                    ! Integrate SOILICE over all layers
                    !
                    allocate(indat3a(nlons,nlats,nlevs),work3da(nlons,nlats,nlevs))
                    allocate(cmordat2d(nlons,nlats))
                    call read_var(ncid(1),var_info(var_found(1))%name,indat3a)
                    allmax(1) = max(allmax(1),maxval(indat3a)) ; allmin(1) = min(allmin(1),minval(indat3a))
                    work3da = 0.
                    do k = 1,nlevs
                       do j = 1,nlats
                          do i = 1,nlons
                             if (indat3a(i,j,k) /= spval) work3da(i,j,k) = indat3a(i,j,k)*lnd_dzsoi(i,j,k)
                          enddo
                       enddo
                    enddo
                    cmordat2d = sum(work3da,dim=3)/sum(lnd_dzsoi,dim=3)
                    write(*,*) 'mrfso at ',it,' X, N: ',maxval(cmordat2d,mask=cmordat2d/=spval),minval(cmordat2d,mask=cmordat2d/=spval)
                 case ('tsl')
                    !
                    ! Pass TSOI straight through
                    !
                    allocate(indat3a(nlons,nlats,nlevs),cmordat3d(nlons,nlats,nlevs))
                    call read_var(ncid(1),var_info(var_found(1))%name,indat3a)
                    cmordat3d = indat3a
                    write(*,*) 'tsl at ',it,' X, N: ',maxval(cmordat3d,mask=cmordat3d/=spval),minval(cmordat3d,mask=cmordat3d/=spval)
                 end select
                 !
                 tval(1)   = time(it)
                 tbnd(1,1) = time_bnds(1,it)
                 tbnd(2,1) = time_bnds(2,it)
                 !
                 ! Pass data to cmor_write
                 !
                 select case (xw(ixw)%entry)
                 case ('evspsblveg','evspsblsoi','tran','mrros','prveg','lai','cVeg','cLitter','cSoil','cProduct','gpp','burntArea','mrro','mrsos','mrso','mrfso')
                    error_flag = cmor_write(      &
                         var_id        = var_ids, &
                         data          = cmordat2d, &
                         ntimes_passed = 1,       &
                         time_vals     = tval,    &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,*) 'Error writing ',xw(ixw)%entry, ', which I call ', xw(ixw)%cesm_vars
                       write(*,*) 'Processing time sample: ', time
                       stop
                    endif
                 case ('tsl')
                    error_flag = cmor_write(      &
                         var_id        = var_ids, &
                         data          = cmordat3d,&
                         ntimes_passed = 1,       &
                         time_vals     = tval,    &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,*) 'Error writing ',xw(ixw)%entry, ', which I call ', xw(ixw)%cesm_vars
                       write(*,*) 'Processing time sample: ', time
                       stop
                    endif
                 end select
                 if (allocated(indat2a))   deallocate(indat2a)
                 if (allocated(indat2b))   deallocate(indat2b)
                 if (allocated(indat2c))   deallocate(indat2c)
                 if (allocated(cmordat2d)) deallocate(cmordat2d)
                 if (allocated(indat3a))   deallocate(indat3a)
                 if (allocated(indat3b))   deallocate(indat3b)
                 if (allocated(work3da))   deallocate(work3da)
                 if (allocated(work3db))   deallocate(work3db)
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
              if (allocated(time))      deallocate(time)
              if (allocated(time_bnds)) deallocate(time_bnds)
           endif
        endif
     enddo xwalk_loop
  enddo table_loop
end program Lmon_CMOR
