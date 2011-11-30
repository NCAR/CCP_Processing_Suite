program Omon_CMOR
  ! Convert CCSM4 ocn monthly (pop.h) data from single-field format
  ! to CMOR-compliant format
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
  use output_times_info
  use mycmor_info
  !
  implicit none
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  integer::error_flag,cmor_var_id
  real,dimension(:,:)  ,allocatable::indat2a,indat2b,indat2c,cmordat2d
  real,dimension(:,:,:),allocatable::indat3a,indat3b,indat3c,cmordat3d,work3da,work3db
  real,dimension(:,:,:,:),allocatable::indat4a ! MOC
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::time_bnds
  double precision,dimension(1)  ::tval
  double precision,dimension(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile,cmor_filename
  integer::i,j,k,m,n,it,ivar,jvar,length,iexp,jexp,itab,ixw,ifile,ilev,ic
  real::spval
  !
  ! Initialize time indices
  ! 
  tidx1 = -999 ; tidx2 = -999
  !
  ! GO!
  !
  mycmor%table_file = 'Tables/CMIP5_Omon'
  call load_table_info
  !
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
  xwalk_file = 'xwalk_Omon.txt'
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
        xw_found     = 0
        time_units   = ' '
        original_name= ' '
        !
        ! The meaty part
        !
        if (xw(ixw)%entry == table(itab)%variable_entry) then
           write(*,'(''MATCH; CMIP5: '',a,'' CESM: '',5(a))') trim(xw(ixw)%entry),(trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
           do ivar = 1,xw(ixw)%ncesm_vars
              if ((trim(xw(ixw)%cesm_vars(ivar)) == 'UNKNOWN').or.(trim(xw(ixw)%cesm_vars(ivar)) == 'UNAVAILABLE')) then
                 write(*,'(''UNAVAILABLE/UNKNOWN: '',a,'' == '',a)') trim(xw(ixw)%entry),trim(table(itab)%variable_entry)
              else
                 ncfile(:,:)(1:) = ' '
                 nc_nfiles(:)    = 0
                 call build_filenames(ixw,ivar,exp(exp_found)%begyr,exp(exp_found)%endyr)
              endif
           enddo
           !
           ! Open CESM file(s) and get information(s)
           !
           do ivar = 1,xw(ixw)%ncesm_vars
              do ifile = 1,nc_nfiles(ivar)
                 call open_cdf(ncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                 call get_dims(ncid(ifile,ivar))
                 call get_vars(ncid(ifile,ivar))
                 !
                 do n=1,dim_counter
                    length = len_trim(dim_info(n)%name)
                    if(dim_info(n)%name(:length).eq.'time') then
                       ntimes(ifile,ivar) = dim_info(n)%length
                    endif
                 enddo
                 call read_att_text(ncid(ifile,ivar),'time','units',time_units)
                 !
                 do n=1,var_counter
                    if (trim(var_info(n)%name) == trim(xw(ixw)%cesm_vars(ivar))) then
                       var_found(ifile,ivar) = n
                       xw_found = ixw
                    endif
                 enddo
                 if (var_found(ifile,ivar) == 0) then
                    write(*,'(''NEVER FOUND: '',a,'' STOP. '')') trim(xw(ixw)%cesm_vars(ivar))
                    stop
                 endif
                 call close_cdf(ncid(ifile,ivar))
                 !
              enddo
           enddo
           ncid = 0
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
              write(*,*) 'ERROR on cmor_dataset!'
              write(*,*) 'outpath               = ',mycmor%outpath
              write(*,*) 'experiment_id         = ',mycmor%experiment_id
              write(*,*) 'institution           = ',mycmor%institution
              write(*,*) 'source                = ',mycmor%source
              write(*,*) 'calendar              = ',mycmor%calendar
              write(*,*) 'realization           = ',mycmor%realization
              write(*,*) 'contact               = ',mycmor%contact
              write(*,*) 'history               = ',mycmor%history
              write(*,*) 'comment               = ',mycmor%comment
              write(*,*) 'references            = ',mycmor%references
              write(*,*) 'model_id              = ',mycmor%model_id
              write(*,*) 'forcing               = ',mycmor%forcing
              write(*,*) 'initialization_method = ',mycmor%initialization_method
              write(*,*) 'physics_version       = ',mycmor%physics_version
              write(*,*) 'institute_id          = ',mycmor%institute_id
              write(*,*) 'parent_experiment_id  = ',mycmor%parent_experiment_id
              write(*,*) 'parent_experiment_rip = ',mycmor%parent_experiment_rip
              write(*,*) 'branch_time           = ',mycmor%branch_time
           endif
           !
           ! Add global metadata
           !
           call add_global_metadata
           !
           ! Define axes via 'cmor_axis'
           !
           table_ids(2) = cmor_load_table('Tables/CMIP5_grids')
           table_ids(3) = cmor_load_table('Tables/CMIP5_fx')
           call cmor_set_table(table_ids(2))
           call define_ocn_axes(table(itab)%dimensions)
           call cmor_set_table(table_ids(1))
           ! 
           ! Make manual alterations so that CMOR works. Silly code!
           !
           if (xw(ixw)%ncesm_vars == 1) then
              write(original_name,'(a)') xw(ixw)%cesm_vars(1)
           endif
           if (xw(ixw)%ncesm_vars == 2) then
              write(original_name,'(a,'','',a)') (trim(xw(ixw)%cesm_vars(jvar)),jvar=1,xw(ixw)%ncesm_vars)
           endif
           if (xw(ixw)%ncesm_vars == 3) then
              write(original_name,'(a,'','',a,'','',a)') (trim(xw(ixw)%cesm_vars(jvar)),jvar=1,xw(ixw)%ncesm_vars)
           endif
           !
           ! Modify units as necessary to accomodate udunits' inability to convert 
           !
           select case (xw(ixw)%entry)
           case ('msftmyz')
              var_info(var_found(1,1))%units = 'kg s-1'
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
           case ('thetao')
              cmor_var_id = cmor_variable(                            &
                   table=mycmor%table_file,                           &
                   table_entry=xw(ixw)%entry,                         &
                   units=var_info(var_found(1,1))%units,              &
                   axis_ids=(/grid_id(1),axis_ids(3),axis_ids(4)/),   &
                   missing_value=var_info(var_found(1,1))%missing_value,&
                   positive=mycmor%positive,                          &
                   original_name=original_name,                       &
                   comment=xw(ixw)%comment)
           case ('msftmyz')
              cmor_var_id = cmor_variable(                            &
                   table=mycmor%table_file,                           &
                   table_entry=xw(ixw)%entry,                         &
                   units=var_info(var_found(1,1))%units,              &
                   axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3),axis_ids(4)/),   &
                   missing_value=var_info(var_found(1,1))%missing_value,&
                   positive=mycmor%positive,                          &
                   original_name=original_name,                       &
                   comment=xw(ixw)%comment)
           case default
              cmor_var_id = cmor_variable(                            &
                   table=mycmor%table_file,                           &
                   table_entry=xw(ixw)%entry,                         &
                   units=var_info(var_found(1,1))%units,              &
                   axis_ids=(/grid_id(1),axis_ids(3)/),               &
                   missing_value=var_info(var_found(1,1))%missing_value,&
                   positive=mycmor%positive,                          &
                   original_name=original_name,                       &
                   comment=xw(ixw)%comment)
           end select
           write(*,'(''called cmor_variable, table_entry, varid: '',a,2x,i10)') trim(xw(ixw)%entry),cmor_var_id
           !
           ! Perform derivations and cycle through time, writing data too
           !
           select case (xw(ixw)%entry)
           case ('tos')
              !
              ! tos: TEMP at k=1
              !
              allocate(indat3a(nlons,nlats,nlevs),cmordat2d(nlons,nlats))
              do ivar = 1,xw(ixw)%ncesm_vars
                 do ifile = 1,nc_nfiles(ivar)
                    call open_cdf(ncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                    call get_dims(ncid(ifile,ivar))
                    call get_vars(ncid(ifile,ivar))
                    !
                    if (.not.(allocated(time)))      allocate(time(ntimes(ifile,ivar)))
                    if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,ivar)))
                    !
                    write(*,*) 'tbnds loop: ',ivar,ifile,ncid(ifile,ivar)
                    do n = 1,ntimes(ifile,ivar)
                       time_counter = n
                       call read_var(ncid(ifile,ivar),'time_bound',time_bnds(:,n))
                    enddo
                    !
                    time_bnds(1,1) = int(time_bnds(1,1))-1
                    time = (time_bnds(1,:)+time_bnds(2,:))/2.
                    nchunks(ifile)   = 1
                    tidx1(1:nchunks(ifile)) = 1
                    tidx2(1:nchunks(ifile)) = ntimes(ifile,ivar)
                    do ic = 1,nchunks(ifile)
                       do it = tidx1(ic),tidx2(ic)
                          time_counter = it
                          !
                          indat3a = spval
                          call read_var(ncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat3a)
                          cmordat2d = indat3a(:,:,1)
                          !
                          tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                          error_flag = cmor_write(          &
                               var_id        = cmor_var_id, &
                               data          = cmordat2d,   &
                               ntimes_passed = 1,           &
                               time_vals     = tval,        &
                               time_bnds     = tbnd)
                          if (error_flag < 0) then
                             write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                             stop
                          endif
                       enddo
                       write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                    enddo
                    call close_cdf(ncid(ifile,ivar))
                    if (allocated(time))      deallocate(time)
                    if (allocated(time_bnds)) deallocate(time_bnds)
                 enddo
              enddo
              error_flag = cmor_close()
              if (error_flag < 0) then
                 write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
              else
                 write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
              endif
           case ('thetao')
              !
              ! theato: No changes
              !
              allocate(indat3a(nlons,nlats,nlevs))
              do ivar = 1,xw(ixw)%ncesm_vars
                 do ifile = 1,nc_nfiles(ivar)
                    call open_cdf(ncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                    call get_dims(ncid(ifile,ivar))
                    call get_vars(ncid(ifile,ivar))
                    !
                    if (.not.(allocated(time)))      allocate(time(ntimes(ifile,ivar)))
                    if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,ivar)))
                    !
                    do n=1,ntimes(ifile,ivar)
                       time_counter = n
                       call read_var(ncid(ifile,ivar),'time_bound',time_bnds(:,n))
                    enddo
                    !
                    time_bnds(1,1) = int(time_bnds(1,1))-1
                    time = (time_bnds(1,:)+time_bnds(2,:))/2.
                    write(*,*) 'thetao loop: ',ifile,nc_nfiles(ivar),ntimes(ifile,ivar),nlons,nlats,nlevs
                    nchunks(ifile)   = 1
                    tidx1(1:nchunks(ifile)) = 1
                    tidx2(1:nchunks(ifile)) = ntimes(ifile,ivar)
                    do ic = 1,nchunks(ifile)
                       do it = tidx1(ic),tidx2(ic)
                          time_counter = it
                          !
                          indat3a = spval
                          call read_var(ncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat3a)
                          !
                          tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                          error_flag = cmor_write(          &
                               var_id        = cmor_var_id, &
                               data          = indat3a,     &
                               ntimes_passed = 1,           &
                               time_vals     = tval,        &
                               time_bnds     = tbnd)
                          if (error_flag < 0) then
                             write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                             stop
                          endif
                       enddo
                       write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                       !
                       cmor_filename = ' '
                       error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                       if (error_flag < 0) then
                          write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
                          stop
                       else
                          write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
                       endif
                    enddo
                    call close_cdf(ncid(ifile,ivar))
                    if (allocated(time))      deallocate(time)
                    if (allocated(time_bnds)) deallocate(time_bnds)
                 enddo
              enddo
           end select
           if (allocated(indat2a))   deallocate(indat2a)
           if (allocated(indat2b))   deallocate(indat2b)
           if (allocated(indat2c))   deallocate(indat2c)
           if (allocated(cmordat2d)) deallocate(cmordat2d)
           if (allocated(indat3a))   deallocate(indat3a)
           if (allocated(indat3b))   deallocate(indat3b)
           if (allocated(work3da))   deallocate(work3da)
           if (allocated(work3db))   deallocate(work3db)
        endif
        !
        ! Reset
        !
        time_counter = 0
        var_counter  = 0
        error_flag   = 0
        var_found    = 0
        mycmor%positive = ' '
        original_name= ' '
        !
        if (allocated(time))      deallocate(time)
        if (allocated(time_bnds)) deallocate(time_bnds)
        !
!!$        error_flag = cmor_close()
!!$        if (error_flag < 0) then
!!$           write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
!!$        else
!!$           write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
!!$        endif
     enddo xwalk_loop
  enddo table_loop
end program Omon_CMOR
