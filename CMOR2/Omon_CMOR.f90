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
  integer::i,j,k,m,n,it,ivar,jvar,length,iexp,jexp,ixw,ilev,ic
  real::spval
  logical::does_exist
  !
  ! GO!
  !
  mycmor%table_file = 'Omon'
  xwalk_file = 'xwalk_'//trim(mycmor%table_file)//'.txt'
  !
  ! Get table information
  !
  mycmor%table_file = 'Tables/CMIP5_'//trim(mycmor%table_file)
  inquire(file=mycmor%table_file,exist=does_exist)
  if (.not.(does_exist)) then
     write(*,*) 'Cannot find ',trim(mycmor%table_file),'. Dying.'
     stop
  endif
  !
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
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
  xwalk_loop: do ixw = 1,num_xw
     call reset_netcdf_var
     mycmor%positive = ' '
     error_flag      = 0
     time_units      = ' '
     !
     ! The meaty part
     !
     if (xw(ixw)%ncesm_vars == 0) then
        write(*,'(a,'' is UNAVAILABLE.'')') trim(xw(ixw)%entry)
        all_continue = .false.
     endif
     !
     do ivar = 1,xw(ixw)%ncesm_vars
        if (trim(xw(ixw)%cesm_vars(ivar)) == 'UNKNOWN') then
           write(*,'(a,'' has UNKNOWN equivalence.'')') trim(xw(ixw)%entry)
           xw(ixw)%ncesm_vars = 0
           all_continue = .false.
        else
           write(*,'(''CHECKING AVAILABILITY OF: '',a,''.'',a,''.'',a,''.* FILES'')') trim(case_read),trim(comp_read),trim(xw(ixw)%cesm_vars(ivar))
           call build_filenames(case_read,comp_read,xw(ixw)%cesm_vars(ivar),ivar,exp(exp_found)%begyr,exp(exp_found)%endyr,mycmor%table_file)
        endif
     enddo
     !
     ! Open CESM file(s) and get information(s)
     !
     if (all_continue) then
        do ivar = 1,xw(ixw)%ncesm_vars
           do ifile = 1,nc_nfiles(ivar)
              call open_cdf(myncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
              call get_dims(myncid(ifile,ivar))
              call get_vars(myncid(ifile,ivar))
              !
              do n=1,dim_counter
                 length = len_trim(dim_info(n)%name)
                 if(dim_info(n)%name(:length).eq.'time') then
                    ntimes(ifile,ivar) = dim_info(n)%length
                 endif
              enddo
              call read_att_text(myncid(ifile,ivar),'time','units',time_units)
              write(*,'(''read_att_text, time_units: '',a)') trim(time_units)
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
              call close_cdf(myncid(ifile,ivar))
              !
           enddo
        enddo
        !
        myncid = 0
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
        call define_ocn_axes(xw(ixw)%dims)
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
        case ('msftmyz','msftbarot')
           var_info(var_found(1,1))%units = 'kg s-1'
        case ('hfss','rlds','rsntds','rsds','tauuo','tauvo')
           mycmor%positive = 'up'
        end select
        !
        spval=var_info(var_found(1,1))%missing_value
        !
        write(*,*) 'calling cmor_variable:'
        write(*,*) 'table         = ',trim(mycmor%table_file)
        write(*,*) 'table_entry   = ',trim(xw(ixw)%entry)
        write(*,*) 'dimensions    = ',trim(xw(ixw)%dims)
        write(*,*) 'units         = ',trim(var_info(var_found(1,1))%units)
        write(*,*) 'axis_ids      = ',axis_ids(1:naxes)
        write(*,*) 'missing_value = ',var_info(var_found(1,1))%missing_value
        write(*,*) 'positive      = ',trim(mycmor%positive)
        write(*,*) 'original_name = ',trim(original_name)
        !
        select case (xw(ixw)%entry)
        case ('thetao','so','agessc','uo','vo','rhopoto','cfc11')
           ! Full-column fields
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
           ! Different MOC grid
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
           ! Single-level fields
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
           if (.not.(allocated(indat3a)))   allocate(indat3a(nlons,nlats,nlevs))
           if (.not.(allocated(cmordat2d))) allocate(cmordat2d(nlons,nlats))
           do ivar = 1,xw(ixw)%ncesm_vars
              do ifile = 1,nc_nfiles(ivar)
                 call open_cdf(myncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                 call get_dims(myncid(ifile,ivar))
                 call get_vars(myncid(ifile,ivar))
                 !
                 if (.not.(allocated(time)))      allocate(time(ntimes(ifile,ivar)))
                 if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,ivar)))
                 !
                 do n = 1,ntimes(ifile,ivar)
                    time_counter = n
                    call read_var(myncid(ifile,ivar),'time_bound',time_bnds(:,n))
                 enddo
                 !
                 time_bnds(1,1) = int(time_bnds(1,1))-1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 select case (ntimes(ifile,ivar))
                 case ( 1152 ) ! RCP, skip 2005
                    nchunks(ifile) = 1
                    tidx1(1:nchunks(ifile)) = 13
                 case default
                    nchunks(ifile)   = 1
                    tidx1(1:nchunks(ifile)) =  1
                 end select
                 tidx2(nchunks(ifile)) = ntimes(ifile,ivar)
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       !
                       indat3a = spval
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat3a)
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
                 call close_cdf(myncid(ifile,ivar))
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
        case ('sos')
           !
           ! sos: SALT at k=1, multiply by 1000 to maintain identity between g kg-1 and psu
           !
           if (.not.(allocated(indat3a)))   allocate(indat3a(nlons,nlats,nlevs))
           if (.not.(allocated(cmordat2d))) allocate(cmordat2d(nlons,nlats))
           do ivar = 1,xw(ixw)%ncesm_vars
              do ifile = 1,nc_nfiles(ivar)
                 call open_cdf(myncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                 call get_dims(myncid(ifile,ivar))
                 call get_vars(myncid(ifile,ivar))
                 !
                 if (.not.(allocated(time)))      allocate(time(ntimes(ifile,ivar)))
                 if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,ivar)))
                 !
                 do n = 1,ntimes(ifile,ivar)
                    time_counter = n
                    call read_var(myncid(ifile,ivar),'time_bound',time_bnds(:,n))
                 enddo
                 !
                 time_bnds(1,1) = int(time_bnds(1,1))-1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 select case (ntimes(ifile,ivar))
                 case ( 1152 ) ! RCP, skip 2005
                    nchunks(ifile) = 1
                    tidx1(1:nchunks(ifile)) = 13
                 case default
                    nchunks(ifile)   = 1
                    tidx1(1:nchunks(ifile)) =  1
                 end select
                 tidx2(nchunks(ifile)) = ntimes(ifile,ivar)
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       !
                       indat3a = spval
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat3a)
                       where (indat3a(:,:,1) /= spval)
                          cmordat2d = indat3a(:,:,1)*1000.
                       endwhere
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
                 call close_cdf(myncid(ifile,ivar))
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
        case ('tauuo','tauvo','hfss','pr','prsn','rsds','rsntds')
           !
           ! tauuo,tauvo,hfss,pr,prsn,rsds,rsntds: No changes
           !
           allocate(indat2a(nlons,nlats))
           do ivar = 1,xw(ixw)%ncesm_vars
              do ifile = 1,nc_nfiles(ivar)
                 call open_cdf(myncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                 call get_dims(myncid(ifile,ivar))
                 call get_vars(myncid(ifile,ivar))
                 !
                 if (.not.(allocated(time)))      allocate(time(ntimes(ifile,ivar)))
                 if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,ivar)))
                 !
                 do n=1,ntimes(ifile,ivar)
                    time_counter = n
                    call read_var(myncid(ifile,ivar),'time_bound',time_bnds(:,n))
                 enddo
                 !
                 time_bnds(1,1) = int(time_bnds(1,1))-1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 select case (ntimes(ifile,ivar))
                 case ( 1152 ) ! RCP, skip 2005
                    nchunks(ifile) = 1
                    tidx1(1:nchunks(ifile)) = 13
                 case default
                    nchunks(ifile)   = 1
                    tidx1(1:nchunks(ifile)) =  1
                 end select
                 tidx2(nchunks(ifile)) = ntimes(ifile,ivar)
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       !
                       indat2a = spval
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat2a)
                       !
                       tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                       error_flag = cmor_write(          &
                            var_id        = cmor_var_id, &
                            data          = indat2a,     &
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
                 call close_cdf(myncid(ifile,ivar))
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
        case ('msftbarot')
           !
           ! msftbarot: Convert BSF from Sv to kg s-1
           !
           if (.not.(allocated(indat2a))) allocate(indat2a(nlons,nlats))
           do ivar = 1,xw(ixw)%ncesm_vars
              do ifile = 1,nc_nfiles(ivar)
                 call open_cdf(myncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                 call get_dims(myncid(ifile,ivar))
                 call get_vars(myncid(ifile,ivar))
                 !
                 if (.not.(allocated(time)))      allocate(time(ntimes(ifile,ivar)))
                 if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,ivar)))
                 !
                 do n=1,ntimes(ifile,ivar)
                    time_counter = n
                    call read_var(myncid(ifile,ivar),'time_bound',time_bnds(:,n))
                 enddo
                 !
                 time_bnds(1,1) = int(time_bnds(1,1))-1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 select case (ntimes(ifile,ivar))
                 case ( 1152 ) ! RCP, skip 2005
                    nchunks(ifile) = 1
                    tidx1(1:nchunks(ifile)) = 13
                 case default
                    nchunks(ifile)   = 1
                    tidx1(1:nchunks(ifile)) =  1
                 end select
                 tidx2(nchunks(ifile)) = ntimes(ifile,ivar)
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       !
                       indat2a = spval
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat2a)
                       where (indat2a /= spval)
                          indat2a = indat2a * (1.e6 * 1000.) ! 10^6 m3 s-1 to kg s-1
                       endwhere
                       !
                       tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                       error_flag = cmor_write(          &
                            var_id        = cmor_var_id, &
                            data          = indat2a,     &
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
                    if ((it-1) /= tidx2(ic)) then
                       cmor_filename = ' '
                       error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                       if (error_flag < 0) then
                          write(*,'(''ERROR close & preserve: '',a)') cmor_filename(1:128)
                          stop
                       else
                          write(*,'('' GOOD close & preserve: '',a)') cmor_filename(1:128)
                       endif
                    endif
                 enddo
                 call close_cdf(myncid(ifile,ivar))
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
        case ('thetao','agessc','uo','vo','rhopoto')
           !
           ! thetao,agessc,uo,vo,rhopoto: No changes
           !
           allocate(indat3a(nlons,nlats,nlevs))
           do ivar = 1,xw(ixw)%ncesm_vars
              do ifile = 1,nc_nfiles(ivar)
                 call open_cdf(myncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                 call get_dims(myncid(ifile,ivar))
                 call get_vars(myncid(ifile,ivar))
                 !
                 if (.not.(allocated(time)))      allocate(time(ntimes(ifile,ivar)))
                 if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,ivar)))
                 !
                 do n=1,ntimes(ifile,ivar)
                    time_counter = n
                    call read_var(myncid(ifile,ivar),'time_bound',time_bnds(:,n))
                 enddo
                 !
                 time_bnds(1,1) = int(time_bnds(1,1))-1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 select case (ntimes(ifile,ivar))
                 case ( 1152 )               ! RCP all in one file from 2005-2100
                    nchunks(ifile) = 10
                    tidx1(1) =  13
                    tidx2(1) =  60
                    do ic = 2,nchunks(ifile)
                       tidx1(ic) = tidx2(ic-1) + 1
                       tidx2(ic) = tidx1(ic) + 119
                    enddo
                 case ( 1872 ) ! 20C from 1850-2005
                    nchunks(ifile) = 16
                    tidx1(1) =   1
                    tidx2(1) = 120
                    do ic = 2,nchunks(ifile)
                       tidx1(ic) = tidx2(ic-1) + 1
                       tidx2(ic) = tidx1(ic) + 119
                    enddo
                 case default
                    nchunks(ifile)   = 1
                    tidx1(1:nchunks(ifile)) =  1
                 end select
                 !
                 tidx2(nchunks(ifile)) = ntimes(ifile,ivar)
                 !
                 write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       !
                       indat3a = spval
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat3a)
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
                 call close_cdf(myncid(ifile,ivar))
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
        case ('so')
           !
           ! so: SALT*1000 to handle CMOR change to psu bug
           !
           allocate(indat3a(nlons,nlats,nlevs))
           do ivar = 1,xw(ixw)%ncesm_vars
              do ifile = 1,nc_nfiles(ivar)
                 call open_cdf(myncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                 call get_dims(myncid(ifile,ivar))
                 call get_vars(myncid(ifile,ivar))
                 !
                 if (.not.(allocated(time)))      allocate(time(ntimes(ifile,ivar)))
                 if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,ivar)))
                 !
                 do n=1,ntimes(ifile,ivar)
                    time_counter = n
                    call read_var(myncid(ifile,ivar),'time_bound',time_bnds(:,n))
                 enddo
                 !
                 time_bnds(1,1) = int(time_bnds(1,1))-1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 select case (ntimes(ifile,ivar))
                 case ( 1152 )               ! RCP all in one file from 2005-2100
                    nchunks(ifile) = 10
                    tidx1(1) =  13
                    tidx2(1) =  60
                    do ic = 2,nchunks(ifile)
                       tidx1(ic) = tidx2(ic-1) + 1
                       tidx2(ic) = tidx1(ic) + 119
                    enddo
                 case ( 1872 ) ! 20C from 1850-2005
                    nchunks(ifile) = 16
                    tidx1(1) =   1
                    tidx2(1) = 120
                    do ic = 2,nchunks(ifile)
                       tidx1(ic) = tidx2(ic-1) + 1
                       tidx2(ic) = tidx1(ic) + 119
                    enddo
                 case default
                    nchunks(ifile)   = 1
                    tidx1(1:nchunks(ifile)) =  1
                 end select
                 tidx2(nchunks(ifile)) = ntimes(ifile,ivar)
                 !
                 write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       !
                       indat3a = spval
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat3a)
                       where (indat3a /= spval)
                          indat3a = indat3a*1000.
                       endwhere
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
                 call close_cdf(myncid(ifile,ivar))
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
!!$           case ('rlds')
!!$              !
!!$              ! rlds: LWUP_F - LWDN_F
!!$              !
!!$              allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),cmordat2d(nlons,nlats))
!!$              do ivar = 1,xw(ixw)%ncesm_vars
!!$                 do ifile = 1,nc_nfiles(ivar)
!!$                    call open_cdf(myncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
!!$                    call get_dims(myncid(ifile,ivar))
!!$                    call get_vars(myncid(ifile,ivar))
!!$                    !
!!$                    if (.not.(allocated(time)))      allocate(time(ntimes(ifile,ivar)))
!!$                    if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,ivar)))
!!$                    !
!!$                    do n=1,ntimes(ifile,ivar)
!!$                       time_counter = n
!!$                       call read_var(myncid(ifile,ivar),'time_bound',time_bnds(:,n))
!!$                    enddo
!!$                    !
!!$                    time_bnds(1,1) = int(time_bnds(1,1))-1
!!$                    time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$                    write(*,*) 'thetao loop: ',ifile,nc_nfiles(ivar),ntimes(ifile,ivar),nlons,nlats,nlevs
!!$                    nchunks(ifile)   = 1
!!$                    tidx1(1:nchunks(ifile)) = 1
!!$                    tidx2(nchunks(ifile)) = ntimes(ifile,ivar)
!!$                    do ic = 1,nchunks(ifile)
!!$                       do it = tidx1(ic),tidx2(ic)
!!$                          time_counter = it
!!$                          !
!!$                          indat3a = spval
!!$                          call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat3a)
!!$                          !
!!$                          tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                          error_flag = cmor_write(          &
!!$                               var_id        = cmor_var_id, &
!!$                               data          = indat3a,     &
!!$                               ntimes_passed = 1,           &
!!$                               time_vals     = tval,        &
!!$                               time_bnds     = tbnd)
!!$                          if (error_flag < 0) then
!!$                             write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                             stop
!!$                          endif
!!$                       enddo
!!$                       write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$                       !
!!$                       cmor_filename = ' '
!!$                       error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$                       if (error_flag < 0) then
!!$                          write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$                          stop
!!$                       else
!!$                          write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$                       endif
!!$                    enddo
!!$                    call close_cdf(myncid(ifile,ivar))
!!$                    if (allocated(time))      deallocate(time)
!!$                    if (allocated(time_bnds)) deallocate(time_bnds)
!!$                 enddo
!!$              enddo
        case ('msftmyz')
           !
           ! msftmyz: MOC
           !
           ! moc_comp(1)="Eulerian Mean" 
           ! moc_comp(2)="Eddy-Induced (bolus)" 
           ! moc_comp(3)="Submeso" 
           ! 
           ! transport_reg(1):="Global Ocean - Marginal Seas" 
           ! transport_reg(2):="Atlantic Ocean + Mediterranean Sea + Labrador Sea + GIN Sea + Arctic Ocean + Hudson Bay" 
           !
           ! MOC:coordinates = "lat_aux_grid moc_z moc_components transport_region time"
           !
           !                Y           Z      comp      basin
           allocate(indat4a(nlats_trans,nmoc_z,nmoc_comp,ntrans_reg))
           !
           ! basin 1: 'atlantic_arctic_ocean'
           ! basin 2: 'indian_pacific_ocean'
           ! basin 3: 'global_ocean'
           !
           !                  Y           Z      basin
           allocate(cmordat3d(nlats_trans,nmoc_z,3))
           do ivar = 1,xw(ixw)%ncesm_vars
              do ifile = 1,nc_nfiles(ivar)
                 call open_cdf(myncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                 call get_dims(myncid(ifile,ivar))
                 call get_vars(myncid(ifile,ivar))
                 !
                 if (.not.(allocated(time)))      allocate(time(ntimes(ifile,ivar)))
                 if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,ivar)))
                 !
                 write(*,*) 'tbnds loop: ',ivar,ifile,myncid(ifile,ivar)
                 do n = 1,ntimes(ifile,ivar)
                    time_counter = n
                    call read_var(myncid(ifile,ivar),'time_bound',time_bnds(:,n))
                 enddo
                 !
                 time_bnds(1,1) = int(time_bnds(1,1))-1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 select case (ntimes(ifile,ivar))
                 case ( 1140, 1872, 2388, 2400 )
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) =  1
                 case ( 1152 )
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = 13
                 end select
                 tidx2(nchunks(ifile)) = ntimes(ifile,ivar)
                 !
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       !
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat4a)
                       !
                       cmordat3d = spval
                       cmordat3d(:,:,2) = spval ! Indo-Pacific not supplied
                       !
                       ! Convert from Sv to kg s-1
                       !
                       where (indat4a(:,:,1,2) /= 0.)
                          cmordat3d(:,:,1) = indat4a(:,:,1,2)*1.e6*1.e3
                       elsewhere
                          cmordat3d(:,:,1) = spval
                       endwhere
                       where (indat4a(:,:,1,1) /= 0.)
                          cmordat3d(:,:,3) = indat4a(:,:,1,1)*1.e6*1.e3
                       elsewhere
                          cmordat3d(:,:,3) = spval
                       endwhere
                       !
                       tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                       error_flag = cmor_write(          &
                            var_id        = cmor_var_id, &
                            data          = cmordat3d,   &
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
                    if (ic < nchunks(ifile)) then
                       cmor_filename(1:) = ' '
                       error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                       if (error_flag < 0) then
                          write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                          stop
                       else
                          write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                       endif
                    endif
                 enddo
              enddo
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
           endif
        end select
        if (allocated(indat2a))   deallocate(indat2a)
        if (allocated(indat2b))   deallocate(indat2b)
        if (allocated(indat2c))   deallocate(indat2c)
        if (allocated(indat4a))   deallocate(indat4a)
        if (allocated(cmordat2d)) deallocate(cmordat2d)
        if (allocated(indat3a))   deallocate(indat3a)
        if (allocated(indat3b))   deallocate(indat3b)
        if (allocated(work3da))   deallocate(work3da)
        if (allocated(work3db))   deallocate(work3db)
        !
        ! Reset
        !
        error_flag   = 0
        mycmor%positive = ' '
        original_name= ' '
        call reset_netcdf_var
     endif
  enddo xwalk_loop
end program Omon_CMOR
