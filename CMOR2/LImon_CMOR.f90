program LImon_CMOR
  ! Convert CESM data from single-field format to CMOR-compliant format
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
  real,dimension(:,:)  ,allocatable::indat2a,indat2b,indat2c,cmordat2d,psdata
  real,dimension(:,:,:),allocatable::indat3a,indat3b,indat3c,cmordat3d,work3da,work3db
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::time_bnds
  double precision,dimension(1)  ::tval
  double precision,dimension(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile,cmor_filename
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,itab,ixw,ilev,ic,jfile
  logical::does_exist
  !
  ! GO!
  !
  mycmor%table_file = 'LImon'
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
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
  xwalk_file = 'xwalk_'//trim(exp(exp_found)%cmip)//'_'//trim(mycmor%table_file)
  call load_xwalk(xwalk_file)
  !
  ! Get table information
  !
  mycmor%table_file = 'Tables/'//trim(exp(exp_found)%cmip)//'_'//trim(mycmor%table_file)
  inquire(file=mycmor%table_file,exist=does_exist)
  if (.not.(does_exist)) then
     write(*,*) 'Cannot find ',trim(mycmor%table_file),'. Dying.'
     stop
  endif
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
              write(*,'(''time length FROM: '',a,'' myncid: '',i10,'' NT: '',i10)') trim(ncfile(ifile,ivar)),myncid(ifile,ivar),ntimes(ifile,ivar)
              do n=1,var_counter
                 if (trim(var_info(n)%name) == trim(xw(ixw)%cesm_vars(ivar))) then
                    var_found(ifile,ivar) = n
                 endif
              enddo
              if (var_found(ifile,ivar) == 0) then
                 !
                 ! Never found - quit
                 !
                 write(*,'(''NEVER FOUND: '',a,'' STOP. '')') trim(xw(ixw)%cesm_vars(ivar))
                 stop
              endif
           enddo
        enddo
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
           write(*,*) 'ERROR on cmor_dataset!'
           write(*,*) 'outpath               = ',trim(mycmor%outpath)
           write(*,*) 'experiment_id         = ',trim(mycmor%experiment_id)
           write(*,*) 'institution           = ',trim(mycmor%institution)
           write(*,*) 'source                = ',trim(mycmor%source)
           write(*,*) 'calendar              = ',trim(mycmor%calendar)
           write(*,*) 'realization           = ',mycmor%realization
           write(*,*) 'initialization_method = ',mycmor%initialization_method
           write(*,*) 'physics_version       = ',mycmor%physics_version
           write(*,*) 'contact               = ',trim(mycmor%contact)
           write(*,*) 'history               = ',trim(mycmor%history)
           write(*,*) 'comment               = ',trim(mycmor%comment)
           write(*,*) 'references            = ',trim(mycmor%references)
           write(*,*) 'model_id              = ',trim(mycmor%model_id)
           write(*,*) 'forcing               = ',trim(mycmor%forcing)
           write(*,*) 'institute_id          = ',trim(mycmor%institute_id)
           write(*,*) 'parent_experiment_id  = ',trim(mycmor%parent_experiment_id)
           write(*,*) 'parent_experiment_rip = ',trim(mycmor%parent_experiment_rip)
           write(*,*) 'branch_time           = ',mycmor%branch_time
           stop
        endif
        !
        ! Add global metadata
        !
        call add_global_metadata
        !
        ! Define axes via 'cmor_axis'
        !
        call define_lnd_axes(xw(ixw)%dims)
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
        if (xw(ixw)%ncesm_vars == 3) write(original_name,'(a,'','',a,'','',a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        !
        ! Modify units as necessary to accomodate udunits' inability to convert 
        !
        select case (xw(ixw)%entry)
        case ('hfdsn')
           mycmor%positive = 'up'
        case ('snw')
           var_info(var_found(1,1))%units = 'kg m-2'
        case ('snc')
           var_info(var_found(1,1))%units = '%'
        case ('snm')
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
        end select
        !
        write(*,*) 'calling cmor_variable:'
        write(*,*) 'table         = ',trim(mycmor%table_file)
        write(*,*) 'table_entry   = ',trim(xw(ixw)%entry)
        write(*,*) 'dimensions    = ',trim(xw(ixw)%dims)
        write(*,*) 'units         = ',trim(var_info(var_found(1,1))%units)
        write(*,*) 'axis_ids      = ',axis_ids(1:4)
        write(*,*) 'missing_value = ',var_info(var_found(1,1))%missing_value
        write(*,*) 'positive      = ',trim(mycmor%positive)
        write(*,*) 'original_name = ',trim(original_name)
        !
        cmor_var_id = cmor_variable(                            &
             table=mycmor%table_file,                           &
             table_entry=xw(ixw)%entry,                         &
             units=var_info(var_found(1,1))%units,                &
             axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3)/),  &
             missing_value=var_info(var_found(1,1))%missing_value,&
             positive=mycmor%positive,                          &
             original_name=original_name,                       &
             comment=xw(ixw)%comment)
        write(*,'(''called cmor_variable; cmor_var_id:'',i8)') cmor_var_id
        !
        ! Perform derivations and cycle through time, writing data too
        !
        select case (xw(ixw)%entry)
        case ('snd','lwsnl','hfdsn','snw')
           !
           ! snd  : SNOWDP, no change
           ! lwsnl: SNOWLIQ, no numerical change except units
           ! snw  : H2OSNO, no numerical change except units
           !
           allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
           do ifile = 1,nc_nfiles(1)
              !
              if (allocated(time)) deallocate(time)
              allocate(time(ntimes(ifile,1)))
              if (allocated(time_bnds)) deallocate(time_bnds)
              allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n=1,ntimes(1,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bounds',time_bnds(:,n))
!                 if (n == 1) time_bnds(1,n) = 0.
              enddo
              time = (time_bnds(1,:)+time_bnds(2,:))/2.
              !
              select case (ntimes(ifile,1))
              case ( 1872,1140,3612,6012,12012 )  ! All data
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case ( 1152 )  ! RCP, 2005-2100, skip 2006
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 13
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 3613
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(1,1)
              end select
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat2a)
                    cmordat2d = 0.
                    where (indat2a /= var_info(var_found(ifile,1))%missing_value)
                       cmordat2d = indat2a
                    elsewhere
                       cmordat2d = 0.
                    endwhere
                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(      &
                         var_id        = cmor_var_id, &
                         data          = cmordat2d,   &
                         ntimes_passed = 1,       &
                         time_vals     = tval,    &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                       stop
                    endif
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                 if (ic <= nchunks(ifile)) then
                    write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                    cmor_filename = ' '
                    error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                    if (error_flag < 0) then
                       write(*,'(''ERROR close chunk: '',i1,'' of '',a)') ic,trim(cmor_filename)
                       stop
                    else
                       write(*,'(''GOOD close chunk : '',i2,'' of '',a)') ic,trim(cmor_filename)
                    endif
                 else
                    error_flag = cmor_close()
                    if (error_flag < 0) then
                       write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
                    else
                       write(*,'(''GOOD cmor_close of  : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
                    endif
                 endif
              enddo
           enddo
        case ('snm')
           !
           ! snm  : QMELT, spval set to zero
           !
           allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
           do ifile = 1,nc_nfiles(1)
              call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
              call get_dims(myncid(ifile,1))
              call get_vars(myncid(ifile,1))
              !
              if (.not.(allocated(time)))      allocate(time(ntimes(ifile,1)))
              if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n=1,ntimes(ifile,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bounds',time_bnds(:,n))
                 if (n == 1) time_bnds(1,n) = 0.
              enddo
              time = (time_bnds(1,:)+time_bnds(2,:))/2.
              !
              select case (ntimes(ifile,1))
              case ( 1872,1140,3612,6012,12012 )  ! All data
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case ( 1152 )  ! RCP, 2005-2100, skip 2006
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 13
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 3613
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(1,1)
              end select
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                    cmordat2d = 0.
                    where (indat2a /= var_info(var_found(1,1))%missing_value)
                       cmordat2d = indat2a
                    elsewhere
                       cmordat2d = 0.
                    endwhere
                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(      &
                         var_id        = cmor_var_id, &
                         data          = cmordat2d,   &
                         ntimes_passed = 1,       &
                         time_vals     = tval,    &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                       stop
                    endif
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
              enddo
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
           endif
           do ifile = 1,nc_nfiles(1)
              call close_cdf(myncid(ifile,1))
           enddo
        case ('snc')
           !
           ! snc : FSNO, mulitply by 100 to get percentage
           !
           allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
           do ifile = 1,nc_nfiles(1)
              call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
              call get_dims(myncid(ifile,1))
              call get_vars(myncid(ifile,1))
              !
              if (.not.(allocated(time)))      allocate(time(ntimes(ifile,1)))
              if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n=1,ntimes(1,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bounds',time_bnds(:,n))
                 if (n == 1) time_bnds(1,n) = 0.
              enddo
              time = (time_bnds(1,:)+time_bnds(2,:))/2.
              !
              select case (ntimes(ifile,1))
              case ( 1872,1140,3612,6012,12012 )  ! All data
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case ( 1152 )  ! RCP, 2005-2100, skip 2006
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 13
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 3613
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(1,1)
              end select
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                    where (indat2a > (0.1*var_info(var_found(1,1))%missing_value))
                       cmordat2d = var_info(var_found(1,1))%missing_value
                    elsewhere
                       cmordat2d = 100*indat2a
                    endwhere
                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(      &
                         var_id        = cmor_var_id, &
                         data          = cmordat2d,   &
                         ntimes_passed = 1,       &
                         time_vals     = tval,    &
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
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
           endif
           do ifile = 1,nc_nfiles(1)
              call close_cdf(myncid(ifile,1))
           enddo
        case ('sootsn')
           !
           ! sootsn : SNOBCMSL + SNODSTMSL + SNOOCMSL
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
           do ifile = 1,nc_nfiles(1)
              call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
              call get_dims(myncid(ifile,1))
              call get_vars(myncid(ifile,1))
              call open_cdf(myncid(ifile,2),trim(ncfile(ifile,2)),.true.)
              call get_dims(myncid(ifile,2))
              call get_vars(myncid(ifile,2))
              call open_cdf(myncid(ifile,3),trim(ncfile(ifile,3)),.true.)
              call get_dims(myncid(ifile,3))
              call get_vars(myncid(ifile,3))
              !
              if (.not.(allocated(time)))      allocate(time(ntimes(ifile,1)))
              if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n=1,ntimes(1,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bounds',time_bnds(:,n))
                 if (n == 1) time_bnds(1,n) = 0.
              enddo
              time = (time_bnds(1,:)+time_bnds(2,:))/2.
              !
              select case (ntimes(ifile,1))
              case ( 1872,1140,3612,6012,12012 )  ! All data
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case ( 1152 )  ! RCP, 2005-2100, skip 2006
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 13
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 3613
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(1,1)
              end select
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat2a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat2b)
                    call read_var(myncid(ifile,3),var_info(var_found(ifile,3))%name,indat2c)
                    where ((indat2a /= var_info(var_found(ifile,1))%missing_value).and.&
                           (indat2b /= var_info(var_found(ifile,2))%missing_value).and.&
                           (indat2c /= var_info(var_found(ifile,3))%missing_value))
                    cmordat2d = indat2a+indat2b+indat2c
                    elsewhere
                    cmordat2d = var_info(var_found(ifile,1))%missing_value
                    endwhere
                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(      &
                         var_id        = cmor_var_id, &
                         data          = cmordat2d,   &
                         ntimes_passed = 1,       &
                         time_vals     = tval,    &
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
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
           endif
           do ifile = 1,nc_nfiles(1)
              call close_cdf(myncid(ifile,1))
              call close_cdf(myncid(ifile,2))
              call close_cdf(myncid(ifile,3))
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
        call reset_netcdf_var
     endif
  enddo xwalk_loop
end program LImon_CMOR
