program day_CMOR
  ! Convert CCSM4 atm monthly (cam2.h0) data from single-field format
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
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,ixw,ilev,ic
  real::spval
  logical::does_exist
  !
  ! GO!
  !
  mycmor%table_file = 'day'
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
              !                    write(*,'(''OPENING: '',a80,'' myncid: '',i10)') trim(ncfile(ifile,ivar)),myncid(ifile,ivar)
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
              call close_cdf(myncid(ifile,ivar))
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
        call define_atm_axes(xw(ixw)%dims)
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
        case ('tauu','tauv','hfss','rlut','rlutcs','hfls','rlus','rsus','rsuscs','rsut','rsutcs')
           mycmor%positive = 'up'
        case ('rlds','rldscs','rsds','rsdscs','rsdt','rtmt')
           mycmor%positive = 'down'
        case ('clt','ci')
           var_info(var_found(1,1))%units = '1'
        case ('hurs','cl')
           var_info(var_found(1,1))%units = '%'
        case ('prc','pr','prsn')
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
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
        case ('ta','ua','va','hus','hur','wap','zg','tro3','tro3Clim','co2','co2Clim','ch4','ch4Clim','n2o','n2oClim')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3),axis_ids(4)/),  &
                missing_value=var_info(var_found(1,1))%missing_value,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        case ('clw','cli','cl')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3),axis_ids(4)/),  &
                missing_value=var_info(var_found(1,1))%missing_value,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        case default
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3)/),  &
                missing_value=var_info(var_found(1,1))%missing_value,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        end select
        write(*,'(''called cmor_variable; cmor_var_id:'',i8)') cmor_var_id
        !
        ! Perform derivations and cycle through time, writing data too
        !
        select case (xw(ixw)%entry)
        case ('ccb','cct','clivi','clwvi','evspsbl','hfls','hfss','hurs','huss',&
             'prw','psl','ps','rldscs','rlds','rlutcs','rlut','rsdscs','rsds','rsdt',&
             'sci','tas','tasmax','tasmin','tauu','tauv','ts')
           !
           ! No change
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
                 do n = 1,ntimes(ifile,ivar)
                    time_counter = n
                    call read_var(myncid(ifile,ivar),'time_bnds',time_bnds(:,n))
                 enddo
                 time_bnds(1,:) = time_bnds(2,:)
                 time_bnds(2,:) = time_bnds(1,:) + 1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 if (ntimes(ifile,ivar) == 56940) then         ! 20C from 1850-2005, use all times, 4 * 35y + 1 * 16y chunks
                    nchunks(ifile)= 5
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551, 38326, 51101/)      ! 1850, 1885, 1920, 1955, 1990
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 38325, 51100, 56940/)      ! 1884, 1919, 1954, 1989, 2005
                 endif
                 if (ntimes(ifile,ivar) == 35040) then         ! RCP from 2005-2100, use only 2006 onwards, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile)= 3
                    tidx1(1:nchunks(ifile)) = (/  366, 13141, 25916/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/13140, 25915, 35040/)      ! 2040, 2075, 2100
                 endif
                 if (ntimes(ifile,ivar) == 34675) then         ! RCP from 2006-2100, use all times, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile)= 3
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 34675/)      ! 2040, 2075, 2100
                 endif
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       !
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
                    if (ic < nchunks(ifile)) then
                       cmor_filename(1:) = ' '
                       error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                       if (error_flag < 0) then
                          write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                          stop                       else
                          write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                       endif
                    endif
                 enddo
              enddo
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR CMOR close of '',a)') cmor_filename(1:128)
              stop
           else
              write(*,'(''GOOD CMOR close of '',a)') cmor_filename(1:128)
           endif
        case ('pr','prsn')
           !
           ! pr  : Add PRECC + PRECL  , unit change from m s-1 to kg m-2 s-1
           ! prsn: Add PRECSC + PRECSL, unit change from m s-1 to kg m-2 s-1
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
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
                    call read_var(myncid(ifile,ivar),'time_bnds',time_bnds(:,n))
                 enddo
                 time_bnds(1,:) = time_bnds(2,:)
                 time_bnds(2,:) = time_bnds(1,:) + 1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 if (ntimes(ifile,ivar) == 56940) then         ! 20C from 1850-2005, use all times, 4 * 35y + 1 * 16y chunks
                    nchunks(ifile)= 5
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551, 38326, 51101/)      ! 1850, 1885, 1920, 1955, 1990
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 38325, 51100, 56940/)      ! 1884, 1919, 1954, 1989, 2005
                 endif
                 if (ntimes(ifile,ivar) == 35040) then         ! RCP from 2005-2100, use only 2006 onwards, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  366, 13141, 25916/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/13140, 25915, 35040/)      ! 2040, 2075, 2100
                 endif
                 if (ntimes(ifile,ivar) == 34675) then         ! RCP from 2006-2100, use all times, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 34675/)      ! 2040, 2075, 2100
                 endif
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat2a)
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat2b)
                       where ((indat2a /= spval).and.(indat2b /= spval))
                          cmordat2d = (indat2a + indat2b)*1000.
                       elsewhere
                          cmordat2d = spval
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
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
           endif
        case ('prc')
           !
           ! prc : PRECC, unit change from m s-1 to kg m-2 s-1
           !
           allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
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
                    call read_var(myncid(ifile,ivar),'time_bnds',time_bnds(:,n))
                 enddo
                 time_bnds(1,:) = time_bnds(2,:)
                 time_bnds(2,:) = time_bnds(1,:) + 1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 if (ntimes(ifile,ivar) == 56940) then         ! 20C from 1850-2005, use all times, 4 * 35y + 1 * 16y chunks
                    nchunks(ifile) = 5
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551, 38326, 51101/)      ! 1850, 1885, 1920, 1955, 1990
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 38325, 51100, 56940/)      ! 1884, 1919, 1954, 1989, 2005
                 endif
                 if (ntimes(ifile,ivar) == 35040) then         ! RCP from 2005-2100, use only 2006 onwards, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  366, 13141, 25916/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/13140, 25915, 35040/)      ! 2040, 2075, 2100
                 endif
                 if (ntimes(ifile,ivar) == 34675) then         ! RCP from 2006-2100, use all times, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 34675/)      ! 2040, 2075, 2100
                 endif
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat2a)
                       where (indat2a /= spval)
                          cmordat2d = indat2a*1000.
                       elsewhere
                          cmordat2d = spval
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
           enddo
        case ('rlus')
           !
           ! rlus: Add FLDS + FLNS
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
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
                    call read_var(myncid(ifile,ivar),'time_bnds',time_bnds(:,n))
                 enddo
                 time_bnds(1,:) = time_bnds(2,:)
                 time_bnds(2,:) = time_bnds(1,:) + 1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 if (ntimes(ifile,ivar) == 56940) then         ! 20C from 1850-2005, use all times, 4 * 35y + 1 * 16y chunks
                    nchunks(ifile) = 5
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551, 38326, 51101/)      ! 1850, 1885, 1920, 1955, 1990
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 38325, 51100, 56940/)      ! 1884, 1919, 1954, 1989, 2005
                 endif
                 if (ntimes(ifile,ivar) == 35040) then         ! RCP from 2005-2100, use only 2006 onwards, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  366, 13141, 25916/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/13140, 25915, 35040/)      ! 2040, 2075, 2100
                 endif
                 if (ntimes(ifile,ivar) == 34675) then         ! RCP from 2006-2100, use all times, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 34675/)      ! 2040, 2075, 2100
                 endif
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat2a)
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat2b)
                       where ((indat2a /= spval).and.(indat2b /= spval))
                          cmordat2d = (indat2a + indat2b)
                       elsewhere
                          cmordat2d = spval
                       endwhere
                       tval(1)   = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                       error_flag = cmor_write(      &
                            var_id        = cmor_var_id, &
                            data          = cmordat2d, &
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
           enddo
        case ('rsus','rsuscs','rsut','rsutcs','rtmt')
           !
           ! rsus   : FSDS  - FSNS
           ! rsuscs : FSDSC - FSNSC
           ! rsut   : SOLIN - FSNTOA
           ! rsutcs : SOLIN - FSNTOAC
           ! rtmt   : FSNT  - FLNT
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
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
                    call read_var(myncid(ifile,ivar),'time_bnds',time_bnds(:,n))
                 enddo
                 time_bnds(1,:) = time_bnds(2,:)
                 time_bnds(2,:) = time_bnds(1,:) + 1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 if (ntimes(ifile,ivar) == 56940) then         ! 20C from 1850-2005, use all times, 4 * 35y + 1 * 16y chunks
                    nchunks(ifile) = 5
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551, 38326, 51101/)      ! 1850, 1885, 1920, 1955, 1990
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 38325, 51100, 56940/)      ! 1884, 1919, 1954, 1989, 2005
                 endif
                 if (ntimes(ifile,ivar) == 35040) then         ! RCP from 2005-2100, use only 2006 onwards, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  366, 13141, 25916/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/13140, 25915, 35040/)      ! 2040, 2075, 2100
                 endif
                 if (ntimes(ifile,ivar) == 34675) then         ! RCP from 2006-2100, use all times, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 34675/)      ! 2040, 2075, 2100
                 endif
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat2a)
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat2b)
                       where ((indat2a /= spval).and.(indat2b /= spval))
                          cmordat2d = indat2a - indat2b
                       elsewhere
                          cmordat2d = spval
                       endwhere
                       tval(1)   = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                       error_flag = cmor_write(      &
                            var_id        = cmor_var_id, &
                            data          = cmordat2d, &
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
           enddo
        case ('ta','ua','va','hus','hur','wap','zg')
           !
           ! Vertically interpolate to standard pressure levels
           !
           allocate(indat3a(nlons,nlats,nlevs),cmordat3d(nlons,nlats,nplev8))
           allocate(psdata(nlons,nlats))
           write(*,*) 'OPENING: ',trim(ncfile(1,1)),' and ',trim(ncfile(1,2))
           if (nc_nfiles(1) == nc_nfiles(2)) then
              do ifile = 1,nc_nfiles(1)
                 call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
                 call get_dims(myncid(ifile,1))
                 call get_vars(myncid(ifile,1))
                 call open_cdf(myncid(ifile,2),trim(ncfile(ifile,2)),.true.)
                 call get_dims(myncid(ifile,2))
                 call get_vars(myncid(ifile,2))
                 !
                 if (.not.(allocated(time)))      allocate(time(ntimes(ifile,1)))
                 if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,1)))
                 !
                 do n = 1,ntimes(ifile,1)
                    time_counter = n
                    call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
                 enddo
                 time_bnds(1,:) = time_bnds(2,:)
                 time_bnds(2,:) = time_bnds(1,:) + 1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 ! Determine amount of data to write, to keep close to ~2 GB limit
                 !
                 select case (ntimes(ifile,1))
                 case ( 1824, 1825, 2190 )
                    nchunks(ifile) = 1
                    tidx1(1:nchunks(ifile)) = 1
                    tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
                 case ( 6012 )
                    nchunks(ifile) = 10
                    tidx1(1) =   1
                    tidx2(1) = 600
                    do ic = 2,nchunks(ifile)
                       tidx1(ic) = tidx2(ic-1) + 1
                       tidx2(ic) = tidx1(ic) + 599
                    enddo
                 end select
                 tidx2(nchunks(ifile)) = ntimes(ifile,1)
                 write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
                 write(*,*) 'GB size cmordat3d: ',(4*ntimes(ifile,1)*size(cmordat3d))/1.e9
                 write(*,*) 'var_info: ',trim(var_info(var_found(ifile,1))%name)
                 write(*,*) 'var_info: ',trim(var_info(var_found(ifile,2))%name)
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                       call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,psdata)
                       !
                       ! Convert PS to mb from Pa
                       !
                       psdata = psdata * 0.01
                       !
                       cmordat3d = spval
                       !
                       ! Do vertical interpolation to pressure levels
                       !
                       call vertint(indat3a,cmordat3d,atm_levs,atm_plev8*0.01,psdata,spval,nlons,nlats,nlevs,nlevs+1,nplev8)
                       !
                       tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                       error_flag = cmor_write(        &
                            var_id        = cmor_var_id,   &
                            data          = cmordat3d, &
                            ntimes_passed = 1,         &
                            time_vals     = tval,      &
                            time_bnds     = tbnd)
                       if (error_flag < 0) then
                          write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                          stop
                       endif
                    enddo
                    write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                    !
                    if (ic .le. nchunks(ifile)) then
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
                 if (allocated(time))      deallocate(time)
                 if (allocated(time_bnds)) deallocate(time_bnds)
              enddo
              !
              error_flag = cmor_close()
              if (error_flag < 0) then
                 write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
              else
                 write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
              endif
           endif
        case ('clw','cli','cl')
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nlevs),indat2a(nlons,nlats))
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
                    call read_var(myncid(ifile,ivar),'time_bnds',time_bnds(:,n))
                 enddo
                 time_bnds(1,:) = time_bnds(2,:)
                 time_bnds(2,:) = time_bnds(1,:) + 1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 ! Determine amount of data to write, to keep close to ~2 GB limit
                 !
                 if (ntimes(ifile,ivar) == 56940) then         ! 20C from 1850-2005, use all times, 4 * 35y + 1 * 16y chunks
                    nchunks(ifile) = 5
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551, 38326, 51101/)      ! 1850, 1885, 1920, 1955, 1990
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 38325, 51100, 56940/)      ! 1884, 1919, 1954, 1989, 2005
                 endif
                 if (ntimes(ifile,ivar) == 35040) then         ! RCP from 2005-2100, use only 2006 onwards, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  366, 13141, 25916/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/13140, 25915, 35040/)      ! 2040, 2075, 2100
                 endif
                 if (ntimes(ifile,ivar) == 34675) then         ! RCP from 2006-2100, use all times, 2 * 35y + 1 * 25y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/    1, 12776, 25551/)      ! 2006, 2041, 2076
                    tidx2(1:nchunks(ifile)) = (/12775, 25550, 34675/)      ! 2040, 2075, 2100
                 endif
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat3a)
                       call read_var(myncid(ifile,ivar),var_info(var_found(ifile,ivar))%name,indat2a)
                       where (indat3a > 5000.) indat3a = spval
                       tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                       error_flag = cmor_write(        &
                            var_id        = cmor_var_id,   &
                            data          = indat3a,   &
                            ntimes_passed = 1,         &
                            time_vals     = tval,      &
                            time_bnds     = tbnd)
                       if (error_flag < 0) then
                          write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                          stop
                       endif
                       error_flag = cmor_write(        &
                            var_id        = zfactor_id,&
                            data          = indat2a,   &
                            ntimes_passed = 1,         &
                            time_vals     = tval,      &
                            time_bnds     = tbnd,      &
                            store_with    = cmor_var_id)
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
                    write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                 enddo
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
        do ivar = 1,xw(ixw)%ncesm_vars
           do ifile = 1,nc_nfiles(ivar)
              call close_cdf(myncid(ifile,ivar))
           enddo
        enddo
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
!!$              error_flag = cmor_close()
!!$              if (error_flag < 0) then
!!$                 write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
!!$              else
!!$                 write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
!!$              endif
        call reset_netcdf_var
     endif
  enddo xwalk_loop
end program day_CMOR
