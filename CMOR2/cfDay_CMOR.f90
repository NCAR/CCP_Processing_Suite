program cfDay_CMOR
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
  real,dimension(:,:)    ,allocatable::indat2a,indat2b,indat2c,cmordat2d,psdata
  real,dimension(:,:,:)  ,allocatable::indat3a,indat3b,indat3c,cmordat3d,work3da,work3db,parasol
  real,dimension(:,:,:,:),allocatable::indat4a
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
  mycmor%table_file = 'cfDay'
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
        call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
        call get_dims(myncid(1,1))
        call get_vars(myncid(1,1))
        call read_att_text(myncid(1,1),'time','units',time_units)
        write(*,'(''time units in: '',i10,5x,a)') myncid(1,1),trim(time_units)
        call close_cdf(myncid(1,1))
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
!              call close_cdf(myncid(ifile,ivar))
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
        case ('rhs')
           var_info(var_found(1,1))%units = '%'
        case ('prc','pr','prsn')
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
        case ('snw')
           var_info(var_found(1,1))%units = 'kg m-2'
        case ('snc')
           var_info(var_found(1,1))%units = '%'
        case ('parasolRefl')
           var_info(var_found(1,1))%units = '1'
        case ('mrro','mrros')
           ! mm s-1 is the same as kg m-2 s-1
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
        case ('ccb','cct')
           var_info(var_found(1,1))%units = 'Pa'
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
        case ('clw','cli','cl','ta','ua','va','hus','hur','wap','zg','mc',&
              'clcalipso','parasolRefl')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(2),axis_ids(3),axis_ids(4),axis_ids(1)/),  &
                missing_value=var_info(var_found(1,1))%missing_value,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        case ('clisccp')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3),axis_ids(4),axis_ids(5)/),&
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
        case ('rhs','clivi','clwvi','evspsbl','hfls','hfss','hurs','huss',&
              'prw','psl','ps','rldscs','rlds','rlutcs','rlut','rsdscs','rsds','rsdt',&
              'sci','tas','tasmax','tasmin','tauu','tauv','ts','clt',&
              'albisccp','ccb','cct','clhcalipso','cllcalipso','clmcalipso',&
              'cltcalipso','cltisccp','pctisccp','ta700')
           !
           ! No change
           !
           allocate(indat2a(nlons,nlats))
           do ifile = 1,nc_nfiles(1)
              call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)             
              write(*,*) 'open_cdf: ',myncid(ifile,1),trim(ncfile(ifile,1))
              call get_dims(myncid(ifile,1))
              call get_vars(myncid(ifile,1))
              !
              if (allocated(time))      deallocate(time)
              if (allocated(time_bnds)) deallocate(time_bnds)
              allocate(time(ntimes(ifile,1)))
              allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n = 1,ntimes(ifile,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
              enddo
              time_bnds(1,:) = time_bnds(2,:)
              time_bnds(2,:) = time_bnds(1,:) + 1
              time = (time_bnds(1,:)+time_bnds(2,:))/2.
              write(*,'(''time length FROM: '',a,'' myncid: '',i10,'' NT: '',i10)') trim(ncfile(ifile,1)),myncid(ifile,1),ntimes(ifile,1)
              !
              select case (ntimes(ifile,1))
              case ( 1825 ) ! "e" series from 2004-2008
                 if (trim(case_read)=='f40.amip_4k_cosp.cam4.1deg.001e') then ! Use only 2006-2008
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = (/ 731/)
                    tidx2(1:nchunks(ifile)) = (/1825/)
                 else
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = (/   1/)
                    tidx2(1:nchunks(ifile)) = (/1825/)
                 endif
              case default
                 nchunks(ifile)= 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    !
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat2a)
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
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                    stop                       
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                 endif
              enddo
           enddo
           !
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           endif
!           do ifile = 1,nc_nfiles(1)
!              call close_cdf(myncid(ifile,1))
!           enddo
        case ('pr','prsn')
           !
           ! pr  : Add PRECC + PRECL  , unit change from m s-1 to kg m-2 s-1
           ! prsn: Add PRECSC + PRECSL, unit change from m s-1 to kg m-2 s-1
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
           if (nc_nfiles(1) == nc_nfiles(2)) then
              do ifile = 1,nc_nfiles(1)
                 call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
                 call get_dims(myncid(ifile,1))
                 call get_vars(myncid(ifile,1))
                 call open_cdf(myncid(ifile,2),trim(ncfile(ifile,2)),.true.)
                 call get_dims(myncid(ifile,2))
                 call get_vars(myncid(ifile,2))
                 !
                 if (allocated(time))      deallocate(time)
                 if (allocated(time_bnds)) deallocate(time_bnds)
                 allocate(time(ntimes(ifile,1)))
                 allocate(time_bnds(2,ntimes(ifile,1)))
                 !
                 do n = 1,ntimes(ifile,1)
                    time_counter = n
                    call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
                 enddo
                 time_bnds(1,:) = time_bnds(2,:)
                 time_bnds(2,:) = time_bnds(1,:) + 1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 write(*,'(''time length FROM: '',a,'' myncid: '',i10,'' NT: '',i10)') trim(ncfile(ifile,1)),myncid(ifile,1),ntimes(ifile,1)
                 !
                 select case (ntimes(ifile,1))
                 case ( 1825 ) ! "e" series; use only 2008 - maybe
                    if (trim(case_read)=='f40.amip_4k_cosp.cam4.1deg.001e') then ! Use only 2006-2008
                       nchunks(ifile)= 1
                       tidx1(1:nchunks(ifile)) = (/ 731/)
                       tidx2(1:nchunks(ifile)) = (/1825/)
                    else
                       nchunks(ifile)= 1
                       tidx1(1:nchunks(ifile)) = (/   1/)
                       tidx2(1:nchunks(ifile)) = (/1825/)
                    endif
                 case default
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = 1
                    tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
                 end select
                 write(*,'(''# chunks '',i3,'':'',20((i6,''-'',i6),'',''))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat2a)
                       call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat2b)
                       cmordat2d = (indat2a + indat2b)*1000.
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
                    cmor_filename(1:) = ' '
                    error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                    if (error_flag < 0) then
                       write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                       stop
                    else
                       write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                    endif
                 enddo
              enddo
           endif
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           endif
!           do ifile = 1,nc_nfiles(1)
 !             call close_cdf(myncid(ifile,1))
 !          enddo
        case ('prc')
           !
           ! prc : PRECC, unit change from m s-1 to kg m-2 s-1
           !
           allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
           do ifile = 1,nc_nfiles(1)
              call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
              call get_dims(myncid(ifile,1))
              call get_vars(myncid(ifile,1))
              !
              !
              if (allocated(time))      deallocate(time)
              if (allocated(time_bnds)) deallocate(time_bnds)
              allocate(time(ntimes(ifile,1)))
              allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n = 1,ntimes(ifile,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
              enddo
              time_bnds(1,:) = time_bnds(2,:)
              time_bnds(2,:) = time_bnds(1,:) + 1
              time = (time_bnds(1,:)+time_bnds(2,:))/2.
              !
              select case (ntimes(ifile,1))
              case ( 1825 ) ! "e" series; use only 2008 - maybe
                 if (trim(case_read)=='f40.amip_4k_cosp.cam4.1deg.001e') then ! Use only 2006-2008
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = (/ 731/)
                    tidx2(1:nchunks(ifile)) = (/1825/)
                 else
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = (/   1/)
                    tidx2(1:nchunks(ifile)) = (/1825/)
                 endif
              case default
                 nchunks(ifile)= 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',20((i6,''-'',i6),'',''))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat2a)
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
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                 endif
              enddo
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR CMOR close of '',a)') trim(cmor_filename(1:))
              stop
           else
              write(*,'(''GOOD CMOR close of '',a)') trim(cmor_filename(1:))
           endif
!           do ifile = 1,nc_nfiles(1)
 !             call close_cdf(myncid(ifile,1))
 !          enddo
        case ('rlus')
           !
           ! rlus: Add FLDS + FLNS
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
           if (nc_nfiles(1) == nc_nfiles(2)) then
              do ifile = 1,nc_nfiles(1)
                 call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
                 call get_dims(myncid(ifile,1))
                 call get_vars(myncid(ifile,1))
                 call open_cdf(myncid(ifile,2),trim(ncfile(ifile,2)),.true.)
                 call get_dims(myncid(ifile,2))
                 call get_vars(myncid(ifile,2))
                 !
                 if (allocated(time))      deallocate(time)
                 if (allocated(time_bnds)) deallocate(time_bnds)
                 allocate(time(ntimes(ifile,1)))
                 allocate(time_bnds(2,ntimes(ifile,1)))
                 !
                 do n = 1,ntimes(ifile,1)
                    time_counter = n
                    call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
                 enddo
                 time_bnds(1,:) = time_bnds(2,:)
                 time_bnds(2,:) = time_bnds(1,:) + 1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 select case (ntimes(ifile,1))
                 case ( 1825 ) ! "e" series; use only 2008 - maybe
                    if (trim(case_read)=='f40.amip_4k_cosp.cam4.1deg.001e') then ! Use only 2006-2008
                       nchunks(ifile)= 1
                       tidx1(1:nchunks(ifile)) = (/ 731/)
                       tidx2(1:nchunks(ifile)) = (/1825/)
                    else
                       nchunks(ifile)= 1
                       tidx1(1:nchunks(ifile)) = (/   1/)
                       tidx2(1:nchunks(ifile)) = (/1825/)
                    endif
                 case default
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = 1
                    tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
                 end select
                 write(*,'(''# chunks '',i3,'':'',20((i6,''-'',i6),'',''))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat2a)
                       call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat2b)
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
                    if (nchunks(ifile) /= 1) then
                       cmor_filename(1:) = ' '
                       error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                       if (error_flag < 0) then
                          write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                          stop
                       else
                          write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                       endif
                    endif
                 enddo
              enddo
           endif
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           endif
!           do ifile = 1,nc_nfiles(1)
 !             call close_cdf(myncid(ifile,1))
 !          enddo
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
           if (nc_nfiles(1) == nc_nfiles(2)) then
              do ifile = 1,nc_nfiles(1)
                 call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
                 call get_dims(myncid(ifile,1))
                 call get_vars(myncid(ifile,1))
                 call open_cdf(myncid(ifile,2),trim(ncfile(ifile,2)),.true.)
                 call get_dims(myncid(ifile,2))
                 call get_vars(myncid(ifile,2))
                 !
                 if (allocated(time))      deallocate(time)
                 if (allocated(time_bnds)) deallocate(time_bnds)
                 allocate(time(ntimes(ifile,1)))
                 allocate(time_bnds(2,ntimes(ifile,1)))
                 !
                 do n = 1,ntimes(ifile,1)
                    time_counter = n
                    call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
                 enddo
                 time_bnds(1,:) = time_bnds(2,:)
                 time_bnds(2,:) = time_bnds(1,:) + 1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 !
                 select case (ntimes(ifile,1))
                 case ( 1825 ) ! "e" series; use only 2008 - maybe
                    if (trim(case_read)=='f40.amip_4k_cosp.cam4.1deg.001e') then ! Use only 2006-2008
                       nchunks(ifile)= 1
                       tidx1(1:nchunks(ifile)) = (/ 731/)
                       tidx2(1:nchunks(ifile)) = (/1825/)
                    else
                       nchunks(ifile)= 1
                       tidx1(1:nchunks(ifile)) = (/   1/)
                       tidx2(1:nchunks(ifile)) = (/1825/)
                    endif
                 case default
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = 1
                    tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
                 end select
                 write(*,'(''# chunks '',i3,'':'',20((i6,''-'',i6),'',''))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat2a)
                       call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat2b)
                       where ((indat2a /= spval).and.(indat2b /= spval))
                          cmordat2d = indat2a - indat2b
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
                    if (nchunks(ifile) /= 1) then
                       cmor_filename(1:) = ' '
                       error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                       if (error_flag < 0) then
                          write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                          stop
                       else
                          write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                       endif
                    endif
                 enddo
              enddo
           endif
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           endif
!           do ifile = 1,nc_nfiles(1)
 !             call close_cdf(myncid(ifile,1))
 !          enddo
        case ('cl','cli','clw','hur','hus','ta','ua','va','wap','zg')
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nlevs),indat2a(nlons,nlats))
           if (nc_nfiles(1) == nc_nfiles(2)) then
              do ifile = 1,nc_nfiles(1)
                 call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
                 call get_dims(myncid(ifile,1))
                 call get_vars(myncid(ifile,1))
                 call open_cdf(myncid(ifile,2),trim(ncfile(ifile,2)),.true.)
                 call get_dims(myncid(ifile,2))
                 call get_vars(myncid(ifile,2))
                 !
                 if (allocated(time))      deallocate(time)
                 if (allocated(time_bnds)) deallocate(time_bnds)
                 allocate(time(ntimes(ifile,1)))
                 allocate(time_bnds(2,ntimes(ifile,1)))
                 !
                 do n = 1,ntimes(ifile,1)
                    time_counter = n
                    call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
                 enddo
                 time_bnds(1,:) = time_bnds(2,:)
                 time_bnds(2,:) = time_bnds(1,:) + 1
                 time = (time_bnds(1,:)+time_bnds(2,:))/2.
                 write(*,'(''time length FROM: '',a,'' myncid: '',i10,'' NT: '',i10)') trim(ncfile(ifile,1)),myncid(ifile,1),ntimes(ifile,1)
                 !
                 ! Determine amount of data to write, to keep close to ~2 4B limit
                 !
                 select case (ntimes(ifile,1))
                 case ( 9855 )         ! 1979-2005, 1 file per year
                    nchunks(ifile) = 27
                    tidx1(1) =   1
                    tidx2(1) = 365
                    do ic = 2,nchunks(ifile)
                       tidx1(ic) = tidx2(ic-1) + 1
                       tidx2(ic) = tidx1(ic) + 364
                    enddo
                    tidx2(nchunks(ifile)) = ntimes(ifile,1)
                 case ( 1825 ) ! "e" series; use only 2006-2008; case-specific
                    if (trim(case_read)=='f40.amip_4k_cosp.cam4.1deg.001e') then ! Use only 2006-2008
                       nchunks(ifile)= 3
                       tidx1(1:nchunks(ifile)) = (/ 731,1096,1461/)
                       tidx2(1:nchunks(ifile)) = (/1095,1460,1825/)
                    else
                       nchunks(ifile)= 1
                       tidx1(1:nchunks(ifile)) = (/   1/)
                       tidx2(1:nchunks(ifile)) = (/1825/)
                    endif
                 case default
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = 1
                    tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
                 end select
                 write(*,'(''# chunks '',i3,'':'',20((i6,''-'',i6),'',''))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                       call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat2a)
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
                    cmor_filename(1:) = ' '
                    error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                    if (error_flag < 0) then
                       write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                       stop
                    else
                       write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                    endif
                 enddo
              enddo
           endif
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR CMOR close of '',a)') trim(cmor_filename(1:))
              stop
           else
              write(*,'(''GOOD CMOR close of '',a)') trim(cmor_filename(1:))
           endif
!           do ifile = 1,nc_nfiles(1)
 !             call close_cdf(myncid(ifile,1))
 !          enddo
        case ('clcalipso')
           !
           ! clcalipso
           !
           allocate(indat3a(nlons,nlats,ncosp_ht))
           do ifile = 1,nc_nfiles(1)
              call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
              call get_dims(myncid(ifile,1))
              call get_vars(myncid(ifile,1))
              !
              if (allocated(time))      deallocate(time)
              if (allocated(time_bnds)) deallocate(time_bnds)
              allocate(time(ntimes(ifile,1)))
              allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n = 1,ntimes(ifile,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
              enddo
              time_bnds(1,:) = time_bnds(2,:)
              time_bnds(2,:) = time_bnds(1,:) + 1
              time = (time_bnds(1,:)+time_bnds(2,:))/2.
              write(*,'(''time length FROM: '',a,'' myncid: '',i10,'' NT: '',i10)') trim(ncfile(ifile,1)),myncid(ifile,1),ntimes(ifile,1)
              !
              ! Determine amount of data to write, to keep close to ~2 4B limit
              !
              select case (ntimes(ifile,1))
              case ( 9855 )         ! 1979-2005, 1 file per year
                 nchunks(ifile) = 27
                 tidx1(1) =   1
                 tidx2(1) = 365
                 do ic = 2,nchunks(ifile)
                    tidx1(ic) = tidx2(ic-1) + 1
                    tidx2(ic) = tidx1(ic) + 364
                 enddo
                 tidx2(nchunks(ifile)) = ntimes(ifile,1)
              case ( 1825 ) ! "e" series; use only 2008 - maybe
                 if (trim(case_read)=='f40.amip_4k_cosp.cam4.1deg.001e') then ! Use only 2006-2008, one file per year
                    nchunks(ifile)= 3
                    tidx1(1:nchunks(ifile)) = (/ 731,1096,1461/)
                    tidx2(1:nchunks(ifile)) = (/1095,1460,1825/)
                 else
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = (/   1/)
                    tidx2(1:nchunks(ifile)) = (/1825/)
                 endif
              case default
                 nchunks(ifile)= 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',30((i6,''-'',i6),'',''))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
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
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                 endif
              enddo
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR CMOR close of '',a)') trim(cmor_filename(1:))
              stop
           else
              write(*,'(''GOOD CMOR close of '',a)') trim(cmor_filename(1:))
           endif
!           do ifile = 1,nc_nfiles(1)
 !             call close_cdf(myncid(ifile,1))
 !          enddo
        case ('parasolRefl')
           !
           ! parasolRefl
           !
!           allocate(indat3a(nlons,nlats,ncosp_sza),parasol(nlons,nlats,ncosp_sza))
           allocate(indat3a(nlons,nlats,ncosp_sza))
           do ifile = 1,nc_nfiles(1)
              call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
              call get_dims(myncid(ifile,1))
              call get_vars(myncid(ifile,1))
              !
              if (allocated(time))      deallocate(time)
              if (allocated(time_bnds)) deallocate(time_bnds)
              allocate(time(ntimes(ifile,1)))
              allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n = 1,ntimes(ifile,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
              enddo
              time_bnds(1,:) = time_bnds(2,:)
              time_bnds(2,:) = time_bnds(1,:) + 1
              time = (time_bnds(1,:)+time_bnds(2,:))/2.
              write(*,'(''time length FROM: '',a,'' myncid: '',i10,'' NT: '',i10)') trim(ncfile(ifile,1)),myncid(ifile,1),ntimes(ifile,1)
              !
              ! Determine amount of data to write, to keep close to ~2 4B limit
              !
              select case (ntimes(ifile,1))
              case ( 9855 )         ! 1979-2005, 1 file per year
                 nchunks(ifile) = 27
                 tidx1(1) =   1
                 tidx2(1) = 365
                 do ic = 2,nchunks(ifile)
                    tidx1(ic) = tidx2(ic-1) + 1
                    tidx2(ic) = tidx1(ic) + 364
                 enddo
                 tidx2(nchunks(ifile)) = ntimes(ifile,1)
              case ( 1825 ) ! "e" series; use only 2008 - maybe
                 if (trim(case_read)=='f40.amip_4k_cosp.cam4.1deg.001e') then ! Use only 2006-2008
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = (/ 731/)
                    tidx2(1:nchunks(ifile)) = (/1825/)
                 else
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = (/   1/)
                    tidx2(1:nchunks(ifile)) = (/1825/)
                 endif
              case default
                 nchunks(ifile)= 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',20((i6,''-'',i6),'',''))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
!!$ NOT NEEDED
!!$                    !
!!$                    ! Linearly interpolate between given levels (0, 15, 30, 45, 60) 
!!$                    ! and requested levels (0. 20. 40. 60. 80.)
!!$                    parasol = 1.e20
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          ! k=1 is same
!!$                          parasol(i,j,        1) = indat3a(i,j,1)
!!$                          ! k=2 is 2/3 15 (k=2) + 1/3 20 (k=3)
!!$                          parasol(i,j,        2) = ((2./3.)*indat3a(i,j,2))+((1./3.)*indat3a(i,j,3))
!!$                          ! k=3 is 1/3 30 (k=3) + 2/3 45 (k=4)
!!$                          parasol(i,j,        3) = ((1./3.)*indat3a(i,j,3))+((2./3.)*indat3a(i,j,4))
!!$                          ! k=4 is k=5
!!$                          parasol(i,j,        4) = indat3a(i,j,ncosp_sza)
!!$                          ! k=5 is set to missing (80 requires extrpolation, so avoid)
!!$                          parasol(i,j,ncosp_sza) = 1.e20
!!$                       enddo
!!$                    enddo
!!$ NOT NEEDED
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
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
              enddo
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR CMOR close of '',a)') trim(xw(ixw)%entry)
              stop
           else
              write(*,'(''GOOD CMOR close of '',a)')  trim(xw(ixw)%entry)
           endif
!           do ifile = 1,nc_nfiles(1)
 !             call close_cdf(myncid(ifile,1))
 !          enddo
        case ('clisccp')
           !
           ! clisccp
           !
           allocate(indat4a(nlons,nlats,nplev7,ncosp_tau))
           do ifile = 1,nc_nfiles(1)
              call open_cdf(myncid(ifile,1),trim(ncfile(ifile,1)),.true.)
              call get_dims(myncid(ifile,1))
              call get_vars(myncid(ifile,1))
              !
              if (allocated(time))      deallocate(time)
              if (allocated(time_bnds)) deallocate(time_bnds)
              allocate(time(ntimes(ifile,1)))
              allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n = 1,ntimes(ifile,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
              enddo
              time_bnds(1,:) = time_bnds(2,:)
              time_bnds(2,:) = time_bnds(1,:) + 1
              time = (time_bnds(1,:)+time_bnds(2,:))/2.
              write(*,'(''time length FROM: '',a,'' myncid: '',i10,'' NT: '',i10)') trim(ncfile(ifile,1)),myncid(ifile,1),ntimes(ifile,1)
              !
              ! Determine amount of data to write, to keep close to ~2 4B limit
              !
              select case (ntimes(ifile,1))
              case ( 9855 )         ! 1979-2005, 1 file per year
                 nchunks(ifile) = 27
                 tidx1(1) =   1
                 tidx2(1) = 365
                 do ic = 2,nchunks(ifile)
                    tidx1(ic) = tidx2(ic-1) + 1
                    tidx2(ic) = tidx1(ic) + 364
                 enddo
                 tidx2(nchunks(ifile)) = ntimes(ifile,1)
              case ( 1825 ) ! "e" series; use only 2008 - maybe
                 if (trim(case_read)=='f40.amip_4k_cosp.cam4.1deg.001e') then ! Use only 2006-2008
                    nchunks(ifile)= 3
                    tidx1(1:nchunks(ifile)) = (/ 731,1096,1461/)
                    tidx2(1:nchunks(ifile)) = (/1095,1460,1825/)
                 else
                    nchunks(ifile)= 1
                    tidx1(1:nchunks(ifile)) = (/   1/)
                    tidx2(1:nchunks(ifile)) = (/1825/)
                 endif
              case default
                 nchunks(ifile)= 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',20((i6,''-'',i6),'',''))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat4a)
                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,&
                         data          = indat4a,   &
                         ntimes_passed = 1,         &
                         time_vals     = tval,      &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                       stop
                    endif
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,trim(cmor_filename(1:))
                 endif
              enddo
           enddo
           cmor_filename(1:) = ' '
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR CMOR close of '',a)') trim(cmor_filename(1:))
              stop
           else
              write(*,'(''GOOD CMOR close of '',a)') trim(cmor_filename(1:))
           endif
!           do ifile = 1,nc_nfiles(1)
 !             call close_cdf(myncid(ifile,1))
 !          enddo
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
        call reset_netcdf_var
     endif
  enddo xwalk_loop
end program cfDay_CMOR
