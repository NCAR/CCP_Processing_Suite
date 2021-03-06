program Amon_CMOR
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
  use mycmor_info
  use output_times_info
  !
  implicit none
  !
  real,parameter::spval = 1.e20
  real,parameter::co2_scale = 28.966/44.
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
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,ixw,ilev,ic,jfile
  logical::does_exist
  !
  ! GO!
  !
  mycmor%table_file = 'Amon'
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
           call build_filenames(case_read,comp_read,xw(ixw)%cesm_vars(ivar),ivar,exp(exp_found)%runbeg,exp(exp_found)%runend,mycmor%table_file)
        endif
     enddo
     !
     ! Open CESM file(s) and get information(s)
     !
     if (all_continue) then
        do ivar = 1,xw(ixw)%ncesm_vars
!           write(*,'(''AVAILABLE: '',a,''.'',a,''.'',a'')') trim(case_read),trim(comp_read),trim(xw(ixw)%cesm_vars(ivar))
           write(*,*) 'AVAILABLE: ',trim(case_read),trim(comp_read),trim(xw(ixw)%cesm_vars(ivar))
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
              !
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
        call define_atm_axes(xw(ixw)%dims)
        ! 
        ! Make manual alterations so that CMOR works. Silly code!
        !
        if (xw(ixw)%ncesm_vars == 1) write(original_name,'(a)') xw(ixw)%cesm_vars(1)
        if (xw(ixw)%ncesm_vars == 2) write(original_name,'(a,'','',a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars == 3) write(original_name,'(a,'','',a,'','',a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        !
        ! Modify units as necessary to accomodate udunits' inability to convert 
        !
        select case (xw(ixw)%entry)
        case ('tauu','tauv','hfss','rlut','rlutcs','hfls','rlus','rsus','rsuscs','rsut','rsutcs','mc')
           mycmor%positive = 'up'
        case ('rlds','rldscs','rsds','rsdscs','rsdt','rtmt')
           mycmor%positive = 'down'
        case ('clt','ci','sci')
           var_info(var_found(1,1))%units = '1'
        case ('hurs','cl')
           var_info(var_found(1,1))%units = '%'
        case ('prc','pr','prsn')
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
        end select
        !
        write(*,*) 'calling cmor_variable:'
        write(*,*) 'table         = ',trim(mycmor%table_file)
        write(*,*) 'table_entry   = ',trim(xw(ixw)%entry)
        write(*,*) 'dimensions    = ',trim(xw(ixw)%dims)
        write(*,*) 'units         = ',var_info(var_found(1,1))%units(1:20)
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
                missing_value=spval,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        case ('ps')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(2),axis_ids(3),axis_ids(1)/), &
                missing_value=spval,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        case ('clw','cli','cl','mc')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(2),axis_ids(3),axis_ids(4),axis_ids(1)/),  &
                missing_value=spval,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        case default
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3)/),  &
                missing_value=spval,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        end select
        write(*,*) 'called cmor_variable:'
        write(*,*) 'varid         = ',cmor_var_id
        !
        ! Perform derivations and cycle through time, writing data too
        !
        select case (xw(ixw)%entry)
        case ('ccb','cct','clivi','clwvi','evspsbl','hfls','hfss','hurs','huss',&
             'prw','psl','ps','rldscs','rlds','rlutcs','rsdscs','rsds','rsdt',&
             'sci','tas','tasmax','tasmin','tauu','tauv','ts','ci','clt','sfcWind')
           !
           ! No change
           !
           allocate(indat2a(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (allocated(time))       deallocate(time)
           if (allocated(time_bnds))  deallocate(time_bnds)
           allocate(time(ntimes(1,1)))
           allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
!              write(*,'(''t '',i6,'': '',3f12.4)') n,time_bnds(1,n),time(n),time_bnds(2,n)
           enddo
           !
           ! Determine amount of data to write, to keep close to ~4 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872,1865,1860,51140,3612,6012,12012,1140 )  ! All data
              select case(exp(exp_found)%model_id)
              case ('CESM1-WACCM')
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 13
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              case default
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 1
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              end select
           case ( 12000 ) ! BGC controls
              select case(exp(exp_found)%model_id)
              case ('CESM1-BGC')
                 select case(exp(exp_found)%expt_id)
                 case ('piControl') ! b40.prescribed_carb.001, 0101 - 0600
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 1201
                    tidx2(1:nchunks(1)) = 7200
                 case ('esmControl') ! b40.coup_carb.004, 0301 - 0800
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 3601
                    tidx2(1:nchunks(1)) = 9600
                 end select
              end select
           case ( 3228 )  ! Abrupt 4XCO2, use 1850-2000 (151 years)
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = 1812
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 13
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 3613
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 where (abs(indat2a) > spval)
                    indat2a = spval
                 endwhere
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
              if (ic < nchunks(1)) then
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
        case ('prc')
           !
           ! prc : PRECC, unit change from m s-1 to kg m-2 s-1
           !
           allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (allocated(time))       deallocate(time)
           if (allocated(time_bnds))  deallocate(time_bnds)
           allocate(time(ntimes(1,1)))
           allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           !
           ! Determine amount of data to write, to keep close to ~4 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872,1865,1860,1140,3612,6012,12012 )  ! All data
              select case(exp(exp_found)%model_id)
              case ('CESM1-WACCM')
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 13
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              case default
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 1
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              end select
           case ( 12000 ) ! BGC controls
              select case(exp(exp_found)%model_id)
              case ('CESM1-BGC')
                 select case(exp(exp_found)%expt_id)
                 case ('piControl') ! b40.prescribed_carb.001, 0101 - 0600
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 1201
                    tidx2(1:nchunks(1)) = 7200
                 case ('esmControl') ! b40.coup_carb.004, 0301 - 0800
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 3601
                    tidx2(1:nchunks(1)) = 9600
                 end select
              end select
           case ( 3228 )  ! Abrupt 4XCO2, use 1850-2000 (151 years)
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = 1812
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 13
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 3613
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 where (abs(indat2a) > spval)
                    indat2a = spval
                 endwhere
                 ! 
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
              if (ic < nchunks(1)) then
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
        case ('pr','prsn')
           !
           ! pr  : Add PRECC + PRECL  , unit change from m s-1 to kg m-2 s-1
           ! prsn: Add PRECSC + PRECSL, unit change from m s-1 to kg m-2 s-1
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (allocated(time))       deallocate(time)
           if (allocated(time_bnds))  deallocate(time_bnds)
           allocate(time(ntimes(1,1)))
           allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           !
           ! Determine amount of data to write, to keep close to ~4 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872,1865,1860,1140,3612,6012,12012 )  ! All data
              select case(exp(exp_found)%model_id)
              case ('CESM1-WACCM')
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 13
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              case default
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 1
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              end select
           case ( 12000 ) ! BGC controls
              select case(exp(exp_found)%model_id)
              case ('CESM1-BGC')
                 select case(exp(exp_found)%expt_id)
                 case ('piControl') ! b40.prescribed_carb.001, 0101 - 0600
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 1201
                    tidx2(1:nchunks(1)) = 7200
                 case ('esmControl') ! b40.coup_carb.004, 0301 - 0800
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 3601
                    tidx2(1:nchunks(1)) = 9600
                 end select
              end select
           case ( 3228 )  ! Abrupt 4XCO2, use 1850-2000 (151 years)
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = 1812
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 13
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 3613
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 cmordat2d = spval
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
                 where (abs(indat2a) > spval)
                    indat2a = spval
                 endwhere
                 where (abs(indat2b) > spval)
                    indat2b = spval
                 endwhere
                 ! 
                 where ((indat2a /= spval).and.(indat2b /= spval))
                    cmordat2d = (indat2a + indat2b)*1000.
                 elsewhere
                    cmordat2d = spval
                 endwhere
                 write(*,*) 'min: ',minval(cmordat2d,mask=cmordat2d/=0.),' max: ',maxval(cmordat2d,mask=cmordat2d/=0.)
                 where (cmordat2d < 1.449e-06)
                    cmordat2d = 1.449e-06
                 endwhere
                 write(*,*) 'min: ',minval(cmordat2d,mask=cmordat2d/=0.),' max: ',maxval(cmordat2d,mask=cmordat2d/=0.)
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
              if (ic < nchunks(1)) then
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
        case ('rlut')
           !
           ! rlut : FSNTOA-FSNT+FLNT
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (allocated(time))      deallocate(time)
           if (allocated(time_bnds)) deallocate(time_bnds)
           allocate(time(ntimes(1,1)))
           allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           !
           ! Determine amount of data to write, to keep close to ~4 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872,1865,1860,1140,3612,6012,12012 )  ! All data
              select case(exp(exp_found)%model_id)
              case ('CESM1-WACCM')
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 13
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              case default
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 1
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              end select
           case ( 12000 ) ! BGC controls
              select case(exp(exp_found)%model_id)
              case ('CESM1-BGC')
                 select case(exp(exp_found)%expt_id)
                 case ('piControl') ! b40.prescribed_carb.001, 0101 - 0600
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 1201
                    tidx2(1:nchunks(1)) = 7200
                 case ('esmControl') ! b40.coup_carb.004, 0301 - 0800
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 3601
                    tidx2(1:nchunks(1)) = 9600
                 end select
              end select
           case ( 3228 )  ! Abrupt 4XCO2, use 1850-2000 (151 years)
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = 1812
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 13
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 3613
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 cmordat2d = spval
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat2c)
                 where (abs(indat2a) > spval)
                    indat2a = spval
                 endwhere
                 where (abs(indat2b) > spval)
                    indat2b = spval
                 endwhere
                 where (abs(indat2c) > spval)
                    indat2c = spval
                 endwhere
                 ! 
                 where ((indat2a /= spval).and.(indat2b /= spval).and.(indat2b /= spval))
                    cmordat2d = indat2a - indat2b + indat2c
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
              if (ic < nchunks(1)) then
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
        case ('rlus')
           !
           ! rlus: Add FLDS + FLNS
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (allocated(time))       deallocate(time)
           if (allocated(time_bnds))  deallocate(time_bnds)
           allocate(time(ntimes(1,1)))
           allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           !
           ! Determine amount of data to write, to keep close to ~4 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872,1865,1860,1140,3612,6012,12012 )  ! All data
              select case(exp(exp_found)%model_id)
              case ('CESM1-WACCM')
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 13
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              case default
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 1
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              end select
           case ( 12000 ) ! BGC controls
              select case(exp(exp_found)%model_id)
              case ('CESM1-BGC')
                 select case(exp(exp_found)%expt_id)
                 case ('piControl') ! b40.prescribed_carb.001, 0101 - 0600
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 1201
                    tidx2(1:nchunks(1)) = 7200
                 case ('esmControl') ! b40.coup_carb.004, 0301 - 0800
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 3601
                    tidx2(1:nchunks(1)) = 9600
                 end select
              end select
           case ( 3228 )  ! Abrupt 4XCO2, use 1850-2000 (151 years)
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = 1812
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 13
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 3613
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 cmordat2d = spval
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
                 where (abs(indat2a) > spval)
                    indat2a = spval
                 endwhere
                 where (abs(indat2b) > spval)
                    indat2b = spval
                 endwhere
                 ! 
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
              if (ic < nchunks(1)) then
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
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (allocated(time))       deallocate(time)
           if (allocated(time_bnds))  deallocate(time_bnds)
           allocate(time(ntimes(1,1)))
           allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           !
           ! Determine amount of data to write, to keep close to ~4 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872,1865,1860,1140,3612,6012,12012 )  ! All data
              select case(exp(exp_found)%model_id)
              case ('CESM1-WACCM')
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 13
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              case default
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 1
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              end select
           case ( 12000 ) ! BGC controls
              select case(exp(exp_found)%model_id)
              case ('CESM1-BGC')
                 select case(exp(exp_found)%expt_id)
                 case ('piControl') ! b40.prescribed_carb.001, 0101 - 0600
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 1201
                    tidx2(1:nchunks(1)) = 7200
                 case ('esmControl') ! b40.coup_carb.004, 0301 - 0800
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 3601
                    tidx2(1:nchunks(1)) = 9600
                 end select
              end select
           case ( 3228 )  ! Abrupt 4XCO2, use 1850-2000 (151 years)
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = 1812
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 13
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 3613
              tidx2(1:nchunks(1)) = ntimes(1,1)
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 cmordat2d = spval
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
                 where (abs(indat2a) > spval)
                    indat2a = spval
                 endwhere
                 where (abs(indat2b) > spval)
                    indat2b = spval
                 endwhere
                 ! 
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
              if (ic < nchunks(1)) then
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
        case ('ta','ua','va','hus','hur','wap','zg',&
              'tro3','tro3Clim','co2Clim','ch4','ch4Clim','n2o','n2oClim')
           select case(exp(exp_found)%model_id)
           case ('CESM1-WACCM')
              !
              do ivar = 1,xw(ixw)%ncesm_vars
                 do ifile = 1,nc_nfiles(ivar)
                    call open_cdf(myncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                    call get_dims(myncid(ifile,ivar))
                    call get_vars(myncid(ifile,ivar))
                    if (allocated(time))       deallocate(time)
                    if (allocated(time_bnds))  deallocate(time_bnds)
                    allocate(time(ntimes(ifile,1)))
                    allocate(time_bnds(2,ntimes(ifile,1)))
                    !
                    do n=1,ntimes(ifile,ivar)
                       time_counter = n
                       call read_var(myncid(ifile,ivar),'time_bnds',time_bnds(:,n))
                       time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
                    enddo
                 enddo
              enddo
              !
              ! Vertically interpolate to standard pressure levels
              !
              allocate(indat3a(nlons,nlats,nlevs),cmordat3d(nlons,nlats,nplev23))
              allocate(psdata(nlons,nlats))
              !
              ! Determine amount of data to write, to keep close to ~4 GB limit
              !
              if (ntimes(1,1)==1140) then ! RCP 2005-2099, keep only 2006-2099
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 13
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              else
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 1
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              endif
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
              do ic = 1,nchunks(1)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,psdata)
                    where (abs(indat3a) > spval)
                       indat3a = spval
                    endwhere
                    where (abs(psdata) > spval)
                       psdata = spval
                    elsewhere
                       psdata = psdata * 0.01
                    endwhere
                    !
                    cmordat3d = spval
                    !
                    ! Do vertical interpolation to pressure levels
                    !
                    call vertint(indat3a,cmordat3d,atm_levs,atm_plev23*0.01,psdata,spval,nlons,nlats,nlevs,nlevs+1,nplev23)
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
                 if (ic < nchunks(1)) then
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
           case default
              ! 
              ! All other models
              !
              allocate(indat3a(nlons,nlats,nlevs),cmordat3d(nlons,nlats,nplev23))
              allocate(psdata(nlons,nlats))
              do ifile = 1,nc_nfiles(1)
                 if (allocated(time))       deallocate(time)
                 if (allocated(time_bnds))  deallocate(time_bnds)
                 allocate(time(ntimes(ifile,1)))
                 allocate(time_bnds(2,ntimes(ifile,1)))
                 !
                 do n=1,ntimes(ifile,1)
                    time_counter = n
                    call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
                    time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
                 enddo
                 !
                 ! Determine amount of data to write, to keep close to ~4 GB limit
                 !
                 select case(ntimes(ifile,1))
                 case ( 1872,1860 )  ! 20C, 1850-2005, ~50y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  1, 601,1201/) ! 1850, 1900, 1951
                    tidx2(1:nchunks(ifile)) = (/600,1200,ntimes(ifile,1)/) ! 1899, 1950, 2005
                 case ( 1865 )  ! 20C, 1850-2005, ~50y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  1, 601,1198/) ! 1850, 1900, 1951
                    tidx2(1:nchunks(ifile)) = (/600,1197,1865/) ! 1899, 1950, 2005
                 case ( 1212 )  ! Pliocene from 450-550
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 601/) ! 0450, 0500
                    tidx2(1:nchunks(ifile)) = (/600,1212/) ! 0499, 0550
                 case ( 3828 )  ! CAM5 piControl
                    nchunks(ifile) = 32
                    tidx1(1) =   1
                    tidx2(1) = 108
                    do ic = 2,nchunks(ifile)
                       tidx1(ic) = tidx2(ic-1) + 1
                       tidx2(ic) = tidx1(ic) + 119
                    enddo
                    tidx2(nchunks(ifile)) = ntimes(ifile,1)
                 case ( 3228 )  ! Abrupt 4XCO2, use 1850-2000 (151 years)
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  1, 601,1201/) ! 1850, 1900, 1951
                    tidx2(1:nchunks(ifile)) = (/600,1200,1812/) ! 1899, 1950, 2000
                 case ( 1152 )  ! RCP, 2005-2100, skip 2006
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/ 13, 541/)      ! 2006, 2050
                    tidx2(1:nchunks(ifile)) = (/540,1152/)      ! 2049, 2100
                 case ( 1248 )  ! COSP 4XCO2
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 601/)
                    tidx2(1:nchunks(ifile)) = (/600,1248/)
                 case ( 1140 )  ! RCP, 2006-2100
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 529/)      ! 2006, 2050
                    tidx2(1:nchunks(ifile)) = (/528,1140/)      ! 2049, 2100
                 case ( 2664 )  ! FASTCHEM piControl
                    nchunks(ifile) = 5
                    tidx1(1:nchunks(ifile)) = (/  1, 361, 961,1561,2161/) ! 0070,0100,0150,0200,0250
                    tidx2(1:nchunks(ifile)) = (/360, 960,1560,2160,2664/) ! 0099,0149,0199,0249,0291
                 case ( 828 )
                    select case (case_read)
                    case ( 'b40.1850_ramp_solar.beta19.005')
                       nchunks(ifile) = 2
                       tidx1(1:nchunks(ifile)) = (/  1, 589/)      ! 1850, 1900
                       tidx2(1:nchunks(ifile)) = (/588, 828/)      ! 1899, 1919
                    case default
                       nchunks(ifile) = 2
                       tidx1(1:nchunks(ifile)) = (/  1, 601/)      ! 1850, 1900
                       tidx2(1:nchunks(ifile)) = (/600, 828/)      ! 1899, 1918
                    end select
                 case ( 829 )
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 590/)      ! 1850, 1900
                    tidx2(1:nchunks(ifile)) = (/589, 829/)      ! 1899, 1918
                 case ( 900 )
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 601/)      ! 1850, 1900
                    tidx2(1:nchunks(ifile)) = (/600, 900/)      ! 1899, 1924
                 case ( 816 )
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 337/)      ! 1850, 1900
                    tidx2(1:nchunks(ifile)) = (/336, 816/)      ! 1899, 1924
                 case ( 1680,3612,6012,12012 ) ! piControl,past1000,midHolocene: ~50Y chunks
                    nchunks(ifile) = ntimes(1,1)/600
                    tidx1(1) =   1
                    tidx2(1) = 600
                    do ic = 2,nchunks(ifile)
                       tidx1(ic) = tidx2(ic-1) + 1
                       tidx2(ic) = tidx1(ic) + 599
                    enddo
                    tidx2(nchunks(ifile)) = ntimes(1,1)
                 case ( 12000 ) ! BGC controls
                    nchunks(ifile) = 10
                    select case(exp(exp_found)%model_id)
                    case ('CESM1-BGC')
                       select case(exp(exp_found)%expt_id)
                       case ('piControl') ! b40.prescribed_carb.001, 0101 - 0600
                          tidx1(1) = 1201
                          tidx2(1) = 1800
                          do ic = 2,nchunks(ifile)
                             tidx1(ic) = tidx2(ic-1) + 1
                             tidx2(ic) = tidx1(ic) + 599
                          enddo
                       case ('esmControl') ! b40.coup_carb.004, 0301 - 0800
                          tidx1(1) = 3601
                          tidx2(1) = 4200
                          do ic = 2,nchunks(ifile)
                             tidx1(ic) = tidx2(ic-1) + 1
                             tidx2(ic) = tidx1(ic) + 599
                          enddo
                       end select
                    end select
                 case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only, ~50y chunks
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/3613,4213/) ! 1850, 1900, 1951
                    tidx2(1:nchunks(ifile)) = (/4212,4824/) ! 1899, 1950, 2005
                 case ( 2388,2400 )
                    nchunks(ifile) = 4
                    tidx1(1:nchunks(ifile)) = (/   1, 589,1189,1789/)
                    tidx2(1:nchunks(ifile)) = (/ 588,1188,1788,ntimes(ifile,1)/)
                 case default
                    nchunks(ifile) = 1
                    tidx1(1:nchunks(ifile)) = 1
                    tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
                 end select
                 write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                       call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,psdata)
                       where (abs(indat3a) > spval)
                          indat3a = spval
                       endwhere
                       where (abs(psdata) > spval)
                          psdata = spval
                       elsewhere
                          psdata = psdata * 0.01
                       endwhere
                       !
                       cmordat3d = spval
                       !
                       ! Do vertical interpolation to pressure levels
                       !
                       call vertint(indat3a,cmordat3d,atm_levs,atm_plev17*0.01,psdata,spval,nlons,nlats,nlevs,nlevs+1,nplev17)
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
                    cmor_filename(1:) = ' '
                    error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                    if (error_flag < 0) then
                       write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                       stop
                    else
                       write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                    endif
                 enddo
              enddo
           end select
        case ('co2')
           select case(exp(exp_found)%model_id)
           case ('CESM1-WACCM')
              !
              do ivar = 1,xw(ixw)%ncesm_vars
                 do ifile = 1,nc_nfiles(ivar)
                    call open_cdf(myncid(ifile,ivar),trim(ncfile(ifile,ivar)),.true.)
                    call get_dims(myncid(ifile,ivar))
                    call get_vars(myncid(ifile,ivar))
                    if (allocated(time))       deallocate(time)
                    if (allocated(time_bnds))  deallocate(time_bnds)
                    allocate(time(ntimes(ifile,1)))
                    allocate(time_bnds(2,ntimes(ifile,1)))
                    !
                    do n=1,ntimes(ifile,ivar)
                       time_counter = n
                       call read_var(myncid(ifile,ivar),'time_bnds',time_bnds(:,n))
                       time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
                    enddo
                 enddo
              enddo
              !
              ! Vertically interpolate to standard pressure levels
              !
              allocate(indat3a(nlons,nlats,nlevs),cmordat3d(nlons,nlats,nplev23))
              allocate(psdata(nlons,nlats))
              !
              ! Determine amount of data to write, to keep close to ~4 GB limit
              !
              if (ntimes(1,1)==1140) then ! RCP 2005-2099, keep only 2006-2099
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 13
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              else
                 nchunks(1) = 1
                 tidx1(1:nchunks(1)) = 1
                 tidx2(1:nchunks(1)) = ntimes(1,1)
              endif
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
              do ic = 1,nchunks(1)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,psdata)
                    indat3a = indat3a * co2_scale ! Get proper mass units
                    where (abs(indat3a) > spval)
                       indat3a = spval
                    endwhere
                    where (abs(psdata) > spval)
                       psdata = spval
                    elsewhere
                       psdata = psdata * 0.01
                    endwhere
                    !
                    cmordat3d = spval
                    !
                    ! Do vertical interpolation to pressure levels
                    !
                    call vertint(indat3a,cmordat3d,atm_levs,atm_plev23*0.01,psdata,spval,nlons,nlats,nlevs,nlevs+1,nplev23)
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
                 if (ic < nchunks(1)) then
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
           case default
              ! 
              ! All other models
              !
              allocate(indat3a(nlons,nlats,nlevs),cmordat3d(nlons,nlats,nplev23))
              allocate(psdata(nlons,nlats))
              do ifile = 1,nc_nfiles(1)
                 if (allocated(time))       deallocate(time)
                 if (allocated(time_bnds))  deallocate(time_bnds)
                 allocate(time(ntimes(ifile,1)))
                 allocate(time_bnds(2,ntimes(ifile,1)))
                 !
                 do n=1,ntimes(ifile,1)
                    time_counter = n
                    call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
                    time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
                 enddo
                 !
                 ! Determine amount of data to write, to keep close to ~4 GB limit
                 !
                 select case(ntimes(ifile,1))
                 case ( 1872,1860 )  ! 20C, 1850-2005, ~50y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  1, 601,1201/) ! 1850, 1900, 1951
                    tidx2(1:nchunks(ifile)) = (/600,1200,ntimes(ifile,1)/) ! 1899, 1950, 2005
                 case ( 1865 )  ! 20C, 1850-2005, ~50y chunks
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  1, 601,1198/) ! 1850, 1900, 1951
                    tidx2(1:nchunks(ifile)) = (/600,1197,1865/) ! 1899, 1950, 2005
                 case ( 1212 )  ! Pliocene from 450-550
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 601/) ! 0450, 0500
                    tidx2(1:nchunks(ifile)) = (/600,1212/) ! 0499, 0550
                 case ( 3828 )  ! CAM5 piControl
                    nchunks(ifile) = 32
                    tidx1(1) =   1
                    tidx2(1) = 108
                    do ic = 2,nchunks(ifile)
                       tidx1(ic) = tidx2(ic-1) + 1
                       tidx2(ic) = tidx1(ic) + 119
                    enddo
                    tidx2(nchunks(ifile)) = ntimes(ifile,1)
                 case ( 3228 )  ! Abrupt 4XCO2, use 1850-2000 (151 years)
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  1, 601,1201/) ! 1850, 1900, 1951
                    tidx2(1:nchunks(ifile)) = (/600,1200,1812/) ! 1899, 1950, 2000
                 case ( 1152 )  ! RCP, 2005-2100, skip 2006
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/ 13, 541/)      ! 2006, 2050
                    tidx2(1:nchunks(ifile)) = (/540,1152/)      ! 2049, 2100
                 case ( 1248 )  ! COSP 4XCO2
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 601/)
                    tidx2(1:nchunks(ifile)) = (/600,1248/)
                 case ( 1140 )  ! RCP, 2006-2100
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 529/)      ! 2006, 2050
                    tidx2(1:nchunks(ifile)) = (/528,1140/)      ! 2049, 2100
                 case ( 2664 )  ! FASTCHEM piControl
                    nchunks(ifile) = 5
                    tidx1(1:nchunks(ifile)) = (/  1, 361, 961,1561,2161/) ! 0070,0100,0150,0200,0250
                    tidx2(1:nchunks(ifile)) = (/360, 960,1560,2160,2664/) ! 0099,0149,0199,0249,0291
                 case ( 828 )
                    select case (case_read)
                    case ( 'b40.1850_ramp_solar.beta19.005')
                       nchunks(ifile) = 2
                       tidx1(1:nchunks(ifile)) = (/  1, 589/)      ! 1850, 1900
                       tidx2(1:nchunks(ifile)) = (/588, 828/)      ! 1899, 1919
                    case default
                       nchunks(ifile) = 2
                       tidx1(1:nchunks(ifile)) = (/  1, 601/)      ! 1850, 1900
                       tidx2(1:nchunks(ifile)) = (/600, 828/)      ! 1899, 1918
                    end select
                 case ( 829 )
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 590/)      ! 1850, 1900
                    tidx2(1:nchunks(ifile)) = (/589, 829/)      ! 1899, 1918
                 case ( 900 )
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 601/)      ! 1850, 1900
                    tidx2(1:nchunks(ifile)) = (/600, 900/)      ! 1899, 1924
                 case ( 816 )
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 337/)      ! 1850, 1900
                    tidx2(1:nchunks(ifile)) = (/336, 816/)      ! 1899, 1924
                 case ( 1680,3612,6012,12012 ) ! piControl,past1000,midHolocene: ~50Y chunks
                    nchunks(ifile) = ntimes(1,1)/600
                    tidx1(1) =   1
                    tidx2(1) = 600
                    do ic = 2,nchunks(ifile)
                       tidx1(ic) = tidx2(ic-1) + 1
                       tidx2(ic) = tidx1(ic) + 599
                    enddo
                    tidx2(nchunks(ifile)) = ntimes(1,1)
                 case ( 12000 ) ! BGC controls
                    nchunks(ifile) = 10
                    select case(exp(exp_found)%model_id)
                    case ('CESM1-BGC')
                       select case(exp(exp_found)%expt_id)
                       case ('piControl') ! b40.prescribed_carb.001, 0101 - 0600
                          tidx1(1) = 1201
                          tidx2(1) = 1800
                          do ic = 2,nchunks(ifile)
                             tidx1(ic) = tidx2(ic-1) + 1
                             tidx2(ic) = tidx1(ic) + 599
                          enddo
                       case ('esmControl') ! b40.coup_carb.004, 0301 - 0800
                          tidx1(1) = 3601
                          tidx2(1) = 4200
                          do ic = 2,nchunks(ifile)
                             tidx1(ic) = tidx2(ic-1) + 1
                             tidx2(ic) = tidx1(ic) + 599
                          enddo
                       end select
                    end select
                 case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only, ~50y chunks
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/3613,4213/) ! 1850, 1900, 1951
                    tidx2(1:nchunks(ifile)) = (/4212,4824/) ! 1899, 1950, 2005
                 case ( 2388,2400 )
                    nchunks(ifile) = 4
                    tidx1(1:nchunks(ifile)) = (/   1, 589,1189,1789/)
                    tidx2(1:nchunks(ifile)) = (/ 588,1188,1788,ntimes(ifile,1)/)
                 case default
                    nchunks(ifile) = 1
                    tidx1(1:nchunks(ifile)) = 1
                    tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
                 end select
                 write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
                 do ic = 1,nchunks(ifile)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                       call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,psdata)
                       indat3a = indat3a * co2_scale ! Get proper mass units
                       where (abs(indat3a) > spval)
                          indat3a = spval
                       endwhere
                       where (abs(psdata) > spval)
                          psdata = spval
                       elsewhere
                          psdata = psdata * 0.01
                       endwhere
                       !
                       cmordat3d = spval
                       !
                       ! Do vertical interpolation to pressure levels
                       !
                       call vertint(indat3a,cmordat3d,atm_levs,atm_plev17*0.01,psdata,spval,nlons,nlats,nlevs,nlevs+1,nplev17)
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
                    cmor_filename(1:) = ' '
                    error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                    if (error_flag < 0) then
                       write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                       stop
                    else
                       write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                    endif
                 enddo
              enddo
           end select
        case ('clw','cli','cl')
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nlevs),psdata(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (allocated(time))       deallocate(time)
           if (allocated(time_bnds))  deallocate(time_bnds)
           allocate(time(ntimes(1,1)))
           allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           !
           ! Determine amount of data to write, to keep close to ~4 GB limit
           !
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case ( 1865 )  ! 20C, 1850-2005, ~50y chunks
                 nchunks(ifile) = 3
                 tidx1(1:nchunks(ifile)) = (/  1, 601,1198/) ! 1850, 1900, 1951
                 tidx2(1:nchunks(ifile)) = (/600,1197,1865/) ! 1899, 1950, 2005
              case ( 1871 )  ! 20C, 1850-2005, ~50y chunks
                 nchunks(ifile) = 3
                 tidx1(1:nchunks(ifile)) = (/  1, 600,1200/) ! 1850, 1900, 1951
                 tidx2(1:nchunks(ifile)) = (/599,1199,1871/) ! 1899, 1950, 2005
              case ( 1212 )  ! Pliocene from 450-550
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/  1, 601/) ! 0450, 0500
                 tidx2(1:nchunks(ifile)) = (/600,1212/) ! 0499, 0550
              case ( 1872,1860 )  ! 20C, 1850-2005
                 select case(exp(exp_found)%model_id)
                 case ('CESM1-CAM5')
                    nchunks(ifile) = 6
                    tidx1(1:nchunks(ifile)) = (/  1, 301, 601, 901,1201,1501/) ! 1850, 1875, 1900, 1925, 1950, 1975
                    tidx2(1:nchunks(ifile)) = (/300, 600, 900,1200,1500,1872/) ! 1874, 1899, 1924, 1949, 1974, 2005
                 case default
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  1, 601,1201/) ! 1850, 1900, 1951
                    tidx2(1:nchunks(ifile)) = (/600,1200,ntimes(ifile,1)/) ! 1899, 1950, 2005
                 end select
              case ( 3828 )  ! CAM5 piControl
                 nchunks(ifile) = 32
                 tidx1(1) =   1
                 tidx2(1) = 108
                 do ic = 2,nchunks(ifile)
                    tidx1(ic) = tidx2(ic-1) + 1
                    tidx2(ic) = tidx1(ic) + 119
                 enddo
                 tidx2(nchunks(ifile)) = ntimes(ifile,1)
              case ( 3228 )  ! Abrupt 4XCO2, use 1850-2000 (151 years)
                 nchunks(ifile) = 3
                 tidx1(1:nchunks(ifile)) = (/  1, 601,1201/) ! 1850, 1900, 1951
                 tidx2(1:nchunks(ifile)) = (/600,1200,1812/) ! 1899, 1950, 2000
              case ( 1152 )  ! RCP, 2005-2100, skip 2006
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/ 13, 541/)      ! 2006, 2050
                 tidx2(1:nchunks(ifile)) = (/540,1152/)      ! 2049, 2100
              case ( 1248 )
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/  1, 601/)
                 tidx2(1:nchunks(ifile)) = (/600,1248/)
              case ( 1140 )  ! RCP, 2006-2100
                 select case(exp(exp_found)%model_id)
                 case ('CESM1-WACCM')
                    nchunks(ifile) = 1
                    tidx1(1:nchunks(ifile)) = 13
                    tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
                 case ('CESM1-CAM5')
                    nchunks(ifile) = 2 ! 1 deg but 30 levels
                    tidx1(1:nchunks(ifile)) = (/   1, 529/)
                    tidx2(1:nchunks(ifile)) = (/ 528,1140/)
                 case default
                    nchunks(ifile) = 1
                    tidx1(1:nchunks(ifile)) = 1
                    tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
                 end select
              case ( 2664 )  ! FASTCHEM piControl
                 nchunks(ifile) = 5
                 tidx1(1:nchunks(ifile)) = (/  1, 361, 961,1561,2161/) ! 0070,0100,0150,0200,0250
                 tidx2(1:nchunks(ifile)) = (/360, 960,1560,2160,2664/) ! 0099,0149,0199,0249,0291
              case ( 828 )
                 select case (case_read)
                 case ( 'b40.1850_ramp_solar.beta19.005')
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 589/)      ! 1850, 1900
                    tidx2(1:nchunks(ifile)) = (/588, 828/)      ! 1899, 1919
                 case default
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 601/)      ! 1850, 1900
                    tidx2(1:nchunks(ifile)) = (/600, 828/)      ! 1899, 1918
                 end select
              case ( 900 )
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/  1, 601/)      ! 1850, 1900
                 tidx2(1:nchunks(ifile)) = (/600, 900/)      ! 1899, 1924
              case ( 816 )
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/  1, 337/)      ! 1850, 1900
                 tidx2(1:nchunks(ifile)) = (/336, 816/)      ! 1899, 1924
              case ( 840 )
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/  1, 349/)
                 tidx2(1:nchunks(ifile)) = (/348, 840/)
              case ( 1680,3612,6012,12012 ) ! piControl,past1000,midHolocene: ~50Y chunks
                 nchunks(ifile) = int(ntimes(ifile,1)/600)+1
                 tidx1(1) =   1
                 tidx2(1) = 600
                 do ic = 2,nchunks(ifile)
                    tidx1(ic) = tidx2(ic-1) + 1
                    tidx2(ic) = tidx1(ic) + 599
                 enddo
                 tidx2(nchunks(ifile)) = ntimes(ifile,1)
              case ( 12000 ) ! BGC controls
                 nchunks(ifile) = 10
                 select case(exp(exp_found)%model_id)
                 case ('CESM1-BGC')
                    select case(exp(exp_found)%expt_id)
                    case ('piControl') ! b40.prescribed_carb.001, 0101 - 0600
                       tidx1(1) = 1201
                       tidx2(1) = 1800
                       do ic = 2,nchunks(ifile)
                          tidx1(ic) = tidx2(ic-1) + 1
                          tidx2(ic) = tidx1(ic) + 599
                       enddo
                    case ('esmControl') ! b40.coup_carb.004, 0301 - 0800
                       tidx1(1) = 3601
                       tidx2(1) = 4200
                       do ic = 2,nchunks(ifile)
                          tidx1(ic) = tidx2(ic-1) + 1
                          tidx2(ic) = tidx1(ic) + 599
                       enddo
                    end select
                 end select
              case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only, ~50y chunks
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/3613,4213/) ! 1850, 1900, 1951
                 tidx2(1:nchunks(ifile)) = (/4212,4824/) ! 1899, 1950, 2005
              case ( 2388,2400 )
                 select case(exp(exp_found)%model_id)
                 case ('CESM1-WACCM')
                    nchunks(ifile) = 4
                    tidx1(1:nchunks(ifile)) = (/   1, 649,1249,1849/) !   96,  150, 200, 250
                    tidx2(1:nchunks(ifile)) = (/ 648,1248,1848,ntimes(ifile,1)/) !  149,  199, 249, 295
                 case default
                    nchunks(ifile) = 4
                    tidx1(1:nchunks(ifile)) = (/   1, 589,1189,1789/)
                    tidx2(1:nchunks(ifile)) = (/ 588,1188,1788,ntimes(ifile,1)/)
                 end select
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,psdata)
                    where (abs(indat3a) > spval)
                       indat3a = spval
                    endwhere
                    where (abs(psdata) > spval)
                       psdata = spval
                    endwhere
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
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
                         data          = psdata,   &
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
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo
        case ('mc')
           !
           ! mc: CMFMC + CMFMCDZM
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nilevs),indat3b(nlons,nlats,nilevs),psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nilevs))
           !
           do ifile = 1,nc_nfiles(1)
              if (allocated(time))       deallocate(time)
              if (allocated(time_bnds))  deallocate(time_bnds)
              allocate(time(ntimes(ifile,1)))
              allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n=1,ntimes(ifile,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
                 time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
              enddo
              !
              ! Determine amount of data to write, to keep close to ~4 GB limit
              !
              select case(ntimes(ifile,1))
              case ( 1865 )  ! 20C, 1850-2005, ~50y chunks
                 nchunks(ifile) = 3
                 tidx1(1:nchunks(ifile)) = (/  1, 601,1198/) ! 1850, 1900, 1951
                 tidx2(1:nchunks(ifile)) = (/600,1197,1865/) ! 1899, 1950, 2005
              case ( 1871 )  ! 20C, 1850-2005, ~50y chunks
                 nchunks(ifile) = 3
                 tidx1(1:nchunks(ifile)) = (/  1, 600,1200/) ! 1850, 1900, 1951
                 tidx2(1:nchunks(ifile)) = (/599,1199,1871/) ! 1899, 1950, 2005
              case ( 1212 )  ! Pliocene from 450-550
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/  1, 601/) ! 0450, 0500
                 tidx2(1:nchunks(ifile)) = (/600,1212/) ! 0499, 0550
              case ( 1872,1860 )  ! 20C, 1850-2005
                 select case(exp(exp_found)%model_id)
                 case ('CESM1-CAM5')
                    nchunks(ifile) = 6
                    tidx1(1:nchunks(ifile)) = (/  1, 301, 601, 901,1201,1501/) ! 1850, 1875, 1900, 1925, 1950, 1975
                    tidx2(1:nchunks(ifile)) = (/300, 600, 900,1200,1500,1872/) ! 1874, 1899, 1924, 1949, 1974, 2005
                 case default
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  1, 601,1201/) ! 1850, 1900, 1951
                    tidx2(1:nchunks(ifile)) = (/600,1200,ntimes(ifile,1)/) ! 1899, 1950, 2005
                 end select
              case ( 3228 )  ! Abrupt 4XCO2, use 1850-2000 (151 years)
                 nchunks(ifile) = 3
                 tidx1(1:nchunks(ifile)) = (/  1, 601,1201/) ! 1850, 1900, 1951
                 tidx2(1:nchunks(ifile)) = (/600,1200,1812/) ! 1899, 1950, 2000
              case ( 1152 )  ! RCP, 2005-2100, skip 2006
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/ 13, 541/)      ! 2006, 2050
                 tidx2(1:nchunks(ifile)) = (/540,1152/)      ! 2049, 2100
              case ( 1248 )
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/  1, 601/)
                 tidx2(1:nchunks(ifile)) = (/600,1248/)
              case ( 1140 )  ! RCP, 2006-2100
                 select case(exp(exp_found)%model_id)
                 case ('CESM1-CAM5')
                    nchunks(ifile) = 3
                    tidx1(1:nchunks(ifile)) = (/  1, 409, 889/)
                    tidx2(1:nchunks(ifile)) = (/408, 888,1140/)
                 case ('CESM1-WACCM')
                    nchunks(ifile) = 1
                    tidx1(1:nchunks(ifile)) = 13
                    tidx2(1:nchunks(ifile)) = ntimes(1,1)
                 case default
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 529/)      ! 2006, 2050
                    tidx2(1:nchunks(ifile)) = (/528,1140/)      ! 2049, 2100
                 end select
              case ( 2664 )  ! FASTCHEM piControl
                 nchunks(ifile) = 5
                 tidx1(1:nchunks(ifile)) = (/  1, 361, 961,1561,2161/) ! 0070,0100,0150,0200,0250
                 tidx2(1:nchunks(ifile)) = (/360, 960,1560,2160,2664/) ! 0099,0149,0199,0249,0291
              case ( 828 )
                 select case (case_read)
                 case ( 'b40.1850_ramp_solar.beta19.005')
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 589/)      ! 1850, 1900
                    tidx2(1:nchunks(ifile)) = (/588, 828/)      ! 1899, 1919
                 case default
                    nchunks(ifile) = 2
                    tidx1(1:nchunks(ifile)) = (/  1, 601/)      ! 1850, 1900
                    tidx2(1:nchunks(ifile)) = (/600, 828/)      ! 1899, 1918
                 end select
              case ( 840 )
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/  1, 349/)
                 tidx2(1:nchunks(ifile)) = (/348, 840/)
              case ( 829 )
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/  1, 590/)      ! 1850, 1900
                 tidx2(1:nchunks(ifile)) = (/589, 829/)      ! 1899, 1918
              case ( 900 )
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/  1, 601/)      ! 1850, 1900
                 tidx2(1:nchunks(ifile)) = (/600, 900/)      ! 1899, 1924
              case ( 816 )
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/  1, 337/)      ! 1850, 1900
                 tidx2(1:nchunks(ifile)) = (/336, 816/)      ! 1899, 1924
              case ( 3612,6012,12012 ) ! piControl,past1000,midHolocene: ~50Y chunks
                 nchunks(ifile) = int(ntimes(1,1)/600)
                 tidx1(1) =   1
                 tidx2(1) = 600
                 do ic = 2,nchunks(ifile)
                    tidx1(ic) = tidx2(ic-1) + 1
                    tidx2(ic) = tidx1(ic) + 599
                 enddo
                 tidx2(nchunks(ifile)) = ntimes(1,1)
              case ( 12000 ) ! BGC controls
                 nchunks(ifile) = 10
                 select case(exp(exp_found)%model_id)
                 case ('CESM1-BGC')
                    select case(exp(exp_found)%expt_id)
                    case ('piControl') ! b40.prescribed_carb.001, 0101 - 0600
                       tidx1(1) = 1201
                       tidx2(1) = 1800
                       do ic = 2,nchunks(ifile)
                          tidx1(ic) = tidx2(ic-1) + 1
                          tidx2(ic) = tidx1(ic) + 599
                       enddo
                    case ('esmControl') ! b40.coup_carb.004, 0301 - 0800
                       tidx1(1) = 3601
                       tidx2(1) = 4200
                       do ic = 2,nchunks(ifile)
                          tidx1(ic) = tidx2(ic-1) + 1
                          tidx2(ic) = tidx1(ic) + 599
                       enddo
                    end select
                 end select
              case ( 1680 )
                 nchunks(ifile) = 3
                 tidx1(1:nchunks(ifile)) = (/  1, 601, 1201/)
                 tidx2(1:nchunks(ifile)) = (/600,1200, 1680/)
              case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only, ~50y chunks
                 nchunks(ifile) = 2
                 tidx1(1:nchunks(ifile)) = (/3613,4213/) ! 1850, 1900, 1951
                 tidx2(1:nchunks(ifile)) = (/4212,4824/) ! 1899, 1950, 2005
              case ( 2388,2400 )
                 select case(exp(exp_found)%model_id)
                 case ('CESM1-WACCM')
                    nchunks(ifile) = 4 ! 2 deg but 66 levels
                    tidx1(1:nchunks(ifile)) = (/   1, 649,1249,1849/) !   96,  150, 200, 250
                    tidx2(1:nchunks(ifile)) = (/ 648,1248,1848,ntimes(ifile,1)/) !  149,  199, 249, 295
                 case ('CESM1-CAM5')
                    nchunks(ifile) = 5 ! 1 deg but 30 levels, 40 y chunks
                    tidx1(1:nchunks(ifile)) = (/   1,481, 961,1441,1921/)
                    tidx2(1:nchunks(ifile)) = (/ 480,960,1440,1920,ntimes(ifile,1)/) !  149,  199, 249, 295
                 case default
                    nchunks(ifile) = 4
                    tidx1(1:nchunks(ifile)) = (/   1, 589,1189,1789/)
                    tidx2(1:nchunks(ifile)) = (/ 588,1188,1788,ntimes(ifile,1)/)
                 end select
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat3b)
                    call read_var(myncid(ifile,3),var_info(var_found(ifile,3))%name,psdata)
                    where (abs(psdata) > spval)
                       psdata = spval
                    endwhere
                    where (abs(indat3a) > spval)
                       indat3a = spval
                    endwhere
                    where (abs(indat3b) > spval)
                       indat3b = spval
                    endwhere
                    where ((indat3a /= spval).and.(indat3b /= spval))
                       cmordat3d = indat3a + indat3b
                    elsewhere
                       cmordat3d = spval
                    endwhere
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
                    error_flag = cmor_write(        &
                         var_id        = zfactor_id,&
                         data          = psdata,   &
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
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
              if (allocated(time))      deallocate(time)
              if (allocated(time_bnds)) deallocate(time_bnds)
           enddo
        end select
        if (allocated(psdata))    deallocate(psdata)
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
        error_flag   = 0
        mycmor%positive = ' '
        original_name= ' '
        !
        if (allocated(time))      deallocate(time)
        if (allocated(time_bnds)) deallocate(time_bnds)
        !
        error_flag = cmor_close()
        if (error_flag < 0) then
           write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
        else
           write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
        endif
        call reset_netcdf_var
     endif
  enddo xwalk_loop
end program Amon_CMOR
