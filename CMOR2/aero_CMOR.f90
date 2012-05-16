program aero_CMOR
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
  !  uninitialized variables used in communicating with CMOR:
  !
  integer::error_flag,cmor_var_id
  real,dimension(:,:)  ,allocatable::indat2a,indat2b,indat2c,indat2d,indat2e,indat2f,cmordat2d,psdata
  real,dimension(:,:,:),allocatable::indat3a,indat3b,indat3c,cmordat3d
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::time_bnds
  double precision,dimension(1)  ::tval
  double precision,dimension(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile,cmor_filename
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,ixw,ilev,ic,jfile
  real::spval
  logical::does_exist
  !
  ! GO!
  !
  mycmor%table_file = 'aero'
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
!!$        select case (xw(ixw)%entry)
!!$        case ('tauu','tauv','hfss','rlut','rlutcs','hfls','rlus','rsus','rsuscs','rsut','rsutcs','mc')
!!$           mycmor%positive = 'up'
!!$        case ('rlds','rldscs','rsds','rsdscs','rsdt','rtmt')
!!$           mycmor%positive = 'down'
!!$        case ('clt','ci','sci')
!!$           var_info(var_found(1,1))%units = '1'
!!$        case ('hurs','cl')
!!$           var_info(var_found(1,1))%units = '%'
!!$        case ('prc','pr','prsn')
!!$           var_info(var_found(1,1))%units = 'kg m-2 s-1'
!!$        end select
        !
        spval=var_info(var_found(1,1))%missing_value
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
        case ('ps')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(2),axis_ids(3),axis_ids(1)/), &
                missing_value=var_info(var_found(1,1))%missing_value,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        case ('cdnc','concaerh2o','concbc','conccmcn','conccn','concdms','concdust','concoa','concso2','concso4','concsoa','concss')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(2),axis_ids(3),axis_ids(4),axis_ids(1)/),  &
                missing_value=var_info(var_found(1,1))%missing_value,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        case ('abs550aer','drybc','drydust','dryoa','drypoa','dryso2','dryso4','dryss',&
             'emibc','emidms','emidust','emiso2','emiso4','emiss',&
             'loadbc','loaddust','loadpoa','loadso4','loadsoa','loadss',&
             'od550aer','reffclwtop','wetbc','wetdust','wetoa','wetso2','wetso4','wetss')
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
        write(*,*) 'called cmor_variable:'
        write(*,*) 'varid         = ',cmor_var_id
        !
        ! Perform derivations and cycle through time, writing data too
        !
        select case (xw(ixw)%entry)
        case ('abs550aer','drydust','dryso2','emidms','emiso4','loadbc','loadpoa','od550aer','ps','reffclwtop','wetso2')
           !
           ! No change
           !
           allocate(indat2a(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (.not.(allocated(time)))      allocate(time(ntimes(1,1)))
           if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
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
              if (ic .le. nchunks(1)) then
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
        case ('drybc','drypoa','emibc','emidust','emiso2','loaddust','loadsoa','wetbc','wetoa')
           !
           ! Two fields summed
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (.not.(allocated(time)))      allocate(time(ntimes(1,1)))
           if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 cmordat2d = spval
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
                 ! 
                 where ((indat2a /= spval).and.(indat2b /= spval))
                    cmordat2d = indat2a + indat2b
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
              if (ic .le. nchunks(1)) then
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
        case ('emiss','loadso4','loadss')
           !
           ! Three fields summed
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (.not.(allocated(time)))      allocate(time(ntimes(1,1)))
           if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 cmordat2d = spval
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat2c)
                 ! 
                 where ((indat2a /= spval).and.(indat2b /= spval).and.(indat2c /= spval))
                    cmordat2d = indat2a + indat2b + indat2c
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
              if (ic .le. nchunks(1)) then
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
        case ('wetdust')
           !
           ! Four fields summed
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats),indat2d(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (.not.(allocated(time)))      allocate(time(ntimes(1,1)))
           if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 cmordat2d = spval
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat2c)
                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat2d)
                 ! 
                 where ((indat2a /= spval).and.(indat2b /= spval).and.&
                        (indat2c /= spval).and.(indat2d /= spval))
                    cmordat2d = indat2a + indat2b + indat2c + indat2d
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
              if (ic .le. nchunks(1)) then
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
        case ('dryoa','dryso4','dryss','wetso4','wetss')
           !
           ! Six fields summed
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats))
           allocate(indat2d(nlons,nlats),indat2e(nlons,nlats),indat2f(nlons,nlats))
           allocate(cmordat2d(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (.not.(allocated(time)))      allocate(time(ntimes(1,1)))
           if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 cmordat2d = spval
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat2c)
                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat2d)
                 call read_var(myncid(1,5),var_info(var_found(1,5))%name,indat2e)
                 call read_var(myncid(1,6),var_info(var_found(1,6))%name,indat2f)
                 ! 
                 where ((indat2a /= spval).and.(indat2b /= spval).and.&
                        (indat2c /= spval).and.(indat2d /= spval).and.&
                        (indat2e /= spval).and.(indat2f /= spval))
                    cmordat2d = indat2a + indat2b + indat2c + indat2d + indat2e + indat2f
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
              if (ic .le. nchunks(1)) then
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
        case ('cdnc','concbc','conccmcn','concdms','concso2')
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nlevs),psdata(nlons,nlats))
           !
           call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
           call get_dims(myncid(1,1))
           call get_vars(myncid(1,1))
           if (.not.(allocated(time)))      allocate(time(ntimes(1,1)))
           if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time_bnds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           !
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,psdata)
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
              if (ic .le. nchunks(1)) then
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
        case ('concdust','concsoa')
           !
           ! Two 3-D fields added together
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs),psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nlevs))
           !
           do ifile = 1,nc_nfiles(1)
              if (.not.(allocated(time)))      allocate(time(ntimes(ifile,1)))
              if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n=1,ntimes(ifile,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
                 time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
              enddo
              !
              nchunks(ifile) = 1
              tidx1(1:nchunks(ifile)) = 1
              tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat3b)
                    call read_var(myncid(ifile,3),var_info(var_found(ifile,3))%name,psdata)
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
        case ('concaerh2o','conccn','concoa','concso4','concss')
           !
           ! Three 3-D fields added together
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs),indat3c(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nlevs))
           !
           do ifile = 1,nc_nfiles(1)
              if (.not.(allocated(time)))      allocate(time(ntimes(ifile,1)))
              if (.not.(allocated(time_bnds))) allocate(time_bnds(2,ntimes(ifile,1)))
              !
              do n=1,ntimes(ifile,1)
                 time_counter = n
                 call read_var(myncid(ifile,1),'time_bnds',time_bnds(:,n))
                 time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
              enddo
              !
              nchunks(ifile) = 1
              tidx1(1:nchunks(ifile)) = 1
              tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat3b)
                    call read_var(myncid(ifile,3),var_info(var_found(ifile,3))%name,indat3c)
                    call read_var(myncid(ifile,4),var_info(var_found(ifile,3))%name,psdata)
                    where ((indat3a /= spval).and.(indat3b /= spval).and.(indat3c /= spval))
                       cmordat3d = indat3a + indat3b + indat3c
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
        if (allocated(indat2a))   deallocate(indat2a)
        if (allocated(indat2b))   deallocate(indat2b)
        if (allocated(indat2c))   deallocate(indat2c)
        if (allocated(indat2c))   deallocate(indat2c)
        if (allocated(indat2d))   deallocate(indat2d)
        if (allocated(indat2e))   deallocate(indat2e)
        if (allocated(indat2f))   deallocate(indat2f)
        if (allocated(cmordat2d)) deallocate(cmordat2d)
        if (allocated(psdata))    deallocate(psdata)
        if (allocated(indat3a))   deallocate(indat3a)
        if (allocated(indat3b))   deallocate(indat3b)
        if (allocated(indat3c))   deallocate(indat3c)
        if (allocated(cmordat3d)) deallocate(cmordat3d)
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
end program aero_CMOR
