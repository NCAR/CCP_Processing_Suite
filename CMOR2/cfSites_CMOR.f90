program cfSites_CMOR
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
  real,dimension(:)  ,allocatable::indat1a,indat1b,indat1c,cmordat1d,psdata
  real,dimension(:,:),allocatable::indat2a,indat2b,indat2c,cmordat2d,work2da,work2db
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::time_bnds
  double precision,dimension(1)  ::tval
  double precision,dimension(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile,cmor_filename
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,ixw,ilev,ic,klo,khi
  real::spval,dataval
  logical::does_exist
  !
  ! GO!
  !
  mycmor%table_file = 'cfSites'
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
           call build_filenames(case_read,comp_read,xw(ixw)%cesm_vars(ivar),ivar,exp(exp_found)%runbeg,exp(exp_found)%runend,mycmor%table_file)
        endif
     enddo
     !
     ! Open CESM file(s) and get information(s)
     !
     if (all_continue) then
        do ivar = 1,xw(ixw)%ncesm_vars
           call open_cdf(myncid(1,ivar),trim(ncfile(1,ivar)),.true.)
           write(*,'(''OPENING: '',a,'' myncid: '',i10)') trim(ncfile(1,ivar)),myncid(1,ivar)
           call get_dims(myncid(1,ivar))
           call get_vars(myncid(1,ivar))
           !
           do n=1,dim_counter
              length = len_trim(dim_info(n)%name)
              if(dim_info(n)%name(:length).eq.'time') then
                 ntimes(1,ivar) = dim_info(n)%length
              endif
           enddo
           call read_att_text(myncid(1,ivar),'time','units',time_units)
           write(*,'(''time length FROM: '',a,'' myncid: '',i10,'' NT: '',i10)') trim(ncfile(1,ivar)),myncid(1,ivar),ntimes(1,ivar)
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
        case ('ccb')
           var_info(var_found(1,1))%units = 'Pa'
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
        write(*,*) 'units         = ',var_info(var_found(1,1))%units(1:20)
!        write(*,*) 'axis_ids      = ',axis_ids(1:naxes)
        write(*,*) 'missing_value = ',var_info(var_found(1,1))%missing_value
        write(*,*) 'positive      = ',trim(mycmor%positive)
        write(*,*) 'original_name = ',trim(original_name)
        !
        select case (xw(ixw)%entry)
!!$        case ('ps')
!!$           cmor_var_id = cmor_variable(                            &
!!$                table=mycmor%table_file,                           &
!!$                table_entry=xw(ixw)%entry,                         &
!!$                units=var_info(var_found(1,1))%units,                &
!!$                axis_ids=(/grid_id(1),axis_ids(1)/), &
!!$                missing_value=var_info(var_found(1,1))%missing_value,&
!!$                positive=mycmor%positive,                          &
!!$                original_name=original_name,                       &
!!$                comment=xw(ixw)%comment)
        case ('ta','ua','va','hus','hur','wap','zg','clw','cli','cl','mc')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/grid_id(1),axis_ids(2),axis_ids(1)/), &
!                axis_ids=(/axis_ids(2),axis_ids(3),axis_ids(1)/), &
                missing_value=var_info(var_found(1,1))%missing_value,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
           write(*,*) 'axis_ids: ',grid_id(1),axis_ids(2),axis_ids(1)
        case default
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/grid_id(1),axis_ids(2)/),  &
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
        case ('ccb','cct','clivi','clwvi','evspsbl','hfls','hfss','hurs','huss',&
              'prw','psl','ps','rldscs','rlds','rlutcs','rlut','rsdscs','rsds','rsdt',&
              'sci','tas','tasmax','tasmin','tauu','tauv','ts','ci','clt')
           !
           ! No change
           !
           if (allocated(time)) deallocate(time)
           allocate(time(ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time',time(n))
           enddo
           !
           ! Determine amount of data to write, to keep close to ~2 GB limit
           !
           if (allocated(indat1a)) deallocate(indat1a)
           allocate(indat1a(nsites))
           if (allocated(indat2a)) deallocate(indat2a)
           allocate(indat2a(nsites,ntimes(1,1)))
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat1a)
                 indat2a(:,it) = indat1a
              enddo
           enddo
!           write(*,'(''0 '',12f8.3)') indat2a(1,1:12)
           indat2a(:,2) = indat2a(:,1)
           do k = 3,it-1,2
              klo = k-1
              khi = k+1
              if (khi.gt.tidx2(nchunks(1))) khi=tidx2(nchunks(1))
!              write(*,*) klo,k,khi
              indat2a(:,k) = (indat2a(:,klo)+indat2a(:,khi))/2.
           enddo
!           write(*,'(''1 '',12f8.3)') indat2a(1,1:12)
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 indat1a = indat2a(:,it)
                 tval(1) = time(it)
                 error_flag = cmor_write(          &
                      var_id        = cmor_var_id, &
                      data          = indat1a,     &
                      ntimes_passed = 1,           &
                      time_vals     = tval)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              enddo
              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           endif
        case ('prc')
           !
           ! prc : PRECC, unit change from m s-1 to kg m-2 s-1
           !
           if (allocated(time)) deallocate(time)
           allocate(time(ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time',time(n))
           enddo
           !
           ! Determine amount of data to write, to keep close to ~2 GB limit
           !
           if (allocated(indat1a)) deallocate(indat1a)
           allocate(indat1a(nsites))
           if (allocated(cmordat1d)) deallocate(cmordat1d)
           allocate(cmordat1d(nsites))
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat1a)
                 cmordat1d = spval
                 where (indat1a /= spval)
                    cmordat1d = indat1a*1000.
                 elsewhere
                    cmordat1d = spval
                 endwhere
                 tval(1) = time(it)
                 error_flag = cmor_write(          &
                      var_id        = cmor_var_id, &
                      data          = indat1a,     &
                      ntimes_passed = 1,           &
                      time_vals     = tval)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              enddo
              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           endif
        case ('pr','prsn')
           !
           ! pr  : Add PRECC + PRECL  , unit change from m s-1 to kg m-2 s-1
           ! prsn: Add PRECSC + PRECSL, unit change from m s-1 to kg m-2 s-1
           !
           if (allocated(time)) deallocate(time)
           allocate(time(ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time',time(n))
           enddo
           !
           if (allocated(indat1a)) deallocate(indat1a)
           allocate(indat1a(nsites))
           if (allocated(indat1b)) deallocate(indat1b)
           allocate(indat1b(nsites))
           if (allocated(cmordat1d)) deallocate(cmordat1d)
           allocate(cmordat1d(nsites))
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat1a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat1b)
                 where ((indat1a /= spval).and.(indat1b /= spval))
                    cmordat1d = (indat1a + indat1b)*1000.
                 elsewhere
                    cmordat1d = spval
                 endwhere
                 tval(1) = time(it)
                 error_flag = cmor_write(          &
                      var_id        = cmor_var_id, &
                      data          = indat1a,     &
                      ntimes_passed = 1,           &
                      time_vals     = tval)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              enddo
              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           endif
        case ('rlus')
           !
           ! Add 
           !
           if (allocated(time)) deallocate(time)
           allocate(time(ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time',time(n))
           enddo
           !
           if (allocated(indat1a)) deallocate(indat1a)
           allocate(indat1a(nsites))
           if (allocated(indat1b)) deallocate(indat1b)
           allocate(indat1b(nsites))
           if (allocated(cmordat1d)) deallocate(cmordat1d)
           allocate(cmordat1d(nsites))
           if (allocated(indat2a)) deallocate(indat2a)
           allocate(indat2a(nsites,ntimes(1,1)))
           if (allocated(indat2b)) deallocate(indat2b)
           allocate(indat2b(nsites,ntimes(1,1)))
           !
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat1a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat1b)
                 indat2a(:,it) = indat1a
                 indat2b(:,it) = indat1b
              enddo
           enddo
!           write(*,'(''0 2 '',119f8.3)') indat2a(:,2)
!           write(*,'(''0 3 '',119f8.3)') indat2a(:,3)
           indat2a(:,1) = indat2a(:,2)
           indat2b(:,1) = indat2b(:,2)
           do k = 3,it-1,2
              indat2a(:,k) = (indat2a(:,k-1)+indat2a(:,k+1))/2.
              indat2b(:,k) = (indat2b(:,k-1)+indat2b(:,k+1))/2.
           enddo
!           write(*,'(''1 2 '',119f8.3)') indat2a(:,2)
!           write(*,'(''1 3 '',119f8.3)') indat2a(:,3)
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 cmordat1d = (indat2a(:,it)+indat2b(:,it))
                 tval(1) = time(it)
                 error_flag = cmor_write(          &
                      var_id        = cmor_var_id, &
                      data          = cmordat1d,     &
                      ntimes_passed = 1,           &
                      time_vals     = tval)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              enddo
              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           endif
        case ('rsus','rsuscs','rsut','rsutcs','rtmt')
           !
           ! rsus   : FSDS  - FSNS
           ! rsuscs : FSDSC - FSNSC
           ! rsut   : SOLIN - FSNTOA
           ! rsutcs : SOLIN - FSNTOAC
           ! rtmt   : FSNT  - FLNT
           !
           if (allocated(time)) deallocate(time)
           allocate(time(ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time',time(n))
           enddo
           !
           if (allocated(indat1a)) deallocate(indat1a)
           allocate(indat1a(nsites))
           if (allocated(indat1b)) deallocate(indat1b)
           allocate(indat1b(nsites))
           if (allocated(cmordat1d)) deallocate(cmordat1d)
           allocate(cmordat1d(nsites))
           if (allocated(indat2a)) deallocate(indat2a)
           allocate(indat2a(nsites,ntimes(1,1)))
           if (allocated(indat2b)) deallocate(indat2b)
           allocate(indat2b(nsites,ntimes(1,1)))
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat1a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat1b)
                 indat2a(:,it) = indat1a
                 indat2b(:,it) = indat1b
              enddo
           enddo
!           write(*,'(''0 2 '',119f8.3)') indat2a(:,2)
!           write(*,'(''0 3 '',119f8.3)') indat2a(:,3)
           indat2a(:,1) = indat2a(:,2)
           indat2b(:,1) = indat2b(:,2)
           do k = 3,it-1,2
              indat2a(:,k) = (indat2a(:,k-1)+indat2a(:,k+1))/2.
              indat2b(:,k) = (indat2b(:,k-1)+indat2b(:,k+1))/2.
           enddo
!           write(*,'(''1 2 '',119f8.3)') indat2a(:,2)
!           write(*,'(''1 3 '',119f8.3)') indat2a(:,3)
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 cmordat1d = (indat2a(:,it)-indat2b(:,it))
                 tval(1) = time(it)
                 error_flag = cmor_write(          &
                      var_id        = cmor_var_id, &
                      data          = cmordat1d,     &
                      ntimes_passed = 1,           &
                      time_vals     = tval)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              enddo
              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           endif
        case ('ta','ua','va','hus','hur','wap','zg','cl')
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           if (allocated(time)) deallocate(time)
           allocate(time(ntimes(1,1)))
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncid(1,1),'time',time(n))
           enddo
           !
           if (allocated(indat2a)) deallocate(indat2a)
           allocate(indat2a(nsites,nlevs))
           if (allocated(indat1b)) deallocate(indat1b)
           allocate(indat1b(nsites))
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 1
           tidx2(1:nchunks(1)) = ntimes(1,1)
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat1b)
                 tval(1) = time(it)
                 error_flag = cmor_write(          &
                      var_id        = cmor_var_id, &
                      data          = indat2a,     &
                      ntimes_passed = 1,           &
                      time_vals     = tval)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
                 error_flag = cmor_write(        &
                      var_id        = zfactor_id,&
                      data          = indat1b,   &
                      ntimes_passed = 1,         &
                      time_vals     = tval,      &
                      store_with    = cmor_var_id)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              enddo
              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
           enddo
           error_flag = cmor_close()
           if (error_flag < 0) then
              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           else
              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
           endif
        end select
        if (allocated(indat1a))   deallocate(indat1a)
        if (allocated(indat1b))   deallocate(indat1b)
        if (allocated(indat1c))   deallocate(indat1c)
        if (allocated(cmordat1d)) deallocate(cmordat1d)
        if (allocated(indat2a))   deallocate(indat2a)
        if (allocated(indat2b))   deallocate(indat2b)
        if (allocated(work2da))   deallocate(work2da)
        if (allocated(work2db))   deallocate(work2db)
     endif
  enddo xwalk_loop
end program cfSites_CMOR
