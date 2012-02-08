program OImon_CMOR
  ! Convert CCSM4 sea ice monthly (cice.h) data from single-field format
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
  use files_info
  use mycmor_info
  use output_times_info
  !
  implicit none
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  integer::error_flag,cmor_var_id
  real,dimension(:,:),allocatable::indat2anh,indat2ash,indat2bnh,indat2bsh,cmordat
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::time_bnds
  double precision,dimension(1)  ::tval
  double precision,dimension(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile,cmor_filename
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,itab,ixw,ic
  integer,dimension(:),allocatable::i_indices,j_indices
  real::spval
  logical::does_exist
  !
  ! GO!
  !
  mycmor%table_file = 'OImon'
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
  call get_ice_grid
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
           write(*,'(''CHECKING AVAILABILITY OF: '',a,''.'',a,''.'',a,''_(nh|sh).* FILES'')') trim(case_read),trim(comp_read),trim(xw(ixw)%cesm_vars(ivar))
           call build_filenames(case_read,comp_read,xw(ixw)%cesm_vars(ivar),ivar,exp(exp_found)%begyr,exp(exp_found)%endyr,mycmor%table_file)
        endif
     enddo
     !
     ! Open CESM file(s) and get information(s)
     !
     if (all_continue) then
        do ivar = 1,xw(ixw)%ncesm_vars
           !
           ! Get NH data and time values
           !
           call open_cdf(myncidnh(1,ivar),trim(ncfilenh(1,ivar)),.true.)
           write(*,'(''OPENING: '',a80,'' ncid: '',i10)') trim(ncfilenh(1,ivar)),myncidnh(1,ivar)
           call get_dims(myncidnh(1,ivar))
           call get_vars(myncidnh(1,ivar))
           !
           do n=1,dim_counter
              length = len_trim(dim_info(n)%name)
              if(dim_info(n)%name(:length).eq.'time') then
                 ntimes(1,1) = dim_info(n)%length
              endif
           enddo
           call read_att_text(myncidnh(1,ivar),'time','units',time_units)
           !
           do n=1,var_counter
              if (trim(var_info(n)%name) == trim(xw(ixw)%cesm_vars(ivar))) then
                 var_found(1,ivar) = n
              endif
           enddo
           if (var_found(1,ivar) == 0) then
              write(*,*) trim(xw(ixw)%cesm_vars(ivar)),' NEVER FOUND. STOP.'
              stop
           endif
           !
           if (.not.(allocated(time)))      then
              allocate(time(ntimes(1,1)))
           endif
           if (.not.(allocated(time_bnds))) then
              allocate(time_bnds(2,ntimes(1,1)))
           endif
           !
           do n=1,ntimes(1,1)
              time_counter = n
              call read_var(myncidnh(1,ivar),'time_bounds',time_bnds(:,n))
              time(n) = (time_bnds(1,n)+time_bnds(2,n))/2.
           enddo
           !
           ! Get SH data info
           !
           call open_cdf(myncidsh(1,ivar),trim(ncfilesh(1,ivar)),.true.)
           write(*,'(''OPENING: '',a80,'' ncid: '',i10)') trim(ncfilesh(1,ivar)),myncidsh(1,ivar)
           call get_dims(myncidsh(1,ivar))
           call get_vars(myncidsh(1,ivar))
           !
           do n=1,var_counter
              if (trim(var_info(n)%name) == trim(xw(ixw)%cesm_vars(ivar))) then
                 var_found(1,ivar) = n
              endif
           enddo
           if (var_found(1,ivar) == 0) then
              write(*,*) trim(xw(ixw)%cesm_vars(ivar)),' NEVER FOUND. STOP.'
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
        call define_ice_axes(xw(ixw)%dims)
        call cmor_set_table(table_ids(1))
        !
        ! Make manual alterations so that CMOR works. Silly code!
        !
        allocate(indat2anh(nlons,104),indat2ash(nlons,76),cmordat(nlons,384))
        if (xw(ixw)%ncesm_vars == 1) then
           write(original_name,'(a)') xw(ixw)%cesm_vars(1)
        endif
        if (xw(ixw)%ncesm_vars .ge. 2) then
           allocate(indat2bnh(nlons,104),indat2bsh(nlons,76))
           write(original_name,'(a,'','',a)') (trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
        endif
        !
        ! Define unit renamings and/or "positive" as needed
        !
        select case (xw(ixw)%entry)
        case ('hflssi','hfssi','rldssi','rlussi','rsdssi','rsussi','sblsi','ssi','strairx','strairy')
           mycmor%positive = 'up'
        case ('evap')
           mycmor%positive = 'up'
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
        case ('bmelt','grCongel','grFrazil','pr','prsn','snomelt','snoToIce','tmelt')
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
        cmor_var_id = cmor_variable(                &
             table=mycmor%table_file,               &
             table_entry=xw(ixw)%entry,             &
             units=var_info(var_found(1,1))%units,  &
             axis_ids=(/grid_id(1),axis_ids(3)/),   &
             missing_value=spval,                   &
             positive=mycmor%positive,              &
             original_name=original_name,           &
             comment=xw(ixw)%comment)
        !
        ! Cycle through time
        !
        select case (ntimes(1,1))
        case ( 6012 )  ! pre-industrial control, 501 years, 250 and 251 year chunks
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 13
        case ( 1152 ) ! RCP, skip 2005
           nchunks(1) = 1
           tidx1(1:nchunks(1)) = 13
        case default
           nchunks(1)   = 1
           tidx1(1:nchunks(1)) =  1
        end select
        tidx2(nchunks(1)) = ntimes(1,1)
        do ic = 1,nchunks(1)
           do it = tidx1(ic),tidx2(ic)
              time_counter = it
              if (xw(ixw)%ncesm_vars == 1) then
                 call read_var(myncidnh(1,1),var_info(var_found(1,1))%name,indat2anh)
                 call read_var(myncidsh(1,1),var_info(var_found(1,1))%name,indat2ash)
                 cmordat(:,  1: 76) = indat2ash(:,1: 76)
                 cmordat(:,281:384) = indat2anh(:,1:104)
              endif
              !
              ! Perform necessary derivations
              !
              select case (xw(ixw)%entry)
              case ('bmelt','evap','grCongel','grFrazil','pr','prsn','snomelt','snoToIce','tmelt')
                 ! Convert cm day-1 to kg m-2 s-1
                 cmordat = cmordat * (1./86400.) * (1./100.) * 1000.
              end select
              !
              ! Write 'em, Dano!
              !
              tval(1)   = time(it)
              tbnd(1,1) = time_bnds(1,it)
              tbnd(2,1) = time_bnds(2,it)
              error_flag = cmor_write(      &
                   var_id        = cmor_var_id, &
                   data          = cmordat, &
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
        if (allocated(indat2anh)) deallocate(indat2anh)
        if (allocated(indat2bnh)) deallocate(indat2bnh)
        if (allocated(indat2ash)) deallocate(indat2ash)
        if (allocated(indat2bsh)) deallocate(indat2bsh)
        if (allocated(cmordat))   deallocate(cmordat)
        do ivar = 1,xw(ixw)%ncesm_vars
           call close_cdf(myncidnh(1,ivar))
           call close_cdf(myncidsh(1,ivar))
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
           write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
        else
           write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
        endif
        call reset_netcdf_var
     endif
  enddo xwalk_loop
end program OImon_CMOR
