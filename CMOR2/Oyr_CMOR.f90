program Oyr_CMOR
  ! Convert CCSM4 ocn monthly (pop.h) data from history format to CMOR-compliant format
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
  real,parameter::rho_0 = 1.0 ! 1 g/cm^3 instead of PD since we have a Boussinesq model
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  integer::error_flag
  integer,dimension(100)::cmor_var_id
  real,dimension(:,:,:),allocatable::indat3a,indat3b,indat3c,cmordat3d
  double precision,dimension(1)  ::time,tval
  double precision,dimension(2,1)::time_bnds,tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile,cmor_filename
  character(len=256)::fcase,fcomp,fsvar,ftime,histfile
  integer::i,j,k,m,n,it,ivar,jvar,kvar,length,iexp,jexp,ixw,ilev,ic,tcount,histncid,yrcount,year1,year2
  real::spval
  logical::does_exist
  !
  ! GO!
  !
  mycmor%table_file = 'Oyr'
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
  ! Open CESM file and get information(s)
  !
  year1 = exp(exp_found)%begyr
  year2 = exp(exp_found)%endyr
  if (year1 == 2005) year1 = 2006
  if (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
     year1 = 301
     year2 = 800
  endif
  if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0301-0800                  
     year1 = 101
     year2 = 600
  endif
  !
  do yrcount = year1,year2
     write(histfile,'(''data/'',a,''.'',a,''.subset.'',i4.4,''.nc'')') trim(case_read),trim(comp_read),yrcount
     call open_cdf(histncid,trim(histfile),.true.)
     call get_dims(histncid)
     call get_vars(histncid)
     call read_att_text(histncid,'time','units',time_units)
     do n=1,var_counter
        do ixw = 1,num_xw
           if (trim(var_info(n)%name)==trim(xw(ixw)%cesm_vars(1))) then
              var_found(1,1) = n
              xw_found = ixw
           endif
        enddo
     enddo
     if (var_found(1,1)==0) then
        write(*,'(''NEVER FOUND: '',a,'' STOP. '')') trim(xw(ixw)%cesm_vars(1))
        stop
     endif
     call close_cdf(histncid)
  enddo
  !
  tcount = 1
  xwalk_loop: do ixw = 1,num_xw
     !
     ! Specify path where tables can be found and indicate that existing netCDF files should be overwritten.
     !
     write(logfile,'(''log_cmor.'',a,''.'',a,''_'',a)') &
          trim(mycmor%experiment_id),&
          trim(xw(ixw)%entry),&
          trim(exp(exp_found)%rip_code)
     !
     ! The meaty part
     !
     if (xw(ixw)%ncesm_vars==0) then
        write(*,'(a,'' is UNAVAILABLE.'')') trim(xw(ixw)%entry)
     endif
     !
     do ivar = 1,xw(ixw)%ncesm_vars
        if (trim(xw(ixw)%cesm_vars(ivar))=='UNKNOWN') then
           write(*,'(a,'' has UNKNOWN equivalence.'')') trim(xw(ixw)%entry)
           xw(ixw)%ncesm_vars = 0
        endif
     enddo
     !
     error_flag = cmor_setup(inpath='CMOR',&
          netcdf_file_action=CMOR_APPEND,&
          logfile=logfile)
     !
     ! Define axes via 'cmor_axis'
     !
     table_ids(1) = cmor_load_table('Tables/CMIP5_Oyr')
     table_ids(2) = cmor_load_table('Tables/CMIP5_grids')
     table_ids(3) = cmor_load_table('Tables/CMIP5_fx')
     call cmor_set_table(table_ids(2))
     call define_ocn_axes(xw(ixw)%dims)
     call cmor_set_table(table_ids(1))
     !
     ! Open CESM file and get information(s)
     !
     kvar   = 1
     yrcount_loop: do yrcount = year1,year2
        write(histfile,'(''data/'',a,''.'',a,''.subset.'',i4.4,''.nc'')') trim(case_read),trim(comp_read),yrcount
        call open_cdf(histncid,trim(histfile),.true.)
        call get_dims(histncid)
        call get_vars(histncid)
!        write(*,'(''Reading from '',a)') trim(histfile)
        time_counter = 1
        call read_var(histncid,'time_bound',time_bnds(:,1))
        time = (time_bnds(1,1)+time_bnds(2,1))/2.
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
!!$        table_ids(2) = cmor_load_table('Tables/CMIP5_grids')
!!$        table_ids(3) = cmor_load_table('Tables/CMIP5_fx')
!!$        call cmor_set_table(table_ids(2))
!!$        call define_ocn_axes(xw(ixw)%dims)
!!$        call cmor_set_table(table_ids(1))
        ! 
        ! Make manual alterations so that CMOR works. Silly code!
        !
        if (xw(ixw)%ncesm_vars==1) then
           write(original_name,'(a)') xw(ixw)%cesm_vars(1)
        endif
        if (xw(ixw)%ncesm_vars==2) then
           write(original_name,'(a,'','',a)') (trim(xw(ixw)%cesm_vars(jvar)),jvar=1,xw(ixw)%ncesm_vars)
        endif
        if (xw(ixw)%ncesm_vars==3) then
           write(original_name,'(a,'','',a,'','',a)') (trim(xw(ixw)%cesm_vars(jvar)),jvar=1,xw(ixw)%ncesm_vars)
        endif
        !
        ! Modify units as necessary to accomodate udunits' inability to convert 
        !
        select case (xw(kvar)%entry)
        case ('msftmyz','msftbarot','wmo','umo','vmo')
           var_info(var_found(1,1))%units = 'kg s-1'
        case ('wmosq')
           var_info(var_found(1,1))%units = 'kg2 s-2'
        case ('masso')
           var_info(var_found(1,1))%units = 'kg'
        case ('cfc11')
           var_info(var_found(1,1))%units = 'mol kg-1'
        case ('fbddtalk','fddtalk','intpp','frn','intppico','intpcalc','intpdiat','intpdiaz','intpn2','intpbsi','intpcalcite')
           var_info(var_found(1,1))%units = 'mol m-2 s-1'
        case ('talk')
           var_info(var_found(1,1))%units = 'mol m-3'
        case ('ph')
           var_info(var_found(1,1))%units = '1'
        case ('dpco2','spco2')
           var_info(var_found(1,1))%units = 'Pa'
        case ('fgco2')
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
           mycmor%positive = 'down'
        case ('intdic')
           var_info(var_found(1,1))%units = 'kg m-2'
        case ('hfss','rlds','rsntds','rsds','tauuo','tauvo')
           mycmor%positive = 'up'
        case ('epc100','epcalc100','epfe100','epsi100','fgo2','fsn')
           mycmor%positive = 'down'
        end select
        !
!!$        write(*,*) 'calling cmor_variable:'
!!$        write(*,*) 'table         = ',trim(mycmor%table_file)
!!$        write(*,*) 'table_entry   = ',trim(xw(kvar)%entry)
!!$        write(*,*) 'dimensions    = ',trim(xw(kvar)%dims)
!!$        write(*,*) 'units         = ',trim(var_info(var_found(1,1))%units)
!!$        write(*,*) 'axis_ids      = ',axis_ids(1:naxes)
!!$        write(*,*) 'missing_value = ',var_info(var_found(1,1))%missing_value
!!$        write(*,*) 'positive      = ',trim(mycmor%positive)
!!$        write(*,*) 'original_name = ',trim(original_name)
        !
        ! All fields are full column
        cmor_var_id(kvar) = cmor_variable(                            &
             table=mycmor%table_file,                           &
             table_entry=xw(ixw)%entry,                         &
             units=var_info(var_found(1,1))%units,              &
             axis_ids=(/grid_id(1),axis_ids(3),axis_ids(4)/),   &
             missing_value=var_info(var_found(1,1))%missing_value,&
             positive=mycmor%positive,                          &
             original_name=original_name,                       &
             comment=xw(kvar)%comment)
        write(*,'(''cmor_variable name: '',a,'' ixw '',i10,'' var_id '',i10)') trim(xw(ixw)%entry),ixw,cmor_var_id(kvar)
        if (abs(cmor_var_id(kvar)) .gt. 1000) then
           write(*,'(''Invalid call to cmor_variable, table_entry, varid: '',a,2x,i10)') trim(xw(ixw)%entry),cmor_var_id(kvar)
           cycle xwalk_loop
!!$        else
!!$           write(*,'(''called cmor_variable, table_entry, varid: '',a,2x,i10)') trim(xw(kvar)%entry),cmor_var_id(kvar)
        endif
        !
        ! Perform derivations and cycle through time, writing data too
        !
        select case (xw(ixw)%entry)
        case ('talk')
           !           !
           ! ALK converted from meq/m3 to mol m-3 via * 1.e-3
           !
           if (.not.(allocated(indat3a)))   allocate(indat3a(nlons,nlats,nlevs))
           if (.not.(allocated(cmordat3d))) allocate(cmordat3d(nlons,nlats,nlevs))
           !
           indat3a   = var_info(var_found(1,1))%missing_value
           cmordat3d = var_info(var_found(1,1))%missing_value
           time_counter = 1
           call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
           do k = 1,nlevs
              do j = 1,nlats
                 do i = 1,nlons
                    if (kmt(i,j) .ge. k) then
                       cmordat3d(i,j,k) = indat3a(i,j,k)*1.e-3
                    endif
                 enddo
              enddo
           enddo
           tval(1) = time(time_counter) ; tbnd(1,1) = time_bnds(1,time_counter) ; tbnd(2,1) = time_bnds(2,time_counter)
           error_flag = cmor_write(          &
                var_id        = cmor_var_id(kvar), &
                data          = cmordat3d,   &
                ntimes_passed = 1,           &
                time_vals     = tval,        &
                time_bnds     = tbnd)
           write(*,'(''cmor_write id '',i5,'' # '',i5,'' flag '',i5,'' T '',i5,'' t '',f12.3)') cmor_var_id(kvar),kvar,error_flag,tcount,tval(1)
           if (error_flag < 0) then
              write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
              stop
           endif
        end select
        call close_cdf(histncid)
        tcount = tcount + 1
     enddo yrcount_loop
!     kvar = kvar + 1
  enddo xwalk_loop
  !
  error_flag = cmor_close()
  if (error_flag < 0) then
     write(*,'(''ERROR close'')')
     stop
  else
     write(*,'('' GOOD close'')')
  endif
end program Oyr_CMOR
