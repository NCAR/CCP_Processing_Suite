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
  integer,dimension(100)::found_xw,cmor_var_id
  real,dimension(:,:,:),allocatable::indat3a,indat3b,indat3c,cmordat3d
  double precision,dimension(1)  ::time,tval
  double precision,dimension(2,1)::time_bnds,tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile,cmor_filename
  character(len=256)::fcase,fcomp,fsvar,ftime,histfile
  integer::i,j,k,m,n,it,ivar,jvar,kvar,length,iexp,jexp,ixw,jxw,ilev,ic,tcount,histncid,yrcount,year1,year2
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
  year1 = exp(exp_found)%runbeg
  year2 = exp(exp_found)%runend
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
     kvar = 0
     do n=1,var_counter
        do ixw = 1,num_xw
           if (trim(var_info(n)%name)==trim(xw(ixw)%cesm_vars(1))) then
              kvar = kvar + 1
              var_found(1,ixw) = n
              found_xw(kvar) = ixw
              xw_found = ixw
           endif
        enddo
     enddo
     if (kvar==0) then
        write(*,'(''NEVER FOUND: '',a,'' STOP. '')') trim(xw(ixw)%cesm_vars(1))
        stop
     endif
     call close_cdf(histncid)
  enddo
  !
  tcount = 1
  do ixw = 1,num_xw
     !
     ! Specify path where tables can be found and indicate that existing netCDF files should be overwritten.
     !
     write(logfile,'(''log_cmor.'',a,''.'',a,''_'',a)') &
          trim(mycmor%experiment_id),&
          trim(xw(ixw)%entry),&
          trim(exp(exp_found)%rip_code)
     !
     error_flag = cmor_setup(inpath='CMOR',&
          netcdf_file_action=CMOR_APPEND,&
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
     ! Define axes via 'cmor_axis'
     !
     table_ids(1) = cmor_load_table('Tables/CMIP5_Oyr')
     table_ids(2) = cmor_load_table('Tables/CMIP5_grids')
     table_ids(3) = cmor_load_table('Tables/CMIP5_fx')
     call cmor_set_table(table_ids(2))
     call define_ocn_axes(xw(ixw)%dims)
     call cmor_set_table(table_ids(1))
     !
     ! Add global metadata
     !
     call add_global_metadata
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
     select case (xw(ixw)%entry)
     case ('fbddtalk','fddtalk','intpp','frn','intppico','intpcalc','intpdiat','intpdiaz','intpn2','intpbsi','intpcalcite')
        var_info(var_found(1,ixw))%units = 'mol m-2 s-1'
     case ('talk')
        var_info(var_found(1,ixw))%units = 'mol m-3'
     case ('ph')
        var_info(var_found(1,ixw))%units = '1'
     case ('dpco2','spco2')
        var_info(var_found(1,ixw))%units = 'Pa'
     case ('fgco2')
        var_info(var_found(1,ixw))%units = 'kg m-2 s-1'
        mycmor%positive = 'down'
     case ('intdic')
        var_info(var_found(1,ixw))%units = 'kg m-2'
     case ('hfss','rlds','rsntds','rsds','tauuo','tauvo')
        mycmor%positive = 'up'
     case ('expcalc','expcfe','expc','expsi','fgo2','fsn')
        mycmor%positive = 'down'
     end select
     !
     ! All fields are full column
     cmor_var_id(ixw) = cmor_variable(             &
          table=mycmor%table_file,                           &
          table_entry=xw(ixw)%entry,                         &
          units=var_info(var_found(1,ixw))%units,              &
          axis_ids=(/grid_id(1),axis_ids(3),axis_ids(4)/),   &
          missing_value=var_info(var_found(1,ixw))%missing_value,&
          positive=mycmor%positive,                          &
          original_name=original_name,                       &
          comment=xw(ixw)%comment)
     write(*,'(''cmor_variable name: '',a,'' ixw '',i10,'' var_id '',i10)') trim(xw(ixw)%entry),ixw,cmor_var_id(ixw)
  enddo
  !
  ! Open CESM file and get information(s)
  !
  yrcount_loop: do yrcount = year1,year2
     write(histfile,'(''data/'',a,''.'',a,''.subset.'',i4.4,''.nc'')') trim(case_read),trim(comp_read),yrcount
     call open_cdf(histncid,trim(histfile),.true.)
     call get_dims(histncid)
     call get_vars(histncid)
     time_counter = 1
     call read_var(histncid,'time_bound',time_bnds(:,1))
     time = (time_bnds(1,1)+time_bnds(2,1))/2.
     !
     ! Perform derivations and cycle through time, writing data too
     !
     time_counter = 1
     do ixw = 1,num_xw
        select case (xw(ixw)%entry)
        case ('talk')
           !           !
           ! ALK converted from meq/m3 to mol m-3 via * 1.e-3
           !
           if (allocated(indat3a))   deallocate(indat3a)
           if (allocated(cmordat3d)) deallocate(cmordat3d)
           allocate(indat3a(nlons,nlats,nlevs))
           allocate(cmordat3d(nlons,nlats,nlevs))
           !
           indat3a   = var_info(var_found(1,ixw))%missing_value
           cmordat3d = var_info(var_found(1,ixw))%missing_value
           call read_var(histncid,var_info(var_found(1,ixw))%name,indat3a)
           write(*,'(''read_var : '',a,'' id '',3i10)') trim(var_info(var_found(1,ixw))%name),ixw
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
                var_id        = cmor_var_id(ixw), &
                data          = cmordat3d,   &
                ntimes_passed = 1,           &
                time_vals     = tval,        &
                time_bnds     = tbnd)
           if (error_flag < 0) then
              write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
              stop
           endif
        case ('ph','expcalc','expcfe','expc','expsi')
           !     
           ! pH_3D
           !
           if (allocated(indat3a)) deallocate(indat3a)
           allocate(indat3a(nlons,nlats,nlevs))
           !
           indat3a = var_info(var_found(1,ixw))%missing_value
           call read_var(histncid,var_info(var_found(1,ixw))%name,indat3a)
           write(*,'(''read_var : '',a,'' id '',3i10)') trim(var_info(var_found(1,ixw))%name),ixw
           tval(1) = time(time_counter) ; tbnd(1,1) = time_bnds(1,time_counter) ; tbnd(2,1) = time_bnds(2,time_counter)
           error_flag = cmor_write(          &
                var_id        = cmor_var_id(ixw), &
                data          = indat3a,     &
                ntimes_passed = 1,           &
                time_vals     = tval,        &
                time_bnds     = tbnd)
           if (error_flag < 0) then
              write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
              stop
           endif
        end select
     enddo
     tcount = tcount + 1
     call close_cdf(histncid)
  enddo yrcount_loop
  !
  error_flag = cmor_close()
  if (error_flag < 0) then
     write(*,'(''ERROR close'')')
     stop
  else
     write(*,'('' GOOD close'')')
  endif
end program Oyr_CMOR
