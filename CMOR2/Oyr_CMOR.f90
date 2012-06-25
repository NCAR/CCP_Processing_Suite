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
  integer::error_flag,cmor_var_id
  real,dimension(:,:,:),allocatable::indat3a,indat3b,indat3c,cmordat3d
  double precision,dimension(1)  ::time,tval
  double precision,dimension(2,1)::time_bnds,tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile,cmor_filename
  character(len=256)::fcase,fcomp,fsvar,ftime,histfile
  integer::i,j,k,m,n,it,ivar,jvar,length,iexp,jexp,ixw,ilev,ic,tcount,histncid,yrcount,year1,year2
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
!     write(*,'(''time units in: '',i10,5x,a,5x,a)') histncid,trim(histfile),trim(time_units)
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
  xwalk_loop: do ixw = 1,num_xw
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
     yrcount_loop: do yrcount = year1,year2
        tcount = 1
        write(histfile,'(''data/'',a,''.'',a,''.subset.'',i4.4,''.nc'')') trim(case_read),trim(comp_read),yrcount
        call open_cdf(histncid,trim(histfile),.true.)
        call get_dims(histncid)
        call get_vars(histncid)
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
        write(*,*) 'calendar = ',trim(mycmor%calendar)
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
        select case (xw(ixw)%entry)
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
        ! All fields are full column
        cmor_var_id = cmor_variable(                            &
             table=mycmor%table_file,                           &
             table_entry=xw(ixw)%entry,                         &
             units=var_info(var_found(1,1))%units,              &
             axis_ids=(/grid_id(1),axis_ids(3),axis_ids(4)/),   &
             missing_value=var_info(var_found(1,1))%missing_value,&
             positive=mycmor%positive,                          &
             original_name=original_name,                       &
             comment=xw(ixw)%comment)
!!$        write(*,*) 'cmor_variable: ', &
!!$             trim(mycmor%table_file),'       ',                           &
!!$             trim(xw(ixw)%entry),'       ',                         &
!!$             trim(var_info(var_found(1,1))%units),'       ',              &
!!$             grid_id(1),axis_ids(3),axis_ids(4),   &
!!$             var_info(var_found(1,1))%missing_value,&
!!$             trim(mycmor%positive),'       ',                          &
!!$             trim(original_name),'       ',                       &
!!$             trim(xw(ixw)%comment)
        if (abs(cmor_var_id) .gt. 1000) then
           write(*,'(''Invalid call to cmor_variable, table_entry, varid: '',a,2x,i10)') trim(xw(ixw)%entry),cmor_var_id
           cycle xwalk_loop
        else
           write(*,'(''called cmor_variable, table_entry, varid: '',a,2x,i10)') trim(xw(ixw)%entry),cmor_var_id
        endif
        !
        ! Perform derivations and cycle through time, writing data too
        !
        select case (xw(ixw)%entry)
!!$        case ('chldiat','chldiaz','chlpico','co3','co3satarag','co3satcalc','dfe','tos','ph',&
!!$             'dissic','dissoc','nh4','no3','o2','phycalc','phydiat','phydiaz','phypico',&
!!$             'physi','po4','si','zooc')
!!$           !
!!$           if (.not.(allocated(indat3a)))   allocate(indat3a(nlons,nlats,nlevs))
!!$           if (.not.(allocated(cmordat2d))) allocate(cmordat2d(nlons,nlats))
!!$           tcount = 1
!!$           !
!!$           indat3a = var_info(var_found(ifile,1))%missing_value
!!$           call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$           tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$           error_flag = cmor_write(          &
!!$                var_id        = cmor_var_id, &
!!$                data          = indat3a,   &
!!$                ntimes_passed = 1,           &
!!$                time_vals     = tval,        &
!!$                time_bnds     = tbnd)
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$              stop
!!$           endif
!!$           if (mod(tcount,100) == 0) then
!!$              error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$              if (error_flag < 0) then
!!$                 write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
!!$                 stop
!!$              else
!!$                 write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
!!$              endif
!!$              write(*,'(''DONE writing '',a,'' T# '',2i8,'' chunk# '',i6)') trim(xw(ixw)%entry),tcount,it-1,ic
!!$           enddo
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
                var_id        = cmor_var_id, &
                data          = cmordat3d,   &
                ntimes_passed = 1,           &
                time_vals     = tval,        &
                time_bnds     = tbnd)
           write(*,*) 'cmor_write: ',cmor_var_id,error_flag,tcount
           if (error_flag < 0) then
              write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
              stop
           endif
!!$           if (mod(tcount,100) == 0) then
!!$              error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$              if (error_flag < 0) then
!!$                 write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
!!$                 stop
!!$              else
!!$                 write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
!!$              endif
!!$              write(*,'(''DONE writing '',a,'' T# '',2i8,'' chunk# '',i6)') trim(xw(ixw)%entry),tcount,it-1,ic
!!$           endif
!!$        case ('fbddtalk','fddtalk')
!!$           !
!!$           ! Convert meq/m3 cm/s to mol m-2 s-1 via * 1.e-5
!!$           !
!!$           allocate(indat2a(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60, 1152 ) ! RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 4824 ) ! LGM from 149901 to 190012; want only 1800-1900
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/3613/) ! 1800-01
!!$                 tidx2(1:nchunks(1)) = (/4824/) ! 1900-12
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    indat2a = var_info(var_found(1,1))%missing_value
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat2a)
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (kmt(i,j) .ge. 1) then
!!$                             indat2a(i,j) = indat2a(i,j) * 1.e-5
!!$                          endif
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = indat2a,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$              enddo
!!$           enddo
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('dpco2','spco2')
!!$           !
!!$           ! Convert ppmv to Pa via * 0.101325
!!$           !
!!$           allocate(indat2a(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60, 1152 ) ! RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 4824 ) ! LGM from 149901 to 190012; want only 1800-1900
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/3613/) ! 1800-01
!!$                 tidx2(1:nchunks(1)) = (/4824/) ! 1900-12
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    indat2a = var_info(var_found(1,1))%missing_value
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat2a)
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (kmt(i,j) .ge. 1) then
!!$                             indat2a(i,j) = indat2a(i,j) * 0.101325
!!$                          endif
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = indat2a,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$              enddo
!!$           enddo
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('fgco2')
!!$           !
!!$           ! Convert mmol/m3 cm/s to kg m-2 s-1 via * 12.0e-8
!!$           !
!!$           allocate(indat2a(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60, 1152 ) ! RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 4824 ) ! LGM from 149901 to 190012; want only 1800-1900
!!$                 nchunks(1) = 1
!!$
!!$                 tidx1(1:nchunks(1)) = (/3613/) ! 1800-01
!!$                 tidx2(1:nchunks(1)) = (/4824/) ! 1900-12
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    indat2a = var_info(var_found(1,1))%missing_value
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat2a)
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (kmt(i,j) .ge. 1) then
!!$                             indat2a(i,j) = indat2a(i,j) * 12.0e-8
!!$                          endif
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = indat2a,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$              enddo
!!$           enddo
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('tauuo','tauvo','hfss','pr','prsn','rsds','rsntds','zos','omlmax',&
!!$              'fbddtdic','fbddtdife','fbddtdip','fbddtdisi',&
!!$              'fddtdic','fddtdife','fddtdip','fddtdisi','fgo2',&
!!$              'o2min','zo2min','zsatarag','zsatcalc')
!!$           !
!!$           ! No changes - just pass through
!!$           !
!!$           allocate(indat2a(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60, 1152 ) ! RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 4824 ) ! LGM from 149901 to 190012; want only 1800-1900
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/3613/) ! 1800-01
!!$                 tidx2(1:nchunks(1)) = (/4824/) ! 1900-12
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    indat2a = var_info(var_found(1,1))%missing_value
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat2a)
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = indat2a,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$              enddo
!!$           enddo
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('msftbarot')
!!$           !
!!$           ! msftbarot: Convert BSF from Sv to kg s-1
!!$           !
!!$           if (allocated(indat2a)) deallocate(indat2a(nlons,nlats))
!!$           allocate(indat2a(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60, 1152 ) ! RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 4824 ) ! LGM from 149901 to 190012; want only 1800-1900
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/3613/) ! 1800-01
!!$                 tidx2(1:nchunks(1)) = (/4824/) ! 1900-12
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    indat2a = var_info(var_found(1,1))%missing_value
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat2a)
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (kmt(i,j) .ge. 1) then
!!$                             indat2a(i,j) = indat2a(i,j) * (1.e6 * 1000.) ! 10^6 m3 s-1 to kg s-1
!!$                          endif
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = indat2a,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$                 else
!!$                    write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$                 endif
!!$              enddo
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              dim_counter  = 0
!!$              var_counter  = 0
!!$              time_counter = 0
!!$              file_counter = 0
!!$           enddo
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('thetao','agessc','uo','vo','rhopoto')
!!$           !
!!$           ! No changes
!!$           !
!!$           allocate(indat3a(nlons,nlats,nlevs))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60 ) ! RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 1152 )               ! RCP all in one file from 2005-2100
!!$                 nchunks(1) = 10
!!$                 tidx1(1) =  13
!!$                 tidx2(1) =  60
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$                 tidx2(nchunks(1)) = ntimes(1,1)
!!$              case ( 1872 ) ! 20C from 1850-2005
!!$                 nchunks(1) = 16
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$                 tidx2(nchunks(1)) = ntimes(1,1)
!!$              case ( 2664 ) ! FASTCHEM piControl from 70-291
!!$                 nchunks(1) = 23
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$              case ( 300 ) ! 30 years
!!$                 nchunks(1) = 3
!!$                 tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                 tidx2(1:nchunks(1)) = (/120, 240, 300/)
!!$              case ( 360 ) ! 30 years
!!$                 if ((exp(exp_found)%begyr==2006).and.(exp(exp_found)%endyr==2035)) then
!!$                    nchunks(1) = 4
!!$                    tidx1(1:nchunks(1)) = (/  1,  49, 169, 289/)
!!$                    tidx2(1:nchunks(1)) = (/ 48, 168, 288, 360/)
!!$                 else
!!$                    nchunks(1) = 3
!!$                    tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                    tidx2(1:nchunks(1)) = (/120, 240, 360/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              !
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    indat3a = var_info(var_found(1,1))%missing_value
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = indat3a,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$                 !
!!$                 cmor_filename = ' '
!!$                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$                    stop
!!$                 else
!!$                    write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$                 endif
!!$              enddo
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$           enddo
!!$           dim_counter  = 0
!!$           var_counter  = 0
!!$           time_counter = 0
!!$           file_counter = 0
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('cfc11')
!!$           !
!!$           ! cfc11 - convert fmol/cm^3 to mol kg-1
!!$           !
!!$           allocate(indat3a(nlons,nlats,nlevs))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 60 ) ! RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 1152 )               ! RCP all in one file from 2005-2100
!!$                 nchunks(1) = 10
!!$                 tidx1(1) =  13
!!$                 tidx2(1) =  60
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$                 tidx2(nchunks(1)) = ntimes(1,1)
!!$              case ( 1872 ) ! 20C from 1850-2005
!!$                 nchunks(1) = 16
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$                 tidx2(nchunks(1)) = ntimes(1,1)
!!$              case ( 2664 ) ! FASTCHEM piControl from 70-291
!!$                 nchunks(1) = 23
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$              case ( 300 ) ! 25 years
!!$                 nchunks(1) = 3
!!$                 tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                 tidx2(1:nchunks(1)) = (/120, 240, 300/)
!!$              case ( 360 ) ! 30 years
!!$                 if ((exp(exp_found)%begyr==2006).and.(exp(exp_found)%endyr==2035)) then
!!$                    nchunks(1) = 4
!!$                    tidx1(1:nchunks(1)) = (/  1,  49, 169, 289/)
!!$                    tidx2(1:nchunks(1)) = (/ 48, 168, 288, 360/)
!!$                 else
!!$                    nchunks(1) = 3
!!$                    tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                    tidx2(1:nchunks(1)) = (/120, 240, 360/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    indat3a = var_info(var_found(1,1))%missing_value
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    do k = 1,nlevs
!!$                       do j = 1,nlats
!!$                          do i = 1,nlons
!!$                             if (kmt(i,j) == k) then
!!$                                indat3a(i,j,k) = indat3a(i,j,k) / 1.e12
!!$                             endif
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = indat3a,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$                 !
!!$                 cmor_filename = ' '
!!$                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$                    stop
!!$                 else
!!$                    write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$                 endif
!!$              enddo
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$           enddo
!!$           dim_counter  = 0
!!$           var_counter  = 0
!!$           time_counter = 0
!!$           file_counter = 0
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('wmo','wmosq')
!!$           !
!!$           ! wmo,wmosq: 0.5 * (WVEL(i,j,k) + WVEL(i,j,k+1)) * TAREA(i,j) * rho_0
!!$           !
!!$           allocate(indat3a(nlons,nlats,nlevs),cmordat3d(nlons,nlats,nlevs))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60 )               ! RCP from 2005-2009, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 1152 )               ! RCP all in one file from 2005-2100
!!$                 nchunks(1) = 10
!!$                 tidx1(1) =  13
!!$                 tidx2(1) =  60
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$              case ( 2664 ) ! FASTCHEM piControl from 70-291
!!$                 nchunks(1) = 23
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$                 tidx2(nchunks(1)) = ntimes(1,1)
!!$              case ( 1872 ) ! 20C from 1850-2005
!!$                 nchunks(1) = 16
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$                 tidx2(nchunks(1)) = ntimes(1,1)
!!$              case ( 300 ) ! 30 years
!!$                 nchunks(1) = 3
!!$                 tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                 tidx2(1:nchunks(1)) = (/120, 240, 300/)
!!$              case ( 360 ) ! 30 years
!!$                 if ((exp(exp_found)%begyr==2006).and.(exp(exp_found)%endyr==2035)) then
!!$                    nchunks(1) = 4
!!$                    tidx1(1:nchunks(1)) = (/  1,  49, 169, 289/)
!!$                    tidx2(1:nchunks(1)) = (/ 48, 168, 288, 360/)
!!$                 else
!!$                    nchunks(1) = 3
!!$                    tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                    tidx2(1:nchunks(1)) = (/120, 240, 360/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    indat3a   = var_info(var_found(1,1))%missing_value
!!$                    cmordat3d = var_info(var_found(1,1))%missing_value
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    do k = 1,nlevs-1
!!$                       do j = 1,nlats
!!$                          do i = 1,nlons
!!$                             if (indat3a(i,j,k) /= var_info(var_found(1,1))%missing_value) then
!!$                                cmordat3d(i,j,k) = ((0.5*(indat3a(i,j,k)+indat3a(i,j,k+1)))*ocn_t_area(i,j)*rho_0)/1000. ! g s-1 to kg s-1
!!$                             else
!!$                                cmordat3d(i,j,k) = var_info(var_found(1,1))%missing_value
!!$                             endif
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                    cmordat3d(:,:,nlevs) = 0.
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat3d,   &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$                 !
!!$                 cmor_filename = ' '
!!$                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$                    stop
!!$                 else
!!$                    write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$                 endif
!!$              enddo
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              dim_counter  = 0
!!$              var_counter  = 0
!!$              time_counter = 0
!!$              file_counter = 0
!!$           enddo
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('umo')
!!$           !
!!$           ! umo: 0.5 * (UVEL(i,j) + UVEL(i,j-1)) * HTE(i,j) * dz(k) * rho_0 (rho_0 = 1 g cm-3)
!!$           !
!!$           allocate(indat3a(nlons,nlats,nlevs),cmordat3d(nlons,nlats,nlevs))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60 )               ! RCP from 2005-2009, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 1152 )               ! RCP all in one file from 2005-2100
!!$                 nchunks(1) = 10
!!$                 tidx1(1) =  13
!!$                 tidx2(1) =  60
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$              case ( 1872 ) ! 20C from 1850-2005
!!$                 nchunks(1) = 16
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$                 tidx2(nchunks(1)) = ntimes(1,1)
!!$              case ( 2664 ) ! FASTCHEM piControl from 70-291
!!$                 nchunks(1) = 23
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$              case ( 300 ) ! 30 years
!!$                 nchunks(1) = 3
!!$                 tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                 tidx2(1:nchunks(1)) = (/120, 240, 300/)
!!$              case ( 360 ) ! 30 years
!!$                 if ((exp(exp_found)%begyr==2006).and.(exp(exp_found)%endyr==2035)) then
!!$                    nchunks(1) = 4
!!$                    tidx1(1:nchunks(1)) = (/  1,  49, 169, 289/)
!!$                    tidx2(1:nchunks(1)) = (/ 48, 168, 288, 360/)
!!$                 else
!!$                    nchunks(1) = 3
!!$                    tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                    tidx2(1:nchunks(1)) = (/120, 240, 360/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    indat3a   = var_info(var_found(1,1))%missing_value
!!$                    cmordat3d = var_info(var_found(1,1))%missing_value
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    do k = 1,nlevs
!!$                       do j = 2,nlats
!!$                          do i = 1,nlons
!!$                             if (indat3a(i,j,k) /= var_info(var_found(1,1))%missing_value) then
!!$                                cmordat3d(i,j,k) = ((0.5*(indat3a(i,j,k)+indat3a(i,j-1,k)))*ocn_t_hte(i,j)*ocn_t_dz(k)*rho_0)/1000. ! g s-1 to kg s-1
!!$                             else
!!$                                cmordat3d(i,j,k) = var_info(var_found(1,1))%missing_value
!!$                             endif
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                    cmordat3d(:,1,:) = 0.
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat3d,   &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$                 !
!!$                 cmor_filename = ' '
!!$                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$                    stop
!!$                 else
!!$                    write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$                 endif
!!$              enddo
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              dim_counter  = 0
!!$              var_counter  = 0
!!$              time_counter = 0
!!$              file_counter = 0
!!$           enddo
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('vmo')
!!$           !
!!$           ! vmo: 0.5 * (VVEL(i,j) + VVEL(i-1,j)) * HTN(i,j) * dz(k) * rho_0
!!$           !
!!$           allocate(indat3a(nlons,nlats,nlevs),cmordat3d(nlons,nlats,nlevs))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60 )               ! RCP from 2005-2009, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 1152 )               ! RCP all in one file from 2005-2100
!!$                 nchunks(1) = 10
!!$                 tidx1(1) =  13
!!$                 tidx2(1) =  60
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$                 tidx2(nchunks(1)) = ntimes(1,1)
!!$              case ( 1872 ) ! 20C from 1850-2005
!!$                 nchunks(1) = 16
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$                 tidx2(nchunks(1)) = ntimes(1,1)
!!$              case ( 2664 ) ! FASTCHEM piControl from 70-291
!!$                 nchunks(1) = 23
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$              case ( 300 ) ! 30 years
!!$                 nchunks(1) = 3
!!$                 tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                 tidx2(1:nchunks(1)) = (/120, 240, 300/)
!!$              case ( 360 ) ! 30 years
!!$                 if ((exp(exp_found)%begyr==2006).and.(exp(exp_found)%endyr==2035)) then
!!$                    nchunks(1) = 4
!!$                    tidx1(1:nchunks(1)) = (/  1,  49, 169, 289/)
!!$                    tidx2(1:nchunks(1)) = (/ 48, 168, 288, 360/)
!!$                 else
!!$                    nchunks(1) = 3
!!$                    tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                    tidx2(1:nchunks(1)) = (/120, 240, 360/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    indat3a   = var_info(var_found(1,1))%missing_value
!!$                    cmordat3d = var_info(var_found(1,1))%missing_value
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    do k = 1,nlevs
!!$                       do j = 1,nlats
!!$                          do i = 2,nlons
!!$                             if (indat3a(i,j,k) /= var_info(var_found(1,1))%missing_value) then
!!$                                cmordat3d(i,j,k) = ((0.5*(indat3a(i,j,k)+indat3a(i-1,j,k)))*ocn_t_htn(i,j)*ocn_t_dz(k)*rho_0)/1000. ! g s-1 to kg s-1
!!$                             else
!!$                                cmordat3d(i,j,k) = var_info(var_found(1,1))%missing_value
!!$                             endif
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                    do k = 1,nlevs
!!$                       do j = 1,nlats
!!$                          if (indat3a(1,j,k) /= var_info(var_found(1,1))%missing_value) &
!!$                             cmordat3d(1,j,k) = ((0.5*(indat3a(1,j,k)+indat3a(nlons,j,k)))*ocn_t_htn(1,j)*ocn_t_dz(k)*rho_0)/1000. ! g s-1 to kg s-1
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat3d,   &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$                 !
!!$                 cmor_filename = ' '
!!$                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$                    stop
!!$                 else
!!$                    write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$                 endif
!!$              enddo
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              dim_counter  = 0
!!$              var_counter  = 0
!!$              time_counter = 0
!!$              file_counter = 0
!!$           enddo
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('so')
!!$           !
!!$           ! so: SALT*1000 to handle CMOR change to psu bug
!!$           !
!!$           allocate(indat3a(nlons,nlats,nlevs))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60 )               ! RCP from 2005-2009, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 1152 )               ! RCP all in one file from 2005-2100
!!$                 nchunks(1) = 10
!!$                 tidx1(1) =  13
!!$                 tidx2(1) =  60
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$              case ( 1872 ) ! 20C from 1850-2005
!!$                 nchunks(1) = 16
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$                 tidx2(nchunks(1)) = ntimes(1,1)
!!$              case ( 2664 ) ! FASTCHEM piControl from 70-291
!!$                 nchunks(1) = 23
!!$                 tidx1(1) =   1
!!$                 tidx2(1) = 120
!!$                 do ic = 2,nchunks(1)
!!$                    tidx1(ic) = tidx2(ic-1) + 1
!!$                    tidx2(ic) = tidx1(ic) + 119
!!$                 enddo
!!$              case ( 300 ) ! 15 years
!!$                 nchunks(1) = 3
!!$                 tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                 tidx2(1:nchunks(1)) = (/120, 240, 300/)
!!$              case ( 360 ) ! 30 years
!!$                 if ((exp(exp_found)%begyr==2006).and.(exp(exp_found)%endyr==2035)) then
!!$                    nchunks(1) = 4
!!$                    tidx1(1:nchunks(1)) = (/  1,  49, 169, 289/)
!!$                    tidx2(1:nchunks(1)) = (/ 48, 168, 288, 360/)
!!$                 else
!!$                    nchunks(1) = 3
!!$                    tidx1(1:nchunks(1)) = (/  1, 121, 241/)
!!$                    tidx2(1:nchunks(1)) = (/120, 240, 360/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    indat3a = var_info(var_found(1,1))%missing_value
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    where (indat3a /= var_info(var_found(1,1))%missing_value)
!!$                       indat3a = indat3a*1000.
!!$                    endwhere
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = indat3a,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$                 !
!!$                 cmor_filename = ' '
!!$                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$                    stop
!!$                 else
!!$                    write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$                 endif
!!$              enddo
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$           enddo
!!$           dim_counter  = 0
!!$           var_counter  = 0
!!$           time_counter = 0
!!$           file_counter = 0
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('soga','thetaoga','masso')
!!$           !
!!$           ! soga     - global average salinity
!!$           ! thetaoga - global average temperature
!!$           ! masso    - total mass (kg)
!!$           !
!!$           allocate(indat3a(nlons,nlats,nlevs)) ! SALT or TEMP
!!$           allocate(volume(nlons,nlats,nlevs))  ! Ocean volume (TAREA * dz), cm^3
!!$           volume = 0.
!!$           do j = 1,nlats
!!$              do i = 1,nlons
!!$                 do k = 1,nlevs
!!$                    if (kmt(i,j).ge.k) then
!!$                       volume(i,j,k) = ocn_t_area(i,j)*ocn_t_dz(k)
!!$                    endif
!!$                 enddo
!!$              enddo
!!$           enddo
!!$!           write(*,'(''VOLUME (cm^3): '',3e20.10)') minval(volume),maxval(volume),sum(volume)
!!$!           write(*,'(''VOLUME (km^3): '',3e20.10)') minval(volume)/1.e15,maxval(volume)/1.e15,sum(volume)/1.e15
!!$           !
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$!              write(*,*) 'OPENING: ',histncid,trim(histfile)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              if (allocated(indat1a))    deallocate(indat1a)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              allocate(indat1a(ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 4824 ) ! LGM from 149901 to 190012; want only 1800-1900
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/3613/) ! 1800-01
!!$                 tidx2(1:nchunks(1)) = (/4824/) ! 1900-12
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              indat1a = 0.
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$!                    write(*,'(''TIME: '',i10,10x,f10.3)') it,time(it)
!!$                    time_counter = it
!!$                    indat3a = 1.e30
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    do k = 1,nlevs
!!$                       do j = 1,nlats
!!$                          do i = 1,nlons
!!$                             if (indat3a(i,j,k).lt.1.e30) indat1a(it) = indat1a(it) + (indat3a(i,j,k)*volume(i,j,k))
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                 enddo
!!$              enddo
!!$              !
!!$!              write(*,*) 'MIN: ',minval(indat1a),'MAX: ',maxval(indat1a)
!!$              if (xw(ixw)%entry=='masso') then
!!$                 indat1a = indat1a/1000.
!!$              elseif (xw(ixw)%entry=='thetaoga') then
!!$                 indat1a = indat1a/sum(volume)
!!$              elseif (xw(ixw)%entry=='soga') then
!!$                 indat1a = 1.e6*(indat1a/sum(volume))
!!$              endif
!!$!              write(*,*) 'MIN: ',minval(indat1a),'MAX: ',maxval(indat1a)
!!$              !
!!$              error_flag = cmor_write(                &
!!$                   var_id        = cmor_var_id,                  &
!!$                   data          = indat1a(tidx1(1):tidx2(1)), &
!!$                   ntimes_passed = ((tidx2(1)-tidx1(1))+1),    &
!!$                   time_vals     = time(tidx1(1):tidx2(1)),    &
!!$                   time_bnds     = time_bnds(:,tidx1(1):tidx2(1)))
!!$              if (error_flag < 0) then
!!$                 write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                 stop
!!$              endif
!!$              write(*,'(''DONE writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it-1
!!$              if (allocated(indat1a))   deallocate(indat1a)
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              dim_counter  = 0
!!$              var_counter  = 0
!!$              time_counter = 0
!!$              file_counter = 0
!!$           enddo
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('zosga','zossga','zostoga')
!!$           !
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              write(*,*) 'OPENING: ',histncid,trim(histfile)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              if (allocated(indat1a))    deallocate(indat1a)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              allocate(indat1a(ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 4824 ) ! LGM from 149901 to 190012; want only 1800-1900
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/3613/) ! 1800-01
!!$                 tidx2(1:nchunks(1)) = (/4824/) ! 1900-12
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              !
!!$              indat1a = 0.
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    indat3a = 1.e30
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat1a)
!!$                 enddo
!!$              enddo
!!$              !
!!$              write(*,*) 'MIN: ',minval(indat1a),'MAX: ',maxval(indat1a)
!!$              !
!!$              error_flag = cmor_write(                &
!!$                   var_id        = cmor_var_id,                  &
!!$                   data          = indat1a(tidx1(1):tidx2(1)), &
!!$                   ntimes_passed = ((tidx2(1)-tidx1(1))+1),    &
!!$                   time_vals     = time(tidx1(1):tidx2(1)),    &
!!$                   time_bnds     = time_bnds(:,tidx1(1):tidx2(1)))
!!$              if (error_flag < 0) then
!!$                 write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                 stop
!!$              endif
!!$              write(*,'(''DONE writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it-1
!!$              if (allocated(indat1a))   deallocate(indat1a)
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              dim_counter  = 0
!!$              var_counter  = 0
!!$              time_counter = 0
!!$              file_counter = 0
!!$           enddo
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('rlds')
!!$           !
!!$           ! rlds: LWUP_F + LWDN_F where IFRAC = 0
!!$           !
!!$           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),icefrac(nlons,nlats),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              call open_cdf(myncid(1,2),trim(ncfile(1,2)),.true.)
!!$              call get_dims(myncid(1,2))
!!$              call get_vars(myncid(1,2))
!!$              call open_cdf(myncid(1,3),trim(ncfile(1,3)),.true.)
!!$              call get_dims(myncid(1,3))
!!$              call get_vars(myncid(1,3))
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    cmordat2d = spval
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat2a)
!!$                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
!!$                    call read_var(myncid(1,3),var_info(var_found(1,3))%name,icefrac)
!!$                    where (icefrac == 0)
!!$                       where ((indat2a /= spval).and.(indat2b /= spval))
!!$                          cmordat2d = indat2a + indat2b
!!$                       elsewhere
!!$                          cmordat2d = spval
!!$                       endwhere
!!$                    elsewhere
!!$                       cmordat2d = spval
!!$                    endwhere
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$                 !
!!$                 cmor_filename = ' '
!!$                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$                    stop
!!$                 else
!!$                    write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$                 endif
!!$              enddo
!!$           enddo
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$              call close_cdf(myncid(1,2))
!!$              call close_cdf(myncid(1,3))
!!$           enddo
!!$        case ('fddtdin','fddtdinrlds','fbddtdin')
!!$           !
!!$           ! Add two fields
!!$           !
!!$           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              call open_cdf(myncid(1,2),trim(ncfile(1,2)),.true.)
!!$              call get_dims(myncid(1,2))
!!$              call get_vars(myncid(1,2))
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    cmordat2d = spval
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat2a)
!!$                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (kmt(i,j) .ge. 1) cmordat2d(i,j) = indat2a(i,j) + indat2b(i,j)
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$              enddo
!!$           enddo
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$              call close_cdf(myncid(1,2))
!!$           enddo
!!$        case ('chl','phyc','phyfe')
!!$           !
!!$           ! Add three full-column fields, but use only topmost layer in sum
!!$           !
!!$           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs),indat3c(nlons,nlats,nlevs),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              call open_cdf(myncid(1,2),trim(ncfile(1,2)),.true.)
!!$              call get_dims(myncid(1,2))
!!$              call get_vars(myncid(1,2))
!!$              call open_cdf(myncid(1,3),trim(ncfile(1,3)),.true.)
!!$              call get_dims(myncid(1,3))
!!$              call get_vars(myncid(1,3))
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    cmordat2d = spval
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
!!$                    call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (kmt(i,j) .ge. 1) cmordat2d(i,j) = indat3a(i,j,1) + indat3b(i,j,1) + indat3c(i,j,1)
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$              enddo
!!$           enddo
!!$           !
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$              call close_cdf(myncid(1,2))
!!$              call close_cdf(myncid(1,3))
!!$           enddo
!!$        case ('phyp')
!!$           !
!!$           ! Add three 15-level fields, but use only topmost layer in sum
!!$           ! phyp = (diatC+spC)*.00855+diazC*.002735
!!$           !
!!$           allocate(indat3a(nlons,nlats,15),indat3b(nlons,nlats,15),indat3c(nlons,nlats,15),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              call open_cdf(myncid(1,2),trim(ncfile(1,2)),.true.)
!!$              call get_dims(myncid(1,2))
!!$              call get_vars(myncid(1,2))
!!$              call open_cdf(myncid(1,3),trim(ncfile(1,3)),.true.)
!!$              call get_dims(myncid(1,3))
!!$              call get_vars(myncid(1,3))
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    cmordat2d = spval
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
!!$                    call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (kmt(i,j) .ge. 1) then
!!$                             cmordat2d(i,j) = ((indat3a(i,j,1)+indat3b(i,j,1))*0.00855) + &
!!$                                               (indat3c(i,j,1)*0.002735)
!!$                          endif
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$              enddo
!!$           enddo
!!$           !
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$              call close_cdf(myncid(1,2))
!!$              call close_cdf(myncid(1,3))
!!$           enddo
!!$        case ('phyn')
!!$           !
!!$           ! Add three 15-level fields, but use only topmost layer in sum
!!$           ! phyp = (diatC+spC+diazC)*.137
!!$           !
!!$           allocate(indat3a(nlons,nlats,15),indat3b(nlons,nlats,15),indat3c(nlons,nlats,15),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              call open_cdf(myncid(1,2),trim(ncfile(1,2)),.true.)
!!$              call get_dims(myncid(1,2))
!!$              call get_vars(myncid(1,2))
!!$              call open_cdf(myncid(1,3),trim(ncfile(1,3)),.true.)
!!$              call get_dims(myncid(1,3))
!!$              call get_vars(myncid(1,3))
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    cmordat2d = spval
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
!!$                    call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (kmt(i,j) .ge. 1) then
!!$                             cmordat2d(i,j) = ((indat3a(i,j,1)+indat3b(i,j,1)+indat3c(i,j,1))*0.137)
!!$                          endif
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$              enddo
!!$           enddo
!!$           !
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$              call close_cdf(myncid(1,2))
!!$              call close_cdf(myncid(1,3))
!!$           enddo
!!$        case ('chlcalc')
!!$           !
!!$           ! chlcalc: spChl*(spCaCO3/spC)
!!$           !
!!$           allocate(indat3a(nlons,nlats,15),indat3b(nlons,nlats,15),indat3c(nlons,nlats,15),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              call open_cdf(myncid(1,2),trim(ncfile(1,2)),.true.)
!!$              call get_dims(myncid(1,2))
!!$              call get_vars(myncid(1,2))
!!$              call open_cdf(myncid(1,3),trim(ncfile(1,3)),.true.)
!!$              call get_dims(myncid(1,3))
!!$              call get_vars(myncid(1,3))
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    cmordat2d = merge(0.,spval,kmt.gt.0)
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
!!$                    call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (kmt(i,j) .ge. 1) cmordat2d(i,j) = indat3a(i,j,1) * (indat3b(i,j,1)/indat3c(i,j,1))
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$              enddo
!!$           enddo
!!$           !
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$              call close_cdf(myncid(1,2))
!!$              call close_cdf(myncid(1,3))
!!$           enddo
!!$        case ('epc100','epcalc100','epfe100','epsi100')
!!$           !
!!$           ! Full-column fields, but use only "depth100m" meters; z_t at index 11 - 
!!$           ! "Please sample those fields at 105m. In the model output, these fields are actually 
!!$           ! fluxes into the top of the cell, and the top of the cell whose center is at 105m is at 100m."
!!$           !
!!$           allocate(indat3a(nlons,nlats,nlevs),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    cmordat2d = merge(0.,spval,kmt.gt.0)
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    do k = 1,11
!!$                       do j = 1,nlats
!!$                          do i = 1,nlons
!!$                             if (kmt(i,j).le.k) then
!!$                                cmordat2d(i,j) = cmordat2d(i,j)+(indat3a(i,j,k))
!!$                             endif
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$              enddo
!!$           enddo
!!$           !
!!$           cmor_filename = ' '
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('intdic')
!!$           !
!!$           ! Single full-column fields, integrate over all of z_t
!!$           ! and convert from mmol/m3 to kg m-2 via *12.0e-8
!!$           !
!!$           allocate(indat3a(nlons,nlats,nlevs),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    cmordat2d = merge(0.,spval,kmt.gt.0)
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    do k = 1,nlevs
!!$                       do j = 1,nlats
!!$                          do i = 1,nlons
!!$                             if (kmt(i,j).ge.k) then
!!$                                cmordat2d(i,j) = cmordat2d(i,j) + ((indat3a(i,j,k)*12.0e-8)*ocn_t_dz(k))
!!$                             endif
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$              enddo
!!$              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$           enddo
!!$           !
!!$           cmor_filename = ' '
!!$           error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('frn')
!!$           !
!!$           ! Single full-column field, integrate over all of z_t
!!$           !
!!$           allocate(indat3a(nlons,nlats,nlevs),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    cmordat2d = merge(0.,spval,kmt.gt.0)
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    do k = 1,nlevs
!!$                       do j = 1,nlats
!!$                          do i = 1,nlons
!!$                             if (kmt(i,j).ge.k) then
!!$                                cmordat2d(i,j) = cmordat2d(i,j) + ((-1.0e-5*indat3a(i,j,k))*ocn_t_dz(k))
!!$                             endif
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$              enddo
!!$              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$           enddo
!!$           !
!!$           cmor_filename = ' '
!!$           error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('intpp')
!!$           !
!!$           ! Add three 15-level fields, integrate over z_t_150m, multiply by 1.e-5 to convert units
!!$           !
!!$           allocate(indat3a(nlons,nlats,15),indat3b(nlons,nlats,15),indat3c(nlons,nlats,15),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              call open_cdf(myncid(1,2),trim(ncfile(1,2)),.true.)
!!$              call get_dims(myncid(1,2))
!!$              call get_vars(myncid(1,2))
!!$              call open_cdf(myncid(1,3),trim(ncfile(1,3)),.true.)
!!$              call get_dims(myncid(1,3))
!!$              call get_vars(myncid(1,3))
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    cmordat2d = merge(0.,spval,kmt.gt.0)
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
!!$                    call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
!!$                    do k = 1,15
!!$                       do j = 1,nlats
!!$                          do i = 1,nlons
!!$                             if (kmt(i,j).ge.k) then
!!$                                cmordat2d(i,j) = cmordat2d(i,j) + &
!!$                                     ((((indat3a(i,j,k)*ocn_t_dz(k))+&
!!$                                       (indat3b(i,j,k)*ocn_t_dz(k))+&
!!$                                       (indat3c(i,j,k)*ocn_t_dz(k))))*1.e-5)
!!$                             endif
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$              enddo
!!$              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$           enddo
!!$           !
!!$           cmor_filename = ' '
!!$           error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$              call close_cdf(myncid(1,2))
!!$              call close_cdf(myncid(1,3))
!!$           enddo
!!$        case ('fsn')
!!$           !
!!$           ! Add single-level fields and add another 15-level field, integrate it over z_t_150m,
!!$           ! multiply by 1.e-5 to convert units
!!$           ! NOx_FLUX + NHy_FLUX + (integrate diaz_Nfix over z_t_150m)
!!$           !
!!$           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat3a(nlons,nlats,15),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              call open_cdf(myncid(1,2),trim(ncfile(1,2)),.true.)
!!$              call get_dims(myncid(1,2))
!!$              call get_vars(myncid(1,2))
!!$              call open_cdf(myncid(1,3),trim(ncfile(1,3)),.true.)
!!$              call get_dims(myncid(1,3))
!!$              call get_vars(myncid(1,3))
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              if (1.eq.1) time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    cmordat2d = merge(0.,spval,kmt.gt.0)
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat2a)
!!$                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
!!$                    call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3a)
!!$                    do k = 1,15
!!$                       do j = 1,nlats
!!$                          do i = 1,nlons
!!$                             if (kmt(i,j).ge.k) then
!!$                                cmordat2d(i,j) = cmordat2d(i,j) + ((indat3a(i,j,k)*ocn_t_dz(k))*1.e-5)
!!$                             endif
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (kmt(i,j).ge.1) then
!!$                             cmordat2d(i,j) = cmordat2d(i,j) + ((indat2a(i,j)+indat2b(i,j))*1.e-5)
!!$                          endif
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$              enddo
!!$              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$           enddo
!!$           !
!!$           cmor_filename = ' '
!!$           error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$              call close_cdf(myncid(1,2))
!!$              call close_cdf(myncid(1,3))
!!$           enddo
!!$        case ('intpcalc','intpdiat','intpdiaz','intpn2','intppico')
!!$           !
!!$           ! Integrate over z_t_150m
!!$           !
!!$           if (allocated(indat3a)) deallocate(indat3a)
!!$           allocate(indat3a(nlons,nlats,15))
!!$           if (allocated(cmordat2d)) deallocate(cmordat2d)
!!$           allocate(cmordat2d(nlons,nlats))
!!$           !
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    cmordat2d = merge(0.,1.e20,kmt.gt.0)
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    do k = 1,15
!!$                       do j = 1,nlats
!!$                          do i = 1,nlons
!!$                             if (kmt(i,j).ge.k) then
!!$                                cmordat2d(i,j) = cmordat2d(i,j) + ((1.0e-5*indat3a(i,j,k))*ocn_t_dz(k))
!!$                             endif
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$              enddo
!!$              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$           enddo
!!$           !
!!$           cmor_filename = ' '
!!$           error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('intpcalcite','intpbsi')
!!$           !
!!$           ! Integrate -1*field over z_t_150m
!!$           !
!!$           if (allocated(indat3a)) deallocate(indat3a)
!!$           allocate(indat3a(nlons,nlats,15))
!!$           if (allocated(cmordat2d)) deallocate(cmordat2d)
!!$           allocate(cmordat2d(nlons,nlats))
!!$           !
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    cmordat2d = merge(0.,1.e20,kmt.gt.0)
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat3a)
!!$                    do k = 1,15
!!$                       do j = 1,nlats
!!$                          do i = 1,nlons
!!$                             if (kmt(i,j).ge.k) then
!!$                                cmordat2d(i,j) = cmordat2d(i,j) + ((-1.0e-5*indat3a(i,j,k))*ocn_t_dz(k))
!!$                             endif
!!$                          enddo
!!$                       enddo
!!$                    enddo
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$              enddo
!!$              write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$           enddo
!!$           !
!!$           cmor_filename = ' '
!!$           error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        case ('intpnitrate')
!!$           !
!!$           ! Add three single-level fields
!!$           !
!!$           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats),cmordat2d(nlons,nlats))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              call open_cdf(myncid(1,2),trim(ncfile(1,2)),.true.)
!!$              call get_dims(myncid(1,2))
!!$              call get_vars(myncid(1,2))
!!$              call open_cdf(myncid(1,3),trim(ncfile(1,3)),.true.)
!!$              call get_dims(myncid(1,3))
!!$              call get_vars(myncid(1,3))
!!$              spval=var_info(var_found(1,1))%missing_value
!!$              !
!!$              if (allocated(time))       deallocate(time)
!!$              if (allocated(time_bnds))  deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n=1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              select case (ntimes(1,1))
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)   = 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    cmordat2d = spval
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat2a)
!!$                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
!!$                    call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat2c)
!!$                    where ((indat2a /= spval).and.(indat2b /= spval).and.(indat2c /= spval))
!!$                       cmordat2d = indat2a + indat2b + indat2c
!!$                    elsewhere
!!$                       cmordat2d = spval
!!$                    endwhere
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat2d,     &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$              enddo
!!$           enddo
!!$           !
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR close: '',a)') cmor_filename(1:128)
!!$              stop
!!$           else
!!$              write(*,'('' GOOD close: '',a)') cmor_filename(1:128)
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$              call close_cdf(myncid(1,2))
!!$              call close_cdf(myncid(1,3))
!!$           enddo
!!$        case ('msftmyz')
!!$           !
!!$           ! msftmyz: MOC
!!$           !
!!$           ! moc_comp(1)="Eulerian Mean" 
!!$           ! moc_comp(2)="Eddy-Induced (bolus)" 
!!$           ! moc_comp(3)="Submeso" 
!!$           ! 
!!$           ! transport_reg(1):="Global Ocean - Marginal Seas" 
!!$           ! transport_reg(2):="Atlantic Ocean + Mediterranean Sea + Labrador Sea + GIN Sea + Arctic Ocean + Hudson Bay" 
!!$           !
!!$           ! MOC:coordinates = "lat_aux_grid moc_z moc_components transport_region time"
!!$           !
!!$           !                Y           Z      comp      basin
!!$           allocate(indat4a(nlats_trans,nmoc_z,nmoc_comp,ntrans_reg))
!!$           !
!!$           ! basin 1: 'atlantic_arctic_ocean'
!!$           ! basin 2: 'indian_pacific_ocean'
!!$           ! basin 3: 'global_ocean'
!!$           !
!!$           !                  Y           Z      basin
!!$           allocate(cmordat3d(nlats_trans,nmoc_z,3))
!!$           do 1 = 1,nc_nfiles(1)
!!$              call open_cdf(histncid,trim(histfile),.true.)
!!$              call get_dims(histncid)
!!$              call get_vars(histncid)
!!$              !
!!$              if (allocated(time))      deallocate(time)
!!$              if (allocated(time_bnds)) deallocate(time_bnds)
!!$              allocate(time(ntimes(1,1)))
!!$              allocate(time_bnds(2,ntimes(1,1)))
!!$              !
!!$              do n = 1,ntimes(1,1)
!!$                 time_counter = n
!!$                 call read_var(histncid,'time_bound',time_bnds(:,n))
!!$              enddo
!!$              !
!!$              time_bnds(1,1) = int(time_bnds(1,1))-1
!!$              time = (time_bnds(1,:)+time_bnds(2,:))/2.
!!$              !
!!$              select case (ntimes(1,1))
!!$              case ( 1140, 1872, 2664, 2388, 2400 )
!!$                 nchunks(1)= 1
!!$                 tidx1(1:nchunks(1)) =  1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              case ( 6192 ) ! midHolocene from 080101-131612; want only 1000-1300
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/2389/) ! 1000
!!$                 tidx2(1:nchunks(1)) = (/6000/) ! 1300
!!$              case ( 4824 ) ! LGM from 149901 to 190012; want only 1800-1900
!!$                 nchunks(1) = 1
!!$                 tidx1(1:nchunks(1)) = (/3613/) ! 1800-01
!!$                 tidx2(1:nchunks(1)) = (/4824/) ! 1900-12
!!$              case ( 12012 )
!!$                 nchunks(1)= 2
!!$                 tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                 tidx2(1:nchunks(1)) = (/6000,12012/)
!!$              case ( 12000 ) ! BGC controls
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/1201,4201/)
!!$                    tidx2(1:nchunks(1)) = (/4200,7200/)
!!$                 endif
!!$                 if (trim(case_read)=='b40.coup_carb.004') then       ! Use only 0301-0800
!!$                    nchunks(1)= 2
!!$                    tidx1(1:nchunks(1)) = (/   1, 6001/)
!!$                    tidx2(1:nchunks(1)) = (/6000,12012/)
!!$                 endif
!!$              case ( 120 ) ! Decade, but may need to subset
!!$                 if (trim(case_read)=='b40.prescribed_carb.001') then ! Use only 0101-0600
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0100') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0600') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 elseif (trim(case_read)=='b40.coup_carb.004') then ! Use only 0301-0800                  
!!$                    call parse_ncfile(trim(histfile),fcase,fcomp,fsvar,ftime)
!!$                    if (ftime(1:4) == '0300') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 13
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    elseif (ftime(1:4) == '0800') then
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = 12
!!$                    else
!!$                       nchunks(1) = 1
!!$                       tidx1(1:nchunks(1)) = 1
!!$                       tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                    endif
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case ( 60,1152 ) ! 2005 -> 2009 or 2100 of RCP, skip 2005
!!$                 if (exp(exp_found)%begyr==2005) then
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 13
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 else
!!$                    nchunks(1) = 1
!!$                    tidx1(1:nchunks(1)) = 1
!!$                    tidx2(1:nchunks(1)) = ntimes(1,1)
!!$                 endif
!!$              case default
!!$                 nchunks(1)= 1
!!$                 tidx1(1:nchunks(1)) = 1
!!$                 tidx2(1:nchunks(1)) = ntimes(1,1)
!!$              end select
!!$              write(*,'(''# chunks '',i3,'':'',10((i4,''-'',i4),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
!!$              do ic = 1,nchunks(1)
!!$                 do it = tidx1(ic),tidx2(ic)
!!$                    time_counter = it
!!$                    !
!!$                    call read_var(histncid,var_info(var_found(1,1))%name,indat4a)
!!$                    !
!!$                    cmordat3d = var_info(var_found(1,1))%missing_value
!!$                    cmordat3d(:,:,2) = var_info(var_found(1,1))%missing_value ! Indo-Pacific not supplied
!!$                    !
!!$                    ! Convert from Sv to kg s-1
!!$                    !
!!$                    where (indat4a(:,:,1,2) /= 0.)
!!$                       cmordat3d(:,:,1) = indat4a(:,:,1,2)*1.e6*1.e3
!!$                    elsewhere
!!$                       cmordat3d(:,:,1) = var_info(var_found(1,1))%missing_value
!!$                    endwhere
!!$                    where (indat4a(:,:,1,1) /= 0.)
!!$                       cmordat3d(:,:,3) = indat4a(:,:,1,1)*1.e6*1.e3
!!$                    elsewhere
!!$                       cmordat3d(:,:,3) = var_info(var_found(1,1))%missing_value
!!$                    endwhere
!!$                    !
!!$                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
!!$                    error_flag = cmor_write(          &
!!$                         var_id        = cmor_var_id, &
!!$                         data          = cmordat3d,   &
!!$                         ntimes_passed = 1,           &
!!$                         time_vals     = tval,        &
!!$                         time_bnds     = tbnd)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
!!$                       stop
!!$                    endif
!!$                 enddo
!!$                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
!!$                 !
!!$                 if (ic < nchunks(1)) then
!!$                    cmor_filename(1:) = ' '
!!$                    error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
!!$                    if (error_flag < 0) then
!!$                       write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
!!$                       stop
!!$                    else
!!$                       write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
!!$                    endif
!!$                 endif
!!$              enddo
!!$           enddo
!!$           dim_counter  = 0
!!$           var_counter  = 0
!!$           time_counter = 0
!!$           file_counter = 0
!!$           error_flag = cmor_close()
!!$           if (error_flag < 0) then
!!$              write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           else
!!$              write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
!!$           endif
!!$           do 1 = 1,nc_nfiles(1)
!!$              call close_cdf(histncid)
!!$           enddo
!!$        end select
!!$        if (allocated(indat2a))   deallocate(indat2a)
!!$        if (allocated(indat2b))   deallocate(indat2b)
!!$        if (allocated(indat2c))   deallocate(indat2c)
!!$        if (allocated(indat4a))   deallocate(indat4a)
!!$        if (allocated(cmordat2d)) deallocate(cmordat2d)
!!$        if (allocated(indat3a))   deallocate(indat3a)
!!$        if (allocated(indat3b))   deallocate(indat3b)
           !
           ! Reset
           !
           error_flag   = 0
           mycmor%positive = ' '
           original_name= ' '
        end select
     enddo yrcount_loop
     tcount = tcount + 1
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
