program A2DM_CMOR
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
  real,parameter::spval  = 1.e20
  real,parameter::pi     = 4.*(atan(1.))
  real,parameter::grav   = 9.81         ! gravity, m s-1
  real,parameter::rearth = 6.3712e6     ! earth radius, m
  real,parameter::rdair  = 287.04       ! Dry air gas constant, J k-1 kg-1
  real,parameter::cscale = 1.4          ! Carbon scaling factor
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  integer::error_flag,cmor_var_id
  real,dimension(:)    ,allocatable::area_wt
  real,dimension(:,:)  ,allocatable::indat2a,indat2b,indat2c,indat2d,indat2e,indat2f,indat2g,cmordat2d
  real,dimension(:,:)  ,allocatable::psdata,tsdata,rho
  real,dimension(:,:,:),allocatable::indat3a,indat3b,indat3c,indat3d,indat3e,indat3f,indat3g,cmordat3d,work3da,work3db
  real,dimension(:,:,:),allocatable::pshybrid,psdelta
  !
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
  mycmor%table_file = '2D-M'
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
  ! Compute area weights
  !
  if (allocated(area_wt)) deallocate(area_wt)
  allocate(area_wt(nlats))
  area_wt = gaussian_wts*(2.*pi/nlons)*(rearth**2)
  write(*,*) 'EARTH AREA   : ',2*pi*rearth**2
  write(*,*) 'AREA         : ',area_wt
  write(*,*) 'TOTAL AREA WT: ',sum(area_wt)
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
        if (xw(ixw)%ncesm_vars == 2) write(original_name,'(a,'','',a)')    (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars == 3) write(original_name,'(2(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars == 4) write(original_name,'(3(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars == 5) write(original_name,'(4(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars == 6) write(original_name,'(5(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars == 7) write(original_name,'(6(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        !
        ! Modify units as necessary to accomodate udunits' inability to convert 
        !
        select case (xw(ixw)%entry)
        case ('tauu','tauv','hfss','rlut','rlutcs','hfls','rlus','rsus','rsuscs','rsut','rsutcs','mc',&
              'emidms','emiso2','emiss','emibc','emioa')
           mycmor%positive = 'up'
        case ('drydms','drydust','drypoa','dryso2','dryso4','drysoa','dryss',&
              'emiso4','wetdust','wetoa','wetso2','wetso4')
           mycmor%positive = 'down'
        case ('loadoa','loadbc','loaddust','loadss','loadso4')
           var_info(var_found(1,1))%units = 'kg m-2'
        case ('wetss')
           mycmor%positive = 'up'
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
        case ('rlds','rldscs','rsds','rsdscs','rsdt','rtmt')
           mycmor%positive = 'down'
        case ('clt','ci','sci')
           var_info(var_found(1,1))%units = '1'
        case ('hurs','cl')
           var_info(var_found(1,1))%units = '%'
        case ('od550ss')
           var_info(var_found(1,1))%units = '1'
        case ('prc','precip','prsn')
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
           mycmor%positive = 'down'
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
        write(*,*) 'varid   = ',cmor_var_id
        write(*,*) 'comment = ',trim(xw(ixw)%comment)
        !
        ! Perform derivations and cycle through time, writing data too
        !
        select case (xw(ixw)%entry)
        case ('clt','od550aer','od550bc','od550oa','od550so4','ps',&
              'rlds','rlutcs','rsdscs','rsds','rsdt','tamax','tamin',&
              'wetso2','wetso4')
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
           enddo
           !
           ! Determine amount of data to write, to keep close to ~4 GB limit
           !
           select case(ntimes(1,1))
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
        case ('drydms','drydust','dryso2','dryso4','dryss','emidms','emiso2','emiso4')
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
           enddo
           !
           ! Determine amount of data to write, to keep close to ~4 GB limit
           !
           select case(ntimes(1,1))
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
        case ('precip')
           !
           ! precip : PRECT, unit change from m s-1 to kg m-2 s-1
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
        case ('emibc','rlus')
           !
           ! Add two fields
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
        case ('drypoa')
           !
           ! Add two fields, 1.4 scale factor too
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
                    cmordat2d = (indat2a + indat2b)*cscale
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
        case ('emioa')
           !
           ! Add two fields, 1.4 scale factor applied to two
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats))
!           allocate(indat2c(nlons,nlats))
!           allocate(indat2d(nlons,nlats),indat2e(nlons,nlats),indat2f(nlons,nlats))
!           allocate(indat2g(nlons,nlats))
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
!                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat2c)
!                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat2d)
!                 call read_var(myncid(1,5),var_info(var_found(1,5))%name,indat2e)
!                 call read_var(myncid(1,6),var_info(var_found(1,6))%name,indat2f)
!                 call read_var(myncid(1,7),var_info(var_found(1,7))%name,indat2g)
                 ! 
!                 cmordat2d = ((indat2a + indat2b)*cscale)+indat2c+indat2d+indat2e+indat2f+indat2g
                 cmordat2d = (indat2a + indat2b)*cscale
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
        case ('od550dust','od550ss','wetdust')
           !
           ! Sum four fields
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats),indat2d(nlons,nlats))
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
                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat2d)
                 !
                 where ((indat2a /= 1.e36).and.(indat2b /= 1.e36).and.(indat2c /= 1.e36).and.(indat2d /= 1.e36))
                    cmordat2d = indat2a + indat2b + indat2c + indat2d
                 elsewhere
                    cmordat2d = 1.e20
                 endwhere
                 !
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
        case ('emiss')
           !
           ! Sum four fields
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats),indat2d(nlons,nlats))
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
                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat2d)
                 !
                 where ((indat2a /= 1.e36).and.(indat2b /= 1.e36).and.(indat2c /= 1.e36).and.(indat2d /= 1.e36))
                    cmordat2d = indat2a + indat2b + indat2c + indat2d
                 elsewhere
                    cmordat2d = 1.e20
                 endwhere
!!$                 !
!!$                 ! Area weighting
!!$                 !
!!$                 write(*,*) 'CM2 0: ',minval(cmordat2d),maxval(cmordat2d)
!!$                 do j = 1,nlats
!!$                    where (cmordat2d(:,j) /= 1.e20)
!!$                       cmordat2d(:,j)=cmordat2d(:,j)*area_wt(j)
!!$                    endwhere
!!$                 enddo
!!$                 write(*,*) 'CM2 1: ',minval(cmordat2d),maxval(cmordat2d)
                 !
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
        case ('drysoa')
           !
           ! Sum five fields
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats),indat2d(nlons,nlats),indat2e(nlons,nlats))
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
                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat2d)
                 call read_var(myncid(1,5),var_info(var_found(1,5))%name,indat2e)
                 !
                 where ((indat2a /= 1.e36).and.(indat2b /= 1.e36).and.(indat2c /= 1.e36).and.(indat2d /= 1.e36).and.(indat2e /= 1.e36))
                    cmordat2d = indat2a + indat2b + indat2c + indat2d + indat2e
                 elsewhere
                    cmordat2d = 1.e20
                 endwhere
                 !
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
        case ('loaddust','loadss','wetss')
           !
           ! Sum four fields integrated over Z
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs))
           allocate(indat3c(nlons,nlats,nlevs),indat3d(nlons,nlats,nlevs))
           allocate(work3da(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats),pshybrid(nlons,nlats,nlevs),psdelta(nlons,nlats,nlevs))
           allocate(tsdata(nlons,nlats),rho(nlons,nlats))
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
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 work3da   = 1.e20
                 cmordat2d = 1.e20
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat3d)
                 call read_var(myncid(1,5),'PS',psdata)
                 call read_var(myncid(1,6),'TS',tsdata)
                 !
                 rho = psdata / (rdair * tsdata)
                 call pres_hybrid_ccm(psdata,pshybrid,nlons,nlats,nlevs)
                 do k = 1,nlevs-1
                    psdelta(:,:,k)=pshybrid(:,:,k+1)-pshybrid(:,:,k)
                 enddo
                 !
                 work3da   = (indat3a + indat3b + indat3c + indat3d)*(psdelta/grav)
                 cmordat2d = sum(work3da,dim=3)
!!$                 !
!!$                 ! Area weighting
!!$                 !
!!$                 write(*,*) 'CM2 0: ',minval(cmordat2d),maxval(cmordat2d)
!!$                 do j = 1,nlats
!!$                    cmordat2d(:,j)=cmordat2d(:,j)*area_wt(j)
!!$                 enddo
!!$                 write(*,*) 'CM2 1: ',minval(cmordat2d),maxval(cmordat2d)
                 !
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
        case ('loadbc')
           !
           ! Sum two fields integrated over Z
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs))
           allocate(work3da(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats),pshybrid(nlons,nlats,nlevs),psdelta(nlons,nlats,nlevs))
           allocate(tsdata(nlons,nlats),rho(nlons,nlats))
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
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 work3da   = 1.e20
                 cmordat2d = 1.e20
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                 call read_var(myncid(1,3),'PS',psdata)
                 call read_var(myncid(1,4),'TS',tsdata)
                 !
                 rho = psdata / (rdair * tsdata)
                 call pres_hybrid_ccm(psdata,pshybrid,nlons,nlats,nlevs)
                 do k = 1,nlevs-1
                    psdelta(:,:,k)=pshybrid(:,:,k+1)-pshybrid(:,:,k)
                 enddo
                 !
                 work3da   = (indat3a + indat3b)*(psdelta/grav)
                 write(*,*) 'CM3 0',minval(work3da),maxval(work3da)
                 cmordat2d = sum(work3da,dim=3)
!!$                 !
!!$                 ! Area weighting
!!$                 !
!!$                 write(*,*) 'CM2 0: ',minval(cmordat2d),maxval(cmordat2d)
!!$                 do j = 1,nlats
!!$                    cmordat2d(:,j)=cmordat2d(:,j)*area_wt(j)
!!$                 enddo
!!$                 write(*,*) 'CM2 1: ',minval(cmordat2d),maxval(cmordat2d)
                 !
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
        case ('loadoa')
           !
           ! Sum seven fields integrated over Z, 1.4X scale factor for OC1 and OC2
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs),indat3c(nlons,nlats,nlevs))
           allocate(indat3d(nlons,nlats,nlevs),indat3e(nlons,nlats,nlevs),indat3f(nlons,nlats,nlevs))
           allocate(indat3g(nlons,nlats,nlevs))
           allocate(work3da(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats),pshybrid(nlons,nlats,nlevs),psdelta(nlons,nlats,nlevs))
           allocate(tsdata(nlons,nlats),rho(nlons,nlats))
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
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 work3da   = 1.e20
                 cmordat2d = 1.e20
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat3d)
                 call read_var(myncid(1,5),var_info(var_found(1,5))%name,indat3e)
!                 call read_var(myncid(1,6),var_info(var_found(1,6))%name,indat3f)
!                 call read_var(myncid(1,7),var_info(var_found(1,7))%name,indat3g)
!                 call read_var(myncid(1,8),'PS',psdata)
!                 call read_var(myncid(1,9),'TS',tsdata)
                 call read_var(myncid(1,6),'PS',psdata)
                 call read_var(myncid(1,7),'TS',tsdata)
                 !
                 rho = psdata / (rdair * tsdata)
                 call pres_hybrid_ccm(psdata,pshybrid,nlons,nlats,nlevs)
                 do k = 1,nlevs-1
                    psdelta(:,:,k)=pshybrid(:,:,k+1)-pshybrid(:,:,k)
                 enddo
                 !
!                 work3da   = (((indat3a+indat3b)*cscale)+indat3c+indat3d+indat3e+indat3f+indat3g)*(psdelta/grav)
                 work3da   = (indat3a+indat3b+indat3c+indat3d+indat3e)*(psdelta/grav)
                 cmordat2d = sum(work3da,dim=3)
!!$                 !
!!$                 ! Area weighting
!!$                 !
!!$                 write(*,*) 'CM2 0: ',minval(cmordat2d),maxval(cmordat2d)
!!$                 do j = 1,nlats
!!$                    cmordat2d(:,j)=cmordat2d(:,j)*area_wt(j)
!!$                 enddo
!!$                 write(*,*) 'CM2 1: ',minval(cmordat2d),maxval(cmordat2d)
                 !
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
        case ('wetoa')
           !
           ! Old: sum six  fields integrated over Z, 1.4X scale factor for the first
           ! New: sum five fields integrated over Z
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats))
           allocate(indat2d(nlons,nlats),indat2e(nlons,nlats),indat2f(nlons,nlats))
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
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 cmordat2d    = 1.e20
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat2c)
                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat2d)
                 call read_var(myncid(1,5),var_info(var_found(1,5))%name,indat2e)
!                 call read_var(myncid(1,6),var_info(var_found(1,6))%name,indat2f)
!                 !
!                 cmordat2d = (indat2a*cscale)+indat2b+indat2c+indat2d+indat2e+indat2f
                 cmordat2d = indat2b+indat2c+indat2d+indat2e
!!$                 !
!!$                 ! Area weighting
!!$                 !
!!$                 write(*,*) 'CM2 0: ',minval(cmordat2d),maxval(cmordat2d)
!!$                 do j = 1,nlats
!!$                    cmordat2d(:,j)=cmordat2d(:,j)*area_wt(j)
!!$                 enddo
!!$                 write(*,*) 'CM2 1: ',minval(cmordat2d),maxval(cmordat2d)
                 !
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
        case ('loadso4')
           !
           ! One field integrated over Z
           !
           allocate(indat3a(nlons,nlats,nlevs))
           allocate(work3da(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats),pshybrid(nlons,nlats,nlevs),psdelta(nlons,nlats,nlevs))
           allocate(tsdata(nlons,nlats),rho(nlons,nlats))
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
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 work3da   = 1.e20
                 cmordat2d = 1.e20
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),'PS',psdata)
                 call read_var(myncid(1,3),'TS',tsdata)
                 !
                 rho = psdata / (rdair * tsdata)
                 call pres_hybrid_ccm(psdata,pshybrid,nlons,nlats,nlevs)
                 do k = 1,nlevs-1
                    psdelta(:,:,k)=pshybrid(:,:,k+1)-pshybrid(:,:,k)
                 enddo
                 !
                 work3da   = indat3a*(psdelta/grav)
                 cmordat2d = sum(indat3a,dim=3)
!!$                 !
!!$                 ! Area weighting
!!$                 !
!!$                 write(*,*) 'CM2 0: ',minval(cmordat2d),maxval(cmordat2d)
!!$                 do j = 1,nlats
!!$                    cmordat2d(:,j)=cmordat2d(:,j)*area_wt(j)
!!$                 enddo
!!$                 write(*,*) 'CM2 1: ',minval(cmordat2d),maxval(cmordat2d)
                 !
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
        end select
        if (allocated(psdata))    deallocate(psdata)
        if (allocated(indat2a))   deallocate(indat2a)
        if (allocated(indat2b))   deallocate(indat2b)
        if (allocated(indat2c))   deallocate(indat2c)
        if (allocated(indat2d))   deallocate(indat2d)
        if (allocated(indat2e))   deallocate(indat2e)
        if (allocated(indat2f))   deallocate(indat2f)
        if (allocated(indat2g))   deallocate(indat2g)
        if (allocated(cmordat2d)) deallocate(cmordat2d)
        if (allocated(indat3a))   deallocate(indat3a)
        if (allocated(indat3b))   deallocate(indat3b)
        if (allocated(indat3c))   deallocate(indat3c)
        if (allocated(indat3d))   deallocate(indat3d)
        if (allocated(indat3e))   deallocate(indat3e)
        if (allocated(indat3f))   deallocate(indat3f)
        if (allocated(indat3g))   deallocate(indat3g)
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
end program A2DM_CMOR
