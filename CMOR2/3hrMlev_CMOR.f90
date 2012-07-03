program Do3hrMlev_CMOR
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
  real,dimension(:,:)  ,allocatable::indat2a,indat2b,indat2c,cmordat2d,psdata
  real,dimension(:,:,:),allocatable::indat3a,indat3b,indat3c,indat3d,cmordat3d
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::time_bnds
  double precision,dimension(1)  ::tval
  double precision,dimension(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile,cmor_filename
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,ixw,ilev,ic
  real::spval,cpair
  logical::does_exist
  !
  cpair = 1004.64
  !
  ! GO!
  !
  mycmor%table_file = '3hrMlev'
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
        do ivar = 1,xw(ixw)%ncesm_vars
           write(*,'(''TO OPEN: '',a)') trim(ncfile(nc_nfiles(ivar),ivar))
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
        case default
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(2),axis_ids(3),axis_ids(4),axis_ids(1)/),  &
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
        case ('ps')
           !
           ! No change
           !
           allocate(psdata(nlons,nlats))
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
           ! Determine amount of data to write, to keep close to ~2 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872,1140,3612,6012,12012 )  ! All data
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
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
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,psdata)
                 tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                 error_flag = cmor_write(          &
                      var_id        = cmor_var_id, &
                      data          = psdata,     &
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
        case ('mc')
           !
           ! Sum two fields
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nilevs),indat3b(nlons,nlats,nilevs),psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nilevs))
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
           ! Determine amount of data to write, to keep close to ~2 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872 )  ! 20C, 1850-2005, ~50y chunks
              nchunks(1) = 3
              tidx1(1:nchunks(1)) = (/  1, 601,1201/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/600,1200,1872/) ! 1899, 1950, 2005
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/ 13, 541/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/540,1152/)      ! 2049, 2100
           case ( 1140 )  ! RCP, 2006-2100
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 529/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/528,1140/)      ! 2049, 2100
           case ( 900 )
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 601/)      ! 1850, 1900
              tidx2(1:nchunks(1)) = (/600, 900/)      ! 1899, 1924
           case ( 3612,6012,12012 ) ! piControl,past1000,midHolocene: ~50Y chunks
              nchunks(1) = int(ntimes(1,1)/600)
              tidx1(1) =   1
              tidx2(1) = 600
              do ic = 2,nchunks(1)
                 tidx1(ic) = tidx2(ic-1) + 1
                 tidx2(ic) = tidx1(ic) + 599
              enddo
              tidx2(nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only, ~50y chunks
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/3613,4213/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/4212,4824/) ! 1899, 1950, 2005
           case ( 2400 )  ! WACCM from 96-295 (200 y) chunks from 96-149, 150-199, 200-249, 250-295
              nchunks(1) = 4
              tidx1(1:nchunks(1)) = (/   1, 649,1249,1849/) !   96,  150, 200, 250
              tidx2(1:nchunks(1)) = (/ 648,1248,1848,2400/) !  149,  199, 249, 295
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,psdata)
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
        case ('tntr','tnthus')
           !
           ! Sum two fields
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs),psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nlevs))
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
           ! Determine amount of data to write, to keep close to ~2 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872 )  ! 20C, 1850-2005, ~50y chunks
              nchunks(1) = 3
              tidx1(1:nchunks(1)) = (/  1, 601,1201/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/600,1200,1872/) ! 1899, 1950, 2005
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/ 13, 541/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/540,1152/)      ! 2049, 2100
           case ( 1140 )  ! RCP, 2006-2100
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 529/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/528,1140/)      ! 2049, 2100
           case ( 900 )
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 601/)      ! 1850, 1900
              tidx2(1:nchunks(1)) = (/600, 900/)      ! 1899, 1924
           case ( 3612,6012,12012 ) ! piControl,past1000,midHolocene: ~50Y chunks
              nchunks(1) = int(ntimes(1,1)/600)
              tidx1(1) =   1
              tidx2(1) = 600
              do ic = 2,nchunks(1)
                 tidx1(ic) = tidx2(ic-1) + 1
                 tidx2(ic) = tidx1(ic) + 599
              enddo
              tidx2(nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only, ~50y chunks
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/3613,4213/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/4212,4824/) ! 1899, 1950, 2005
           case ( 2400 )  ! WACCM from 96-295 (200 y) chunks from 96-149, 150-199, 200-249, 250-295
              nchunks(1) = 4
              tidx1(1:nchunks(1)) = (/   1, 649,1249,1849/) !   96,  150, 200, 250
              tidx2(1:nchunks(1)) = (/ 648,1248,1848,2400/) !  149,  199, 249, 295
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,psdata)
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
        case ('tntc')
           !
           ! ZMDT+EVAPTZM+ZMMTT+CMFDT
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs))
           allocate(indat3c(nlons,nlats,nlevs),indat3d(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nlevs))
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
           ! Determine amount of data to write, to keep close to ~2 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872 )  ! 20C, 1850-2005, ~50y chunks
              nchunks(1) = 3
              tidx1(1:nchunks(1)) = (/  1, 601,1201/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/600,1200,1872/) ! 1899, 1950, 2005
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/ 13, 541/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/540,1152/)      ! 2049, 2100
           case ( 1140 )  ! RCP, 2006-2100
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 529/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/528,1140/)      ! 2049, 2100
           case ( 900 )
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 601/)      ! 1850, 1900
              tidx2(1:nchunks(1)) = (/600, 900/)      ! 1899, 1924
           case ( 3612,6012,12012 ) ! piControl,past1000,midHolocene: ~50Y chunks
              nchunks(1) = int(ntimes(1,1)/600)
              tidx1(1) =   1
              tidx2(1) = 600
              do ic = 2,nchunks(1)
                 tidx1(ic) = tidx2(ic-1) + 1
                 tidx2(ic) = tidx1(ic) + 599
              enddo
              tidx2(nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only, ~50y chunks
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/3613,4213/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/4212,4824/) ! 1899, 1950, 2005
           case ( 2400 )  ! WACCM from 96-295 (200 y) chunks from 96-149, 150-199, 200-249, 250-295
              nchunks(1) = 4
              tidx1(1:nchunks(1)) = (/   1, 649,1249,1849/) !   96,  150, 200, 250
              tidx2(1:nchunks(1)) = (/ 648,1248,1848,2400/) !  149,  199, 249, 295
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat3d)
                 call read_var(myncid(1,5),var_info(var_found(1,5))%name,psdata)
                 where ((indat3a /= spval).and.(indat3b /= spval).and.(indat3c /= spval).and.(indat3d /= spval))
                    cmordat3d = indat3a + indat3b + indat3c + indat3d
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
        case ('tnhusc')
           !
           ! ZMDQ+EVAPQZM+CMFDQ
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs),indat3c(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nlevs))
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
           ! Determine amount of data to write, to keep close to ~2 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872 )  ! 20C, 1850-2005, ~50y chunks
              nchunks(1) = 3
              tidx1(1:nchunks(1)) = (/  1, 601,1201/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/600,1200,1872/) ! 1899, 1950, 2005
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/ 13, 541/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/540,1152/)      ! 2049, 2100
           case ( 1140 )  ! RCP, 2006-2100
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 529/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/528,1140/)      ! 2049, 2100
           case ( 900 )
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 601/)      ! 1850, 1900
              tidx2(1:nchunks(1)) = (/600, 900/)      ! 1899, 1924
           case ( 3612,6012,12012 ) ! piControl,past1000,midHolocene: ~50Y chunks
              nchunks(1) = int(ntimes(1,1)/600)
              tidx1(1) =   1
              tidx2(1) = 600
              do ic = 2,nchunks(1)
                 tidx1(ic) = tidx2(ic-1) + 1
                 tidx2(ic) = tidx1(ic) + 599
              enddo
              tidx2(nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only, ~50y chunks
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/3613,4213/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/4212,4824/) ! 1899, 1950, 2005
           case ( 2400 )  ! WACCM from 96-295 (200 y) chunks from 96-149, 150-199, 200-249, 250-295
              nchunks(1) = 4
              tidx1(1:nchunks(1)) = (/   1, 649,1249,1849/) !   96,  150, 200, 250
              tidx2(1:nchunks(1)) = (/ 648,1248,1848,2400/) !  149,  199, 249, 295
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,psdata)
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
        case ('tntscpbl')
           !
           ! HPROGCLD/CPAIR+HSED/CPAIR (CPAIR = 1004.64)
           !
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs),psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nlevs))
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
           ! Determine amount of data to write, to keep close to ~2 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872 )  ! 20C, 1850-2005, ~50y chunks
              nchunks(1) = 3
              tidx1(1:nchunks(1)) = (/  1, 601,1201/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/600,1200,1872/) ! 1899, 1950, 2005
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/ 13, 541/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/540,1152/)      ! 2049, 2100
           case ( 1140 )  ! RCP, 2006-2100
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 529/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/528,1140/)      ! 2049, 2100
           case ( 900 )
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 601/)      ! 1850, 1900
              tidx2(1:nchunks(1)) = (/600, 900/)      ! 1899, 1924
           case ( 3612,6012,12012 ) ! piControl,past1000,midHolocene: ~50Y chunks
              nchunks(1) = int(ntimes(1,1)/600)
              tidx1(1) =   1
              tidx2(1) = 600
              do ic = 2,nchunks(1)
                 tidx1(ic) = tidx2(ic-1) + 1
                 tidx2(ic) = tidx1(ic) + 599
              enddo
              tidx2(nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only, ~50y chunks
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/3613,4213/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/4212,4824/) ! 1899, 1950, 2005
           case ( 2400 )  ! WACCM from 96-295 (200 y) chunks from 96-149, 150-199, 200-249, 250-295
              nchunks(1) = 4
              tidx1(1:nchunks(1)) = (/   1, 649,1249,1849/) !   96,  150, 200, 250
              tidx2(1:nchunks(1)) = (/ 648,1248,1848,2400/) !  149,  199, 249, 295
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
           write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
           do ic = 1,nchunks(1)
              do it = tidx1(ic),tidx2(ic)
                 time_counter = it
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,psdata)
                 where ((indat3a /= spval).and.(indat3b /= spval))
                    cmordat3d = (indat3a/cpair) + (indat3b/cpair)
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
           ! Determine amount of data to write, to keep close to ~2 GB limit
           !
           select case(ntimes(1,1))
           case ( 1872 )  ! 20C, 1850-2005, ~50y chunks
              nchunks(1) = 3
              tidx1(1:nchunks(1)) = (/  1, 601,1201/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/600,1200,1872/) ! 1899, 1950, 2005
           case ( 1152 )  ! RCP, 2005-2100, skip 2006
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/ 13, 541/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/540,1152/)      ! 2049, 2100
           case ( 1140 )  ! RCP, 2006-2100
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 529/)      ! 2006, 2050
              tidx2(1:nchunks(1)) = (/528,1140/)      ! 2049, 2100
           case ( 900 )
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/  1, 601/)      ! 1850, 1900
              tidx2(1:nchunks(1)) = (/600, 900/)      ! 1899, 1924
           case ( 3612,6012,12012 ) ! piControl,past1000,midHolocene: ~50Y chunks
              nchunks(1) = int(ntimes(1,1)/600)
              tidx1(1) =   1
              tidx2(1) = 600
              do ic = 2,nchunks(1)
                 tidx1(ic) = tidx2(ic-1) + 1
                 tidx2(ic) = tidx1(ic) + 599
              enddo
              tidx2(nchunks(1)) = ntimes(1,1)
           case ( 4824 )  ! LGM from 1499-1900, 1800-1900 (101y) only, ~50y chunks
              nchunks(1) = 2
              tidx1(1:nchunks(1)) = (/3613,4213/) ! 1850, 1900, 1951
              tidx2(1:nchunks(1)) = (/4212,4824/) ! 1899, 1950, 2005
           case ( 2400 )  ! WACCM from 96-295 (200 y) chunks from 96-149, 150-199, 200-249, 250-295
              nchunks(1) = 4
              tidx1(1:nchunks(1)) = (/   1, 649,1249,1849/) !   96,  150, 200, 250
              tidx2(1:nchunks(1)) = (/ 648,1248,1848,2400/) !  149,  199, 249, 295
           case default
              nchunks(1) = 1
              tidx1(1:nchunks(1)) = 1
              tidx2(1:nchunks(1)) = ntimes(1,1)
           end select
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
        if (allocated(indat2a))   deallocate(indat2a)
        if (allocated(psdata))    deallocate(psdata)
        if (allocated(indat2b))   deallocate(indat2b)
        if (allocated(indat2c))   deallocate(indat2c)
        if (allocated(cmordat2d)) deallocate(cmordat2d)
        if (allocated(indat3a))   deallocate(indat3a)
        if (allocated(indat3b))   deallocate(indat3b)
        if (allocated(indat3c))   deallocate(indat3c)
        if (allocated(indat3d))   deallocate(indat3d)
        do ivar = 1,xw(ixw)%ncesm_vars
           call close_cdf(myncid(1,ivar))
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
end program Do3hrMlev_CMOR
