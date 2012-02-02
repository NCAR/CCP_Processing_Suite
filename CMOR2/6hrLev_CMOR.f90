program Do6hrLev_CMOR
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
  real,dimension(:,:,:),allocatable::indat3a,indat3b,indat3c,cmordat3d,work3da,work3db
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::time_bnds
  double precision,dimension(1)  ::tval
  double precision,dimension(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile,cmor_filename
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,itab,ixw,ilev,ic
  real::spval
  !
  ! GO!
  !
  mycmor%table_file = 'CMIP5_6hrLev'
  call load_table_info
  !
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
  xwalk_file = 'xwalk_6hrLev.txt'
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
  table_loop: do itab = 1,num_tab
     xwalk_loop: do ixw = 1,num_xw
        mycmor%positive = ' '
        time_counter    = 0
        var_counter     = 0
        error_flag      = 0
        var_found       = 0
        all_continue    = .false.
        time_units      = ' '
        original_name   = ' '
        ncfile(:,:)(1:) = ' '
        nc_nfiles(:)    = 0
        !
        ! The meaty part
        !
        if (xw(ixw)%entry == table(itab)%variable_entry) then
           write(*,'(''MATCH CMIP5: '',a,'' CESM: '',5(a,'',''))') trim(xw(ixw)%entry),(trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
           do ivar = 1,xw(ixw)%ncesm_vars
              if ((trim(xw(ixw)%cesm_vars(ivar)) == 'UNKNOWN').or.(trim(xw(ixw)%cesm_vars(ivar)) == 'UNAVAILABLE')) then
                 write(*,'(''UNAVAILABLE/UNKNOWN: '',a,'' == '',a)') trim(xw(ixw)%entry),trim(table(itab)%variable_entry)
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
                 call read_att_text(myncid(1,1),'time','units',time_units)
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
                 !
                 if (.not.(allocated(time)))      then
                    allocate(time(ntimes(1,1)))
                 endif
                 !
                 do n=1,ntimes(1,1)
                    time_counter = n
                    call read_var(myncid(1,ivar),'time',time(n))
                 enddo
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
              call define_atm_axes(table(itab)%dimensions)
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
              spval=var_info(var_found(1,1))%missing_value
              !
              write(*,*) 'calling cmor_variable:'
              write(*,*) 'table         = ',trim(mycmor%table_file)
              write(*,*) 'table_entry   = ',trim(xw(ixw)%entry)
              write(*,*) 'dimensions    = ',trim(table(itab)%dimensions)
              write(*,*) 'units         = ',var_info(var_found(1,1))%units(1:20)
              write(*,*) 'axis_ids      = ',axis_ids(1:4)
              write(*,*) 'missing_value = ',var_info(var_found(1,1))%missing_value
              write(*,*) 'positive      = ',trim(mycmor%positive)
              write(*,*) 'original_name = ',trim(original_name)
              !
              select case (xw(ixw)%entry)
              case ('ta','ua','va','hus')
                 cmor_var_id = cmor_variable(                            &
                      table=mycmor%table_file,                           &
                      table_entry=xw(ixw)%entry,                         &
                      units=var_info(var_found(1,1))%units,                &
                      axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3),axis_ids(4)/),  &
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
                 allocate(indat2a(nlons,nlats))
                 !
                 ! Determine amount of data to write, to keep close to ~2 GB limit
                 !
                 select case(ntimes(1,1))
                 case ( 1460 )  ! One year, four pieces, one per calendar quarter 01/01-03/31,04/01-06/30,07/01-09/30,10/01-12/31
                    nchunks(1) = 4
                    tidx1(1:nchunks(1)) = (/  1, 360, 725, 1093/)
                    tidx2(1:nchunks(1)) = (/359, 724,1092, ntimes(1,1)/)
                 case default
                    nchunks(1) = 1
                    tidx1(1:nchunks(1)) = 1
                    tidx2(1:nchunks(1)) = ntimes(1,1)
                 end select
                 write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
                 do ic = 1,nchunks(1)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
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
              case ('ta','ua','va','hus')
                 !
                 ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
                 ! break up into nicely-sized chunks along time
                 !
                 allocate(indat3a(nlons,nlats,nlevs),indat2a(nlons,nlats))
                 !
                 ! Determine amount of data to write, to keep close to ~2 GB limit
                 !
                 select case(ntimes(1,1))
                 case ( 1460 )  ! One year, four pieces, one per calendar quarter 01/01-03/31,04/01-06/30,07/01-09/30,10/01-12/31
                    nchunks(1) = 4
                    tidx1(1:nchunks(1)) = (/  1, 360, 725, 1093/)
                    tidx2(1:nchunks(1)) = (/359, 724,1092, ntimes(1,1)/)
                 case default
                    nchunks(1) = 1
                    tidx1(nchunks(1)) = 1
                    tidx2(nchunks(1)) = ntimes(1,1)
                 end select
                 write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),'',''))') nchunks(1),(tidx1(ic),tidx2(ic),ic=1,nchunks(1))
                 do ic = 1,nchunks(1)
                    do it = tidx1(ic),tidx2(ic)
                       time_counter = it
                       call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                       call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat2a)
                       where (indat3a > 1.e6) indat3a = spval
                       tval(1) = time(it)
                       error_flag = cmor_write(        &
                            var_id        = cmor_var_id,   &
                            data          = indat3a,   &
                            ntimes_passed = 1,         &
                            time_vals     = tval)
                       if (error_flag < 0) then
                          write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                          stop
                       endif
                       error_flag = cmor_write(        &
                            var_id        = zfactor_id,&
                            data          = indat2a,   &
                            ntimes_passed = 1,         &
                            time_vals     = tval,      &
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
              if (allocated(indat2b))   deallocate(indat2b)
              if (allocated(indat2c))   deallocate(indat2c)
              if (allocated(cmordat2d)) deallocate(cmordat2d)
              if (allocated(indat3a))   deallocate(indat3a)
              if (allocated(indat3b))   deallocate(indat3b)
              if (allocated(work3da))   deallocate(work3da)
              if (allocated(work3db))   deallocate(work3db)
              do ivar = 1,xw(ixw)%ncesm_vars
                 call close_cdf(myncid(1,ivar))
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
              !
              error_flag = cmor_close()
              if (error_flag < 0) then
                 write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
              else
                 write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
              endif
           endif
        endif
     enddo xwalk_loop
  enddo table_loop
end program Do6hrLev_CMOR
