program Amon_CMOR
  ! Convert CCSM4 atmospheric monthly (cam2.h0) data from single-field format
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
  use atm_grid_info
  use mycmor_info
  !
  implicit none
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  INTEGER::error_flag,var_ids
  REAL,DIMENSION(:,:),ALLOCATABLE::data2d
  REAL,DIMENSION(:,:,:),ALLOCATABLE::data3d
  real::allmax,allmin,scale_factor
  double precision,dimension(:)  ,allocatable::time
  double precision,dimension(:,:),allocatable::bnds_time
  DOUBLE PRECISION,DIMENSION(1)  ::tval
  DOUBLE PRECISION,DIMENSION(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::time_units,exp_file,xwalk_file,ncfile,table_file
  character(len=256)::case_read,comp_read,svar,tstr
  integer::i,j,m,n,tcount,it,ivar,length,iexp,jexp,var_found
  integer::ncid,ntimes,nsamps_read,itab,ixw
  integer::ilon,ilat,ipres,ilev,itim,itim2,ilon2,ilat2
  logical::continue
  !
  ! GO!
  !
  mycmor%table_file = 'Tables/CMIP5_Amon'
  call load_table
  !
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
  xwalk_file = 'xwalk_Amon.txt'
  call load_xwalk(xwalk_file)
  !
  ! Get experiment information
  !
  exp_file = 'experiments.txt'
  call load_exp(exp_file)
  !
  read(*,*) case_read
  read(*,*) comp_read
  read(*,*) nsamps_read
  !
  ! Get experiment metadata from exp table and input case information
  !
  call get_exp_metadata(case_read,comp_read,nsamps_read)
  !
  ! Get grid information
  !
  call get_atm_grid
  write(*,*) 'GRID: ',size(alats),size(alons)
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
  do itab = 1,num_tab
     do ixw = 1,num_xw
        mycmor%positive = ' '
        time_counter = 0
        var_counter  = 0
        error_flag   = 0
        var_found    = 0
        scale_factor = 1.
        allmax       = -1.e36
        allmin       =  1.e36
        continue     = .false.
        if (xw(ixw)%ncesm_vars == 1) then
           if (xw(ixw)%entry == table(itab)%variable_entry) then
              if (xw(ixw)%cesm_vars(1)(1:7) /= 'UNKNOWN') then
                 if (xw(ixw)%cesm_vars(1)(1:12) /= 'UNAVAILABLE') then
                    write(ncfile,'(''data/'',a,''.'',a,''.'',a,''.'',a,''01-'',a,''12.nc'')') &
                         trim(case_read),&
                         trim(comp_read),&
                         trim(xw(ixw)%cesm_vars(1)),&
                         exp(exp_found)%begin_end(1:4),&
                         exp(exp_found)%begin_end(6:9)
                    inquire(file=trim(ncfile),exist=continue)
                    if (.not.(continue)) then
                       write(ncfile,'(''data/'',a,''.'',a,''.'',a,''.'',a,''-01_cat_'',a,''-12.nc'')') &
                            trim(case_read),&
                            trim(comp_read),&
                            trim(xw(ixw)%cesm_vars(1)),&
                            exp(exp_found)%begin_end(1:4),&
                            exp(exp_found)%begin_end(6:9)
                       inquire(file=trim(ncfile),exist=continue)
                    endif
                    if (continue) then
                       write(*,'('' GOOD TO GO : '',a,'' == '',a,'' from CESM file: '',a)') &
                            trim(xw(ixw)%entry),&
                            trim(table(itab)%variable_entry),&
                         trim(ncfile)
                    else
                       write(*,*) trim(ncfile(1:)),' NOT FOUND.'
                    endif
                 else
                    write(*,'('' UNAVAILABLE: '',a,'' == '',a)') &
                         trim(xw(ixw)%entry),&
                         trim(table(itab)%variable_entry)
                 endif
              else
                 write(*,'('' UNKNOWN    : '',a,'' == '',a)') &
                      trim(xw(ixw)%entry),&
                      trim(table(itab)%variable_entry)
              endif
           endif
           !
           ! Open CESM file and get information
           !
           if (continue) then
              call open_cdf(ncid,trim(ncfile),.true.)
              write(*,'(''OPENING: '',a80,'' ncid: '',i10)') trim(ncfile),ncid
              call get_dims(ncid)
              call get_vars(ncid)
              !
              do n=1,dim_counter
                 length = len_trim(dim_info(n)%name)
                 if(dim_info(n)%name(:length).eq.'time') then
                    ntimes = dim_info(n)%length
                 endif
              enddo
              time_units(1:) = ' '
              call read_att_text(ncid,'time','units',time_units)
              !
              do n=1,var_counter
                 if (trim(var_info(n)%name) == trim(xw(ixw)%cesm_vars(1))) then
                    var_found = n
                 endif
              enddo
              if (var_found == 0) then
                 write(*,*) trim(xw(ixw)%cesm_vars(1)),' NEVER FOUND. STOP.'
                 stop
              endif
              !
              ALLOCATE(time(ntimes),bnds_time(2,ntimes))
              !
              do n=1,ntimes
                 time_counter = n
                 call read_var(ncid,'time_bnds',bnds_time(:,n))
                 time(n) = (bnds_time(1,n)+bnds_time(2,n))/2.
              enddo
              !
              ! Specify path where tables can be found and indicate that existing netCDF files should be overwritten.
              !
              error_flag = cmor_setup(inpath='CMOR',netcdf_file_action=CMOR_REPLACE,logfile='LOG_CMOR.'//trim(xw(ixw)%entry))
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
              !
              ! Add acknowledgements
              !
              if (exp(exp_found)%loc(1:2) == 'NC') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(mycmor%ack_NC))
              if (exp(exp_found)%loc(1:2) == 'NE') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(mycmor%ack_NE))
              if (exp(exp_found)%loc(1:2) == 'OR') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(mycmor%ack_OR))
              !
              ! Add grid information
              !
              if (exp(exp_found)%grid(1:1) /= ' ') error_flag = cmor_set_cur_dataset_attribute("resolution",trim(adjustl(exp(exp_found)%grid)))
              !
              ! Add additional forcing information
              !
              if (mycmor%forcing_note(1:1) /= ' ') error_flag = cmor_set_cur_dataset_attribute("forcing_note",trim(adjustl(mycmor%forcing_note)))
              !
              ! Define all axes that will be needed
              !
              ilat = cmor_axis(                  &
                   table=mycmor%table_file,  &
                   table_entry='latitude',       &
                   units='degrees_north',        &
                   length=SIZE(alats),           &
                   coord_vals=alats,             &
                   cell_bounds=bnds_lat)

              ilon = cmor_axis(  &
                   table=mycmor%table_file,  &
                   table_entry='longitude',      &
                   length=SIZE(alons),           &
                   units='degrees_east',         &
                   coord_vals=alons,             &
                   cell_bounds=bnds_lon)
              !
              ! Note that the time axis is defined next, but the time coordinate
              ! values and bounds will be passed to cmor through function
              ! cmor_write (later, below).
              !
              itim = cmor_axis(  &
                   table=mycmor%table_file,      &
                   table_entry='time',           &
                   units=time_units,             &
                   length=ntimes,                &
                   interval='30 days')

              write(*,*) 'CMOR axes defined'
              ! 
              ! Make manual alterations so that CMOR works. Silly code!
              !
              if (xw(ixw)%entry == 'hurs')   var_info(var_found)%units = '1'
              if (xw(ixw)%entry == 'clt')    var_info(var_found)%units = '1'
              if (xw(ixw)%entry == 'ci')     var_info(var_found)%units = '1'
              if (xw(ixw)%entry == 'tauu')   mycmor%positive = 'up'
              if (xw(ixw)%entry == 'tauv')   mycmor%positive = 'up'
              if (xw(ixw)%entry == 'hfss')   mycmor%positive = 'up'
              if (xw(ixw)%entry == 'rlds')   mycmor%positive = 'down'
              if (xw(ixw)%entry == 'rlut')   mycmor%positive = 'up'
              if (xw(ixw)%entry == 'rlutcs') mycmor%positive = 'up'
              if (xw(ixw)%entry == 'rldscs') mycmor%positive = 'down'
              if (xw(ixw)%entry == 'rsds')   mycmor%positive = 'down'
              if (xw(ixw)%entry == 'rsdscs') mycmor%positive = 'down'
              if (xw(ixw)%entry == 'rsdt')   mycmor%positive = 'down'
              if (xw(ixw)%entry == 'hfls')   mycmor%positive = 'up'
              if (xw(ixw)%entry == 'prc') then
                 scale_factor = 1000.
                 var_info(var_found)%units = 'kg m-2 s-1'
              endif
              !
              write(*,*) &
                   'table=',trim(mycmor%table_file),' ',                  &
                   'table_entry=',trim(xw(ixw)%entry),' ',                &
                   'units=',trim(var_info(var_found)%units),' ',          &
                   'missing_value=',var_info(var_found)%missing_value,' ',&
                   'positive=',trim(mycmor%positive),' ',                 &
                   'original_name=',xw(ixw)%cesm_vars(1)
              !
              var_ids = cmor_variable(                     &
                   table=mycmor%table_file,                  &
                   table_entry=xw(ixw)%entry,               &
                   units=var_info(var_found)%units,              &
                   axis_ids=(/ ilon, ilat, itim /),              &
                   missing_value=var_info(var_found)%missing_value,&
                   positive=mycmor%positive,                       &
                   original_name=xw(ixw)%cesm_vars(1))
              !
              allocate(data2d(nlons,nlats))
              ! Cycle through time
              time_loop: DO it=1, ntimes
                 ! append each to the appropriate netCDF file.
                 time_counter = it
!                 write(*,'(''Reading '',a20,'' T= '',i10)') trim(var_info(var_found)%name),it
                 call read_var(ncid,var_info(var_found)%name,data2d)
                 tval(1)   = time(it)
                 tbnd(1,1) = bnds_time(1,it)
                 tbnd(2,1) = bnds_time(2,it)
                 m = 1
                 error_flag = cmor_write(      &
                      var_id        = var_ids, &
                      data          = data2d*scale_factor,  &
                      ntimes_passed = 1,       &
                      time_vals     = tval,    &
                      time_bnds     = tbnd)
                 !     write(*,'(i6,2g12.4)') it,minval(data2d),maxval(data2d)
                 allmax = max(allmax,maxval(data2d))
                 allmin = min(allmin,minval(data2d))

                 if (error_flag < 0) then
                    write(*,*) 'Error writing ',xw(ixw)%entry, ', which I call ', xw(ixw)%cesm_vars
                    write(*,*) 'Processing time sample: ', time
                    stop
                 end if

              end do time_loop
              ivar = ivar + 1
              !
              ! Close all files opened by CMOR.
              !
              error_flag = cmor_close()
              write(*,'(''********************************************************************************'')')
              write(*,'(''********************************************************************************'')')
              write(*,'(''CMOR executed to completion; T#: '',i5,'' X#: '',i5,'' EXT: '',2g10.4)') itab,ixw,allmin,allmax
              write(*,'(''********************************************************************************'')')
              write(*,'(''********************************************************************************'')')
              call close_cdf(ncid)
              time_counter = 0
              var_counter  = 0
              error_flag   = 0
              var_found    = 0
              scale_factor = 1.
              allmax       = -1.e36
              allmin       =  1.e36
              continue     = .false.
              mycmor%positive = ' '
           endif
        endif
        continue = .false.
     end do
  end do
end program Amon_CMOR
