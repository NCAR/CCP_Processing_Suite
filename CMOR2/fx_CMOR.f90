program fx_CMOR
  ! Convert CCSM4 data from single-field format to CMOR-compliant format
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
!  integer,parameter::nexp = 56
  integer,parameter::nexp = 1
  real,parameter::shr_const_g  = 9.80616 ! acceleration of gravity ~ m/s^2
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  integer::error_flag,cmor_var_id
  real,dimension(:)    ,allocatable::dx,dy,dlat
  integer,dimension(:,:)  ,allocatable::region_mask,basin
  real,dimension(:,:)  ,allocatable::indat2a,cmordat2d
  real,dimension(:,:,:),allocatable::indat3a,work3da,cmordat3d
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile
  character(len=512)::acknowledgement
  integer::i,j,k,m,n,tcount,ivar,length,iexp,jexp,ixw
  character(len=15),dimension(nexp)::expids
  character(len=256)::whoami,prochost,ccps_rev,ccps_date,ccps_uuid,info_file
  character(len=10)::pdate,ptime
  logical::exists
  !
!!$  expids = (/&
!!$       ! Long-term fully-coupled
!!$       'piControl     ','1pctCO2       ','abrupt4xCO2   ',&
!!$       'historical    ','historicalExt ','historicalGHG ','historicalMisc','historicalNat ',&
!!$       'rcp26         ','rcp45         ','rcp60         ','rcp85         ',&
!!$       ! Paleo runs
!!$       'lgm           ','midHolocene   ','past1000      ',&
!!$       ! AMIP, aquaplanet and other atm-only runs
!!$       'amip          ','amip4K        ','amip4xCO2     ','amipFuture    ',&
!!$       'aqua4K        ','aqua4xCO2     ','aquaControl   ',&
!!$       'sst2030       ','sstClim       ','sstClim4xCO2  ','sstClimAerosol','sstClimSulfate',&
!!$       ! ESM runs
!!$       'esmControl    ','esmFdbk1      ','esmFdbk2      ','esmFixClim1   ','esmFixClim2   ','esmHistorical ','esmrcp85      ',&
!!$       ! Decadal predictions
!!$       'decadal1961   ','decadal1966   ',&
!!$       'decadal1971   ','decadal1975   ','decadal1976   ',&
!!$       'decadal1980   ','decadal1981   ','decadal1985   ','decadal1986   ',&
!!$       'decadal1990   ','decadal1991   ','decadal1995   ','decadal1996   ',&
!!$       'decadal2000   ','decadal2001   ','decadal2002   ','decadal2003   ','decadal2004   ','decadal2005   ','decadal2006   ',&
!!$       'noVolcXXXX    ','volcIn2010    '/)
  !
  ! GO!
  !
  mycmor%table_file = 'PMIP3_fx'
  expids = ('PlioExp2a      ')
  !
  ! Get table information
  !
  mycmor%table_file = 'Tables/'//trim(mycmor%table_file)
  inquire(file=mycmor%table_file,exist=exists)
  if (.not.(exists)) then
     write(*,*) 'Cannot find ',trim(mycmor%table_file),'. Dying.'
     stop
  endif
  !
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
  xwalk_file = 'xwalk_PMIP3_fx'
  call load_xwalk(xwalk_file)
  !
  ! Get grid information
  !
  mycmor%model_id      = 'CCSM4'
  exp_found = 1
  exp(exp_found)%model_id = 'CCSM4'
  call get_atm_grid
!  call get_lnd_grid
!  call get_ocn_grid
  !
  ! Set up CMOR subroutine arguments
  !
  acknowledgement = 'The CESM project is supported by the National Science Foundation and the Office of Science (BER) of the U.S. Department of Energy. '//&
       'NCAR is sponsored by the National Science Foundation. '//&
       'Computing resources were provided by the Climate Simulation Laboratory at the NCAR Computational and Information Systems Laboratory (CISL), '//&
       'sponsored by the National Science Foundation and other agencies.'//&
       'This research used resources of the National Energy Research Scientific Computing Center, which is supported by the Office of Science (BER) '//&
       'of the U.S. Department of Energy under Contract No. DE-AC02-05CH11231'//&
       'This research used resources of the Oak Ridge Leadership Computing Facility, located in the National Center for Computational Sciences '//&
       'at Oak Ridge National Laboratory, which is supported by the Office of Science (BER) of the Department of Energy under Contract DE-AC05-00OR22725.'
  !
  mycmor%forcing_note = 'Additional information on the external forcings used in this experiment can be found at '//&
       'http://www.cesm.ucar.edu/CMIP5/forcing_information'
  !
  ! Define arguments to 'cmor_dataset' - set by load_exp and init routines 
  !
  mycmor%outpath       = 'CMOR'
  mycmor%source        = 'CCSM4'
  mycmor%calendar      = 'noleap'
  mycmor%contact       = 'cesm_data@ucar.edu'
  mycmor%history       = ' '
  mycmor%comment       = ' '
  mycmor%positive      = ' '
  !
  ! References
  !
  select case (mycmor%model_id)
  case ('CCSM4')
     mycmor%references    = 'Gent P. R., et.al. 2011: The Community Climate System Model version 4. J. Climate, doi: 10.1175/2011JCLI4083.1'
     mycmor%institute_id  = 'NCAR'
     mycmor%institution   = 'NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-CAM5')
     mycmor%references    = 'TBD'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-BGC')
     mycmor%references    = 'TBD'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-FASTCHEM')
     mycmor%references    = 'TBD'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case ('CESM1-WACCM')
     mycmor%references    = 'TBD'
     mycmor%institute_id  = 'NSF-DOE-NCAR'
     mycmor%institution   = 'NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA'
  case default
     write(*,*) 'Unknown model_id: ',trim(adjustl(exp(exp_found)%model_id)),' Stopping.'
     stop
  end select
  mycmor%forcing         = 'N/A'
  !
  ! Step through CMOR table entries to see what CESM fields we can read and in process, and if so, do it!
  !
  expt_loop: do iexp = 1,nexp
     xwalk_loop: do ixw = 1,num_xw
        call reset_netcdf_var
        mycmor%positive = ' '
        error_flag   = 0
        original_name= ' '
        !
        ! The meaty part
        !
        xw(ixw)%entry = 'orog'
        select case (xw(ixw)%entry)
        case ( 'areacella','sftlf','orog' )
           ncfile(1,1) = 'atm_grid_cam4_f09.nc'
           xw(ixw)%ncesm_vars = 1
        case ( 'mrsofc','sftgif','rootd' )
           ncfile(1,1) = 'lnd_grid_f09.nc'
           xw(ixw)%ncesm_vars = 1
        case ( 'areacello','basin','deptho','hfgeou','sftof','thkcello','volcello')
           ncfile(1,1) = 'ocn_grid_gx1.nc'
           xw(ixw)%ncesm_vars = 1
        end select
        call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
        write(*,'(''OPENING: '',a,'' myncid: '',i10)') trim(ncfile(1,1)),myncid(1,1)
        call get_dims(myncid(1,1))
        call get_vars(myncid(1,1))
        !
        var_loop: do ivar = 1,xw(ixw)%ncesm_vars
           if (trim(xw(ixw)%cesm_vars(ivar)) == 'UNKNOWN') then
              write(*,'(a,'' IS UNKNOWN EQUIVALENCE.'')') trim(xw(ixw)%cesm_vars(ivar))
           else
              do n=1,var_counter
                 if (trim(var_info(n)%name) == trim(xw(ixw)%cesm_vars(ivar))) then
                    write(*,*) 'var_found n: ',n,' name: ',trim(var_info(n)%name),' ixw ',ixw,' entry: ',trim(xw(ixw)%cesm_vars(ivar))
                    var_found(1,ivar) = n
                    xw_found = ixw
                 endif
              enddo
              if (var_found(1,ivar) == 0) then
                 write(*,'(''NEVER FOUND: '',a,'' STOP. '')') trim(xw(ixw)%cesm_vars(ivar))
                 stop
              endif
              !
              ! Specify path where tables can be found and indicate that existing netCDF files should be overwritten.
              !
              write(logfile,'(''log_cmor.'',a,''.'',a,''_'',a)') &
                   trim(expids(iexp)),&
                   'r0i0p0',&
                   trim(xw(ixw)%entry)
              error_flag = cmor_setup(inpath='CMOR',&
                   netcdf_file_action=CMOR_REPLACE,&
                   logfile=logfile)
              !
              error_flag = cmor_dataset(                              &
                   outpath=mycmor%outpath,                            &
                   experiment_id=trim(expids(iexp)),                  &
                   institution=mycmor%institution,                    &
                   source=mycmor%source,                              &
                   calendar=mycmor%calendar,                          &
                   realization=0,                                     &
                   contact=mycmor%contact,                            &
                   history=mycmor%history,                            &
                   comment=mycmor%comment,                            &
                   references=mycmor%references,                      &
                   model_id=mycmor%model_id,                          &
                   forcing=mycmor%forcing,                            &
                   initialization_method=0,                           &
                   physics_version=0,                                 &
                   institute_id=mycmor%institute_id,                  &
                   parent_experiment_id='N/A',                &
                   parent_experiment_rip='N/A',               &
                   branch_time=0.d0)
              if (error_flag < 0) then
                 write(*,*) 'Error on cmor_dataset!'
                 write(*,*) 'outpath=',mycmor%outpath
                 write(*,*) 'experiment_id=',trim(expids(iexp))
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
                 write(*,*) 'parent_experiment_id=N/A'
                 write(*,*) 'parent_experiment_rip=N/A'
                 write(*,*) 'branch_time=0'
              endif
              !
              ! Add acknowledgements
              !
              error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(acknowledgement))
              !
              ! Add case name, repo tag, and compset
              !
              error_flag = cmor_set_cur_dataset_attribute("cesm_casename",'not applicable')
              error_flag = cmor_set_cur_dataset_attribute("cesm_repotag" ,'not applicable')
              error_flag = cmor_set_cur_dataset_attribute("cesm_compset" ,'not applicable')
              !
              ! Add grid information
              !
              error_flag = cmor_set_cur_dataset_attribute("resolution",'f09_g16 (0.9x1.25_gx1v6)')
              !
              ! Add additional forcing information
              !
              if (mycmor%forcing_note(1:1) /= ' ') error_flag = cmor_set_cur_dataset_attribute("forcing_note",trim(adjustl(mycmor%forcing_note)))
              !
              info_file = 'Info_in.fx'
              inquire(file=trim(info_file),exist=exists)
              if (exists) then
                 call date_and_time(date=pdate,time=ptime)
                 open(30,file=trim(info_file),form='formatted')
                 read(30,'(a)') whoami
                 read(30,'(a)') prochost
                 read(30,'(a)') ccps_rev
                 read(30,'(a)') ccps_date
                 read(30,'(a)') ccps_uuid
                 close(30)
                 error_flag = cmor_set_cur_dataset_attribute("processed_by ",trim(whoami)//" on "//trim(prochost)//" at "//pdate//"-"//ptime)
                 if (error_flag /= 0) then
                    write(*,*) "Error globalMD: ",error_flag
                 endif
                 error_flag = cmor_set_cur_dataset_attribute("processing_code_information ",trim(ccps_rev)//" "//trim(ccps_date)//" "//trim(ccps_uuid))
                 if (error_flag /= 0) then
                    write(*,*) "Error globalMD: ",error_flag
                 endif
              else
                 write(*,*) "Information file: ",trim(info_file)," missing. Stop."
              endif

              ! Define axes via 'cmor_axis'
              !
              table_ids(2) = cmor_load_table('Tables/CMIP5_grids')
              select case (xw(ixw)%entry)
              case ( 'areacello','basin','deptho','hfgeou','sftof','thkcello','volcello')
                 call get_ocn_grid
                 call cmor_set_table(table_ids(2))
                 call define_ocn_axes(xw(ixw)%dims)
              case ( 'mrsofc','sftgif','rootd' )
                 call get_lnd_grid
                 call cmor_set_table(table_ids(2))
                 call define_lnd_axes(xw(ixw)%dims)
              case ( 'areacella','sftlf','orog')
                 call get_atm_grid
                 call cmor_set_table(table_ids(2))
                 call define_atm_axes(xw(ixw)%dims)
              end select
              call cmor_set_table(table_ids(1))
              write(*,*) 'cmor_set_table: ',table_ids(1),table_ids(2),table_ids(3)
              !
              if (xw(ixw)%ncesm_vars == 1) write(original_name,'(a)') xw(ixw)%cesm_vars(1)
              if (xw(ixw)%ncesm_vars == 2) write(original_name,'(a,'','',a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
              if (xw(ixw)%ncesm_vars == 3) write(original_name,'(a,'','',a,'','',a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
              ! 
              ! Make manual alterations so that CMOR works. Silly code!
              !
              select case (xw(ixw)%entry)
              case ('sftlf')
                 var_info(var_found(1,1))%units = '1'
              case ('mrsofc')
                 var_info(var_found(1,1))%units = 'kg m-2'
              case ('areacella')
                 var_info(var_found(1,1))%units = 'm2'
              case ('orog')
                 var_info(var_found(1,1))%units = 'm'
              case ('sftgif')
                 var_info(var_found(1,1))%units = '%'
              case ('volcello')
                 var_info(var_found(1,1))%units = 'm3'
              end select
              !
              write(*,*) 'calling cmor_variable:'
              write(*,*) 'table         = ',trim(mycmor%table_file)
              write(*,*) 'table_entry   = ',trim(xw(ixw)%entry)
              write(*,*) 'dimensions    = ',trim(xw(ixw)%dims)
              write(*,*) 'units         = ',trim(var_info(var_found(1,1))%units)
              write(*,*) 'axis_ids      = ',axis_ids(1:4)
              !              write(*,*) 'missing_value = ',var_info(var_found(1,1))%missing_value
              !              write(*,*) 'missing_value = ',var_info(var_found(1,1))%missing_value
              write(*,*) 'positive      = ',trim(mycmor%positive)
              write(*,*) 'original_name = ',trim(original_name)
              !
              select case (xw(ixw)%entry)
              case ('deptho','areacello')
                 cmor_var_id = cmor_variable(                &
                      table=mycmor%table_file,               &
                      table_entry=xw(ixw)%entry,             &
                      units=var_info(var_found(1,1))%units,  &
                      axis_ids=(/grid_id(1)/),               &
                      missing_value=var_info(var_found(1,1))%missing_value,                   &
                      positive=mycmor%positive,              &
                      original_name=original_name,           &
                      comment=xw(ixw)%comment)
              case ('basin')
                 cmor_var_id = cmor_variable(                &
                      table=mycmor%table_file,               &
                      table_entry=xw(ixw)%entry,             &
                      units=var_info(var_found(1,1))%units,  &
                      axis_ids=(/grid_id(1)/),               &
                      missing_value=var_info(var_found(1,1))%int_missing_value,                &
                      positive=mycmor%positive,              &
                      original_name=original_name,           &
                      comment=xw(ixw)%comment)
              case ('volcello','thkcello')
                 cmor_var_id = cmor_variable(                &
                      table=mycmor%table_file,               &
                      table_entry=xw(ixw)%entry,             &
                      units=var_info(var_found(1,1))%units,  &
                      axis_ids=(/grid_id(1),axis_ids(3)/),   &
                      missing_value=var_info(var_found(1,1))%missing_value,                   &
                      positive=mycmor%positive,              &
                      original_name=original_name,           &
                      comment=xw(ixw)%comment)
              case default
                 cmor_var_id = cmor_variable(                &
                      table=mycmor%table_file,               &
                      table_entry=xw(ixw)%entry,             &
                      units=var_info(var_found(1,1))%units,  &
                      axis_ids=(/axis_ids(1),axis_ids(2)/),   &
                      missing_value=var_info(var_found(1,1))%missing_value,                   &
                      positive=mycmor%positive,              &
                      original_name=original_name,           &
                      comment=xw(ixw)%comment)
                 write(*,*) 'called cmor_variable'
                 write(*,*)
                 write(*,*) 'cmor_var_id   = ',cmor_var_id
                 write(*,*)
                 write(*,*) 'table         = ',trim(mycmor%table_file)
                 write(*,*) 'table_entry   = ',trim(xw(ixw)%entry)
                 write(*,*) 'units         = ',trim(var_info(var_found(1,1))%units)
                 write(*,*) 'axis_ids      = ',axis_ids(1),axis_ids(2)
                 write(*,*) 'missing_value = ',var_info(var_found(1,1))%missing_value
                 write(*,*) 'positive      = ',trim(mycmor%positive)
                 write(*,*) 'original_name = ',trim(original_name)
                 write(*,*) 'comment       = ',trim(xw(ixw)%comment)
              end select
              !
              select case (xw(ixw)%entry)
!!$              case ('deptho','areacello')
!!$                 call get_ocn_grid
!!$                 allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
!!$                 write(*,*) 'TO READ: ',trim(var_info(var_found(1,1))%name)
!!$                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
!!$                 write(*,*) 'READ: ',trim(var_info(var_found(1,1))%name),maxval(indat2a),minval(indat2a)
!!$                 where (kmt == 0)
!!$                    cmordat2d = var_info(var_found(1,1))%missing_value
!!$                 elsewhere
!!$                    cmordat2d = indat2a
!!$                 endwhere
!!$                 error_flag = cmor_write(      &
!!$                      var_id        = cmor_var_id, &
!!$                      data          = cmordat2d)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR writing '',a)') trim(xw(ixw)%entry)
!!$                    stop
!!$                 endif
!!$              case ('thkcello')
!!$                 call get_ocn_grid
!!$                 allocate (cmordat3d(nlons,nlats,nlevs))
!!$                 do k = 1,nlevs
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          cmordat3d(i,j,k) = ocn_t_dz(k)
!!$                       enddo
!!$                    enddo
!!$                 enddo
!!$                 write(*,*) 'COMPUTED thkcello: ',minval(cmordat3d),maxval(cmordat3d)
!!$                 error_flag = cmor_write(      &
!!$                      var_id        = cmor_var_id, &
!!$                      data          = cmordat3d)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR writing '',a)') trim(xw(ixw)%entry)
!!$                    stop
!!$                 endif
!!$                 !
!!$                 !    CMIP5 basins        POP gx1v6 REGION_MASK equivalent values
!!$                 !  0 : global_land           0
!!$                 !  1 : southern_ocean        1
!!$                 !  2 : atlantic_ocean        6
!!$                 !  3 : pacific_ocean         2
!!$                 !  4 : arctic_ocean         10
!!$                 !  5 : indian_ocean          3
!!$                 !  6 : mediterranean_sea     7
!!$                 !  7 : black_sea           -13
!!$                 !  8 : hudson_bay           11
!!$                 !  9 : baltic_sea          -12
!!$                 ! 10 : red_sea               5
!!$                 !
!!$              case ('basin')
!!$                 call get_ocn_grid
!!$                 allocate(region_mask(nlons,nlats),basin(nlons,nlats))
!!$                 write(*,*) 'TO READ: ',trim(var_info(var_found(1,1))%name)
!!$                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,region_mask)
!!$                 basin = var_info(var_found(1,1))%int_missing_value
!!$                 where (region_mask ==   0) basin =  0
!!$                 where (region_mask ==   1) basin =  1
!!$                 where (region_mask ==   6) basin =  2
!!$                 where (region_mask ==   2) basin =  3
!!$                 where (region_mask ==  10) basin =  4
!!$                 where (region_mask ==   3) basin =  5
!!$                 where (region_mask ==   7) basin =  6
!!$                 where (region_mask == -13) basin =  7
!!$                 where (region_mask ==  11) basin =  8
!!$                 where (region_mask == -12) basin =  9
!!$                 where (region_mask ==  -5) basin = 10
!!$                 error_flag = cmor_write(      &
!!$                      var_id        = cmor_var_id, &
!!$                      data          = basin)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR writing '',a)') trim(xw(ixw)%entry)
!!$                    stop
!!$                 endif
!!$              case ('sftgif')
!!$                 call get_lnd_grid
!!$                 allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
!!$                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
!!$                 cmordat2d = indat2a
!!$                 error_flag = cmor_write(      &
!!$                      var_id        = cmor_var_id, &
!!$                      data          = cmordat2d)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR writing '',a)') trim(xw(ixw)%entry)
!!$                    stop
!!$                 endif
!!$              case ('mrsofc')
!!$                 call get_lnd_grid
!!$                 allocate(indat3a(nlons,nlats,nlevs),work3da(nlons,nlats,nlevs),cmordat2d(nlons,nlats))
!!$                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
!!$                 ! Integrate WATSAT over Z
!!$                 do k = 1,nlevs
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (indat3a(i,j,k) /= var_info(var_found(1,1))%missing_value) &
!!$                               work3da(i,j,k) = indat3a(i,j,k)*lnd_dzsoi(i,j,k)
!!$                       enddo
!!$                    enddo
!!$                 enddo
!!$                 cmordat2d = sum(work3da,dim=3)/sum(lnd_levs)
!!$                 error_flag = cmor_write(      &
!!$                      var_id        = cmor_var_id, &
!!$                      data          = cmordat2d)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR writing '',a)') trim(xw(ixw)%entry)
!!$                    stop
!!$                 endif
!!$              case ('sftlf')
!!$                 call get_atm_grid
!!$                 allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
!!$                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
!!$                 cmordat2d = indat2a
!!$                 error_flag = cmor_write(      &
!!$                      var_id        = cmor_var_id, &
!!$                      data          = cmordat2d)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR writing '',a)') trim(xw(ixw)%entry)
!!$                    stop
!!$                 endif
!!$              case ('areacella')
!!$                 call get_atm_grid
!!$                 allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
!!$                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
!!$                 cmordat2d = indat2a*1000.*1000.
!!$                 error_flag = cmor_write(      &
!!$                      var_id        = cmor_var_id, &
!!$                      data          = cmordat2d)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR writing '',a)') trim(xw(ixw)%entry)
!!$                    stop
!!$                 endif
              case ('orog')
                 call get_atm_grid
                 allocate(cmordat2d(nlons,nlats))
                 cmordat2d = 0.
                 where (landfrac .ne. 0)
                    cmordat2d = phis/shr_const_g
                 endwhere
                 where (cmordat2d .lt. 0.)
                    cmordat2d = 0.
                 endwhere
                 error_flag = cmor_write(      &
                      var_id        = cmor_var_id, &
                      data          = cmordat2d)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a)') trim(xw(ixw)%entry)
                    stop
                 endif
!!$              case ('volcello')
!!$                 call get_ocn_grid
!!$                 write(*,*) 'volcello: ',nlons,nlats
!!$                 allocate(cmordat3d(nlons,nlats,nlevs))
!!$                 cmordat3d = var_info(var_found(1,1))%missing_value
!!$                 do k = 1,nlevs
!!$                    do j = 1,nlats
!!$                       do i = 1,nlons
!!$                          if (kmt(i,j).ge.k) then
!!$                             cmordat3d(i,j,k) = (ocn_t_area(i,j)*ocn_t_levs(k))/1.e6
!!$                          else
!!$                             cmordat3d(i,j,k) = var_info(var_found(1,1))%missing_value
!!$                          endif
!!$                       enddo
!!$                    enddo
!!$                 enddo
!!$                 !
!!$                 error_flag = cmor_write(      &
!!$                      var_id        = cmor_var_id, &
!!$                      data          = cmordat3d)
!!$                 if (error_flag < 0) then
!!$                    write(*,'(''ERROR writing '',a)') trim(xw(ixw)%entry)
!!$                    stop
!!$                 endif
              end select
              !
              ! Close all files opened by CMOR.
              !
              error_flag = cmor_close()
              if (error_flag < 0) then
                 write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
              else
                 write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') trim(xw(ixw)%entry),error_flag
              endif
              error_flag   = 0
              var_found    = 0
              original_name= ' '
              !
              if (allocated(indat2a))   deallocate(indat2a)
              if (allocated(cmordat2d)) deallocate(cmordat2d)
              if (allocated(cmordat3d)) deallocate(cmordat3d)
           endif
        enddo var_loop
        call reset_netcdf_var
     enddo xwalk_loop
  enddo expt_loop
end program fx_CMOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! define module with maximum dims
!
module max_parms
  !
  integer,parameter::max_entries   = 150000
  integer,parameter::max_exprmnt   =   1000
  integer,parameter::max_cesm_vars =      6
  integer,parameter::max_ncfiles   =    110
  integer,parameter::max_nchunks   =   1000
!  integer,parameter::max_ntimes    = 400000
  !
end module max_parms
!
! CMOR "table" information
!
module table_info
  use max_parms
  type TableInfo
     character(len=256)::approx_interval,axis,axis_entry,baseURL,bounds_values,cell_measures,cell_methods,cf_version,climatology,cmor_version
     character(len=256)::comment,coords_attrib,dimensions,expt_id_ok,forcings,formula,frequency,generic_levels,index_only,long_name,missing_value
     character(len=256)::modeling_realm,must_call_cmor_grid,must_have_bounds,ok_max_mean_abs,ok_min_mean_abs,out_name,positive,product,project_id
     character(len=256)::required_global_attributes,standard_name,stored_direction,table_date,table_id,tolerance,type,units,valid_max,valid_min,value
     character(len=256)::requested,requested_bounds,variable_entry,z_bounds_factors,z_factors
  end type TableInfo
  integer::num_tab,tab_found
  !
  ! table_ids(1) = Amon, day, OImon, Omon, etc.' table_ids(2) = CMIP5_grids. Always.
  !
  integer,dimension(20)::table_ids
  type(TableInfo),dimension(0:max_entries)::table
end module table_info
!
! CESM experiment information
!
module exp_info
  !
  use max_parms
  !
  type SimInfo
     character(len=256)::case,loc,model_id,expt_id,rip_code,cmip,run_refcase,run_refdate
     character(len=256)::runbegend,mipbegend,grid,compset,repotag,start_fin,mach,dout,forcing
     integer::runbeg,runend,mipbeg,mipend,runlen,miplen
     character(len=256),dimension(10)::icase
  end type SimInfo
  integer::num_exp,exp_found,parent_found
  character(len=256)::case_read,comp_read
  !
  type(SimInfo),dimension(max_entries)::exp
  !
end module exp_info
!
! CMOR arguments information
!
module mycmor_info
  !
  use max_parms
  !
  type CMORInfo
     character(len=256)::table_file,outpath,experiment_id,institution,source,calendar,contact,history
     character(len=256)::references,model_id,forcing,institute_id
     character(len=256)::parent_experiment_id,parent_experiment_rip,positive
     character(len=256)::forcing_note
     character(len=512)::ack_NC,ack_OR,ack_NE
!     character(len=1024)::comment
     character(len=256)::comment
     integer::realization,initialization_method,physics_version
     double precision::branch_time
  end type CMORInfo
  !
  type(CMORInfo)::mycmor
  !
end module mycmor_info
!
! Crosswalk (xwalk) information
!
module xwalk_info
  !
  use max_parms
  !
  type XWInfo
     character(len=256)::table,entry,realm,sname,dims,comment
     character(len=256),dimension(10)::cesm_vars
     integer::ncesm_vars
  end type XWInfo
  integer::num_xw,xw_found
  !
  type(XWInfo),dimension(max_entries)::xw
  !
end module xwalk_info
!
! Input netCDF files information
!
module files_info
  use max_parms
  !
  character(len=256),dimension(max_ncfiles,max_cesm_vars)::ncfile
  character(len=256),dimension(max_ncfiles,max_cesm_vars)::ncfile_nh,ncfile_sh ! For OImon processing
  integer,dimension(max_cesm_vars)::nc_nfiles
  integer,dimension(max_cesm_vars)::nc_nfiles_nh,nc_nfiles_sh ! For OImon processing
  integer,dimension(max_ncfiles,max_cesm_vars)::myncid,var_found,ntimes
  integer,dimension(max_ncfiles,max_cesm_vars)::myncid_nh,myncid_sh ! For OImon processing
!  real   ,dimension(max_ncfiles,max_cesm_vars,max_ntimes)::timevals
  integer::ifile
  logical::all_continue
  !
end module files_info
!
! Grid information
!
module grid_info
  real,dimension(:),    allocatable::atm_lats,atm_lons,atm_levs,atm_levs_bnds,atm_ilevs,atm_ilevs_bnds
  real,dimension(:),    allocatable::cosp_tau,cosp_prs,cosp_ht,cosp_dbze,cosp_sza
  real,dimension(:),    allocatable::atm_plev23,atm_plev17,atm_plev8,atm_plev7,atm_plev7_bnds,atm_plev3
  real,dimension(:,:),  allocatable::atm_lats_bnds,atm_lons_bnds,landfrac,phis
  real,dimension(:,:),  allocatable::cosp_tau_bnds,cosp_prs_bnds,cosp_ht_bnds,cosp_dbze_bnds
  real,dimension(:),    allocatable::a_coeff,b_coeff,a_coeff_bnds,b_coeff_bnds
  real,dimension(:),    allocatable::ah_coeff,bh_coeff
  real,dimension(:,:),  allocatable::ice_t_lats,ice_t_lons,ice_u_lats,ice_u_lons
  real,dimension(:,:,:),allocatable::ice_t_lats_bnds,ice_t_lons_bnds,ice_u_lats_bnds,ice_u_lons_bnds
  real,dimension(:),    allocatable::ocn_t_levs,ocn_trans_lats,ocn_trans_levs,ocn_t_dz
  real,dimension(:,:),  allocatable::ocn_t_levs_bnds,ocn_trans_lats_bnds,ocn_trans_levs_bnds
  real,dimension(:,:),  allocatable::ocn_t_lats,ocn_t_lons,ocn_t_area,ocn_t_hte,ocn_t_htn
  real,dimension(:,:,:),allocatable::ocn_t_lats_bnds,ocn_t_lons_bnds
  real,dimension(:,:),  allocatable::ocn_u_lats,ocn_u_lons
  real,dimension(:,:,:),allocatable::ocn_u_lats_bnds,ocn_u_lons_bnds
  integer,dimension(:,:),allocatable::kmt
  real,dimension(:),    allocatable::lnd_lats,lnd_lons,lnd_levs,lnd_levs_bnds,atm_sites
  real,dimension(:,:,:),allocatable::lnd_zsoi,lnd_dzsoi  ! CLM soil depth (m), CLM soil layer thickness (m)
  real,dimension(:,:),  allocatable::lnd_lats_bnds,lnd_lons_bnds
  real::p0
  integer::nlons,nlats,nlevs,nilevs,nplev23,nplev17,nplev8,nplev7,nplev3
  integer::ncosp_tau,ncosp_prs,ncosp_ht,ncosp_dbze,nsites,ncosp_sza
  integer::naxes,zfactor_id
  integer::nlats_trans,nmoc_z,ntrans_reg,nmoc_comp,ntrans_comp
  integer,dimension(1)::grid_id
  integer,dimension(10)::axis_ids
  character(len=256),dimension(10)::dimnames,dimunits
  character(len=256)::time_units
end module grid_info
!
! Output time chunks information
!
module output_times_info
  use max_parms
  integer,dimension(max_nchunks)::nchunks,tidx1,tidx2
end module output_times_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_atm_grid
  !
  ! Read coordinate information from model into arrays that will be passed to CMOR.
  !
  use counters_netcdf_jfl
  use interfaces_netcdf_jfl
  use definitions_netcdf_jfl
  use grid_info
  use exp_info
  use mycmor_info
  !
  implicit none
  !
  double precision,allocatable,dimension(:)::slon,slat
  character(len=256)::gridfile
  integer::gridid,i,j,k,n,length
  logical::exists
  !
  select case (exp(exp_found)%model_id)
  case ('CESM1-WACCM')
     gridfile = 'atm_grid_f19.nc'
  case ('CESM1-CAM5')
     gridfile = 'atm_grid_cam5_f09.nc'
  case ('CESM1-CAM5.1-FV2')
     gridfile = 'atm_grid_cam5_f19.nc'
  case default
     gridfile = 'atm_grid_cam4_f09.nc'
  end select
  !
  inquire(file=trim(gridfile),exist=exists)
  if (.not.(exists)) then
     write(*,*) 'Cannot find ',trim(gridfile),' - STOPPING.'
     stop
  endif
  !
  call open_cdf(gridid,trim(gridfile),.true.)
  !
  ! Read time-invariant dimensions and variables from gridfile
  !
  call get_dims(gridid)
  !
  do n=1,dim_counter
     length = len_trim(dim_info(n)%name)
     select case (dim_info(n)%name)
     case('lat')
        nlats = dim_info(n)%length
     case('lon') 
        nlons = dim_info(n)%length
     case('lev') 
        nlevs = dim_info(n)%length
     case('ilev') 
        nilevs = dim_info(n)%length
     case('plev23') 
        nplev23 = dim_info(n)%length
     case('plev17') 
        nplev17 = dim_info(n)%length
     case('plev8') 
        nplev8 = dim_info(n)%length
     case('plev7') 
        nplev7 = dim_info(n)%length
     case('plev3') 
        nplev3 = dim_info(n)%length
     case('cosp_tau') 
        ncosp_tau  = dim_info(n)%length
     case('cosp_prs') 
        ncosp_prs  = dim_info(n)%length
     case('cosp_ht') 
        ncosp_ht   = dim_info(n)%length
     case('cosp_dbze') 
        ncosp_dbze = dim_info(n)%length
     case('cosp_sza') 
        ncosp_sza = dim_info(n)%length
     case('ncol') 
        nsites = dim_info(n)%length
     end select
  enddo
  allocate(atm_lons(nlons),atm_lats(nlats),slon(nlons),slat(nlats))
  allocate(atm_lons_bnds(2,nlons),atm_lats_bnds(2,nlats))
  allocate(landfrac(nlons,nlats),phis(nlons,nlats))
  allocate(atm_levs(nlevs),atm_levs_bnds(nlevs+1),atm_sites(nsites))
  allocate(atm_ilevs(nilevs),atm_ilevs_bnds(nilevs+1))
  allocate(atm_plev23(nplev23),atm_plev17(nplev17),atm_plev8(nplev8),atm_plev7(nplev7),atm_plev3(nplev3))
  allocate(a_coeff(nlevs),b_coeff(nlevs),a_coeff_bnds(nlevs+1),b_coeff_bnds(nlevs+1))
  allocate(ah_coeff(nilevs),bh_coeff(nilevs))
  allocate(cosp_tau(ncosp_tau),cosp_tau_bnds(2,ncosp_tau))
  allocate(cosp_sza(ncosp_sza))
  allocate(cosp_prs(ncosp_prs),cosp_prs_bnds(2,ncosp_prs))
  allocate(cosp_ht(ncosp_ht),cosp_ht_bnds(2,ncosp_ht))
  allocate(cosp_dbze(ncosp_dbze),cosp_dbze_bnds(2,ncosp_dbze))
  !
  do i = 1,nsites
     atm_sites(i) = i
  enddo
  !
  call get_vars(gridid)
  call read_var(gridid,'lon'   ,atm_lons)
  call read_var(gridid,'lat'   ,atm_lats)
  call read_var(gridid,'plev23',atm_plev23)
  call read_var(gridid,'plev17',atm_plev17)
  call read_var(gridid,'plev8' ,atm_plev8)
  call read_var(gridid,'plev7' ,atm_plev7)
  call read_var(gridid,'plev3' ,atm_plev3)
  call read_var(gridid,'lev'   ,atm_levs)
  call read_var(gridid,'ilev'  ,atm_levs_bnds)
  call read_var(gridid,'ilev'  ,atm_ilevs)
  call read_var(gridid,'hyam'  ,a_coeff)
  call read_var(gridid,'hyai'  ,a_coeff_bnds)
  call read_var(gridid,'hybm'  ,b_coeff)
  call read_var(gridid,'hybi'  ,b_coeff_bnds)
  call read_var(gridid,'hyai'  ,ah_coeff)
  call read_var(gridid,'hybi'  ,bh_coeff)
  call read_var(gridid,'P0'    ,p0)
  call read_var(gridid,'LANDFRAC',landfrac)
  call read_var(gridid,'PHIS'  ,phis)
  call read_var(gridid,'cosp_tau',cosp_tau)
  call read_var(gridid,'cosp_tau_bnds',cosp_tau_bnds)
  call read_var(gridid,'cosp_sza',cosp_sza)
  call read_var(gridid,'cosp_prs',cosp_prs)
  call read_var(gridid,'cosp_prs_bnds',cosp_prs_bnds)
  cosp_prs      = cosp_prs      * 100
  cosp_prs_bnds = cosp_prs_bnds * 100
  call read_var(gridid,'cosp_ht',cosp_ht)
  call read_var(gridid,'cosp_ht_bnds',cosp_ht_bnds)
  call read_var(gridid,'cosp_dbze',cosp_dbze)
  call read_var(gridid,'cosp_dbze_bnds',cosp_dbze_bnds)
  !
  ! Replace given sza5 values ( 0 15 30 45 60 ) with correct values ( 0 20 40 60 80 )
  !
  cosp_sza = (/0.,20.,40.,60.,80./)
  !
  ! Create atm_ilev_bnds from regular levs
  !
  atm_ilevs_bnds(1) = atm_levs(1)-((atm_levs(2)-atm_levs(1))*0.5)
  do k = 2,nilevs
     atm_ilevs_bnds(k) = atm_levs(k-1)
  enddo
  atm_ilevs_bnds(nilevs+1) = atm_levs(nlevs)+((atm_levs(nlevs)-atm_levs(nlevs-1))*0.5)
!  do k = 1,nlevs
!     write(*,*) atm_levs_bnds(k),atm_levs(k),atm_levs_bnds(k)
!  enddo
  !
  ! Convert Pa values to mb (for vertint)
  !
  p0 = p0 * 0.01
  !
  ! Transfer bounds for lons and lats
  !
  call read_var(gridid,'slon',slon)
  call read_var(gridid,'slat',slat)
  atm_lats_bnds(1,1)     = -90.
  atm_lats_bnds(2,nlats) =  90.
  do j = 1,nlats-1
     atm_lats_bnds(2,j) = slat(j)
  end do
  do j = 2,nlats
     atm_lats_bnds(1,j) =  atm_lats_bnds(2,j-1)
  end do
  !
  atm_lons_bnds(1,    1) = slon(1)
  atm_lons_bnds(2,nlons) = atm_lons(nlons) + ((slon(nlons)-slon(nlons-1))/2.)
  do i = 1,nlons-1
     atm_lons_bnds(2,i) = slon(i+1)
  end do
  do i = 2,nlons
     atm_lons_bnds(1,i) =  atm_lons_bnds(2,i-1)
  end do
  !
  call close_cdf(gridid)
  write(*,'(''ATM grid loaded'')')
  !
!!$  do j = 1,nlats
!!$     write(*,*) atm_lats_bnds(1,j),atm_lats(j),atm_lats_bnds(2,j)
!!$  enddo
!!$  do i = 1,nlons
!!$     write(*,*) atm_lons_bnds(1,i),atm_lons(i),atm_lons_bnds(2,i)
!!$  enddo
end subroutine get_atm_grid
