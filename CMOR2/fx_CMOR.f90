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
  integer,parameter::nexp = 37
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  integer::error_flag,cmor_var_id
  real,dimension(:)    ,allocatable::dx,dy,dlat
  real,dimension(:,:)  ,allocatable::indat2a,cmordat2d
  real,dimension(:,:,:),allocatable::cmordat3d
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile
  character(len=512)::acknowledgement
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,itab,ixw
  real::spval,pi,rad,rr,rearth,dlon
  character(len=15),dimension(nexp)::expids
  character(len=256)::whoami,prochost,ccps_rev,ccps_date,ccps_uuid,info_file
  character(len=10)::pdate,ptime
  logical::exists
  !
  expids = (/&
       ! Long-term fully-coupled
       'piControl     ','1pctCO2       ','abrupt4xCO2   ',&
       'historical    ','historicalExt ','historicalGHG ','historicalMisc','historicalNat ',&
       'rcp26         ','rcp45         ','rcp60         ','rcp85         ',&
       ! Paleo runs
       'lgm           ','midHolocene   ','past1000      ',&
       ! AMIP, aquaplanet and other atm-only runs
       'amip          ','amip4K        ','amip4xCO2     ','amipFuture    ',&
       'aqua4K        ','aqua4xCO2     ','aquaControl   ',&
       'sst2030       ','sstClim       ','sstClim4xCO2  ','sstClimAerosol','sstClimSulfate',&
       ! ESM runs
       'esmControl    ','esmFdbk1      ','esmFdbk2      ','esmFixClim1   ','esmFixClim2   ','esmHistorical ','esmrcp85      ',&
       ! Decadal predictions
       'decadalXXXX   ','noVolcXXXX    ','volcIn2010    '/)
  !
  ! GO!
  !
  mycmor%table_file = 'CMIP5_fx'
  call load_table_info
  !
  ! Get "crossxwalk" (xwalk) information
  !   Provides information on relationship between CMOR variables and
  !   model variables
  !
  xwalk_file = 'xwalk_fx.txt'
  call load_xwalk(xwalk_file)
  !
  ! Get grid information
  !
  call get_atm_grid
  call get_lnd_grid
  call get_ocn_grid
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
  mycmor%model_id      = 'CCSM4'
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
     table_loop: do itab = 1,num_tab
        xwalk_loop: do ixw = 1,num_xw
           mycmor%positive = ' '
           time_counter = 0
           var_counter  = 0
           error_flag   = 0
           var_found    = 0
           original_name= ' '
           !
           ! The meaty part
           !
           if (xw(ixw)%entry == table(itab)%variable_entry) then
              select case (xw(ixw)%entry)
              case ( 'areacella','sftlf')
                 ncfile(1,1) = 'atm_grid_f09.nc'
                 xw(ixw)%ncesm_vars = 1
              case ( 'mrsofc','sftgif','rootd' )
                 ncfile(1,1) = 'lnd_grid_f09.nc'
                 xw(ixw)%ncesm_vars = 1
              case ( 'areacello','basin','deptho','hfgeou','sftof','thkcello','volcello')
                 ncfile(1,1) = 'ocn_grid_gx1.nc'
                 xw(ixw)%ncesm_vars = 1
              end select
              write(*,'(''OPENING: '',a,'' myncid: '',i10)') trim(ncfile(1,1)),myncid(1,1)
              call open_cdf(myncid(1,1),trim(ncfile(1,1)),.true.)
              call get_dims(myncid(1,1))
              call get_vars(myncid(1,1))
              !
              do ivar = 1,xw(ixw)%ncesm_vars
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
              enddo
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
              !
              ! Define axes via 'cmor_axis'
              !
              table_ids(2) = cmor_load_table('Tables/CMIP5_grids')
              select case (xw(ixw)%entry)
              case ( 'areacello','basin','deptho','hfgeou','sftof','thkcello','volcello')
                 call get_ocn_grid
                 call cmor_set_table(table_ids(2))
                 call define_ocn_axes(table(itab)%dimensions)
                 write(*,*) 'define_ocn_axes: ',trim(table(itab)%dimensions)
              case ( 'mrsofc','sftgif','rootd' )
                 call get_lnd_grid
                 call cmor_set_table(table_ids(2))
                 call define_lnd_axes(table(itab)%dimensions)
              write(*,*) 'define_lnd_axes: ',trim(table(itab)%dimensions)
              case ( 'areacella','sftlf')
                 call get_atm_grid
                 call cmor_set_table(table_ids(2))
                 call define_atm_axes(table(itab)%dimensions)
                 write(*,*) 'define_atm_axes: ',trim(table(itab)%dimensions)
              end select
              call cmor_set_table(table_ids(1))
              write(*,*) 'cmor_set_table: ',table_ids(1),table_ids(2),table_ids(3)
              !
              if (xw(ixw)%ncesm_vars == 1) write(original_name,'(a)') xw(ixw)%cesm_vars(1)
              if (xw(ixw)%ncesm_vars == 2) write(original_name,'(a,'','',a)') (trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
              if (xw(ixw)%ncesm_vars == 3) write(original_name,'(a,'','',a,'','',a)') (trim(xw(ixw)%cesm_vars(ivar)),ivar=1,xw(ixw)%ncesm_vars)
              ! 
              ! Make manual alterations so that CMOR works. Silly code!
              !
              select case (xw(ixw)%entry)
              case ('sftlf')
                 var_info(var_found(1,1))%units = '1'
              case ('sftgif')
                 var_info(var_found(1,1))%units = '%'
              case ('volcello')
                 var_info(var_found(1,1))%units = 'm3'
              end select
              !
              spval=var_info(var_found(1,1))%missing_value
              !
              write(*,*) 'calling cmor_variable:'
              write(*,*) 'table         = ',trim(mycmor%table_file)
              write(*,*) 'table_entry   = ',trim(xw(ixw)%entry)
              write(*,*) 'dimensions    = ',trim(table(itab)%dimensions)
              write(*,*) 'units         = ',trim(var_info(var_found(1,1))%units)
              write(*,*) 'axis_ids      = ',axis_ids(1:4)
              write(*,*) 'missing_value = ',var_info(var_found(1,1))%missing_value
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
                      missing_value=spval,                   &
                      positive=mycmor%positive,              &
                      original_name=original_name,           &
                      comment=xw(ixw)%comment)
              case ('volcello')
                 cmor_var_id = cmor_variable(                &
                      table=mycmor%table_file,               &
                      table_entry=xw(ixw)%entry,             &
                      units=var_info(var_found(1,1))%units,  &
                      axis_ids=(/grid_id(1),axis_ids(3)/),   &
                      missing_value=spval,                   &
                      positive=mycmor%positive,              &
                      original_name=original_name,           &
                      comment=xw(ixw)%comment)
              case default
                 cmor_var_id = cmor_variable(                &
                      table=mycmor%table_file,               &
                      table_entry=xw(ixw)%entry,             &
                      units=var_info(var_found(1,1))%units,  &
                      axis_ids=(/axis_ids(1),axis_ids(2)/),   &
                      missing_value=spval,                   &
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
                 write(*,*) 'missing_value = ',spval
                 write(*,*) 'positive      = ',trim(mycmor%positive)
                 write(*,*) 'original_name = ',trim(original_name)
                 write(*,*) 'comment       = ',trim(xw(ixw)%comment)
              end select
              !
              select case (xw(ixw)%entry)
              case ('deptho','areacello')
                 call get_ocn_grid
                 allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
                 write(*,*) 'TO READ: ',trim(var_info(var_found(1,1))%name)
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 write(*,*) 'READ: ',trim(var_info(var_found(1,1))%name),maxval(indat2a),minval(indat2a)
                 where (kmt == 0)
                    cmordat2d = spval
                 elsewhere
                    cmordat2d = indat2a
                 endwhere
                 error_flag = cmor_write(      &
                      var_id        = cmor_var_id, &
                      data          = cmordat2d)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              case ('sftgif')
                 call get_lnd_grid
                 allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 cmordat2d = indat2a
                 error_flag = cmor_write(      &
                      var_id        = cmor_var_id, &
                      data          = cmordat2d)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              case ('sftlf')
                 call get_atm_grid
                 allocate(indat2a(nlons,nlats),cmordat2d(nlons,nlats))
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 cmordat2d = indat2a
                 error_flag = cmor_write(      &
                      var_id        = cmor_var_id, &
                      data          = cmordat2d)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              case ('areacella')
                 call get_atm_grid
!                 allocate(cmordat2d(nlons,nlats),dx(nlons),dy(nlats),dlat(nlats))
                 allocate(cmordat2d(nlons,nlats),dx(nlats),dy(nlats),dlat(nlats))
                 !
                 rad    = 4.*atan(1.)/180.0
                 rearth = 6.37122e3
                 rr     = rearth*rad
                 !
                 dlon = rr*(atm_lons(2) - atm_lons(1))
                 write(*,*) 'dlon*cos(atm_lats*rad): ',size(dlon*cos(atm_lats*rad))
                 do j = 1,nlats
                    dx(j) = dlon*cos((atm_lats(j)*rad))
                 enddo
                 do j = 2,nlats
                    dlat(i) = atm_lats(i) - atm_lats(i-1)
                 enddo
                 dlat(1) = dlat(nlats)
                 dy      = dlat*rr
                 do j = 1,nlats
                    do i = 1,nlons
                       cmordat2d(i,j) = dx(i) * dy(j)
                    enddo
                 enddo
                 write(*,*) 'A N,X,S: ',minval(cmordat2d),maxval(cmordat2d),sum(cmordat2d)
                 error_flag = cmor_write(      &
                      var_id        = cmor_var_id, &
                      data          = cmordat2d)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              case ('volcello')
                 call get_ocn_grid
                 write(*,*) 'volcello: ',nlons,nlats
                 allocate(indat2a(nlons,nlats),cmordat3d(nlons,nlats,nlevs))
                 write(*,*) 'TO READ: ',trim(var_info(var_found(1,1))%name)
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat2a)
                 cmordat3d = spval
                 do k = 1,nlevs
                    do j = 1,nlats
                       do i = 1,nlons
                          if (kmt(i,j).ge.k) then
                             cmordat3d(i,j,k) = (indat2a(i,j)*ocn_t_levs(k))/(100. * 100. * 100.)
                          else
                             cmordat3d(i,j,k) = spval
                          endif
                       enddo
                    enddo
                 enddo
                 !
                 error_flag = cmor_write(      &
                      var_id        = cmor_var_id, &
                      data          = cmordat3d)
                 if (error_flag < 0) then
                    write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                    stop
                 endif
              end select
              !
              ! Close all files opened by CMOR.
              !
              error_flag = cmor_close()
              if (error_flag < 0) then
                 write(*,'(''ERROR cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
              else
                 write(*,'('' GOOD cmor_close of : '',a,'' flag: '',i6)') ,trim(xw(ixw)%entry),error_flag
              endif
              time_counter = 0
              var_counter  = 0
              error_flag   = 0
              var_found    = 0
              original_name= ' '
              !
              if (allocated(indat2a))   deallocate(indat2a)
              if (allocated(cmordat2d)) deallocate(cmordat2d)
              if (allocated(cmordat3d)) deallocate(cmordat3d)
           endif
        enddo xwalk_loop
     enddo table_loop
  enddo expt_loop
end program fx_CMOR
