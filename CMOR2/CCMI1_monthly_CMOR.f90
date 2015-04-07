program CCMI_monthly_CMOR
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
  real,parameter::spval  = 1.e20
  real,parameter::pi     = 4.*(atan(1.))
  real,parameter::grav   = 9.81         ! gravity, m s-1
  real,parameter::rearth = 6.37122e6    ! earth radius, m
  real,parameter::rdair  = 287.04       ! Dry air gas constant, J k-1 kg-1
  real,parameter::cscale = 1.4          ! Carbon scaling factor
  real,parameter::iscale = 60./68.1          ! Isoprene scaleing 
  real,parameter::avogn  = 6.022e23
  !
  ! Molecular weights, kg/mole
  !
  real,parameter::mw_dryair = 28.97e-3
!  real,parameter::mw_so4    = 96.06
  real,parameter::mw_ch4  =  16.
  real,parameter::mw_co   =  28.
  real,parameter::mw_h2o  =  18.
  real,parameter::mw_nox  =  46.0055
  real,parameter::mw_o3   =  47.9982
  real,parameter::mw_pan  = 121.
  real,parameter::mw_h2o2 =  34.0147
  real,parameter::mw_hno3 =  63.0128
  real,parameter::mw_so2 =  64.066
  real,parameter::mw_nh3 =  17.031
  real,parameter::mw_oh   =  17.0073
  real,parameter::mw_SO4  =  97.078
  ! SOA components molecular weights
  ! varsoa  = (/"SOAI","SOAT","SOAB","SOAX","SOAM"/)
  ! mwsoa_c = (/  60.0550  ,    84.077   ,    72.0660  ,    96.0880  , 120.1100   /) WRONG
  ! mwsoa   = (/ 136.141400,   141.141800,   127.116000,   155.167600, 200.226000 /) CORRECT
  !  real,parameter::mw_soai =  60.055
  !  real,parameter::mw_soat =  84.077
  !  real,parameter::mw_soab =  72.066
  !  real,parameter::mw_soax =  96.088
  !  real,parameter::mw_soam = 120.110
  real,parameter::mw_soai = 136.141400
  real,parameter::mw_soat = 141.141800
  real,parameter::mw_soab = 127.116000
  real,parameter::mw_soax = 155.167670
  real,parameter::mw_soam = 200.226000
  !
  !  uninitialized variables used in communicating with CMOR:
  !
  integer::error_flag,cmor_var_id
  real,dimension(:) ,allocatable::area_wt
  real,dimension(:,:),allocatable::indat2a,indat2b,indat2c,indat2d,indat2e,indat2f,indat2g,cmordat2d
  real,dimension(:,:),allocatable::indat2_01,indat2_02,indat2_03,indat2_04,indat2_05,indat2_06,indat2_07,indat2_08,indat2_09,indat2_10
  real,dimension(:,:),allocatable::indat2_11,indat2_12,indat2_13,indat2_14,indat2_15,indat2_16,indat2_17,indat2_18,indat2_19,indat2_20
  real,dimension(:,:) ,allocatable::psdata
  real,dimension(:,:,:),allocatable::indat3a,indat3b,indat3c,indat3d,indat3e,indat3f,indat3g,cmordat3d,work3da,work3db,zonave
  real,dimension(:,:,:),allocatable::pshybrid,psdelta,pshybrid_mid,tdata,rho
  !
  double precision,dimension(:) ,allocatable::time
  double precision,dimension(:,:),allocatable::time_bnds
  double precision,dimension(1) ::tval
  double precision,dimension(2,1)::tbnd
  !
  ! Other variables
  !
  character(len=256)::exp_file,xwalk_file,table_file,svar,tstr,original_name,logfile,cmor_filename
  integer::i,j,k,m,n,tcount,it,ivar,length,iexp,jexp,ixw,ilev,ic,jfile,ii,ij,ik,lon_count
  logical::does_exist
  !
  ! GO
  !
  mycmor%table_file = 'monthly'
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
  ! Provides information on relationship between CMOR variables and
  ! model variables
  !
  xwalk_file = 'xwalk_'//trim(exp(exp_found)%cmip)//'_'//trim(mycmor%table_file)
  call load_xwalk(xwalk_file)
  !
  ! Get table information
  !
  mycmor%table_file = 'Tables/'//trim(exp(exp_found)%cmip)//'_'//trim(mycmor%table_file)
  inquire(file=mycmor%table_file,exist=does_exist)
  if (.not.(does_exist)) then
     write(*,*) 'Cannot find',trim(mycmor%table_file),'. Dying.'
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
  ! Step through CMOR table entries to see what CESM fields we can read and in process, and if so, do it
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
        if (xw(ixw)%ncesm_vars ==  1) write(original_name,'(a)') xw(ixw)%cesm_vars(1)
        if (xw(ixw)%ncesm_vars ==  2) write(original_name,'(a,'','',a)')    (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars ==  3) write(original_name,'(2(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars ==  4) write(original_name,'(3(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars ==  5) write(original_name,'(4(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars ==  6) write(original_name,'(5(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars ==  7) write(original_name,'(6(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars ==  8) write(original_name,'(7(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars ==  9) write(original_name,'(8(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars == 10) write(original_name,'(9(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        if (xw(ixw)%ncesm_vars == 18) write(original_name,'(18(a,'',''),a)') (trim(xw(ixw)%cesm_vars(i)),i=1,xw(ixw)%ncesm_vars)
        !
        ! Modify units as necessary to accomodate udunits' inability to convert
        !
        select case (xw(ixw)%entry)
        case('tauu','tauv','hfss','rlut','rlutcs','hfls','rlus','rsus','rsuscs','rsut','rsutcs','mc',&
             'emidms','emiso2','emiss','emibc','emidust','emiso4',&
             'emico','emibisop','eminh3','emipom','emivoc','eminox')
           mycmor%positive = 'up'
        case('drydms','drydust','drypoa','dryso2','dryso4','drysoa','dryss','dryoa',&
             'wetdust','wetoa','wetso4','wetbc','wetnh4',&
             'dryhno3','drynh4','dryno2','drynoy','dryo3','drybc')
           mycmor%positive = 'down'
        case ('loadoa','loadbc','loaddust','loadss','loadso4','airmass')
           var_info(var_found(1,1))%units = 'kg m-2'
        case ('wetss')
           mycmor%positive = 'up'
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
        case ('rlds','rldscs','rsds','rsdscs','rsdt','rtmt')
           mycmor%positive = 'down'
        case ('clt','ci','sci','cod')
           var_info(var_found(1,1))%units = '1'
        case ('hurs','cl','sic')
           var_info(var_found(1,1))%units = '%'
        case('od550ss')
           var_info(var_found(1,1))%units = '1'
        case('prc','pr','prsn')
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
!           mycmor%positive = 'down'
        case('wethno3','wetnh3','wetso2')
           mycmor%positive = 'down'
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
        case('emiisop')
           var_info(var_found(1,1))%units = 'kg (C) m-2 s-1'
           mycmor%positive = 'up'
        case ('chegpso4','do3chm','chepsoa','lso3chm','tpo3chm')
           var_info(var_found(1,1))%units = 'kg m-2 s-1'
        case ('photo1d','jno2')
           var_info(var_found(1,1))%units = 's-1'
        case('lossch4','lossco','o3loss','o3prod','ohloss',&
             'losso1dviah2o','losso3viaho2','losso3viaoh',&
             'lossrcoo2viano2','lossro2viaho2',&
             'lossro2viano','lossro2viano3','lossro2viaro2',&
              'prodh2o2viaho2','prodhno3viano2oh','prodhpx',&
              'prodo1d','prodo3viaho2','prodo3viaro2','prodoh')
           var_info(var_found(1,1))%units = 'mole m-3 s-1'
        case ('vmraoanh')
           var_info(var_found(1,1))%units = 'year'
        case ('ptp')
           var_info(var_found(1,1))%units = 'Pa'
        case ('reffclwtop','ztp','dh')
           var_info(var_found(1,1))%units = 'm'
        case ('tatp')
           var_info(var_found(1,1))%units = 'K'
        case ('toz')
           var_info(var_found(1,1))%units = 'DU'
        case ('acceldivf','accelgw','accelnogw','accelogw')
           var_info(var_found(1,1))%units = 'm s-2'
        end select
!
        write(*,*) 'calling cmor_variable:'
        write(*,*) 'table =',trim(mycmor%table_file)
        write(*,*) 'table_entry =',trim(xw(ixw)%entry)
        write(*,*) 'dimensions = ',trim(xw(ixw)%dims)
        write(*,*) 'units = ',var_info(var_found(1,1))%units(1:20)
        write(*,*) 'axis_ids      = ',axis_ids(1:naxes)
        write(*,*) 'missing_value = ',var_info(var_found(1,1))%missing_value
        write(*,*) 'positive      = ',trim(mycmor%positive)
        write(*,*) 'original_name = ',trim(original_name)
        !
        select case (xw(ixw)%entry)
        case ('aoa','vmraoanh','chegpso4','cl','cli','clw',&
             'dh','do3chm','jno2','lossch4','lossco',&
             'mcu','mmrbc','mmrdust','mmroa','mmrsoa','mmrss',&
             'o3loss','o3prod','ohloss','photo1d','pilev','pmlev',&
             'prodh2o2viaho2','ta','ua','va','hus','zg',&
             'vmrc2h2','vmrc2h6','vmrch2o','vmrch3ccl3','vmrch3cn','vmrch4','vmrco25','vmrco50',&
             'vmrco','vmrdms','vmre90','vmre90n','vmre90s','vmrh2o','vmrhcl','vmrhcn','vmrhno3',&
             'vmrisop','vmrn2o','vmrnh50','vmrnh50w','vmrnh5','vmrno2','vmrno','vmro3','vmro3s',&
             'vmroh','vmrpan','vmrsf6','vmrso2','vmrso2t','vmrst8025',&
             'losso1dviah2o','losso3viaho2','losso3viaoh','lossrcoo2viano2','lossro2viaho2',&
              'lossro2viano','lossro2viano3','lossro2viaro2',&
              'prodhno3viano2oh','prodhpx','prodo1d','prodo3viaho2','prodo3viaro2','prodoh',&
              'mmraernh4','mmraerno3','mmraerso4','emilnox')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(2),axis_ids(3),axis_ids(4),axis_ids(1)/),  &
                missing_value=spval,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        case ('zmbrcl','zmbr','zmbro','zmbrono2','zmbry','zmc2br2f4','zmcbrclf2','zmcbrf3','zmccl4',&
              'zmcf2cl2','zmcfcl3','zmch2br2','zmch2o','zmch3br','zmch3ccl2f','zmch3ccl3','zmch3cclf2',&
              'zmch3cl','zmch3ooh','zmch4','zmchbr3','zmchclf2','zmcl2o2','zmcl','zmclo','zmclono2',&
              'zmcly','zmco','zmh2','zmh2o2','zmh2o','zmhbr','zmhcl','zmhno3','zmhno4','zmho2','zmhobr',&
              'zmhocl','zmmnstrage','zmn2o5','zmn2o','zmn','zmno2','zmno','zmnoy','zmo3','zmoclo','zmoh',&
              'zmta','zmtnt','zmua','zmva','zmzg','acceldivf','accelgw','accelnogw','accelogw','airmass')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(5),axis_ids(3),axis_ids(4),axis_ids(1)/),  &
                missing_value=spval,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
        case ('ps')
           cmor_var_id = cmor_variable(                            &
                table=mycmor%table_file,                           &
                table_entry=xw(ixw)%entry,                         &
                units=var_info(var_found(1,1))%units,                &
                axis_ids=(/axis_ids(2),axis_ids(3),axis_ids(1)/), &
                missing_value=spval,&
                positive=mycmor%positive,                          &
                original_name=original_name,                       &
                comment=xw(ixw)%comment)
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
        case ('clt','od550aer','od550bc','od550oa','ps',&
               'rlds','rlutcs','rsdscs','rsds','rsdt','clivi',&
               'abs550aer','dryhno3','drynh3','drynh4','dryno2','drynoy','dryo3',&
               'emico','emiisop','eminh3','eminox','emibisop',&
              'hfls','hfss','od550so4','reffclwtop',&
                'wetso4','wetbc','wetnh4','ztp',&
              'drydms','dryso2','dryso4','emidms','emiso2','emiso4',&
              'dryss','wetss','drydust','wetdust',&
              'ptp','tatp','flashrate','lwp','swclrefc')
           !
           ! No change or multiply by a factor
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
                
                 if (var_info(var_found(1,1))%name.eq.'SFISOP') then 
                   indat2a = indat2a*iscale
                 endif
                 if (var_info(var_found(1,1))%name.eq.'FLASHFRQ') then 
                   indat2a = indat2a/60./sum(area_wt)
                 endif
                 if (var_info(var_found(1,1))%name.eq.'TGCLDLWP') then 
                   indat2a = indat2a*1000.
                 endif
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
        case ('tos')
           !
           ! Use 'landfrac' to mask 'TS'
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
                 cmordat2d = indat2a
                 where (landfrac /= 0.)
                    cmordat2d = spval
                 endwhere
                 tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                 error_flag = cmor_write(          &
                      var_id        = cmor_var_id, &
                      data          = cmordat2d,     &
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
        case ('sic')
           !
           ! Convert ICEFRAC from 0..1 to percent
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
                 cmordat2d = indat2a*100.
                 tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                 error_flag = cmor_write(          &
                      var_id        = cmor_var_id, &
                      data          = cmordat2d,     &
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
        case ('pr')
           !
           ! PRECT, unit change from m s-1 to kg m-2 s-1
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
        case ('emibc','rlus','drybc')
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
        case ('drypoa','emioa')
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
        case ('emipom')
           !
           ! Add two fields, 1.4 scale factor applied
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
        case ('od550dust','od550ss','emidust','emiss')
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
                    cmordat2d = spval
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
        case ('emivoc','emibvoc')
           !
           ! Sum 19 fields
           !
           allocate(indat2_01(nlons,nlats),indat2_02(nlons,nlats),indat2_03(nlons,nlats),indat2_04(nlons,nlats),indat2_05(nlons,nlats))
           allocate(indat2_06(nlons,nlats),indat2_07(nlons,nlats),indat2_08(nlons,nlats),indat2_09(nlons,nlats),indat2_10(nlons,nlats))
           allocate(indat2_11(nlons,nlats),indat2_12(nlons,nlats),indat2_13(nlons,nlats),indat2_14(nlons,nlats),indat2_15(nlons,nlats),indat2_16(nlons,nlats),indat2_17(nlons,nlats),indat2_18(nlons,nlats),indat2_19(nlons,nlats))
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
                 call read_var(myncid(1, 1),var_info(var_found(1, 1))%name,indat2_01)
                 call read_var(myncid(1, 2),var_info(var_found(1, 2))%name,indat2_02)
                 call read_var(myncid(1, 3),var_info(var_found(1, 3))%name,indat2_03)
                 call read_var(myncid(1, 4),var_info(var_found(1, 4))%name,indat2_04)
                 call read_var(myncid(1, 5),var_info(var_found(1, 5))%name,indat2_05)
                 call read_var(myncid(1, 6),var_info(var_found(1, 6))%name,indat2_06)
                 call read_var(myncid(1, 7),var_info(var_found(1, 7))%name,indat2_07)
                 call read_var(myncid(1, 8),var_info(var_found(1, 8))%name,indat2_08)
                 call read_var(myncid(1, 9),var_info(var_found(1, 9))%name,indat2_09)
                 call read_var(myncid(1,10),var_info(var_found(1,10))%name,indat2_10)
                 call read_var(myncid(1,11),var_info(var_found(1,11))%name,indat2_11)
                 call read_var(myncid(1,12),var_info(var_found(1,12))%name,indat2_12)
                 call read_var(myncid(1,13),var_info(var_found(1,13))%name,indat2_13)
                 call read_var(myncid(1,14),var_info(var_found(1,14))%name,indat2_14)
                 call read_var(myncid(1,15),var_info(var_found(1,15))%name,indat2_15)
                 call read_var(myncid(1,16),var_info(var_found(1,16))%name,indat2_16)
                 call read_var(myncid(1,17),var_info(var_found(1,17))%name,indat2_17)
                 call read_var(myncid(1,18),var_info(var_found(1,18))%name,indat2_18)
                 !
                 where ((indat2_01 /= 1.e36).and.(indat2_02 /= 1.e36).and.(indat2_03 /= 1.e36).and.(indat2_04 /= 1.e36).and.&
                        (indat2_05 /= 1.e36).and.(indat2_06 /= 1.e36).and.(indat2_07 /= 1.e36).and.(indat2_08 /= 1.e36).and.&
                        (indat2_09 /= 1.e36).and.(indat2_10 /= 1.e36).and.(indat2_11 /= 1.e36).and.(indat2_12 /= 1.e36).and.&
                        (indat2_13 /= 1.e36).and.(indat2_14 /= 1.e36).and.(indat2_15 /= 1.e36).and.(indat2_16 /= 1.e36).and.&
                        (indat2_17 /= 1.e36).and.(indat2_18 /= 1.e36))
                   indat2_01 = indat2_01 * (12./ 32.)
                   indat2_02 = indat2_02 * (3*12./ 58.)
                   indat2_03 = indat2_03 * (2* 12./ 44.)
                   indat2_04 = indat2_04 * (12./ 30.)
                   indat2_05 = indat2_05 * (12./ 28.)
                   indat2_06 = indat2_06 * (2*12./ 30.)
                   indat2_07 = indat2_07 * (3*12./ 44.)
                   indat2_08 = indat2_08 * (2*12./ 28.)
                   indat2_09 = indat2_09 * (3*12./ 42.)
                   indat2_10 = indat2_10 * (2*12./ 46.)
                   indat2_11 = indat2_11 * (10*12./ 136.)
                   indat2_12 = indat2_12 * (5*12./ 72.)
                   indat2_13 = indat2_13 * (4*12./ 56.)
                   indat2_14 = indat2_14 * (6*12./ 92.)
                   indat2_15 = indat2_15 * (4*12./ 72.)
                   indat2_16 = indat2_16 * (12./ 46.)
                   indat2_17 = indat2_17 * (2*12./ 60.)
                   indat2_18 = indat2_18 * (5*12./ 68.1)

                   cmordat2d = indat2_01+indat2_02+indat2_03+indat2_04+indat2_05+indat2_06+indat2_07+indat2_08+indat2_09+indat2_10+indat2_11+indat2_12+indat2_13+indat2_14+indat2_15+indat2_16+indat2_17+indat2_18
                 elsewhere
                    cmordat2d = spval
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
           deallocate(indat2_01,indat2_02,indat2_03,indat2_04,indat2_05)
           deallocate(indat2_06,indat2_07,indat2_08,indat2_09,indat2_10)
           deallocate(indat2_11,indat2_12,indat2_13,indat2_14,indat2_15,indat2_16,indat2_17,indat2_18)

        case ('wetoa')
           !
           ! Sum six fields
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats),indat2d(nlons,nlats),indat2e(nlons,nlats),indat2f(nlons,nlats))
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
                 call read_var(myncid(1,6),var_info(var_found(1,6))%name,indat2f)
                 !
                 where ((indat2a /= 1.e36).and.(indat2b /= 1.e36).and.(indat2c /= 1.e36).and.(indat2d /= 1.e36).and.(indat2e /= 1.e36).and.(indat2f /= 1.e36).and.(indat2g /= 1.e36))
                    cmordat2d = 1.4*(indat2a) + indat2b + indat2c + indat2d + indat2e + indat2f + indat2g
                 elsewhere
                    cmordat2d = spval
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

        case ('dryoa')
           !
           ! Sum seven fields
           !
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats),indat2c(nlons,nlats),indat2d(nlons,nlats),indat2e(nlons,nlats),indat2f(nlons,nlats),indat2g(nlons,nlats))
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
                 call read_var(myncid(1,6),var_info(var_found(1,6))%name,indat2f)
                 call read_var(myncid(1,7),var_info(var_found(1,7))%name,indat2g)
                 !
                 where ((indat2a /= 1.e36).and.(indat2b /= 1.e36).and.(indat2c /= 1.e36).and.(indat2d /= 1.e36).and.(indat2e /= 1.e36).and.(indat2f /= 1.e36).and.(indat2g /= 1.e36))
                    cmordat2d = 1.4*(indat2a + indat2b) + indat2c + indat2d + indat2e + indat2f + indat2g 
                 elsewhere
                    cmordat2d = spval
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
                    cmordat2d = spval
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
        case ('loaddust','loadss')
           !
           ! Sum four fields integrated over Z
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs))
           allocate(indat3c(nlons,nlats,nlevs),indat3d(nlons,nlats,nlevs))
           allocate(work3da(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats),pshybrid(nlons,nlats,nlevs),psdelta(nlons,nlats,nlevs))
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
                 work3da   = spval
                 cmordat2d = spval
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                 call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
                 call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat3d)
                 call read_var(myncid(1,5),'PS',psdata)
                 !
                 call pres_hybrid_ccm(psdata,pshybrid,nlons,nlats,nlevs)
                 do k = 1,nlevs-1
                    psdelta(:,:,k)=pshybrid(:,:,k+1)-pshybrid(:,:,k)
                 enddo
                 !
                 work3da   = (indat3a + indat3b + indat3c + indat3d)*(psdelta/grav)
                 cmordat2d = sum(work3da,dim=3)
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
                 work3da   = spval
                 cmordat2d = spval
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                 call read_var(myncid(1,3),'PS',psdata)
                 !
                 call pres_hybrid_ccm(psdata,pshybrid,nlons,nlats,nlevs)
                 do k = 1,nlevs-1
                    psdelta(:,:,k)=pshybrid(:,:,k+1)-pshybrid(:,:,k)
                 enddo
                 !
                 work3da   = (indat3a + indat3b)*(psdelta/grav)
                 write(*,*) 'CM3 0',minval(work3da),maxval(work3da)
                 cmordat2d = sum(work3da,dim=3)
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
        case ('mmroa')
           !
           ! Sum seven fields, leave on model levels
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs),indat3c(nlons,nlats,nlevs))
           allocate(indat3d(nlons,nlats,nlevs),indat3e(nlons,nlats,nlevs),indat3f(nlons,nlats,nlevs))
           allocate(indat3g(nlons,nlats,nlevs),cmordat3d(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    cmordat3d = spval
                    call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                    call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
                    call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat3d)
                    call read_var(myncid(1,5),var_info(var_found(1,5))%name,indat3e)
                    call read_var(myncid(1,6),var_info(var_found(1,6))%name,indat3f)
                    call read_var(myncid(1,7),var_info(var_found(1,7))%name,indat3g)
                    call read_var(myncid(1,8),'PS',psdata)
                    !
                    cmordat3d = indat3a+indat3b+indat3c+indat3d+indat3e+indat3f+indat3g
                    !
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,   &
                         data          = cmordat3d,     &
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
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo
        case ('mmrsoa')
           !
           ! Sum five fields, leave on model levels
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs),indat3c(nlons,nlats,nlevs))
           allocate(indat3d(nlons,nlats,nlevs),indat3e(nlons,nlats,nlevs))
           allocate(cmordat3d(nlons,nlats,nlevs),psdata(nlons,nlats))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    cmordat3d = spval
                    call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                    call read_var(myncid(1,2),var_info(var_found(1,2))%name,indat3b)
                    call read_var(myncid(1,3),var_info(var_found(1,3))%name,indat3c)
                    call read_var(myncid(1,4),var_info(var_found(1,4))%name,indat3d)
                    call read_var(myncid(1,5),var_info(var_found(1,5))%name,indat3e)
                    call read_var(myncid(1,6),'PS',psdata)
                    !
                    cmordat3d = indat3a+indat3b+indat3c+indat3d+indat3e
                    !
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,   &
                         data          = cmordat3d,     &
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
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo

        case ('toz','lso3chm','tpo3chm','cod')
           !
           ! One field integrated over Z and scaled
           !
           allocate(indat3a(nlons,nlats,nlevs))
           allocate(work3da(nlons,nlats,nlevs))
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
                 work3da   = spval
                 cmordat2d = spval
                 call read_var(myncid(1,1),var_info(var_found(1,1))%name,indat3a)
                 !
                 !
                  cmordat2d = sum(indat3a,dim=3)
                  if (var_info(var_found(ifile,1))%name.eq.'O3') then
                      cmordat2d = cmordat2d * 2.1e+22 / 2.69e16 
                    endif
                  if (var_info(var_found(ifile,1))%name.eq.'DO3CHM_LMS') then
                      cmordat2d = cmordat2d / sum(area_wt) 
                    endif
                  if (var_info(var_found(ifile,1))%name.eq.'DO3CHM_TRP') then
                      cmordat2d = cmordat2d / sum(area_wt) 
                    endif
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
        case ('zmbrcl','zmbr','zmbro','zmbrono2','zmbry','zmc2br2f4','zmcbrclf2','zmcbrf3','zmccl4',&
              'zmcf2cl2','zmcfcl3','zmch2br2','zmch2o','zmch3br','zmch3ccl2f','zmch3ccl3','zmch3cclf2',&
              'zmch3cl','zmch3ooh','zmch4','zmchbr3','zmchclf2','zmcl2o2','zmcl','zmclo','zmclono2',&
              'zmcly','zmco','zmh2','zmh2o2','zmh2o','zmhbr','zmhcl','zmhno3','zmhno4','zmho2','zmhobr',&
              'zmhocl','zmmnstrage','zmn2o5','zmn2o','zmn','zmno2','zmno','zmnoy','zmo3','zmoclo','zmoh',&
              'zmta','zmtnt','zmua','zmva','zmzg','airmass','accelogw')
           !
           ! Just one field, interpolated to plevs then zonally averaged
           !
           allocate(indat3a(nlons,nlats,nlevs))
           allocate(cmordat3d(nlons,nlats,nplev31))
           allocate(psdata(nlons,nlats),zonave(1,nlats,nplev31))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,psdata)
                    where (abs(indat3a) > spval)
                       indat3a = spval
                    endwhere
                    where (abs(psdata) > spval)
                       psdata = spval
                    elsewhere
                       psdata = psdata * 0.01
                    endwhere
                    !
                    cmordat3d = spval
                    zonave    = 0.
                    !
                    ! Do vertical interpolation to pressure levels
                    !
                    call vertint(indat3a,cmordat3d,atm_levs,atm_plev31*0.01,psdata,spval,nlons,nlats,nlevs,nlevs+1,nplev31)
                    !
                    ! Zonal average
                    !
                    do ik = 1,nplev31
                       do ij = 1,nlats
                          lon_count = 0
                          do ii = 1,nlons
                             if (cmordat3d(ii,ij,ik) /= spval) then
                                zonave(1,ij,ik) = zonave(1,ij,ik) + cmordat3d(ii,ij,ik)
                                lon_count = lon_count + 1
                             endif
                          enddo
                          if (lon_count == 0) then
                             zonave(1,ij,ik) = spval
                          else
                             zonave(1,ij,ik) = zonave(1,ij,ik)/lon_count
                          endif
                       enddo
                    enddo
!                    write(*,*) minval(cmordat3d),maxval(cmordat3d,mask=cmordat3d/=spval),minval(zonave),maxval(zonave,mask=zonave/=spval)
                    !
                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id, &
                         data          = zonave,      &
                         ntimes_passed = 1,         &
                         time_vals     = tval,      &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                       stop
                    endif
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                 endif
              enddo
           enddo

        case ('acceldivf')
           !
           ! Just one field, interpolated to plevs 
           !
           allocate(indat3a(nlons,nlats,nlevs))
           allocate(cmordat3d(nlons,nlats,nplev31))
           allocate(psdata(nlons,nlats))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,psdata)
                    where (abs(indat3a) > spval)
                       indat3a = spval
                    endwhere
                    where (abs(psdata) > spval)
                       psdata = spval
                    elsewhere
                       psdata = psdata * 0.01
                    endwhere
                    !
                    cmordat3d = spval
                    !
                    ! Do vertical interpolation to pressure levels
                    !
                    call vertint(indat3a,cmordat3d,atm_levs,atm_plev31*0.01,psdata,spval,nlons,nlats,nlevs,nlevs+1,nplev31)
                    !
                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id, &
                         data          = cmordat3d,      &
                         ntimes_passed = 1,         &
                         time_vals     = tval,      &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                       stop
                    endif
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                 endif
              enddo
           enddo

        case ('accelnogw')
           !
           ! Two fields, interpolated to plevs 
           !
           allocate(indat3a(nlons,nlats,nlevs))
           allocate(indat3b(nlons,nlats,nlevs))
           allocate(indat3c(nlons,nlats,nlevs))
           allocate(cmordat3d(nlons,nlats,nplev31))
           allocate(psdata(nlons,nlats))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,1))%name,indat3b)
                    call read_var(myncid(ifile,3),var_info(var_found(ifile,2))%name,psdata)
                    indat3c = indat3a + indat3b
                    where (abs(psdata) > spval)
                       psdata = spval
                    elsewhere
                       psdata = psdata * 0.01
                    endwhere
                    !
                    cmordat3d = spval
                    !
                    ! Do vertical interpolation to pressure levels
                    !
                    call vertint(indat3c,cmordat3d,atm_levs,atm_plev31*0.01,psdata,spval,nlons,nlats,nlevs,nlevs+1,nplev31)
                    !
                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id, &
                         data          = cmordat3d,      &
                         ntimes_passed = 1,         &
                         time_vals     = tval,      &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                       stop
                    endif
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                 endif
              enddo
           enddo
        case ('accelgw')
           !
           ! Three fields, interpolated to plevs and zonal average
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs),indat3c(nlons,nlats,nlevs))
           allocate(cmordat3d(nlons,nlats,nplev31))
           allocate(psdata(nlons,nlats),zonave(1,nlats,nplev31))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat3b)
                    call read_var(myncid(ifile,3),var_info(var_found(ifile,3))%name,indat3c)
                    call read_var(myncid(ifile,4),var_info(var_found(ifile,4))%name,psdata)
                    indat3d = indat3a + indat3b + indat3c
                    where (abs(psdata) > spval)
                       psdata = spval
                    elsewhere
                       psdata = psdata * 0.01
                    endwhere
                    !
                    cmordat3d = spval
                    zonave    = 0.
                    !
                    ! Do vertical interpolation to pressure levels
                    !
                    call vertint(indat3c,cmordat3d,atm_levs,atm_plev31*0.01,psdata,spval,nlons,nlats,nlevs,nlevs+1,nplev31)
                    !
                    ! Zonal average
                    !
                    do ik = 1,nplev31
                       do ij = 1,nlats
                          lon_count = 0
                          do ii = 1,nlons
                             if (cmordat3d(ii,ij,ik) /= spval) then
                                zonave(1,ij,ik) = zonave(1,ij,ik) + cmordat3d(ii,ij,ik)
                                lon_count = lon_count + 1
                             endif
                          enddo
                          if (lon_count == 0) then
                             zonave(1,ij,ik) = spval
                          else
                             zonave(1,ij,ik) = zonave(1,ij,ik)/lon_count
                          endif
                       enddo
                    enddo
                    !
                    tval(1) = time(it) ; tbnd(1,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id, &
                         data          = zonave,      &
                         ntimes_passed = 1,         &
                         time_vals     = tval,      &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                       stop
                    endif
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic-1,cmor_filename(1:128)
                 endif
              enddo
           enddo


        case ('aoa','vmraoanh','cl','cli','clw','hus','jno2','mcu','photo1d','pilev','pmlev',&
              'ta','ua','va','vmrc2h2','vmrc2h6','vmrch2o','vmrch3ccl3','vmrch3cn',&
              'vmrch4','vmrco','vmrco25','vmrco50','vmrdms','vmre90','vmre90n','vmre90s','vmrh2o',&
              'vmrhcl','vmrhcn','vmrhno3','vmrisop','vmrn2o','vmrnh5','vmrnh50','vmrnh50w','vmrno',&
              'vmrno2','vmro3','vmro3s','vmroh','vmrpan','vmrsf6','vmrso2','vmrso2t','vmrst8025','zg',&
              'mmraernh4','mmraerno3','mmraerso4','dh')
           !
           ! Just one field, leave on model levels
           !
           allocate(indat3a(nlons,nlats,nlevs))
           allocate(cmordat3d(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),'PS',psdata)
                    !
                    if (var_info(var_found(ifile,1))%name.eq.'NH4') then
                      indat3a = 18./28.97*indat3a   
                    endif
                    if (var_info(var_found(ifile,1))%name.eq.'AOA_NH') then
                      indat3a = 1.e07*indat3a   
                    endif
                    if (var_info(var_found(ifile,1))%name.eq.'Z3') then
                     do k = 1,nlevs-1
                       indat3a(:,:,k)=indat3a(:,:,k+1)-indat3a(:,:,k)
                     enddo
                    endif
                    cmordat3d = indat3a
                    !
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,   &
                         data          = cmordat3d,     &
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
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo

        case ('do3chm')
           !
           ! Just one field, leave on model levels
           !
           allocate(indat3a(nlons,nlats,nlevs))
           allocate(cmordat3d(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),'PS',psdata)
                    !
                    indat3a = indat3a / sum(area_wt)
                    cmordat3d = indat3a
                    !
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,   &
                         data          = cmordat3d,     &
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
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo
         
        case ('chegpso4')
           !
           ! Just one field, scaled and converted
           !
           allocate(indat3a(nlons,nlats,nlevs),tdata(nlons,nlats,nlevs),rho(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nlevs),pshybrid(nlons,nlats,nlevs),pshybrid_mid(nlons,nlats,nlevs),psdelta(nlons,nlats,nlevs))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),'PS',psdata)
                    call read_var(myncid(ifile,3),'T',tdata)
                    ! Convert to kg m-2 s-1
                    indat3a = indat3a * (mw_dryair * 1.e-3 / avogn )
                    !
                    call pres_hybrid_ccm(psdata,pshybrid,nlons,nlats,nlevs)
                    do k = 1,nlevs-1
                       psdelta(:,:,k)=pshybrid(:,:,k+1)-pshybrid(:,:,k)
                    enddo
                    !
                    call pres_hybrid_mid_ccm(psdata,pshybrid_mid,nlons,nlats,nlevs)
                    rho = pshybrid_mid/(287.04*tdata)
                    cmordat3d = indat3a*mw_so4 * psdelta/grav/rho
                    !
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,   &
                         data          = cmordat3d,     &
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
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo
        case ('lossch4','lossco','o3loss',&
              'losso1dviah2o','losso3viaho2','losso3viaoh','lossrcoo2viano2','lossro2viaho2',&
              'lossro2viano','lossro2viano3','lossro2viaro2', 'o3prod','ohloss',&
              'prodh2o2viaho2','prodhno3viano2oh','prodhpx','prodo1d','prodo3viaho2','prodo3viaro2','prodoh')
           !
           ! Just one field, scaled and converted
           !
           allocate(indat3a(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nlevs),pshybrid(nlons,nlats,nlevs),psdelta(nlons,nlats,nlevs))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),'PS',psdata)
                    ! Convert to kg m-2 s-1
                    indat3a = indat3a * 1.e-6 / avogn 
                    !
                    call pres_hybrid_ccm(psdata,pshybrid,nlons,nlats,nlevs)
                    do k = 1,nlevs-1
                       psdelta(:,:,k)=pshybrid(:,:,k+1)-pshybrid(:,:,k)
                    enddo
                    !
                    cmordat3d = indat3a*(psdelta/grav)
                    !
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,   &
                         data          = cmordat3d,     &
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
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo
        case ('emilnox')
          !
           ! scale LNOx 
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           ! nc1:LNO_PROD 
           !
           allocate(indat3a(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
           allocate(tdata(nlons,nlats,nlevs))
           allocate(pshybrid(nlons,nlats,nlevs),pshybrid_mid(nlons,nlats,nlevs),psdelta(nlons,nlats,nlevs))
           allocate(cmordat3d(nlons,nlats,nlevs))
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
            do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),'PS',psdata)
                    call read_var(myncid(ifile,3),'T',tdata)
                    ! Convert to kg m-2 s-1
                    indat3a = indat3a * (14. / 30.) * 1.e-03/avogn
                    !
                    call pres_hybrid_ccm(psdata,pshybrid,nlons,nlats,nlevs)
                    do k = 1,nlevs-1
                       psdelta(:,:,k)=pshybrid(:,:,k+1)-pshybrid(:,:,k)
                    enddo
                    !
                     call pres_hybrid_mid_ccm(psdata,pshybrid_mid,nlons,nlats,nlevs)
                     rho = pshybrid_mid/(287.04*tdata)
                    !
                    cmordat3d = (indat3a*(psdelta/grav))/rho
                    !
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,   &
                         data          = cmordat3d,     &
                         ntimes_passed = 1,         &
                         time_vals     = tval,      &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                       stop
                    endif
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                 !
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo


        case ('wethno3','wetnh3','wetso2')
           !
           ! Integrate over Z  and scale
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           allocate(indat3a(nlons,nlats,nlevs),tdata(nlons,nlats,nlevs),rho(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
           allocate(work3da(nlons,nlats,nlevs),pshybrid(nlons,nlats,nlevs),pshybrid_mid(nlons,nlats,nlevs),psdelta(nlons,nlats,nlevs))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    cmordat2d = spval
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),'PS',psdata)
                    call read_var(myncid(ifile,3),'T',tdata)
                    ! Convert to kg m-2 s-1
                    !
                    call pres_hybrid_ccm(psdata,pshybrid,nlons,nlats,nlevs)
                    do k = 1,nlevs-1
                       psdelta(:,:,k)=pshybrid(:,:,k+1)-pshybrid(:,:,k)
                    enddo
                    !
                    call pres_hybrid_mid_ccm(psdata,pshybrid_mid,nlons,nlats,nlevs)
                    rho = pshybrid_mid/(287.04*tdata)

                    if (var_info(var_found(ifile,1))%name.eq.'DTWR_HNO3') then
                     indat3a = indat3a*mw_hno3/mw_dryair
                    endif
                    if (var_info(var_found(ifile,1))%name.eq.'DTWR_NH3') then
                     indat3a = indat3a*mw_nh3/mw_dryair
                    endif
                    if (var_info(var_found(ifile,1))%name.eq.'DTWR_SO2') then
                     indat3a = indat3a*mw_so2/mw_dryair
                    endif

                    work3da = (indat3a*(psdelta/grav))/rho
                    cmordat2d = sum(work3da,dim=3)
                    !
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,   &
                         data          = cmordat2d,     &
                         ntimes_passed = 1,         &
                         time_vals     = tval,      &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                       stop
                    endif
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                 !
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo

        case ('chepsoa')
           !
           ! Sum five fields, integrate over Z  and scale
           ! Non-vertically interpolated data; pass straight through, but include 'PS' as required, and
           ! break up into nicely-sized chunks along time
           !
           ! nc1:SOAM nc2:SOAI nc3:SOAT nc4:SOAB nc5:SOAX
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs),indat3c(nlons,nlats,nlevs))
           allocate(indat3d(nlons,nlats,nlevs),indat3e(nlons,nlats,nlevs),tdata(nlons,nlats,nlevs),rho(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
           allocate(work3da(nlons,nlats,nlevs),pshybrid(nlons,nlats,nlevs),pshybrid_mid(nlons,nlats,nlevs),psdelta(nlons,nlats,nlevs))
           allocate(cmordat2d(nlons,nlats))
           allocate(indat2a(nlons,nlats),indat2b(nlons,nlats))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    cmordat2d = spval
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat3b)
                    call read_var(myncid(ifile,3),var_info(var_found(ifile,3))%name,indat3c)
                    call read_var(myncid(ifile,4),var_info(var_found(ifile,4))%name,indat3d)
                    call read_var(myncid(ifile,5),var_info(var_found(ifile,5))%name,indat3e)
                    call read_var(myncid(ifile,6),var_info(var_found(ifile,6))%name,indat2a)
                    call read_var(myncid(ifile,7),var_info(var_found(ifile,7))%name,indat2b)
                    call read_var(myncid(ifile,8),'PS',psdata)
                    call read_var(myncid(ifile,9),'T',tdata)
                    ! Convert to kg m-2 s-1
                    indat3a = indat3a * (mw_soam / mw_dryair)
                    indat3b = indat3b * (mw_soai / mw_dryair)
                    indat3c = indat3c * (mw_soat / mw_dryair)
                    indat3d = indat3d * (mw_soab / mw_dryair)
                    indat3e = indat3e * (mw_soax / mw_dryair)
                    !
                    call pres_hybrid_ccm(psdata,pshybrid,nlons,nlats,nlevs)
                    do k = 1,nlevs-1
                       psdelta(:,:,k)=pshybrid(:,:,k+1)-pshybrid(:,:,k)
                    enddo
                    !
                    call pres_hybrid_mid_ccm(psdata,pshybrid_mid,nlons,nlats,nlevs)
                    rho = pshybrid_mid/(287.04*tdata)
                    !
                    work3da = ((indat3a + indat3b + indat3c + indat3d + indat3e)*(psdelta/grav))/rho
                    cmordat2d = sum(work3da,dim=3)
                    cmordat2d = cmordat2d + 1.4 * (indat2a + indat2b)
                    !
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,   &
                         data          = cmordat2d,     &
                         ntimes_passed = 1,         &
                         time_vals     = tval,      &
                         time_bnds     = tbnd)
                    if (error_flag < 0) then
                       write(*,'(''ERROR writing '',a,'' T# '',i6)') trim(xw(ixw)%entry),it
                       stop
                    endif
                 enddo
                 write(*,'(''DONE writing '',a,'' T# '',i6,'' chunk# '',i6)') trim(xw(ixw)%entry),it-1,ic
                 !
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo
        case ('mmrdust','mmrss')
           !
           ! Sum four fields, leave on model levels
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs))
           allocate(indat3c(nlons,nlats,nlevs),indat3d(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nlevs))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat3b)
                    call read_var(myncid(ifile,3),var_info(var_found(ifile,3))%name,indat3c)
                    call read_var(myncid(ifile,4),var_info(var_found(ifile,4))%name,indat3d)
                    call read_var(myncid(ifile,5),'PS',psdata)
                    !
                    cmordat3d = indat3a + indat3b + indat3c + indat3d
                    !
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,   &
                         data          = cmordat3d,     &
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
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo
        case ('mmrbc')
           !
           ! Sum two fields, leave on model levels
           !
           allocate(indat3a(nlons,nlats,nlevs),indat3b(nlons,nlats,nlevs))
           allocate(psdata(nlons,nlats))
           allocate(cmordat3d(nlons,nlats,nlevs))
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
           do ifile = 1,nc_nfiles(1)
              select case(ntimes(ifile,1))
              case default
                 nchunks(ifile) = 1
                 tidx1(1:nchunks(ifile)) = 1
                 tidx2(1:nchunks(ifile)) = ntimes(ifile,1)
              end select
              write(*,'(''# chunks '',i3,'':'',10((i6,''-'',i6),1x))') nchunks(ifile),(tidx1(ic),tidx2(ic),ic=1,nchunks(ifile))
              do ic = 1,nchunks(ifile)
                 do it = tidx1(ic),tidx2(ic)
                    time_counter = it
                    call read_var(myncid(ifile,1),var_info(var_found(ifile,1))%name,indat3a)
                    call read_var(myncid(ifile,2),var_info(var_found(ifile,2))%name,indat3b)
                    call read_var(myncid(ifile,3),'PS',psdata)
                    !
                    cmordat3d = indat3a + indat3b
                    !
                    tval(1) = time(it) ; tbnd(ifile,1) = time_bnds(1,it) ; tbnd(2,1) = time_bnds(2,it)
                    error_flag = cmor_write(        &
                         var_id        = cmor_var_id,   &
                         data          = cmordat3d,     &
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
                 cmor_filename(1:) = ' '
                 error_flag = cmor_close(var_id=cmor_var_id,file_name=cmor_filename,preserve=1)
                 if (error_flag < 0) then
                    write(*,'(''ERROR close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                    stop
                 else
                    write(*,'(''GOOD close chunk: '',i6,'' of '',a)') ic,cmor_filename(1:128)
                 endif
              enddo
           enddo
        end select
        !
        if (allocated(cmordat2d)) deallocate(cmordat2d)
        if (allocated(indat2a))   deallocate(indat2a)
        if (allocated(indat2b))   deallocate(indat2b)
        if (allocated(indat2c))   deallocate(indat2c)
        if (allocated(indat2d))   deallocate(indat2d)
        if (allocated(indat2e))   deallocate(indat2e)
        if (allocated(indat2f))   deallocate(indat2f)
        if (allocated(indat2g))   deallocate(indat2g)
        if (allocated(indat3a))   deallocate(indat3a)
        if (allocated(indat3b))   deallocate(indat3b)
        if (allocated(indat3c))   deallocate(indat3c)
        if (allocated(indat3d))   deallocate(indat3d)
        if (allocated(indat3e))   deallocate(indat3e)
        if (allocated(indat3f))   deallocate(indat3f)
        if (allocated(indat3g))   deallocate(indat3g)
        if (allocated(psdata))    deallocate(psdata)
        if (allocated(psdelta))   deallocate(psdelta)
        if (allocated(pshybrid))  deallocate(pshybrid)
        if (allocated(rho))       deallocate(rho)
        if (allocated(tdata))     deallocate(tdata)
        if (allocated(work3da))   deallocate(work3da)
        if (allocated(work3db))   deallocate(work3db)
        if (allocated(pshybrid_mid))   deallocate(pshybrid_mid)
        !
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
end program CCMI_monthly_CMOR
