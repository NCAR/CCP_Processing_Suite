subroutine define_atm_axes(dimensions)
  !
  ! Define all axes that will be needed, as defined by the table entry
  !
  ! 'time' must be last
  !
  use cmor_users_functions
  use files_info
  use table_info
  use grid_info
  use exp_info
  use mycmor_info
  !
  implicit none
  integer,dimension(:),allocatable::location
  character(len=256),intent(in)::dimensions
  !
  integer::i,j,idim,error_flag,ilev
  integer,dimension(0:10)::idxb
  !
  allocate(location(nlons*nlats))
  do i = 1,nlons*nlats
     location(i) = i
  enddo
  !
  ! Parse "dimensions" to find names and how many
  !
  naxes = 1 ; idxb = 1 ; dimnames(:)(1:) = ' '
  do i = 1,len_trim(adjustl(dimensions))
     if (dimensions(i:i) == ' ') then
        idxb(naxes) = i
        naxes = naxes + 1
     endif
  enddo
  !
  dimunits(:)(1:) = ' '
  idxb(0)     = 0
  idxb(naxes) = len_trim(dimensions)+1
  do i = 1,naxes
     dimnames(i)(1:) = dimensions(idxb(i-1)+1:idxb(i)-1)
     select case (dimnames(i))
     case ( 'basin','location','oline','site','typebare','typec3pft','typec4pft','typepdec','typepever','typesdec','typesever','vegtype')
        dimunits(i) = ' '
     case ( 'rho' )
        dimunits(i) = 'kg m-3'
     case ( 'alt40','depth0m','depth100m','height10m','height2m','olayer100m','sdepth','sdepth1','olevel')
        dimunits(i) = 'm'
     case ( 'wv550nm','wv440nm','wv870nm')
        dimunits(i) = 'nm'
     case ( 'p220','p500','p560','p700','p840','plev3','plev7','plev8','plevs')
        dimunits(i) = 'Pa'
     case ( 'longitude')
        dimunits(i) = 'degrees_east'
     case ( 'latitude')
        dimunits(i) = 'degrees_north'
     case ( 'time','time1','time2')
        dimunits(i) = time_units
     case ( 'dbze')
        dimunits(i) = 'dBZ'
     case ( 'sza5')
        dimunits(i) = 'degree'
     case ( 'scatratio','tau','alevel','alev1','alevhalf')
        dimunits(i) = '1'
     case default
        write(*,*) 'Unknown dimension: ',trim(dimnames(i)),' Stopping.'
        stop
     end select
  enddo
  !
  do i = 1,naxes
     write(*,'(''0 define_atm_axes #:'',i4,'' name: '',a32,'' units: '',a32)') i,dimnames(i)(1:32),dimunits(i)(1:32)
  enddo
  axis_ids = 0 ; idim = 1
  !
  ! Define axes
  !
  do i = 1,naxes
     select case(dimnames(i))
!     write(*,'(''L define_atm_axes #:'',i4,'' name: '',a32)') i,dimnames(i)(1:32)
     case ('time','time1','time2')
        select case (mycmor%table_file)
        case ('Tables/CMIP5_Amon','Tables/GeoMIP_Amon','Tables/CMIP5_aero','Tables/CMIP5_cfMon','Tables/PMIP3_Amon',&
              'Tables/AEROCOM-ACC_2D-M','Tables/AEROCOM-ACC_3D-M','Tables/CCMI1_monthly','Tables/HTAP2-monthly')
           axis_ids(idim) = cmor_axis(  &
                table=mycmor%table_file,&
                table_entry=dimnames(i),&
                units=time_units,       &
                interval='30 days')
           write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
           idim = idim + 1
        case ('Tables/CMIP5_day','Tables/GeoMIP_day','Tables/CMIP5_cfDay','Tables/PMIP3_day')
           axis_ids(idim) = cmor_axis(  &
                table=mycmor%table_file,&
                table_entry=dimnames(i),&
                units=time_units,       &
                interval='1 day')
           write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
           idim = idim + 1
        case ('Tables/CMIP5_6hrLev','Tables/CMIP5_6hrPlev','Tables/GeoMIP_6hrLev','Tables/GeoMIP_6hrPlev')
           axis_ids(idim) = cmor_axis(  &
                table=mycmor%table_file,&
                table_entry=dimnames(i),&
                units=time_units,       &
                interval='6 hours')
           write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
           idim = idim + 1
        case ('Tables/CMIP5_3hr','Tables/CMIP5_cf3hr')
           axis_ids(idim) = cmor_axis(  &
                table=mycmor%table_file,&
                table_entry=dimnames(i),&
                units=time_units,       &
                interval='3 hours')
           write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
           idim = idim + 1
        case ('Tables/TAMIP_3hrCurt','Tables/TAMIP_3hrMlev','Tables/TAMIP_3hrPlev','Tables/TAMIP_3hrSlev')
           axis_ids(idim) = cmor_axis(  &
                table=mycmor%table_file,&
                table_entry=dimnames(i),&
                units=time_units,       &
                interval='3 hours')
           write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
           idim = idim + 1
        case ('Tables/TAMIP_sites','Tables/CMIP5_cfSites')
           axis_ids(idim) = cmor_axis(  &
                table=mycmor%table_file,&
                table_entry=dimnames(i),&
                units=time_units,       &
                interval='30 minutes')
           write(*,'('' dimension: '',a,'' idim: '',i4,'' axis_id: '',i4)') trim(dimnames(i)),idim,axis_ids(idim)
           idim = idim + 1
        end select
     end select
  enddo
  !
  do i = 1,naxes
     select case(dimnames(i))
     case ('location')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=SIZE(location),        &
             coord_vals=location)
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('p220')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=1, &
             coord_vals=(/22000/))
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('p560')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=1,                     &
             coord_vals=(/56000/))
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('p840')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=1,                     &
             coord_vals=(/84000/))
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('site')
        write(*,'('' SITE '')')
        select case (mycmor%table_file)
        case ('Tables/CMIP5_cfSites')
!           grid_id(2) = cmor_grid(axis_ids=(/axis_ids(i)/))
!           write(*,*) 'CMOR GRID defined, dimnames,axis_ids,grid_id: ',trim(dimnames(i)),i,grid_id(1)
           zfactor_id = cmor_zfactor(   &
                zaxis_id=ilev,          &
                zfactor_name='ps',      &
!                axis_ids=(/axis_ids(2),axis_ids(1)/), &
                axis_ids=(/grid_id(1),axis_ids(1)/), &
                units='Pa')
        case default
           axis_ids(idim) = cmor_axis(        &
                table=mycmor%table_file,      &
                table_entry=dimnames(i),      &
                units=dimunits(i),            &
                length=SIZE(atm_sites),       &
                coord_vals=atm_sites)
           write(*,'('' dimension: '',a,'' idim: '',i4,'' axis_id: '',i4)') trim(dimnames(i)),idim,axis_ids(idim)
           grid_id(1) = cmor_grid(axis_ids=(/axis_ids(1)/))
           write(*,'('' grid_id defined: '',2i8)') grid_id(1),idim
           idim = idim + 1
        end select
     case ('alt40')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=SIZE(cosp_ht),         &
             coord_vals=cosp_ht,           &
             cell_bounds=cosp_ht_bnds)
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('dbze')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=SIZE(cosp_dbze),       &
             coord_vals=cosp_dbze,         &
             cell_bounds=cosp_dbze_bnds)
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('latitude')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=SIZE(atm_lats),        &
             coord_vals=atm_lats,          &
             cell_bounds=atm_lats_bnds)
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('longitude')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(atm_lons),        &
             units=dimunits(i),            &
             coord_vals=atm_lons,          &
             cell_bounds=atm_lons_bnds)
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('tau')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(cosp_tau),        &
             units=dimunits(i),            &
             coord_vals=cosp_tau,          &
             cell_bounds=cosp_tau_bnds)
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('sza5')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(cosp_sza),        &
             units=dimunits(i),            &
             coord_vals=cosp_sza)
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('plevs')
        select case(exp(exp_found)%model_id)
        case ('CESM1-WACCM')
           select case (mycmor%table_file)
           case ('Tables/CMIP5_Amon')
             axis_ids(idim) = cmor_axis(        &
                  table=mycmor%table_file,      &
                  table_entry=dimnames(i),      &
                  length=SIZE(atm_plev23),       &
                  units=dimunits(i),            &
                  coord_vals=atm_plev23)
             write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
           case ('Tables/CCMI1_monthly')
             axis_ids(idim) = cmor_axis(        &
                  table=mycmor%table_file,      &
                  table_entry=dimnames(i),      &
                  length=SIZE(atm_plev31),       &
                  units=dimunits(i),            &
                  coord_vals=atm_plev31)
             write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
           !
           ! KLUDGE in degenerate longitude as axis_ids(5); for zonal means - Geez
           !
             axis_ids(5) = cmor_axis(           &
                table=mycmor%table_file,      &
                table_entry='longitude',      &
                length=1,                     &
                units='degrees_east',         &
                coord_vals=(/180.0/),         &
                cell_bounds=(/0.0,360.0/))
           end select
        case ('CESM1-CAM4Chem','CESM1-CAM4ChemSD','CESM1-EASALL','CESM1-EURALL','CESM1-GLOALL', &
              'CESM1-NDEALL','CESM1-NAMALL','CESM1-RBUALL','CESM1-SASALL','CESM1-WACCMSD')
           axis_ids(idim) = cmor_axis(        &
                table=mycmor%table_file,      &
                table_entry=dimnames(i),      &
                length=SIZE(atm_plev31),       &
                units=dimunits(i),            &
                coord_vals=atm_plev31)
           write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
           !
           ! KLUDGE in degenerate longitude as axis_ids(5); for zonal means - Geez
           !
           axis_ids(5) = cmor_axis(           &
                table=mycmor%table_file,      &
                table_entry='longitude',      &
                length=1,                     &
                units='degrees_east',         &
                coord_vals=(/180.0/),         &
                cell_bounds=(/0.0,360.0/))
        case default
           axis_ids(idim) = cmor_axis(        &
                table=mycmor%table_file,      &
                table_entry=dimnames(i),      &
                length=SIZE(atm_plev17),       &
                units=dimunits(i),            &
                coord_vals=atm_plev17)
           write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        end select
        idim = idim + 1
     case ('plev8')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(atm_plev8),       &
             units=dimunits(i),            &
             coord_vals=atm_plev8)
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('plev7')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(cosp_prs),        &
             units=dimunits(i),            &
             coord_vals=cosp_prs,          &
             cell_bounds=cosp_prs_bnds)
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('plev3')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(atm_plev3),       &
             units=dimunits(i),            &
             coord_vals=atm_plev3)
        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('alevel')
        ilev = cmor_axis(                        &
             table=mycmor%table_file,            &
             table_entry='standard_hybrid_sigma',&
             length=SIZE(atm_levs),              &
             units='1',                          &
             coord_vals=atm_levs/1000.,          &
             cell_bounds=atm_levs_bnds/1000.)
        error_flag = cmor_zfactor(    &
             zaxis_id=ilev,     &
             zfactor_name='p0',       &
             units='Pa',              &
             zfactor_values=p0Pa)
        error_flag = cmor_zfactor(       &
             zaxis_id=ilev,        &
             zfactor_name='b',           &
             axis_ids= (/ilev/),   &
             units=' ',                  &
             zfactor_values=b_coeff,     &
             zfactor_bounds=b_coeff_bnds)
        error_flag = cmor_zfactor(       &
             zaxis_id=ilev,        &
             zfactor_name='a',           &
             axis_ids= (/ilev/),   &
             units=' ',                  &
             zfactor_values=a_coeff,     &
             zfactor_bounds=a_coeff_bnds)
        select case (mycmor%table_file)
        case ('Tables/TAMIP_sites')
           zfactor_id = cmor_zfactor(       &
                zaxis_id=ilev,        &
                zfactor_name='ps',          &
                axis_ids=(/axis_ids(2),axis_ids(1)/),&
                units='Pa')
        case ('Tables/CMIP5_cfSites')
           zfactor_id = cmor_zfactor(   &
                zaxis_id=ilev,          &
                zfactor_name='ps',      &
                axis_ids=(/grid_id(1)/))
!                axis_ids=(/axis_ids(2)/),&
!                units='Pa')
        case default
           zfactor_id = cmor_zfactor(       &
                zaxis_id=ilev,        &
                zfactor_name='ps',          &
                axis_ids=(/axis_ids(2),axis_ids(3),axis_ids(1)/), &
                units='Pa')
        end select
        axis_ids(idim) = ilev
        write(*,'('' dimension #:'',i4,'' name: '',a,'' defined, id:'',i4)') i,dimnames(i),axis_ids(idim)
        idim = idim + 1
     case ('alevhalf')
        ilev = cmor_axis(                        &
             table=mycmor%table_file,            &
             table_entry='standard_hybrid_sigma',&
             length=SIZE(atm_ilevs),             &
             units='1',                          &
             coord_vals=atm_ilevs/1000.,         &
             cell_bounds=atm_ilevs_bnds/1000.)
        error_flag = cmor_zfactor(    &
             zaxis_id=ilev,     &
             zfactor_name='p0',       &
             units='Pa',              &
             zfactor_values=p0Pa)
        error_flag = cmor_zfactor(       &
             zaxis_id=ilev,        &
             zfactor_name='b',           &
             axis_ids= (/ilev/),   &
             units=' ',                  &
             zfactor_values=b_coeff,     &
             zfactor_bounds=b_coeff_bnds)
        error_flag = cmor_zfactor(       &
             zaxis_id=ilev,        &
             zfactor_name='a',           &
             axis_ids= (/ilev/),   &
             units=' ',                  &
             zfactor_values=a_coeff,     &
             zfactor_bounds=a_coeff_bnds)
        zfactor_id = cmor_zfactor(       &
             zaxis_id=ilev,        &
             zfactor_name='ps',          &
             axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3)/), &
             units='Pa')
        axis_ids(idim) = ilev
        write(*,'('' dimension: '',a,'' defined: '',i4)') 'standard_hybrid_sigma',axis_ids(idim)
        idim = idim + 1
     end select
  enddo
  write(*,'(''CMOR axes defined, axis_ids: '',5i5)') (axis_ids(i),i=1,naxes)
end subroutine define_atm_axes
