subroutine define_ocn_axes(dimensions)
  !
  ! Define all axes that will be needed, as defined by the table entry
  !
  ! 'time' must be last
  !
  use cmor_users_functions
  use table_info
  use grid_info
  use files_info
  use mycmor_info
  use xwalk_info
  !
  implicit none
  real,dimension(:),allocatable::i_indices,j_indices
  character(len=256),intent(in)::dimensions
  !
  integer::i,j,idim,table_id
  integer,dimension(0:10)::idxb
  !
  allocate(i_indices(nlons),j_indices(nlats))
  do i = 1,nlons
     i_indices(i) = i
  enddo
  do j = 1,nlats
     j_indices(j) = j
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
!  do i = 1,naxes
!     write(*,*) 'DIMS: ',dimnames(i)(1:32),' UNITS: ',dimunits(i)(1:32)
!  enddo
  axis_ids = 0 ; idim = 1
  !
  ! Define axes depending on CMIP5 field. Complicated
  !
  select case (xw(xw_found)%entry)
  case ('msftmyz')
     ! Meridional overturning circulation
     do i = 1,naxes
        select case(dimnames(i))
        case ('latitude')
           axis_ids(idim) = cmor_axis(        &
                table=mycmor%table_file,      &
                table_entry=dimnames(i),      &
                units=dimunits(i),            &
                length=nlats_trans,           &
                coord_vals=ocn_trans_lats,    &
                cell_bounds=ocn_trans_lats_bnds)
           write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
           idim = idim + 1
        case ('olevel')
           call cmor_set_table(table_ids(1))
           axis_ids(idim) = cmor_axis(        &
                table=mycmor%table_file,      &
                table_entry='depth_coord',    &
                length=nmoc_z,                &
                units='m',                    &
                coord_vals=ocn_trans_levs,    &
                cell_bounds=ocn_trans_levs_bnds)
           write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
           idim = idim + 1
        case ('basin')
           call cmor_set_table(table_ids(1))
           axis_ids(idim) = cmor_axis(        &
                table_entry=dimnames(i),      &
                coord_vals=(/'atlantic_arctic_ocean','indian_pacific_ocean','global_ocean'/),&
                units='')
           write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
           idim = idim + 1
        case ('time')
           call cmor_set_table(table_ids(1))
           axis_ids(idim) = cmor_axis(        &
                table_entry=dimnames(i),      &
                units=dimunits(i),            &
                interval='30 days')
           write(*,*) 'time defined, axis_id: ',idim,axis_ids(idim)
           idim = idim + 1
        end select ! dimnames(i)
     enddo
  case ('msftbarot','tauuo','tauvo','uo','vo')
     ! U-grid fields: BSF, TAUX, TAUY, UVEL, VVEL
     do i = 1,naxes
        select case(dimnames(i))
        case ('longitude')
           call cmor_set_table(table_ids(2))
           axis_ids(idim) = cmor_axis(        &
                table_entry='i_index',      &
                units='1',&
                coord_vals=i_indices)
           write(*,*) 'longitude defined, axis_id: ',idim,axis_ids(idim)
           idim = idim + 1
        case ('latitude')
           call cmor_set_table(table_ids(2))
           axis_ids(idim) = cmor_axis(        &
                table_entry='j_index',       &
                units='1',&
                coord_vals=j_indices)
           write(*,*) 'latitude defined, axis_id: ',idim,axis_ids(idim)
           grid_id(1) = cmor_grid(                    &
                axis_ids=(/axis_ids(1),axis_ids(2)/), &
                latitude=ocn_u_lats,                    &
                longitude=ocn_u_lons,                   &
                latitude_vertices=ocn_u_lats_bnds,      &
                longitude_vertices=ocn_u_lons_bnds)
           write(*,*) 'CMOR GRID defined, grid_id: ',grid_id(1)
           idim = idim + 1
        case ('olevel')
           call cmor_set_table(table_ids(1))
           axis_ids(idim) = cmor_axis(        &
                table=mycmor%table_file,      &
                table_entry='depth_coord',    &
                length=nlevs,                 &
                units='m',                    &
                coord_vals=ocn_t_levs,          &
                cell_bounds=ocn_t_levs_bnds)
           idim = idim + 1
        case ('time')
           call cmor_set_table(table_ids(1))
           axis_ids(idim) = cmor_axis(        &
                table_entry=dimnames(i),      &
                units=dimunits(i),            &
                interval='30 days')
           write(*,*) 'time defined, axis_id: ',idim,axis_ids(idim)
           idim = idim + 1
        end select ! dimnames(i)
     enddo
  case ( 'so','thetao','tos','sos','volo','hfss','pr','prsn','rlds','rsds','rsntds','agessc','rhopoto','tossq')
     ! T-grid fields
     do i = 1,naxes
        select case(dimnames(i))
        case ('longitude')
           call cmor_set_table(table_ids(2))
           axis_ids(idim) = cmor_axis(        &
                table_entry='i_index',      &
                units='1',&
                coord_vals=i_indices)
           write(*,*) 'longitude defined, axis_id: ',idim,axis_ids(idim)
           idim = idim + 1
        case ('latitude')
           call cmor_set_table(table_ids(2))
           axis_ids(idim) = cmor_axis(        &
                table_entry='j_index',       &
                units='1',&
                coord_vals=j_indices)
           write(*,*) 'latitude  defined, axis_id: ',idim,axis_ids(idim)
           idim = idim + 1
           grid_id(1) = cmor_grid(                    &
                axis_ids=(/axis_ids(1),axis_ids(2)/), &
                latitude=ocn_t_lats,                    &
                longitude=ocn_t_lons,                   &
                latitude_vertices=ocn_t_lats_bnds,      &
                longitude_vertices=ocn_t_lons_bnds)
           write(*,'('' grid defined: '',i4,'' axis_ids: '',2i10)') grid_id(1),axis_ids(1),axis_ids(2)
        case ('olevel')
           call cmor_set_table(table_ids(1))
           axis_ids(idim) = cmor_axis(        &
                table=mycmor%table_file,      &
                table_entry='depth_coord',    &
                length=nlevs,                 &
                units=dimunits(i),            &
                coord_vals=ocn_t_levs,        &
                cell_bounds=ocn_t_levs_bnds)
           write(*,*) 'olevel   defined, axis_id: ',idim,axis_ids(idim)
           idim = idim + 1
        case ('time')
           call cmor_set_table(table_ids(1))
           axis_ids(idim) = cmor_axis(        &
                table_entry=dimnames(i),      &
                units=dimunits(i),            &
                interval='30 days')
           write(*,'('' dimension: '',a,'' defined: '',i4,'' units: '',a)') 'time',axis_ids(idim),trim(dimunits(i))
           idim = idim + 1
        end select ! dimnames(i)
     enddo
  case ( 'areacello','basin','deptho','hfgeou','sftof','thkcello','volcello')
     ! T-grid fields
     do i = 1,naxes
        select case(dimnames(i))
        case ('longitude')
           call cmor_set_table(table_ids(2))
           axis_ids(idim) = cmor_axis(        &
                table_entry='i_index',      &
                units='1',&
                coord_vals=i_indices)
           write(*,*) 'longitude defined, axis_id: ',idim,axis_ids(idim)
           idim = idim + 1
        case ('latitude')
           call cmor_set_table(table_ids(2))
           axis_ids(idim) = cmor_axis(        &
                table_entry='j_index',       &
                units='1',&
                coord_vals=j_indices)
           write(*,*) 'latitude  defined, axis_id: ',idim,axis_ids(idim)
           idim = idim + 1
           grid_id(1) = cmor_grid(                    &
                axis_ids=(/axis_ids(1),axis_ids(2)/), &
                latitude=ocn_t_lats,                    &
                longitude=ocn_t_lons,                   &
                latitude_vertices=ocn_t_lats_bnds,      &
                longitude_vertices=ocn_t_lons_bnds)
           write(*,'('' grid defined: '',i4,'' axis_ids: '',2i10)') grid_id(1),axis_ids(1),axis_ids(2)
        case ('olevel')
           call cmor_set_table(table_ids(1))
           axis_ids(idim) = cmor_axis(        &
                table=mycmor%table_file,      &
                table_entry='depth_coord',    &
                length=nlevs,                 &
                units=dimunits(i),            &
                coord_vals=ocn_t_levs,        &
                cell_bounds=ocn_t_levs_bnds)
           write(*,*) 'olevel   defined, axis_id: ',idim,axis_ids(idim)
           idim = idim + 1
        end select ! dimnames(i)
     enddo
  end select ! xw(ixw)%entry
  write(*,*) 'CMOR axes defined, axis_ids: ',(axis_ids(i),i=1,naxes)
end subroutine define_ocn_axes
