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
  idxb(0)     = 0
  idxb(naxes) = len_trim(dimensions)+1
  do i = 1,naxes
     dimnames(i)(1:) = dimensions(idxb(i-1)+1:idxb(i)-1)
     do j = 1,num_tab
        if (dimnames(i) == table(j)%axis_entry) dimunits(i) = table(j)%units
     enddo
     if (dimnames(i) == 'time') dimunits(i) = time_units
  enddo
  !
  do i = 1,naxes
     write(*,*) 'DIMS: ',dimnames(i)(1:32),' UNITS: ',dimunits(i)(1:32)
  enddo
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
  case default
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
                latitude=ocn_lats,                    &
                longitude=ocn_lons,                   &
                latitude_vertices=ocn_lats_bnds,      &
                longitude_vertices=ocn_lons_bnds)
           write(*,*) 'CMOR GRID defined, grid_id: ',grid_id(1)
           idim = idim + 1
        case ('olevel')
           call cmor_set_table(table_ids(1))
           axis_ids(idim) = cmor_axis(        &
                table=mycmor%table_file,      &
                table_entry='depth_coord',    &
                length=nlevs,                 &
                units='m',                    &
                coord_vals=ocn_levs,          &
                cell_bounds=ocn_levs_bnds)
           idim = idim + 1
        case ('time')
           call cmor_set_table(table_ids(1))
           axis_ids(idim) = cmor_axis(        &
                table_entry=dimnames(i),      &
                units=dimunits(i),            &
!                length=ntimes(1,1),           &
                interval='30 days')
           write(*,*) 'time defined, axis_id: ',idim,axis_ids(idim)
           idim = idim + 1
        end select ! dimnames(i)
     enddo
  end select ! xw(ixw)%entry
  write(*,*) 'CMOR axes defined, axis_ids: ',(axis_ids(i),i=1,naxes)
end subroutine define_ocn_axes
