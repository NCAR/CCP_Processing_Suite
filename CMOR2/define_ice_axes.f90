subroutine define_ice_axes(dimensions)
  !
  ! Define all axes that will be needed, as defined by the table entry
  !
  ! 'time' must be last
  !
  use cmor_users_functions
  use table_info
  use grid_info
  use mycmor_info
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
             latitude=ice_lats,                    &
             longitude=ice_lons,                   &
             latitude_vertices=ice_lats_bnds,      &
             longitude_vertices=ice_lons_bnds)
        write(*,*) 'CMOR GRID defined, grid_id: ',grid_id(1)
        idim = idim + 1
     case ('time')
        call cmor_set_table(table_ids(1))
        axis_ids(idim) = cmor_axis(        &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=ntimes,                &
             interval='30 days')
        write(*,*) 'time defined, axis_id: ',idim,axis_ids(idim)
        idim = idim + 1
     end select
  enddo
  write(*,*) 'CMOR axes defined, axis_ids: ',(axis_ids(i),i=1,naxes)
end subroutine define_ice_axes
