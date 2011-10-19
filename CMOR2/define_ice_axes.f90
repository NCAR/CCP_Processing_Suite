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
  character(len=256),intent(in)::dimensions
  !
  integer::i,j,idim,status,table_id,grid_id
  integer,dimension(0:10)::idxb
  !
  character(len=256),dimension(10)::dimnames,dimunits
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
  dimids = 0 ; idim = 1
  !
  table_id = cmor_load_table('Tables/CMIP5_grids')
  call cmor_set_table(table_id)
  !
  do i = 1,naxes
     select case(dimnames(i))
     case ('latitude')
        dimids(idim) = cmor_axis(          &
             table_entry='grid_latitude',  &
             units=dimunits(i))
!             coord_vals=ice_lats)
        idim = idim + 1
     case ('longitude')
        dimids(idim) = cmor_axis(          &
             table_entry='grid_longitude', &
             units=dimunits(i))
!             coord_vals=ice_lons)
        idim = idim + 1
     case ('time')
        dimids(idim) = cmor_axis(          &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=ntimes,                &
             interval='30 days')
        idim = idim + 1
     end select
  enddo
  write(*,*) 'CMOR axes defined, dimids: ',(dimids(i),i=1,naxes)
  grid_id  = cmor_grid(                  &
       axis_ids=(/dimids(1),dimids(2)/), &
       latitude=ice_lats,                &
       longitude=ice_lons,               &
       latitude_vertices=ice_lats_bnds,  &
       longitude_vertices=ice_lons_bnds, &
       4)
end subroutine define_ice_axes
