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
  axis_ids = 0 ; idim = 1
  !
  table_id = cmor_load_table('Tables/CMIP5_grids')
  call cmor_set_table(table_id)
  !
  do i = 1,naxes
     select case(dimnames(i))
     case ('latitude')
        axis_ids(idim) = cmor_axis(       &
             table_entry='j',             &
             length=size(ice_lats),       &
             units='1')
        idim = idim + 1
     case ('longitude')
        axis_ids(idim) = cmor_axis(       &
             table_entry='i',             &
             length=size(ice_lons),       &
             units='1')
        idim = idim + 1
     case ('time')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=ntimes,                &
             interval='30 days')
        idim = idim + 1
     end select
  enddo
  write(*,*) 'CMOR axes defined, axis_ids: ',(axis_ids(i),i=1,naxes)
  grid_id = cmor_grid_2d_r(                  &
       axis_ids=(/axis_ids(1),axis_ids(2)/), &
       latitude=ice_lats,                    &
       longitude=ice_lons,                   &
       latitude_vertices=ice_lats_bnds,      &
       longitude_vertices=ice_lons_bnds)
end subroutine define_ice_axes
