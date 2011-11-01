subroutine define_lnd_axes(dimensions)
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
  integer::i,j,idim
  integer,dimension(0:10)::idxb
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
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(lnd_lons),        &
             units=dimunits(i),            &
             coord_vals=lnd_lons,          &
             cell_bounds=lnd_lons_bnds)
        idim = idim + 1
     case ('latitude')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=SIZE(lnd_lats),        &
             coord_vals=lnd_lats,          &
             cell_bounds=lnd_lats_bnds)
        idim = idim + 1
     case ('sdepth')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(lnd_levs),        &
             units=dimunits(i),            &
             coord_vals=lnd_levs,          &
             cell_bounds=lnd_levs_bnds)
        idim = idim + 1
        write(*,*) 'lnd_levs: ',lnd_levs
     case ('time')
        axis_ids(idim) = cmor_axis(          &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
!             length=ntimes,                &
             interval='30 days')
        idim = idim + 1
     end select
  enddo
  write(*,*) 'lnd axes defined, axis_ids: ',(axis_ids(i),i=1,naxes)
end subroutine define_lnd_axes
