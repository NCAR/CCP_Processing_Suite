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
             interval='30 days')
        idim = idim + 1
     end select
  enddo
  write(*,*) 'lnd axes defined, axis_ids: ',(axis_ids(i),i=1,naxes)
end subroutine define_lnd_axes
