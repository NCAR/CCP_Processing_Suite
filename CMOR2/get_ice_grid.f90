!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_ice_grid
  !
  ! Read coordinate information from model into arrays that will be passed to CMOR.
  !
  use counters_netcdf_jfl
  use interfaces_netcdf_jfl
  use definitions_netcdf_jfl
  use grid_info
  !
  implicit none
  !
  integer::gridid,i,j,n,length
  logical::exists
  !
  inquire(file='ice_grid_gx1.nc',exist=exists)
  if (.not.(exists)) then
     write(*,*) 'Cannot find ice_grid_gx1.nc - STOPPING.'
     stop
  endif
  !
  call open_cdf(gridid,'ice_grid_gx1.nc',.true.)
  !
  ! Read time-invariant dimensions and variables from 'ice_grid_gx1.nc'
  !
  call get_dims(gridid)
  !
  do n=1,dim_counter
     length = len_trim(dim_info(n)%name)
     if (dim_info(n)%name(:length).eq.'nj') nlats = dim_info(n)%length
     if (dim_info(n)%name(:length).eq.'ni') nlons = dim_info(n)%length
  enddo
  allocate(ice_t_lons(nlons,nlats),ice_t_lats(nlons,nlats))
  allocate(ice_u_lons(nlons,nlats),ice_u_lats(nlons,nlats))
  !
  call get_vars(gridid)
  call read_var(gridid,'TLON',ice_t_lons)
  call read_var(gridid,'TLAT',ice_t_lats)
  call read_var(gridid,'ULON',ice_u_lons)
  call read_var(gridid,'ULAT',ice_u_lats)
  !
  ! Transfer bounds for lons and lats
  !
  allocate(ice_t_lons_bnds(4,nlons,nlats),ice_t_lats_bnds(4,nlons,nlats))
  allocate(ice_u_lons_bnds(4,nlons,nlats),ice_u_lats_bnds(4,nlons,nlats))
  call read_var(gridid,'lont_bounds',ice_t_lons_bnds)
  call read_var(gridid,'latt_bounds',ice_t_lats_bnds)
  call read_var(gridid,'lonu_bounds',ice_u_lons_bnds)
  call read_var(gridid,'latu_bounds',ice_u_lats_bnds)
  !
  call close_cdf(gridid)
  write(*,'(''ICE grid loaded'')')
end subroutine get_ice_grid
