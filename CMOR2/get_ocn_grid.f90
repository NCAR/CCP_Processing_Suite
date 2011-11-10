!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_ocn_grid
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
  !
  call open_cdf(gridid,'ocn_grid_gx1.nc',.true.)
  !
  ! Read time-invariant dimensions and variables from 'ocn_grid_gx1.nc'
  !
  call get_dims(gridid)
  !
  do n=1,dim_counter
     length = len_trim(dim_info(n)%name)
     if(dim_info(n)%name(:length).eq.'nlat') then
        nlats = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'nlon') then
        nlons = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'z_t') then
        nlevs = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'lat_aux_grid') then ! heat, salt and mass transports on this grid
        nlats_trans = dim_info(n)%length
     endif
  enddo
  allocate(ocn_lons(nlons,nlats),ocn_lats(nlons,nlats),ocn_levs(nlevs),ocn_levs_bnds(nlevs+1))
  allocate(ocn_trans_lats(nlats_trans))
  !
  call get_vars(gridid)
  call read_var(gridid,'TLONG',ocn_lons)
  call read_var(gridid,'TLAT',ocn_lats)
  call read_var(gridid,'z_t',ocn_levs)
  call read_var(gridid,'z_w',ocn_levs_bnds(1:nlevs))
  call read_var(gridid,'lat_aux_grid',ocn_trans_lats)
  write(*,*) 'Got lat_aux_grid: ',ocn_trans_lats(1),ocn_trans_lats(nlats_trans)
  !
  !  Fill last ocn_bnds level
  !
  ocn_levs_bnds(nlevs+1) = ocn_levs(nlevs)+(0.5*ocn_levs_bnds(nlevs))
  !
  ! Convert ocn_levs and ocn_levs_bnds from cm to m
  !
  ocn_levs      = ocn_levs      / 100.
  ocn_levs_bnds = ocn_levs_bnds / 100.
  !
  ! Transfer bounds for lons and lats
  !
  allocate(ocn_lons_bnds(4,nlons,nlats),ocn_lats_bnds(4,nlons,nlats))
  call read_var(gridid,'lont_bounds',ocn_lons_bnds)
  call read_var(gridid,'latt_bounds',ocn_lats_bnds)
  call close_cdf(gridid)
  write(*,'(''OCN grid loaded'')')
end subroutine get_ocn_grid
