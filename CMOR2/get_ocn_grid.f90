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
     elseif(dim_info(n)%name(:length).eq.'moc_z') then ! heat, salt and mass transports on this grid
        nmoc_z = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'transport_reg') then ! heat, salt and mass transports on this grid
        ntrans_reg = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'moc_comp') then ! heat, salt and mass transports on this grid
        nmoc_comp = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'transport_comp') then ! heat, salt and mass transports on this grid
        ntrans_comp = dim_info(n)%length
     endif
  enddo
  !
  allocate(ocn_t_lons(nlons,nlats),ocn_t_lats(nlons,nlats),kmt(nlons,nlats))
  allocate(ocn_t_lons_bnds(4,nlons,nlats),ocn_t_lats_bnds(4,nlons,nlats))
  allocate(ocn_t_levs(nlevs),ocn_t_levs_bnds(2,nlevs))
  !
  allocate(ocn_u_lons(nlons,nlats),ocn_u_lats(nlons,nlats),kmt(nlons,nlats))
  allocate(ocn_u_lons_bnds(4,nlons,nlats),ocn_u_lats_bnds(4,nlons,nlats))
  !
  allocate(ocn_trans_lats(nlats_trans),ocn_trans_lats_bnds(2,nlats_trans))
  allocate(ocn_trans_levs(nmoc_z),ocn_trans_levs_bnds(2,nmoc_z))
  !
  call get_vars(gridid)
  call read_var(gridid,'TLONG',ocn_t_lons)
  call read_var(gridid,'TLAT',ocn_t_lats)
  call read_var(gridid,'ULONG',ocn_u_lons)
  call read_var(gridid,'ULAT',ocn_u_lats)
  call read_var(gridid,'KMT',kmt)
  call read_var(gridid,'z_t',ocn_t_levs)
  call read_var(gridid,'lat_aux_grid',ocn_trans_lats)
  call read_var(gridid,'moc_z',ocn_trans_levs)
  !
  ! Create bounds for ocn_trans_lats, ocn_trans_levs, and ocn_levs
  !
  ocn_trans_lats_bnds(1,1)           = ocn_trans_lats(1)
  ocn_trans_lats_bnds(2,nlats_trans) = ocn_trans_lats(nlats_trans)
  do j = 1,nlats_trans-1
     ocn_trans_lats_bnds(2,j) = ocn_trans_lats(j)+((ocn_trans_lats(j+1)-ocn_trans_lats(j))/2.)
  end do
  do j = 2,nlats_trans
     ocn_trans_lats_bnds(1,j) =  ocn_trans_lats_bnds(2,j-1)
  end do
  !
  ocn_trans_levs_bnds(1,    1) = 0.
  ocn_trans_levs_bnds(2,nmoc_z) = ocn_trans_levs(nmoc_z) + ((ocn_trans_levs(nmoc_z)-ocn_trans_levs(nmoc_z-1))/2.)
  do i = 1,nmoc_z-1
     ocn_trans_levs_bnds(2,i) = ocn_trans_levs(i)+((ocn_trans_levs(i+1)-ocn_trans_levs(i))/2.)
  end do
  do i = 2,nmoc_z
     ocn_trans_levs_bnds(1,i) =  ocn_trans_levs_bnds(2,i-1)
  end do
  !
  ocn_t_levs_bnds(1,    1) = 0.
  ocn_t_levs_bnds(2,nlevs) = ocn_t_levs(nlevs) + ((ocn_t_levs(nlevs)-ocn_t_levs(nlevs-1))/2.)
  do i = 1,nlevs-1
     ocn_t_levs_bnds(2,i) = ocn_t_levs(i)+((ocn_t_levs(i+1)-ocn_t_levs(i))/2.)
  end do
  do i = 2,nlevs
     ocn_t_levs_bnds(1,i) =  ocn_t_levs_bnds(2,i-1)
  end do
  !
  ! Convert ocn_t_levs and ocn_t_levs_bnds from cm to m
  !
  ocn_t_levs            = ocn_t_levs            / 100.
  ocn_t_levs_bnds       = ocn_t_levs_bnds       / 100.
  ocn_trans_levs      = ocn_trans_levs      / 100.
  ocn_trans_levs_bnds = ocn_trans_levs_bnds / 100.
  !
  ! Transfer bounds for lons and lats
  !
  call read_var(gridid,'lont_bounds',ocn_t_lons_bnds)
  call read_var(gridid,'latt_bounds',ocn_t_lats_bnds)
  call read_var(gridid,'lonu_bounds',ocn_u_lons_bnds)
  call read_var(gridid,'latu_bounds',ocn_u_lats_bnds)
  call close_cdf(gridid)
  write(*,'(''OCN grid loaded'')')
end subroutine get_ocn_grid
