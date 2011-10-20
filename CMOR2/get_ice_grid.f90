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
  !
  call open_cdf(gridid,'ice_grid_gx1.nc',.true.)
  !
  ! Read time-invariant dimensions and variables from 'ice_grid_gx1.nc'
  !
  call get_dims(gridid)
  !
  do n=1,dim_counter
     length = len_trim(dim_info(n)%name)
     if(dim_info(n)%name(:length).eq.'nj') then
        nlats = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'ni') then
        nlons = dim_info(n)%length
     endif
  enddo
  allocate(ice_lons(nlons,nlats),ice_lats(nlons,nlats))
  !
  call get_vars(gridid)
  call read_var(gridid,'TLON',ice_lons)
  call read_var(gridid,'TLAT',ice_lats)
   !
  ! Transfer bounds for lons and lats
  !
  allocate(ice_lons_bnds(4,nlons,nlats),ice_lats_bnds(4,nlons,nlats))
  call read_var(gridid,'lont_bounds',ice_lons_bnds)
  call read_var(gridid,'latt_bounds',ice_lats_bnds)
!!$  bnds_lat(1,1)     = -90.
!!$  bnds_lat(2,nlats) =  90.
!!$  do j = 1,nlats-1
!!$     bnds_lat(2,j) = slat(j)
!!$  end do
!!$  do j = 2,nlats
!!$     bnds_lat(1,j) =  bnds_lat(2,j-1)
!!$  end do
!!$  !
!!$  bnds_lon(1,    1) = slon(1)
!!$  bnds_lon(2,nlons) = alons(nlons) + ((slon(nlons)-slon(nlons-1))/2.)
!!$  do i = 1,nlons-1
!!$     bnds_lon(2,i) = slon(i+1)
!!$  end do
!!$  do i = 2,nlons
!!$     bnds_lon(1,i) =  bnds_lon(2,i-1)
!!$  end do
!!$  !
  call close_cdf(gridid)
  write(*,'(''ICE grid loaded'')')
!!$  !
!!$  do j = 1,nlats
!!$     write(*,*) bnds_lat(1,j),alats(j),bnds_lat(2,j)
!!$  enddo
!!$  do i = 1,nlons
!!$     write(*,*) bnds_lon(1,i),alons(i),bnds_lon(2,i)
!!$  enddo
end subroutine get_ice_grid
