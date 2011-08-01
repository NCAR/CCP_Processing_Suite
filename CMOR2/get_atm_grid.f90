!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_atm_grid
  !
  ! Read coordinate information from model into arrays that will be passed
  !   to CMOR.
  ! Read latitude, longitude, and pressure coordinate values into
  !   alats, alons, and plevs, respectively.  Also generate latitude and
  !   longitude bounds, and store in bnds_lat and bnds_lon, respectively.
  !
  use counters_netcdf_jfl
  use interfaces_netcdf_jfl
  use definitions_netcdf_jfl
  use atm_grid_info
  !
  implicit none
  !
  integer::gridid,i,j,n,length
  !
  call open_cdf(gridid,'atm_grid_f09.nc',.true.)
  !
  ! Read time-invariant dimensions and variables from 'atm_grid_f09.nc'
  !
  call get_dims(gridid)
  !
  do n=1,dim_counter
     length = len_trim(dim_info(n)%name)
     if(dim_info(n)%name(:length).eq.'lat') then
        nlats = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'lon') then
        nlons = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'lev') then
        nlevs = dim_info(n)%length
     endif
  enddo
  ALLOCATE(alons(nlons),alats(nlats),slon(nlons),slat(nlats-1))
  ALLOCATE(zlevs(nlevs),zlev_bnds(nlevs+1))
  ALLOCATE(a_coeff(nlevs),b_coeff(nlevs),a_coeff_bnds(nlevs+1),b_coeff_bnds(nlevs+1))
  ALLOCATE(bnds_lon(2,nlons),bnds_lat(2,nlats),plevs(17))
  !
  call get_vars(gridid)
  call read_var(gridid,'lon',alons)
  call read_var(gridid,'lat',alats)
  call read_var(gridid,'lev',zlevs)
  !
  ! Convert zlevs from mb to hPa
  !
  zlevs = zlevs * 100
  plevs = (/100000., 92500., 85000., 70000.,60000., 50000., 40000., 30000., 25000., &
       20000., 15000., 10000.,  7000., 5000.,  3000.,  2000.,  1000./)
  !
  ! Transfer bounds for lons and lats
  !
  call read_var(gridid,'slon',slon)
  call read_var(gridid,'slat',slat)
  bnds_lat(1,1)     = -90.
  bnds_lat(2,nlats) =  90.
  do j = 1,nlats-1
     bnds_lat(2,j) = slat(j)
  end do
  do j = 2,nlats
     bnds_lat(1,j) =  bnds_lat(2,j-1)
  end do
  !
  bnds_lon(1,    1) = slon(1)
  bnds_lon(2,nlons) = alons(nlons) + ((slon(nlons)-slon(nlons-1))/2.)
  do i = 1,nlons-1
     bnds_lon(2,i) = slon(i+1)
  end do
  do i = 2,nlons
     bnds_lon(1,i) =  bnds_lon(2,i-1)
  end do
  !
  call close_cdf(gridid)
  write(*,'(''ATM grid loaded'')')
  !
!!$  do j = 1,nlats
!!$     write(*,*) bnds_lat(1,j),alats(j),bnds_lat(2,j)
!!$  enddo
!!$  do i = 1,nlons
!!$     write(*,*) bnds_lon(1,i),alons(i),bnds_lon(2,i)
!!$  enddo
end subroutine get_atm_grid
