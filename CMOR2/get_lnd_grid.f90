!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_lnd_grid
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
  integer::gridid,i,j,k,n,length
  !
  call open_cdf(gridid,'lnd_grid_f09.nc',.true.)
  !
  ! Read time-invariant dimensions and variables from 'lnd_grid_f09.nc'
  !
  call get_dims(gridid)
  !
  do n=1,dim_counter
     length = len_trim(dim_info(n)%name)
     if(dim_info(n)%name(:length).eq.'lat') then
        nlats = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'lon') then
        nlons = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'levgrnd') then
        nlevs = dim_info(n)%length
     endif
  enddo
  ALLOCATE(lnd_lons(nlons),lnd_lats(nlats),lnd_levs(nlevs))
  ALLOCATE(lnd_lons_bnds(2,nlons),lnd_lats_bnds(2,nlats))
  !
  call get_vars(gridid)
  call read_var(gridid,'lon',lnd_lons)
  call read_var(gridid,'lat',lnd_lats)
  !
  ! Create bounds for lons and lats
  !
  lnd_lats_bnds(1,1)     = -90.
  lnd_lats_bnds(2,nlats) =  90.
  do j = 1,nlats-1
     lnd_lats_bnds(2,j) = lnd_lats(j)+((lnd_lats(j+1)-lnd_lats(j))/2.)
  end do
  do j = 2,nlats
     lnd_lats_bnds(1,j) =  lnd_lats_bnds(2,j-1)
  end do
  !
  lnd_lons_bnds(1,    1) = lnd_lons(1)
  lnd_lons_bnds(2,nlons) = lnd_lons(nlons) + ((lnd_lons(nlons)-lnd_lons(nlons-1))/2.)
  do i = 1,nlons-1
     lnd_lons_bnds(2,i) = lnd_lons(i)+((lnd_lons(i+1)-lnd_lons(i))/2.)
  end do
  do i = 2,nlons
     lnd_lons_bnds(1,i) =  lnd_lons_bnds(2,i-1)
  end do
  !
  ! Get soil depths
  !
  call read_var(gridid,'levgrnd',lnd_levs)
  write(*,'(''LND levs: '',15f9.4)') (lnd_levs(k),k=1,nlevs)
  call close_cdf(gridid)
  write(*,'(''LND grid loaded'')')
  !
!!$  do j = 1,nlats
!!$     write(*,*) lnd_lats_bnds(1,j),lnd_lats(j),lnd_lats_bnds(2,j)
!!$  enddo
!!$  do i = 1,nlons
!!$     write(*,*) lnd_lons_bnds(1,i),lnd_lons(i),lnd_lons_bnds(2,i)
!!$  enddo
end subroutine get_lnd_grid
