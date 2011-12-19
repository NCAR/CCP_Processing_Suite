!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_atm_grid
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
  double precision,allocatable,dimension(:)::slon,slat
  integer::gridid,i,j,k,n,length
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
     elseif(dim_info(n)%name(:length).eq.'ilev') then
        nilevs = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'plevs') then
        nplevs = dim_info(n)%length
     elseif(dim_info(n)%name(:length).eq.'plev8') then
        nplev8 = dim_info(n)%length
     endif
  enddo
  allocate(atm_lons(nlons),atm_lats(nlats),slon(nlons),slat(nlats))
  allocate(atm_lons_bnds(2,nlons),atm_lats_bnds(2,nlats))
  allocate(atm_levs(nlevs),atm_levs_bnds(nlevs+1))
  allocate(atm_ilevs(nilevs),atm_ilevs_bnds(nilevs+1))
  allocate(atm_plevs(nplevs),atm_plev8(nplev8))
  allocate(a_coeff(nlevs),b_coeff(nlevs),a_coeff_bnds(nlevs+1),b_coeff_bnds(nlevs+1))
  !
  call get_vars(gridid)
  call read_var(gridid,'lon'  ,atm_lons)
  call read_var(gridid,'lat'  ,atm_lats)
  call read_var(gridid,'plevs',atm_plevs)
  call read_var(gridid,'plev8',atm_plev8)
  call read_var(gridid,'lev'  ,atm_levs)
  call read_var(gridid,'ilev' ,atm_levs_bnds)
  call read_var(gridid,'ilev' ,atm_ilevs)
  call read_var(gridid,'hyam' ,a_coeff)
  call read_var(gridid,'hyai' ,a_coeff_bnds)
  call read_var(gridid,'hybm' ,b_coeff)
  call read_var(gridid,'hybi' ,b_coeff_bnds)
  call read_var(gridid,'P0'   ,p0)
  !
  !
  ! Create atm_ilev_bnds from regular levs
  !
  atm_ilevs_bnds(1) = atm_levs(1)-((atm_levs(2)-atm_levs(1))*0.5)
  do k = 2,nilevs
     atm_ilevs_bnds(k) = atm_levs(k-1)
  enddo
  atm_ilevs_bnds(nilevs+1) = atm_levs(nlevs)+((atm_levs(nlevs)-atm_levs(nlevs-1))*0.5)
  !
  ! Convert Pa values to mb (for vertint)
  !
  p0 = p0 * 0.01
  !
  ! Transfer bounds for lons and lats
  !
  call read_var(gridid,'slon',slon)
  call read_var(gridid,'slat',slat)
  atm_lats_bnds(1,1)     = -90.
  atm_lats_bnds(2,nlats) =  90.
  do j = 1,nlats-1
     atm_lats_bnds(2,j) = slat(j)
  end do
  do j = 2,nlats
     atm_lats_bnds(1,j) =  atm_lats_bnds(2,j-1)
  end do
  !
  atm_lons_bnds(1,    1) = slon(1)
  atm_lons_bnds(2,nlons) = atm_lons(nlons) + ((slon(nlons)-slon(nlons-1))/2.)
  do i = 1,nlons-1
     atm_lons_bnds(2,i) = slon(i+1)
  end do
  do i = 2,nlons
     atm_lons_bnds(1,i) =  atm_lons_bnds(2,i-1)
  end do
  !
  call close_cdf(gridid)
  write(*,'(''ATM grid loaded'')')
  !
!!$  do j = 1,nlats
!!$     write(*,*) atm_lats_bnds(1,j),atm_lats(j),atm_lats_bnds(2,j)
!!$  enddo
!!$  do i = 1,nlons
!!$     write(*,*) atm_lons_bnds(1,i),atm_lons(i),atm_lons_bnds(2,i)
!!$  enddo
end subroutine get_atm_grid
