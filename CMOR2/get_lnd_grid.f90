!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_lnd_grid
  !
  ! Read coordinate information from model into arrays that will be passed to CMOR.
  !
  use counters_netcdf_jfl
  use interfaces_netcdf_jfl
  use definitions_netcdf_jfl
  use grid_info
  use exp_info
  use mycmor_info
  !
  implicit none
  !
  character(len=256)::gridfile
  integer::gridid,i,j,k,n,length
  logical::exists
  !
  select case (exp(exp_found)%model_id)
  case ('CESM1-WACCM')
     gridfile = 'lnd_grid_f19.nc'
  case default
     gridfile = 'lnd_grid_f09.nc'
  end select
  !
  inquire(file=trim(gridfile),exist=exists)
  if (.not.(exists)) then
     write(*,*) 'Cannot find ',trim(gridfile),' - STOPPING.'
     stop
  endif
  !
  call open_cdf(gridid,trim(gridfile),.true.)
  !
  ! Read time-invariant dimensions and variables from gridfile
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
  allocate(lnd_lons(nlons),lnd_lats(nlats),lnd_levs(nlevs),lnd_levs_bnds(nlevs+1))
  allocate(lnd_lons_bnds(2,nlons),lnd_lats_bnds(2,nlats))
  allocate(lnd_zsoi(nlons,nlats,nlevs),lnd_dzsoi(nlons,nlats,nlevs))
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
  ! Get soil layer thickness (1d)
  !
  call read_var(gridid,'levgrnd',lnd_levs)
  !
  ! Create bounds for levs
  !
  lnd_levs_bnds(1) = 0.
  do k = 1,nlevs-1
     lnd_levs_bnds(k+1) = 0.5*(lnd_levs(k)+lnd_levs(k+1))
  enddo
  lnd_levs_bnds(nlevs+1) = lnd_levs(nlevs)+(0.5*lnd_levs_bnds(nlevs-1))
  !
!  write(*,'(''LND levs: '',15(''Z '',f9.4))') (lnd_levs(k),k=1,nlevs)
!  write(*,'(''LND bnds: '',16(''B '',f9.4))') (lnd_levs_bnds(k),k=1,nlevs+1)
  !
  ! Get soil depths and soil layer thickness (3d)
  !
  call read_var(gridid,'ZSOI',lnd_zsoi)
  call read_var(gridid,'DZSOI',lnd_dzsoi)
  call close_cdf(gridid)
  write(*,'(''LND grid loaded'')')
!!$  !
!!$  do j = 1,nlats
!!$     write(*,*) lnd_lats_bnds(1,j),lnd_lats(j),lnd_lats_bnds(2,j)
!!$  enddo
!!$  do i = 1,nlons
!!$     write(*,*) lnd_lons_bnds(1,i),lnd_lons(i),lnd_lons_bnds(2,i)
!!$  enddo
end subroutine get_lnd_grid
