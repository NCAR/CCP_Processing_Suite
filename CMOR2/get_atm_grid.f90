!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_atm_grid
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
  double precision,allocatable,dimension(:)::slon,slat
  character(len=256)::gridfile
  integer::gridid,i,j,k,n,length
  logical::exists
  !
  select case (exp(exp_found)%model_id)
  case ('CESM1-WACCM')
     gridfile = 'atm_grid_f19_L66.nc'
  case ('CESM1-WACCMSD')
     gridfile = 'atm_grid_f19_L88.nc'
  case ('CESM1-CAM4Chem' )
     gridfile = 'atm_grid_cam4_f19_26L.nc'
  case ('CESM1-CAM4ChemSD' )
     gridfile = 'atm_grid_cam4_f19_56L.nc'
  case ('CESM1-BASE' )
     gridfile = 'atm_grid_cam4_f19_56L.nc'
  case ('CESM1-EASALL' )
     gridfile = 'atm_grid_cam4_f19_56L.nc'
  case ('CESM1-EURALL' )
     gridfile = 'atm_grid_cam4_f19_56L.nc'
  case ('CESM1-GLOALL' )
     gridfile = 'atm_grid_cam4_f19_56L.nc'
  case ('CESM1-MDEALL' )
     gridfile = 'atm_grid_cam4_f19_56L.nc'
  case ('CESM1-NAMALL' )
     gridfile = 'atm_grid_cam4_f19_56L.nc'
  case ('CESM1-RBUALL' )
     gridfile = 'atm_grid_cam4_f19_56L.nc'
  case ('CESM1-SASALL' )
     gridfile = 'atm_grid_cam4_f19_56L.nc'
  case ('CESM1-CAM5')
     gridfile = 'atm_grid_cam5_f09.nc'
  case ('CESM1-CAM5.1-FV2')
     gridfile = 'atm_grid_cam5_f19.nc'
  case ('CCSM4')
     gridfile = 'atm_grid_cam4_f09.nc'
     if (exp(exp_found)%expt_id == 'AEROCOM-A2-CTRL') then
        gridfile = 'atm_grid_cam4_f19_56L.nc'
     endif
  case default
     gridfile = 'atm_grid_cam4_f09.nc'
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
     select case (dim_info(n)%name)
     case('lat')
        nlats = dim_info(n)%length
     case('lon') 
        nlons = dim_info(n)%length
     case('lev') 
        nlevs = dim_info(n)%length
     case('ilev') 
        nilevs = dim_info(n)%length
     case('plev31') 
        nplev31 = dim_info(n)%length
     case('plev23') 
        nplev23 = dim_info(n)%length
     case('plev17') 
        nplev17 = dim_info(n)%length
     case('plev8') 
        nplev8 = dim_info(n)%length
     case('plev7') 
        nplev7 = dim_info(n)%length
     case('plev3') 
        nplev3 = dim_info(n)%length
     case('cosp_tau') 
        ncosp_tau  = dim_info(n)%length
     case('cosp_prs') 
        ncosp_prs  = dim_info(n)%length
     case('cosp_ht') 
        ncosp_ht   = dim_info(n)%length
     case('cosp_dbze') 
        ncosp_dbze = dim_info(n)%length
     case('cosp_sza') 
        ncosp_sza = dim_info(n)%length
     case('ncol') 
        nsites = dim_info(n)%length
     end select
  enddo
  allocate(atm_lons(nlons),atm_lats(nlats),slon(nlons),slat(nlats))
  allocate(atm_lons_bnds(2,nlons),atm_lats_bnds(2,nlats))
  allocate(landfrac(nlons,nlats),phis(nlons,nlats),gaussian_wts(nlats))
  allocate(area(nlons,nlats))
  allocate(atm_levs(nlevs),atm_levs_bnds(nlevs+1),atm_sites(nsites))
  allocate(atm_ilevs(nilevs),atm_ilevs_bnds(nilevs+1))
  allocate(atm_plev23(nplev23),atm_plev17(nplev17),atm_plev8(nplev8),atm_plev7(nplev7),atm_plev3(nplev3),atm_plev31(nplev31))
  allocate(a_coeff(nlevs),b_coeff(nlevs),a_coeff_bnds(nlevs+1),b_coeff_bnds(nlevs+1))
  allocate(ah_coeff(nilevs),bh_coeff(nilevs))
  allocate(cosp_tau(ncosp_tau),cosp_tau_bnds(2,ncosp_tau))
  allocate(cosp_sza(ncosp_sza))
  allocate(cosp_prs(ncosp_prs),cosp_prs_bnds(2,ncosp_prs))
  allocate(cosp_ht(ncosp_ht),cosp_ht_bnds(2,ncosp_ht))
  allocate(cosp_dbze(ncosp_dbze),cosp_dbze_bnds(2,ncosp_dbze))
  !
  do i = 1,nsites
     atm_sites(i) = i
  enddo
  !
  call get_vars(gridid)
  call read_var(gridid,'lon'   ,atm_lons)
  call read_var(gridid,'lat'   ,atm_lats)
  call read_var(gridid,'plev31',atm_plev31)
  call read_var(gridid,'plev23',atm_plev23)
  call read_var(gridid,'plev17',atm_plev17)
  call read_var(gridid,'plev8' ,atm_plev8)
  call read_var(gridid,'plev7' ,atm_plev7)
  call read_var(gridid,'plev3' ,atm_plev3)
  call read_var(gridid,'lev'   ,atm_levs)
  call read_var(gridid,'ilev'  ,atm_levs_bnds)
  call read_var(gridid,'ilev'  ,atm_ilevs)
  call read_var(gridid,'hyam'  ,a_coeff)
  call read_var(gridid,'hyai'  ,a_coeff_bnds)
  call read_var(gridid,'hybm'  ,b_coeff)
  call read_var(gridid,'hybi'  ,b_coeff_bnds)
  call read_var(gridid,'hyai'  ,ah_coeff)
  call read_var(gridid,'hybi'  ,bh_coeff)
  call read_var(gridid,'P0'    ,p0Pa)
  call read_var(gridid,'LANDFRAC',landfrac)
  call read_var(gridid,'PHIS'  ,phis)
  call read_var(gridid,'AREA'  ,area)
  call read_var(gridid,'gw'    ,gaussian_wts)
  call read_var(gridid,'cosp_tau',cosp_tau)
  call read_var(gridid,'cosp_tau_bnds',cosp_tau_bnds)
  call read_var(gridid,'cosp_sza',cosp_sza)
  call read_var(gridid,'cosp_prs',cosp_prs)
  call read_var(gridid,'cosp_prs_bnds',cosp_prs_bnds)
  cosp_prs      = cosp_prs      * 100
  cosp_prs_bnds = cosp_prs_bnds * 100
  call read_var(gridid,'cosp_ht',cosp_ht)
  call read_var(gridid,'cosp_ht_bnds',cosp_ht_bnds)
  call read_var(gridid,'cosp_dbze',cosp_dbze)
  call read_var(gridid,'cosp_dbze_bnds',cosp_dbze_bnds)
  !
  ! Replace given sza5 values ( 0 15 30 45 60 ) with correct values ( 0 20 40 60 80 )
  !
  cosp_sza = (/0.,20.,40.,60.,80./)
  !
  ! Create atm_ilev_bnds from regular levs
  !
  atm_ilevs_bnds(1) = atm_levs(1)-((atm_levs(2)-atm_levs(1))*0.5)
  do k = 2,nilevs
     atm_ilevs_bnds(k) = atm_levs(k-1)
  enddo
  atm_ilevs_bnds(nilevs+1) = atm_levs(nlevs)+((atm_levs(nlevs)-atm_levs(nlevs-1))*0.5)
!  do k = 1,nlevs
!     write(*,*) atm_levs_bnds(k),atm_levs(k),atm_levs_bnds(k)
!  enddo
  !
  ! Convert Pa values to mb (for vertint)
  !
  p0mb = p0Pa * 0.01
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
