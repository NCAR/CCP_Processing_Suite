!!$program DRIVE_file_output_times
!!$  use output_times_info
!!$  implicit none
!!$  character(len=256)::model_id,realm,dims
!!$  real::size_bytes,max_size,tot_size
!!$  integer,dimension(21)::itimes
!!$  integer::ntimes,i,nmax
!!$  !
!!$  itimes = (/60,120,828,900,1140,1152,1460,1872,2388,2400,3612,4824,6012,11315,12012,14600,29200,32120,34675,35040,56940/)
!!$  !
!!$  max_size = (2.**32.)-1
!!$  write(*,*) 'ms: ',max_size
!!$  !
!!$  model_id = 'CCSM4'
!!$!  realm    = 'atmos'
!!$  realm    = 'atmos'
!!$  dims     = 'longitude latitude alevel time'
!!$  call calc_size(size_bytes,model_id,realm,dims)
!!$  do i = 1,21
!!$     ntimes = itimes(i)
!!$     call file_output_times(ntimes)
!!$     tot_size  = size_bytes*ntimes
!!$     nmax = (tot_size/max_size)+1
!!$     write(*,'(''GB: '',f10.3,'' nm: '',i10,'' nt: '',i10,'' t0: '',i10,'' t1: '',i10)') tot_size/1.e9,nmax,ntimes,tidx1(1),tidx2(1)
!!$  enddo
!!$end program DRIVE_file_output_times
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine file_output_times(ntimes)
  use output_times_info
  implicit none
  integer,intent(in)::ntimes
  !
  select case (ntimes)
  case ( 1152 )           ! RCP from 2005-2100, skip 2005
     nchunks(1) = 1
     tidx1(1:nchunks(1)) = (/    13/)
     tidx2(1:nchunks(1)) = (/ntimes/)
  case default
     nchunks(1) = 1
     tidx1(1:nchunks(1)) = (/     1/)
     tidx2(1:nchunks(1)) = (/ntimes/)
  end select
end subroutine file_output_times
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_size(size_bytes,model_id,realm,dims)
  implicit none
  real,intent(out)::size_bytes
  character(len=256),intent(in)::model_id,realm,dims
  integer::nxa,nya,nza,nzai,nxl,nyl,nzl,nxo,nyo,nzo,nyto
  integer::alt40,dbze,scatratio,sites,sza5,oline,basin
  integer::plev17,plev8,plev7,plev3,plev23
  !
  alt40     =  40
  dbze      =  15
  scatratio =  15
  sza5      =   5
  oline     =  15
  basin     =   3
  sites     = 119
  plev23    =  23
  plev17    =  17
  plev8     =   8
  plev7     =   7
  plev3     =   3
  !
  size_bytes = 0.
  !
  select case (trim(model_id))
  case ('CESM1-WACCM')
     nxa = 144 ; nya =  96 ; nza = 66 ; nzai =  67 ; nzl = 15 
     nxo = 384 ; nyo = 320 ; nzo = 60 ; nyto = 395
  case ('CESM1-CAM5')
     nxa = 288 ; nya = 192 ; nza = 30 ; nzai =  31 ; nzl = 15 
     nxo = 384 ; nyo = 320 ; nzo = 60 ; nyto = 395
  case default
     nxa = 288 ; nya = 192 ; nza = 26 ; nzai =  27 ; nzl = 15 
     nxo = 384 ; nyo = 320 ; nzo = 60 ; nyto = 395
  end select
  !
  select case (trim(realm))
  case ( 'aerosol','aerosol land')
     select case (trim(dims))
     case ( 'longitude latitude alev1 time','longitude latitude time')
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude alevel time')
        size_bytes = float((nxa * nya * nza * 4) + (nxa * nya * 4))
     end select
  case ('atmos','atmos atmosChem','atmos land')
     select case (trim(dims))
     case ( 'alevel site time1' )
        size_bytes = float(sites * nza * 4)
     case ( 'alevhalf site time1' )
        size_bytes = float(sites * nzai * 4)
     case ( 'location alt40 dbze time1' )
        size_bytes = float(sites * alt40 * dbze * 4)
     case ( 'location alt40 scatratio time1' )
        size_bytes = float(sites * alt40 * scatratio * 4)
     case ( 'location alt40 time1' )
        size_bytes = float(sites * alt40 * 4)
     case ( 'location sza5 time1' )
        size_bytes = float(sites * sza5 * 4)
     case ( 'location time1','location time1 p220','location time1 p560','location time1 p840' )
        size_bytes = float(sites * 4)
     case ( 'longitude latitude alevel time','longitude latitude alevel time1' )
        size_bytes = float((nxa * nya * nza * 4) + (nxa * nya * 4))
     case ( 'longitude latitude alevhalf time','longitude latitude alevhalf time1' )
        size_bytes = float(nxa * nya * nzai * 4)
     case ( 'longitude latitude alt40 dbze time' )
        size_bytes = float(nxa * nya * alt40 * dbze * 4)
     case ( 'longitude latitude alt40 scatratio time' )
        size_bytes = float(nxa * nya * alt40 * scatratio * 4)
     case ( 'longitude latitude alt40 time' )
        size_bytes = float(nxa * nya * alt40 * 4)
     case ( 'longitude latitude plev3 time1' )
        size_bytes = float(nxa * nya * plev3 * 4)
     case ( 'longitude latitude plev7 tau time' )
        size_bytes = float(nxa * nya * plev7 * 4)
     case ( 'longitude latitude plev8 time' )
        size_bytes = float(nxa * nya * plev8 * 4)
     case ( 'longitude latitude plevs time','longitude latitude plevs time2' )
        size_bytes = float(nxa * nya * plev17 * 4)
     case ( 'longitude latitude sza5 time' )
        size_bytes = float(nxa * nya * sza5 * 4)
     case ( 'longitude latitude time','longitude latitude time1','longitude latitude time1 height10m','longitude latitude time1 height2m',&
            'longitude latitude time height10m','longitude latitude time height2m',&
            'longitude latitude time p220','longitude latitude time p500','longitude latitude time p560','longitude latitude time p700','longitude latitude time p840' )
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude' )
        size_bytes = float(nxa * nya * 4)
     case ( 'site time1','site time1 height10m','site time1 height2m' )
        size_bytes = float(sites * 4)
     case ( 'time','time2' )
        size_bytes = 4.
     end select
  case ( 'land','landIce land','land landIce' )
     select case (trim(dims))
     case ( 'longitude latitude sdepth time' )
        size_bytes = float(nxa * nya * nzl * 4)
     case ( 'longitude latitude time','longitude latitude time1' )
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude time1 sdepth1','longitude latitude time sdepth1' )
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude time typebare' )
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude time typec3pft' )
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude time typec4pft' )
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude time typepdec' )
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude time typepever' )
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude time typesdec' )
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude time typesever' )
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude vegtype time' )
        size_bytes = float(nxa * nya * 4)
     case ( 'longitude latitude' )
        size_bytes = float(nxa * nya * 4)
     end select
  case ('ocean','ocnBgchem')
     select case (trim(dims))
     case ('latitude basin time')
        size_bytes = float(nyto * basin * 4)
     case ('latitude olevel basin time')
        size_bytes = float(nyto * nzo * basin * 4)
     case ('latitude rho basin time')
        size_bytes = float(nyto * basin * 4)
     case ('longitude latitude olevel time','longitude latitude olevel time2')
        size_bytes = float(nxo * nyo * nzo * 4)
     case ('longitude latitude olevel')
        size_bytes = float(nxo * nyo * nzo * 4) + float(nxo * nyo * 4 * 11)
     case ('longitude latitude time','longitude latitude time1','longitude latitude time2',&
           'longitude latitude time depth0m','longitude latitude time depth100m','longitude latitude time olayer100m')
        size_bytes = float(nxo * nyo * 4)
     case ('longitude latitude')
        size_bytes = float(nxo * nyo * 4) + float(nxo * nyo * 4 * 11)
     case ('oline time')
        size_bytes = float(oline * 4)
     case ('time')
        size_bytes = 4.
     end select
  case ('ocean seaIce','seaIce','seaIce ocean')
     select case (trim(dims))
     case ('longitude latitude olevel time')
        size_bytes = float(nxo * nyo * nzo * 4)
     case ('longitude latitude time')
        size_bytes = float(nxo * nyo * 4)
     case ('time')
        size_bytes = 4.
     end select
  end select
end subroutine calc_size
