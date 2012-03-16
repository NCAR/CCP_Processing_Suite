!!$program DRIVE_build_filenames
!!$  use exp_info
!!$  use xwalk_info
!!$  use files_info
!!$  implicit none
!!$  integer::i,j,ixw
!!$  !
!!$  case_read = 'b40.1850.track1.1deg.006'
!!$  comp_read = 'pop.h'
!!$  xw_found  = 1
!!$  exp_found = 1
!!$  ixw       = 1
!!$  xw(ixw)%cesm_vars(1) = 'SALT'
!!$  xw(ixw)%cesm_vars(2) = 'TEMP'
!!$  xw(ixw)%ncesm_vars   = 2
!!$  exp(exp_found)%begyr =  800
!!$  exp(exp_found)%endyr = 1300
!!$  !
!!$  ncfile(:,:)(1:) = ' '
!!$  nc_nfiles(:)    = 0
!!$  !
!!$  do j = 1,xw(ixw)%ncesm_vars
!!$     do i = 1,1
!!$        call build_filenames(i,j,exp(exp_found)%begyr,exp(exp_found)%endyr)
!!$     enddo
!!$     do i = 1,nc_nfiles(j)
!!$        write(*,'(''I: '',2i8,5x,a)') i,j,trim(ncfile(i,j))
!!$     enddo
!!$  enddo
!!$  write(*,*) nc_nfiles(1:10)
!!$  do j = 1,xw(ixw)%ncesm_vars
!!$     do i = 1,nc_nfiles(j)
!!$        write(*,'(''O: '',2i8,5x,a)') i,j,trim(ncfile(i,j))
!!$     enddo
!!$  enddo
!!$end program DRIVE_build_filenames
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_filenames(case,comp,cesm_var,ivar,begyr,endyr,table)
  use files_info
  !
  implicit none
  character(len=256),intent(in)::case,comp,cesm_var,table
  integer,intent(in)::ivar,begyr,endyr
  integer::i,j,year1,year2,ndt,idt
  character(len=256)::checkname
  character(len=256),dimension(2)::dtbeg,dtend
  logical::exists
  !
!  write(*,*) 'Entering build_filenames: ',trim(case),' ',trim(comp),' ',trim(cesm_var),ivar,begyr,endyr,trim(table)
  !
  select case (table)
  case ('Tables/CMIP5_Amon','Tables/CMIP5_Lmon','Tables/CMIP5_LImon','Tables/CMIP5_Omon','Tables/CMIP5_OImon','Tables/CMIP5_aero','Tables/CMIP5_cfMon',&
        'Tables/GeoMIP_Amon','Tables/GeoMIP_Lmon','Tables/GeoMIP_LImon','Tables/GeoMIP_Omon','Tables/GeoMIP_OImon','Tables/GeoMIP_aero','Tables/GeoMIP_cfMon')
     ndt = 1 ; dtbeg(1) = '01' ; dtend(1) = '12'
  case ('Tables/CMIP5_day','Tables/CMIP5_cfDay','Tables/GeoMIP_day','Tables/GeoMIP_cfDay')
     ndt = 2
     dtbeg(1) = '0101' ; dtend(1) = '1231'
     dtbeg(2) = '0102' ; dtend(2) = '1231'
  case ('Tables/CMIP5_6hrLev','Tables/CMIP5_6hrPlev','Tables/GeoMIP_6hrLev','Tables/GeoMIP_6hrPlev')
     ndt = 1 ; dtbeg(1) = '010100Z' ; dtend(1) = '123118Z'
  case ('Tables/CMIP5_3hr','Tables/CMIP5_cf3hr','Tables/GeoMIP_3hr','Tables/GeoMIP_cf3hr')
     ndt = 1 ; dtbeg(1) = '010100Z' ; dtend(1) = '123121Z'
  end select
  !
  select case (table)
  case ('Tables/CMIP5_OImon','Tables/GeoMIP_OImon')
     exists = .false.
     all_continue = .true.
     do year1 = begyr,endyr
        do year2 = endyr,begyr,-1
           do idt = 1,ndt
              write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,a,''-'',i4.4,a,''.nc'')') &
                   trim(case),&
                   trim(comp),&
                   trim(cesm_var)//'_nh',&
                   year1,trim(dtbeg(idt)),&
                   year2,trim(dtend(idt))
              inquire(file=checkname,exist=exists)
              if (exists) then
                 nc_nfilesnh(ivar) = nc_nfilesnh(ivar) + 1
                 ncfilenh(nc_nfilesnh(ivar),ivar) = checkname
              endif
              write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,a,''-'',i4.4,a,''.nc'')') &
                   trim(case),&
                   trim(comp),&
                   trim(cesm_var)//'_sh',&
                   year1,trim(dtbeg(idt)),&
                   year2,trim(dtend(idt))
              inquire(file=checkname,exist=exists)
              if (exists) then
                 nc_nfilessh(ivar) = nc_nfilessh(ivar) + 1
                 ncfilesh(nc_nfilessh(ivar),ivar) = checkname
              endif
           enddo
        enddo
     enddo
     !
     all_continue = all_continue.and.(nc_nfilesnh(ivar) /= 0)
     all_continue = all_continue.and.(nc_nfilessh(ivar) /= 0)
     write(*,*) 'build_filenames all_continue: ',all_continue
     if (all_continue) write(*,'(''nfiles NH: '',20i5)') nc_nfilesnh(1:ivar)
     if (all_continue) write(*,'(''nfiles SH: '',20i5)') nc_nfilessh(1:ivar)
!     if (all_continue) write(*,'(''NH  files: '',10(a))') (trim(ncfilenh(i,ivar)),i=1,nc_nfilesnh(ivar))
!     if (all_continue) write(*,'(''SH  files: '',10(a))') (trim(ncfilenh(i,ivar)),i=1,nc_nfilessh(ivar))
  case default
     exists = .false.
     do year1 = begyr,endyr
        do year2 = endyr,begyr,-1
           do idt = 1,ndt
!              write(*,'(''dtbeg, dtend: '',i4,5x,a,10x,a)') idt,trim(dtbeg(idt)),trim(dtend(idt))
              write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,a,''-'',i4.4,a,''.nc'')') &
                   trim(case),&
                   trim(comp),&
                   trim(cesm_var),&
                   year1,trim(dtbeg(idt)),&
                   year2,trim(dtend(idt))
              inquire(file=checkname,exist=exists)
              if (exists) then
                 nc_nfiles(ivar) = nc_nfiles(ivar) + 1
                 ncfile(nc_nfiles(ivar),ivar) = checkname
              endif
           enddo
        enddo
     enddo
     !
     all_continue = all_continue.and.(nc_nfiles(ivar) /= 0)
     write(*,*) 'build_filenames all_continue: ',all_continue
     if (all_continue) write(*,'(''nfiles: '',100i5)') nc_nfiles(1:ivar)
!     if (all_continue) write(*,'('' files: '',100(a))') (trim(ncfile(i,ivar)),i=1,nc_nfiles(ivar))
  end select
end subroutine build_filenames
