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
  integer::i,j,year1,year2
  character(len=256)::checkname,dtbeg,dtend
  logical::exists
  !
  write(*,*) 'Entering build_filenames: ',trim(case),' ',trim(comp),' ',trim(cesm_var),ivar,begyr,endyr,trim(table)
  !
  select case (table)
  case ('Tables/CMIP5_Amon','Tables/CMIP5_Lmon','Tables/CMIP5_LImon','Tables/CMIP5_Omon','Tables/CMIP5_OImon','Tables/CMIP5_aero','Tables/CMIP5_cfMon')
     dtbeg = '01' ; dtend = '12'
  case ('Tables/CMIP5_day','Tables/CMIP5_cfDay')
     dtbeg = '0101' ; dtend = '1231'
  case ('Tables/CMIP5_6hrLev','Tables/CMIP5_6hrPlev')
     dtbeg = '010100Z' ; dtend = '123118Z'
  case ('Tables/CMIP5_3hr','Tables/CMIP5_cf3hr')
     dtbeg = '010100Z' ; dtend = '123121Z'
  end select
  !
  exists = .false.
  do year1 = begyr,endyr
     do year2 = endyr,begyr,-1
        write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,a,''-'',i4.4,a,''.nc'')') &
             trim(case),&
             trim(comp),&
             trim(cesm_var),&
             year1,trim(dtbeg),&
             year2,trim(dtend)
        inquire(file=checkname,exist=exists)
        all_continue = all_continue.or.exists
        if (exists) then
           nc_nfiles(ivar) = nc_nfiles(ivar) + 1
           ncfile(nc_nfiles(ivar),ivar) = checkname
        endif
     enddo
  enddo
  !
  write(*,*) 'build_filenames all_continue: ',all_continue
  if (all_continue) write(*,'(''nfiles: '',10i5)') nc_nfiles(1:ivar)
  if (all_continue) write(*,'('' files: '',10(a))') (trim(ncfile(i,ivar)),i=1,nc_nfiles(ivar))
end subroutine build_filenames
