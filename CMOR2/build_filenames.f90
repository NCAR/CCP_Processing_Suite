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
subroutine build_filenames(ixw,ivar,year1,year2)
  use exp_info
  use xwalk_info
  use files_info
  !
  implicit none
  integer,intent(in)::ivar,ixw,year1,year2
  integer::i,j,ifile
  character(len=512)::checkname
  logical::exists
  !
  ifile = 1
  do j = year1,year2
     do i = year1,year2
        exists = .false.
        write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,''01-'',i4.4,''12.nc'')') &
             trim(case_read),&
             trim(comp_read),&
             trim(xw(ixw)%cesm_vars(ivar)),i,j
        inquire(file=checkname,exist=exists)
        if (exists) then
           ncfile(ifile,ivar) = checkname
           nc_nfiles(ivar) = ifile
           ifile = ifile + 1
        endif
     enddo
  enddo
end subroutine build_filenames
