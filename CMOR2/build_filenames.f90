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
subroutine build_filenames(case,comp,cesm_var,ivar,begyr,endyr)
  use files_info
  !
  implicit none
  character(len=256),intent(in)::case,comp,cesm_var
  integer,intent(in)::ivar,begyr,endyr
  integer::i,j,year1,year2
  character(len=256)::checkname
  logical::exists
  !
  write(*,*) 'Entering build_filenames: ',trim(case),' ',trim(comp),' ',trim(cesm_var),ivar,begyr,endyr
  exists = .false.
  do year1 = begyr,endyr
     do year2 = endyr,begyr,-1
        write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,''01-'',i4.4,''12.nc'')') &
             trim(case),&
             trim(comp),&
             trim(cesm_var),&
             year1,year2
        inquire(file=checkname,exist=exists)
        all_continue = all_continue.or.exists
        if (exists) then
           nc_nfiles(ivar) = nc_nfiles(ivar) + 1
           ncfile(nc_nfiles(ivar),ivar) = checkname
        endif
     enddo
  enddo
  if (nc_nfiles(ivar) == 0) then
     exists = .false.
     do year1 = begyr,endyr
        do year2 = endyr,begyr,-1
           write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,''0102-'',i4.4,''1231.nc'')') &
                trim(case),&
                trim(comp),&
                trim(cesm_var),&
                year1,year2
           inquire(file=checkname,exist=exists)
           all_continue = all_continue.or.exists
           if (exists) then
              nc_nfiles(ivar) = nc_nfiles(ivar) + 1
              ncfile(nc_nfiles(ivar),ivar) = checkname
           endif
        enddo
     enddo
     exists = .false.
     do year1 = begyr,endyr
        do year2 = endyr,begyr,-1
           write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,''0101-'',i4.4,''1231.nc'')') &
                trim(case),&
                trim(comp),&
                trim(cesm_var),&
                year1,year2
           inquire(file=checkname,exist=exists)
           all_continue = all_continue.or.exists
           if (exists) then
              nc_nfiles(ivar) = nc_nfiles(ivar) + 1
              ncfile(nc_nfiles(ivar),ivar) = checkname
           endif
        enddo
     enddo
  endif
  !
  write(*,*) 'build_filenames all_continue: ',all_continue
  if (all_continue) write(*,'(''nfiles: '',10i5)') nc_nfiles(1:ivar)
  if (all_continue) write(*,'('' files: '',10(a))') (trim(ncfile(i,ivar)),i=1,nc_nfiles(ivar))
end subroutine build_filenames
