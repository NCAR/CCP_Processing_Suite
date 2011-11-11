!!$program DRIVE_build_filenames
!!$  use exp_info
!!$  use xwalk_info
!!$  use files_info
!!$  implicit none
!!$  integer::i,j,ixw,ifile
!!$  !
!!$  case_read = 'b40.20th.track1.1deg.008'
!!$  comp_read = 'pop.h'
!!$  xw_found  = 1
!!$  exp_found = 1
!!$  ixw       = 1
!!$  xw(ixw)%cesm_vars(1) = 'TEMP'
!!$  xw(ixw)%ncesm_vars   = 1
!!$  exp(exp_found)%begyr = 1850
!!$  exp(exp_found)%endyr = 2005
!!$  !
!!$  do j = 1,1
!!$     do i = 1,1
!!$        call build_filenames(i,j)
!!$     enddo
!!$  enddo
!!$  do i = 1,xw(ixw)%ncesm_vars
!!$     do ifile = 1,nc_nfiles(i)
!!$        if (exists(ifile,i)) write(*,'(a,'' GOOD'')') trim(ncfile(ifile,i))
!!$     enddo
!!$  enddo
!!$end program DRIVE_build_filenames
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_filenames(ivar,ixw)
  use exp_info
  use xwalk_info
  use files_info
  !
  implicit none
  integer,intent(in)::ivar,ixw
  integer::dec_b,dec_e,idec,ifile,year1,year2
  !
  ncfile(:,:)(1:) = ' '
  exists(:,:)     = .false.
  nc_nfiles(:)    = 0
  !
  write(ncfile(1,ivar),'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,''01-'',i4.4,''12.nc'')') &
       trim(case_read),&
       trim(comp_read),&
       trim(xw(ixw)%cesm_vars(ivar)),&
       exp(exp_found)%begyr,&
       exp(exp_found)%endyr
  inquire(file=trim(ncfile(1,ivar)),exist=exists(1,ivar))
  if (.not.(exists(1,ivar))) then
     dec_b = exp(exp_found)%begyr/10
     dec_e = exp(exp_found)%endyr/10
     ifile = 1
     do idec = dec_b,dec_e
        year1 = idec*10
        year2 = (idec*10) + 9
        write(ncfile(ifile,ivar),'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,''01-'',i4.4,''12.nc'')') &
             trim(case_read),&
             trim(comp_read),&
             trim(xw(ixw)%cesm_vars(ivar)),&
             year1,year2
        inquire(file=trim(ncfile(ifile,ivar)),exist=exists(ifile,ivar))
        nc_nfiles(ivar) = ifile
        ifile = ifile + 1
     enddo
  endif
  if (exp(exp_found)%endyr == 2005) then
     ifile = ifile + 1
     write(ncfile(ifile,ivar),'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,''01-'',i4.4,''12.nc'')') &
          trim(case_read),&
          trim(comp_read),&
          trim(xw(ixw)%cesm_vars(ivar)),&
          year1,exp(exp_found)%endyr
     inquire(file=trim(ncfile(ifile,ivar)),exist=exists(ifile,ivar))
     nc_nfiles(ivar) = ifile
  endif
end subroutine build_filenames
