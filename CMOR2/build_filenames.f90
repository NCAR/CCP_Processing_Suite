!!$program DRIVE_build_filenames
!!$  implicit none
!!$  character(len=256)::case_read,comp_read
!!$  character(len=256),dimension(10)::cesm_vars
!!$  integer::begin,end,ivar,nvars
!!$  !
!!$  case_read = 'b40.rcp8_5.1deg.005'
!!$  comp_read = 'cam2.h1'
!!$  cesm_vars(1) = 'PRECC'
!!$  cesm_vars(2) = 'PRECL'
!!$  begin = 2006
!!$  end   = 2100
!!$  nvars = 2
!!$  !
!!$  call build_filenames(case_read,comp_read,cesm_vars,nvars,begin,end)
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
  exists(:,:) = .false.
  !
  !  write(*,'('' build_filenames: '',2i5,'' case '',a,'' comp '',a,'' cesm_var '',a,'' b,e '',2i4)') &
  !       ivar,ixw,trim(case_read),trim(comp_read),trim(xw(ixw)%cesm_vars(ivar)),exp(exp_found)%begyr,exp(exp_found)%endyr
  write(ncfile(1,ivar),'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,''0101-'',i4.4,''1231.nc'')') &
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
        if (year1 == 2000) year1 = 2006
        if ((year1 == 2090).and.(year2 == 2099)) year2 = 2100
        write(ncfile(ifile,ivar),'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,''0101-'',i4.4,''1231.nc'')') &
             trim(case_read),&
             trim(comp_read),&
             trim(xw(ixw)%cesm_vars(ivar)),&
             year1,year2
!           write(*,*) idec,idec*10,(idec*10)+9
!           write(*,'('' ncfile : '',a)') trim(ncfile(ifile,ivar))
        inquire(file=trim(ncfile(ifile,ivar)),exist=exists(ifile,ivar))
        nc_nfiles(ivar) = ifile
        ifile = ifile + 1
     enddo
  endif
!  do ivar = 1,xw(ixw)%ncesm_vars
!  do ifile = 1,nc_nfiles(ivar)
!     if (exists(ifile,ivar)) write(*,'(a,'' GOOD'')') trim(ncfile(ifile,ivar))
!  enddo
!  enddo
end subroutine build_filenames
