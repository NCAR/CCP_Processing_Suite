subroutine build_filenames(ivar)
  use exp_info
  use xwalk_info
  use files_info
  !
  implicit none
  integer,intent(in)::ivar
  !
  files%ncfile(:,:)(1:) = ' '
  !
  write(ncfile(ivar),'(''data/'',a,''.'',a,''.'',a,''.'',a,''0101-'',a,''1231.nc'')') &
       trim(case_read),&
       trim(comp_read),&
       trim(xw(ixw)%cesm_vars(ivar)),&
       exp(exp_found)%begin_end(1:4),&
       exp(exp_found)%begin_end(6:9)
  inquire(file=trim(ncfile(ivar)),exist=continue(ivar))
