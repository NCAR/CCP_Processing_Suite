!!$program DRIVE_file_output_times
!!$  use output_times_info
!!$  implicit none
!!$  integer::ntimes
!!$  !
!!$  ntimes = 2388 ; call file_output_times(ntimes) ; write(*,*) tidx1(1),tidx2(1)
!!$  ntimes = 1872 ; call file_output_times(ntimes) ; write(*,*) tidx1(1),tidx2(1)
!!$  ntimes = 1140 ; call file_output_times(ntimes) ; write(*,*) tidx1(1),tidx2(1)
!!$  ntimes = 1152 ; call file_output_times(ntimes) ; write(*,*) tidx1(1),tidx2(1)
!!$end program DRIVE_file_output_times
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine file_output_times(ntimes)
  use output_times_info
  implicit none
  integer,intent(in)::ntimes
  !
  if (ntimes == 2388) then          ! RCP from 2101-2299, use all times
     nchunks(1) = 1
     tidx1(1:nchunks(1)) = (/   1/)
     tidx2(1:nchunks(1)) = (/2388/)
  endif
  if (ntimes == 1872) then          ! 20C from 1850-2005, use all times
     nchunks(1) = 1
     tidx1(1:nchunks(1)) = (/   1/)
     tidx2(1:nchunks(1)) = (/1872/)
  endif
  if (ntimes == 1140) then          ! RCP from 2006-2100, use all times
     nchunks(1) = 1
     tidx1(1:nchunks(1)) = (/   1/)
     tidx2(1:nchunks(1)) = (/1140/)
  endif
  if (ntimes == 1152) then          ! RCP from 2005-2100, skip 2005
     nchunks(1) = 1
     tidx1(1:nchunks(1)) = (/  13/)
     tidx2(1:nchunks(1)) = (/1152/)
  endif
end subroutine file_output_times
