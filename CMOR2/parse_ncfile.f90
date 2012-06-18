!!$program driver
!!$  implicit none
!!$  character(len=256)::ncfile
!!$  character(len=256)::case,comp,svar,time
!!$  integer::i,j,k
!!$  !
!!$  read(*,*) ncfile
!!$  call parse_ncfile(ncfile,case,comp,svar,time)
!!$  write(*,'(''case: '',a)') case(1:len_trim(case))
!!$  write(*,'(''comp: '',a)') comp(1:len_trim(comp))
!!$  write(*,'(''svar: '',a)') svar(1:len_trim(svar))
!!$  write(*,'(''time: '',a)') time(1:len_trim(time))
!!$  !
!!$end program driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parse_ncfile(ncfile,case,comp,svar,time)
  ! 
  ! Parse input CCSM4/CESM1 single-field format file into its components
  ! case name, component, variable name, time
  !
  implicit none
  !
  character(len=256),intent(in)::ncfile
  character(len=256),intent(out)::case,comp,svar,time
  !
  integer,parameter::maxdot = 50
  character(len=256)::work
  integer,dimension(maxdot)::idot
  integer::i,len,c1,c2,p1,p2,v1,v2,t1,t2,ndot,nslash
  !
  !          1         2         3         4         5
  ! 123456789012345678901234567890123456789012345678901234567890
  ! data/b40.20th.track1.1deg.012.cam2.h0.TS.185001-200512.nc
  !
  work(1:) = ' '
  work(1:) = adjustl(ncfile)
  len      = len_trim(work)
  !
  nslash = 0
  write(*,'(''parse_ncfile: '',a)') work(1:len)
  do i = 1,len
     if (work(i:i) == "/") then
        nslash = i + 1
     endif
  end do
  if (nslash.ne.0) then
     work(1:) = ncfile(nslash:len)
  else
     work(1:) = adjustl(ncfile)
  endif
  !
  write(*,'(''parse_ncfile: '',a)') work(1:len)
  !
  ndot   = 1
  len = len_trim(work)
  do i = 1,len
     if (work(i:i) == ".") then
        idot(ndot) = i
        ndot       = ndot + 1
     endif
  enddo
  ndot = ndot - 1
  c1 = 1
  do i = 2,ndot
     if (work(idot(i-1):idot(i)) == '.cam2.') then
        c2 = idot(i-1)-1
        p1 = idot(i-1)+1
        p2 = idot(i+1)-1
        v1 = idot(i+1)+1
        v2 = idot(i+2)-1
        t1 = idot(i+2)+1
        t2 = idot(i+3)-1
     endif
     if (work(idot(i-1):idot(i)) == '.pop.') then
        c2 = idot(i-1)-1
        p1 = idot(i-1)+1
        p2 = idot(i+1)-1
        v1 = idot(i+1)+1
        v2 = idot(i+2)-1
        t1 = idot(i+2)+1
        t2 = idot(i+3)-1
     endif
  enddo
  case(1:) = work(c1:c2)
  comp(1:) = work(p1:p2)
  svar(1:) = work(v1:v2)
  time(1:) = work(t1:t2)
  write(*,'(''parse_ncfile case '',a)') trim(case)
  write(*,'(''parse_ncfile comp '',a)') trim(comp)
  write(*,'(''parse_ncfile svar '',a)') trim(svar)
  write(*,'(''parse_ncfile time '',a)') trim(time)
end subroutine parse_ncfile
