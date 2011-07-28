!!$program driver
!!$  implicit none
!!$  character(len=256)::rip
!!$  integer::realization,initialization_method,physics_version
!!$  integer::i,j,k
!!$  !
!!$  do i = 1,9
!!$     do j = 1,9
!!$        do k = 1,9
!!$           write(rip,'(''r'',i1,''i'',i1,''p'',i1)') i,j,k
!!$           call parse_rip(rip,realization,initialization_method,physics_version)
!!$           write(*,'(''RIP: '',a,'' '',3i3)') rip(1:len_trim(rip)),realization,initialization_method,physics_version
!!$        enddo
!!$     enddo
!!$  enddo
!!$  !
!!$end program driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parse_rip
  !
  ! Parse RIP code into its components
  ! 
  use cmor_info
  use exp_info
  implicit none
  !
  ! Accepts one- or two-digit R, I, and P, i.e., one of:
  ! 123456789
  ! r1i1p1
  ! r1i1p99
  ! r1i99p1
  ! r1i99p99
  ! r99i1p1
  ! r99i1p99
  ! r99i99p1
  ! r99i99p99
  !
  character(len=512)::work
  integer::i,len
  !
  work(1:) = ' '
  work(1:) = adjustl(exp(exp_found)%rip_code(1:len_trim(exp(exp_found)%rip_code)))
  len      = len_trim(work)
  !
  do i = 1,len
     if ((work(i:i) == 'r').or.(work(i:i) == 'i').or.(work(i:i) == 'p')) work(i:i) = ' '
  enddo
  !
  write(*,*) 'WORK: ',work(1:len_trim(work))
  read(work,*) cmor%realization,cmor%initialization_method,cmor%physics_version
end subroutine parse_rip
