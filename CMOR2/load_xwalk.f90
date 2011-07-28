!!$program driver
!!$  use xwalk_info
!!$  !
!!$  implicit none
!!$  character(len=256)::xw_file
!!$  !
!!$  xw(:)%varin2d(1:) = ' '
!!$  xw(:)%units2d(1:) = ' '
!!$  xw(:)%positive2d(1:) = ' '
!!$  xw(:)%entry2d(1:) = ' '
!!$  !
!!$  read(*,*) xw_file
!!$  !
!!$  call load_xwalk(xw_file)
!!$  !
!!$end program driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine load_xwalk(xw_file)
  use xwalk_info
  implicit none
  !
  character(len=256),intent(in)::xw_file
  !
  integer::iostat,i,j,ixw,icol,ibng
  character(len=256)::instring,csv_file
  !
  !  csv_file = trim(xw_file)//'.csv'
  !
  ! Get table information
  !
  xw(:)%varin2d(1:) = ' '
  xw(:)%units2d(1:) = ' '
  xw(:)%entry2d(1:) = ' '
  open(20,file=trim(adjustl(xw_file)),form='formatted')
  iostat = 0 ; ixw = 0
  do while (iostat == 0)
     instring(1:)  = ' '
     read(20,'(a)',iostat=iostat) instring
     if (iostat == 0) then
        !         1         2         3
        !123456789012345678901234567890
        !ci:FREQZM
        !clivi:TGCLDIWP
        !clt:CLDTOT
        !clwvi:TGCLDLWP
        !evspsbl:QFLX
        !hfls:LHFLX
        !hfss:SHFLX
        !hurs:RHREFHT
        !huss:QREFHT
        !prc:PRECC
        !
        icol = index(instring,':')
        ixw  = ixw + 1
        xw(ixw)%entry2d(1:) = instring(1:icol-1)
        xw(ixw)%varin2d(1:) = instring(icol+1:len_trim(instring))
     endif
  enddo
!  ixw = ixw - 1
!!$  select case (trim(adjustl(var_info(n)%name)))
!!$  case ('FSDTOA','FSN200','FSN200C','FSNIRTOA','FSNRTOAC','FSNRTOAS','FLUT','FLUTC','FSNS','FSNSC','FSNT','FSNTC','FSNTOA','FSNTOAC','FSUTOA')
!!$     positive = 'up'
!!$     write(*,*) trim(adjustl(var_info(n)%name)),' up'
!!$  case ('FSDS','FSDSC','FDL','FDLC','FLDS','FLDSC','FLN200','FLN200C','FLNS','FLNSC','FLNT','FLNTC','FUL','FULC')
!!$     positive = 'down'
!!$     write(*,*) trim(adjustl(var_info(n)%name)),' down'
!!$  case default
!!$  end select
  num_xw = ixw
  close(20)
!!$  do i = 1,ixw
!!$     write(*,'(3(a,'',''),a)') &
!!$             'XW',&
!!$             trim(adjustl(xw(i)%entry2d)),&
!!$             trim(adjustl(xw(i)%varin2d))
!!$  enddo
end subroutine load_xwalk
