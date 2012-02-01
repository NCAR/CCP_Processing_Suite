!!$program driver
!!$  use xwalk_info
!!$  !
!!$  implicit none
!!$  character(len=256)::xw_file
!!$  !
!!$  xw_file = 'xwalk_Amon.txt'
!!$  !  read(*,*) xw_file
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
  integer::iostat,i,j,ixw
  integer,dimension(20)::icol,iblk
  logical::does_exist
  character(len=256)::instring
  !
  ! Get crosswalk information
  !
  do i = 1,max_entries
     xw(i)%table(1:)        = ' '
     xw(i)%entry(1:)        = ' '
     xw(i)%sname(1:)        = ' '
     xw(i)%dims(1:)         = ' '
     xw(i)%realm(1:)        = ' '
     xw(i)%comment(1:)      = ' '
     xw(i)%cesm_vars(:)(1:) = ' '
     xw(i)%ncesm_vars       = 0
  enddo
  !
  inquire(file=xw_file,exist=does_exist)
  if (.not.(does_exist)) then
     write(*,*) 'Cannot find ',trim(xw_file),'. Dying.'
     stop
  endif
  open(20,file=xw_file,form='formatted')
  iostat = 0 ; ixw = 0
  do while (iostat == 0)
     instring(1:)  = ' '
     read(20,'(a)',iostat=iostat) instring
     if (iostat == 0) then
        icol = 0 ; iblk = 0
        ixw  = ixw + 1
        !         1         2         3
        !123456789012345678901234567890
        !
        !6HrLev:hus:specific_humidity:longitude latitude alevel time1:atmos:2:Q PS:Q interpolated to standard plevs
        !6HrLev:ps:surface_air_pressure:longitude latitude time1:atmos:1:PS:PS no change
        !6HrLev:ta:air_temperature:longitude latitude alevel time1:atmos:2:T PS:T interpolated to standard plevs
        !6HrLev:ua:eastward_wind:longitude latitude alevel time1:atmos:2:U PS:U interpolated to standard plevs
        !6HrLev:va:northward_wind:longitude latitude alevel time1:atmos:2:V PS:V interpolated to standard plevs
        !
        j = 1
        do i = 1,len_trim(adjustl(instring))
           if (instring(i:i) == ':') then
              icol(j) = i
              j = j + 1
           endif
        enddo
        xw(ixw)%table(1:)         = instring(1        :icol(1)-1)
        xw(ixw)%entry(1:)         = instring(icol(1)+1:icol(2)-1)
        xw(ixw)%sname(1:)         = instring(icol(2)+1:icol(3)-1)
        xw(ixw)%dims(1:)          = instring(icol(3)+1:icol(4)-1)
        xw(ixw)%realm(1:)         = instring(icol(4)+1:icol(5)-1)
        read(instring(icol(5)+1:icol(6)-1),*) xw(ixw)%ncesm_vars
        xw(ixw)%comment(1:)       = instring(icol(7)+1:len_trim(instring))
        !
        if (xw(ixw)%ncesm_vars .gt. 1) then
           iblk(1) = icol(6)
           j = 2
           do i = icol(6)+1,icol(7)-1
              if (instring(i:i) == ' ') then
                 iblk(j) = i
                 j = j + 1
              endif
           enddo
           iblk(xw(ixw)%ncesm_vars+1) = icol(7)
           do i = 1,xw(ixw)%ncesm_vars
              xw(ixw)%cesm_vars(i) = instring(iblk(i)+1:iblk(i+1)-1)
           enddo
        else
           xw(ixw)%cesm_vars(1) = instring(icol(6)+1:icol(7)-1)
        endif
     endif
  enddo
  close(20)
  num_xw = ixw
  write(*,'(''xwalk file loaded : '',i5,'' entries.'')') num_xw
!!$  do i = 1,ixw
!!$     select case (xw(i)%ncesm_vars)
!!$        case ( 1 )
!!$           write(*,'(a8,''|'',,a20,''|'',i2,''|'',a12,''|'',a50)') &
!!$                xw(i)%table(1:),&
!!$                xw(i)%entry(1:),&
!!$                xw(i)%ncesm_vars,(xw(i)%cesm_vars(j)(1:),j=1,xw(i)%ncesm_vars),&
!!$                xw(i)%comment
!!$        case ( 2 )
!!$           write(*,'(a8,''|'',,a20,''|'',i2,''|'',a12,'':'',a12,''|'',a50)') &
!!$                xw(i)%table(1:),&
!!$                xw(i)%entry(1:),&
!!$                xw(i)%ncesm_vars,(xw(i)%cesm_vars(j)(1:),j=1,xw(i)%ncesm_vars),&
!!$                xw(i)%comment
!!$        case ( 3 )
!!$           write(*,'(a8,''|'',,a20,''|'',i2,''|'',2(a12,'':''),a12,''|'',a50)') &
!!$                xw(i)%table(1:),&
!!$                xw(i)%entry(1:),&
!!$                xw(i)%ncesm_vars,(xw(i)%cesm_vars(j)(1:),j=1,xw(i)%ncesm_vars),&
!!$                xw(i)%comment
!!$        case ( 4 )
!!$           write(*,'(a8,''|'',,a20,''|'',i2,''|'',3(a12,'':''),a12,''|'',a50)') &
!!$                xw(i)%table(1:),&
!!$                xw(i)%entry(1:),&
!!$                xw(i)%ncesm_vars,(xw(i)%cesm_vars(j)(1:),j=1,xw(i)%ncesm_vars),&
!!$                xw(i)%comment
!!$        case ( 5 )
!!$           write(*,'(a8,''|'',,a20,''|'',i2,''|'',4(a12,'':''),a12,''|'',a50)') &
!!$                xw(i)%table(1:),&
!!$                xw(i)%entry(1:),&
!!$                xw(i)%ncesm_vars,(xw(i)%cesm_vars(j)(1:),j=1,xw(i)%ncesm_vars),&
!!$                xw(i)%comment
!!$        end select
!!$  enddo
end subroutine load_xwalk
