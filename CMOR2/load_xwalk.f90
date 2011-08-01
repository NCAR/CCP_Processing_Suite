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
  integer,dimension(10)::icol,icom
  logical::does_exist
  character(len=256)::instring
  !
  ! Get crosswalk information
  !
  xw(:)%table(1:)         = ' '
  xw(:)%entry(1:)         = ' '
  xw(:)%standard_name(1:) = ' '
  xw(:)%realm(1:)         = ' '
  xw(:)%priority          = -9999.
  xw(:)%cesm_vars(:)(1:)  = ' '
  xw(:)%ncesm_vars        = 0
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
        icol = 0 ; icom = 0
        ixw  = ixw + 1
        !         1         2         3
        !123456789012345678901234567890
        !
        !Amon:1.0:ccb:air_pressure_at_convective_cloud_base:atmos:1:CLDBOT
        !Amon:1.0:cct:air_pressure_at_convective_cloud_top:atmos:1:CLDTOP
        ![...]
        !Amon:1.0:mc:atmosphere_net_upward_convective_mass_flux:atmos:2:CMFMC,CMFMCDZM
        !Amon:1.0:pr:precipitation_flux:atmos:2:PRECC,PRECL
        !Amon:1.0:prsn:snowfall_flux:atmos:2:PRECSC,PRECSL
        !
        j = 1
        do i = 1,len_trim(adjustl(instring))
           if (instring(i:i) == ':') then
              icol(j) = i
              j = j + 1
           endif
        enddo
        xw(ixw)%table(1:)         = instring(1:icol(1)-1)
        xw(ixw)%entry(1:)         = instring(icol(2)+1:icol(3)-1)
        xw(ixw)%standard_name(1:) = instring(icol(3)+1:icol(4)-1)
        xw(ixw)%realm(1:)         = instring(icol(4)+1:icol(5)-1)
        read(instring(icol(1)+1:icol(2)-1),*) xw(ixw)%priority
        read(instring(icol(5)+1:icol(6)-1),*) xw(ixw)%ncesm_vars
        !
        if (xw(ixw)%ncesm_vars .gt. 1) then
           icom(1) = icol(6)
           j = 2
           do i = icol(6)+1,len_trim(instring)
              if (instring(i:i) == ',') then
                 icom(j) = i
                 j = j + 1
              endif
           enddo
           icom(xw(ixw)%ncesm_vars+1) = len_trim(instring)+1
           do i = 1,xw(ixw)%ncesm_vars
              xw(ixw)%cesm_vars(i) = instring(icom(i)+1:icom(i+1)-1)
           enddo
        else
           xw(ixw)%cesm_vars(1) = instring(icol(6)+1:len_trim(instring))
        endif
     endif
  enddo
  close(20)
  num_xw = ixw
  write(*,'(''xwalk file loaded : '',i5,'' entries.'')') num_xw
!!$  do i = 1,ixw
!!$     write(*,'(a8,''|'',,a20,''|'',f3.1,''|'',i2,''|'',5(a12,'', ''))') &
!!$          xw(i)%table(1:),&
!!$          xw(i)%entry(1:),&
!!$          xw(i)%priority,&
!!$          xw(i)%ncesm_vars,(xw(i)%cesm_vars(j)(1:),j=1,xw(i)%ncesm_vars)
!!$  enddo
end subroutine load_xwalk
