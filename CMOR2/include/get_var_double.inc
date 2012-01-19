!
! local
!
integer i,j,length,index_info
integer status,start(10),count(10)
logical time_dependent
!
! find info
!
index_info = -1
do j=1,file_counter
  do i=1,var_counter
     length = len_trim(var_info(i)%name)
     if(var_info(i)%name(:length).eq.var_name.and.var_info(i)%file_name.eq.nc_filename(j)) then
        index_info = i
        exit
     endif
  end do
end do
if(index_info.eq.-1) then
   write(*,*) 'ERROR: ',trim(var_name),' is not in ',trim(nc_filename(1))
   write(*,*) 'STOPPING'
   stop
endif
!
! check if the variable is time-dependent
!
time_dependent = .false.
if(index(var_info(index_info)%format,'time').ne.0)&
  time_dependent = .true.
!
if(time_dependent) then
   do i=1,var_info(index_info)%nvdims-1
      start(i) = 1
      count(i) = var_info(index_info)%vdims(i)
   end do
   start(var_info(index_info)%nvdims) = time_counter
   count(var_info(index_info)%nvdims) = 1
   status = nf_get_vara_double(ncid&
                                ,var_info(index_info)%id&
                                ,start&
                                ,count&
                                ,array)
else
   status = nf_get_var_double (ncid&
                                ,var_info(index_info)%id&
                                ,array)
endif
!
if(status.ne.nf_noerr) call handle_err(status)
!
