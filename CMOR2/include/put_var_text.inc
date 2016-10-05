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
do i=1,var_counter
   if ((var_info(i)%name.eq.var_name).and.(var_info(i)%ncid.eq.ncid)) then
      index_info = i
      exit
   endif
end do
if(index_info.eq.-1) then
   write(*,*) 'ERROR: ',trim(var_name),' NOT FOUND in NCID: ',ncid
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
   status = nf_put_vara_text(ncid&
        ,var_info(index_info)%id&
        ,start&
        ,count&
        ,array)
else
   status = nf_put_var_text (ncid&
        ,var_info(index_info)%id&
        ,array)
endif
!
if(status.ne.nf_noerr) call handle_err(status)
!
