!
! local
!
integer i,length,index_info
integer status,start(10),count(10)
logical time_dependent
!
! find info
!
index_info = -1
do i=1,var_counter
   length = len_trim(var_info(i)%name)
   if(var_info(i)%name(:length).eq.var_name.and.&
           var_info(i)%file_name.eq.netcdf_file_name(ncid)) then
      index_info = i
      exit
   endif
end do
if(index_info.eq.-1) then
   print *,' '
   print *,var_name,' is not in the list of defined variables'
   print *,'stopping'
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
   status = nf_get_vara_int(ncid&
                             ,var_info(index_info)%id&
                             ,start&
                             ,count&
                             ,array)
else
   status = nf_get_var_int (ncid&
                             ,var_info(index_info)%id&
                             ,array)
endif
!
if(status.ne.nf_noerr) call handle_err(status)

