!
! local
!
integer::i,length,index_info,status
!
! find info
!
if (trim(var_name) /= 'global') then
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
   status = nf_get_att_int(ncid, &
        var_info(index_info)%id, &
        att_name,&
        NF_INT, &
        att_val)
   if(status.ne.nf_noerr) call handle_err(status)
else
   status = nf_put_att_int(ncid, &
        NF_GLOBAL, &
        att_name,&
        NF_INT, &
        att_val)
   if(status.ne.nf_noerr) call handle_err(status)
endif

