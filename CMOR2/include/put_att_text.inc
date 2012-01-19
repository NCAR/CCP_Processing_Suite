!
! local
!
integer::i,j,length,index_info,status,att_len
!
! find info
!
if (trim(var_name) /= 'global') then
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
   att_len = len_trim(att_val)
   status = nf_put_att_text(ncid, &
        var_info(index_info)%id, &
        att_name,&
        att_len, &
        att_val)
   if(status.ne.nf_noerr) call handle_err(status)
else
   att_len = len_trim(att_val)
   status = nf_put_att_text(ncid, &
        nf_global, &
        att_name,&
        att_len, &
        att_val)
   if(status.ne.nf_noerr) call handle_err(status)
endif
