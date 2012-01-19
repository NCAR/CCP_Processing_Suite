!
! local
!
integer::i,j,length,index_info,status
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
   status = nf_put_att_int(ncid, &
        var_info(index_info)%id, &
        att_name,&
        NF_INT, &
        1, &
        att_val)
   if(status.ne.nf_noerr) call handle_err(status)
else
   status = nf_put_att_int(ncid, &
        NF_GLOBAL, &
        att_name,&
        NF_INT, &
        1, &
        att_val)
   if(status.ne.nf_noerr) call handle_err(status)
endif

