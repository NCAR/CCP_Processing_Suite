!
! local
!
integer::i,j,length,index_info,status
!
! find info
!
if (trim(var_name) /= 'global') then
  index_info = -1
  do i=1,var_counter
     if(var_info(i)%name.eq.var_name.and.var_info(i)%ncid.eq.ncid) then
        index_info = i
        exit
     endif
   end do
   if(index_info.eq.-1) then
      write(*,*) 'ERROR: ',trim(var_name),' is not in ncid',ncid
      write(*,*) 'STOPPING'
      stop
   endif
   status = nf_get_att_real(ncid, &
        var_info(index_info)%id, &
        att_name,&
        NF_FLOAT, &
        att_val)
   if(status.ne.nf_noerr) call handle_err(status)
else
   status = nf_get_att_real(ncid, &
        NF_GLOBAL, &
        att_name,&
        NF_FLOAT, &
        att_val)
   if(status.ne.nf_noerr) call handle_err(status)
endif

