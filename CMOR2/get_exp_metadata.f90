!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_exp_metadata
  use exp_info
  use mycmor_info
  implicit none
  !
  integer::iexp,refyear,refmon,refday
  !
  exp_found = 0
  do iexp = 1,num_exp
     if (trim(exp(iexp)%case) == trim(case_read)) then
        exp_found = iexp
     endif
  enddo
  if (exp_found /= 0) then
     write(*,*) 'MATCH experiment: ',trim(exp(exp_found)%case),' ',trim(case_read),' ',trim(exp(exp_found)%model_id)
     write(*,*) trim(exp(exp_found)%run_refcase),' ',trim(exp(exp_found)%repotag),' ',trim(exp(exp_found)%compset),' ',trim(exp(exp_found)%rip_code)
  else
     write(*,*) 'case ',case_read(1:len_trim(case_read)),' NOT FOUND. Dying.'
     stop
  endif
  if (exp(exp_found)%forcing(1:) == 'unknown forcings') then
     write(*,*) 'FORCINGS UNDEFINED. PLEASE FIX. STOPPING.'
     stop
  endif
  !
  parent_found = 0
  do iexp = 1,num_exp
     if (trim(exp(exp_found)%run_refcase) == trim(exp(iexp)%case)) then
        parent_found = iexp
     endif
  enddo
  if (parent_found /= 0) then
     write(*,*) 'MATCH parent: ',trim(case_read),' ',trim(exp(parent_found)%case),' ',trim(exp(parent_found)%rip_code)
     mycmor%parent_experiment_id  = exp(parent_found)%expt_id(1:)
     mycmor%parent_experiment_rip = exp(parent_found)%rip_code(1:)
     read(exp(exp_found)%run_refdate(1: 4),'(i2.2)') refyear
     mycmor%branch_time = refyear
     read(exp(exp_found)%run_refdate(6: 7),'(i2.2)') refmon
     read(exp(exp_found)%run_refdate(9:10),'(f2.2)') refday
  else
     if (trim(exp(exp_found)%run_refcase) == "N/A") then
        mycmor%parent_experiment_id  = "N/A"
        mycmor%parent_experiment_rip = "N/A"
        mycmor%branch_time = 0.
     else
        write(*,*) 'parent ',case_read(1:len_trim(case_read)),' NOT FOUND. Dying.'
        stop
     endif
  endif
  write(*,'(''Experiment metadata loaded'')')
end subroutine get_exp_metadata
