!!$program drive_system
!!$  implicit none
!!$  call add_global_metadata
!!$end program drive_system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_global_metadata
  use cmor_users_functions
  use exp_info
  use mycmor_info
  !
  implicit none
  integer::error_flag
  character(len=256)::whoami,prochost,ccps_rev,ccps_date,ccps_uuid,info_file
  character(len=10)::pdate,ptime
  logical::exists
  !
  ! Add acknowledgements
  !
  if (exp(exp_found)%loc(1:2) == 'NC') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(mycmor%ack_NC))
  if (exp(exp_found)%loc(1:2) == 'NE') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(mycmor%ack_NE))
  if (exp(exp_found)%loc(1:2) == 'OR') error_flag = cmor_set_cur_dataset_attribute("acknowledgements",trim(mycmor%ack_OR))
  !
  ! Add case name, repo tag, and compset
  !
  error_flag = cmor_set_cur_dataset_attribute("cesm_casename",trim(adjustl(exp(exp_found)%case)))
  error_flag = cmor_set_cur_dataset_attribute("cesm_repotag" ,trim(adjustl(exp(exp_found)%repotag)))
  error_flag = cmor_set_cur_dataset_attribute("cesm_compset" ,trim(adjustl(exp(exp_found)%compset)))
  !
  ! Add grid information
  !
  if (exp(exp_found)%grid(1:1) /= ' ') error_flag = cmor_set_cur_dataset_attribute("resolution",trim(adjustl(exp(exp_found)%grid)))
  !
  ! Add additional forcing information
  !
  if (mycmor%forcing_note(1:1) /= ' ') error_flag = cmor_set_cur_dataset_attribute("forcing_note",trim(adjustl(mycmor%forcing_note)))
  !
  ! Add cmor_version for HTAP2 (odd)
  !
  if (index(mycmor%table_file,'CCMI1') /= 0) then
     error_flag = cmor_set_cur_dataset_attribute("initialization_description"," ")
     error_flag = cmor_set_cur_dataset_attribute("physics_description"," ")
  endif
  !
  info_file = 'Info_in.'//trim(case_read)//'.'//trim(comp_read)
  inquire(file=trim(info_file),exist=exists)
  if (exists) then
     call date_and_time(date=pdate,time=ptime)
     open(30,file=trim(info_file),form='formatted')
     read(30,'(a)') whoami
     read(30,'(a)') prochost
     read(30,'(a)') ccps_rev
     read(30,'(a)') ccps_date
     read(30,'(a)') ccps_uuid
     close(30)
     error_flag = cmor_set_cur_dataset_attribute("processed_by ",trim(whoami)//" on "//trim(prochost)//" at "//pdate//"-"//ptime)
     if (error_flag /= 0) then
        write(*,*) "Error globalMD: ",error_flag
     endif
     error_flag = cmor_set_cur_dataset_attribute("processing_code_information ",trim(ccps_rev)//" "//trim(ccps_date)//" "//trim(ccps_uuid))
     if (error_flag /= 0) then
        write(*,*) "Error globalMD: ",error_flag
     endif
  else
     write(*,*) "Information file: ",trim(info_file)," missing. Stop."
  endif
end subroutine add_global_metadata
