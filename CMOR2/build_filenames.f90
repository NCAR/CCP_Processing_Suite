!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_filenames(case,comp,cesm_var,ivar,runbeg,runend,table)
  use files_info
  !
  implicit none
  character(len=256),intent(in)::case,comp,cesm_var,table
  integer,intent(in)::ivar,runbeg,runend
  integer::i,j,year1,year2,ndt,idt
  character(len=512)::checkname
  character(len=256),dimension(2)::dtbeg,dtend
  logical::exists
  !
   write(*,*) 'Entering build_filenames: ',trim(case),' ',trim(comp),' ',trim(cesm_var),ivar,runbeg,runend,trim(table)
  !
  select case (table)
  case ('Tables/CMIP5_Amon','Tables/CMIP5_Lmon','Tables/CMIP5_LImon','Tables/CMIP5_Omon','Tables/CMIP5_OImon','Tables/CMIP5_aero','Tables/CMIP5_cfMon',&
        'Tables/PMIP3_Amon','Tables/PMIP3_Lmon','Tables/PMIP3_LImon','Tables/PMIP3_Omon','Tables/PMIP3_OImon',&
        'Tables/GeoMIP_Amon','Tables/GeoMIP_Lmon','Tables/GeoMIP_LImon','Tables/GeoMIP_Omon','Tables/GeoMIP_OImon','Tables/GeoMIP_aero','Tables/GeoMIP_cfMon',&
        'Tables/AEROCOM-ACC_2D-M','Tables/AEROCOM-ACC_3D-M','Tables/AEROCOM-ACC_2D-I',&
        'Tables/CCMI1_monthly','Tables/HTAP2_monthly')
     ndt = 2
     dtbeg(1) = '01' ; dtend(1) = '12'
     dtbeg(2) = '01' ; dtend(2) = '11'
  case ('Tables/CMIP5_day','Tables/CMIP5_cfDay','Tables/GeoMIP_day','Tables/GeoMIP_cfDay','Tables/PMIP3_day')
     ndt = 2
     dtbeg(1) = '0101' ; dtend(1) = '1231'
     dtbeg(2) = '0102' ; dtend(2) = '1231'
  case ('Tables/CMIP5_6hrLev','Tables/CMIP5_6hrPlev','Tables/GeoMIP_6hrLev','Tables/GeoMIP_6hrPlev')
     ndt = 2
     dtbeg(1) = '010100Z' ; dtend(1) = '123118Z'
     dtbeg(2) = '010106Z' ; dtend(2) = '123118Z'
  case ('Tables/CMIP5_3hr','Tables/CMIP5_cf3hr','Tables/GeoMIP_3hr','Tables/GeoMIP_cf3hr')
     ndt = 1 ; dtbeg(1) = '010100Z' ; dtend(1) = '123121Z'
  case ('Tables/CMIP5_cfSites')
     ndt = 1 ; dtbeg(1) = '001010000Z' ; dtend(1) = '12312330Z'
  end select
  !
  select case (table)
!  case ('Tables/CMIP5_OImon','Tables/GeoMIP_OImon','Tables/PMIP3_OImon')
  case ('Tables/CMIP5_OImon','Tables/GeoMIP_OImon','Tables/PMIP3_OImon','Tables/CMIP5_day')
     exists = .false.
     do year1 = runbeg,runend
        do year2 = runend,year1,-1
           do idt = 1,ndt
              write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,a,''-'',i4.4,a,''.nc'')') &
                   trim(case),&
                   trim(comp),&
                   trim(cesm_var)//'_nh',&
                   year1,trim(dtbeg(idt)),&
                   year2,trim(dtend(idt))
              inquire(file=checkname,exist=exists)
              if (exists) then
                 nc_nfiles_nh(ivar) = nc_nfiles_nh(ivar) + 1
                 ncfile_nh(nc_nfiles_nh(ivar),ivar) = checkname
              endif
              write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,a,''-'',i4.4,a,''.nc'')') &
                   trim(case),&
                   trim(comp),&
                   trim(cesm_var)//'_sh',&
                   year1,trim(dtbeg(idt)),&
                   year2,trim(dtend(idt))
              inquire(file=checkname,exist=exists)
              if (exists) then
                 nc_nfiles_sh(ivar) = nc_nfiles_sh(ivar) + 1
                 ncfile_sh(nc_nfiles_sh(ivar),ivar) = checkname
              endif
           enddo
        enddo
     enddo
     !
     all_continue = all_continue.and.(nc_nfiles_nh(ivar) /= 0)
     all_continue = all_continue.and.(nc_nfiles_sh(ivar) /= 0)
     write(*,*) 'build_filenames all_continue: ',all_continue
!     if (all_continue) write(*,'(''nfiles NH: '',20i5)') nc_nfiles_nh(1:ivar)
!     if (all_continue) write(*,'(''nfiles SH: '',20i5)') nc_nfiles_sh(1:ivar)
!     if (all_continue) write(*,'(''NH  files: '',10(a))') (trim(ncfile_nh(i,ivar)),i=1,nc_nfiles_nh(ivar))
!     if (all_continue) write(*,'(''SH  files: '',10(a))') (trim(ncfile_sh(i,ivar)),i=1,nc_nfiles_sh(ivar))
  case ('Tables/TAMIP_3hrCurt','Tables/TAMIP_3hrMlev','Tables/TAMIP_3hrPlev','Tables/TAMIP_3hrSlev','Tables/TAMIP_sites')
     exists = .false.
     write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,''.nc'')') &
          trim(case),&
          trim(comp),&
          trim(cesm_var),&
          runbeg
     inquire(file=checkname,exist=exists)
!     write(*,'(''build_filenames: '',a,5x,i4,5x)') trim(checkname),year1
     if (exists) then
        nc_nfiles(ivar) = nc_nfiles(ivar) + 1
        ncfile(nc_nfiles(ivar),ivar) = checkname
     endif
     !
     all_continue = all_continue.and.(nc_nfiles(ivar) /= 0)
     write(*,*) 'build_filenames all_continue: ',all_continue
!     if (all_continue) write(*,'(''nfiles: '',100i5)') nc_nfiles(1:ivar)
!      if (all_continue) write(*,'('' files: '',100(a))') (trim(ncfile(i,ivar)),i=1,nc_nfiles(ivar))
  case default
     exists = .false.
     if (runbeg .lt. runend) then
        do year1 = runbeg,runend
           do year2 = runend,year1,-1
              do idt = 1,ndt
!                  write(*,'(''dtbeg, dtend: '',i4,5x,a,10x,a)') idt,trim(dtbeg(idt)),trim(dtend(idt))
                  write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,a,''-'',i4.4,a,''.nc'')') &
                      trim(case),&
                      trim(comp),&
                      trim(cesm_var),&
                      year1,trim(dtbeg(idt)),&
                      year2,trim(dtend(idt))
                 inquire(file=checkname,exist=exists)
!                     write(*,'(''build_filenames: '',a,5x,i4,5x)') trim(checkname),year1
                 if (exists) then
                    nc_nfiles(ivar) = nc_nfiles(ivar) + 1
                    ncfile(nc_nfiles(ivar),ivar) = checkname
                 endif
              enddo
           enddo
        enddo
     else if (runbeg .eq. runend) then
        write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,a,''-'',i4.4,a,''.nc'')') &
             trim(case),&
             trim(comp),&
             trim(cesm_var),&
             runbeg,trim(dtbeg(1)),&
             runend,trim(dtend(1))
        inquire(file=checkname,exist=exists)
        write(*,'(''build_filenames: '',a,5x,i4,5x)') trim(checkname),year1
        if (exists) then
           nc_nfiles(ivar) = nc_nfiles(ivar) + 1
           ncfile(nc_nfiles(ivar),ivar) = checkname
        endif
     else
        write(checkname,'(''data/'',a,''.'',a,''.'',a,''.'',i4.4,a,''-'',i4.4,a,''.nc'')') &
             trim(case),&
             trim(comp),&
             trim(cesm_var),&
             runbeg,trim(dtbeg(1)),&
             runend,trim(dtend(1))
        inquire(file=checkname,exist=exists)
        write(*,'(''build_filenames: '',a,5x,i4,5x)') trim(checkname),year1
        if (exists) then
           nc_nfiles(ivar) = nc_nfiles(ivar) + 1
           ncfile(nc_nfiles(ivar),ivar) = checkname
        endif
     endif
     !
     all_continue = all_continue.and.(nc_nfiles(ivar) /= 0)
!     write(*,*) 'build_filenames all_continue: ',all_continue
     if (all_continue) write(*,'(''nfiles: '',100i5)') nc_nfiles(1:ivar)
     if (all_continue) write(*,'('' files: '',100(a))') (trim(ncfile(i,ivar)),i=1,nc_nfiles(ivar))
  end select
end subroutine build_filenames
