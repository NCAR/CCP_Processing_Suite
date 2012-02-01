subroutine define_atm_axes(dimensions)
  !
  ! Define all axes that will be needed, as defined by the table entry
  !
  ! 'time' must be last
  !
  use cmor_users_functions
  use files_info
  use table_info
  use grid_info
  use mycmor_info
  !
  implicit none
  character(len=256),intent(in)::dimensions
  !
  integer::i,j,idim,error_flag,ilev
  integer,dimension(0:10)::idxb
  !
  ! Parse "dimensions" to find names and how many
  !
  naxes = 1 ; idxb = 1 ; dimnames(:)(1:) = ' '
  do i = 1,len_trim(adjustl(dimensions))
     if (dimensions(i:i) == ' ') then
        idxb(naxes) = i
        naxes = naxes + 1
     endif
  enddo
  !
  dimunits(:)(1:) = ' '
  idxb(0)     = 0
  idxb(naxes) = len_trim(dimensions)+1
  do i = 1,naxes
     dimnames(i)(1:) = dimensions(idxb(i-1)+1:idxb(i)-1)
     do j = 1,num_tab
        if (dimnames(i) == table(j)%axis_entry) dimunits(i) = table(j)%units
     enddo
     if (dimnames(i) == 'time') dimunits(i) = time_units
  enddo
  !
!  do i = 1,naxes
!     write(*,*) 'DIMS: ',dimnames(i)(1:32),' UNITS: ',dimunits(i)(1:32)
!  enddo
  axis_ids = 0 ; idim = 1
  !
  ! Define 'time' first
  !
  select case (mycmor%table_file)
  case ('Tables/CMIP5_Amon')
     axis_ids(idim) = cmor_axis(  &
          table=mycmor%table_file,&
          table_entry='time',     &
          units=time_units,       &
          interval='30 days')
  case ('Tables/CMIP5_day')
     axis_ids(idim) = cmor_axis(  &
          table=mycmor%table_file,&
          table_entry='time',     &
          units=time_units,       &
          interval='1 day')
  end select
!  write(*,'('' dimension: '',a,'' defined: '',i4)') 'time',axis_ids(idim)
  idim = idim + 1
  !
  do i = 1,naxes
     select case(dimnames(i))
     case ('latitude')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             units=dimunits(i),            &
             length=SIZE(atm_lats),        &
             coord_vals=atm_lats,          &
             cell_bounds=atm_lats_bnds)
!        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('longitude')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(atm_lons),        &
             units=dimunits(i),            &
             coord_vals=atm_lons,          &
             cell_bounds=atm_lons_bnds)
!        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('plevs')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(atm_plev17),       &
             units=dimunits(i),            &
             coord_vals=atm_plev17)
!        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('plev8')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(atm_plev8),       &
             units=dimunits(i),            &
             coord_vals=atm_plev8)
!        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('plev7')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(atm_plev7),       &
             units=dimunits(i),            &
             coord_vals=atm_plev7)
!        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('plev3')
        axis_ids(idim) = cmor_axis(        &
             table=mycmor%table_file,      &
             table_entry=dimnames(i),      &
             length=SIZE(atm_plev3),       &
             units=dimunits(i),            &
             coord_vals=atm_plev3)
!        write(*,'('' dimension: '',a,'' defined: '',i4)') trim(dimnames(i)),axis_ids(idim)
        idim = idim + 1
     case ('alevel')
        ilev = cmor_axis(                        &
             table=mycmor%table_file,            &
             table_entry='standard_hybrid_sigma',&
             length=SIZE(atm_levs),              &
             units='1',                          &
             coord_vals=atm_levs/1000.,          &
             cell_bounds=atm_levs_bnds/1000.)
        error_flag = cmor_zfactor(    &
             zaxis_id=ilev,     &
             zfactor_name='p0',       &
             units='Pa',              &
             zfactor_values=p0)
        error_flag = cmor_zfactor(       &
             zaxis_id=ilev,        &
             zfactor_name='b',           &
             axis_ids= (/ilev/),   &
             units=' ',                  &
             zfactor_values=b_coeff,     &
             zfactor_bounds=b_coeff_bnds)
        error_flag = cmor_zfactor(       &
             zaxis_id=ilev,        &
             zfactor_name='a',           &
             axis_ids= (/ilev/),   &
             units=' ',                  &
             zfactor_values=a_coeff,     &
             zfactor_bounds=a_coeff_bnds)
        zfactor_id = cmor_zfactor(       &
             zaxis_id=ilev,        &
             zfactor_name='ps',          &
             axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3)/), &
             units='Pa')
        axis_ids(idim) = ilev
!!        write(*,'('' dimension: '',a,'' defined: '',i4)') 'standard_hybrid_sigma',axis_ids(idim)
        idim = idim + 1
     case ('alevhalf')
        ilev = cmor_axis(                        &
             table=mycmor%table_file,            &
             table_entry='standard_hybrid_sigma',&
             length=SIZE(atm_ilevs),             &
             units='1',                          &
             coord_vals=atm_ilevs/1000.,         &
             cell_bounds=atm_ilevs_bnds/1000.)
        error_flag = cmor_zfactor(    &
             zaxis_id=ilev,     &
             zfactor_name='p0',       &
             units='Pa',              &
             zfactor_values=p0)
        error_flag = cmor_zfactor(       &
             zaxis_id=ilev,        &
             zfactor_name='b',           &
             axis_ids= (/ilev/),   &
             units=' ',                  &
             zfactor_values=b_coeff,     &
             zfactor_bounds=b_coeff_bnds)
        error_flag = cmor_zfactor(       &
             zaxis_id=ilev,        &
             zfactor_name='a',           &
             axis_ids= (/ilev/),   &
             units=' ',                  &
             zfactor_values=a_coeff,     &
             zfactor_bounds=a_coeff_bnds)
        zfactor_id = cmor_zfactor(       &
             zaxis_id=ilev,        &
             zfactor_name='ps',          &
             axis_ids=(/axis_ids(1),axis_ids(2),axis_ids(3)/), &
             units='Pa')
        axis_ids(idim) = ilev
!!        write(*,'('' dimension: '',a,'' defined: '',i4)') 'standard_hybrid_sigma',axis_ids(idim)
        idim = idim + 1
     end select
  enddo
!  write(*,'(''CMOR axes defined, axis_ids: '',5i5)') (axis_ids(i),i=1,naxes)
end subroutine define_atm_axes
