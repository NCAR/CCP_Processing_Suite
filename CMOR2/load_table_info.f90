!!$program driver
!!$  use mycmor_info
!!$  use table_info
!!$  !
!!$  implicit none
!!$  !
!!$  mycmor%table_file = 'CMIP5_Amon'
!!$  !  read(*,*) table_file
!!$  !
!!$  call load_table_info
!!$  !
!!$end program driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine load_table_info
  use mycmor_info
  use table_info
  implicit none
  !
  integer::iostat,i,j,id(10),ixw,ixt,icol,ibng
  logical::does_exist
  character(len=256)::instring
  !
  table(:)%variable_entry(1:) = ' '
  table(:)%axis_entry(1:) = ' '
  !
  table(:)%approx_interval(1:) = ' '
  table(:)%axis(1:) = ' '
  table(:)%baseURL(1:) = ' '
  table(:)%bounds_values(1:) = ' '
  table(:)%cell_measures(1:) = ' '
  table(:)%cell_methods(1:) = ' '
  table(:)%cf_version(1:) = ' '
  table(:)%climatology(1:) = ' '
  table(:)%cmor_version(1:) = ' '
  table(:)%comment(1:) = ' '
  table(:)%coords_attrib(1:) = ' '
  table(:)%dimensions(1:) = ' '
  table(:)%expt_id_ok(1:) = ' '
  table(:)%forcings(1:) = ' '
  table(:)%formula(1:) = ' '
  table(:)%frequency(1:) = ' '
  table(:)%generic_levels(1:) = ' '
  table(:)%index_only(1:) = ' '
  table(:)%long_name(1:) = ' '
  table(:)%missing_value(1:) = ' '
  table(:)%modeling_realm(1:) = ' '
  table(:)%must_call_cmor_grid(1:) = ' '
  table(:)%must_have_bounds(1:) = ' '
  table(:)%ok_max_mean_abs(1:) = ' '
  table(:)%ok_min_mean_abs(1:) = ' '
  table(:)%out_name(1:) = ' '
  table(:)%positive(1:) = ' '
  table(:)%product(1:) = ' '
  table(:)%project_id(1:) = ' '
  table(:)%requested(1:) = ' '
  table(:)%requested_bounds(1:) = ' '
  table(:)%required_global_attributes(1:) = ' '
  table(:)%standard_name(1:) = ' '
  table(:)%stored_direction(1:) = ' '
  table(:)%table_date(1:) = ' '
  table(:)%table_id(1:) = ' '
  table(:)%tolerance(1:) = ' '
  table(:)%type(1:) = ' '
  table(:)%units(1:) = ' '
  table(:)%valid_max(1:) = ' '
  table(:)%valid_min(1:) = ' '
  table(:)%value(1:) = ' '
  table(:)%z_bounds_factors(1:) = ' '
  table(:)%z_factors(1:) = ' '
  !
  ! Get table information
  !
  mycmor%table_file = 'Tables/'//trim(mycmor%table_file)
  inquire(file=mycmor%table_file,exist=does_exist)
  if (.not.(does_exist)) then
     write(*,*) 'Cannot find ',trim(mycmor%table_file),'. Dying.'
     stop
  endif
  !
  open(20,file=mycmor%table_file,form='formatted')
  iostat = 0 ; ixt = 0
  do while (iostat == 0)
     instring(1:)  = ' '
     read(20,'(a)',iostat=iostat) instring
     if (iostat == 0) then
        !         1         2         3
        !123456789012345678901234567890
        !table_id: Table 6hrLev
        !modeling_realm: atmos
        !frequency: 6hr
        !variable_entry:    ta
        !modeling_realm:    atmos
        !standard_name:     air_temperature
        !dimensions:        longitude latitude alevel time1
        !units:           Pa
        !valid_min:         -700
        !valid_max:         1.00E+04
        !cell_measures:      area: areacella
        !type:              real
        !out_name:          ta
        !
        icol = index(instring,':')
        if (icol /= 0) then
           icol = icol + 1
           ibng = index(instring,'!')
           if (ibng /= 0) ibng = ibng - 1
           !
           select case (instring(1:icol-1))
           case ('variable_entry:')
              ixt = ixt + 1
              if (ibng /= 0) then
                 table(ixt)%variable_entry(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%variable_entry(1:) = adjustl(instring(icol:))
              endif
           case ('axis_entry:')
              ixt = ixt + 1
              if (ibng /= 0) then
                 table(ixt)%axis_entry(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%axis_entry(1:) = adjustl(instring(icol:))
              endif
           case ('axis:')
              if (ibng /= 0) then
                 table(ixt)%axis(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%axis(1:) = adjustl(instring(icol:))
              endif
           case ('bounds_values:')
              if (ibng /= 0) then
                 table(ixt)%bounds_values(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%bounds_values(1:) = adjustl(instring(icol:))
              endif
           case ('cell_measures:')
              if (ibng /= 0) then
                 table(ixt)%cell_measures(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%cell_measures(1:) = adjustl(instring(icol:))
              endif
           case ('cell_methods:')
              if (ibng /= 0) then
                 table(ixt)%cell_methods(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%cell_methods(1:) = adjustl(instring(icol:))
              endif
           case ('climatology:')
              if (ibng /= 0) then
                 table(ixt)%climatology(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%climatology(1:) = adjustl(instring(icol:))
              endif
           case ('comment:')
              if (ibng /= 0) then
                 table(ixt)%comment(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%comment(1:) = adjustl(instring(icol:))
              endif
           case ('coords_attrib:')
              if (ibng /= 0) then
                 table(ixt)%coords_attrib(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%coords_attrib(1:) = adjustl(instring(icol:))
              endif
           case ('dimensions:')
              if (ibng /= 0) then
                 table(ixt)%dimensions(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%dimensions(1:) = adjustl(instring(icol:))
              endif
           case ('expt_id_ok:')
              if (ibng /= 0) then
                 table(ixt)%expt_id_ok(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%expt_id_ok(1:) = adjustl(instring(icol:))
              endif
           case ('formula:')
              if (ibng /= 0) then
                 table(ixt)%formula(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%formula(1:) = adjustl(instring(icol:))
              endif
           case ('generic_levels:')
              if (ibng /= 0) then
                 table(ixt)%generic_levels(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%generic_levels(1:) = adjustl(instring(icol:))
              endif
           case ('index_only:')
              if (ibng /= 0) then
                 table(ixt)%index_only(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%index_only(1:) = adjustl(instring(icol:))
              endif
           case ('long_name:')
              if (ibng /= 0) then
                 table(ixt)%long_name(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%long_name(1:) = adjustl(instring(icol:))
              endif
           case ('modeling_realm:')
              if (ibng /= 0) then
                 table(ixt)%modeling_realm(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%modeling_realm(1:) = adjustl(instring(icol:))
              endif
           case ('must_call_cmor_grid:')
              if (ibng /= 0) then
                 table(ixt)%must_call_cmor_grid(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%must_call_cmor_grid(1:) = adjustl(instring(icol:))
              endif
           case ('must_have_bounds:')
              if (ibng /= 0) then
                 table(ixt)%must_have_bounds(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%must_have_bounds(1:) = adjustl(instring(icol:))
              endif
           case ('ok_max_mean_abs:')
              if (ibng /= 0) then
                 table(ixt)%ok_max_mean_abs(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%ok_max_mean_abs(1:) = adjustl(instring(icol:))
              endif
           case ('ok_min_mean_abs:')
              if (ibng /= 0) then
                 table(ixt)%ok_min_mean_abs(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%ok_min_mean_abs(1:) = adjustl(instring(icol:))
              endif
           case ('out_name:')
              if (ibng /= 0) then
                 table(ixt)%out_name(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%out_name(1:) = adjustl(instring(icol:))
              endif
           case ('positive:')
              if (ibng /= 0) then
                 table(ixt)%positive(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%positive(1:) = adjustl(instring(icol:))
              endif
           case ('requested:')
              if (ibng /= 0) then
                 table(ixt)%requested(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%requested(1:) = adjustl(instring(icol:))
              endif
           case ('requested_bounds:')
              if (ibng /= 0) then
                 table(ixt)%requested_bounds(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%requested_bounds(1:) = adjustl(instring(icol:))
              endif
           case ('standard_name:')
              if (ibng /= 0) then
                 table(ixt)%standard_name(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%standard_name(1:) = adjustl(instring(icol:))
              endif
           case ('stored_direction:')
              if (ibng /= 0) then
                 table(ixt)%stored_direction(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%stored_direction(1:) = adjustl(instring(icol:))
              endif
           case ('tolerance:')
              if (ibng /= 0) then
                 table(ixt)%tolerance(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%tolerance(1:) = adjustl(instring(icol:))
              endif
           case ('type:')
              if (ibng /= 0) then
                 table(ixt)%type(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%type(1:) = adjustl(instring(icol:))
              endif
           case ('units:')
              if (ibng /= 0) then
                 table(ixt)%units(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%units(1:) = adjustl(instring(icol:))
              endif
           case ('valid_max:')
              if (ibng /= 0) then
                 table(ixt)%valid_max(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%valid_max(1:) = adjustl(instring(icol:))
              endif
           case ('valid_min:')
              if (ibng /= 0) then
                 table(ixt)%valid_min(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%valid_min(1:) = adjustl(instring(icol:))
              endif
           case ('value:')
              if (ibng /= 0) then
                 table(ixt)%value(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%value(1:) = adjustl(instring(icol:))
              endif
           case ('z_bounds_factors:')
              if (ibng /= 0) then
                 table(ixt)%z_bounds_factors(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%z_bounds_factors(1:) = adjustl(instring(icol:))
              endif
           case ('z_factors:')
              if (ibng /= 0) then
                 table(ixt)%z_factors(1:) = adjustl(instring(icol:ibng))
              else
                 table(ixt)%z_factors(1:) = adjustl(instring(icol:))
              endif
           case ('approx_interval:')
              if (ibng /= 0) then
                 table(ixt)%approx_interval = adjustl(instring(icol:ibng))
              else
                 table(ixt)%approx_interval = adjustl(instring(icol:))
              endif
           case ('baseURL:')
              if (ibng /= 0) then
                 table(ixt)%baseURL = adjustl(instring(icol:ibng))
              else
                 table(ixt)%baseURL = adjustl(instring(icol:))
              endif
           case ('cf_version:')
              if (ibng /= 0) then
                 table(ixt)%cf_version = adjustl(instring(icol:ibng))
              else
                 table(ixt)%cf_version = adjustl(instring(icol:))
              endif
           case ('cmor_version:')
              if (ibng /= 0) then
                 table(ixt)%cmor_version = adjustl(instring(icol:ibng))
              else
                 table(ixt)%cmor_version = adjustl(instring(icol:))
              endif
           case ('forcings:')
              if (ibng /= 0) then
                 table(ixt)%forcings = adjustl(instring(icol:ibng))
              else
                 table(ixt)%forcings = adjustl(instring(icol:))
              endif
           case ('frequency:')
              if (ibng /= 0) then
                 table(ixt)%frequency = adjustl(instring(icol:ibng))
              else
                 table(ixt)%frequency = adjustl(instring(icol:))
              endif
           case ('missing_value:')
              if (ibng /= 0) then
                 table(ixt)%missing_value = adjustl(instring(icol:ibng))
              else
                 table(ixt)%missing_value = adjustl(instring(icol:))
              endif
           case ('product:')
              if (ibng /= 0) then
                 table(ixt)%product = adjustl(instring(icol:ibng))
              else
                 table(ixt)%product = adjustl(instring(icol:))
              endif
           case ('project_id:')
              if (ibng /= 0) then
                 table(ixt)%project_id = adjustl(instring(icol:ibng))
              else
                 table(ixt)%project_id = adjustl(instring(icol:))
              endif
           case ('required_global_attributes:')
              if (ibng /= 0) then
                 table(ixt)%required_global_attributes = adjustl(instring(icol:ibng))
              else
                 table(ixt)%required_global_attributes = adjustl(instring(icol:))
              endif
           case ('table_date:')
              if (ibng /= 0) then
                 table(ixt)%table_date = adjustl(instring(icol:ibng))
              else
                 table(ixt)%table_date = adjustl(instring(icol:))
              endif
           case ('table_id:')
              if (ibng /= 0) then
                 table(ixt)%table_id = adjustl(instring(icol:ibng))
              else
                 table(ixt)%table_id = adjustl(instring(icol:))
              endif
           case default
           end select
        endif
     endif
  enddo
  num_tab = ixt
  write(*,'(''CMOR table loaded : '',i5,'' entries.'')') num_tab
  close(20)
!!$  open(21,file=mycmor%table_file(1:len_trim(mycmor%table_file))//'.txt')
!!$  write(21,*) 'variable_entry,out_name,standard_name,long_name,dimensions,type,units,cell_measures,cell_methods,modeling_realm,must_have_bounds,ok_max_mean_abs,ok_min_mean_abs,positive,stored_direction,valid_max,valid_min'
!!$  do i = 1,ixt
!!$     if (table(i)%variable_entry(1:) == ' ') then
!!$        write(20,'(50(a,'',''),a)') &
!!$             'axis',&
!!$             trim(adjustl(table(i)%axis_entry)),&
!!$             trim(adjustl(table(i)%out_name)),&
!!$             trim(adjustl(table(i)%standard_name)),&
!!$             trim(adjustl(table(i)%long_name)),&
!!$             trim(adjustl(table(i)%dimensions)),&
!!$             trim(adjustl(table(i)%type)),&
!!$             trim(adjustl(table(i)%units)),&
!!$             trim(adjustl(table(i)%cell_measures)),&
!!$             trim(adjustl(table(i)%cell_methods)),&
!!$             trim(adjustl(table(i)%modeling_realm)),&
!!$             trim(adjustl(table(i)%must_have_bounds)),&
!!$             trim(adjustl(table(i)%ok_max_mean_abs)),&
!!$             trim(adjustl(table(i)%ok_min_mean_abs)),&
!!$             trim(adjustl(table(i)%positive)),&
!!$             trim(adjustl(table(i)%stored_direction)),&
!!$             trim(adjustl(table(i)%valid_max)),&
!!$             trim(adjustl(table(i)%valid_min))
!!$     else
!!$        write(21,'(50(a,'':''),a)') &
!!$             'variable',&
!!$             trim(adjustl(table(i)%variable_entry)),&
!!$             trim(adjustl(table(i)%out_name)),&
!!$             trim(adjustl(table(i)%standard_name)),&
!!$             trim(adjustl(table(i)%long_name)),&
!!$             trim(adjustl(table(i)%dimensions)),&
!!$             trim(adjustl(table(i)%type)),&
!!$             trim(adjustl(table(i)%units)),&
!!$             trim(adjustl(table(i)%cell_measures)),&
!!$             trim(adjustl(table(i)%cell_methods)),&
!!$             trim(adjustl(table(i)%modeling_realm)),&
!!$             trim(adjustl(table(i)%must_have_bounds)),&
!!$             trim(adjustl(table(i)%ok_max_mean_abs)),&
!!$             trim(adjustl(table(i)%ok_min_mean_abs)),&
!!$             trim(adjustl(table(i)%positive)),&
!!$             trim(adjustl(table(i)%stored_direction)),&
!!$             trim(adjustl(table(i)%valid_max)),&
!!$             trim(adjustl(table(i)%valid_min))
!!$     endif
!!$  enddo
end subroutine load_table_info
