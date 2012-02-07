!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine reset_netcdf_var
  !
  ! Reset JFL netcdf interface 'var_info' information
  !
  use counters_netcdf_jfl
  use interfaces_netcdf_jfl
  use definitions_netcdf_jfl
  !
  integer::i
  !
  dim_counter  = 0
  var_counter  = 0
  time_counter = 0
  file_counter = 0
  !
  dim_info(:)%name(1:)          = ' '
  dim_info(:)%file_name(1:)     = ' '
  dim_info(:)%id                = 0
  dim_info(:)%length            = 0
  !
  var_info(:)%name(1:)          = ' '
  var_info(:)%format(1:)        = ' '
  var_info(:)%type(1:)          = ' '
  var_info(:)%file_name(1:)     = ' '
  var_info(:)%units(1:)         = ' '
  var_info(:)%missing_value     = 0.
  var_info(:)%int_missing_value = 0
  var_info(:)%FillValue         = 0.
  var_info(:)%int_FillValue     = 0
  var_info(:)%id                = 0
  var_info(:)%nvdims            = 0
  var_info(:)%units(1:)         = ' '
  do i = 1,10
     var_info(:)%vdims(i)       = 0
  enddo
  !
  nc_filename(:)(1:)            = ' '
end subroutine reset_netcdf_var
