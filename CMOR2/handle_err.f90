subroutine handle_err(status)
  !
  implicit none
  include 'netcdf.inc'
  !
  integer::status
  if (status /= nf_noerr) then
     print *, nf_strerror(status)
     call abort()
  endif
  !
  return
end subroutine handle_err
