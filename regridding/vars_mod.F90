MODULE vars_mod
  !---------------------------------------------------------------------------
  !   This module declares 'global' variables.
  !
  !   CVS:$Id: vars_mod.F90,v 1.2 2005/03/08 16:15:04 strandwg Exp $
  !---------------------------------------------------------------------------

  USE kinds_mod

  !***************************************************************************

  IMPLICIT NONE
  SAVE

  !***************************************************************************
  INTEGER (KIND=int_kind),PARAMETER::max_vars = 10
  CHARACTER(LEN=char_len):: src_file,out_form
  CHARACTER(LEN=char_len),DIMENSION(max_vars):: src_var, src_var2
  CHARACTER(LEN=char_len),DIMENSION(max_vars):: dst_file
  CHARACTER(LEN=char_len),DIMENSION(max_vars):: dst_var
  CHARACTER(LEN=char_len):: grid_var
  INTEGER (KIND=int_kind):: nk,nvars
  INTEGER(KIND=int_kind),DIMENSION(max_vars)::src_km
  CHARACTER(LEN=char_len):: map_file
  CHARACTER(LEN=char_len):: z_file
  LOGICAL (KIND=log_kind):: dst_lclobber,dst_z_standard

  INTEGER (KIND=int_kind):: dst_imt, dst_jmt, dst_km
  !*****************************************************************************

END MODULE vars_mod
