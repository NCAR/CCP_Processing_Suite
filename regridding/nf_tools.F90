MODULE nf_tools

  !-----------------------------------------------------------------------------
  !   CVS:$Id: nf_tools.F90,v 1.1 2004/11/10 21:23:25 strandwg Exp $
  !-----------------------------------------------------------------------------

  USE kinds_mod
  USE array_tools
  USE nf_wrap

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  !   module parameters
  !-----------------------------------------------------------------------------

  INTEGER(KIND=int_kind), PARAMETER :: &
       max_dims = 7

  !-----------------------------------------------------------------------------
  !   derived types
  !-----------------------------------------------------------------------------

  TYPE nf_dim
     CHARACTER(LEN=char_len) :: name
     INTEGER(KIND=int_kind) :: len, id
  END TYPE nf_dim

  TYPE nf_char_att
     CHARACTER(LEN=char_len) :: name
     INTEGER(KIND=int_kind) :: len
     CHARACTER(LEN=char_len) :: text
  END TYPE nf_char_att

  TYPE nf_real_att
     CHARACTER(LEN=char_len) :: name
     INTEGER(KIND=int_kind) :: len
     REAL(KIND=real_kind) :: rval
     REAL(KIND=real_kind), DIMENSION(:), POINTER :: rvals
  END TYPE nf_real_att

  TYPE nf_dbl_att
     CHARACTER(LEN=char_len) :: name
     INTEGER(KIND=int_kind) :: len
     REAL(KIND=dbl_kind) :: dval
     REAL(KIND=dbl_kind), DIMENSION(:), POINTER :: dvals
  END TYPE nf_dbl_att

  TYPE nf_var
     CHARACTER(LEN=char_len) :: name
     INTEGER(KIND=int_kind) :: type, ndims, id
     TYPE(nf_dim), DIMENSION(max_dims) :: dimids
     TYPE(nf_char_att), DIMENSION(:), POINTER :: char_atts
     TYPE(nf_real_att), DIMENSION(:), POINTER :: real_atts
     TYPE(nf_dbl_att),  DIMENSION(:), POINTER :: dbl_atts
  END TYPE nf_var

  !-----------------------------------------------------------------------------
  !   generic interfaces
  !-----------------------------------------------------------------------------

  INTERFACE nf_dimlookup
     MODULE PROCEDURE nf_dimlookup_real, nf_dimlookup_dbl
  END INTERFACE

  INTERFACE nf_varlookup
     MODULE PROCEDURE &
          nf_varlookup_real_1D, nf_varlookup_dbl_1D, &
          nf_varlookup_real_2D, nf_varlookup_dbl_2D
  END INTERFACE

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  SUBROUTINE nf_dimlookup_real(NCID, VALS, dimname, dimid, MSG)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: NCID
    REAL(KIND=real_kind), DIMENSION(:), INTENT(IN) :: VALS
    CHARACTER(LEN=char_len), INTENT(OUT) :: dimname
    INTEGER(KIND=int_kind), INTENT(OUT) :: dimid
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: MSG

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=char_len) :: wrap_msg
    INTEGER(KIND=int_kind) :: valslen, ndims, dimlen, varid, varndims
    INTEGER(KIND=int_kind), DIMENSION(1) :: vardimids
    INTEGER(KIND=int_kind) :: stat
    REAL(KIND=real_kind), DIMENSION(SIZE(VALS)) :: dimvals

    IF (PRESENT(MSG)) THEN
       wrap_msg = 'nf_dimlookup_real :: ' // TRIM(MSG)
    ELSE
       wrap_msg = 'nf_dimlookup_real :: '
    END IF

    valslen = SIZE(VALS)

    CALL nf_inq_ndims_wrap(NCID, ndims, wrap_msg)

    DO dimid = 1,ndims
       CALL nf_inq_dim_wrap(NCID, dimid, dimname, dimlen, wrap_msg)
       IF (dimlen /= valslen) CYCLE

       CALL nf_inq_varid_wrap(NCID, dimname, varid, wrap_msg, NF_ENOTVAR, stat)
       IF (stat == NF_ENOTVAR) CYCLE

       CALL nf_inq_varndims_wrap(NCID, varid, varndims, wrap_msg)
       IF (varndims /= 1) CYCLE

       CALL nf_inq_vardimid_wrap(NCID, varid, vardimids, wrap_msg)
       IF (vardimids(1) /= dimid) CYCLE

       CALL nf_get_var_wrap(NCID, varid, dimvals, wrap_msg)
       IF (equal(VALS, dimvals)) RETURN
    END DO

    dimname = 'notfound'
    dimid = -1

  END SUBROUTINE nf_dimlookup_real

  !*****************************************************************************

  SUBROUTINE nf_dimlookup_dbl(NCID, VALS, dimname, dimid, MSG)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: NCID
    REAL(KIND=dbl_kind), DIMENSION(:), INTENT(IN) :: VALS
    CHARACTER(LEN=char_len), INTENT(OUT) :: dimname
    INTEGER(KIND=int_kind), INTENT(OUT) :: dimid
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: MSG

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=char_len) :: wrap_msg
    INTEGER(KIND=int_kind) :: valslen, ndims, dimlen, varid, varndims
    INTEGER(KIND=int_kind), DIMENSION(1) :: vardimids
    INTEGER(KIND=int_kind) :: stat
    REAL(KIND=dbl_kind), DIMENSION(SIZE(VALS)) :: dimvals

    IF (PRESENT(MSG)) THEN
       wrap_msg = 'nf_dimlookup_dbl :: ' // TRIM(MSG)
    ELSE
       wrap_msg = 'nf_dimlookup_dbl :: '
    END IF

    valslen = SIZE(VALS)

    CALL nf_inq_ndims_wrap(NCID, ndims, wrap_msg)

    DO dimid = 1,ndims
       CALL nf_inq_dim_wrap(NCID, dimid, dimname, dimlen, wrap_msg)
       IF (dimlen /= valslen) CYCLE

       CALL nf_inq_varid_wrap(NCID, dimname, varid, wrap_msg, NF_ENOTVAR, stat)
       IF (stat == NF_ENOTVAR) CYCLE

       CALL nf_inq_varndims_wrap(NCID, varid, varndims, wrap_msg)
       IF (varndims /= 1) CYCLE

       CALL nf_inq_vardimid_wrap(NCID, varid, vardimids, wrap_msg)
       IF (vardimids(1) /= dimid) CYCLE

       CALL nf_get_var_wrap(NCID, varid, dimvals, wrap_msg)
       IF (equal(VALS, dimvals)) RETURN
    END DO

    dimname = 'notfound'
    dimid = -1

  END SUBROUTINE nf_dimlookup_dbl

  !*****************************************************************************

  SUBROUTINE nf_varlookup_real_1D(NCID, VALS, varname, varid, MSG)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: NCID
    REAL(KIND=real_kind), DIMENSION(:), INTENT(IN) :: VALS
    CHARACTER(LEN=char_len), INTENT(OUT) :: varname
    INTEGER(KIND=int_kind), INTENT(OUT) :: varid
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: MSG

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=char_len) :: wrap_msg
    INTEGER(KIND=int_kind) :: valslen, nvars, dimlen, varndims
    INTEGER(KIND=int_kind), DIMENSION(1) :: vardimids
    REAL(KIND=real_kind), DIMENSION(SIZE(VALS)) :: varvals

    IF (PRESENT(MSG)) THEN
       wrap_msg = 'nf_varlookup_real_1D :: ' // TRIM(MSG)
    ELSE
       wrap_msg = 'nf_varlookup_real_1D :: '
    END IF

    valslen = SIZE(VALS)

    CALL nf_inq_nvars_wrap(NCID, nvars, wrap_msg)

    DO varid = 1,nvars
       CALL nf_inq_varndims_wrap(NCID, varid, varndims, wrap_msg)
       IF (varndims /= 1) CYCLE

       CALL nf_inq_vardimid_wrap(NCID, varid, vardimids, wrap_msg)
       CALL nf_inq_dimlen_wrap(NCID, vardimids(1), dimlen, wrap_msg)
       IF (dimlen /= valslen) CYCLE

       CALL nf_get_var_wrap(NCID, varid, varvals, wrap_msg)
       IF (equal(VALS, varvals)) THEN
          CALL nf_inq_varname_wrap(NCID, varid, varname, wrap_msg)
          RETURN
       END IF
    END DO

    varname = 'notfound'
    varid = -1

  END SUBROUTINE nf_varlookup_real_1D

  !*****************************************************************************

  SUBROUTINE nf_varlookup_dbl_1D(NCID, VALS, varname, varid, MSG)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: NCID
    REAL(KIND=dbl_kind), DIMENSION(:), INTENT(IN) :: VALS
    CHARACTER(LEN=char_len), INTENT(OUT) :: varname
    INTEGER(KIND=int_kind), INTENT(OUT) :: varid
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: MSG

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=char_len) :: wrap_msg
    INTEGER(KIND=int_kind) :: valslen, nvars, dimlen, varndims
    INTEGER(KIND=int_kind), DIMENSION(1) :: vardimids
    REAL(KIND=dbl_kind), DIMENSION(SIZE(VALS)) :: varvals

    IF (PRESENT(MSG)) THEN
       wrap_msg = 'nf_varlookup_dbl_1D :: ' // TRIM(MSG)
    ELSE
       wrap_msg = 'nf_varlookup_dbl_1D :: '
    END IF

    valslen = SIZE(VALS)

    CALL nf_inq_nvars_wrap(NCID, nvars, wrap_msg)

    DO varid = 1,nvars
       CALL nf_inq_varndims_wrap(NCID, varid, varndims, wrap_msg)
       IF (varndims /= 1) CYCLE

       CALL nf_inq_vardimid_wrap(NCID, varid, vardimids, wrap_msg)
       CALL nf_inq_dimlen_wrap(NCID, vardimids(1), dimlen, wrap_msg)
       IF (dimlen /= valslen) CYCLE

       CALL nf_get_var_wrap(NCID, varid, varvals, wrap_msg)
       IF (equal(VALS, varvals)) THEN
          CALL nf_inq_varname_wrap(NCID, varid, varname, wrap_msg)
          RETURN
       END IF
    END DO

    varname = 'notfound'
    varid = -1

  END SUBROUTINE nf_varlookup_dbl_1D

  !*****************************************************************************

  SUBROUTINE nf_varlookup_real_2D(NCID, VALS, varname, varid, MSG)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: NCID
    REAL(KIND=real_kind), DIMENSION(:,:), INTENT(IN) :: VALS
    CHARACTER(LEN=char_len), INTENT(OUT) :: varname
    INTEGER(KIND=int_kind), INTENT(OUT) :: varid
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: MSG

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=char_len) :: wrap_msg
    INTEGER(KIND=int_kind) :: nvars, dimlen, varndims
    INTEGER(KIND=int_kind), DIMENSION(2) :: valslen, vardimids
    REAL(KIND=real_kind), DIMENSION(SIZE(VALS,1),SIZE(VALS,2)) :: varvals

    IF (PRESENT(MSG)) THEN
       wrap_msg = 'nf_varlookup_real_2D :: ' // TRIM(MSG)
    ELSE
       wrap_msg = 'nf_varlookup_real_2D :: '
    END IF

    valslen = SHAPE(VALS)

    CALL nf_inq_nvars_wrap(NCID, nvars, wrap_msg)

    DO varid = 1,nvars
       CALL nf_inq_varndims_wrap(NCID, varid, varndims, wrap_msg)
       IF (varndims /= 2) CYCLE

       CALL nf_inq_vardimid_wrap(NCID, varid, vardimids, wrap_msg)
       CALL nf_inq_dimlen_wrap(NCID, vardimids(1), dimlen, wrap_msg)
       IF (dimlen /= valslen(1)) CYCLE
       CALL nf_inq_dimlen_wrap(NCID, vardimids(2), dimlen, wrap_msg)
       IF (dimlen /= valslen(2)) CYCLE

       CALL nf_get_var_wrap(NCID, varid, varvals, wrap_msg)
       IF (equal(VALS, varvals)) THEN
          CALL nf_inq_varname_wrap(NCID, varid, varname, wrap_msg)
          RETURN
       END IF
    END DO

    varname = 'notfound'
    varid = -1

  END SUBROUTINE nf_varlookup_real_2D

  !*****************************************************************************

  SUBROUTINE nf_varlookup_dbl_2D(NCID, VALS, varname, varid, MSG)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: NCID
    REAL(KIND=dbl_kind), DIMENSION(:,:), INTENT(IN) :: VALS
    CHARACTER(LEN=char_len), INTENT(OUT) :: varname
    INTEGER(KIND=int_kind), INTENT(OUT) :: varid
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: MSG

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=char_len) :: wrap_msg
    INTEGER(KIND=int_kind) :: nvars, dimlen, varndims
    INTEGER(KIND=int_kind), DIMENSION(2) :: valslen, vardimids
    REAL(KIND=dbl_kind), DIMENSION(SIZE(VALS,1),SIZE(VALS,2)) :: varvals

    IF (PRESENT(MSG)) THEN
       wrap_msg = 'nf_varlookup_dbl_2D :: ' // TRIM(MSG)
    ELSE
       wrap_msg = 'nf_varlookup_dbl_2D :: '
    END IF

    valslen = SHAPE(VALS)

    CALL nf_inq_nvars_wrap(NCID, nvars, wrap_msg)

    DO varid = 1,nvars
       CALL nf_inq_varndims_wrap(NCID, varid, varndims, wrap_msg)
       IF (varndims /= 2) CYCLE

       CALL nf_inq_vardimid_wrap(NCID, varid, vardimids, wrap_msg)
       CALL nf_inq_dimlen_wrap(NCID, vardimids(1), dimlen, wrap_msg)
       IF (dimlen /= valslen(1)) CYCLE
       CALL nf_inq_dimlen_wrap(NCID, vardimids(2), dimlen, wrap_msg)
       IF (dimlen /= valslen(2)) CYCLE

       CALL nf_get_var_wrap(NCID, varid, varvals, wrap_msg)
       IF (equal(VALS, varvals)) THEN
          CALL nf_inq_varname_wrap(NCID, varid, varname, wrap_msg)
          RETURN
       END IF
    END DO

    varname = 'notfound'
    varid = -1

  END SUBROUTINE nf_varlookup_dbl_2D

  !*****************************************************************************

END MODULE nf_tools
