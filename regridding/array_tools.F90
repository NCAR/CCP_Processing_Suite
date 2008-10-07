MODULE array_tools

  !-----------------------------------------------------------------------------
  !   CVS:$Id: array_tools.F90,v 1.1 2004/11/10 21:23:25 strandwg Exp $
  !-----------------------------------------------------------------------------

  USE kinds_mod
  USE msg_mod

  IMPLICIT NONE
  PRIVATE

  !-----------------------------------------------------------------------------
  !   public functions
  !-----------------------------------------------------------------------------

  PUBLIC :: rms, equal

  !-----------------------------------------------------------------------------
  !   generic interfaces
  !-----------------------------------------------------------------------------

  INTERFACE rms
     MODULE PROCEDURE &
          rms_1D_s, rms_1D_d, &
          rms_2D_s, rms_2D_d
  END INTERFACE

  INTERFACE equal
     MODULE PROCEDURE &
          equal_0D_ss, equal_0D_sd, equal_0D_ds, equal_0D_dd, &
          equal_1D_ss, equal_1D_sd, equal_1D_ds, equal_1D_dd, &
          equal_2D_ss, equal_2D_sd, equal_2D_ds, equal_2D_dd
  END INTERFACE

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  REAL(KIND=real_kind) FUNCTION rms_1D_s(A, WEIGHT, MASK)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=real_kind), DIMENSION(:), INTENT(IN) :: A
    REAL(KIND=real_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: MASK

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'rms_1D_s'

    INTEGER(KIND=int_kind) :: imax, i, cnt
    REAL(KIND=real_kind) :: work

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A) /= SIZE(WEIGHT)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A) /= SIZE(MASK)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    imax = SIZE(A)

    work = 0.0_real_kind

    IF (PRESENT(MASK)) THEN
       cnt = 0
       IF (PRESENT(WEIGHT)) THEN
          DO i = 1,imax
             IF (MASK(i)) THEN
                work = work + WEIGHT(i) * A(i) ** 2
                cnt = cnt + 1
             END IF
          END DO
       ELSE
          DO i = 1,imax
             IF (MASK(i)) THEN
                work = work + WEIGHT(i) * A(i) ** 2
                cnt = cnt + 1
             END IF
          END DO
       END IF
       IF (cnt /= 0) work = work / REAL(cnt,real_kind)
    ELSE
       IF (PRESENT(WEIGHT)) THEN
          DO i = 1,imax
             work = work + WEIGHT(i) * A(i) ** 2
          END DO
       ELSE
          DO i = 1,imax
             work = work + A(i) ** 2
          END DO
       END IF
       IF (imax /= 0) work = work / REAL(imax,real_kind)
    END IF

    rms_1D_s = SQRT(work)

  END FUNCTION rms_1D_s

  !*****************************************************************************

  REAL(KIND=dbl_kind) FUNCTION rms_1D_d(A, WEIGHT, MASK)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=dbl_kind), DIMENSION(:), INTENT(IN) :: A
    REAL(KIND=dbl_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: MASK

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'rms_1D_d'

    INTEGER(KIND=int_kind) :: imax, i, cnt
    REAL(KIND=dbl_kind) :: work

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A) /= SIZE(WEIGHT)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A) /= SIZE(MASK)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    imax = SIZE(A)

    work = 0.0_dbl_kind

    IF (PRESENT(MASK)) THEN
       cnt = 0
       IF (PRESENT(WEIGHT)) THEN
          DO i = 1,imax
             IF (MASK(i)) THEN
                work = work + WEIGHT(i) * A(i) ** 2
                cnt = cnt + 1
             END IF
          END DO
       ELSE
          DO i = 1,imax
             IF (MASK(i)) THEN
                work = work + WEIGHT(i) * A(i) ** 2
                cnt = cnt + 1
             END IF
          END DO
       END IF
       IF (cnt /= 0) work = work / REAL(cnt,dbl_kind)
    ELSE
       IF (PRESENT(WEIGHT)) THEN
          DO i = 1,imax
             work = work + WEIGHT(i) * A(i) ** 2
          END DO
       ELSE
          DO i = 1,imax
             work = work + A(i) ** 2
          END DO
       END IF
       IF (imax /= 0) work = work / REAL(imax,dbl_kind)
    END IF

    rms_1D_d = SQRT(work)

  END FUNCTION rms_1D_d

  !*****************************************************************************

  REAL(KIND=real_kind) FUNCTION rms_2D_s(A, WEIGHT, MASK)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=real_kind), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=real_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'rms_2D_s'

    INTEGER(KIND=int_kind) :: imax, jmax, i, j, cnt
    REAL(KIND=real_kind) :: work

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A,1) /= SIZE(WEIGHT,1) .OR. SIZE(A,2) /= SIZE(WEIGHT,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A,1) /= SIZE(MASK,1) .OR. SIZE(A,2) /= SIZE(MASK,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    imax = SIZE(A,1)
    jmax = SIZE(A,2)

    work = 0.0_real_kind

    IF (PRESENT(MASK)) THEN
       cnt = 0
       IF (PRESENT(WEIGHT)) THEN
          DO j = 1,jmax
             DO i = 1,imax
                IF (MASK(i,j)) THEN
                   work = work + WEIGHT(i,j) * A(i,j) ** 2
                   cnt = cnt + 1
                END IF
             END DO
          END DO
       ELSE
          DO j = 1,jmax
             DO i = 1,imax
                IF (MASK(i,j)) THEN
                   work = work + A(i,j) ** 2
                   cnt = cnt + 1
                END IF
             END DO
          END DO
       END IF
       IF (cnt /= 0) work = work / REAL(cnt,real_kind)
    ELSE
       IF (PRESENT(WEIGHT)) THEN
          DO j = 1,jmax
             DO i = 1,imax
                work = work + WEIGHT(i,j) * A(i,j) ** 2
             END DO
          END DO
       ELSE
          DO j = 1,jmax
             DO i = 1,imax
                work = work + A(i,j) ** 2
             END DO
          END DO
       END IF
       IF (imax /= 0 .AND. jmax /= 0) work = work / REAL(imax*jmax,real_kind)
    END IF

    rms_2D_s = SQRT(work)

  END FUNCTION rms_2D_s

  !*****************************************************************************

  REAL(KIND=dbl_kind) FUNCTION rms_2D_d(A, WEIGHT, MASK)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=dbl_kind), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=dbl_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'rms_2D_d'

    INTEGER(KIND=int_kind) :: imax, jmax, i, j, cnt
    REAL(KIND=dbl_kind) :: work

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A,1) /= SIZE(WEIGHT,1) .OR. SIZE(A,2) /= SIZE(WEIGHT,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A,1) /= SIZE(MASK,1) .OR. SIZE(A,2) /= SIZE(MASK,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    imax = SIZE(A,1)
    jmax = SIZE(A,2)

    work = 0.0_dbl_kind

    IF (PRESENT(MASK)) THEN
       cnt = 0
       IF (PRESENT(WEIGHT)) THEN
          DO j = 1,jmax
             DO i = 1,imax
                IF (MASK(i,j)) THEN
                   work = work + WEIGHT(i,j) * A(i,j) ** 2
                   cnt = cnt + 1
                END IF
             END DO
          END DO
       ELSE
          DO j = 1,jmax
             DO i = 1,imax
                IF (MASK(i,j)) THEN
                   work = work + A(i,j) ** 2
                   cnt = cnt + 1
                END IF
             END DO
          END DO
       END IF
       IF (cnt /= 0) work = work / REAL(cnt,dbl_kind)
    ELSE
       IF (PRESENT(WEIGHT)) THEN
          DO j = 1,jmax
             DO i = 1,imax
                work = work + WEIGHT(i,j) * A(i,j) ** 2
             END DO
          END DO
       ELSE
          DO j = 1,jmax
             DO i = 1,imax
                work = work + A(i,j) ** 2
             END DO
          END DO
       END IF
       IF (imax /= 0 .AND. jmax /= 0) work = work / REAL(imax*jmax,dbl_kind)
    END IF

    rms_2D_d = SQRT(work)

  END FUNCTION rms_2D_d

  !*****************************************************************************

  LOGICAL FUNCTION equal_0D_ss(A, B, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=real_kind), INTENT(IN) :: A
    REAL(KIND=real_kind), INTENT(IN) :: B
    REAL(KIND=real_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_0D_ss'

    REAL(KIND=real_kind) :: reltol_actual, abstol_actual
    REAL(KIND=real_kind) :: diff, mag

    IF (PRESENT(RELTOL)) THEN
       reltol_actual = RELTOL
    ELSE
       reltol_actual = 1.0e-4_real_kind
    END IF

    IF (PRESENT(ABSTOL)) THEN
       abstol_actual = ABSTOL
    ELSE
       abstol_actual = 0.0_real_kind
    END IF

    equal_0D_ss = .FALSE.

    diff = ABS(A - B)
    IF (diff < abstol_actual) equal_0D_ss = .TRUE.

    mag = MAX(ABS(A), ABS(B))
    IF (diff < mag * reltol_actual) equal_0D_ss = .TRUE.

  END FUNCTION equal_0D_ss

  !*****************************************************************************

  LOGICAL FUNCTION equal_0D_sd(A, B, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=real_kind), INTENT(IN) :: A
    REAL(KIND=dbl_kind),  INTENT(IN) :: B
    REAL(KIND=real_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_0D_sd'

    REAL(KIND=real_kind) :: B_real

    B_real = REAL(B, real_kind)
    equal_0D_sd = equal_0D_ss(A, B_real, RELTOL, ABSTOL)

  END FUNCTION equal_0D_sd

  !*****************************************************************************

  LOGICAL FUNCTION equal_0D_ds(A, B, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=dbl_kind),  INTENT(IN) :: A
    REAL(KIND=real_kind), INTENT(IN) :: B
    REAL(KIND=real_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_0D_ds'

    REAL(KIND=real_kind) :: A_real

    A_real = REAL(A, real_kind)
    equal_0D_ds = equal_0D_ss(A_real, B, RELTOL, ABSTOL)

  END FUNCTION equal_0D_ds

  !*****************************************************************************

  LOGICAL FUNCTION equal_0D_dd(A, B, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=dbl_kind), INTENT(IN) :: A
    REAL(KIND=dbl_kind), INTENT(IN) :: B
    REAL(KIND=dbl_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_0D_dd'

    REAL(KIND=dbl_kind) :: reltol_actual, abstol_actual
    REAL(KIND=dbl_kind) :: diff, mag

    IF (PRESENT(RELTOL)) THEN
       reltol_actual = RELTOL
    ELSE
       reltol_actual = 1.0e-8_dbl_kind
    END IF

    IF (PRESENT(ABSTOL)) THEN
       abstol_actual = ABSTOL
    ELSE
       abstol_actual = 0.0_dbl_kind
    END IF

    equal_0D_dd = .FALSE.

    diff = ABS(A - B)
    IF (diff < abstol_actual) equal_0D_dd = .TRUE.

    mag = MAX(ABS(A), ABS(B))
    IF (diff < mag * reltol_actual) equal_0D_dd = .TRUE.

  END FUNCTION equal_0D_dd

  !*****************************************************************************

  LOGICAL FUNCTION equal_1D_ss(A, B, WEIGHT, MASK, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=real_kind), DIMENSION(:), INTENT(IN) :: A
    REAL(KIND=real_kind), DIMENSION(:), INTENT(IN) :: B
    REAL(KIND=real_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: MASK
    REAL(KIND=real_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_1D_ss'

    REAL(KIND=real_kind) :: reltol_actual, abstol_actual
    REAL(KIND=real_kind) :: diff, mag

    IF (SIZE(A) /= SIZE(B)) THEN
       equal_1D_ss = .FALSE.
       RETURN
    END IF

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A) /= SIZE(WEIGHT)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A) /= SIZE(MASK)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    IF (PRESENT(RELTOL)) THEN
       reltol_actual = RELTOL
    ELSE
       reltol_actual = 1.0e-4_real_kind
    END IF

    IF (PRESENT(ABSTOL)) THEN
       abstol_actual = ABSTOL
    ELSE
       abstol_actual = 0.0_real_kind
    END IF

    equal_1D_ss = .FALSE.

    diff = rms(A - B, WEIGHT, MASK)
    IF (diff < abstol_actual) equal_1D_ss = .TRUE.

    mag = MAX(rms(A, WEIGHT, MASK), rms(B, WEIGHT, MASK))
    IF (diff < mag * reltol_actual) equal_1D_ss = .TRUE.

  END FUNCTION equal_1D_ss

  !*****************************************************************************

  LOGICAL FUNCTION equal_1D_sd(A, B, WEIGHT, MASK, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=real_kind), DIMENSION(:), INTENT(IN) :: A
    REAL(KIND=dbl_kind), DIMENSION(:), INTENT(IN) :: B
    REAL(KIND=real_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: MASK
    REAL(KIND=real_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_1D_sd'

    REAL(KIND=real_kind), DIMENSION(SIZE(B)) :: B_real

    IF (SIZE(A) /= SIZE(B)) THEN
       equal_1D_sd = .FALSE.
       RETURN
    END IF

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A) /= SIZE(WEIGHT)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A) /= SIZE(MASK)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    B_real = REAL(B, real_kind)
    equal_1D_sd = equal_1D_ss(A, B_real, WEIGHT, MASK, RELTOL, ABSTOL)

  END FUNCTION equal_1D_sd

  !*****************************************************************************

  LOGICAL FUNCTION equal_1D_ds(A, B, WEIGHT, MASK, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=dbl_kind), DIMENSION(:), INTENT(IN) :: A
    REAL(KIND=real_kind), DIMENSION(:), INTENT(IN) :: B
    REAL(KIND=real_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: MASK
    REAL(KIND=real_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_1D_ds'

    REAL(KIND=real_kind), DIMENSION(SIZE(A)) :: A_real

    IF (SIZE(A) /= SIZE(B)) THEN
       equal_1D_ds = .FALSE.
       RETURN
    END IF

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A) /= SIZE(WEIGHT)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A) /= SIZE(MASK)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    A_real = REAL(A, real_kind)
    equal_1D_ds = equal_1D_ss(A_real, B, WEIGHT, MASK, RELTOL, ABSTOL)

  END FUNCTION equal_1D_ds

  !*****************************************************************************

  LOGICAL FUNCTION equal_1D_dd(A, B, WEIGHT, MASK, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=dbl_kind), DIMENSION(:), INTENT(IN) :: A
    REAL(KIND=dbl_kind), DIMENSION(:), INTENT(IN) :: B
    REAL(KIND=dbl_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:), INTENT(IN), OPTIONAL :: MASK
    REAL(KIND=dbl_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_1D_dd'

    REAL(KIND=dbl_kind) :: reltol_actual, abstol_actual
    REAL(KIND=dbl_kind) :: diff, mag

    IF (SIZE(A) /= SIZE(B)) THEN
       equal_1D_dd = .FALSE.
       RETURN
    END IF

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A) /= SIZE(WEIGHT)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A) /= SIZE(MASK)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    IF (PRESENT(RELTOL)) THEN
       reltol_actual = RELTOL
    ELSE
       reltol_actual = 1.0e-8_dbl_kind
    END IF

    IF (PRESENT(ABSTOL)) THEN
       abstol_actual = ABSTOL
    ELSE
       abstol_actual = 0.0_dbl_kind
    END IF

    equal_1D_dd = .FALSE.

    diff = rms(A - B, WEIGHT, MASK)
    IF (diff < abstol_actual) equal_1D_dd = .TRUE.

    mag = MAX(rms(A, WEIGHT, MASK), rms(B, WEIGHT, MASK))
    IF (diff < mag * reltol_actual) equal_1D_dd = .TRUE.

  END FUNCTION equal_1D_dd

  !*****************************************************************************

  LOGICAL FUNCTION equal_2D_ss(A, B, WEIGHT, MASK, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=real_kind), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=real_kind), DIMENSION(:,:), INTENT(IN) :: B
    REAL(KIND=real_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK
    REAL(KIND=real_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_2D_ss'

    REAL(KIND=real_kind) :: reltol_actual, abstol_actual
    REAL(KIND=real_kind) :: diff, mag

    IF (SIZE(A,1) /= SIZE(B,1) .OR. SIZE(A,2) /= SIZE(B,2)) THEN
       equal_2D_ss = .FALSE.
       RETURN
    END IF

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A,1) /= SIZE(WEIGHT,1) .OR. SIZE(A,2) /= SIZE(WEIGHT,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A,1) /= SIZE(MASK,1) .OR. SIZE(A,2) /= SIZE(MASK,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    IF (PRESENT(RELTOL)) THEN
       reltol_actual = RELTOL
    ELSE
       reltol_actual = 1.0e-4_real_kind
    END IF

    IF (PRESENT(ABSTOL)) THEN
       abstol_actual = ABSTOL
    ELSE
       abstol_actual = 0.0_real_kind
    END IF

    equal_2D_ss = .FALSE.

    diff = rms(A - B, WEIGHT, MASK)
    IF (diff < abstol_actual) equal_2D_ss = .TRUE.

    mag = MAX(rms(A, WEIGHT, MASK), rms(B, WEIGHT, MASK))
    IF (diff < mag * reltol_actual) equal_2D_ss = .TRUE.

  END FUNCTION equal_2D_ss

  !*****************************************************************************

  LOGICAL FUNCTION equal_2D_sd(A, B, WEIGHT, MASK, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=real_kind), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=dbl_kind), DIMENSION(:,:), INTENT(IN) :: B
    REAL(KIND=real_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK
    REAL(KIND=real_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_2D_sd'

    REAL(KIND=real_kind), DIMENSION(SIZE(B,1), SIZE(B,2)) :: B_real

    IF (SIZE(A,1) /= SIZE(B,1) .OR. SIZE(A,2) /= SIZE(B,2)) THEN
       equal_2D_sd = .FALSE.
       RETURN
    END IF

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A,1) /= SIZE(WEIGHT,1) .OR. SIZE(A,2) /= SIZE(WEIGHT,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A,1) /= SIZE(MASK,1) .OR. SIZE(A,2) /= SIZE(MASK,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    B_real = REAL(B, real_kind)
    equal_2D_sd = equal_2D_ss(A, B_real, WEIGHT, MASK, RELTOL, ABSTOL)

  END FUNCTION equal_2D_sd

  !*****************************************************************************

  LOGICAL FUNCTION equal_2D_ds(A, B, WEIGHT, MASK, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=dbl_kind), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=real_kind), DIMENSION(:,:), INTENT(IN) :: B
    REAL(KIND=real_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK
    REAL(KIND=real_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_2D_ds'

    REAL(KIND=real_kind), DIMENSION(SIZE(A,1), SIZE(A,2)) :: A_real

    IF (SIZE(A,1) /= SIZE(B,1) .OR. SIZE(A,2) /= SIZE(B,2)) THEN
       equal_2D_ds = .FALSE.
       RETURN
    END IF

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A,1) /= SIZE(WEIGHT,1) .OR. SIZE(A,2) /= SIZE(WEIGHT,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A,1) /= SIZE(MASK,1) .OR. SIZE(A,2) /= SIZE(MASK,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    A_real = REAL(A, real_kind)
    equal_2D_ds = equal_2D_ss(A_real, B, WEIGHT, MASK, RELTOL, ABSTOL)

  END FUNCTION equal_2D_ds

  !*****************************************************************************

  LOGICAL FUNCTION equal_2D_dd(A, B, WEIGHT, MASK, RELTOL, ABSTOL)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=dbl_kind), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=dbl_kind), DIMENSION(:,:), INTENT(IN) :: B
    REAL(KIND=dbl_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: WEIGHT
    LOGICAL(KIND=log_kind), DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK
    REAL(KIND=dbl_kind), INTENT(IN), OPTIONAL :: RELTOL, ABSTOL

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'equal_2D_dd'

    REAL(KIND=dbl_kind) :: reltol_actual, abstol_actual
    REAL(KIND=dbl_kind) :: diff, mag

    IF (SIZE(A,1) /= SIZE(B,1) .OR. SIZE(A,2) /= SIZE(B,2)) THEN
       equal_2D_dd = .FALSE.
       RETURN
    END IF

    IF (PRESENT(WEIGHT)) THEN
       IF (SIZE(A,1) /= SIZE(WEIGHT,1) .OR. SIZE(A,2) /= SIZE(WEIGHT,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & WEIGHT')
       END IF
    END IF

    IF (PRESENT(MASK)) THEN
       IF (SIZE(A,1) /= SIZE(MASK,1) .OR. SIZE(A,2) /= SIZE(MASK,2)) THEN
          CALL msg_write(sub_name, 'size mismatch between A & MASK')
       END IF
    END IF

    IF (PRESENT(RELTOL)) THEN
       reltol_actual = RELTOL
    ELSE
       reltol_actual = 1.0e-8_dbl_kind
    END IF

    IF (PRESENT(ABSTOL)) THEN
       abstol_actual = ABSTOL
    ELSE
       abstol_actual = 0.0_dbl_kind
    END IF

    equal_2D_dd = .FALSE.

    diff = rms(A - B, WEIGHT, MASK)
    IF (diff < abstol_actual) equal_2D_dd = .TRUE.

    mag = MAX(rms(A, WEIGHT, MASK), rms(B, WEIGHT, MASK))
    IF (diff < mag * reltol_actual) equal_2D_dd = .TRUE.

  END FUNCTION equal_2D_dd

  !*****************************************************************************

END MODULE array_tools
