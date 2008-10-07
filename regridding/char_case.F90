MODULE char_case

  !-----------------------------------------------------------------------------
  !   convert character and strings to upper/lower case
  !   CVS:$Id: char_case.F90,v 1.1 2004/11/10 21:23:25 strandwg Exp $
  !-----------------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

  !-----------------------------------------------------------------------------
  !   module variables
  !-----------------------------------------------------------------------------

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
       lower_letters = 'abcdefghijklmnopqrstuvwxyz', &
       upper_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  FUNCTION upper(STR)

    CHARACTER(LEN=*), INTENT(IN) :: STR
    CHARACTER(LEN=LEN_TRIM(STR)) :: upper
    INTEGER :: i, j

    DO i = 1,LEN_TRIM(STR)
       j = SCAN(lower_letters, STR(i:i))
       IF (j /= 0) THEN
          upper(i:i) = upper_letters(j:j)
       ELSE
          upper(i:i) = STR(i:i)
       END IF
    END DO

  END FUNCTION upper

  !*****************************************************************************

  FUNCTION lower(STR)

    CHARACTER(LEN=*), INTENT(IN) :: STR
    CHARACTER(LEN=LEN_TRIM(STR)) :: lower
    INTEGER :: i, j

    DO i = 1,LEN_TRIM(STR)
       j = SCAN(upper_letters, STR(i:i))
       IF (j /= 0) THEN
          lower(i:i) = lower_letters(j:j)
       ELSE
          lower(i:i) = STR(i:i)
       END IF
    END DO

  END FUNCTION lower

  !*****************************************************************************

END MODULE char_case
