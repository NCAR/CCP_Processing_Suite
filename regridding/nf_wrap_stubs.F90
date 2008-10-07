MODULE nf_wrap_stubs

  USE kinds_mod

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  SUBROUTINE nf_wrap_stop

    CALL ABORT

  END SUBROUTINE nf_wrap_stop

  !*****************************************************************************

  SUBROUTINE nf_wrap_broadcast(STAT)

    INTEGER(kind=int_kind), INTENT(IN) :: STAT

  END SUBROUTINE nf_wrap_broadcast

  !*****************************************************************************

END MODULE nf_wrap_stubs
