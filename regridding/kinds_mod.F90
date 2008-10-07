MODULE kinds_mod

  !-----------------------------------------------------------------------------
  !   This module defines kinds for common data types.
  !   CVS:$Id: kinds_mod.F90,v 1.1 2004/11/10 21:23:25 strandwg Exp $
  !-----------------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

  !*****************************************************************************

  INTEGER, PARAMETER :: char_len  = 256
  INTEGER, PARAMETER :: int_kind  = KIND(1)
  INTEGER, PARAMETER :: log_kind  = KIND(.TRUE.)
  INTEGER, PARAMETER :: real_kind = SELECTED_REAL_KIND(6)
  INTEGER, PARAMETER :: dbl_kind  = SELECTED_REAL_KIND(13)

  !*****************************************************************************

END MODULE kinds_mod
