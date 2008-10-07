MODULE const_mod

  !-----------------------------------------------------------------------------
  !   CVS:$Id: const_mod.F90,v 1.1 2004/11/10 21:23:25 strandwg Exp $
  !-----------------------------------------------------------------------------

  USE kinds_mod

  REAL (KIND=dbl_kind), PARAMETER :: &
       tiny  = 1e-6_dbl_kind, &
       c0    = 0.0_dbl_kind, &
       p125  = 0.125_dbl_kind, &
       p25   = 0.25_dbl_kind, &
       p5    = 0.5_dbl_kind, &
       p1    = 0.1_dbl_kind, &
       p01   = 0.01_dbl_kind, &
       p001  = 0.001_dbl_kind, &
       p0001 = 0.0001_dbl_kind, &
       c1    = 1.0_dbl_kind, &
       c2    = 2.0_dbl_kind, &
       c3    = 3.0_dbl_kind, &
       c4    = 4.0_dbl_kind, &
       c10   = 10.0_dbl_kind, &
       c90   = 90.0_dbl_kind, &
       c91   = 91.0_dbl_kind, &
       c100  = 100.0_dbl_kind, &
       c1000 = 1000.0_dbl_kind, &
       pi    = 3.1415926535897932385_dbl_kind, &
       pi2   = c2*pi, &
       pih   = p5*pi

END MODULE const_mod
