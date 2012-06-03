! NCLFORTSTART
subroutine dphybrid (p0,hyai,hybi,psfc,mlon,nlat,klev,dphy)
  implicit none

  ! input
  integer::mlon, nlat, klev
  real   ::p0
  real,dimension(klev)::hyai,hybi
  real,dimension(mlon,nlat)::psfc
  real,dimension(mlon,nlat,klev-1)::dphy
  ! NCLEND

  ! NCL: dphy = dpres_hybrid_ccm (ps,p0,hyai,hybi)
  ! Note: dphy will have one less vertical level

  ! this routine interploates ccm2/3 hybrid coordinate data
  !     to pressure coordinates. then it computes

  !     formula for the pressure of a hybrid surface is;
  !          phy(k) = hyai(k)*p0 + hybi(k)*psfc

  !     input
  !          hyai   - is the "a" or pressure hybrid coef
  !          hybi   - is the "b" or sigma coeficient
  !          p0     - is the base pressure in pascals [Pa]
  !          psfc   - is the surface pressure in pascals [Pa]
  !          mlon   - longitude dimension
  !          nlat   - latitude  dimension
  !          klev   - number of levels
  !     output
  !          dphy    - delta pressure [Pa; always positive]

  ! local
  integer:: kl, nl, ml
  real   :: pa, pb

  do nl=1,nlat
     do ml=1,mlon
        do kl=1,klev-1
           pa = p0*hyai(kl)   + hybi(kl  )*psfc(ml,nl)
           pb = p0*hyai(kl+1) + hybi(kl+1)*psfc(ml,nl)
           dphy(ml,nl,kl) = abs(pb-pa)
        end do
     end do
  end do

  return
end subroutine dphybrid
