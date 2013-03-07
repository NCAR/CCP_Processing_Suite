subroutine pres_hybrid_ccsm(ps_in,ps_out,nx,ny,nz)
  use grid_info
  implicit none
  !
  integer::nx,ny,nz
  real,dimension(nx,ny)::ps_in,ps_out
  !
  integer::i,j,k
  !
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx
           ps_out(i,j) = a_coeff_bnds(k)*p0 + b_coeff_bnds(k)*ps_in(i,j)
        enddo
     enddo
  enddo
end subroutine pres_hybrid_ccsm
