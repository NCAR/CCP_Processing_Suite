subroutine pres_hybrid_ccm(ps_in,ps_out,nx,ny,nz)
  use grid_info
  implicit none
  !
  integer::nx,ny,nz
  real::pzero
  real,dimension(nx,ny)::ps_in
  real,dimension(nx,ny,nz)::ps_out
  !
  integer::i,j,k
  !
  pzero = p0 * 100.
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx
           ps_out(i,j,k) = a_coeff_bnds(k)*pzero + b_coeff_bnds(k)*ps_in(i,j)
        enddo
     enddo
  enddo
!!$  write(*,'(''pres_hybrid_ccsm PS: '',2f16.4)') minval(ps_in),maxval(ps_in)
!!$  do k = 1,nz
!!$     write(*,'(''pres_hybrid_ccsm PI'',i2.2,'' : '',2f16.4)') k,minval(ps_out(:,:,k)),maxval(ps_out(:,:,k))
!!$  enddo
end subroutine pres_hybrid_ccm
