! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_fillps
  use mod_types
  implicit none
  private
  public fillps
  contains
  subroutine fillps(n,dxfi,dyfi,dzfi,dti,u,v,w,p)
    !
    !  fill the right-hand side of the Poisson equation for the correction pressure.
    !
    !  the discrete divergence is:
    !
    !  w(i,j,k)-w(i,j,k-1)   v(i,j,k)-v(i,j-1,k)   u(i,j,k)-u(i-1,j,k)
    !  ------------------- + ------------------- + -------------------  = div
    !          dz                    dy                    dx
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(0:) :: dxfi,dyfi,dzfi
    real(rp), intent(in ) :: dti
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(0:,0:,0:) :: p
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          p(i,j,k) = ( &
                      (w(i,j,k)-w(i,j,k-1))*dti*dzfi(k) + &
                      (v(i,j,k)-v(i,j-1,k))*dti*dyfi(j) + &
                      (u(i,j,k)-u(i-1,j,k))*dti*dxfi(i) &
                     )
        end do
      end do
    end do
  end subroutine fillps
end module mod_fillps
