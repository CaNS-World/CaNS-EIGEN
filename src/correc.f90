! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_correc
  use mod_types
  implicit none
  private
  public correc
  contains
  subroutine correc(n,dxci,dyci,dzci,dt,p,u,v,w)
    !
    ! corrects the velocity so that it is divergence free
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: dxci,dyci,dzci
    real(rp), intent(in) :: dt
    real(rp), intent(in   ), dimension(0:,0:,0:) :: p
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)
          u(i,j,k) = u(i,j,k) - dt*dxci(i)*(p(i+1,j,k)-p(i,j,k))
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)+1
      do j=0,n(2)
        do i=0,n(1)+1
          v(i,j,k) = v(i,j,k) - dt*dyci(j)*(p(i,j+1,k)-p(i,j,k))
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)
      do j=0,n(2)+1
        do i=0,n(1)+1
          w(i,j,k) = w(i,j,k) - dt*dzci(k)*(p(i,j,k+1)-p(i,j,k))
        end do
      end do
    end do
  end subroutine correc
end module mod_correc
