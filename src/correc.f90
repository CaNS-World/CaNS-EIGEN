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
  subroutine correc(n,nh_v,nh_p,nh_g,dxci,dyci,dzci,dt,p,u,v,w)
    !
    ! corrects the velocity so that it is divergence free
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh_v,nh_p,nh_g
    real(rp), intent(in), dimension(1-nh_g:) :: dxci,dyci,dzci
    real(rp), intent(in) :: dt
    real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    integer :: i,j,k,m
    !
    ! the velocity is also corrected in the halo layers covered by `p`, since the
    ! predicted velocity halos are updated before; only the upper halo layers along
    ! the normal direction are left for the halo exchange that follows
    !
    m = min(nh_v,nh_p)
    !
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1-m,n(3)+m
      do j=1-m,n(2)+m
        do i=1-m,n(1)
          u(i,j,k) = u(i,j,k) - dt*dxci(i)*(p(i+1,j,k)-p(i,j,k))
        end do
      end do
    end do
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1-m,n(3)+m
      do j=1-m,n(2)
        do i=1-m,n(1)+m
          v(i,j,k) = v(i,j,k) - dt*dyci(j)*(p(i,j+1,k)-p(i,j,k))
        end do
      end do
    end do
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1-m,n(3)
      do j=1-m,n(2)+m
        do i=1-m,n(1)+m
          w(i,j,k) = w(i,j,k) - dt*dzci(k)*(p(i,j,k+1)-p(i,j,k))
        end do
      end do
    end do
  end subroutine correc
end module mod_correc
