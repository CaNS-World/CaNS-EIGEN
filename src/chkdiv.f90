! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_chkdiv
  use mpi
  use mod_common_mpi, only: ierr
  use mod_types
  implicit none
  private
  public chkdiv
  contains
  subroutine chkdiv(lo,hi,l,dxfi,dyfi,dzfi,u,v,w,divtot,divmax)
    !
    ! checks the divergence of the velocity field
    !
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(in), dimension(lo(1)-1:) :: dxfi
    real(rp), intent(in), dimension(lo(2)-1:) :: dyfi
    real(rp), intent(in), dimension(lo(3)-1:) :: dzfi
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(rp), intent(out) :: divtot,divmax
    real(rp) :: div
    integer :: i,j,k

    divtot = 0.
    divmax = 0.
    !$acc data copy(divtot,divmax) async(1)
    !$acc parallel loop collapse(3) default(present) private(div) reduction(+:divtot) reduction(max:divmax) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(div) REDUCTION(+:divtot) REDUCTION(max:divmax)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          div = (w(i,j,k)-w(i,j,k-1))*dzfi(k) + &
                (v(i,j,k)-v(i,j-1,k))*dyfi(j) + &
                (u(i,j,k)-u(i-1,j,k))*dxfi(i)
          divmax = max(divmax,abs(div))
          divtot = divtot + abs(div)/(dxfi(i)*dyfi(j)*dzfi(k)) ! abs(div)*cell_volume
          !if(abs(div) >= 1.e-12) print*,div,'Large divergence at grid cell: ',i,j,k,div
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,divtot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,divmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    divtot = divtot/product(l(:))
  end subroutine chkdiv
end module mod_chkdiv
