! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_debug
  use mpi
  use mod_common_mpi, only: ierr
  use mod_types
  implicit none
  private
  public chk_helmholtz
  contains
  !
  subroutine chk_helmholtz(lo,hi,l,dxci,dxfi,dyci,dyfi,dzci,dzfi,alphai,fp,fpp,bc,is_bound,c_or_f,difftot,diffmax)
    !
    ! this subroutine checks the correctness of the solution of a Helmholtz equation
    ! with collocated or staggered boundary conditions
    !
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(in), dimension(lo(1)-1:) :: dxci,dxfi
    real(rp), intent(in), dimension(lo(2)-1:) :: dyci,dyfi
    real(rp), intent(in), dimension(lo(3)-1:) :: dzci,dzfi
    real(rp), intent(in) :: alphai
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: fp,fpp
    character(len=1), intent(in), dimension(0:1,3) :: bc
    character(len=1), intent(in), dimension(3) :: c_or_f
    logical         , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(out) :: difftot,diffmax
    real(rp) :: val,res
    integer :: i,j,k
    integer :: idir
    integer, dimension(3) :: q
    q(:) = 0
    do idir = 1,3
      if(bc(1,idir) /= 'P'.and.c_or_f(idir) == 'f'.and.is_bound(1,idir)) q(idir) = 1
    end do
    difftot = 0.
    diffmax = 0.
    !$acc wait
    !$acc data copy(difftot,diffmax,q)
    select case(c_or_f(1)//c_or_f(2)//c_or_f(3))
    case('ccc')
      !$acc parallel loop collapse(3) default(present) private(val,res) reduction(+:difftot) reduction(max:diffmax)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared ) PRIVATE(val,res) REDUCTION(+:difftot) REDUCTION(max:diffmax)
      do k=lo(3),hi(3)-q(3)
        do j=lo(2),hi(2)-q(2)
          do i=lo(1),hi(1)-q(1)
            val = alphai*fpp(i,j,k) + &
                  ((fpp(i+1,j,k)-fpp(i  ,j,k))*dxci(i  ) - &
                   (fpp(i  ,j,k)-fpp(i-1,j,k))*dxci(i-1))*dxfi(i) + &
                  ((fpp(i,j+1,k)-fpp(i,j  ,k))*dyci(j  ) - &
                   (fpp(i,j  ,k)-fpp(i,j-1,k))*dyci(j-1))*dyfi(j) + &
                  ((fpp(i,j,k+1)-fpp(i,j,k  ))*dzci(k  ) - &
                   (fpp(i,j,k  )-fpp(i,j,k-1))*dzci(k-1))*dzfi(k)
            res = abs(val-fp(i,j,k)*alphai)
            difftot = difftot + res/(dxfi(i)*dyfi(j)*dzfi(k)) ! abs(res)*cell_volume
            diffmax = max(diffmax,res)
          end do
        end do
      end do
    case('fcc')
      !$acc parallel loop collapse(3) default(present) private(val,res) reduction(+:difftot) reduction(max:diffmax)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared ) PRIVATE(val,res) REDUCTION(+:difftot) REDUCTION(max:diffmax)
      do k=lo(3),hi(3)-q(3)
        do j=lo(2),hi(2)-q(2)
          do i=lo(1),hi(1)-q(1)
            val = alphai*fpp(i,j,k) + &
                  ((fpp(i+1,j,k)-fpp(i  ,j,k))*dxfi(i+1) - &
                   (fpp(i  ,j,k)-fpp(i-1,j,k))*dxfi(i  ))*dxci(i) + &
                  ((fpp(i,j+1,k)-fpp(i,j  ,k))*dyci(j  ) - &
                   (fpp(i,j  ,k)-fpp(i,j-1,k))*dyci(j-1))*dyfi(j) + &
                  ((fpp(i,j,k+1)-fpp(i,j,k  ))*dzci(k  ) - &
                   (fpp(i,j,k  )-fpp(i,j,k-1))*dzci(k-1))*dzfi(k)
            res = abs(val-fp(i,j,k)*alphai)
            difftot = difftot + res/(dxci(i)*dyfi(j)*dzfi(k)) ! abs(res)*cell_volume
            diffmax = max(diffmax,res)
          end do
        end do
      end do
    case('cfc')
      !$acc parallel loop collapse(3) default(present) private(val,res) reduction(+:difftot) reduction(max:diffmax)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared ) PRIVATE(val,res) REDUCTION(+:difftot) REDUCTION(max:diffmax)
      do k=lo(3),hi(3)-q(3)
        do j=lo(2),hi(2)-q(2)
          do i=lo(1),hi(1)-q(1)
            val = alphai*fpp(i,j,k) + &
                  ((fpp(i+1,j,k)-fpp(i  ,j,k))*dxci(i  ) - &
                   (fpp(i  ,j,k)-fpp(i-1,j,k))*dxci(i-1))*dxfi(i) + &
                  ((fpp(i,j+1,k)-fpp(i,j  ,k))*dyfi(j+1) - &
                   (fpp(i,j  ,k)-fpp(i,j-1,k))*dyfi(j  ))*dyci(j) + &
                  ((fpp(i,j,k+1)-fpp(i,j,k  ))*dzci(k  ) - &
                   (fpp(i,j,k  )-fpp(i,j,k-1))*dzci(k-1))*dzfi(k)
            res = abs(val-fp(i,j,k)*alphai)
            difftot = difftot + res/(dxfi(i)*dyci(j)*dzfi(k)) ! abs(res)*cell_volume
            diffmax = max(diffmax,res)
          end do
        end do
      end do
    case('ccf')
      !$acc parallel loop collapse(3) default(present) private(val,res) reduction(+:difftot) reduction(max:diffmax)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared ) PRIVATE(val,res) REDUCTION(+:difftot) REDUCTION(max:diffmax)
      do k=lo(3),hi(3)-q(3)
        do j=lo(2),hi(2)-q(2)
          do i=lo(1),hi(1)-q(1)
            val = alphai*fpp(i,j,k) + &
                  ((fpp(i+1,j,k)-fpp(i  ,j,k))*dxci(i  ) - &
                   (fpp(i  ,j,k)-fpp(i-1,j,k))*dxci(i-1))*dxfi(i) + &
                  ((fpp(i,j+1,k)-fpp(i,j  ,k))*dyci(j  ) - &
                   (fpp(i,j  ,k)-fpp(i,j-1,k))*dyci(j-1))*dyfi(j) + &
                  ((fpp(i,j,k+1)-fpp(i,j,k  ))*dzfi(k+1) - &
                   (fpp(i,j,k  )-fpp(i,j,k-1))*dzfi(k  ))*dzci(k)
            res = abs(val-fp(i,j,k)*alphai)
            difftot = difftot + res/(dxfi(i)*dyfi(j)*dzci(k)) ! abs(res)*cell_volume
            diffmax = max(diffmax,res)
          end do
        end do
      end do
    end select
    !$acc end data
    call MPI_ALLREDUCE(MPI_IN_PLACE,difftot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,diffmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    difftot = difftot/product(l(:))
  end subroutine chk_helmholtz
  !
  subroutine chk_poisson(lo,hi,l,dxci,dxfi,dyci,dyfi,dzci,dzfi,fp,fpp,difftot,diffmax)
    !
    ! this subroutine checks if the Poisson equation is correctly solved;
    ! the inputs should be downcasted to precision `gp` to ensure proper
    ! testing when single precision is used
    !
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(in), dimension(lo(1)-1:) :: dxci,dxfi
    real(rp), intent(in), dimension(lo(2)-1:) :: dyci,dyfi
    real(rp), intent(in), dimension(lo(3)-1:) :: dzci,dzfi
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: fp,fpp
    real(rp), intent(out) :: difftot,diffmax
    real(rp) :: val,res
    integer :: i,j,k
    difftot = 0.
    diffmax = 0.
    !$acc wait
    !$acc data copy(difftot,diffmax)
    !$acc parallel loop collapse(3) default(present) private(val,res) reduction(+:difftot) reduction(max:diffmax)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared ) PRIVATE(val,res) REDUCTION(+:difftot) REDUCTION(max:diffmax)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          val = ((fpp(i+1,j,k)-fpp(i  ,j,k))*dxci(i  ) - &
                 (fpp(i  ,j,k)-fpp(i-1,j,k))*dxci(i-1))*dxfi(i) + &
                ((fpp(i,j+1,k)-fpp(i,j  ,k))*dyci(j  ) - &
                 (fpp(i,j  ,k)-fpp(i,j-1,k))*dyci(j-1))*dyfi(j) + &
                ((fpp(i,j,k+1)-fpp(i,j,k  ))*dzci(k  ) - &
                 (fpp(i,j,k  )-fpp(i,j,k-1))*dzci(k-1))*dzfi(k)
          res = abs(val-fp(i,j,k))
          difftot = difftot + res/(dxfi(i)*dyfi(j)*dzfi(k)) ! abs(res)*cell_volume
          diffmax = max(diffmax,res)
          !if(abs(val-fp(i,j,k)) > 1.e-8) print*, 'Large difference : ', val-fp(i,j,k),i,j,k
        end do
      end do
    end do
    !$acc end data
    call MPI_ALLREDUCE(MPI_IN_PLACE,difftot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,diffmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    difftot = difftot/product(l(:))
  end subroutine chk_poisson
end module mod_debug
