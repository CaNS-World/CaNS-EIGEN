! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_chkdt
  use mpi
  use mod_common_mpi, only:ierr
  use mod_param     , only:is_impdiff,is_impdiff_1d
  use mod_types
  implicit none
  private
  public chkdt
  contains
  subroutine chkdt(n,dxci,dxfi,dyci,dyfi,dzci,dzfi,visc,alpha,u,v,w,dtmax)
    !
    ! computes maximum allowed time step
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: dxci,dxfi,dyci,dyfi,dzci,dzfi
    real(rp), intent(in) :: visc,alpha
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out) :: dtmax
    real(rp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(rp) :: dtix,dtiy,dtiz,dti
    integer :: i,j,k
    real(rp), save :: dlmin
    logical , save :: is_first = .true.
    !
    if(is_first) then ! calculate dlmin only once
      is_first = .false.
      dlmin = min(minval(1./dxfi),minval(1./dyfi))
      if(.not.is_impdiff_1d) then
        dlmin = min(dlmin,minval(1./dzfi))
      end if
      call MPI_ALLREDUCE(MPI_IN_PLACE,dlmin,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    end if
    !
    dti = 0.
    !$acc data copy(dti) async(1)
    !$acc parallel loop collapse(3) default(present) private(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) reduction(max:dti) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) REDUCTION(max:dti)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ux = abs(u(i,j,k))
          vx = 0.25*abs( v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k) )
          wx = 0.25*abs( w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1) )
          dtix = ux*dxci(i)+vx*dyfi(j)+wx*dzfi(k)
          uy = 0.25*abs( u(i,j,k)+u(i,j+1,k)+u(i-1,j+1,k)+u(i-1,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25*abs( w(i,j,k)+w(i,j+1,k)+w(i,j+1,k-1)+w(i,j,k-1) )
          dtiy = uy*dxfi(i)+vy*dyci(j)+wy*dzfi(k)
          uz = 0.25*abs( u(i,j,k)+u(i-1,j,k)+u(i-1,j,k+1)+u(i,j,k+1) )
          vz = 0.25*abs( v(i,j,k)+v(i,j-1,k)+v(i,j-1,k+1)+v(i,j,k+1) )
          wz = abs(w(i,j,k))
          dtiz = uz*dxfi(i)+vz*dyfi(j)+wz*dzci(k)
          dti = max(dti,dtix,dtiy,dtiz)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,dti,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(dti < epsilon(0._rp)) dti = 1.
    if(is_impdiff .and. .not.is_impdiff_1d) then
      dtmax = sqrt(3.)/dti
    else
      dtmax = min(1.65/12./max(visc,alpha)*dlmin**2,sqrt(3.)/dti)
    end if
  end subroutine chkdt
end module mod_chkdt
