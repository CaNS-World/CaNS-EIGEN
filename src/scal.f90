! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_scal
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mod_types
  use mod_solve_helmholtz, only: rhs_bound
  implicit none
  private
  public scal,scalar,cmpt_scalflux,bulk_forcing_s,initialize_scalars
  !
  ! scalar derived type
  !
  type scalar
    real(rp), allocatable, dimension(:,:,:) :: val
    real(rp), allocatable, dimension(:,:,:) :: dsdtrko
    real(rp) :: alpha
    character(len=100) :: ini
    character(len=1), dimension(0:1,3) :: cbc
    real(rp),         dimension(0:1,3) :: bc
    real(rp) :: source
    logical  :: is_forced
    real(rp) :: scalf
    real(rp) :: f
    real(rp), dimension(0:1,3) :: fluxo
#if !defined(_OPENACC) || defined(_USE_HIP)
    type(C_PTR), dimension(2,2) :: arrplan
#else
    integer    , dimension(2,2) :: arrplan
#endif
    real(rp), allocatable, dimension(:) :: lambdax_g,lambday_g
    real(rp), allocatable, dimension(:,:) :: lambdaxy,eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd
    real(rp), allocatable, dimension(:) :: a,b,c
    real(rp), dimension(2) :: normfft
    type(rhs_bound) :: rhsb
  end type scalar
  !
  contains
  subroutine scal(nx,ny,nz,dxci,dxfi,dyci,dyfi,dzci,dzfi,visc,u,v,w,s,dsdt,dsdtd)
    use mod_param, only: impdiff_mode,impdiff_explicit,impdiff_z,impdiff_yz,impdiff_xyz
    !
    ! computes convective and diffusive fluxes
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(0:) :: dxci,dxfi,dyci,dyfi,dzci,dzfi
    real(rp), intent(in) :: visc
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w,s
    real(rp), dimension(:,:,:), intent(out) :: dsdt
    real(rp), dimension(:,:,:), intent(out), optional :: dsdtd
    integer :: i,j,k
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: dsdtd_x,dsdtd_y,dsdtd_z,dsdtd_xy
    !
#if !defined(_LOOP_UNSWITCHING)
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_x,dsdtd_y,dsdtd_z) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) &
    !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_x,dsdtd_y,dsdtd_z)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          usim  = 0.5*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
          usip  = 0.5*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
          vsjm  = 0.5*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
          vsjp  = 0.5*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
          wskm  = 0.5*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
          wskp  = 0.5*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
          dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxci(i  )
          dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxci(i-1)
          dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyci(j  )
          dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyci(j-1)
          dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
          dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
          !
          dsdt(i,j,k) = dxfi(i)*( -usip + usim ) + &
                        dyfi(j)*( -vsjp + vsjm ) + &
                        dzfi(k)*( -wskp + wskm )
          dsdtd_x = (dsdxp-dsdxm)*visc*dxfi(i)
          dsdtd_y = (dsdyp-dsdym)*visc*dyfi(j)
          dsdtd_z = (dsdzp-dsdzm)*visc*dzfi(k)
          select case(impdiff_mode)
          case(impdiff_explicit)
            dsdt(i,j,k)  = dsdt(i,j,k) + dsdtd_x + dsdtd_y + dsdtd_z
          case(impdiff_z)
            dsdt(i,j,k)  = dsdt(i,j,k) + dsdtd_x + dsdtd_y
            dsdtd(i,j,k) = dsdtd_z
          case(impdiff_yz)
            dsdt(i,j,k)  = dsdt(i,j,k) + dsdtd_x
            dsdtd(i,j,k) = dsdtd_y + dsdtd_z
          case(impdiff_xyz)
            dsdtd(i,j,k) = dsdtd_x + dsdtd_y + dsdtd_z
          end select
        end do
      end do
    end do
#else
    if(impdiff_mode == impdiff_explicit) then
      !$acc parallel loop collapse(3) default(present) &
      !$acc private(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) &
      !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            usim  = 0.5*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
            usip  = 0.5*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
            vsjm  = 0.5*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
            vsjp  = 0.5*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
            wskm  = 0.5*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
            wskp  = 0.5*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
            dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxci(i  )
            dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxci(i-1)
            dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyci(j  )
            dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyci(j-1)
            dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
            dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
            dsdt(i,j,k) = dxfi(i)*( -usip + usim ) + &
                          dyfi(j)*( -vsjp + vsjm ) + &
                          dzfi(k)*( -wskp + wskm )
            dsdtd_xy = (dsdxp-dsdxm)*visc*dxfi(i) + &
                       (dsdyp-dsdym)*visc*dyfi(j)
            dsdtd_z  = (dsdzp-dsdzm)*visc*dzfi(k)
            dsdt(i,j,k) = dsdt(i,j,k) + dsdtd_xy + dsdtd_z
          end do
        end do
      end do
    else if(impdiff_mode == impdiff_z) then
      !$acc parallel loop collapse(3) default(present) &
      !$acc private(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) &
      !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            usim  = 0.5*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
            usip  = 0.5*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
            vsjm  = 0.5*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
            vsjp  = 0.5*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
            wskm  = 0.5*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
            wskp  = 0.5*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
            dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxci(i  )
            dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxci(i-1)
            dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyci(j  )
            dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyci(j-1)
            dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
            dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
            dsdt(i,j,k) = dxfi(i)*( -usip + usim ) + &
                          dyfi(j)*( -vsjp + vsjm ) + &
                          dzfi(k)*( -wskp + wskm )
            dsdtd_xy = (dsdxp-dsdxm)*visc*dxfi(i) + &
                       (dsdyp-dsdym)*visc*dyfi(j)
            dsdtd_z  = (dsdzp-dsdzm)*visc*dzfi(k)
            dsdt(i,j,k)  = dsdt(i,j,k) + dsdtd_xy
            dsdtd(i,j,k) = dsdtd_z
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(3) default(present) &
      !$acc private(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_x,dsdtd_y,dsdtd_z) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) &
      !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_x,dsdtd_y,dsdtd_z)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            usim  = 0.5*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
            usip  = 0.5*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
            vsjm  = 0.5*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
            vsjp  = 0.5*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
            wskm  = 0.5*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
            wskp  = 0.5*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
            dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxci(i  )
            dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxci(i-1)
            dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyci(j  )
            dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyci(j-1)
            dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
            dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
            dsdt(i,j,k) = dxfi(i)*( -usip + usim ) + &
                          dyfi(j)*( -vsjp + vsjm ) + &
                          dzfi(k)*( -wskp + wskm )
            dsdtd_x = (dsdxp-dsdxm)*visc*dxfi(i)
            dsdtd_y = (dsdyp-dsdym)*visc*dyfi(j)
            dsdtd_z = (dsdzp-dsdzm)*visc*dzfi(k)
            if(impdiff_mode == impdiff_yz) then
              dsdt( i,j,k) = dsdt(i,j,k) + dsdtd_x
              dsdtd(i,j,k) = dsdtd_y + dsdtd_z
            else
              dsdtd(i,j,k) = dsdtd_x + dsdtd_y + dsdtd_z
            end if
          end do
        end do
      end do
    end if
#endif
  end subroutine scal
  !
  subroutine cmpt_scalflux(n,is_bound,l,dxci,dxfi,dyci,dyfi,dzci,dzfi,alpha,s,flux)
    use mpi
    use mod_param, only: cbcpre ! it is fine to use the pressure BC to check for periodicity
    implicit none
    integer , intent(in ), dimension(3) :: n
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l
    real(rp), intent(in ), dimension(0:)    :: dxci,dxfi,dyci,dyfi,dzci,dzfi
    real(rp), intent(in )                   :: alpha
    real(rp), intent(in ), dimension(0:,0:,0:) :: s
    real(rp), intent(out), dimension(0:1,1:3) :: flux
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: flux_xp,flux_xm,flux_yp,flux_ym,flux_zp,flux_zm
    integer :: ierr
    !
    integer :: i,j,k,nx,ny,nz
    real(rp) :: lx,ly,lz
    !
    nx = n(1); ny = n(2); nz = n(3)
    lx = l(1); ly = l(2); lz = l(3)
    flux_xp = 0.
    flux_xm = 0.
    !$acc data copy(flux_xp,flux_xm) async(1)
    if(is_bound(0,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdxp) reduction(+:flux_xp) async(1)
      !$OMP PARALLEL DO   COLLAPSE(2) DEFAULT(shared ) private(dsdxp) reduction(+:flux_xp)
      do k=1,nz
        do j=1,ny
          dsdxp   = (s(1 ,j,k)-s(0   ,j,k))*dxci(0)*alpha
          flux_xp = flux_xp + dsdxp/(dyfi(j)*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_bound(1,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdxm) reduction(+:flux_xm) async(1)
      !$OMP PARALLEL DO   collapse(2) default(shared ) private(dsdxm) reduction(+:flux_xm)
      do k=1,nz
        do j=1,ny
          dsdxm   = (s(nx,j,k)-s(nx+1,j,k))*dxci(nx)*alpha
          flux_xm = flux_xm + dsdxm/(dyfi(j)*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$acc end data
    flux_yp = 0.
    flux_ym = 0.
    !$acc data copy(flux_yp,flux_ym) async(1)
    if(is_bound(0,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdyp) reduction(+:flux_yp) async(1)
      !$OMP PARALLEL DO   collapse(2) default(shared ) private(dsdyp) reduction(+:flux_yp)
      do k=1,nz
        do i=1,nx
          dsdyp   = (s(i,1 ,k)-s(i,0   ,k))*dyci(0)*alpha
          flux_yp = flux_yp + dsdyp/(dxfi(i)*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_bound(1,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdym) reduction(+:flux_ym) async(1)
      !$OMP PARALLEL DO   collapse(2) default(shared ) private(dsdym) reduction(+:flux_ym)
      do k=1,nz
        do i=1,nx
          dsdym   = (s(i,ny,k)-s(i,ny+1,k))*dyci(ny)*alpha
          flux_ym = flux_ym + dsdym/(dxfi(i)*dzfi(k)*lx*lz)
        end do
      end do
    end if
    !$acc end data
    flux_zp = 0.
    flux_zm = 0.
    !$acc data copy(flux_zp,flux_zm) async(1)
    if(is_bound(0,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdzp) reduction(+:flux_zp) async(1)
      !$OMP PARALLEL DO   collapse(2) default(shared ) private(dsdzp) reduction(+:flux_zp)
      do j=1,ny
        do i=1,nx
          dsdzp   = (s(i,j,1 )-s(i,j,0   ))*dzci(0)*alpha
          flux_zp = flux_zp + dsdzp/(dxfi(i)*dyfi(j)*lx*ly)
        end do
      end do
    end if
    if(is_bound(1,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdzm) reduction(+:flux_zm) async(1)
      !$OMP PARALLEL DO   collapse(2) default(shared ) private(dsdzm) reduction(+:flux_zm)
      do j=1,ny
        do i=1,nx
          dsdzm   = (s(i,j,nz)-s(i,j,nz+1))*dzci(nz)*alpha
          flux_zm = flux_zm + dsdzm/(dxfi(i)*dyfi(j)*lx*ly)
        end do
      end do
    end if
    !$acc end data
    !$acc wait(1)
    !
    !flux(:) = [flux_x,flux_y,flux_z]
    !
    ! discern upper and lower boundary fluxes
    !
    flux(:,1) = [flux_xp,flux_xm]
    flux(:,2) = [flux_yp,flux_ym]
    flux(:,3) = [flux_zp,flux_zm]
    !
    call MPI_ALLREDUCE(MPI_IN_PLACE,flux,6,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine cmpt_scalflux
  !
  subroutine bulk_forcing_s(n,is_forced,ff,p)
    implicit none
    integer , intent(in   ), dimension(3) :: n
    logical , intent(in   )               :: is_forced
    real(rp), intent(in   )               :: ff
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer :: i,j,k
    if(is_forced) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,k) = p(i,j,k) + ff
          end do
        end do
      end do
    end if
  end subroutine bulk_forcing_s
  !
  subroutine initialize_scalars(scalars,nscal,ng,n,n_z)
    use mod_param, only: alphai,iniscal,cbcscal,bcscal,ssource,is_sforced,scalf, &
                         impdiff_mode,impdiff_explicit
    !
    ! initializes/allocates members of an array of `nscal` scalar derived types
    !
    implicit none
    type(scalar), intent(inout), dimension(:) :: scalars
    integer     , intent(in   ) :: nscal,ng(3),n(3),n_z(3)
    integer :: iscal
    do iscal=1,nscal
      allocate(scalars(iscal)%val(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      allocate(scalars(iscal)%dsdtrko(n(1),n(2),n(3)))
      if(impdiff_mode /= impdiff_explicit) then
        allocate(scalars(iscal)%lambdax_g(ng(1)),scalars(iscal)%lambday_g(ng(2)))
        allocate(scalars(iscal)%lambdaxy(n_z(1),n_z(2)))
        allocate(scalars(iscal)%eigvecx_fwd(ng(1),ng(1)),scalars(iscal)%eigvecx_bwd(ng(1),ng(1)), &
                 scalars(iscal)%eigvecy_fwd(ng(2),ng(2)),scalars(iscal)%eigvecy_bwd(ng(2),ng(2)))
        allocate(scalars(iscal)%a(n_z(3)), &
                 scalars(iscal)%b(n_z(3)), &
                 scalars(iscal)%c(n_z(3)))
        allocate(scalars(iscal)%rhsb%x(n(2),n(3),0:1), &
                 scalars(iscal)%rhsb%y(n(1),n(3),0:1), &
                 scalars(iscal)%rhsb%z(n(1),n(2),0:1))
      end if
      scalars(iscal)%alpha     = alphai(iscal)**(-1)
      scalars(iscal)%ini       = iniscal(iscal)
      scalars(iscal)%cbc       = cbcscal(:,:,iscal)
      scalars(iscal)%bc        = bcscal(:,:,iscal)
      scalars(iscal)%source    = ssource(iscal)
      scalars(iscal)%is_forced = is_sforced(iscal)
      scalars(iscal)%scalf     = scalf(iscal)
    end do
  end subroutine initialize_scalars
end module mod_scal
