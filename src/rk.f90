! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_rk
  use mod_mom  , only: momx_a => momx_a_vv, &
                       momy_a => momy_a_vv, &
                       momz_a => momz_a_vv, &
                       momx_d,momy_d,momz_d, &
                       momx_p,momy_p,momz_p, &
                       cmpt_wallshear,bulk_forcing, &
                       momx_d_x ,momy_d_x ,momz_d_x, &
                       momx_d_yz,momy_d_yz,momz_d_yz, &
                       momx_d_xy,momy_d_xy,momz_d_xy, &
                       momx_d_z ,momy_d_z ,momz_d_z, &
                       mom_xyz_ad
  use mod_param, only: impdiff_mode,impdiff_explicit,impdiff_z,impdiff_yz,impdiff_xyz, &
                       is_boussinesq_buoyancy,is_fast_mom_kernels
  use mod_scal , only: scal,cmpt_scalflux,scalar
  use mod_utils, only: bulk_mean
  use mod_types
  implicit none
  public rk,rk_scal,bulk_forcing_wall_shear
  contains
  subroutine rk(rkpar,n,nh_v,nh_p,nh_g,l,dxci,dxfi,dyci,dyfi,dzci,dzfi,dxc,dxf,dyc,dyf,dzc,dzf,dt,visc,p, &
                is_forced,velf,bforce,gacc,beta,scalars,dudtrko,dvdtrko,dwdtrko,u,v,w,f)
#if defined(_OPENACC)
    use mod_common_cudecomp, only: dudtrk_t => work, &
                                   dvdtrk_t => solver_buf_0, &
                                   dwdtrk_t => solver_buf_1
#endif
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the momentum equations.
    !
    implicit none
    real(rp), intent(in   ), dimension(2)        :: rkpar
    integer , intent(in   ), dimension(3)        :: n
    integer , intent(in   )                      :: nh_v,nh_p,nh_g
    real(rp), intent(in   ), dimension(3)        :: l
    real(rp), intent(in   ), dimension(1-nh_g:)  :: dxci,dxfi,dyci,dyfi,dzci,dzfi
    real(rp), intent(in   ), dimension(1-nh_g:)  :: dxc,dxf,dyc,dyf,dzc,dzf
    real(rp), intent(in   )                      :: dt,visc
    real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    logical , intent(in   ), dimension(3)        :: is_forced
    real(rp), intent(in   ), dimension(3)        :: velf,bforce
    real(rp), intent(in   ), dimension(3)        :: gacc
    real(rp), intent(in   )                      :: beta
    type(scalar), intent(in   ), target, dimension(:), optional :: scalars
    real(rp), intent(inout), dimension(1:,1:,1:) :: dudtrko,dvdtrko,dwdtrko
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(out  ), dimension(3)        :: f
    real(rp), pointer, contiguous, dimension(:,:,:) :: s
#if !defined(_OPENACC)
    real(rp), target       , allocatable, dimension(:,:,:), save :: dudtrk_t,dvdtrk_t,dwdtrk_t
#endif
    real(rp), pointer      , contiguous , dimension(:,:,:), save :: dudtrk  ,dvdtrk  ,dwdtrk
    real(rp),                allocatable, dimension(:,:,:), save :: dudtrkd ,dvdtrkd ,dwdtrkd
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12,forcing_u,forcing_v,forcing_w
    logical  :: is_buoyancy
    integer  :: i,j,k
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !
    ! is_boussinesq_buoyancy = T will always imply present(scalars), as we do not allow for
    ! nscal = 0 with buoyancy on (see `param.f90`)
    !
    is_buoyancy = present(scalars).and.is_boussinesq_buoyancy
    if(is_buoyancy) then
      s => scalars(1)%val
    end if
    !
    ! initialization
    !
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
#if !defined(_OPENACC)
      allocate(dudtrk_t(n(1),n(2),n(3)),dvdtrk_t(n(1),n(2),n(3)),dwdtrk_t(n(1),n(2),n(3)))
#endif
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dudtrko(i,j,k) = 0._rp
            dvdtrko(i,j,k) = 0._rp
            dwdtrko(i,j,k) = 0._rp
          end do
        end do
      end do
      if(impdiff_mode /= impdiff_explicit) then
        allocate(dudtrkd(n(1),n(2),n(3)),dvdtrkd(n(1),n(2),n(3)),dwdtrkd(n(1),n(2),n(3)))
        !$acc enter data create(dudtrkd,dvdtrkd,dwdtrkd) async(1)
      end if
#if defined(_OPENACC)
      dudtrk(1:n(1),1:n(2),1:n(3)) => dudtrk_t(1:product(int(n(:),i8)))
      dvdtrk(1:n(1),1:n(2),1:n(3)) => dvdtrk_t(1:product(int(n(:),i8)))
      dwdtrk(1:n(1),1:n(2),1:n(3)) => dwdtrk_t(1:product(int(n(:),i8)))
#else
      dudtrk  => dudtrk_t
      dvdtrk  => dvdtrk_t
      dwdtrk  => dwdtrk_t
#endif
    end if
    !
    if(is_fast_mom_kernels) then
      call mom_xyz_ad(n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,dyci,dyfi,dzci,dzfi,visc,u,v,w,dudtrk,dvdtrk,dwdtrk, &
                      dudtrkd,dvdtrkd,dwdtrkd)
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dudtrk(i,j,k) = 0._rp
            dvdtrk(i,j,k) = 0._rp
            dwdtrk(i,j,k) = 0._rp
          end do
        end do
      end do
      if(impdiff_mode == impdiff_explicit) then
        call momx_d(n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,dyci,dyfi,dzci,dzfi,visc,u,dudtrk)
        call momy_d(n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,dyci,dyfi,dzci,dzfi,visc,v,dvdtrk)
        call momz_d(n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,dyci,dyfi,dzci,dzfi,visc,w,dwdtrk)
      else
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              dudtrkd(i,j,k) = 0._rp
              dvdtrkd(i,j,k) = 0._rp
              dwdtrkd(i,j,k) = 0._rp
            end do
          end do
        end do
        select case(impdiff_mode)
        case(impdiff_z)
          call momx_d_xy(n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,dyci,dyfi,visc,u,dudtrk)
          call momy_d_xy(n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,dyci,dyfi,visc,v,dvdtrk)
          call momz_d_xy(n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,dyci,dyfi,visc,w,dwdtrk)
          call momx_d_z( n(1),n(2),n(3),nh_v,nh_g,dzci,dzfi,visc,u,dudtrkd)
          call momy_d_z( n(1),n(2),n(3),nh_v,nh_g,dzci,dzfi,visc,v,dvdtrkd)
          call momz_d_z( n(1),n(2),n(3),nh_v,nh_g,dzci,dzfi,visc,w,dwdtrkd)
        case(impdiff_yz)
          call momx_d_x( n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,visc,u,dudtrk)
          call momy_d_x( n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,visc,v,dvdtrk)
          call momz_d_x( n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,visc,w,dwdtrk)
          call momx_d_yz(n(1),n(2),n(3),nh_v,nh_g,dyci,dyfi,dzci,dzfi,visc,u,dudtrkd)
          call momy_d_yz(n(1),n(2),n(3),nh_v,nh_g,dyci,dyfi,dzci,dzfi,visc,v,dvdtrkd)
          call momz_d_yz(n(1),n(2),n(3),nh_v,nh_g,dyci,dyfi,dzci,dzfi,visc,w,dwdtrkd)
        case(impdiff_xyz)
          call momx_d(n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,dyci,dyfi,dzci,dzfi,visc,u,dudtrkd)
          call momy_d(n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,dyci,dyfi,dzci,dzfi,visc,v,dvdtrkd)
          call momz_d(n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,dyci,dyfi,dzci,dzfi,visc,w,dwdtrkd)
        end select
      end if
      call momx_a(n(1),n(2),n(3),nh_v,nh_g,dxci,dxfi,dyfi,dzfi,u,v,w,dudtrk)
      call momy_a(n(1),n(2),n(3),nh_v,nh_g,dxfi,dyci,dyfi,dzfi,u,v,w,dvdtrk)
      call momz_a(n(1),n(2),n(3),nh_v,nh_g,dxfi,dyfi,dzci,dzfi,u,v,w,dwdtrk)
    end if
    !
#if !defined(_LOOP_UNSWITCHING)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                factor12*(bforce(1) - dxci(i)*( p(i+1,j,k)-p(i,j,k)))
          if(is_buoyancy) then
            u(i,j,k) = u(i,j,k) - factor12*gacc(1)*beta*0.5*(s(i+1,j,k)+s(i,j,k))
          end if
          !
          v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                factor12*(bforce(2) - dyci(j)*( p(i,j+1,k)-p(i,j,k)))
          if(is_buoyancy) then
            v(i,j,k) = v(i,j,k) - factor12*gacc(2)*beta*0.5*(s(i,j+1,k)+s(i,j,k))
          end if
          !
          w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)))
          if(is_buoyancy) then
            w(i,j,k) = w(i,j,k) - factor12*gacc(3)*beta*0.5*(s(i,j,k+1)+s(i,j,k))
          end if
          !
          if(impdiff_mode /= impdiff_explicit) then
            u(i,j,k) = u(i,j,k) + factor12*dudtrkd(i,j,k)
            v(i,j,k) = v(i,j,k) + factor12*dvdtrkd(i,j,k)
            w(i,j,k) = w(i,j,k) + factor12*dwdtrkd(i,j,k)
          end if
        end do
      end do
    end do
#else
    if(impdiff_mode == impdiff_explicit .and. .not.is_buoyancy) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dxci(i)*(p(i+1,j,k)-p(i,j,k)))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dyci(j)*(p(i,j+1,k)-p(i,j,k)))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)))
          end do
        end do
      end do
    else if(impdiff_mode /= impdiff_explicit .and. .not.is_buoyancy) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dxci(i)*(p(i+1,j,k)-p(i,j,k)) + &
                                            dudtrkd(i,j,k))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dyci(j)*(p(i,j+1,k)-p(i,j,k)) + &
                                            dvdtrkd(i,j,k))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)) + &
                                            dwdtrkd(i,j,k))
          end do
        end do
      end do
    else if(impdiff_mode == impdiff_explicit .and. is_buoyancy) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dxci(i)*(p(i+1,j,k)-p(i,j,k)) &
                                            - gacc(1)*beta*0.5*(s(i+1,j,k)+s(i,j,k)))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dyci(j)*(p(i,j+1,k)-p(i,j,k)) &
                                            - gacc(2)*beta*0.5*(s(i,j+1,k)+s(i,j,k)))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)) &
                                            - gacc(3)*beta*0.5*(s(i,j,k+1)+s(i,j,k)))
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dxci(i)*(p(i+1,j,k)-p(i,j,k)) &
                                            - gacc(1)*beta*0.5*(s(i+1,j,k)+s(i,j,k)) + &
                                            dudtrkd(i,j,k))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dyci(j)*(p(i,j+1,k)-p(i,j,k)) &
                                            - gacc(2)*beta*0.5*(s(i,j+1,k)+s(i,j,k)) + &
                                            dvdtrkd(i,j,k))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)) &
                                            - gacc(3)*beta*0.5*(s(i,j,k+1)+s(i,j,k)) + &
                                            dwdtrkd(i,j,k))
          end do
        end do
      end do
    end if
#endif
    !
    ! replaced previous pointer swap to save memory on GPUs by using already allocated
    ! buffers
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
        end do
      end do
    end do
!#if 0 /*pressure gradient term treated explicitly above */
!    !$acc parallel loop collapse(3) default(present) async(1)
!    !$OMP parallel do   collapse(3) DEFAULT(shared)
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!          dudtrk(i,j,k) = 0._rp
!          dvdtrk(i,j,k) = 0._rp
!          dwdtrk(i,j,k) = 0._rp
!        end do
!      end do
!    end do
!    call momx_p(n(1),n(2),n(3),dxci,bforce(1),p,dudtrk)
!    call momy_p(n(1),n(2),n(3),dyci,bforce(2),p,dvdtrk)
!    call momz_p(n(1),n(2),n(3),dzci,bforce(3),p,dwdtrk)
!    !$acc parallel loop collapse(3)
!    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!          u(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
!          v(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
!          w(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
!        end do
!      end do
!    end do
!#endif
    !
    ! compute bulk velocity forcing
    !
    call cmpt_bulk_forcing(n,nh_v,nh_g,l,is_forced,velf,dxc,dxf,dyc,dyf,dzc,dzf,u,v,w,f)
    forcing_u = f(1)
    forcing_v = f(2)
    forcing_w = f(3)
    !
    if(impdiff_mode /= impdiff_explicit) then
      !
      ! apply bulk forcing and compute rhs of the Helmholtz equation
      !
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = (u(i,j,k) - .5_rp*factor12*dudtrkd(i,j,k)) + forcing_u
            v(i,j,k) = (v(i,j,k) - .5_rp*factor12*dvdtrkd(i,j,k)) + forcing_v
            w(i,j,k) = (w(i,j,k) - .5_rp*factor12*dwdtrkd(i,j,k)) + forcing_w
          end do
        end do
      end do
    else if(any(is_forced(:))) then
      !
      ! apply bulk forcing
      !
      call bulk_forcing(n,nh_v,is_forced,f,u,v,w)
    end if
  end subroutine rk
  !
  subroutine rk_scal(rkpar,n,nh_v,nh_s,nh_g,l,dxci,dxfi,dyci,dyfi,dzci,dzfi,dxf,dyf,dzf,alpha,dt,is_bound,u,v,w, &
                     is_forced,scalf,ssource,fluxo,dsdtrko,s,f)
#if defined(_OPENACC)
    use mod_common_cudecomp, only: dsdtrk_t => work
#endif
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the scalar field.
    !
    implicit none
    logical , parameter :: is_cmpt_wallflux = .false.
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_s,nh_g
    real(rp), intent(in   ), dimension(3) :: l
    real(rp), intent(in   ), dimension(1-nh_g:) :: dxci,dxfi,dyci,dyfi,dzci,dzfi
    real(rp), intent(in   ), dimension(1-nh_g:) :: dxf,dyf,dzf
    real(rp), intent(in   ) :: alpha,dt
    logical , intent(in   ), dimension(0:1,3)    :: is_bound
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    logical , intent(in   ) :: is_forced
    real(rp), intent(in   ) :: scalf,ssource
    real(rp), intent(inout), dimension(0:1,3) :: fluxo
    real(rp), intent(inout), dimension(1:,1:,1:) :: dsdtrko
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(out  ) :: f
    !
#if !defined(_OPENACC)
    real(rp), target       , allocatable, dimension(:,:,:), save :: dsdtrk_t
#endif
    real(rp), pointer      , contiguous , dimension(:,:,:), save :: dsdtrk
    real(rp), target       , allocatable, dimension(:,:,:), save :: dsdtrkd
    logical, save :: is_first = .true.
    !
    real(rp) :: factor1,factor2,factor12
    real(rp), dimension(0:1,3) :: flux
    integer :: i,j,k
    real(rp) :: mean
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
#if !defined(_OPENACC)
      allocate(dsdtrk_t(1:n(1),1:n(2),1:n(3)))
#endif
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dsdtrko(i,j,k) = 0._rp
          end do
        end do
      end do
      if(impdiff_mode /= impdiff_explicit) then
        allocate(dsdtrkd(n(1),n(2),n(3)))
        !$acc enter data create(dsdtrkd) async(1)
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              dsdtrkd(i,j,k) = 0._rp
            end do
          end do
        end do
      end if
    end if
#if defined(_OPENACC)
    dsdtrk(1:n(1),1:n(2),1:n(3)) => dsdtrk_t(1:product(int(n(:),i8)))
#else
    dsdtrk => dsdtrk_t
#endif
    !
    call scal(n(1),n(2),n(3),nh_v,nh_s,nh_g,dxci,dxfi,dyci,dyfi,dzci,dzfi,alpha,u,v,w,s,dsdtrk,dsdtrkd)
#if !defined(_LOOP_UNSWITCHING)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource
          if(impdiff_mode /= impdiff_explicit) then
            s(i,j,k) = s(i,j,k) + factor12*dsdtrkd(i,j,k)
          end if
        end do
      end do
    end do
#else
    if(impdiff_mode == impdiff_explicit) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + &
                                  factor12*(ssource + dsdtrkd(i,j,k))
          end do
        end do
      end do
    end if
#endif
    !
    ! compute wall scalar flux
    !
    if(is_cmpt_wallflux) then
      call cmpt_scalflux(n,nh_s,nh_g,is_bound,l,dxci,dxfi,dyci,dyfi,dzci,dzfi,alpha,s(:,:,:),flux)
      f = (factor1*sum((flux( 0,:)+flux( 1,:))/l(:)) + &
           factor2*sum((fluxo(0,:)+fluxo(1,:))/l(:)))
      fluxo(:,:) = flux(:,:)
    end if
    !
    ! bulk scalar forcing
    !
    if(is_forced) then
      call bulk_mean(n,nh_s,nh_g,l,dxf,dyf,dzf,s(:,:,:),mean)
      f = scalf - mean
    end if
    if(impdiff_mode /= impdiff_explicit) then
      !
      ! compute rhs of Helmholtz equation
      !
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            s(i,j,k) = s(i,j,k) - .5_rp*factor12*dsdtrkd(i,j,k)
          end do
        end do
      end do
    end if
    !
    ! replaced previous pointer swap to save memory on GPUs by using already allocated
    ! buffers
    !
    !$acc parallel loop collapse(3) default(present) async(1) ! not really necessary
    !$OMP parallel do   collapse(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dsdtrko(i,j,k) = dsdtrk(i,j,k)
        end do
      end do
    end do
  end subroutine rk_scal
  !
  subroutine cmpt_bulk_forcing(n,nh,nh_g,l,is_forced,velf,dxc,dxf,dyc,dyf,dzc,dzf,u,v,w,f)
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh,nh_g
    real(rp), intent(in   ), dimension(3) :: l
    logical , intent(in   ), dimension(3) :: is_forced
    real(rp), intent(in   ), dimension(3) :: velf
    real(rp), intent(in   ), dimension(1-nh_g:) :: dxc,dxf,dyc,dyf,dzc,dzf
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: u,v,w
    real(rp), intent(out  ), dimension(3) :: f
    real(rp) :: mean
    !
    ! bulk velocity forcing
    !
    f(:) = 0.
    if(is_forced(1)) then
      call bulk_mean(n,nh,nh_g,l,dxc,dyf,dzf,u,mean)
      f(1) = velf(1) - mean
    end if
    if(is_forced(2)) then
      call bulk_mean(n,nh,nh_g,l,dxf,dyc,dzf,v,mean)
      f(2) = velf(2) - mean
    end if
    if(is_forced(3)) then
      call bulk_mean(n,nh,nh_g,l,dxf,dyf,dzc,w,mean)
      f(3) = velf(3) - mean
    end if
  end subroutine cmpt_bulk_forcing
  !
  subroutine bulk_forcing_wall_shear(irk,rkpar,dt,n,nh,nh_g,ng,lo,l,dxci,dxfi,dyci,dyfi,dzci,dzfi, &
                                     dxc,dxf,dyc,dyf,dzc,dzf,dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g, &
                                     visc,is_bound,is_forced,velf,bforce,tau,tauo,u,v,w,f)
    !
    ! bulk forcing that balances exactly the wall shear over a Runge-Kutta substage
    ! (zero net acceleration)
    !
    ! walls with explicit diffusion use the stresses of this and of the previous substage with the
    ! Runge-Kutta weights, and walls with implicit diffusion the stresses before and after the
    ! implicit step, with the forcing applied as f*g, (1 + alpha*lap_W) g = 1 along the wall-bounded
    ! directions W with implicit diffusion, i.e., as if it were added before the implicit step
    !
    ! for expert use, `relax_time > 0` makes the bulk velocity relax towards `velf` over this time scale,
    ! which can counter a slow drift of the bulk velocity due to round-off errors, at the cost of one
    ! volume average per time step
    !
    use mod_bound, only: set_bc,bc_dist
    use mod_param, only: cbcvel,bcvel,cbcpre
    implicit none
    real(rp), parameter :: relax_time = 0.
    integer , intent(in   )                      :: irk
    real(rp), intent(in   ), dimension(2)        :: rkpar
    real(rp), intent(in   )                      :: dt
    integer , intent(in   ), dimension(3)        :: n,ng,lo
    integer , intent(in   )                      :: nh,nh_g
    real(rp), intent(in   ), dimension(3)        :: l
    real(rp), intent(in   ), dimension(1-nh_g:)  :: dxci,dxfi,dyci,dyfi,dzci,dzfi
    real(rp), intent(in   ), dimension(1-nh_g:)  :: dxc,dxf,dyc,dyf,dzc,dzf
    real(rp), intent(in   ), dimension(1-nh_g:)  :: dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g
    real(rp), intent(in   )                      :: visc
    logical , intent(in   ), dimension(0:1,3)    :: is_bound
    logical , intent(in   ), dimension(3)        :: is_forced
    real(rp), intent(in   ), dimension(3)        :: velf,bforce
    real(rp), intent(in   ), dimension(0:1,3,3)  :: tau,tauo
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: u,v,w
    real(rp), intent(out  ), dimension(3)        :: f
    real(rp), dimension(3), save :: ubulk
    logical , save :: is_ubulk = .false.
    real(rp), allocatable, dimension(:,:,:,:), save :: gg
    real(rp), dimension(3), save :: alpha_s = huge(1._rp)
    real(rp), dimension(3), save :: gmean_s
    real(rp), dimension(0:1,3,3) :: tau_i
    real(rp), dimension(3,3) :: tw,two
    real(rp), dimension(3) :: db
    real(rp), allocatable, dimension(:,:,:), save :: dr
    real(rp) :: factor1,factor2,factor12,alpha
    logical , dimension(3) :: is_w
    integer :: idir,ib
    f(:) = 0.
    if(.not.any(is_forced(:))) return ! nothing to force
    !
    ! time-step factors as in `rk`
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    db(:) = 0.
    if(relax_time > 0..and.is_ubulk) db(:) = factor12*(velf(:)-ubulk(:))/relax_time
    do idir=1,3
      tw(:,idir) = sum(tau(:,:,idir),1)/l(:); two(:,idir) = sum(tauo(:,:,idir),1)/l(:)
    end do
    !
    ! wall-bounded directions with implicit diffusion
    !
    is_w(:) = [impdiff_mode == impdiff_xyz, &
               impdiff_mode == impdiff_yz.or.impdiff_mode == impdiff_xyz, &
               impdiff_mode /= impdiff_explicit]
    do idir=1,3
      is_w(idir) = is_w(idir).and.cbcpre(0,idir)//cbcpre(1,idir) /= 'PP'
    end do
    if(.not.any(is_w)) then
      do idir=1,3
        if(.not.is_forced(idir)) cycle
        f(idir) = db(idir) - factor12*bforce(idir) + (factor1*sum(tw(:,idir)) + factor2*sum(two(:,idir)))
      end do
      call bulk_forcing(n,nh,is_forced,f,u,v,w)
    else
      if(.not.allocated(gg)) then
        allocate(gg(merge(n(1),1,is_w(1)),merge(n(2),1,is_w(2)),merge(n(3),1,is_w(3)),3))
        !$acc enter data create(gg)
      end if
      alpha = -.5*visc*(sum(rkpar(:))*dt) ! as in the main program
      if(abs(alpha-alpha_s(irk)) > 0.) then
        alpha_s(irk) = alpha
        call cmpt_unit_response(alpha,n,nh_g,ng,lo,l,is_w,dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g, &
                                gg(:,:,:,irk),gmean_s(irk))
        !$acc update device(gg(:,:,:,irk)) async(1)
      end if
      !
      ! wall shear of the unforced solution, with the ghost values set from the boundary conditions
      ! (`dr` -> mirror distances, computed once)
      !
      if(.not.allocated(dr)) then
        allocate(dr(0:nh,0:1,3))
        do ib=0,1
          dr(:,ib,1) = bc_dist(n(1),nh,nh_g,ib,.true.,dxc)
          dr(:,ib,2) = bc_dist(n(2),nh,nh_g,ib,.true.,dyc)
          dr(:,ib,3) = bc_dist(n(3),nh,nh_g,ib,.true.,dzc)
        end do
      end if
      do idir=1,3
        if(.not.is_w(idir)) cycle
        do ib=0,1
          if(.not.is_bound(ib,idir)) cycle
          if(is_forced(1)) call set_bc(cbcvel(ib,idir,1),ib,idir,nh,.true.,bcvel(ib,idir,1),dr(:,ib,idir),u)
          if(is_forced(2)) call set_bc(cbcvel(ib,idir,2),ib,idir,nh,.true.,bcvel(ib,idir,2),dr(:,ib,idir),v)
          if(is_forced(3)) call set_bc(cbcvel(ib,idir,3),ib,idir,nh,.true.,bcvel(ib,idir,3),dr(:,ib,idir),w)
        end do
      end do
      call cmpt_wallshear(n,nh,nh_g,is_forced,is_bound,l,dxci,dxfi,dyci,dyfi,dzci,dzfi,visc,u,v,w, &
                          tau_i(:,:,1),tau_i(:,:,2),tau_i(:,:,3))
      do idir=1,3
        if(.not.is_forced(idir)) cycle
        f(idir) = (db(idir) - factor12*bforce(idir) + &
                   sum(merge(.5*factor12*(tw(:,idir) + sum(tau_i(:,:,idir),1)/l(:)), &
                             factor1*tw(:,idir) + factor2*two(:,idir),is_w(:))))/gmean_s(irk)
      end do
      if(is_forced(1)) call add_response(n,nh,f(1),gg(:,:,:,irk),u)
      if(is_forced(2)) call add_response(n,nh,f(2),gg(:,:,:,irk),v)
      if(is_forced(3)) call add_response(n,nh,f(3),gg(:,:,:,irk),w)
    end if
    if(relax_time > 0..and.irk == 3) then
      if(is_forced(1)) call bulk_mean(n,nh,nh_g,l,dxc,dyf,dzf,u,ubulk(1))
      if(is_forced(2)) call bulk_mean(n,nh,nh_g,l,dxf,dyc,dzf,v,ubulk(2))
      if(is_forced(3)) call bulk_mean(n,nh,nh_g,l,dxf,dyf,dzc,w,ubulk(3))
      is_ubulk = .true.
    end if
  contains
    subroutine add_response(n,nh,ff,gg,p)
      implicit none
      integer , intent(in   ), dimension(3)        :: n
      integer , intent(in   )                      :: nh
      real(rp), intent(in   )                      :: ff
      real(rp), intent(in   ), dimension(:,:,:)    :: gg
      real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
      integer :: i,j,k,ni,nj,nk
      !
      ! `gg` has extent 1 along the directions where the response is uniform
      !
      ni = size(gg,1); nj = size(gg,2); nk = size(gg,3)
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,k) = p(i,j,k) + ff*gg(min(i,ni),min(j,nj),min(k,nk))
          end do
        end do
      end do
    end subroutine add_response
  end subroutine bulk_forcing_wall_shear
  !
  subroutine cmpt_unit_response(alpha,n,nh,ng,lo,l,is_w,dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g,g,gmean)
    !
    ! response to a uniform unit forcing, (1 + alpha*lap_W) g = 1, with homogeneous Dirichlet conditions
    ! along the directions W where `is_w` is true (one or two, as a forced direction is periodic), solved
    ! on the global grid; returns the local part of g, with extent 1 along the other directions, and its
    ! volume average `gmean`
    !
    ! with two directions, the first is diagonalized (once, as the operator does not depend on alpha)
    ! and a tridiagonal system is solved along the second for each eigenvalue
    !
    use mod_initsolver, only: init_eigendecomp_axis,tridmatrix
    use mod_linalg    , only: gemm
    implicit none
    real(rp), intent(in ) :: alpha
    integer , intent(in ), dimension(3)     :: n,ng,lo
    integer , intent(in )                   :: nh
    real(rp), intent(in ), dimension(3)     :: l
    logical , intent(in ), dimension(3)     :: is_w
    real(rp), intent(in ), dimension(1-nh:) :: dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g
    real(rp), intent(out), dimension(:,:,:) :: g
    real(rp), intent(out)                   :: gmean
    real(rp), allocatable, dimension(:)  , save :: lambda_a,r_a,s_a
    real(rp), allocatable, dimension(:,:), save :: vbwd_a
    real(rp), allocatable, dimension(:)   :: dci_a,dfi_a,dci_b,dfi_b,a,b,c,gl
    real(rp), allocatable, dimension(:,:) :: vfwd,vbwd,gh,gp
    integer , allocatable, dimension(:)   :: iw
    integer :: ia,ib,na,nb,m
    !
    iw = pack([1,2,3],is_w)
    ia = iw(1); na = ng(ia)
    call grid_g(ia,dci_a,dfi_a)
    if(size(iw) == 1) then
      allocate(a(na),b(na),c(na),gl(na))
      call tridmatrix(['D','D'],na,nh,dci_a,dfi_a,'c',.false.,a,b,c)
      gl(:) = 1.
      call solve_tridiag(alpha*a,1.+alpha*b,alpha*c,gl)
      gmean = sum(gl(:)/dfi_a(1:na))/l(ia)
      g(:,:,:) = reshape(gl(lo(ia):lo(ia)-1+n(ia)),shape(g))
    else
      ib = iw(2); nb = ng(ib)
      call grid_g(ib,dci_b,dfi_b)
      if(.not.allocated(lambda_a)) then
        allocate(lambda_a(na),r_a(na),s_a(na),vfwd(na,na),vbwd(na,na))
        call init_eigendecomp_axis(na,nh,dci_a,dfi_a,1./dci_a,1./dfi_a,['D','D'],'c',lambda_a,vfwd,vbwd)
        r_a(:) = sum(vfwd(:,:),2)
        s_a(:) = matmul(1./dfi_a(1:na),vbwd(:,:))
        vbwd_a = vbwd(lo(ia):lo(ia)-1+n(ia),:)
      end if
      allocate(a(nb),b(nb),c(nb),gh(nb,na),gp(n(ia),n(ib)))
      call tridmatrix(['D','D'],nb,nh,dci_b,dfi_b,'c',.false.,a,b,c)
      do m=1,na
        gh(:,m) = r_a(m)
        call solve_tridiag(alpha*a,1.+alpha*(b+lambda_a(m)),alpha*c,gh(:,m))
      end do
      gmean = sum(s_a(:)*matmul(1./dfi_b(1:nb),gh(:,:)))/(l(ia)*l(ib))
      call gemm('N','T',n(ia),n(ib),na,1._rp,vbwd_a,n(ia),gh(lo(ib):lo(ib)-1+n(ib),:),n(ib),0._rp,gp,n(ia))
      g(:,:,:) = reshape(gp,shape(g))
    end if
  contains
    subroutine grid_g(idir,dci,dfi)
      implicit none
      integer , intent(in) :: idir
      real(rp), intent(out), allocatable, dimension(:) :: dci,dfi
      select case(idir)
      case(1)
        dci = dxci_g; dfi = dxfi_g
      case(2)
        dci = dyci_g; dfi = dyfi_g
      case(3)
        dci = dzci_g; dfi = dzfi_g
      end select
    end subroutine grid_g
    subroutine solve_tridiag(a,b,c,p)
      implicit none
      real(rp), intent(in   ), dimension(:) :: a,b,c
      real(rp), intent(inout), dimension(:) :: p
      real(rp), dimension(size(p)) :: d
      real(rp) :: z
      integer :: k,nn
      nn = size(p)
      z = 1._rp/b(1)
      d(1) = c(1)*z
      p(1) = p(1)*z
      do k=2,nn
        z = 1._rp/(b(k) - a(k)*d(k-1))
        d(k) = c(k)*z
        p(k) = (p(k) - a(k)*p(k-1))*z
      end do
      do k=nn-1,1,-1
        p(k) = p(k) - d(k)*p(k+1)
      end do
    end subroutine solve_tridiag
  end subroutine cmpt_unit_response
end module mod_rk
