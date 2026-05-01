! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_solver
  use, intrinsic :: iso_c_binding, only: C_PTR,c_f_pointer,c_loc
  use decomp_2d
  use mod_fft   , only: fft
  use mod_param , only: ipencil_axis,is_poisson_pcr_tdma
  use mod_linalg, only: gemm
  use mod_types
  implicit none
  private
  public solver,solver_gaussel_z,solver_gaussel_yz
  contains
  subroutine solver(n,ng,is_fft,arrplan,normfft,lambdaxy,eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd,a,b,c, &
                    bc,c_or_f,p,is_ptdma_update,aa_z,cc_z)
    !
    ! note: some of the transposes below are suboptimal in slab decompositions,
    ! as they would be a no-op if done in-place (e.g., `px = py` for xy slabs)
    !
    implicit none
    integer , intent(in), dimension(3) :: n,ng
    logical , intent(in), dimension(2) :: is_fft
    type(C_PTR), intent(in), dimension(2,2) :: arrplan
    real(rp), intent(in) :: normfft
    real(rp), intent(in), dimension(:,:) :: lambdaxy,eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd
    real(rp), intent(in), dimension(:) :: a,b,c
    character(len=1), dimension(0:1,3), intent(in) :: bc
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    logical , intent(inout), optional :: is_ptdma_update
    real(rp), intent(inout), dimension(:,:,:), optional :: aa_z,cc_z
    real(rp), target, allocatable, dimension(:,:,:) :: px,py,pz
    real(rp), target, allocatable, dimension(:,:,:) :: px_aux,py_aux_1,py_aux_2
    real(rp), pointer, contiguous, dimension(:,:)   :: px_2d,px_aux_2d,py_aux_1_2d,py_aux_2_2d
    integer :: q,i,j,k
    logical :: is_periodic_z
    integer, dimension(3) :: n_z,hi_z
    logical :: is_ptdma_update_
    real(rp) :: norm
    !
    norm = normfft
    !
    is_ptdma_update_ = .true.
    if(present(is_ptdma_update)) is_ptdma_update_ = is_ptdma_update
    n_z(:)  = zsize(:)
    hi_z(:) = zend(:)
    if(is_poisson_pcr_tdma) then
      n_z(:)  = ysize(:)
      hi_z(:) = yend(:)
    end if
    allocate(px(xsize(1),xsize(2),xsize(3)))
    allocate(py(ysize(1),ysize(2),ysize(3)))
    allocate(pz(zsize(1),zsize(2),zsize(3)))
    if(.not.is_fft(1)) then
      allocate(px_aux(xsize(1),xsize(2),xsize(3)))
      call c_f_pointer(c_loc(px    ),px_2d    ,[xsize(1),xsize(2)*xsize(3)])
      call c_f_pointer(c_loc(px_aux),px_aux_2d,[xsize(1),xsize(2)*xsize(3)])
    end if
    !
    if(.not.is_fft(2)) then
      allocate(py_aux_1(ysize(2),ysize(3),ysize(1)))
      allocate(py_aux_2(ysize(2),ysize(3),ysize(1)))
      call c_f_pointer(c_loc(py_aux_1),py_aux_1_2d,[ysize(2),ysize(3)*ysize(1)])
      call c_f_pointer(c_loc(py_aux_2),py_aux_2_2d,[ysize(2),ysize(3)*ysize(1)])
    end if
    select case(ipencil_axis)
    case(1)
      !$OMP PARALLEL WORKSHARE
      px(:,:,:) = p(1:n(1),1:n(2),1:n(3))
      !$OMP END PARALLEL WORKSHARE
    case(2)
      !$OMP PARALLEL WORKSHARE
      py(:,:,:) = p(1:n(1),1:n(2),1:n(3))
      !$OMP END PARALLEL WORKSHARE
      call transpose_y_to_x(py,px)
    case(3)
      !$OMP PARALLEL WORKSHARE
      pz(:,:,:) = p(1:n(1),1:n(2),1:n(3))
      !$OMP END PARALLEL WORKSHARE
      !call transpose_z_to_x(pz,px)
      call transpose_z_to_y(pz,py)
      call transpose_y_to_x(py,px)
    end select
    !
    if(is_fft(1)) then
      call fft(arrplan(1,1),px) ! fwd transform in x
    else
      call gemm('N','N',ng(1),xsize(2)*xsize(3),ng(1),1._rp, &
                eigvecx_fwd,ng(1),px_2d,ng(1),0._rp,px_aux_2d,ng(1))
      !$OMP PARALLEL WORKSHARE
      px(:,:,:) = px_aux(:,:,:)
      !$OMP END PARALLEL WORKSHARE
    end if
    !
    call transpose_x_to_y(px,py)
    if(is_fft(2)) then
      call fft(arrplan(1,2),py) ! fwd transform in y
    else
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do k=1,ysize(3)
        do j=1,ysize(2)
          do i=1,ysize(1)
            py_aux_1(j,k,i) = py(i,j,k)
          end do
        end do
      end do
      call gemm('N','N',ng(2),ysize(1)*ysize(3),ng(2),1._rp, &
                eigvecy_fwd,ng(2),py_aux_1_2d,ng(2),0._rp,py_aux_2_2d,ng(2))
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do k=1,ysize(3)
        do j=1,ysize(2)
          do i=1,ysize(1)
            py(i,j,k) = py_aux_2(j,k,i)
          end do
        end do
      end do
    end if
    !
    q = merge(1,0,c_or_f(3) == 'f'.and.bc(1,3) == 'D'.and.hi_z(3) == ng(3))
    is_periodic_z = bc(0,3)//bc(1,3) == 'PP'
    if(.not.is_poisson_pcr_tdma) then
      call transpose_y_to_z(py,pz)
      !
      call gaussel(n_z(1),n_z(2),n_z(3)-q,0,a,b,c,is_periodic_z,norm,pz,lambdaxy)
      !
      call transpose_z_to_y(pz,py)
    else
      call gaussel_ptdma(n_z(1),n_z(2),n_z(3)-q,0,a,b,c,is_periodic_z,norm,py,lambdaxy,is_ptdma_update_,aa_z,cc_z)
      if(present(is_ptdma_update)) is_ptdma_update = is_ptdma_update_
    end if
    if(is_fft(2)) then
      call fft(arrplan(2,2),py) ! bwd transform in y
    else
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do k=1,ysize(3)
        do j=1,ysize(2)
          do i=1,ysize(1)
            py_aux_1(j,k,i) = py(i,j,k)
          end do
        end do
      end do
      call gemm('N','N',ng(2),ysize(1)*ysize(3),ng(2),1._rp, &
                eigvecy_bwd,ng(2),py_aux_1_2d,ng(2),0._rp,py_aux_2_2d,ng(2))
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do k=1,ysize(3)
        do j=1,ysize(2)
          do i=1,ysize(1)
            py(i,j,k) = py_aux_2(j,k,i)
          end do
        end do
      end do
    end if
    !
    call transpose_y_to_x(py,px)
    if(is_fft(1)) then
      call fft(arrplan(2,1),px) ! bwd transform in x
    else
      call gemm('N','N',ng(1),xsize(2)*xsize(3),ng(1),1._rp, &
                eigvecx_bwd,ng(1),px_2d,ng(1),0._rp,px_aux_2d,ng(1))
      !$OMP PARALLEL WORKSHARE
      px(:,:,:) = px_aux(:,:,:)
      !$OMP END PARALLEL WORKSHARE
    end if
    !
    select case(ipencil_axis)
    case(1)
      !$OMP PARALLEL WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = px(:,:,:)
      !$OMP END PARALLEL WORKSHARE
    case(2)
      call transpose_x_to_y(px,py)
      !$OMP PARALLEL WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = py(:,:,:)
      !$OMP END PARALLEL WORKSHARE
    case(3)
      !call transpose_x_to_z(px,pz)
      call transpose_x_to_y(px,py)
      call transpose_y_to_z(py,pz)
      !$OMP PARALLEL WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = pz(:,:,:)
      !$OMP END PARALLEL WORKSHARE
    end select
  end subroutine solver
  !
  subroutine gaussel(nx,ny,n,nh,a,b,c,is_periodic,norm,p,lambdaxy)
    use mod_param, only: eps
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(rp), intent(in), dimension(:) :: a,b,c
    logical , intent(in) :: is_periodic
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp), intent(in), dimension(nx,ny), optional :: lambdaxy
    real(rp), allocatable, dimension(:,:,:) :: d,p2
    real(rp) :: z
    integer :: i,j,k,nn
    !
    ! solve tridiagonal system
    !
    nn = n
    if(is_periodic) nn = n-1
    !
    ! allocate work arrays
    !
    allocate(d(nx,ny,nn))
    !
    ! forward elimination
    !
    if(present(lambdaxy)) then
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(z)
      do j=1,ny
        do i=1,nx
          z = 1./(b(1) + lambdaxy(i,j) + eps)
          d(i,j,1) = c(1)*z
          p(i,j,1) = p(i,j,1)*norm*z
        end do
      end do
      !
      do k=2,nn
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(z)
        do j=1,ny
          do i=1,nx
            z = 1./(b(k) + lambdaxy(i,j) - a(k)*d(i,j,k-1) + eps)
            d(i,j,k) = c(k)*z
            p(i,j,k) = (p(i,j,k)*norm - a(k)*p(i,j,k-1))*z
          end do
        end do
      end do
    else
      z = 1./(b(1) + eps)
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do j=1,ny
        do i=1,nx
          d(i,j,1) = c(1)*z
          p(i,j,1) = p(i,j,1)*norm*z
        end do
      end do
      !
      do k=2,nn
        z = 1./(b(k) - a(k)*d(1,1,k-1) + eps)
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
        do j=1,ny
          do i=1,nx
            d(i,j,k) = c(k)*z
            p(i,j,k) = (p(i,j,k)*norm - a(k)*p(i,j,k-1))*z
          end do
        end do
      end do
    end if
    !
    ! backward substitution
    !
    do k=nn-1,1,-1
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do j=1,ny
        do i=1,nx
          p(i,j,k) = p(i,j,k) - d(i,j,k)*p(i,j,k+1)
        end do
      end do
    end do
    !
    ! handle periodic closure with an auxiliary tridiagonal solve
    !
    if(is_periodic) then
      allocate(p2(nx,ny,nn))
      !
      ! initialize the auxiliary right-hand side for the periodic correction
      !
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do j=1,ny
        do i=1,nx
          p2(i,j,1:nn) = 0.
          p2(i,j,1)  = -a(1)
          p2(i,j,nn) = p2(i,j,nn) - c(nn)
        end do
      end do
      !
      ! forward elimination for the auxiliary system
      !
      if(present(lambdaxy)) then
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(z)
        do j=1,ny
          do i=1,nx
            z = 1./(b(1) + lambdaxy(i,j) + eps)
            d(i,j,1) = c(1)*z
            p2(i,j,1) = p2(i,j,1)*z
          end do
        end do
        do k=2,nn
          !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(z)
          do j=1,ny
            do i=1,nx
              z = 1./(b(k) + lambdaxy(i,j) - a(k)*d(i,j,k-1) + eps)
              d(i,j,k) = c(k)*z
              p2(i,j,k) = (p2(i,j,k) - a(k)*p2(i,j,k-1))*z
            end do
          end do
        end do
      else
        z = 1./(b(1) + eps)
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
        do j=1,ny
          do i=1,nx
            d(i,j,1) = c(1)*z
            p2(i,j,1) = p2(i,j,1)*z
          end do
        end do
        do k=2,nn
          z = 1./(b(k) - a(k)*d(1,1,k-1) + eps)
          !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
          do j=1,ny
            do i=1,nx
              d(i,j,k) = c(k)*z
              p2(i,j,k) = (p2(i,j,k) - a(k)*p2(i,j,k-1))*z
            end do
          end do
        end do
      end if
      !
      ! backward substitution for the auxiliary system
      !
      do k=nn-1,1,-1
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
        do j=1,ny
          do i=1,nx
            p2(i,j,k) = p2(i,j,k) - d(i,j,k)*p2(i,j,k+1)
          end do
        end do
      end do
      !
      ! solve for the periodic closure value and correct the interior solution
      !
      if(present(lambdaxy)) then
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
        do j=1,ny
          do i=1,nx
            p(i,j,nn+1) = (p(i,j,nn+1)*norm        - c(nn+1)*p( i,j,1) - a(nn+1)*p(i,j,nn)) / &
                          (b(nn+1) + lambdaxy(i,j) + c(nn+1)*p2(i,j,1) + a(nn+1)*p2(i,j,nn) + eps)
          end do
        end do
      else
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
        do j=1,ny
          do i=1,nx
            p(i,j,nn+1) = (p(i,j,nn+1)*norm - c(nn+1)*p( i,j,1) - a(nn+1)*p( i,j,nn)) / &
                          (b(nn+1)          + c(nn+1)*p2(i,j,1) + a(nn+1)*p2(i,j,nn) + eps)
          end do
        end do
      end if
      !
      ! apply the (Sherman-Morrison) periodic correction to all interior points
      !
      do k=1,nn
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
        do j=1,ny
          do i=1,nx
            p(i,j,k) = p(i,j,k) + p2(i,j,k)*p(i,j,nn+1)
          end do
        end do
      end do
    end if
  end subroutine gaussel
  !
  subroutine gaussel_ptdma(nx,ny,n,nh,a,b,c,is_periodic,norm,p,lambdaxy,is_update,aa_z_save,cc_z_save)
    !
    ! distributed TDMA solver
    !
    use mod_common_mpi, only: dinfo_ptdma
    use mod_param     , only: eps
    !
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(rp), intent(in), dimension(:) :: a,b,c
    logical , intent(in) :: is_periodic
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp), intent(in), dimension(:,:), optional :: lambdaxy
    logical , intent(inout), optional :: is_update
    real(rp), intent(inout), dimension(:,:,:), optional :: aa_z_save,cc_z_save
    real(rp),              dimension(nx,ny,n) :: aa,cc
    real(rp), allocatable, dimension(: ,: ,:) :: aa_y,cc_y,pp_y,aa_z,cc_z,pp_z
    real(rp), allocatable, dimension(: ,: ,:) :: pp_z_2,cc_z_0
    real(rp) :: z,zz(2),bb(n)
    integer :: i,j,k
    integer , dimension(3) :: nr_z
    integer :: nx_r,ny_r,nn
    !
    nr_z(:) = dinfo_ptdma%zsz(:)
    allocate(aa_y(nx,ny,2), &
             cc_y(nx,ny,2), &
             pp_y(nx,ny,2), &
             aa_z(nr_z(1),nr_z(2),nr_z(3)), &
             cc_z(nr_z(1),nr_z(2),nr_z(3)), &
             pp_z(nr_z(1),nr_z(2),nr_z(3)))
    if(is_periodic) then
      allocate(cc_z_0(nr_z(1),nr_z(2),nr_z(3)), &
               pp_z_2(nr_z(1),nr_z(2),nr_z(3)))
    end if
    !
    if(present(lambdaxy)) then
      !
      ! factor inner rows of z-distributed systems so that they are only coupled to the boundaries:
      !
      !$OMP PARALLEL DEFAULT(shared) PRIVATE(zz,bb,z)
      !$OMP DO COLLAPSE(2)
      do j=1,ny
        do i=1,nx
          !
          bb(:) = b(1:n) + lambdaxy(i,j)
          zz(:) = 1./(bb(1:2)+eps)
          aa(i,j,1:2) = a(1:2)*zz(:)
          cc(i,j,1:2) = c(1:2)*zz(:)
          p( i,j,1:2) = p(i,j,1:2)*norm*zz(:)
          !
          ! elimination of lower diagonals
          !
          do k=3,n
            z = 1./(bb(k)-a(k)*cc(i,j,k-1)+eps)
            p(i,j,k) = (p(i,j,k)*norm-a(k)*p(i,j,k-1))*z
            aa(i,j,k) = -a(k)*aa(i,j,k-1)*z
            cc(i,j,k) = c(k)*z
          end do
          !
          ! elimination of upper diagonals
          !
          do k=n-2,2,-1
            p(i,j,k)  = p(i,j,k) - cc(i,j,k)*p(i,j,k+1)
            aa(i,j,k) =  aa(i,j,k)-cc(i,j,k)*aa(i,j,k+1)
            cc(i,j,k) = -cc(i,j,k)*cc(i,j,k+1)
          end do
          z = 1./(1.-aa(i,j,2)*cc(i,j,1)+eps)
          p(i,j,1) = (p(i,j,1)-cc(i,j,1)*p(i,j,2))*z
          aa(i,j,1) = aa(i,j,1)*z
          cc(i,j,1) = -cc(i,j,1)*cc(i,j,2)*z
          !
          ! gather reduced systems
          !
          aa_y(i,j,1) = aa(i,j,1); aa_y(i,j,2) = aa(i,j,n)
          cc_y(i,j,1) = cc(i,j,1); cc_y(i,j,2) = cc(i,j,n)
          pp_y(i,j,1) = p(i,j,1) ; pp_y(i,j,2) = p(i,j,n)
        end do
      end do
      !$OMP END PARALLEL
    else
      !$OMP PARALLEL DEFAULT(shared) PRIVATE(zz,z)
      !$OMP DO COLLAPSE(2)
      do j=1,ny
        do i=1,nx
          zz(:) = 1./(b(1:2)+eps)
          aa(i,j,1:2) = a(1:2)*zz(:)
          cc(i,j,1:2) = c(1:2)*zz(:)
          p( i,j,1:2) = p(i,j,1:2)*norm*zz(:)
          !
          ! elimination of lower diagonals
          !
          do k=3,n
            z = 1./(b(k)-a(k)*cc(i,j,k-1)+eps)
            p(i,j,k) = (p(i,j,k)*norm-a(k)*p(i,j,k-1))*z
            aa(i,j,k) = -a(k)*aa(i,j,k-1)*z
            cc(i,j,k) = c(k)*z
          end do
          !
          ! elimination of upper diagonals
          !
          do k=n-2,2,-1
            p(i,j,k)  = p(i,j,k) - cc(i,j,k)*p(i,j,k+1)
            aa(i,j,k) =  aa(i,j,k)-cc(i,j,k)*aa(i,j,k+1)
            cc(i,j,k) = -cc(i,j,k)*cc(i,j,k+1)
          end do
          z = 1./(1.-aa(i,j,2)*cc(i,j,1)+eps)
          p(i,j,1) = (p(i,j,1)-cc(i,j,1)*p(i,j,2))*z
          aa(i,j,1) = aa(i,j,1)*z
          cc(i,j,1) = -cc(i,j,1)*cc(i,j,2)*z
          !
          ! gather reduced systems
          !
          aa_y(i,j,1) = aa(i,j,1); aa_y(i,j,2) = aa(i,j,n)
          cc_y(i,j,1) = cc(i,j,1); cc_y(i,j,2) = cc(i,j,n)
          pp_y(i,j,1) = p(i,j,1) ; pp_y(i,j,2) = p(i,j,n)
        end do
      end do
      !$OMP END PARALLEL
    end if
    !
    ! transpose to gather reduced subdomain boundary systems along z
    !
    if(present(is_update) .and. present(aa_z_save) .and. present(cc_z_save)) then
      if(is_update) then
        is_update = .false.
        call transpose_y_to_z(aa_y,aa_z_save,dinfo_ptdma)
        call transpose_y_to_z(cc_y,cc_z_save,dinfo_ptdma)
      end if
      aa_z(:,:,:) = aa_z_save(:,:,:)
      cc_z(:,:,:) = cc_z_save(:,:,:)
    else
      call transpose_y_to_z(aa_y,aa_z,dinfo_ptdma)
      call transpose_y_to_z(cc_y,cc_z,dinfo_ptdma)
    end if
    call transpose_y_to_z(pp_y,pp_z,dinfo_ptdma)
    !
    ! solve reduced systems
    !
    nn   = nr_z(3)
    ny_r = nr_z(2)
    nx_r = nr_z(1)
    if(is_periodic) then
      nn = nn-1
      cc_z_0(:,:,:) = cc_z(:,:,:)
    end if
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(z)
    !$OMP DO COLLAPSE(2)
    do j=1,ny_r
      do i=1,nx_r
        do k=2,nn
          z = 1./(1.-aa_z(i,j,k)*cc_z(i,j,k-1)+eps)
          pp_z(i,j,k) = (pp_z(i,j,k)-aa_z(i,j,k)*pp_z(i,j,k-1))*z
          cc_z(i,j,k) = cc_z(i,j,k)*z
        end do
        do k=nn-1,1,-1
          pp_z(i,j,k) = pp_z(i,j,k) - cc_z(i,j,k)*pp_z(i,j,k+1)
        end do
      end do
    end do
    !$OMP END PARALLEL
    if(is_periodic) then
      associate(cc_z => cc_z_0)
      !$OMP PARALLEL DEFAULT(shared) PRIVATE(z)
      !$OMP DO COLLAPSE(2)
      do j=1,ny_r
        do i=1,nx_r
          pp_z_2(i,j,1:nn) = 0.
          pp_z_2(i,j,1 ) = -aa_z(i,j,1 )
          pp_z_2(i,j,nn) = pp_z_2(i,j,nn) - cc_z(i,j,nn)
          !
          do k=2,nn
            z = 1./(1.-aa_z(i,j,k)*cc_z(i,j,k-1)+eps)
            pp_z_2(i,j,k) = (pp_z_2(i,j,k)-aa_z(i,j,k)*pp_z_2(i,j,k-1))*z
            cc_z(i,j,k) = cc_z(i,j,k)*z
          end do
          !
          do k=nn-1,1,-1
            pp_z_2(i,j,k) = pp_z_2(i,j,k) - cc_z(i,j,k)*pp_z_2(i,j,k+1)
          end do
          pp_z(i,j,nn+1) = (pp_z(i,j,nn+1) - cc_z(i,j,nn+1)*pp_z(  i,j,1) - aa_z(i,j,nn+1)*pp_z(  i,j,nn)) / &
                           (1.             + cc_z(i,j,nn+1)*pp_z_2(i,j,1) + aa_z(i,j,nn+1)*pp_z_2(i,j,nn)+eps)
          do k=1,nn
            pp_z(i,j,k) = pp_z(i,j,k) + pp_z_2(i,j,k)*pp_z(i,j,nn+1)
          end do
        end do
      end do
      !$OMP END PARALLEL
      end associate
    end if
    !
    ! transpose solution to the original z-distributed form
    !
    call transpose_z_to_y(pp_z,pp_y,dinfo_ptdma)
    !
    ! obtain final solution on the inner points
    !
    !$OMP PARALLEL DEFAULT(shared)
    !$OMP DO COLLAPSE(2)
    do j=1,ny
      do i=1,nx
        p(i,j,1) = pp_y(i,j,1)
        p(i,j,n) = pp_y(i,j,2)
        do k=2,n-1
          p(i,j,k) = p(i,j,k) - aa(i,j,k)*p(i,j,1) - cc(i,j,k)*p(i,j,n)
        end do
      end do
    end do
    !$OMP END PARALLEL
  end subroutine gaussel_ptdma
  !
  subroutine dgtsv_homebrewed(n,a,b,c,norm,p)
    use mod_param, only: eps
    implicit none
    integer , intent(in) :: n
    real(rp), intent(in   ), dimension(:) :: a,b,c
    real(rp), intent(in   )               :: norm
    real(rp), intent(inout), dimension(:) :: p
    real(rp), dimension(n) :: d
    real(rp) :: z
    integer :: l
    !
    ! Gauss elimination
    !
    z = 1./(b(1)+eps)
    d(1) = c(1)*z
    p(1) = p(1)*norm*z
    do l=2,n
      z    = 1./(b(l)-a(l)*d(l-1)+eps)
      d(l) = c(l)*z
      p(l) = (p(l)*norm-a(l)*p(l-1))*z
    end do
    !
    ! backward substitution
    !
    do l=n-1,1,-1
      p(l) = p(l) - d(l)*p(l+1)
    end do
  end subroutine dgtsv_homebrewed
  !
  subroutine solver_gaussel_z(n,ng,hi,a,b,c,bcz,c_or_f,norm,p)
    implicit none
    integer , intent(in), dimension(3) :: n,ng,hi
    real(rp), intent(in), dimension(:) :: a,b,c
    character(len=1), dimension(0:1), intent(in) :: bcz
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), allocatable, dimension(:,:,:) :: px,py,pz
    integer :: q
    logical :: is_periodic_z
    integer, dimension(3) :: n_z,hi_z
    logical :: is_no_decomp_z
    !
    n_z(:)  = zsize(:)
    hi_z(:) = zend(:)
    if(is_poisson_pcr_tdma) then
      n_z(:)  = ysize(:)
      hi_z(:) = yend(:)
    end if
    is_no_decomp_z = xsize(3) == n_z(3).or.ipencil_axis == 3 ! not decomposed along z: xsize(3) == ysize(3) == ng(3) when dims(2) = 1
    if(.not.is_poisson_pcr_tdma .and. .not.is_no_decomp_z) then
      allocate(px(xsize(1),xsize(2),xsize(3)))
      allocate(py(ysize(1),ysize(2),ysize(3)))
      allocate(pz(zsize(1),zsize(2),zsize(3)))
      select case(ipencil_axis)
      case(1)
        !$OMP PARALLEL WORKSHARE
        px(:,:,:) = p(1:n(1),1:n(2),1:n(3))
        !$OMP END PARALLEL WORKSHARE
        !call transpose_x_to_z(px,pz)
        call transpose_x_to_y(px,py)
        call transpose_y_to_z(py,pz)
      case(2)
        !$OMP PARALLEL WORKSHARE
        py(:,:,:) = p(1:n(1),1:n(2),1:n(3))
        !$OMP END PARALLEL WORKSHARE
        call transpose_y_to_z(py,pz)
      end select
    end if
    !
    q = merge(1,0,c_or_f(3) == 'f'.and.bcz(1) == 'D'.and.hi_z(3) == ng(3))
    is_periodic_z = bcz(0)//bcz(1) == 'PP'
    if(.not.is_no_decomp_z) then
      if(.not.is_poisson_pcr_tdma) then
        call gaussel(      n_z(1),n_z(2),n_z(3)-q,0,a,b,c,is_periodic_z,norm,pz)
      else
        call gaussel_ptdma(n_z(1),n_z(2),n_z(3)-q,1,a,b,c,is_periodic_z,norm,p)
      end if
    else
      call gaussel(n(1),n(2),n(3)-q,1,a,b,c,is_periodic_z,norm,p)
    end if
    !
    if(.not.is_poisson_pcr_tdma .and. .not.is_no_decomp_z) then
      select case(ipencil_axis)
      case(1)
        !call transpose_z_to_x(pz,px)
        call transpose_z_to_y(pz,py)
        call transpose_y_to_x(py,px)
        !$OMP PARALLEL WORKSHARE
        p(1:n(1),1:n(2),1:n(3)) = px(:,:,:)
        !$OMP END PARALLEL WORKSHARE
      case(2)
        call transpose_z_to_y(pz,py)
        !$OMP PARALLEL WORKSHARE
        p(1:n(1),1:n(2),1:n(3)) = py(:,:,:)
        !$OMP END PARALLEL WORKSHARE
      end select
    end if
  end subroutine solver_gaussel_z
  !
  subroutine solver_gaussel_yz(n,ng,is_fft,arrplan,normfft,lambday,eigvecy_fwd,eigvecy_bwd,a,b,c,bc,c_or_f,p)
    implicit none
    integer , intent(in), dimension(3) :: n,ng
    logical , intent(in), dimension(2) :: is_fft
    type(C_PTR), intent(in), dimension(2,2) :: arrplan
    real(rp), intent(in) :: normfft
    real(rp), intent(in), dimension(:) :: lambday
    real(rp), intent(in), dimension(:,:) :: eigvecy_fwd,eigvecy_bwd
    real(rp), intent(in), dimension(:) :: a,b,c
    character(len=1), dimension(0:1,3), intent(in) :: bc
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), target, allocatable, dimension(:,:,:) :: py,pz
    real(rp), target, allocatable, dimension(:,:,:) :: py_aux_1,py_aux_2
    real(rp), pointer, contiguous, dimension(:,:) :: py_aux_1_2d,py_aux_2_2d
    integer :: i,j,k,q
    logical :: is_periodic_z
    !
    allocate(py(ysize(1),ysize(2),ysize(3)))
    if(ipencil_axis == 3) allocate(pz(zsize(1),zsize(2),zsize(3)))
    if(.not.is_fft(2)) then
      allocate(py_aux_1(ysize(2),ysize(3),ysize(1)))
      allocate(py_aux_2(ysize(2),ysize(3),ysize(1)))
      call c_f_pointer(c_loc(py_aux_1),py_aux_1_2d,[ysize(2),ysize(3)*ysize(1)])
      call c_f_pointer(c_loc(py_aux_2),py_aux_2_2d,[ysize(2),ysize(3)*ysize(1)])
    end if
    select case(ipencil_axis)
    case(2)
      !$OMP PARALLEL WORKSHARE
      py(:,:,:) = p(1:n(1),1:n(2),1:n(3))
      !$OMP END PARALLEL WORKSHARE
    case(3)
      !$OMP PARALLEL WORKSHARE
      pz(:,:,:) = p(1:n(1),1:n(2),1:n(3))
      !$OMP END PARALLEL WORKSHARE
      call transpose_z_to_y(pz,py)
    end select
    !
    if(is_fft(2)) then
      call fft(arrplan(1,2),py)
    else
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do k=1,ysize(3)
        do j=1,ysize(2)
          do i=1,ysize(1)
            py_aux_1(j,k,i) = py(i,j,k)
          end do
        end do
      end do
      call gemm('N','N',ng(2),ysize(1)*ysize(3),ng(2),1._rp, &
                eigvecy_fwd,ng(2),py_aux_1_2d,ng(2),0._rp,py_aux_2_2d,ng(2))
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do k=1,ysize(3)
        do j=1,ysize(2)
          do i=1,ysize(1)
            py(i,j,k) = py_aux_2(j,k,i)
          end do
        end do
      end do
    end if
    !
    q = merge(1,0,c_or_f(3) == 'f'.and.bc(1,3) == 'D'.and.yend(3) == ng(3))
    is_periodic_z = bc(0,3)//bc(1,3) == 'PP'
    call gaussel_yz(ysize(1),ysize(2),ysize(3)-q,0,a,b,c,is_periodic_z,normfft,py,lambday)
    !
    if(is_fft(2)) then
      call fft(arrplan(2,2),py)
    else
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do k=1,ysize(3)
        do j=1,ysize(2)
          do i=1,ysize(1)
            py_aux_1(j,k,i) = py(i,j,k)
          end do
        end do
      end do
      call gemm('N','N',ng(2),ysize(1)*ysize(3),ng(2),1._rp, &
                eigvecy_bwd,ng(2),py_aux_1_2d,ng(2),0._rp,py_aux_2_2d,ng(2))
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do k=1,ysize(3)
        do j=1,ysize(2)
          do i=1,ysize(1)
            py(i,j,k) = py_aux_2(j,k,i)
          end do
        end do
      end do
    end if
    !
    select case(ipencil_axis)
    case(2)
      !$OMP PARALLEL WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = py(:,:,:)
      !$OMP END PARALLEL WORKSHARE
    case(3)
      call transpose_y_to_z(py,pz)
      !$OMP PARALLEL WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = pz(:,:,:)
      !$OMP END PARALLEL WORKSHARE
    end select
  end subroutine solver_gaussel_yz
  !
  subroutine gaussel_yz(nx,ny,n,nh,a,b,c,is_periodic,norm,p,lambday)
    use mod_param, only: eps
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(rp), intent(in), dimension(:) :: a,b,c,lambday
    logical , intent(in) :: is_periodic
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp), allocatable, dimension(:,:,:) :: d,p2
    real(rp) :: z,lambda
    integer :: i,j,k,nn
    !
    nn = n
    if(is_periodic) nn = n-1
    allocate(d(nx,ny,nn))
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(z,lambda)
    do j=1,ny
      do i=1,nx
        lambda = lambday(j)
        z = 1./(b(1) + lambda + eps)
        d(i,j,1) = c(1)*z
        p(i,j,1) = p(i,j,1)*norm*z
      end do
    end do
    do k=2,nn
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(z,lambda)
      do j=1,ny
        do i=1,nx
          lambda = lambday(j)
          z = 1./(b(k) + lambda - a(k)*d(i,j,k-1) + eps)
          d(i,j,k) = c(k)*z
          p(i,j,k) = (p(i,j,k)*norm - a(k)*p(i,j,k-1))*z
        end do
      end do
    end do
    do k=nn-1,1,-1
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do j=1,ny
        do i=1,nx
          p(i,j,k) = p(i,j,k) - d(i,j,k)*p(i,j,k+1)
        end do
      end do
    end do
    if(is_periodic) then
      allocate(p2(nx,ny,nn))
      !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(shared)
      do k=1,nn
        do j=1,ny
          do i=1,nx
            p2(i,j,k) = 0.
          end do
        end do
      end do
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
      do j=1,ny
        do i=1,nx
          p2(i,j,1 ) = -a(1 )
          p2(i,j,nn) = -c(nn)
        end do
      end do
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(z,lambda)
      do j=1,ny
        do i=1,nx
          lambda = lambday(j)
          z = 1./(b(1) + lambda + eps)
          d( i,j,1) = c(1)*z
          p2(i,j,1) = p2(i,j,1)*z
        end do
      end do
      do k=2,nn
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(z,lambda)
        do j=1,ny
          do i=1,nx
            lambda = lambday(j)
            z = 1./(b(k) + lambda - a(k)*d(i,j,k-1) + eps)
            d( i,j,k) = c(k)*z
            p2(i,j,k) = (p2(i,j,k) - a(k)*p2(i,j,k-1))*z
          end do
        end do
      end do
      do k=nn-1,1,-1
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
        do j=1,ny
          do i=1,nx
            p2(i,j,k) = p2(i,j,k) - d(i,j,k)*p2(i,j,k+1)
          end do
        end do
      end do
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(lambda)
      do j=1,ny
        do i=1,nx
          lambda = lambday(j)
          p(i,j,nn+1) = (p(i,j,nn+1)*norm - c(nn+1)*p(i,j,1) - a(nn+1)*p(i,j,nn)) / &
                        (b(nn+1) + lambda + c(nn+1)*p2(i,j,1) + a(nn+1)*p2(i,j,nn) + eps)
        end do
      end do
      do k=1,nn
        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared)
        do j=1,ny
          do i=1,nx
            p(i,j,k) = p(i,j,k) + p2(i,j,k)*p(i,j,nn+1)
          end do
        end do
      end do
    end if
  end subroutine gaussel_yz
  !
#if 0
  subroutine gaussel_lapack(nx,ny,n,a,b,c,p)
    implicit none
#if !defined(_SINGLE_PRECISION)
    external :: dgttrf,dgttrs
    procedure(), pointer :: gttrf => dgttrf, gttrs => dgttrs
#else
    external :: sgttrf,sgttrs
    procedure(), pointer :: gttrf => sgttrf, gttrs => sgttrs
#endif
    integer , intent(in) :: nx,ny,n
    real(rp), intent(in), dimension(:) :: a,b,c
    real(rp), intent(inout), dimension(:,:,:) :: p
    real(rp), allocatable, dimension(:) :: aa,bb,cc,ccc
    integer , allocatable, dimension(:) :: ipiv
    integer :: i,j,info
    !real(rp), dimension(n,nx,ny) :: p_t
    !
    allocate(aa,source=a(2:n  ))
    allocate(bb,source=b(1:n  ))
    allocate(cc,source=c(1:n-1))
    allocate(ccc(n-2),ipiv(n))
    call gttrf(n,aa,bb,cc,ccc,ipiv,info)
    do j=1,ny
      do i=1,nx
        call gttrs('N',n,1,aa,bb,cc,ccc,ipiv,p(i,j,1:n),n,info)
      end do
    end do
    !p_t = reshape(p(1:nx,1:ny,1:n),shape(p_t),order=[2,3,1])
    !call gttrs('N',n,nx*ny,aa,bb,cc,ccc,ipiv,p_t(1:n,:,:),n,info)
    !p(1:nx,1:ny,1:n) = reshape(p_t,shape(p(1:nx,1:ny,1:n)),order=[3,1,2])
  end subroutine gaussel_lapack
#endif
end module mod_solver
