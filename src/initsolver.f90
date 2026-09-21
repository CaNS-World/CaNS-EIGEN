! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_initsolver
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mod_fft   , only: fftini
  use mod_linalg, only: stedc,syevd
  use mod_types
  implicit none
  private
  public initsolver
  contains
  subroutine initsolver(is_poisson_fft,ng,n_x_fft,n_y_fft,lo_z,hi_z,dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g,cbc,bc, &
                        lambdax_g,lambday_g,lambdaxy,eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd,c_or_f,a,b,c,arrplan,normfft, &
                        rhsbx,rhsby,rhsbz)
    !
    ! initializes the Poisson/Helmholtz solver
    !
    implicit none
    logical , intent(in), dimension(2) :: is_poisson_fft
    integer , intent(in), dimension(3) :: ng,n_x_fft,n_y_fft,lo_z,hi_z
    real(rp), intent(in), dimension(0:) :: dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(rp)        , intent(in), dimension(0:1,3) :: bc
    real(rp), intent(out), dimension(:) :: lambdax_g,lambday_g
    real(rp), intent(out), dimension(lo_z(1):,lo_z(2):) :: lambdaxy
    real(rp), intent(out), dimension(:,:) :: eigvecx_fwd,eigvecx_bwd, &
                                             eigvecy_fwd,eigvecy_bwd
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(out), dimension(lo_z(3):) :: a,b,c
#if !(defined(_OPENACC) || defined(_OPENMP)) || defined(_USE_HIP)
    type(C_PTR), intent(out), dimension(2,2) :: arrplan
#else
    integer    , intent(out), dimension(2,2) :: arrplan
#endif
    real(rp), intent(out), dimension(:,:,0:) :: rhsbx
    real(rp), intent(out), dimension(:,:,0:) :: rhsby
    real(rp), intent(out), dimension(:,:,0:) :: rhsbz
    real(rp), intent(out), dimension(2) :: normfft
    real(rp), dimension(2)         :: dl,dli
    real(rp), dimension(0:ng(1)+1) :: dxc_g,dxf_g
    real(rp), dimension(0:ng(2)+1) :: dyc_g,dyf_g
    real(rp), dimension(0:ng(3)+1) :: dzc_g,dzf_g
    integer :: i,j
    real(rp), dimension(ng(3))      :: az_g,bz_g,cz_g
    !
    dli(1) = dxfi_g(0)
    dli(2) = dyfi_g(0)
    dl(:) = dli(:)**(-1)
    dxc_g(:) = dxci_g(:)**(-1)
    dxf_g(:) = dxfi_g(:)**(-1)
    dyc_g(:) = dyci_g(:)**(-1)
    dyf_g(:) = dyfi_g(:)**(-1)
    dzc_g(:) = dzci_g(:)**(-1)
    dzf_g(:) = dzfi_g(:)**(-1)
    !
    ! Generating eigenvalues/eigenvectors consistent with the BCs:
    !
    !  - For     uniform grid spacing, use the standard FFT-based approach.
    !  - For non-uniform grid spacing, use the generalization of the former approach by performing a numerical eigendecomposition:
    !
    !    1. Symmetrize tri-diagonal matrix using a similarity transformation:
    !
    !         T^{symm} = M^{1/2} T M^{-1/2},
    !
    !       where M = diag(delta(1),...,delta(n)) is a 'weight/mass' matrix with the cell widths (see tridmatrix for more details).
    !    2. Compute eigendecomposition using LAPACK's `stedc`:
    !
    !         T^{symm} = U Lambda U^t,
    !
    !       with Lambda the diagonal matrix of eigenvalues and U the matrix of eigenvectors.
    !    3. Re-scale eigenvectors with similarity transformation above, to obtain the generalized eigenvectors:
    !
    !         Q = M^{-1/2} U; Q^{-1} = U^t M^{1/2} => T = Q Lambda Q^{-1},
    !
    !       where, in the Poisson solver:
    !         - Q^{-1} provides the forward  transforms to be applied along x/y (`eigvecx/y_fwd` below), and
    !         - Q      provides the backward transforms to be applied along x/y (`eigvecx/y_bwd` below).
    !
    if(is_poisson_fft(1)) then
      call eigenvalues(ng(1),cbc(:,1),c_or_f(1),lambdax_g)
      lambdax_g(:) = lambdax_g(:)*dli(1)**2
    else
      call init_eigendecomp_axis(ng(1),dxci_g,dxfi_g,dxc_g,dxf_g,cbc(:,1),c_or_f(1),lambdax_g,eigvecx_fwd,eigvecx_bwd)
    end if
    if(is_poisson_fft(2)) then
      call eigenvalues(ng(2),cbc(:,2),c_or_f(2),lambday_g)
      lambday_g(:) = lambday_g(:)*dli(2)**2
    else
      call init_eigendecomp_axis(ng(2),dyci_g,dyfi_g,dyc_g,dyf_g,cbc(:,2),c_or_f(2),lambday_g,eigvecy_fwd,eigvecy_bwd)
    end if
    !
    ! add eigenvalues
    !
    do j=lo_z(2),hi_z(2)
      do i=lo_z(1),hi_z(1)
        lambdaxy(i,j) = lambdax_g(i)+lambday_g(j)
      end do
    end do
    !
    ! compute and distribute coefficients for tridiagonal solver
    !
    call tridmatrix(cbc(:,3),ng(3),dzci_g,dzfi_g,c_or_f(3),.false.,az_g,bz_g,cz_g)
    a(:) = az_g(lo_z(3):hi_z(3))
    b(:) = bz_g(lo_z(3):hi_z(3))
    c(:) = cz_g(lo_z(3):hi_z(3))
    !
    ! compute values to be added to the right hand side
    !
    if(     c_or_f(1) == 'c') then
      call bc_rhs(cbc(:,1),bc(:,1),[dxc_g(0),dxc_g(ng(1)  )],[dxf_g(1),dxf_g(ng(1))],c_or_f(1),rhsbx)
    else if(c_or_f(1) == 'f') then
      call bc_rhs(cbc(:,1),bc(:,1),[dxc_g(1),dxc_g(ng(1)-1)],[dxf_g(1),dxf_g(ng(1))],c_or_f(1),rhsbx)
    end if
    if(     c_or_f(2) == 'c') then
      call bc_rhs(cbc(:,2),bc(:,2),[dyc_g(0),dyc_g(ng(2)  )],[dyf_g(1),dyf_g(ng(2))],c_or_f(2),rhsby)
    else if(c_or_f(2) == 'f') then
      call bc_rhs(cbc(:,2),bc(:,2),[dyc_g(1),dyc_g(ng(2)-1)],[dyf_g(1),dyf_g(ng(2))],c_or_f(2),rhsby)
    end if
    if(     c_or_f(3) == 'c') then
      call bc_rhs(cbc(:,3),bc(:,3),[dzc_g(0),dzc_g(ng(3)  )],[dzf_g(1),dzf_g(ng(3))],c_or_f(3),rhsbz)
    else if(c_or_f(3) == 'f') then
      call bc_rhs(cbc(:,3),bc(:,3),[dzc_g(1),dzc_g(ng(3)-1)],[dzf_g(1),dzf_g(ng(3))],c_or_f(3),rhsbz)
    end if
    !
    ! prepare ffts
    !
    call fftini(ng,is_poisson_fft,n_x_fft,n_y_fft,cbc(:,1:2),c_or_f(1:2),arrplan,normfft)
  end subroutine initsolver
  !
  subroutine init_eigendecomp_axis(n,dci_g,dfi_g,dc_g,df_g,cbc,c_or_f,lambda_g,eigvec_fwd,eigvec_bwd)
    !
    ! numerical eigendecomposition along one direction
    !
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: ieee_exceptions, only: ieee_status_type,ieee_get_status,ieee_set_status, &
                                             ieee_set_halting_mode,ieee_all
    implicit none
    integer , intent(in ) :: n
    real(rp), intent(in ), dimension(0:n+1) :: dci_g,dfi_g,dc_g,df_g
    character(len=1), intent(in), dimension(0:1) :: cbc
    character(len=1), intent(in) :: c_or_f
    real(rp), intent(out), dimension(n) :: lambda_g
    real(rp), intent(out), dimension(n,n) :: eigvec_fwd,eigvec_bwd
    integer :: q,i,j
    real(rp), dimension(n) :: a_g,b_g,c_g
    real(rp), allocatable, dimension(:,:) :: eigvecs
    real(rp), allocatable, dimension(:) :: work
    integer , allocatable, dimension(:) :: iwork
    integer :: wsize,iwsize,info
    type(ieee_status_type) :: fp_status
    !
    q = merge(1,0,c_or_f == 'f'.and.cbc(1) /= 'P')
    call tridmatrix(cbc,n,dci_g,dfi_g,c_or_f,.true.,a_g,b_g,c_g)
    allocate(eigvecs(n,n),work(1),iwork(1))
    !
    ! preserve the caller's IEEE state while LAPACK uses nonhalting arithmetic
    !
    call ieee_get_status(fp_status)
    call ieee_set_halting_mode(ieee_all,.false.)
    if(cbc(0)//cbc(1) /= 'PP') then
      !
      ! non-periodic BCs: simple symmetric tridiagonal matrix
      !
      call stedc('I',n-q,b_g,c_g,eigvecs(1:n-q,1:n-q),n-q,work,-1,iwork,-1,info) ! workspace size query
      if(info /= 0) error stop 'ERROR: LAPACK workspace query failed.'
      wsize = int(work(1),kind(wsize)); iwsize = iwork(1)
      deallocate(work,iwork)
      allocate(work(wsize),iwork(iwsize))
      call stedc('I',n-q,b_g,c_g,eigvecs(1:n-q,1:n-q),n-q,work,wsize,iwork,iwsize,info)
      lambda_g(:) = b_g(:)
    else
      !
      ! periodic BCs: define full cyclic symmetric tridiagonal matrix (upper diagonal)
      !
      eigvecs(:,:) = 0.
      do i=1,n-q
        eigvecs(i,  i) = b_g(i)
      end do
      do i=1,n-q-1
        eigvecs(i,i+1) = c_g(i)
      end do
      !
      ! the cyclic contribution shares an entry with the regular stencil for n <= 2
      !
      eigvecs(1,n-q) = eigvecs(1,n-q) + c_g(n-q)
      call syevd('V','U',n-q,eigvecs(1:n-q,1:n-q),n-q,lambda_g(1:n-q),work,-1,iwork,-1,info) ! workspace size query
      if(info /= 0) error stop 'ERROR: LAPACK workspace query failed.'
      wsize = int(work(1),kind(wsize)); iwsize = iwork(1)
      deallocate(work,iwork)
      allocate(work(wsize),iwork(iwsize))
      call syevd('V','U',n-q,eigvecs(1:n-q,1:n-q),n-q,lambda_g(1:n-q),work,wsize,iwork,iwsize,info)
    end if
    call ieee_set_status(fp_status)
    if(info /= 0) error stop 'ERROR: LAPACK eigendecomposition failed.'
    if(.not.all(ieee_is_finite(lambda_g(1:n-q)))) error stop 'ERROR: nonfinite eigenvalues.'
    !
    ! set the constant eigenvalue to zero, independently of the mode ordering
    !
    if(cbc(0)//cbc(1) == 'PP'.or.cbc(0)//cbc(1) == 'NN') &
      lambda_g(minloc(abs(lambda_g(1:n-q)),dim=1)) = 0.
    !
    ! compute generalized eigenvectors
    !
    select case(c_or_f)
    case('c')
      do j=1,n
        do i=1,n
          eigvec_fwd(i,j) = eigvecs(j,i)*sqrt(df_g(j))
          eigvec_bwd(i,j) = sqrt(df_g(i))**(-1)*eigvecs(i,j)
        end do
      end do
    case('f')
      if(q == 1) then ! set trivial equation for the boundary point
        lambda_g(n) = 0.
        eigvecs(n,:) = 0.
        eigvecs(:,n) = 0.
        eigvecs(n,n) = 1.
      end if
      do j=1,n
        do i=1,n
          eigvec_fwd(i,j) = eigvecs(j,i)*sqrt(dc_g(j))
          eigvec_bwd(i,j) = sqrt(dc_g(i))**(-1)*eigvecs(i,j)
        end do
      end do
    end select
    deallocate(work,iwork,eigvecs)
  end subroutine init_eigendecomp_axis
  !
  subroutine eigenvalues(n,bc,c_or_f,lambda)
    use mod_param, only: pi
    implicit none
    integer , intent(in ) :: n
    character(len=1), intent(in), dimension(0:1) :: bc
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f -> face-centered
    real(rp), intent(out), dimension(n) :: lambda
    integer :: l
    select case(bc(0)//bc(1))
    case('PP')
      do l=1,n
        lambda(l  )   = -2.*(1.-cos((2*(l-1))*pi/(1.*n)))
      end do
#if defined(_OPENACC) || defined(_OPENMP)
      block
        !
        ! new format: (r[0],r[n],r[1],i[1],...,r[n-1],i[n-1])
        ! note that i[0] = i[n] = 0 in a R2C DFT
        !
        integer :: nh,iswap(n)
        nh = (n+1)/2
        iswap(1) = 1
        if(n > 1) iswap(2) = nh+(1-mod(n,2))
        do l=2,n-1
          if(l <= nh) then ! real eigenvalue
            iswap(2*l-1                  ) = l
          else             ! imaginary eigenvalue
            iswap(n-2*(l-(nh+1))-mod(n,2)) = l+1
          end if
        end do
        lambda(:) = lambda(iswap(:))
      end block
#endif
    case('NN')
      if(     c_or_f == 'c') then
        do l=1,n
          lambda(l)   = -2.*(1.-cos((l-1  )*pi/(1.*n)))
        end do
      else if(c_or_f == 'f') then
        do l=1,n-1 ! point at n is a dependent boundary value
          lambda(l)   = -2.*(1.-cos((l-1  )*pi/(1.*(n-1))))
        end do
        lambda(n) = 0.
      end if
    case('DD')
      if(     c_or_f == 'c') then
        do l=1,n
          lambda(l)   = -2.*(1.-cos((l    )*pi/(1.*n)))
        end do
      else if(c_or_f == 'f') then
        do l=1,n-1 ! point at n is a boundary and is excluded here
          lambda(l)   = -2.*(1.-cos((l    )*pi/(1.*(n+1-1))))
        end do
        lambda(n) = 0.
      end if
    case('ND','DN')
      if(     c_or_f == 'c') then
        do l=1,n
          lambda(l)   = -2.*(1.-cos((2*l-1)*pi/(2.*n)))
        end do
      else if(c_or_f == 'f') then
        do l=1,n-1 ! point at n is prescribed by the boundary condition
          lambda(l)   = -2.*(1.-cos((2*l-1)*pi/(2.*n-1.)))
        end do
        lambda(n) = 0.
      end if
    end select
  end subroutine eigenvalues
  !
  subroutine tridmatrix(bc,n,dzci,dzfi,c_or_f,is_symm,a,b,c)
    implicit none
    character(len=1), intent(in), dimension(0:1) :: bc
    integer , intent(in) :: n
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f-face-centered
    logical , intent(in) :: is_symm
    real(rp), intent(out), dimension(n) :: a,b,c
    integer :: k
    integer :: ibound
    real(rp), dimension(0:1) :: factor
    !
    if(n == 1.and.bc(0)//bc(1) == 'PP') then
      a(:) = 0.; b(:) = 0.; c(:) = 0.
      return
    end if
    !
    select case(c_or_f)
    case('c')
      do k=1,n
        a(k) = dzfi(k)*dzci(k-1)
        c(k) = dzfi(k)*dzci(k)
      end do
    case('f')
      do k = 1,n
        a(k) = dzci(k)*dzfi(k)
        c(k) = dzci(k)*dzfi(k+1)
      end do
    end select
    b(:) = -(a(:)+c(:))
    do ibound = 0,1
      select case(bc(ibound))
      case('P')
        factor(ibound) = 0.
      case('D')
        factor(ibound) = -1.
      case('N')
        factor(ibound) = 1.
      end select
    end do
    select case(c_or_f)
    !
    ! If the matrix needs to be symmetrized for the numerical eigendecomposition,
    ! one can set `is_symm = .true.`
    !
    ! Note:
    !
    !   Left-multiplying the discrete 1D operator by the local cell widths M = diag(dz(1),...,dz(n)) makes the matrix symmetric:
    !
    !   (M T)^t = (T^t) M = M T (i.e., T is M-symmetric; M is often referred to as a weight/mass matrix).
    !
    !   However, in multiple dimensions, this left-multiply breaks the Kronecker-product (separable) structure
    !   and introduces cross-directional couplings that prevent solving the problem dimension-by-dimension.
    !
    !   Instead, to preserve separability while obtaining a symmetric tridiagonal matrix,
    !   one can use the following similarity transform:
    !
    !   T^{symm} = M^{1/2} T M^{-1/2}.
    !
    case('c')
      b(1) = b(1) + factor(0)*a(1)
      b(n) = b(n) + factor(1)*c(n)
      if(is_symm) then ! similarity transform
        !a(1:n) = a(1:n)*sqrt(dzfi(0:n-1)/dzfi(1:n))
        c(1:n) = c(1:n)*sqrt(dzfi(2:n+1)/dzfi(1:n)) ! include grid BCs to handle cyclic matrices
      end if
    case('f')
      if(bc(0) == 'N') b(1) = b(1) + factor(0)*a(1)
      if(bc(1) == 'N') b(n-1) = b(n-1) + factor(1)*c(n-1)
      if(is_symm) then ! similarity transform
        !a(1:n) = a(1:n)*sqrt(dzci(0:n-1)/dzfi(1:n))
        c(1:n) = c(1:n)*sqrt(dzci(2:n+1)/dzci(1:n)) ! include grid BCs to handle cyclic matrices
      end if
    end select
  end subroutine tridmatrix
  !
  subroutine bc_rhs(cbc,bc,dlc,dlf,c_or_f,rhs)
    implicit none
    character(len=1), intent(in), dimension(0:1) :: cbc
    real(rp), intent(in), dimension(0:1) :: bc
    real(rp), intent(in), dimension(0:1) :: dlc,dlf
    real(rp), intent(out), dimension(:,:,0:) :: rhs
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f -> face-centered
    real(rp), dimension(0:1) :: factor
    real(rp) :: sgn
    integer :: ibound
    !
    select case(c_or_f)
    case('c')
      do ibound = 0,1
        select case(cbc(ibound))
        case('P')
          factor(ibound) = 0.
        case('D')
          factor(ibound) = -2.*bc(ibound)
        case('N')
          if(ibound == 0) sgn =  1.
          if(ibound == 1) sgn = -1.
          factor(ibound) = sgn*dlc(ibound)*bc(ibound)
        end select
        rhs(:,:,ibound) = factor(ibound)/dlc(ibound)/dlf(ibound)
      end do
    case('f')
      do ibound = 0,1
        select case(cbc(ibound))
        case('P')
          factor(ibound) = 0.
        case('D')
          factor(ibound) = -bc(ibound)
        case('N')
          if(ibound == 0) sgn =  1.
          if(ibound == 1) sgn = -1.
          factor(ibound) = sgn*dlf(ibound)*bc(ibound)
        end select
        rhs(:,:,ibound) = factor(ibound)/dlc(ibound)/dlf(ibound)
      end do
    end select
  end subroutine bc_rhs
end module mod_initsolver
