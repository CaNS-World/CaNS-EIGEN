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
#if !defined(_OPENACC) || defined(_USE_HIP)
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
      call init_eigen_axis(ng(1),cbc(:,1),c_or_f(1),dxci_g,dxfi_g, &
                           lambdax_g,eigvecx_fwd,eigvecx_bwd)
    end if
    if(is_poisson_fft(2)) then
      call eigenvalues(ng(2),cbc(:,2),c_or_f(2),lambday_g)
      lambday_g(:) = lambday_g(:)*dli(2)**2
    else
      call init_eigen_axis(ng(2),cbc(:,2),c_or_f(2),dyci_g,dyfi_g, &
                           lambday_g,eigvecy_fwd,eigvecy_bwd)
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
  subroutine init_eigen_axis(n,bc,c_or_f,dci,dfi,lambda,fwd,bwd)
    ! Diagonalize the same active operator used by the staggered FD kernels.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: ieee_exceptions, only: ieee_status_type,ieee_get_status,ieee_set_status, &
                                            ieee_set_halting_mode,ieee_all
    implicit none
    integer, intent(in) :: n
    character(len=1), intent(in) :: bc(0:1),c_or_f
    real(rp), intent(in) :: dci(0:),dfi(0:)
    real(rp), intent(out) :: lambda(n),fwd(n,n),bwd(n,n)
    real(rp) :: a(n),b(n),c(n),mass(n),off(n),scale,tol
    real(rp), allocatable :: vectors(:,:),work(:)
    integer, allocatable :: iwork(:)
    integer :: m,i,j,info,lwork,liwork,izero
    type(ieee_status_type) :: fp_status
    ! Nonperiodic upper faces are prescribed or reconstructed by bounduvw.
    m = n-merge(1,0,c_or_f == 'f'.and.bc(1) /= 'P')
    if(m < 1) error stop 'ERROR: a nonperiodic face transform needs at least two grid points.'
    call tridmatrix(bc,n,dci,dfi,c_or_f,.false.,a,b,c)
    if(c_or_f == 'c') then
      mass = 1._rp/dfi(1:n)
    else
      mass = 1._rp/dci(1:n)
    end if
    if(any(mass <= 0._rp)) error stop 'ERROR: eigenproblem has nonpositive grid weights.'
    allocate(vectors(m,m),work(1),iwork(1))
    scale = maxval(abs(a)+abs(b)+abs(c))
    if(bc(0)//bc(1) == 'PP') then
      ! Accumulate, rather than overwrite: cyclic neighbors coincide at N=1/2.
      vectors = 0.
      do i=1,m
        vectors(i,i) = vectors(i,i)+b(i)
        j = modulo(i-2,m)+1
        vectors(i,j) = vectors(i,j)+a(i)*sqrt(mass(i)/mass(j))
        j = modulo(i,m)+1
        vectors(i,j) = vectors(i,j)+c(i)*sqrt(mass(i)/mass(j))
      end do
      ! LAPACK may use nonhalting IEEE arithmetic internally (also in workspace
      ! queries). Restore the caller's flags and traps after each library call.
      call ieee_get_status(fp_status)
      call ieee_set_halting_mode(ieee_all,.false.)
      call syevd('V','U',m,vectors,m,lambda,work,-1,iwork,-1,info)
      call ieee_set_status(fp_status)
      call check_lapack(info,'syevd workspace query')
      lwork = max(1,ceiling(work(1))); liwork = max(1,iwork(1))
      deallocate(work,iwork)
      allocate(work(lwork),iwork(liwork))
      call ieee_set_halting_mode(ieee_all,.false.)
      call syevd('V','U',m,vectors,m,lambda,work,lwork,iwork,liwork,info)
      call ieee_set_status(fp_status)
      call check_lapack(info,'syevd')
    else
      off = 0.
      do i=1,m-1
        off(i) = c(i)*sqrt(mass(i)/mass(i+1))
      end do
      call ieee_get_status(fp_status)
      call ieee_set_halting_mode(ieee_all,.false.)
      call stedc('I',m,b,off,vectors,m,work,-1,iwork,-1,info)
      call ieee_set_status(fp_status)
      call check_lapack(info,'stedc workspace query')
      lwork = max(1,ceiling(work(1))); liwork = max(1,iwork(1))
      deallocate(work,iwork)
      allocate(work(lwork),iwork(liwork))
      call ieee_set_halting_mode(ieee_all,.false.)
      call stedc('I',m,b,off,vectors,m,work,lwork,iwork,liwork,info)
      call ieee_set_status(fp_status)
      call check_lapack(info,'stedc')
      lambda(1:m) = b(1:m)
    end if
    if(.not.all(ieee_is_finite(lambda(1:m))).or..not.all(ieee_is_finite(vectors))) &
      error stop 'ERROR: LAPACK returned a nonfinite eigendecomposition.'
    if(bc(0)//bc(1) == 'PP'.or.bc(0)//bc(1) == 'NN') then
      ! LAPACK sorts modes differently from FFTs; identify the null mode by value.
      izero = minloc(abs(lambda(1:m)),dim=1)
      tol = 64._rp*epsilon(1._rp)*max(1._rp,scale)*m
      if(abs(lambda(izero)) > tol) error stop 'ERROR: eigenproblem lost its constant mode.'
      lambda(izero) = 0.
      vectors(:,izero) = sqrt(mass(1:m)/sum(mass(1:m)))
    end if
    fwd = 0.; bwd = 0.
    do j=1,m
      do i=1,m
        fwd(i,j) = vectors(j,i)*sqrt(mass(j))
        bwd(i,j) = vectors(i,j)/sqrt(mass(i))
      end do
    end do
    if(m < n) then
      ! Keep the allocated pencil extent; this decoupled slot is later overwritten by its BC.
      lambda(n) = 0.
      fwd(n,n) = 1.
      bwd(n,n) = 1.
    end if
  end subroutine init_eigen_axis
  !
  subroutine check_lapack(info,operation)
    integer, intent(in) :: info
    character(len=*), intent(in) :: operation
    if(info /= 0) then
      print*, 'ERROR: LAPACK ',operation,' failed, INFO = ',info
      error stop 'Eigenproblem initialization failed'
    end if
  end subroutine check_lapack
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
#if defined(_OPENACC)
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
    if(n == 1.and.bc(0)//bc(1) == 'PP') then
      a = 0.; b = 0.; c = 0.
      return
    end if
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
