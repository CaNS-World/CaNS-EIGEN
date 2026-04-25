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
                        lambdaxy,eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd,c_or_f,a,b,c,arrplan,normfft, &
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
    real(rp), intent(out) :: normfft
    real(rp), dimension(2)         :: dl,dli
    real(rp), dimension(0:ng(1)+1) :: dxc_g,dxf_g
    real(rp), dimension(0:ng(2)+1) :: dyc_g,dyf_g
    real(rp), dimension(0:ng(3)+1) :: dzc_g,dzf_g
    integer :: q,i,j
    real(rp), dimension(ng(1))      :: lambdax
    real(rp), dimension(ng(2))      :: lambday
    real(rp), dimension(ng(1))      :: ax_g,bx_g,cx_g
    real(rp), dimension(ng(2))      :: ay_g,by_g,cy_g
    real(rp), dimension(ng(3))      :: az_g,bz_g,cz_g
    real(rp), allocatable, dimension(:,:) :: eigvecs
    real(rp), allocatable, dimension(:) :: work
    integer , allocatable, dimension(:) :: iwork
    integer :: wsize,iwsize,info
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
    !         Q = M^{1/2} U; Q^{-1} = M^{-1/2} U^t => T = Q Lambda Q^{-1},
    !
    !       where, in the Poisson solver:
    !         - Q^{-1} provides the forward  transforms to be applied along x/y (`eigvecx/y_fwd` below), and
    !         - Q      provides the backward transforms to be applied along x/y (`eigvecx/y_bwd` below).
    !
    if(is_poisson_fft(1)) then
      call eigenvalues(ng(1),cbc(:,1),c_or_f(1),lambdax)
      lambdax(:) = lambdax(:)*dli(1)**2
    else
      !
      ! numerical eigendecomposition
      !
      q = merge(1,0,c_or_f(1) == 'f'.and.cbc(1,1) == 'D')
      call tridmatrix(cbc(:,1),ng(1),dxci_g,dxfi_g,c_or_f(1),.true.,ax_g,bx_g,cx_g)
      allocate(eigvecs(ng(1),ng(1)),work(1),iwork(1))
      if(cbc(0,1)//cbc(1,1) /= 'PP') then
        !
        ! non-periodic BCs: simple symmetric tridiagonal matrix
        !
        call stedc('I',ng(1)-q,bx_g,cx_g,eigvecs(1:ng(1)-q,1:ng(1)-q),ng(1)-q,work,-1,iwork,-1,info) ! workspace size query
        wsize = int(work(1),kind(wsize)); iwsize = iwork(1)
        deallocate(work,iwork)
        allocate(work(wsize),iwork(iwsize))
        call stedc('I',ng(1)-q,bx_g,cx_g,eigvecs(1:ng(1)-q,1:ng(1)-q),ng(1)-q,work,wsize,iwork,iwsize,info)
        lambdax(:) = bx_g(:)
      else
        !
        ! periodic BCs: define full cyclic symmetric tridiagonal matrix (upper diagonal)
        !
        eigvecs(:,:) = 0.
        do i=1,ng(1)-q
          eigvecs(i,  i) = bx_g(i)
        end do
        do i=1,ng(1)-q-1
          eigvecs(i,i+1) = cx_g(i)
        end do
        eigvecs(1,ng(1)-q) = cx_g(ng(1)-q) ! == ax_g(1) cyclic BC (but note that ax_g array is not symmetrized in `tridmatrix`)
        call syevd('V','U',ng(1)-q,eigvecs(1:ng(1)-q,1:ng(1)-q),ng(1)-q,lambdax(1:ng(1)-q),work,-1,iwork,-1,info) ! workspace size query
        wsize = int(work(1),kind(wsize)); iwsize = iwork(1)
        deallocate(work,iwork)
        allocate(work(wsize),iwork(iwsize))
        call syevd('V','U',ng(1)-q,eigvecs(1:ng(1)-q,1:ng(1)-q),ng(1)-q,lambdax(1:ng(1)-q),work,wsize,iwork,iwsize,info)
      end if
      !
      ! compute generalized eigenvectors
      !
      select case(c_or_f(1))
      case('c')
        do j=1,ng(1)
          do i=1,ng(1)
            eigvecx_fwd(i,j) = eigvecs(j,i)*sqrt(dxf_g(j))
            eigvecx_bwd(i,j) = sqrt(dxf_g(i))**(-1)*eigvecs(i,j)
          end do
        end do
      case('f')
        if(q == 1) then ! set trivial equation for the boundary point
          lambdax(ng(1)) = 0.
          eigvecs(ng(1),:    ) = 0.
          eigvecs(:    ,ng(1)) = 0.
          eigvecs(ng(1),ng(1)) = 1.
        end if
        do j=1,ng(1)
          do i=1,ng(1)
            eigvecx_fwd(i,j) = eigvecs(j,i)*sqrt(dxc_g(j))
            eigvecx_bwd(i,j) = sqrt(dxc_g(i))**(-1)*eigvecs(i,j)
          end do
        end do
      end select
      deallocate(work,iwork,eigvecs)
    end if
    if(is_poisson_fft(2)) then
      call eigenvalues(ng(2),cbc(:,2),c_or_f(2),lambday)
      lambday(:) = lambday(:)*dli(2)**2
    else
      !
      ! numerical eigendecomposition
      !
      q = merge(1,0,c_or_f(2) == 'f'.and.cbc(1,2) == 'D')
      call tridmatrix(cbc(:,2),ng(2),dyci_g,dyfi_g,c_or_f(2),.true.,ay_g,by_g,cy_g)
      allocate(eigvecs(ng(2),ng(2)),work(1),iwork(1))
      if(cbc(0,2)//cbc(1,2) /= 'PP') then
        !
        ! non-periodic BCs: simple symmetric tridiagonal matrix
        !
        call stedc('I',ng(2)-q,by_g,cy_g,eigvecs(1:ng(2)-q,1:ng(2)-q),ng(2)-q,work,-1,iwork,-1,info) ! workspace size query
        wsize = int(work(1),kind(wsize)); iwsize = iwork(1)
        deallocate(work,iwork)
        allocate(work(wsize),iwork(iwsize))
        call stedc('I',ng(2)-q,by_g,cy_g,eigvecs(1:ng(2)-q,1:ng(2)-q),ng(2)-q,work,wsize,iwork,iwsize,info)
        lambday(:) = by_g(:)
      else
        !
        ! periodic BCs: define full cyclic symmetric tridiagonal matrix (upper diagonal)
        !
        eigvecs(:,:) = 0.
        do j=1,ng(2)-q
          eigvecs(j,  j) = by_g(j)
        end do
        do j=1,ng(2)-q-1
          eigvecs(j,j+1) = cy_g(j)
        end do
        eigvecs(1,ng(2)-q) = cy_g(ng(2)-q) ! == ay_g(1) cyclic BC (but note that ay_g array is not symmetrized in `tridmatrix`)
        call syevd('V','U',ng(2)-q,eigvecs(1:ng(2)-q,1:ng(2)-q),ng(2)-q,lambday(1:ng(2)-q),work,-1,iwork,-1,info) ! workspace size query
        wsize = int(work(1),kind(wsize)); iwsize = iwork(1)
        deallocate(work,iwork)
        allocate(work(wsize),iwork(iwsize))
        call syevd('V','U',ng(2)-q,eigvecs(1:ng(2)-q,1:ng(2)-q),ng(2)-q,lambday(1:ng(2)-q),work,wsize,iwork,iwsize,info)
      end if
      !
      ! compute generalized eigenvectors
      !
      select case(c_or_f(2))
      case('c')
        do j=1,ng(2)
          do i=1,ng(2)
            eigvecy_fwd(i,j) = eigvecs(j,i)*sqrt(dyf_g(j))
            eigvecy_bwd(i,j) = sqrt(dyf_g(i))**(-1)*eigvecs(i,j)
          end do
        end do
      case('f')
        if(q == 1) then ! set trivial equation for the boundary point
          lambday(ng(2)) = 0.
          eigvecs(ng(2),:    ) = 0.
          eigvecs(:    ,ng(2)) = 0.
          eigvecs(ng(2),ng(2)) = 1.
        end if
        do j=1,ng(2)
          do i=1,ng(2)
            eigvecy_fwd(i,j) = eigvecs(j,i)*sqrt(dyc_g(j))
            eigvecy_bwd(i,j) = sqrt(dyc_g(i))**(-1)*eigvecs(i,j)
          end do
        end do
      end select
      deallocate(work,iwork,eigvecs)
    end if
    !
    ! add eigenvalues
    !
    do j=lo_z(2),hi_z(2)
      do i=lo_z(1),hi_z(1)
        lambdaxy(i,j) = lambdax(i)+lambday(j)
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
        iswap(2) = nh+(1-mod(n,2))
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
        do l=1,n
          lambda(l)   = -2.*(1.-cos((l-1  )*pi/(1.*(n-1+1))))
        end do
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
      do l=1,n
        lambda(l)     = -2.*(1.-cos((2*l-1)*pi/(2.*n)))
      end do
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
      if(bc(1) == 'N') b(n) = b(n) + factor(1)*c(n)
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
