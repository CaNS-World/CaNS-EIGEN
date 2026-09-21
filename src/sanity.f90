! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_sanity
  use, intrinsic :: iso_c_binding, only: C_PTR
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi
  use decomp_2d
  use mod_bound          , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv         , only: chkdiv
  use mod_common_mpi     , only: myid,ierr
  use mod_correc         , only: correc
  use mod_debug          , only: chk_helmholtz
  use mod_fft            , only: fftend
  use mod_fillps         , only: fillps
  use mod_initflow       , only: add_noise
  use mod_initmpi        , only: initmpi
  use mod_initsolver     , only: initsolver
  use mod_solve_helmholtz, only: solve_helmholtz
  use mod_param          , only: ipencil_axis,impdiff_mode,impdiff_z,impdiff_yz,impdiff_xyz,is_poisson_dtdma,small
  use mod_param          , only: is_poisson_fft_param => is_poisson_fft
#if !defined(_OPENACC)
  use mod_solver         , only: solver
#else
  use mod_solver_gpu     , only: solver => solver_gpu
#endif
  use mod_types
  implicit none
  private
  public test_sanity_input,test_sanity_grid,test_sanity_solver
  contains
  subroutine test_sanity_input(ng,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_forced)
    !
    ! performs some a priori checks of the input files before the calculation starts
    !
    implicit none
    integer , intent(in), dimension(3) :: ng
    integer , intent(in), dimension(2) :: dims
    logical , intent(in), dimension(3) :: stop_type
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
    real(rp)        , intent(in), dimension(0:1,3,3) :: bcvel
    real(rp)        , intent(in), dimension(0:1,3)   :: bcpre
    logical         , intent(in), dimension(3)       :: is_forced
    logical :: passed,passed_loc
    !
    call chk_dims(ng,dims,cbcvel,cbcpre,passed);                 if(.not.passed) call abortit
    call chk_stop_type(stop_type,passed);          if(.not.passed) call abortit
    call chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed); if(.not.passed) call abortit
    call chk_forcing(cbcpre,is_forced,passed);     if(.not.passed) call abortit
    if(is_poisson_dtdma) then
      passed_loc = ng(3)/dims(2) >= 2
      if(impdiff_mode /= 0 .and. (impdiff_mode /= impdiff_z .or. dims(2) > 1)) then
        if(cbcvel(1,3,3) /= 'P') passed_loc = passed_loc.and.(ng(3)/dims(2)-1 >= 2)
      end if
      if(myid == 0.and.(.not.passed_loc)) &
        print*, 'ERROR: DTDMA requires at least two active points per Z slab.'
      if(.not.passed_loc) call abortit
    end if
    if(impdiff_mode == impdiff_z .and. .not.(ipencil_axis == 3) .and. .not.is_poisson_dtdma) then
      if(dims(2) > 1) then
        if(myid == 0)  print*, 'Warning: a run with implicit Z diffusion (`impdiff_mode = 1`) is much more efficient &
                                       & when the flow is not decomposed along the Z direction.'
      end if
    end if
    if(impdiff_mode == impdiff_yz) then
      if(is_poisson_dtdma) then
        if(myid == 0) print*, 'ERROR: `impdiff_mode = 2` does not support `is_poisson_dtdma = T`.'; call abortit
      end if
      if(dims(2) /= 1 .or. .not.any(ipencil_axis == [2,3])) then
        if(myid == 0) print*, 'ERROR: `impdiff_mode = 2` requires `ipencil_axis = 2` or `3` and `dims(2) = 1`.'
        call abortit
      end if
      if(ipencil_axis == 3) then
        if(myid == 0) print*, 'Warning: `impdiff_mode = 2` may be more efficient with `ipencil_axis = 2`.'
      end if
    end if
    if(is_poisson_dtdma .and. (ipencil_axis == 3)) then
      if(myid == 0)  print*, 'ERROR: `is_poisson_dtdma = T` requires X/Y-aligned pencils.'; call abortit
    end if
  end subroutine test_sanity_input
  !
  subroutine test_sanity_grid(is_poisson_fft,dxf,dyf)
    implicit none
    logical, intent(in), dimension(2) :: is_poisson_fft
    real(rp), intent(in), dimension(0:) :: dxf,dyf
    logical, dimension(2) :: is_non_uniform_grid
    logical :: passed
    integer :: idir
    passed = .true.
    is_non_uniform_grid(:) = [any(dxf /= dxf(1)), any(dyf /= dyf(1))]
    do idir=1,2
      if(is_non_uniform_grid(idir) .and. is_poisson_fft(idir)) then
        if(myid == 0) print*, 'ERROR: FFT-based synthesis cannot be used with non-uniform grid.'
        if(myid == 0) print*, 'Check grid along direction: ', idir, '.'
        passed = .false.
      end if
    end do
    if(.not.passed) call abortit
  end subroutine test_sanity_grid
  !
  subroutine chk_stop_type(stop_type,passed)
    implicit none
    logical, intent(in), dimension(3) :: stop_type
    logical, intent(out) :: passed
    passed = .true.
    if(.not.any(stop_type(:))) then
      if(myid == 0) print*, 'ERROR: stopping criterion not chosen.'
      passed = .false.
    end if
  end subroutine chk_stop_type
  !
  subroutine chk_dims(ng,dims,cbcvel,cbcpre,passed)
    use mod_param, only: nscal,cbcscal
    implicit none
    integer, intent(in) :: ng(3),dims(2)
    character(len=1), intent(in) :: cbcvel(0:1,3,3),cbcpre(0:1,3)
    logical, intent(out) :: passed
    integer :: ii(2),idir,ivel,iscal
    logical :: periodic
    passed = all(ng >= 1).and.any(ipencil_axis == [1,2,3])
    if(.not.passed) then
      if(myid == 0) print*, 'ERROR: positive grid sizes and a valid pencil axis are required.'
      return
    end if
    ii = pack([1,2,3],[1,2,3] /= ipencil_axis)
    passed = all(dims <= ng(ii)).and.all(dims >= 1)
    if(myid == 0.and..not.passed) print*, 'ERROR: process grid exceeds the physical grid.'
    do idir=1,3
      if(ng(idir) /= 1) cycle
      periodic = cbcpre(0,idir)//cbcpre(1,idir) == 'PP'
      do ivel=1,3
        periodic = periodic.and.(cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel) == 'PP')
      end do
      do iscal=1,nscal
        periodic = periodic.and.(cbcscal(0,idir,iscal)//cbcscal(1,idir,iscal) == 'PP')
      end do
      passed = passed.and.periodic
      if(myid == 0.and..not.periodic) print*, 'ERROR: singleton directions require periodic boundary conditions.'
    end do
  end subroutine chk_dims
  !
  subroutine chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed)
    use mod_param, only: nscal,cbcscal
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3  ) :: cbcpre
    real(rp)        , intent(in), dimension(0:1,3,3) :: bcvel
    real(rp)        , intent(in), dimension(0:1,3  ) :: bcpre
    logical         , intent(out) :: passed
    character(len=2) :: bc01v,bc01p
    integer :: ivel,idir,iscal
    logical :: passed_loc
    passed = .true.
    !
    ! check validity of pressure and velocity BCs
    !
    passed_loc = .true.
    do ivel = 1,3
      do idir=1,3
        bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
        passed_loc = passed_loc.and.( (bc01v == 'PP').or. &
                                      (bc01v == 'ND').or. &
                                      (bc01v == 'DN').or. &
                                      (bc01v == 'NN').or. &
                                      (bc01v == 'DD') )
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity BCs not valid.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01p == 'PP').or. &
                                    (bc01p == 'ND').or. &
                                    (bc01p == 'DN').or. &
                                    (bc01p == 'NN').or. &
                                    (bc01p == 'DD') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: pressure BCs not valid.'
    passed = passed.and.passed_loc
    !
    ! check that all variables have the same periodic directions
    !
    passed_loc = .true.
    do idir=1,3
      do ivel=1,3
        passed_loc = passed_loc.and.all((cbcvel(:,idir,ivel) == 'P').eqv.(cbcpre(:,idir) == 'P'))
      end do
      do iscal=1,nscal
        passed_loc = passed_loc.and.all((cbcscal(:,idir,iscal) == 'P').eqv.(cbcpre(:,idir) == 'P'))
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: velocity and scalar periodic BCs must match pressure BCs in every direction.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      ivel = idir
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01v == 'PP'.and.bc01p == 'PP').or. &
                                    (bc01v == 'ND'.and.bc01p == 'DN').or. &
                                    (bc01v == 'DN'.and.bc01p == 'ND').or. &
                                    (bc01v == 'DD'.and.bc01p == 'NN').or. &
                                    (bc01v == 'NN'.and.bc01p == 'DD') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity and pressure BCs not compatible.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,2
      passed_loc = passed_loc.and.((bcpre(0,idir) == 0.).and.(bcpre(1,idir) == 0.))
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: pressure BCs in directions x and y must be homogeneous (value = 0.).'
    passed = passed.and.passed_loc
    if(impdiff_mode == impdiff_yz .or. impdiff_mode == impdiff_xyz) then
      passed_loc = .true.
      do ivel=1,3
        do idir=1,2
          if((impdiff_mode == impdiff_xyz.or.idir == 2).and.is_poisson_fft_param(idir)) then
            passed_loc = passed_loc.and.(bcvel(0,idir,ivel) == 0.).and.(bcvel(1,idir,ivel) == 0.)
          end if
        end do
      end do
      if(myid == 0.and..not.passed_loc) &
        print*, 'ERROR: velocity BCs with FFT-based implicit diffusion in directions x/y must be homogeneous (value = 0.).'
      passed = passed.and.passed_loc
    end if
    block
      character(len=2) :: pair
      passed_loc = .true.
      do iscal=1,nscal
        do idir=1,3
          pair = cbcscal(0,idir,iscal)//cbcscal(1,idir,iscal)
          passed_loc = passed_loc.and.any(pair == ['PP','DD','NN','DN','ND'])
        end do
      end do
      if(myid == 0.and..not.passed_loc) print*, 'ERROR: scalar boundary conditions are not valid.'
      passed = passed.and.passed_loc
    end block
  end subroutine chk_bc
  !
  subroutine chk_forcing(cbcpre,is_forced,passed)
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre
    logical         , intent(in), dimension(3) :: is_forced
    logical         , intent(out) :: passed
    integer :: idir
    passed = .true.
    !
    ! 1) check for compatibility between pressure BCs and flow forcing
    !
    do idir=1,3
      if(is_forced(idir)) then
        passed = passed.and.(cbcpre(0,idir)//cbcpre(1,idir) == 'PP')
      end if
    end do
    if(myid == 0.and.(.not.passed)) &
    print*, 'ERROR: Flow cannot be forced in a non-periodic direction; check the BCs and is_forced in `input.nml`.'
  end subroutine chk_forcing
  !
  subroutine test_sanity_solver(ng,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,is_poisson_fft,l, &
                                dxc,dxf,dyc,dyf,dzc,dzf,dxci,dxfi,dyci,dyfi,dzci,dzfi,dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g, &
                                nb,is_bound,cbcvel,cbcpre,bcvel,bcpre)
#if defined(_OPENACC)
    use mod_workspaces     , only: set_cufft_wspace
    use mod_common_cudecomp, only: istream_acc_queue_1
#endif
    use mod_param, only: nscal,cbcscal,bcscal
    implicit none
    integer , intent(in), dimension(3) :: ng,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
    logical , intent(in), dimension(2) :: is_poisson_fft
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(in), dimension(0:) :: dxc,dxf,dyc,dyf,dzc,dzf, &
                                           dxci,dxfi,dyci,dyfi,dzci,dzfi, &
                                           dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
    real(rp), intent(in), dimension(0:1,3,3)         :: bcvel
    real(rp), intent(in), dimension(0:1,3)           :: bcpre
    real(rp), allocatable, target, dimension(:,:,:) :: u,v,w,p,phi
    real(rp), pointer, contiguous :: field(:,:,:)
#if !defined(_OPENACC) || defined(_USE_HIP)
    type(C_PTR), dimension(2,2) :: arrplan
#else
    integer    , dimension(2,2) :: arrplan
#endif
    real(rp), dimension(2) :: normfft
    real(rp), allocatable, dimension(:) :: lambdax_g,lambday_g
    real(rp), allocatable, dimension(:,:) :: lambdaxy,eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd
    real(rp), allocatable, dimension(:) :: a,b,c
    real(rp), allocatable, dimension(:,:,:) :: rhsbx,rhsby,rhsbz
    real(rp) :: dt,dti,alpha,alphai,div_initial,relative_residual,op_norm,field_norm(2)
    character(len=1) :: field_bc(0:1,3),center(3)
    real(rp) :: field_values(0:1,3)
    logical :: implicit_dir(3)
    integer :: icomponent
    real(rp) :: divtot,divmax,restot,resmax
    integer :: i,j,k
    logical :: passed,passed_loc
    passed = .true.
    !$acc wait
    allocate(u(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
             v(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
             w(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
             p(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
             phi(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
             lambdax_g(ng(1)),lambday_g(ng(2)),lambdaxy(n_z(1),n_z(2)), &
             eigvecx_fwd(ng(1),ng(1)),eigvecx_bwd(ng(1),ng(1)), &
             eigvecy_fwd(ng(2),ng(2)),eigvecy_bwd(ng(2),ng(2)), &
             a(n_z(3)),b(n_z(3)),c(n_z(3)), &
             rhsbx(n(2),n(3),0:1), &
             rhsby(n(1),n(3),0:1), &
             rhsbz(n(1),n(2),0:1))
    !$acc enter data copyin(n,n_z)
    !
    ! initialize velocity below with some random noise
    !
    u(:,:,:) = 0.
    v(:,:,:) = 0.
    w(:,:,:) = 0.
    p(:,:,:) = 0.
    phi(:,:,:) = 0.
    call add_noise(ng,lo,123,.5_rp,u(1:n(1),1:n(2),1:n(3)))
    call add_noise(ng,lo,456,.5_rp,v(1:n(1),1:n(2),1:n(3)))
    call add_noise(ng,lo,789,.5_rp,w(1:n(1),1:n(2),1:n(3)))
    !$acc enter data copyin(u,v,w,p,phi)
    !
    ! test pressure correction
    !
    call initsolver(is_poisson_fft,ng,n_x_fft,n_y_fft,lo_z,hi_z,dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g, &
                    cbcpre,bcpre(:,:),lambdax_g,lambday_g,lambdaxy,eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd, &
                    ['c','c','c'],a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
    !$acc enter data copyin(lambday_g,lambdaxy,eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd,a,b,c,rhsbx,rhsby,rhsbz)
#if defined(_OPENACC)
    call set_cufft_wspace(pack(arrplan,.true.),istream_acc_queue_1)
#endif
    dt  = acos(-1.) ! value is irrelevant
    dti = dt**(-1)
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,dxc,dxf,dyc,dyf,dzc,dzf,u,v,w,.false.,-1)
    call chkdiv(lo,hi,l,dxfi,dyfi,dzfi,u,v,w,divtot,div_initial)
    call fillps(n,dxfi,dyfi,dzfi,dti,u,v,w,p)
    call updt_rhs_b(['c','c','c'],cbcpre,n,is_bound,rhsbx,rhsby,rhsbz,p)
    call solver(n,ng,is_poisson_fft,arrplan,product(normfft(:)),lambdaxy,eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd,a,b,c, &
                cbcpre,['c','c','c'],p)
    call boundp(cbcpre,n,bcpre,nb,is_bound,dxc,dyc,dzc,p)
    call correc(n,dxci,dyci,dzci,dt,p,u,v,w)
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,dxc,dxf,dyc,dyf,dzc,dzf,u,v,w,.true.,1,0)
    call chkdiv(lo,hi,l,dxfi,dyfi,dzfi,u,v,w,divtot,divmax)
    !$acc update self(u,v,w,p)
    passed_loc = all(ieee_is_finite(u)).and.all(ieee_is_finite(v)).and. &
                 all(ieee_is_finite(w)).and.all(ieee_is_finite(p))
    call MPI_ALLREDUCE(MPI_IN_PLACE,passed_loc,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
    if(.not.passed_loc) error stop 'ERROR: nonfinite pressure-correction result.'
    relative_residual = divmax/max(1._rp,div_initial)
    if(myid == 0) print*, 'CHECK pressure relative divergence:',relative_residual
    passed_loc = relative_residual < small
    if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: Pressure correction: Divergence is too large, with maximum = ', divmax
    passed = passed.and.passed_loc
    call fftend(arrplan)
    !
    implicit_dir = [ impdiff_mode == impdiff_xyz, &
                    (impdiff_mode == impdiff_yz).or.(impdiff_mode == impdiff_xyz), &
                     impdiff_mode /= 0 ]
    if(impdiff_mode /= 0) then
      alpha = -acos(-1._rp)
      alphai = 1._rp/alpha
      do icomponent=1,3+nscal
        center = 'c'
        if(icomponent <= 3) then
          center(icomponent) = 'f'
          field_bc     = cbcvel(:,:,icomponent)
          field_values =  bcvel(:,:,icomponent)
        else
          field_bc     = cbcscal(:,:,icomponent-3)
          field_values =  bcscal(:,:,icomponent-3)
        end if
        select case(icomponent)
        case(1)
          field => u
        case(2)
          field => v
        case(3)
          field => w
        case default
          field => phi
        end select
        field(:,:,:) = 0.
        call add_noise(ng,lo,123*icomponent,.5_rp,field(1:n(1),1:n(2),1:n(3)))
        !$acc update device(field)
        call initsolver(is_poisson_fft,ng,n_x_fft,n_y_fft,lo_z,hi_z,dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g, &
                        field_bc,field_values,lambdax_g,lambday_g,lambdaxy, &
                        eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd,center,a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
        !$acc update device(lambday_g,lambdaxy,eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd,a,b,c,rhsbx,rhsby,rhsbz)
        if(impdiff_mode /= impdiff_xyz) call fftend(arrplan,1)
        if(impdiff_mode == impdiff_z  ) call fftend(arrplan,2)
#if defined(_OPENACC)
        call set_cufft_wspace(pack(arrplan,.true.),istream_acc_queue_1)
#endif
        if(icomponent <= 3) then
          call bounduvw(cbcvel,n,bcvel,nb,is_bound,dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
        else
          call boundp(field_bc,n,field_values,nb,is_bound,dxc,dyc,dzc,field)
        end if
        !$acc wait
        !$acc parallel loop collapse(3) default(present)
        !$OMP parallel do   collapse(3) default(shared)
        do k=0,n(3)+1
          do j=0,n(2)+1
            do i=0,n(1)+1
              p(i,j,k) = field(i,j,k)
            end do
          end do
        end do
        call solve_helmholtz(n,ng,hi,is_poisson_fft,arrplan,normfft,alpha,lambday_g,lambdaxy, &
                             eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd, &
                             a,b,c,rhsbx,rhsby,rhsbz,is_bound,field_bc,center,field)
        !$acc wait
        call fftend(arrplan)
        if(icomponent <= 3) then
          call bounduvw(cbcvel,n,bcvel,nb,is_bound,dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
        else
          call boundp(field_bc,n,field_values,nb,is_bound,dxc,dyc,dzc,field)
        end if
        call chk_helmholtz(lo,hi,l,dxci,dxfi,dyci,dyfi,dzci,dzfi,alphai,p,field,field_bc,is_bound,center, &
                           restot,resmax,implicit_dir)
        !$acc update self(field,p)
        passed_loc = all(ieee_is_finite(field)).and.all(ieee_is_finite(p))
        call MPI_ALLREDUCE(MPI_IN_PLACE,passed_loc,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        if(.not.passed_loc) error stop 'ERROR: nonfinite Helmholtz result.'
        field_norm = [maxval(abs(p)),maxval(abs(field))]
        call MPI_ALLREDUCE(MPI_IN_PLACE,field_norm,2,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
        op_norm = abs(alphai)
        if(implicit_dir(1)) op_norm = op_norm+4._rp*maxval(dxci_g)*maxval(dxfi_g)
        if(implicit_dir(2)) op_norm = op_norm+4._rp*maxval(dyci_g)*maxval(dyfi_g)
        if(implicit_dir(3)) op_norm = op_norm+4._rp*maxval(dzci_g)*maxval(dzfi_g)
        relative_residual = resmax/max(1._rp,abs(alphai)*field_norm(1)+op_norm*field_norm(2))
        if(myid == 0) print*, 'CHECK Helmholtz mode/field/residual:',impdiff_mode,icomponent,relative_residual
        passed_loc = relative_residual < small
        if(myid == 0.and..not.passed_loc) print*, 'ERROR: wrong solution of Helmholtz equation.'
        passed = passed.and.passed_loc
      end do
    end if
    !$acc exit data delete(u,v,w,p,phi,lambday_g,lambdaxy,eigvecx_fwd,eigvecx_bwd,eigvecy_fwd,eigvecy_bwd,a,b,c,rhsbx,rhsby,rhsbz)
    if(.not.passed) then
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      error stop
    end if
  end subroutine test_sanity_solver
  !
  subroutine abortit
    implicit none
    if(myid == 0) print*, ''
    if(myid == 0) print*, '*** Simulation aborted due to errors in the input file ***'
    if(myid == 0) print*, '    check `input.nml`.'
    call decomp_2d_finalize
    call MPI_FINALIZE(ierr)
    error stop
  end subroutine abortit
end module mod_sanity
