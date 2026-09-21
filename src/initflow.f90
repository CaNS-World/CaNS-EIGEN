! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_initflow
  use mpi
  use mod_common_mpi, only: ierr,myid
  use mod_param     , only: pi
  use mod_types
  implicit none
  private
  public initflow,initscal,add_noise
  contains
  subroutine initflow(inivel,cbcvel,bcvel,ng,lo,l,xc,xf,yc,yf,zc,zf,dxc,dxf,dyc,dyf,dzc,dzf, &
                      visc,is_forced,velf,bforce,is_wallturb,u,v,w,p)
    !
    ! computes initial conditions for the velocity field
    !
    implicit none
    character(len=*), intent(in) :: inivel
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    real(rp), intent(in), dimension(0:1,3,3) :: bcvel
    integer , intent(in), dimension(3) :: ng,lo
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(in), dimension(0:) :: xc,xf,yc,yf,zc,zf,dxc,dxf,dyc,dyf,dzc,dzf
    real(rp), intent(in)               :: visc
    logical , intent(in), dimension(3) :: is_forced
    real(rp), intent(in), dimension(3) :: velf,bforce
    logical , intent(in)               :: is_wallturb
    real(rp), target, dimension(0:,0:,0:), intent(inout) :: u,v,w,p
    real(rp), pointer, dimension(:,:,:) :: us,ut,un
    real(rp), allocatable, dimension(:) :: u1d_n,u1d_t,rn,rt,rs,rft,rfn
    character(len=len(inivel)) :: ini
    integer :: i,j,k,idir_s,idir_t,idir_n,ih,ia,ib,isuffix
    logical :: is_noise,is_mean,is_pair,is_duct
    real(rp) :: xxc,yyc,zzc,xxf,yyf,zzf
    real(rp) :: uref,lref,lprof
    real(rp) :: ubulk,reb,retau
    real(rp), dimension(3) :: xyz_c,xyz_f,vel
    integer, dimension(3) :: n,ii
    !
    n(:) = shape(p) - 2*1
    ini = inivel
    isuffix = len_trim(ini)
    idir_s = 1
    if(isuffix == 5) then
      if((ini(isuffix-1:isuffix-1) == '-').and.(.not.any(ini(1:3) == ['tgv','ant']))) then
        select case(ini(isuffix:isuffix))
        case('x','y','z')
          idir_s = index('xyz',ini(isuffix:isuffix))
          ini = ini(:isuffix-2)
        end select
      end if
    end if
    !
    ! streamwise, transverse and normal directions for the profile
    !
    select case(idir_s)
    case(1)
      idir_t = 2; idir_n = 3
      us => u; ut => v; un => w
    case(2)
      idir_t = 1; idir_n = 3
      us => v; ut => u; un => w
    case(3)
      idir_t = 2; idir_n = 1
      us => w; ut => v; un => u
    end select
    allocate(u1d_n(n(idir_n)),u1d_t(n(idir_t)),rn(0:n(idir_n)+1),rt(0:n(idir_t)+1))
    select case(idir_n)
    case(1)
      rn(:) = xc(0:n(1)+1)
    case(2)
      rn(:) = yc(0:n(2)+1)
    case(3)
      rn(:) = zc(0:n(3)+1)
    end select
    select case(idir_t)
    case(1)
      rt(:) = xc(0:n(1)+1)
    case(2)
      rt(:) = yc(0:n(2)+1)
    case(3)
      rt(:) = zc(0:n(3)+1)
    end select
    u1d_t(:) = 1.
    lprof = l(idir_n)
    if(any(trim(ini) == ['hcp','hcl','hdc'])) lprof = 2.*lprof
    is_duct = any(trim(ini) == ['poi','log']) .and. &
              (cbcvel(0,idir_t,idir_s)//cbcvel(1,idir_t,idir_s) == 'DD')
    is_noise = .false.
    is_mean  = .false.
    is_pair  = .false.
    uref  = 1.
    ubulk = uref
    if(is_forced(idir_s)) ubulk = velf(idir_s)
    select case(trim(ini))
    case('cou')
      call couette(   n(idir_n),rn/lprof,1._rp,u1d_n)
      u1d_n(:) = u1d_n(:) + 0.5 ! from 1 to 0
      u1d_n(:) = bcvel(0,idir_n,idir_s)*(u1d_n(:)) + bcvel(1,idir_n,idir_s)*(1.-u1d_n(:))
      uref = abs(bcvel(1,idir_n,idir_s)-bcvel(0,idir_n,idir_s))
    case('poi','hcp')
      call poiseuille(n(idir_n),rn/lprof,ubulk,u1d_n)
      if(is_duct) call poiseuille(n(idir_t),rt/l(idir_t),1._rp,u1d_t)
      is_mean = .true.
    case('tbl')
      call temporal_bl(n(idir_n),rn,1._rp,visc,uref,u1d_n)
      is_noise = .true.
    case('iop') ! reversed 'poi'
      !
      ! convective reference frame moving with velocity `ubulk`;
      ! walls have negative velocity equal to `ubulk` in the laboratory frame
      !
      ubulk = 0.5*abs(bcvel(0,idir_n,idir_s)+bcvel(1,idir_n,idir_s))
      call poiseuille(n(idir_n),rn/lprof,ubulk,u1d_n)
      u1d_n(:) = u1d_n(:) - ubulk
      is_mean = .true.
    case('zer')
      u1d_n(:) = 0.
    case('uni')
      u1d_n(:) = uref
    case('log','hcl')
      reb = ubulk*lprof/visc
      call log_profile(n(idir_n),rn/lprof,reb,u1d_n)
      if(is_duct) then
        reb = ubulk*l(idir_t)/visc
        call log_profile(n(idir_t),rt/l(idir_t),reb,u1d_t)
      end if
      is_noise = .true.
      is_mean = .true.
    case('tgv')
      do k=1,n(3)
        zzc = zc(k)/l(3)*2.*pi
        do j=1,n(2)
          yyc = yc(j)/l(2)*2.*pi
          yyf = yf(j)/l(2)*2.*pi
          do i=1,n(1)
            xxc = xc(i)/l(1)*2.*pi
            xxf = xf(i)/l(1)*2.*pi
            u(i,j,k) =  sin(xxf)*cos(yyc)*cos(zzc)*uref
            v(i,j,k) = -cos(xxc)*sin(yyf)*cos(zzc)*uref
            w(i,j,k) = 0.
            p(i,j,k) = 0.!(cos(2.*xxc)+cos(2.*yyc))*(cos(2.*zzc)+2.)/16.*uref**2
          end do
        end do
      end do
    case('tgv-2d-x','tgv-2d-y','tgv-2d-z')
      ih = index('xyz',ini(isuffix:isuffix))
      ia = mod(ih  ,3)+1
      ib = mod(ih+1,3)+1
      do k=1,n(3)
        xyz_c(3) = zc(k)
        xyz_f(3) = zf(k)
        do j=1,n(2)
          xyz_c(2) = yc(j)
          xyz_f(2) = yf(j)
          do i=1,n(1)
            xyz_c(1) = xc(i)
            xyz_f(1) = xf(i)
            vel(:) = 0.
            vel(ia) =  cos(xyz_f(ia))*sin(xyz_c(ib))*uref
            vel(ib) = -sin(xyz_c(ia))*cos(xyz_f(ib))*uref
            u(i,j,k) = vel(1)
            v(i,j,k) = vel(2)
            w(i,j,k) = vel(3)
            p(i,j,k) = -(cos(2.*xyz_c(ia))+cos(2.*xyz_c(ib)))/4.*uref**2
          end do
        end do
      end do
    case('ant')
      !
      ! see M. Antuono, JFM 890, A23 (2020)
      !
      do k=1,n(3)
        zzc = zc(k)/l(3)*2.*pi+0.5*pi
        zzf = zf(k)/l(3)*2.*pi+0.5*pi
        do j=1,n(2)
          yyc = yc(j)/l(2)*2.*pi+0.5*pi
          yyf = yf(j)/l(2)*2.*pi+0.5*pi
          do i=1,n(1)
            xxc = xc(i)/l(1)*2.*pi+0.5*pi
            xxf = xf(i)/l(1)*2.*pi+0.5*pi
            u(i,j,k) = (4.*sqrt(2.)/3./sqrt(3.))*(sin(xxf-5.*pi/6.)*cos(yyc-1.*pi/6.)*sin(zzc         ) - &
                                                  sin(xxf-1.*pi/6.)*sin(yyc         )*cos(zzc-5.*pi/6.))*uref
            v(i,j,k) = (4.*sqrt(2.)/3./sqrt(3.))*(sin(xxc         )*sin(yyf-5.*pi/6.)*sin(zzc-1.*pi/6.) - &
                                                  cos(xxc-5.*pi/6.)*sin(yyf-1.*pi/6.)*sin(zzc         ))*uref
            w(i,j,k) = (4.*sqrt(2.)/3./sqrt(3.))*(cos(xxc-1.*pi/6.)*sin(yyc         )*sin(zzf-5.*pi/6.) - &
                                                  sin(xxc         )*cos(yyc-5.*pi/6.)*sin(zzf-1.*pi/6.))*uref
            p(i,j,k) = -(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)/2.
          end do
        end do
      end do
    case('pdc','hdc')
      lref = lprof/2.
      if(is_wallturb) then ! turbulent flow
        uref  = (bforce(idir_s)*lref)**(0.5) ! utau = sqrt(-dpdx*h)
        retau = uref*lref/visc
        reb   = (retau/.09)**(1./.88)
        ubulk = reb*visc/(2*lref)
      else                 ! laminar flow
        ubulk = (bforce(idir_s)*lref**2/(3.*visc))
      end if
      call poiseuille(n(idir_n),rn/lprof,ubulk,u1d_n)
      is_mean = .true.
    case default
      if(myid == 0) print*, 'ERROR: invalid name for initial velocity field'
      if(myid == 0) print*, ''
      if(myid == 0) print*, '*** Simulation aborted due to errors in the case file ***'
      if(myid == 0) print*, '    check INFO_INPUT.md'
      call MPI_FINALIZE(ierr)
      error stop
    end select
    if(.not.any(trim(ini(1:3)) == ['tgv','ant'])) then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            ii = [i,j,k]
            us(i,j,k) = u1d_t(ii(idir_t))*u1d_n(ii(idir_n))
            ut(i,j,k) = 0.
            un(i,j,k) = 0.
            p(i,j,k) = 0.
          end do
        end do
      end do
    end if
    if(is_noise) then
      call add_noise(ng,lo,123,.05_rp,us(1:n(1),1:n(2),1:n(3)))
      call add_noise(ng,lo,456,.05_rp,ut(1:n(1),1:n(2),1:n(3)))
      call add_noise(ng,lo,789,.05_rp,un(1:n(1),1:n(2),1:n(3)))
    end if
    if(is_mean) then
      if(trim(ini) /= 'iop') then
        select case(idir_s)
        case(1)
          call set_mean(n,ubulk,dxc,dyf,dzf,l,us(1:n(1),1:n(2),1:n(3)))
        case(2)
          call set_mean(n,ubulk,dxf,dyc,dzf,l,us(1:n(1),1:n(2),1:n(3)))
        case(3)
          call set_mean(n,ubulk,dxf,dyf,dzc,l,us(1:n(1),1:n(2),1:n(3)))
        end select
      end if
    end if
    if(is_wallturb) is_pair = .true.
    if(is_pair) then
      allocate(rs(0:n(idir_s)+1),rft(0:n(idir_t)+1),rfn(0:n(idir_n)+1))
      select case(idir_s)
      case(1)
        rs(:) = xc(0:n(1)+1)
      case(2)
        rs(:) = yc(0:n(2)+1)
      case(3)
        rs(:) = zc(0:n(3)+1)
      end select
      select case(idir_t)
      case(1)
        rft(:) = xf(0:n(1)+1)
      case(2)
        rft(:) = yf(0:n(2)+1)
      case(3)
        rft(:) = zf(0:n(3)+1)
      end select
      select case(idir_n)
      case(1)
        rfn(:) = xf(0:n(1)+1)
      case(2)
        rfn(:) = yf(0:n(2)+1)
      case(3)
        rfn(:) = zf(0:n(3)+1)
      end select
      if(.false.) then
        !
        ! initialize a streamwise vortex pair for a fast transition
        ! to turbulence in a pressure-driven channel:
        !        psi(x,y,z)  = f(z)*g(x,y), with
        !        f(z)        = (1-z**2)**2, and
        !        g(x,y)      = y*exp[-(16x**2-4y**2)]
        ! (x,y,z) --> (streamwise, spanwise, wall-normal) directions
        !
        ! see Henningson and Kim, JFM 1991
        !
        rfn(:) = 2.*rfn(:)/l(idir_n) - 1.
        rn(:) = 2.*rn(:)/l(idir_n) - 1. ! normal coordinate between -1 and +1
        rs(:)  = (rs(:) -.5*l(idir_s))*2./l(idir_n)
        rt(:)  = (rt(:) -.5*l(idir_t))*2./l(idir_n)
        rft(:) = (rft(:)-.5*l(idir_t))*2./l(idir_n)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              ii = [i,j,k]
              xxc = rs(ii(idir_s))
              yyc = rt(ii(idir_t)); yyf = rft(ii(idir_t))
              zzc = rn(ii(idir_n)); zzf = rfn(ii(idir_n))
              ut(i,j,k) = -1.*gxy(yyf,xxc)*dfz(zzc)*ubulk*1.5
              un(i,j,k) =  1.*fz(zzf)*dgxy(yyc,xxc)*ubulk*1.5
              p(i,j,k) = 0.
            end do
          end do
        end do
      else
        !
        ! alternatively, using a Taylor-Green vortex
        ! for the cross-stream velocity components
        !
        rfn(:) = rfn(:)/l(idir_n)*2.*pi
        rn(:)  = rn(:) /l(idir_n)*2.*pi
        rs(:)  = rs(:) /l(idir_s)*2.*pi
        rt(:)  = rt(:) /l(idir_t)*2.*pi
        rft(:) = rft(:)/l(idir_t)*2.*pi
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              ii = [i,j,k]
              xxc = rs(ii(idir_s))
              yyc = rt(ii(idir_t)); yyf = rft(ii(idir_t))
              zzc = rn(ii(idir_n)); zzf = rfn(ii(idir_n))
              ut(i,j,k) =  sin(xxc)*cos(yyf)*cos(zzc)*ubulk
              un(i,j,k) = -cos(xxc)*sin(yyc)*cos(zzf)*ubulk
              p(i,j,k) = 0.!(cos(2.*xxc)+cos(2.*yyc))*(cos(2.*zzc)+2.)/16.
            end do
          end do
        end do
      end if
    end if
  end subroutine initflow
  !
  subroutine initscal(iniscal,bcscal,ng,lo,l,xc,xf,yc,yf,zc,zf,dxc,dxf,dyc,dyf,dzc,dzf, &
                      salpha,is_sforced,scalf,s)
    !
    ! computes initial conditions for the scalar field
    !
    implicit none
    character(len=3), intent(in)                 :: iniscal
    real(rp), intent(in   ), dimension(0:1,3)    :: bcscal
    integer , intent(in   ), dimension(3)        :: ng,lo
    real(rp), intent(in   ), dimension(3)        :: l
    real(rp), intent(in   ), dimension(0:)       :: xc,xf,yc,yf,zc,zf,dxc,dxf,dyc,dyf,dzc,dzf
    real(rp), intent(in   )                      :: salpha
    logical , intent(in   )                      :: is_sforced
    real(rp), intent(in   )                      :: scalf
    real(rp), intent(inout), dimension(0:,0:,0:) :: s
    real(rp), allocatable, dimension(:) :: s1d
    integer :: i,j,k
    logical :: is_noise,is_mean
    real(rp) :: sref
    integer, dimension(3) :: n
    real(rp) :: xx
    !
    n(:) = shape(s) - 2*1
    allocate(s1d(n(3)))
    is_noise = .false.
    is_mean  = .false.
    !sref = 0.
    sref = 0.5*(bcscal(0,3)+bcscal(1,3)) ! bottom and top bcs
    if(is_sforced) sref = scalf
    select case(trim(iniscal))
    case('zer')
      s1d(:) = 0._rp
    case('uni')
      s1d(:) = sref
    case('cou')
      call couette(   n(3),zc/l(3),1._rp,s1d)
      s1d(:) = s1d(:) + 0.5 ! from 1 to 0
      s1d(:) = bcscal(0,3)*(s1d(:)) + bcscal(1,3)*(1.-s1d(:))
      sref = abs(bcscal(1,3) - bcscal(0,3))
    case('dhc')
      s1d(:) = 0._rp
    case('tbl')
      sref = 1.
      call temporal_bl(n(3),zc,1._rp,salpha,sref,s1d)
      is_noise = .true.
    case default
      if(myid == 0) print*, 'ERROR: invalid name for initial scalar field'
      if(myid == 0) print*, ''
      if(myid == 0) print*, '*** Simulation aborted due to errors in the case file ***'
      if(myid == 0) print*, '    check INFO_INPUT.md'
      call MPI_FINALIZE(ierr)
      error stop
    end select
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s1d(k)
        end do
      end do
    end do
    !
    if(trim(iniscal) == 'dhc') then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            xx = xc(i)/l(1)
            s(i,j,k) = ((1.-xx)*bcscal(0,1) + xx*bcscal(1,1))
          end do
        end do
      end do
    end if
    !
    if(is_noise) then
      call add_noise(ng,lo,123,.05_rp,s(1:n(1),1:n(2),1:n(3)))
    end if
    if(is_mean) then
      call set_mean(n,sref,dxf,dyf,dzf,l,s(1:n(1),1:n(2),1:n(3)))
    end if
  end subroutine initscal
  !
  subroutine add_noise(ng,lo,iseed,norm,p)
    implicit none
    integer , intent(in), dimension(3) :: ng,lo
    integer , intent(in) :: iseed
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(:,:,:) :: p
    integer(4), allocatable, dimension(:) :: seed
    real(rp) :: rn
    integer, dimension(3) :: n
    integer :: i,j,k,ii,jj,kk
    !
    n(:) = shape(p)
    allocate(seed(64))
    seed(:) = iseed
    call random_seed( put = seed )
    do k=1,ng(3)
      kk = k-(lo(3)-1)
      do j=1,ng(2)
        jj = j-(lo(2)-1)
        do i=1,ng(1)
          ii = i-(lo(1)-1)
          call random_number(rn)
          if(ii >= 1.and.ii <= n(1) .and. &
             jj >= 1.and.jj <= n(2) .and. &
             kk >= 1.and.kk <= n(3) ) then
             p(ii,jj,kk) = p(ii,jj,kk) + 2.*(rn-.5)*norm
          end if
        end do
      end do
    end do
  end subroutine add_noise
  !
  subroutine set_mean(n,mean,dx,dy,dz,l,p)
  implicit none
  integer , intent(in), dimension(3) :: n
  real(rp), intent(in) :: mean
  real(rp), intent(in), dimension(0:) :: dx,dy,dz
  real(rp), intent(in), dimension(3) :: l
  real(rp), intent(inout), dimension(:,:,:) :: p
  real(rp) :: vol,meanold
  integer :: i,j,k
  vol = product(l(:))
  meanold = 0.
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(shared) REDUCTION(+:meanold)
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        meanold = meanold + p(i,j,k)*dx(i)*dy(j)*dz(k)
      end do
    end do
  end do
  meanold = meanold/vol
  call MPI_ALLREDUCE(MPI_IN_PLACE,meanold,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  if(abs(meanold) > epsilon(0._rp)) then
    !$OMP PARALLEL WORKSHARE
    p(:,:,:) = p(:,:,:)/meanold*mean
    !$OMP END PARALLEL WORKSHARE
  end if
  end subroutine set_mean
  !
  subroutine couette(n,zc,norm,p)
    !
    ! plane couette profile normalized by the wall velocity difference
    !
    implicit none
    integer , intent(in)   :: n
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in)   :: norm
    real(rp), intent(out), dimension(n) :: p
    integer :: k
    real(rp) :: z
    do k=1,n
      z    = zc(k)
      p(k) = .5*(1.-2.*z)*norm
    end do
  end subroutine couette
  !
  subroutine poiseuille(n,zc,norm,p)
    implicit none
    integer , intent(in)   :: n
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in)   :: norm
    real(rp), intent(out), dimension(n) :: p
    integer :: k
    real(rp) :: z
    !
    ! plane poiseuille profile normalized by the bulk velocity
    !
    do k=1,n
      z    = zc(k)
      p(k) = 6.*z*(1.-z)*norm
    end do
  end subroutine poiseuille
  !
  subroutine temporal_bl(n,zc,d,nu,norm,p)
    implicit none
    integer , intent(in )   :: n
    real(rp), intent(in ), dimension(0:) :: zc
    real(rp), intent(in )   :: d,nu,norm
    real(rp), intent(out), dimension(n) :: p
    integer  :: k
    real(rp) :: theta
    !
    ! temporal boundary layer profile
    ! with thickness d, viscosity nu, and wall velocity norm (at z=0)
    !
    theta = 54.*nu/norm
    do k=1,n
      p(k)=(0.5+(0.5)*tanh((d/(2.*theta))*(1.-zc(k)/d)))*norm
    end do
  end subroutine temporal_bl
  !
  subroutine log_profile(n,zc,reb,p)
    implicit none
    integer , intent(in)   :: n
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in)   :: reb
    real(rp), intent(out), dimension(n) :: p
    integer :: k
    real(rp) :: z,retau ! z/lz and bulk Reynolds number
    retau = 0.09*reb**(0.88) ! from Pope's book
    do k=1,n
      z = zc(k)*2.*retau
      if(z >= retau) z = 2.*retau-z
      p(k) = 2.5*log(z) + 5.5
      if(z <= 11.6 ) p(k)=z
    end do
  end subroutine log_profile
  !
  ! functions to initialize the streamwise vortex pair
  ! (explained above)
  !
  function fz(zc)
  real(rp), intent(in) :: zc
  real(rp) :: fz
    fz = ((1.-zc**2)**2)
  end function
  !
  function dfz(zc)
  real(rp), intent(in) :: zc
  real(rp) :: dfz
    dfz = -4.*zc*(1.-zc**2)
  end function
  !
  function gxy(xc,yc)
  real(rp), intent(in) :: xc,yc
  real(rp) :: gxy
    gxy = yc*exp(-4.*(4.*xc**2+yc**2))
  end function
  !
  function dgxy(xc,yc)
  real(rp), intent(in) :: xc,yc
  real(rp) :: dgxy
    dgxy = exp(-4.*(4.*xc**2+yc**2))*(1.-8.*yc**2)
  end function
end module mod_initflow
