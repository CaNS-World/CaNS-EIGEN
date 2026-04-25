! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_linalg
  use mod_types, only: sp,dp
  implicit none
  private
  public :: gemm,stedc,syevd
  !
  ! explicit interfaces for external BLAS/LAPACK routines
  !
  interface
    subroutine sstedc(compz,n,d,e,z,ldz,work,lwork,iwork,liwork,info)
      import :: sp
      implicit none
      character, intent(in) :: compz
      integer, intent(in) :: n,ldz,lwork,liwork
      real(sp), intent(inout) :: d(*),e(*),z(ldz,*),work(*)
      integer, intent(inout) :: iwork(*)
      integer, intent(out) :: info
    end subroutine sstedc
    subroutine dstedc(compz,n,d,e,z,ldz,work,lwork,iwork,liwork,info)
      import :: dp
      implicit none
      character, intent(in) :: compz
      integer, intent(in) :: n,ldz,lwork,liwork
      real(dp), intent(inout) :: d(*),e(*),z(ldz,*),work(*)
      integer, intent(inout) :: iwork(*)
      integer, intent(out) :: info
    end subroutine dstedc
  end interface
  interface stedc
    procedure :: sstedc,dstedc
  end interface
  !
  interface
    subroutine ssyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
      import :: sp
      implicit none
      character, intent(in) :: jobz,uplo
      integer, intent(in) :: n,lda,lwork,liwork
      real(sp), intent(inout) :: a(lda,*)
      real(sp), intent(out) :: w(*)
      real(sp), intent(inout) :: work(*)
      integer, intent(inout) :: iwork(*)
      integer, intent(out) :: info
    end subroutine ssyevd
    subroutine dsyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
      import :: dp
      implicit none
      character, intent(in) :: jobz,uplo
      integer, intent(in) :: n,lda,lwork,liwork
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(out) :: w(*)
      real(dp), intent(inout) :: work(*)
      integer, intent(inout) :: iwork(*)
      integer, intent(out) :: info
    end subroutine dsyevd
  end interface
  interface syevd
    procedure :: ssyevd,dsyevd
  end interface
  !
  interface
    subroutine sgemm(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
      import :: sp
      implicit none
      character, intent(in) :: transa,transb
      integer, intent(in) :: m,n,k,lda,ldb,ldc
      real(sp), intent(in) :: alpha,beta
      real(sp), intent(in) :: A(lda,*),B(ldb,*)
      real(sp), intent(inout) :: C(ldc,*)
    end subroutine sgemm
    subroutine dgemm(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
      import :: dp
      implicit none
      character, intent(in) :: transa,transb
      integer, intent(in) :: m,n,k,lda,ldb,ldc
      real(dp), intent(in) :: alpha,beta
      real(dp), intent(in) :: A(lda,*),B(ldb,*)
      real(dp), intent(inout) :: C(ldc,*)
    end subroutine dgemm
  end interface
  interface gemm
    procedure :: sgemm,dgemm
  end interface
end module mod_linalg
