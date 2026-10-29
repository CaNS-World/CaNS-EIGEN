! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_linalg
  use mod_types, only: sp,dp,rp
#if defined(_OPENACC)
#if !defined(_USE_HIP)
#if !defined(_SINGLE_PRECISION)
  use cublas         , only: gemm_gpu => cublasDgemm_v2, &
#else
  use cublas         , only: gemm_gpu => cublasSgemm_v2, &
#endif
                             gemm_op_n => CUBLAS_OP_N
#else
#if !defined(_SINGLE_PRECISION)
  use hipfort_hipblas, only: hipblas_gemm => hipblasDgemm, &
#else
  use hipfort_hipblas, only: hipblas_gemm => hipblasSgemm, &
#endif
                             gemm_op_n => HIPBLAS_OP_N
#endif
#endif
  implicit none
  private
  public :: gemm,stedc,syevd
#if defined(_OPENACC)
  public :: gemm_gpu,gemm_op_n
#endif
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
#if defined(_OPENACC) && defined(_USE_HIP)
  contains
  function gemm_gpu(handle,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) result(istat)
    !
    ! wrapper around the hipfort C binding of hipblas{D,S}gemm, which takes device addresses
    ! as `type(c_ptr)`; to be called within an OpenACC `host_data use_device` region
    !
    use, intrinsic :: iso_c_binding, only: c_ptr,c_loc
    implicit none
    type(c_ptr), intent(in) :: handle
    integer(kind(gemm_op_n)), intent(in) :: transa,transb
    integer , intent(in) :: m,n,k,lda,ldb,ldc
    real(rp), intent(in) :: alpha,beta
    real(rp), intent(in   ), target, dimension(*) :: a,b
    real(rp), intent(inout), target, dimension(*) :: c
    integer :: istat
    istat = hipblas_gemm(handle,transa,transb,m,n,k,alpha,c_loc(a),lda,c_loc(b),ldb,beta,c_loc(c),ldc)
  end function gemm_gpu
#endif
end module mod_linalg
