module acc_linalg_m
  use ISO_C_BINDING
  use algos_common_m
  use blas_m
  use message_m, only: die
  use nrtype_m
  implicit none
  private
  public :: acc_dgemm
  public :: acc_xgemm
  public :: acc_zgemm
  ! WPH: Implementation of dgemm as a generic subroutine
  interface acc_xgemm
    module procedure acc_zgemm
    module procedure acc_dgemm
  end interface
contains
  subroutine acc_zgemm(transa, transb, &
                       m, n, k, &
                       alpha, a, lda, b, ldb, &
                       beta, c, ldc, &
                       algo)
    implicit none
    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double complex :: alpha, beta
    double complex :: a(lda, *), b(ldb, *), c(ldc, *)
    integer(kind(CPU_ALGO)) :: algo
    select case (algo)
    case (OPENACC_ALGO)
      call die("OpenACC version of acc_zgemm requested, but OpenACC not "&
            &"compiled into this executable", only_root_writes = .true.)
    case (CPU_ALGO)
      call zgemm(transa, transb, &
                 m, n, k, &
                 alpha, a, lda, b, ldb, &
                 beta, c, ldc)
    case default
      call die("Invald algorithm for acc_zgemm", only_root_writes = .true.)
    end select
  end subroutine acc_zgemm
  subroutine acc_dgemm(transa, transb, &
                       m, n, k, &
                       alpha, a, lda, b, ldb, &
                       beta, c, ldc, &
                       algo)
    implicit none
    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double precision :: alpha, beta
    double precision :: a(lda, *), b(ldb, *), c(ldc, *)
    integer(kind(CPU_ALGO)) :: algo
    select case (algo)
    case (OPENACC_ALGO)
      call die("OpenACC version of acc_dgemm requested, but OpenACC not "&
            &"compiled into this executable", only_root_writes = .true.)
    case (CPU_ALGO)
      call dgemm(transa, transb, &
                 m, n, k, &
                 alpha, a, lda, b, ldb, &
                 beta, c, ldc)
    case default
      call die("Invald algorithm for acc_dgemm", only_root_writes = .true.)
    end select
  end subroutine acc_dgemm
end module acc_linalg_m
