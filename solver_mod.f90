module solver_mod
  use data_mod, only: dp
  implicit none

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      use, intrinsic :: iso_fortran_env, only: real64
      integer, intent(in) :: n, nrhs, lda, ldb
      real(real64), intent(inout) :: a(lda, *), b(ldb, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine dgesv
  end interface

contains

  subroutine solve_linear_system(A_mat, B_vec, n_size)
    implicit none
    integer, intent(in) :: n_size
    real(dp), intent(inout) :: A_mat(n_size, n_size)
    real(dp), intent(inout) :: B_vec(n_size)

    integer, allocatable :: ipiv(:)
    integer :: info

    allocate(ipiv(n_size))

    call dgesv(n_size, 1, A_mat, n_size, ipiv, B_vec, n_size, info)

    if (info /= 0) then
      write(*, '(A, I6)') ' *** ERROR in LAPACK DGESV: info = ', info
      stop 'Solver failure: Singular matrix or bad arguments in LAPACK DGESV.'
    end if

    deallocate(ipiv)
  end subroutine solve_linear_system

end module solver_mod
