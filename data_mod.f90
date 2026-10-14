module data_mod
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  ! Double precision definition
  integer, parameter :: dp = real64
  real(dp), parameter :: pi = 3.14159265358979323846_dp

  ! Problem dimensions
  integer :: n = 0        ! Number of boundary nodes
  integer :: ne = 0       ! Number of boundary elements
  integer :: l = 0        ! Number of internal points
  integer :: ntrack = 1   ! Tracked node for disp_trac.dat output

  ! File unit numbers and names
  integer, parameter :: inp = 15
  integer, parameter :: ipr = 16
  integer, parameter :: iout_track = 25

  character(len=64) :: filein = 'cube8.dat'
  character(len=64) :: fileout = 'cube8_out.dat'
  character(len=64) :: filetrack = 'disp_trac.dat'
  character(len=80) :: title = ''

  ! Material properties (Elastostatic 3D)
  real(dp) :: young = 0.0_dp  ! Young's modulus (E)
  real(dp) :: rge = 0.0_dp    ! Shear modulus (mu or G)
  real(dp) :: xnu = 0.0_dp    ! Poisson's ratio (nu)

  ! Kelvin fundamental solution constants (3D elastostatics)
  real(dp) :: c3 = 0.0_dp     ! 1 / (16 * pi * (1 - nu) * mu)
  real(dp) :: c4 = 0.0_dp     ! 3 - 4 * nu
  real(dp) :: c5 = 0.0_dp     ! -1 / (8 * pi * (1 - nu))
  real(dp) :: c6 = 0.0_dp     ! 1 - 2 * nu
  real(dp) :: delta(3, 3) = 0.0_dp

  ! Dynamic Arrays (Real double precision for elastostatics)
  real(dp), allocatable :: x(:), y(:), z(:)
  real(dp), allocatable :: cx(:), cy(:), cz(:)

  integer, allocatable :: kode(:)
  real(dp), allocatable :: bc(:,:)

  integer, allocatable :: nconec(:)
  integer, allocatable :: kconec(:,:)
  integer, allocatable :: nesim(:)
  integer, allocatable :: ke(:,:)

  real(dp), allocatable :: A(:,:)
  real(dp), allocatable :: FI(:)
  real(dp), allocatable :: DFI(:,:)
  real(dp), allocatable :: dsol(:)

contains

  subroutine allocate_bem_data(n_nodes, n_elems, n_int_pts)
    implicit none
    integer, intent(in) :: n_nodes, n_elems, n_int_pts

    n = n_nodes
    ne = n_elems
    l = n_int_pts

    if (allocated(x)) deallocate(x, y, z)
    allocate(x(n), y(n), z(n))

    if (l > 0) then
      if (allocated(cx)) deallocate(cx, cy, cz, dsol)
      allocate(cx(l), cy(l), cz(l), dsol(3*l))
    end if

    if (allocated(kode)) deallocate(kode, bc)
    allocate(kode(3*ne), bc(27, ne))

    if (allocated(nconec)) deallocate(nconec, kconec, nesim, ke)
    allocate(nconec(ne), kconec(9, ne), nesim(n), ke(8, n))

    if (allocated(A)) deallocate(A, FI, DFI)
    allocate(A(3*n, 3*n), FI(3*n), DFI(27, ne))

    ! Initialize Kronecker delta
    delta = 0.0_dp
    delta(1, 1) = 1.0_dp
    delta(2, 2) = 1.0_dp
    delta(3, 3) = 1.0_dp
  end subroutine allocate_bem_data

  subroutine deallocate_bem_data()
    implicit none
    if (allocated(x)) deallocate(x, y, z)
    if (allocated(cx)) deallocate(cx, cy, cz, dsol)
    if (allocated(kode)) deallocate(kode, bc)
    if (allocated(nconec)) deallocate(nconec, kconec, nesim, ke)
    if (allocated(A)) deallocate(A, FI, DFI)
  end subroutine deallocate_bem_data

  subroutine update_material_constants()
    implicit none
    if (young > 0.0_dp .and. rge == 0.0_dp) then
      rge = young / (2.0_dp * (1.0_dp + xnu))
    else if (rge > 0.0_dp .and. young == 0.0_dp) then
      young = 2.0_dp * rge * (1.0_dp + xnu)
    end if
    c3 = 1.0_dp / (16.0_dp * pi * (1.0_dp - xnu) * rge)
    c4 = 3.0_dp - 4.0_dp * xnu
    c5 = -1.0_dp / (8.0_dp * pi * (1.0_dp - xnu))
    c6 = 1.0_dp - 2.0_dp * xnu
  end subroutine update_material_constants

end module data_mod
