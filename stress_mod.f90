module stress_mod
  use data_mod
  use assembly_mod, only: all_eta1, all_eta2, all_eta3, all_xja, &
                          all_xg1, all_xg2, all_xg3, all_ff, &
                          ome_quad, omei_tri, omej_tri
  implicit none

  ! Dynamic Arrays for boundary stresses
  ! Cauchy stress ordering: 1=xx, 2=yy, 3=zz, 4=xy, 5=yz, 6=zx
  real(dp), allocatable, public :: sigma_elem(:, :, :)   ! (6, 9, ne): element-nodal stresses
  real(dp), allocatable, public :: sigma_node(:, :)      ! (6, n): nodal averaged stresses
  real(dp), allocatable, public :: sigma_principal(:, :) ! (3, n): principal stresses (s1 >= s2 >= s3)
  real(dp), allocatable, public :: sigma_vm(:)           ! (n): Von Mises equivalent stress

  ! Dynamic Arrays for internal points
  real(dp), allocatable, public :: u_int(:, :)          ! (3, l): internal displacements
  real(dp), allocatable, public :: sigma_int(:, :)      ! (6, l): internal Cauchy stresses
  real(dp), allocatable, public :: sigma_int_prin(:, :) ! (3, l): internal principal stresses
  real(dp), allocatable, public :: sigma_int_vm(:)      ! (l): internal Von Mises stress

  character(len=64), public :: filestress = 'stresses.dat'
  character(len=64), public :: fileintpts = 'internal_points.dat'

contains

  !=============================================================================
  ! shape_deriv_lagrange
  ! Computes shape functions and intrinsic derivatives for:
  ! - 9-node Lagrangian quadrilateral elements (ncone = 9)
  ! - 6-node quadratic triangular elements (ncone = 6)
  !=============================================================================
  subroutine shape_deriv_lagrange(ncone, g1, g2, ff, fd1, fd2)
    implicit none
    integer, intent(in) :: ncone
    real(dp), intent(in) :: g1, g2
    real(dp), intent(out) :: ff(9), fd1(9), fd2(9)

    real(dp) :: g3

    ff = 0.0_dp
    fd1 = 0.0_dp
    fd2 = 0.0_dp

    if (ncone == 9) then
      ! 9-node Lagrangian quadrilateral element
      ! Intrinsic coordinates (g1, g2) in [-1, 1] x [-1, 1]
      ! Node ordering:
      ! 1: (-1, -1), 2: ( 0, -1), 3: ( 1, -1)
      ! 4: ( 1,  0), 5: ( 1,  1), 6: ( 0,  1)
      ! 7: (-1,  1), 8: (-1,  0), 9: ( 0,  0)
      ff(1) = 0.25_dp * g1 * g2 * (g1 - 1.0_dp) * (g2 - 1.0_dp)
      ff(2) = 0.5_dp  * (1.0_dp - g1**2) * g2 * (g2 - 1.0_dp)
      ff(3) = 0.25_dp * g1 * g2 * (g1 + 1.0_dp) * (g2 - 1.0_dp)
      ff(4) = 0.5_dp  * g1 * (g1 + 1.0_dp) * (1.0_dp - g2**2)
      ff(5) = 0.25_dp * g1 * g2 * (g1 + 1.0_dp) * (g2 + 1.0_dp)
      ff(6) = 0.5_dp  * (1.0_dp - g1**2) * g2 * (g2 + 1.0_dp)
      ff(7) = 0.25_dp * g1 * g2 * (g1 - 1.0_dp) * (g2 + 1.0_dp)
      ff(8) = 0.5_dp  * g1 * (g1 - 1.0_dp) * (1.0_dp - g2**2)
      ff(9) = (1.0_dp - g1**2) * (1.0_dp - g2**2)

      ! Intrinsic derivatives with respect to g1 (xi_1)
      fd1(1) = 0.25_dp * g2 * (g2 - 1.0_dp) * (2.0_dp * g1 - 1.0_dp)
      fd1(2) = -g1 * g2 * (g2 - 1.0_dp)
      fd1(3) = 0.25_dp * g2 * (g2 - 1.0_dp) * (2.0_dp * g1 + 1.0_dp)
      fd1(4) = 0.5_dp  * (1.0_dp - g2**2) * (2.0_dp * g1 + 1.0_dp)
      fd1(5) = 0.25_dp * g2 * (g2 + 1.0_dp) * (2.0_dp * g1 + 1.0_dp)
      fd1(6) = -g1 * g2 * (g2 + 1.0_dp)
      fd1(7) = 0.25_dp * g2 * (g2 + 1.0_dp) * (2.0_dp * g1 - 1.0_dp)
      fd1(8) = 0.5_dp  * (1.0_dp - g2**2) * (2.0_dp * g1 - 1.0_dp)
      fd1(9) = -2.0_dp * g1 * (1.0_dp - g2**2)

      ! Intrinsic derivatives with respect to g2 (xi_2)
      fd2(1) = 0.25_dp * g1 * (g1 - 1.0_dp) * (2.0_dp * g2 - 1.0_dp)
      fd2(2) = 0.5_dp  * (1.0_dp - g1**2) * (2.0_dp * g2 - 1.0_dp)
      fd2(3) = 0.25_dp * g1 * (g1 + 1.0_dp) * (2.0_dp * g2 - 1.0_dp)
      fd2(4) = -g1 * g2 * (g1 + 1.0_dp)
      fd2(5) = 0.25_dp * g1 * (g1 + 1.0_dp) * (2.0_dp * g2 + 1.0_dp)
      fd2(6) = 0.5_dp  * (1.0_dp - g1**2) * (2.0_dp * g2 + 1.0_dp)
      fd2(7) = 0.25_dp * g1 * (g1 - 1.0_dp) * (2.0_dp * g2 + 1.0_dp)
      fd2(8) = -g1 * g2 * (g1 - 1.0_dp)
      fd2(9) = -2.0_dp * g2 * (1.0_dp - g1**2)
    else
      ! 6-node quadratic triangular element
      ! Area coordinates (g1, g2, g3) with g3 = 1 - g1 - g2
      ! Node ordering:
      ! 1: (1, 0, 0), 2: (0, 1, 0), 3: (0, 0, 1)
      ! 4: (1/2, 1/2, 0), 5: (0, 1/2, 1/2), 6: (1/2, 0, 1/2)
      g3 = 1.0_dp - g1 - g2
      ff(1) = g1 * (2.0_dp * g1 - 1.0_dp)
      ff(2) = g2 * (2.0_dp * g2 - 1.0_dp)
      ff(3) = g3 * (2.0_dp * g3 - 1.0_dp)
      ff(4) = 4.0_dp * g1 * g2
      ff(5) = 4.0_dp * g2 * g3
      ff(6) = 4.0_dp * g3 * g1

      ! Derivatives with respect to g1
      fd1(1) = 4.0_dp * g1 - 1.0_dp
      fd1(2) = 0.0_dp
      fd1(3) = 1.0_dp - 4.0_dp * g3
      fd1(4) = 4.0_dp * g2
      fd1(5) = -4.0_dp * g2
      fd1(6) = 4.0_dp * (g3 - g1)

      ! Derivatives with respect to g2
      fd2(1) = 0.0_dp
      fd2(2) = 4.0_dp * g2 - 1.0_dp
      fd2(3) = 1.0_dp - 4.0_dp * g3
      fd2(4) = 4.0_dp * g1
      fd2(5) = 4.0_dp * (g3 - g2)
      fd2(6) = -4.0_dp * g1
    end if
  end subroutine shape_deriv_lagrange

  !=============================================================================
  ! axes_cos_lagrange
  ! Computes the local orthonormal coordinate triad (e1, e2, e3) and the
  ! surface transformation matrix CL (Gao & Davies Eq. 4.71) at a given point.
  !=============================================================================
  subroutine axes_cos_lagrange(ncone, fd1, fd2, coords, L, CL, v1, v2)
    implicit none
    integer, intent(in) :: ncone
    real(dp), intent(in) :: fd1(9), fd2(9), coords(3, 9)
    real(dp), intent(out) :: L(3, 3), CL(2, 2), v1, v2

    real(dp) :: g1(3), g2(3), g3(3), jcb, cos_th, sin_th
    integer :: a, i

    g1 = 0.0_dp
    g2 = 0.0_dp
    do a = 1, ncone
      do i = 1, 3
        g1(i) = g1(i) + fd1(a) * coords(i, a)
        g2(i) = g2(i) + fd2(a) * coords(i, a)
      end do
    end do

    v1 = sqrt(dot_product(g1, g1))
    v2 = sqrt(dot_product(g2, g2))

    ! Cross product g3 = g1 x g2
    g3(1) = g1(2) * g2(3) - g1(3) * g2(2)
    g3(2) = g1(3) * g2(1) - g1(1) * g2(3)
    g3(3) = g1(1) * g2(2) - g1(2) * g2(1)

    jcb = sqrt(dot_product(g3, g3))

    ! Row 3 of L: unit normal e3 (normal vector)
    L(3, 1:3) = g3(1:3) / jcb

    ! Row 1 of L: unit tangent e1
    L(1, 1:3) = g1(1:3) / v1

    ! Row 2 of L: unit tangent e2 = e3 x e1 (completing orthonormal triad)
    L(2, 1) = L(3, 2) * L(1, 3) - L(3, 3) * L(1, 2)
    L(2, 2) = L(3, 3) * L(1, 1) - L(3, 1) * L(1, 3)
    L(2, 3) = L(3, 1) * L(1, 2) - L(3, 2) * L(1, 1)

    ! Transformation matrix CL relating intrinsic to local Cartesian coordinates
    cos_th = dot_product(g1, g2) / (v1 * v2)
    if (cos_th > 1.0_dp) cos_th = 1.0_dp
    if (cos_th < -1.0_dp) cos_th = -1.0_dp

    sin_th = sqrt(max(0.0_dp, 1.0_dp - cos_th**2))
    if (sin_th < 1.0e-12_dp) sin_th = 1.0e-12_dp

    CL(1, 1) = 1.0_dp / v1
    CL(1, 2) = -cos_th / (v1 * sin_th)
    CL(2, 1) = 0.0_dp
    CL(2, 2) = 1.0_dp / (v2 * sin_th)
  end subroutine axes_cos_lagrange

  !=============================================================================
  ! boundary_stresses_elem
  ! Evaluates full Cauchy stress tensor at all nodes of boundary element ie.
  ! Traction-recovery method based on local elasticity Hooke's Law:
  ! - Out-of-plane stresses obtained directly from surface tractions
  ! - In-plane stresses obtained from surface displacement gradients
  !=============================================================================
  subroutine boundary_stresses_elem(ie, sig_elem)
    implicit none
    integer, intent(in) :: ie
    real(dp), intent(out) :: sig_elem(6, 9)

    integer :: ncone, k, a, ik
    real(dp) :: coords(3, 9), u_elem(3, 9), t_elem(3, 9)
    real(dp) :: gi1(9), gi2(9), g1, g2
    real(dp) :: ff(9), fd1(9), fd2(9)
    real(dp) :: L(3, 3), CL(2, 2), v1, v2
    real(dp) :: dN_dx1(9), dN_dx2(9)
    real(dp) :: u_loc(3, 9), t_loc(3)
    real(dp) :: u1_1, u2_2, u1_2, u2_1, eps11, eps22, eps12
    real(dp) :: S_loc(3, 3), S_glob(3, 3)
    real(dp) :: c_eps, c_nu, c_sig3

    ncone = nconec(ie)
    sig_elem = 0.0_dp

    ! Retrieve nodal coordinates, displacements and tractions for this element
    do a = 1, ncone
      ik = kconec(a, ie)
      coords(1, a) = x(ik)
      coords(2, a) = y(ik)
      coords(3, a) = z(ik)

      u_elem(1, a) = FI(3*ik - 2)
      u_elem(2, a) = FI(3*ik - 1)
      u_elem(3, a) = FI(3*ik)

      t_elem(1, a) = DFI(3*a - 2, ie)
      t_elem(2, a) = DFI(3*a - 1, ie)
      t_elem(3, a) = DFI(3*a,     ie)
    end do

    ! Set node parametric coordinates
    if (ncone == 9) then
      gi1 = [-1.0_dp,  0.0_dp,  1.0_dp,  1.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, 0.0_dp]
      gi2 = [-1.0_dp, -1.0_dp, -1.0_dp,  0.0_dp, 1.0_dp, 1.0_dp,  1.0_dp,  0.0_dp, 0.0_dp]
    else
      gi1(1:6) = [1.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp]
      gi2(1:6) = [0.0_dp, 1.0_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp]
    end if

    ! Elastic constants for Hooke's Law in local tangent plane
    c_eps  = 2.0_dp * rge / (1.0_dp - xnu)
    c_nu   = 2.0_dp * rge * xnu / (1.0_dp - xnu)
    c_sig3 = xnu / (1.0_dp - xnu)

    ! Compute stresses at each node k of element ie
    do k = 1, ncone
      g1 = gi1(k)
      g2 = gi2(k)

      ! 1. Evaluate shape functions and intrinsic derivatives
      call shape_deriv_lagrange(ncone, g1, g2, ff, fd1, fd2)

      ! 2. Form orthonormal triad L and surface transformation CL
      call axes_cos_lagrange(ncone, fd1, fd2, coords, L, CL, v1, v2)

      ! 3. Local Cartesian derivatives of shape functions on tangent plane
      do a = 1, ncone
        dN_dx1(a) = CL(1, 1) * fd1(a) + CL(2, 1) * fd2(a)
        dN_dx2(a) = CL(1, 2) * fd1(a) + CL(2, 2) * fd2(a)
      end do

      ! 4. Local displacement of all element nodes in the local basis of node k
      do a = 1, ncone
        u_loc(1, a) = L(1, 1)*u_elem(1, a) + L(1, 2)*u_elem(2, a) + L(1, 3)*u_elem(3, a)
        u_loc(2, a) = L(2, 1)*u_elem(1, a) + L(2, 2)*u_elem(2, a) + L(2, 3)*u_elem(3, a)
        u_loc(3, a) = L(3, 1)*u_elem(1, a) + L(3, 2)*u_elem(2, a) + L(3, 3)*u_elem(3, a)
      end do

      ! 5. Local displacement gradients and in-plane strains at node k
      u1_1 = 0.0_dp; u2_2 = 0.0_dp; u1_2 = 0.0_dp; u2_1 = 0.0_dp
      do a = 1, ncone
        u1_1 = u1_1 + dN_dx1(a) * u_loc(1, a)
        u2_2 = u2_2 + dN_dx2(a) * u_loc(2, a)
        u1_2 = u1_2 + dN_dx2(a) * u_loc(1, a)
        u2_1 = u2_1 + dN_dx1(a) * u_loc(2, a)
      end do

      eps11 = u1_1
      eps22 = u2_2
      eps12 = 0.5_dp * (u1_2 + u2_1)

      ! 6. Local traction at node k: t_loc = L * t_elem
      t_loc(1) = L(1, 1)*t_elem(1, k) + L(1, 2)*t_elem(2, k) + L(1, 3)*t_elem(3, k)
      t_loc(2) = L(2, 1)*t_elem(1, k) + L(2, 2)*t_elem(2, k) + L(2, 3)*t_elem(3, k)
      t_loc(3) = L(3, 1)*t_elem(1, k) + L(3, 2)*t_elem(2, k) + L(3, 3)*t_elem(3, k)

      ! 7. Out-of-plane stresses from surface traction (Cauchy's formula)
      S_loc(1, 3) = t_loc(1); S_loc(3, 1) = S_loc(1, 3)
      S_loc(2, 3) = t_loc(2); S_loc(3, 2) = S_loc(2, 3)
      S_loc(3, 3) = t_loc(3)

      ! 8. In-plane stresses from isotropic Hooke's Law
      S_loc(1, 1) = c_eps * eps11 + c_nu * eps22 + c_sig3 * S_loc(3, 3)
      S_loc(2, 2) = c_eps * eps22 + c_nu * eps11 + c_sig3 * S_loc(3, 3)
      S_loc(1, 2) = 2.0_dp * rge * eps12
      S_loc(2, 1) = S_loc(1, 2)

      ! 9. Transform local stress tensor back to global coordinates: S_glob = L^T * S_loc * L
      S_glob = matmul(transpose(L), matmul(S_loc, L))

      ! 10. Store Cauchy stress components (xx, yy, zz, xy, yz, zx)
      sig_elem(1, k) = S_glob(1, 1)
      sig_elem(2, k) = S_glob(2, 2)
      sig_elem(3, k) = S_glob(3, 3)
      sig_elem(4, k) = 0.5_dp * (S_glob(1, 2) + S_glob(2, 1))
      sig_elem(5, k) = 0.5_dp * (S_glob(2, 3) + S_glob(3, 2))
      sig_elem(6, k) = 0.5_dp * (S_glob(3, 1) + S_glob(1, 3))
    end do
  end subroutine boundary_stresses_elem

  !=============================================================================
  ! calc_principal_and_vm
  ! Computes Von Mises equivalent stress and the 3 principal stresses
  ! sorted descending (s1 >= s2 >= s3) using exact analytical Cardano solution.
  !=============================================================================
  subroutine calc_principal_and_vm(sig, sig_prin, sig_vm)
    implicit none
    real(dp), intent(in) :: sig(6)
    real(dp), intent(out) :: sig_prin(3), sig_vm

    real(dp) :: sxx, syy, szz, sxy, syz, szx
    real(dp) :: i1, q, p_xx, p_yy, p_zz, j2, j3, r, phi, cos_3phi, tmp
    real(dp), parameter :: tol = 1.0e-20_dp

    sxx = sig(1); syy = sig(2); szz = sig(3)
    sxy = sig(4); syz = sig(5); szx = sig(6)

    ! Von Mises equivalent stress
    sig_vm = sqrt(0.5_dp * ((sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2 + &
                            6.0_dp * (sxy**2 + syz**2 + szx**2)))

    ! Mean stress / hydrostatic pressure
    i1 = sxx + syy + szz
    q = i1 / 3.0_dp

    ! Deviatoric stress tensor components
    p_xx = sxx - q
    p_yy = syy - q
    p_zz = szz - q

    ! Second invariant of deviatoric stress J2
    j2 = 0.5_dp * (p_xx**2 + p_yy**2 + p_zz**2) + sxy**2 + syz**2 + szx**2

    if (j2 <= tol * (i1**2 + 1.0_dp)) then
      sig_prin(1) = q
      sig_prin(2) = q
      sig_prin(3) = q
    else
      ! Third invariant of deviatoric stress J3 = det(s)
      j3 = p_xx * (p_yy * p_zz - syz**2) - &
           sxy  * (sxy  * p_zz - syz * szx) + &
           szx  * (sxy  * syz  - p_yy * szx)

      r = sqrt(4.0_dp * j2 / 3.0_dp)
      cos_3phi = 4.0_dp * j3 / (r**3)
      if (cos_3phi > 1.0_dp)  cos_3phi = 1.0_dp
      if (cos_3phi < -1.0_dp) cos_3phi = -1.0_dp

      phi = acos(cos_3phi) / 3.0_dp

      sig_prin(1) = q + r * cos(phi)
      sig_prin(3) = q + r * cos(phi + 2.0_dp * pi / 3.0_dp)
      sig_prin(2) = 3.0_dp * q - sig_prin(1) - sig_prin(3)

      ! Sort descending: sig_prin(1) >= sig_prin(2) >= sig_prin(3)
      if (sig_prin(1) < sig_prin(2)) then
        tmp = sig_prin(1); sig_prin(1) = sig_prin(2); sig_prin(2) = tmp
      end if
      if (sig_prin(2) < sig_prin(3)) then
        tmp = sig_prin(2); sig_prin(2) = sig_prin(3); sig_prin(3) = tmp
      end if
      if (sig_prin(1) < sig_prin(2)) then
        tmp = sig_prin(1); sig_prin(1) = sig_prin(2); sig_prin(2) = tmp
      end if
    end if
  end subroutine calc_principal_and_vm

  !=============================================================================
  ! compute_boundary_stresses
  ! Master routine to compute element-nodal stresses, nodal averaged stresses,
  ! principal stresses and Von Mises stress across the entire boundary.
  !=============================================================================
  subroutine compute_boundary_stresses()
    implicit none
    integer :: ie, k, in_idx, j, n_shared

    if (allocated(sigma_elem)) deallocate(sigma_elem)
    allocate(sigma_elem(6, 9, ne))
    sigma_elem = 0.0_dp

    if (allocated(sigma_node)) deallocate(sigma_node, sigma_principal, sigma_vm)
    allocate(sigma_node(6, n), sigma_principal(3, n), sigma_vm(n))
    sigma_node = 0.0_dp
    sigma_principal = 0.0_dp
    sigma_vm = 0.0_dp

    ! 1. Calculate element-nodal stresses for all elements (parallelized with OpenMP)
    !$omp parallel do private(ie) schedule(dynamic)
    do ie = 1, ne
      call boundary_stresses_elem(ie, sigma_elem(:, :, ie))
    end do
    !$omp end parallel do

    ! 2. Compute nodal averaged stresses
    do in_idx = 1, n
      n_shared = min(nesim(in_idx), 8)
      if (n_shared > 0) then
        do j = 1, n_shared
          ie = ke(j, in_idx)
          do k = 1, nconec(ie)
            if (kconec(k, ie) == in_idx) then
              sigma_node(:, in_idx) = sigma_node(:, in_idx) + sigma_elem(:, k, ie)
              exit
            end if
          end do
        end do
        sigma_node(:, in_idx) = sigma_node(:, in_idx) / real(n_shared, dp)
      end if

      ! 3. Calculate Von Mises and Principal stresses for this node
      call calc_principal_and_vm(sigma_node(:, in_idx), sigma_principal(:, in_idx), sigma_vm(in_idx))
    end do
  end subroutine compute_boundary_stresses

  !=============================================================================
  ! output_boundary_stresses
  ! Outputs boundary stresses to the main output file (ipr) and to stresses.dat
  !=============================================================================
  subroutine output_boundary_stresses()
    implicit none
    integer :: i, ie, k, ik, iunit_str
    integer, parameter :: iunit_file = 26

    ! Write Nodal Averaged Stresses to ipr (cube8_out.dat)
    write(ipr, '(A)') ' ' // repeat('*', 79)
    write(ipr, '(A)') ' BOUNDARY NODES STRESSES (NODAL AVERAGED)'
    write(ipr, '(A)') ' ' // repeat('*', 79)
    write(ipr, 100)
100 format(1X, 'NODE', 6X, 'SIGMA_XX', 8X, 'SIGMA_YY', 8X, 'SIGMA_ZZ', 8X, &
           'SIGMA_XY', 8X, 'SIGMA_YZ', 8X, 'SIGMA_ZX', 7X, 'VON_MISES', 9X, &
           'SIGMA_1', 9X, 'SIGMA_2', 9X, 'SIGMA_3'/)

    do i = 1, n
      write(ipr, 110) i, (sigma_node(k, i), k = 1, 6), sigma_vm(i), &
                         (sigma_principal(k, i), k = 1, 3)
    end do
110 format(1X, I5, 10(1X, E15.7))

    ! Write Element-Nodal Stresses to ipr
    write(ipr, '(//A)') ' ' // repeat('*', 79)
    write(ipr, '(A)') ' BOUNDARY ELEMENTS NODAL STRESSES'
    write(ipr, '(A)') ' ' // repeat('*', 79)
    write(ipr, 200)
200 format(1X, 'ELEMENT', 2X, 'NODE', 6X, 'SIGMA_XX', 8X, 'SIGMA_YY', 8X, &
           'SIGMA_ZZ', 8X, 'SIGMA_XY', 8X, 'SIGMA_YZ', 8X, 'SIGMA_ZX'/)

    do ie = 1, ne
      write(ipr, '(1X, I5)') ie
      do k = 1, nconec(ie)
        ik = kconec(k, ie)
        write(ipr, 210) ik, (sigma_elem(i, k, ie), i = 1, 6)
      end do
    end do
210 format(10X, I5, 6(1X, E15.7))

    ! Write dedicated stresses.dat file for post-processing and plotting
    open(unit=iunit_file, file=trim(filestress), status='replace', action='write')
    write(iunit_file, '(A)') '# BOUNDARY NODES STRESSES (QUADESTH 3D ELASTOSTATICS)'
    write(iunit_file, '(A)') '# Col 1: Node ID'
    write(iunit_file, '(A)') '# Col 2-4: X, Y, Z coordinates'
    write(iunit_file, '(A)') '# Col 5-10: Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_yz, Sigma_zx'
    write(iunit_file, '(A)') '# Col 11: Von Mises Equivalent Stress'
    write(iunit_file, '(A)') '# Col 12-14: Principal Stresses (Sigma_1 >= Sigma_2 >= Sigma_3)'
    write(iunit_file, '(A)') '# ' // repeat('-', 140)

    do i = 1, n
      write(iunit_file, 300) i, x(i), y(i), z(i), &
                            (sigma_node(k, i), k = 1, 6), sigma_vm(i), &
                            (sigma_principal(k, i), k = 1, 3)
    end do
300 format(I6, 3(1X, E15.7), 10(1X, E15.7))

    close(iunit_file)
    write(*, '(A)') ' Boundary stresses exported to: ' // trim(filestress)
  end subroutine output_boundary_stresses

  !=============================================================================
  ! compute_internal_points
  ! Evaluates displacements and Cauchy stresses at all internal points
  ! using Somigliana displacement and stress identities with OpenMP parallelization
  !=============================================================================
  subroutine compute_internal_points()
    implicit none
    integer :: ip, ie, k, m, ik, ncone, npg, i, j, igp, tot_gp
    real(dp) :: wg, ff_val
    real(dp) :: ug(3), tg(3)

    ! Flattened Gauss point arrays
    real(dp), allocatable :: gp_x(:), gp_y(:), gp_z(:)
    real(dp), allocatable :: gp_nx(:), gp_ny(:), gp_nz(:)
    real(dp), allocatable :: gp_u(:,:), gp_t(:,:)

    ! Loop variables for internal points
    real(dp) :: xi_x, xi_y, xi_z
    real(dp) :: rx, ry, rz, r2, r, inv_r, inv_r2, inv_r3
    real(dp) :: rd1, rd2, rd3, n1, n2, n3, rdn
    real(dp) :: tu1, tu2, tu3, tt1, tt2, tt3
    real(dp) :: rdt, rdu, nu_val
    real(dp) :: fac_U, fac_T, fac_D, fac_S
    real(dp) :: u_tmp(3), sig_tmp(6)
    real(dp) :: D_xx, D_yy, D_zz, D_xy, D_yz, D_zx
    real(dp) :: S_xx, S_yy, S_zz, S_xy, S_yz, S_zx
    real(dp) :: c4, c6, pr4, nu, const_U, const_T, const_D, const_S

    if (l <= 0) return

    ! Allocate internal point result arrays
    if (allocated(u_int)) deallocate(u_int)
    if (allocated(sigma_int)) deallocate(sigma_int)
    if (allocated(sigma_int_prin)) deallocate(sigma_int_prin)
    if (allocated(sigma_int_vm)) deallocate(sigma_int_vm)

    allocate(u_int(3, l), sigma_int(6, l), sigma_int_prin(3, l), sigma_int_vm(l))
    u_int = 0.0_dp
    sigma_int = 0.0_dp
    sigma_int_prin = 0.0_dp
    sigma_int_vm = 0.0_dp

    ! 1. Count total Gauss points across all boundary elements
    tot_gp = 0
    do ie = 1, ne
      if (nconec(ie) == 9) then
        tot_gp = tot_gp + 36
      else
        tot_gp = tot_gp + 16
      end if
    end do

    allocate(gp_x(tot_gp), gp_y(tot_gp), gp_z(tot_gp))
    allocate(gp_nx(tot_gp), gp_ny(tot_gp), gp_nz(tot_gp))
    allocate(gp_u(3, tot_gp), gp_t(3, tot_gp))

    ! 2. Precompute Gauss point positions, normals, and weighted field values (u * dGamma, t * dGamma)
    igp = 0
    do ie = 1, ne
      ncone = nconec(ie)
      if (ncone == 9) then
        npg = 6
      else
        npg = 4
      end if

      do i = 1, npg
        do j = 1, npg
          k = npg * (i - 1) + j
          igp = igp + 1

          if (ncone == 9) then
            wg = ome_quad(i) * ome_quad(j) * all_xja(k, ie)
          else
            wg = omei_tri(i) * omej_tri(j) * all_xja(k, ie)
          end if

          gp_x(igp) = all_xg1(k, ie)
          gp_y(igp) = all_xg2(k, ie)
          gp_z(igp) = all_xg3(k, ie)

          gp_nx(igp) = all_eta1(k, ie)
          gp_ny(igp) = all_eta2(k, ie)
          gp_nz(igp) = all_eta3(k, ie)

          ug = 0.0_dp
          tg = 0.0_dp
          do m = 1, ncone
            ik = kconec(m, ie)
            ff_val = all_ff(m, k, ie)
            ug(1) = ug(1) + ff_val * FI(3*ik - 2)
            ug(2) = ug(2) + ff_val * FI(3*ik - 1)
            ug(3) = ug(3) + ff_val * FI(3*ik)

            tg(1) = tg(1) + ff_val * DFI(3*m - 2, ie)
            tg(2) = tg(2) + ff_val * DFI(3*m - 1, ie)
            tg(3) = tg(3) + ff_val * DFI(3*m,     ie)
          end do

          gp_u(:, igp) = ug(:) * wg
          gp_t(:, igp) = tg(:) * wg
        end do
      end do
    end do

    ! Fundamental constants for Kelvin kernels
    nu = xnu
    c4 = 3.0_dp - 4.0_dp * nu
    c6 = 1.0_dp - 2.0_dp * nu
    pr4 = 1.0_dp - 4.0_dp * nu

    const_U = 1.0_dp / (16.0_dp * pi * (1.0_dp - nu) * rge)
    const_T = -1.0_dp / (8.0_dp * pi * (1.0_dp - nu))
    const_D = 1.0_dp / (8.0_dp * pi * (1.0_dp - nu))
    const_S = rge / (4.0_dp * pi * (1.0_dp - nu))

    ! 3. Compute internal points displacements and stresses in parallel
    !$omp parallel do private(ip, xi_x, xi_y, xi_z, igp, rx, ry, rz, r2, r, inv_r, inv_r2, inv_r3, &
    !$omp& rd1, rd2, rd3, n1, n2, n3, rdn, tu1, tu2, tu3, tt1, tt2, tt3, rdt, rdu, nu_val, &
    !$omp& fac_U, fac_T, fac_D, fac_S, u_tmp, sig_tmp, D_xx, D_yy, D_zz, D_xy, D_yz, D_zx, &
    !$omp& S_xx, S_yy, S_zz, S_xy, S_yz, S_zx) schedule(guided)
    do ip = 1, l
      xi_x = cx(ip)
      xi_y = cy(ip)
      xi_z = cz(ip)

      u_tmp = 0.0_dp
      sig_tmp = 0.0_dp

      do igp = 1, tot_gp
        rx = gp_x(igp) - xi_x
        ry = gp_y(igp) - xi_y
        rz = gp_z(igp) - xi_z

        r2 = rx*rx + ry*ry + rz*rz
        r = sqrt(r2)
        inv_r = 1.0_dp / r
        inv_r2 = inv_r * inv_r
        inv_r3 = inv_r2 * inv_r

        rd1 = rx * inv_r
        rd2 = ry * inv_r
        rd3 = rz * inv_r

        n1 = gp_nx(igp)
        n2 = gp_ny(igp)
        n3 = gp_nz(igp)

        rdn = rd1*n1 + rd2*n2 + rd3*n3

        tu1 = gp_u(1, igp); tu2 = gp_u(2, igp); tu3 = gp_u(3, igp)
        tt1 = gp_t(1, igp); tt2 = gp_t(2, igp); tt3 = gp_t(3, igp)

        rdt = rd1*tt1 + rd2*tt2 + rd3*tt3
        rdu = rd1*tu1 + rd2*tu2 + rd3*tu3
        nu_val = n1*tu1 + n2*tu2 + n3*tu3

        ! Displacement kernels: U* and T*
        fac_U = const_U * inv_r
        fac_T = const_T * inv_r2

        ! u_i += U_ik * t_k - T_ik * u_k
        u_tmp(1) = u_tmp(1) + fac_U * (c4 * tt1 + rd1 * rdt) - &
                              fac_T * (rdn * (c6 * tu1 + 3.0_dp * rd1 * rdu) + c6 * (nu_val * rd1 - n1 * rdu))
        u_tmp(2) = u_tmp(2) + fac_U * (c4 * tt2 + rd2 * rdt) - &
                              fac_T * (rdn * (c6 * tu2 + 3.0_dp * rd2 * rdu) + c6 * (nu_val * rd2 - n2 * rdu))
        u_tmp(3) = u_tmp(3) + fac_U * (c4 * tt3 + rd3 * rdt) - &
                              fac_T * (rdn * (c6 * tu3 + 3.0_dp * rd3 * rdu) + c6 * (nu_val * rd3 - n3 * rdu))

        ! Stress kernels: D* and S*
        fac_D = const_D * inv_r2
        fac_S = const_S * inv_r3

        D_xx = fac_D * (c6 * (2.0_dp * tt1 * rd1 - rdt) + 3.0_dp * rd1 * rd1 * rdt)
        D_yy = fac_D * (c6 * (2.0_dp * tt2 * rd2 - rdt) + 3.0_dp * rd2 * rd2 * rdt)
        D_zz = fac_D * (c6 * (2.0_dp * tt3 * rd3 - rdt) + 3.0_dp * rd3 * rd3 * rdt)
        D_xy = fac_D * (c6 * (tt1 * rd2 + tt2 * rd1)    + 3.0_dp * rd1 * rd2 * rdt)
        D_yz = fac_D * (c6 * (tt2 * rd3 + tt3 * rd2)    + 3.0_dp * rd2 * rd3 * rdt)
        D_zx = fac_D * (c6 * (tt3 * rd1 + tt1 * rd3)    + 3.0_dp * rd3 * rd1 * rdt)

        S_xx = fac_S * (3.0_dp * rdn * (c6 * rdu + 2.0_dp * nu * tu1 * rd1 - 5.0_dp * rd1 * rd1 * rdu) + &
                        6.0_dp * nu * rdu * n1 * rd1 + c6 * (3.0_dp * nu_val * rd1 * rd1 + 2.0_dp * n1 * tu1) - &
                        pr4 * nu_val)

        S_yy = fac_S * (3.0_dp * rdn * (c6 * rdu + 2.0_dp * nu * tu2 * rd2 - 5.0_dp * rd2 * rd2 * rdu) + &
                        6.0_dp * nu * rdu * n2 * rd2 + c6 * (3.0_dp * nu_val * rd2 * rd2 + 2.0_dp * n2 * tu2) - &
                        pr4 * nu_val)

        S_zz = fac_S * (3.0_dp * rdn * (c6 * rdu + 2.0_dp * nu * tu3 * rd3 - 5.0_dp * rd3 * rd3 * rdu) + &
                        6.0_dp * nu * rdu * n3 * rd3 + c6 * (3.0_dp * nu_val * rd3 * rd3 + 2.0_dp * n3 * tu3) - &
                        pr4 * nu_val)

        S_xy = fac_S * (3.0_dp * rdn * (nu * (tu1 * rd2 + tu2 * rd1) - 5.0_dp * rd1 * rd2 * rdu) + &
                        3.0_dp * nu * rdu * (n1 * rd2 + n2 * rd1) + &
                        c6 * (3.0_dp * nu_val * rd1 * rd2 + n1 * tu2 + n2 * tu1))

        S_yz = fac_S * (3.0_dp * rdn * (nu * (tu2 * rd3 + tu3 * rd2) - 5.0_dp * rd2 * rd3 * rdu) + &
                        3.0_dp * nu * rdu * (n2 * rd3 + n3 * rd2) + &
                        c6 * (3.0_dp * nu_val * rd2 * rd3 + n2 * tu3 + n3 * tu2))

        S_zx = fac_S * (3.0_dp * rdn * (nu * (tu3 * rd1 + tu1 * rd3) - 5.0_dp * rd3 * rd1 * rdu) + &
                        3.0_dp * nu * rdu * (n3 * rd1 + n1 * rd3) + &
                        c6 * (3.0_dp * nu_val * rd3 * rd1 + n3 * tu1 + n1 * tu3))

        ! sigma_ij += D_kij * t_k - S_kij * u_k
        sig_tmp(1) = sig_tmp(1) + (D_xx - S_xx)
        sig_tmp(2) = sig_tmp(2) + (D_yy - S_yy)
        sig_tmp(3) = sig_tmp(3) + (D_zz - S_zz)
        sig_tmp(4) = sig_tmp(4) + (D_xy - S_xy)
        sig_tmp(5) = sig_tmp(5) + (D_yz - S_yz)
        sig_tmp(6) = sig_tmp(6) + (D_zx - S_zx)
      end do

      u_int(1:3, ip) = u_tmp(1:3)
      sigma_int(1:6, ip) = sig_tmp(1:6)

      ! Synchronize with legacy dsol array
      dsol(3*ip - 2) = u_tmp(1)
      dsol(3*ip - 1) = u_tmp(2)
      dsol(3*ip    ) = u_tmp(3)

      call calc_principal_and_vm(sigma_int(:, ip), sigma_int_prin(:, ip), sigma_int_vm(ip))
    end do
    !$omp end parallel do

    deallocate(gp_x, gp_y, gp_z, gp_nx, gp_ny, gp_nz, gp_u, gp_t)
  end subroutine compute_internal_points

  !=============================================================================
  ! output_internal_points
  ! Outputs internal points displacements and stresses to ipr and internal_points.dat
  !=============================================================================
  subroutine output_internal_points()
    implicit none
    integer :: ip, k
    integer, parameter :: iunit_int = 27

    if (l <= 0) return

    ! Write Internal Points Stresses to main output file (ipr)
    write(ipr, '(//A)') ' ' // repeat('*', 79)
    write(ipr, '(A)') ' INTERNAL POINTS STRESSES'
    write(ipr, '(A)') ' ' // repeat('*', 79)
    write(ipr, 100)
100 format(1X, 'POINT', 5X, 'SIGMA_XX', 8X, 'SIGMA_YY', 8X, 'SIGMA_ZZ', 8X, &
           'SIGMA_XY', 8X, 'SIGMA_YZ', 8X, 'SIGMA_ZX', 7X, 'VON_MISES', 9X, &
           'SIGMA_1', 9X, 'SIGMA_2', 9X, 'SIGMA_3'/)

    do ip = 1, l
      write(ipr, 110) ip, (sigma_int(k, ip), k = 1, 6), sigma_int_vm(ip), &
                          (sigma_int_prin(k, ip), k = 1, 3)
    end do
110 format(1X, I5, 10(1X, E15.7))

    ! Write dedicated internal_points.dat file for post-processing and plotting
    open(unit=iunit_int, file=trim(fileintpts), status='replace', action='write')
    write(iunit_int, '(A)') '# INTERNAL POINTS DISPLACEMENTS AND STRESSES (QUADESTH 3D ELASTOSTATICS)'
    write(iunit_int, '(A, I6)') '# Total Internal Points: ', l
    write(iunit_int, '(A)') '# Col 1: Point ID'
    write(iunit_int, '(A)') '# Col 2-4: X, Y, Z coordinates'
    write(iunit_int, '(A)') '# Col 5-7: u_x, u_y, u_z displacements'
    write(iunit_int, '(A)') '# Col 8-13: Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_yz, Sigma_zx'
    write(iunit_int, '(A)') '# Col 14: Von Mises Equivalent Stress'
    write(iunit_int, '(A)') '# Col 15-17: Principal Stresses (Sigma_1 >= Sigma_2 >= Sigma_3)'
    write(iunit_int, '(A)') '# ' // repeat('-', 140)

    do ip = 1, l
      write(iunit_int, 300) ip, cx(ip), cy(ip), cz(ip), &
                           (u_int(k, ip), k = 1, 3), &
                           (sigma_int(k, ip), k = 1, 6), sigma_int_vm(ip), &
                           (sigma_int_prin(k, ip), k = 1, 3)
    end do
300 format(I6, 3(1X, E15.7), 3(1X, E15.7), 6(1X, E15.7), 4(1X, E15.7))

    close(iunit_int)
    write(*, '(A)') ' Internal points exported to: ' // trim(fileintpts)
  end subroutine output_internal_points

  !=============================================================================
  ! deallocate_stress_data
  ! Frees all allocated memory in stress_mod
  !=============================================================================
  subroutine deallocate_stress_data()
    implicit none
    if (allocated(sigma_elem)) deallocate(sigma_elem)
    if (allocated(sigma_node)) deallocate(sigma_node)
    if (allocated(sigma_principal)) deallocate(sigma_principal)
    if (allocated(sigma_vm)) deallocate(sigma_vm)
    if (allocated(u_int)) deallocate(u_int)
    if (allocated(sigma_int)) deallocate(sigma_int)
    if (allocated(sigma_int_prin)) deallocate(sigma_int_prin)
    if (allocated(sigma_int_vm)) deallocate(sigma_int_vm)
  end subroutine deallocate_stress_data

end module stress_mod
