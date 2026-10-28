module stress_mod
  use data_mod
  implicit none

  ! Dynamic Arrays for boundary stresses
  ! Cauchy stress ordering: 1=xx, 2=yy, 3=zz, 4=xy, 5=yz, 6=zx
  real(dp), allocatable, public :: sigma_elem(:, :, :)   ! (6, 9, ne): element-nodal stresses
  real(dp), allocatable, public :: sigma_node(:, :)      ! (6, n): nodal averaged stresses
  real(dp), allocatable, public :: sigma_principal(:, :) ! (3, n): principal stresses (s1 >= s2 >= s3)
  real(dp), allocatable, public :: sigma_vm(:)           ! (n): Von Mises equivalent stress

  character(len=64), public :: filestress = 'stresses.dat'

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
  ! deallocate_stress_data
  ! Frees all allocated memory in stress_mod
  !=============================================================================
  subroutine deallocate_stress_data()
    implicit none
    if (allocated(sigma_elem)) deallocate(sigma_elem)
    if (allocated(sigma_node)) deallocate(sigma_node)
    if (allocated(sigma_principal)) deallocate(sigma_principal)
    if (allocated(sigma_vm)) deallocate(sigma_vm)
  end subroutine deallocate_stress_data

end module stress_mod
