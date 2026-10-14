module assembly_mod
  use data_mod
  implicit none

  ! Gauss-Legendre quadrature constants
  real(dp), parameter :: gi_quad(6) = [ &
       0.238619186083197_dp, -0.238619186083197_dp, &
       0.661209386466265_dp, -0.661209386466265_dp, &
       0.932469514203152_dp, -0.932469514203152_dp ]

  real(dp), parameter :: ome_quad(6) = [ &
       0.467913934572691_dp, 0.467913934572691_dp, &
       0.360761573048139_dp, 0.360761573048139_dp, &
       0.171324492379170_dp, 0.171324492379170_dp ]

  real(dp), parameter :: rr_tri(4) = [ 0.0694318422_dp, 0.3300094782_dp, 0.6699905218_dp, 0.9305681558_dp ]
  real(dp), parameter :: ss_tri(4) = [ 0.0571041961_dp, 0.2768430136_dp, 0.5835904324_dp, 0.8602401357_dp ]
  real(dp), parameter :: omei_tri(4) = [ 0.1739274226_dp, 0.3260725774_dp, 0.3260725774_dp, 0.1739274226_dp ]
  real(dp), parameter :: omej_tri(4) = [ 0.1355069134_dp, 0.2034645680_dp, 0.1298475476_dp, 0.0311809709_dp ]

  ! Precomputed Gauss geometry data for elements
  real(dp), allocatable :: all_eta1(:,:), all_eta2(:,:), all_eta3(:,:), all_xja(:,:)
  real(dp), allocatable :: all_xg1(:,:), all_xg2(:,:), all_xg3(:,:)
  real(dp), allocatable :: all_ff(:,:,:)

contains

  subroutine precompute_element_geometry()
    implicit none
    integer :: iej, m, k, npts
    real(dp) :: xj1(9), xj2(9), xj3(9)

    npts = 36
    if (allocated(all_eta1)) deallocate(all_eta1, all_eta2, all_eta3, all_xja, all_xg1, all_xg2, all_xg3, all_ff)
    allocate(all_eta1(npts, ne), all_eta2(npts, ne), all_eta3(npts, ne), all_xja(npts, ne))
    allocate(all_xg1(npts, ne), all_xg2(npts, ne), all_xg3(npts, ne))
    allocate(all_ff(9, npts, ne))

    do iej = 1, ne
      do m = 1, nconec(iej)
        k = kconec(m, iej)
        xj1(m) = x(k)
        xj2(m) = y(k)
        xj3(m) = z(k)
      end do

      call etaxja(xj1, xj2, xj3, all_eta1(:, iej), all_eta2(:, iej), all_eta3(:, iej), &
                  all_xja(:, iej), all_xg1(:, iej), all_xg2(:, iej), all_xg3(:, iej), &
                  all_ff(:, :, iej), nconec(iej))
    end do
  end subroutine precompute_element_geometry

  subroutine ghmat8()
    implicit none
    integer :: in, iej, m, k, nodo, i, ik, ind, ld, kode0, j
    real(dp) :: xi(3), xj1(9), xj2(9), xj3(9)
    real(dp) :: h_elem(3, 27), g_elem(3, 27)
    real(dp) :: d1(3, 3), d2(3, 3), vp

    ! Reset system matrix and RHS vector
    A = 0.0_dp
    FI = 0.0_dp

    ! Parallelize assembly over collocation points using OpenMP
    !$omp parallel do default(none) &
    !$omp& shared(n, ne, x, y, z, nconec, kconec, nesim, ke, kode, bc, rge, A, FI) &
    !$omp& shared(all_eta1, all_eta2, all_eta3, all_xja, all_xg1, all_xg2, all_xg3, all_ff) &
    !$omp& private(in, xi, iej, xj1, xj2, xj3, nodo, i, ik, ind, k, ld, kode0, m, j, vp) &
    !$omp& private(h_elem, g_elem, d1, d2) &
    !$omp& schedule(dynamic)
    do in = 1, n
      xi(1) = x(in)
      xi(2) = y(in)
      xi(3) = z(in)

      do iej = 1, ne
        do m = 1, nconec(iej)
          k = kconec(m, iej)
          xj1(m) = x(k)
          xj2(m) = y(k)
          xj3(m) = z(k)
        end do

        nodo = 0
        do i = 1, nconec(iej)
          if (in == kconec(i, iej)) then
            nodo = i
            exit
          end if
        end do

        if (nodo > 0) then
          call locin8(xj1, xj2, xj3, h_elem, g_elem, nodo, nconec(iej))
        else
          call extin8(all_eta1(:, iej), all_eta2(:, iej), all_eta3(:, iej), all_xja(:, iej), &
                      all_xg1(:, iej), all_xg2(:, iej), all_xg3(:, iej), all_ff(:, :, iej), &
                      xi, h_elem, g_elem, nconec(iej))
        end if

        d1 = 0.0_dp
        d2 = 0.0_dp

        if (nodo == 0) then
          do j = 1, 3
            do i = 1, 3
              do k = 1, nconec(iej)
                ind = 3*k - 3 + j
                d1(i, j) = d1(i, j) - h_elem(i, ind)
              end do
            end do
          end do
        else
          do j = 1, 3
            do i = 1, 3
              do k = 1, nconec(iej)
                if (k /= nodo) then
                  ind = 3*k - 3 + j
                  d2(i, j) = d2(i, j) - h_elem(i, ind)
                end if
              end do
            end do
          end do
        end if

        ! Form system equations according to KODE
        do ld = 1, 3
          if (kode(3*iej - 3 + ld) == 0) then
            ! Displacement prescribed on component ld of element iej
            do k = 1, nconec(iej)
              ik = kconec(k, iej)
              A(3*in-2, 3*ik-3+ld) = A(3*in-2, 3*ik-3+ld) - g_elem(1, 3*k-3+ld) * rge
              A(3*in-1, 3*ik-3+ld) = A(3*in-1, 3*ik-3+ld) - g_elem(2, 3*k-3+ld) * rge
              A(3*in  , 3*ik-3+ld) = A(3*in  , 3*ik-3+ld) - g_elem(3, 3*k-3+ld) * rge

              if (in /= ik) then
                FI(3*in-2) = FI(3*in-2) - h_elem(1, 3*k-3+ld) * bc(3*k-3+ld, iej)
                FI(3*in-1) = FI(3*in-1) - h_elem(2, 3*k-3+ld) * bc(3*k-3+ld, iej)
                FI(3*in  ) = FI(3*in  ) - h_elem(3, 3*k-3+ld) * bc(3*k-3+ld, iej)
              else
                FI(3*in-2) = FI(3*in-2) - d2(1, ld) * bc(3*k-3+ld, iej)
                FI(3*in-1) = FI(3*in-1) - d2(2, ld) * bc(3*k-3+ld, iej)
                FI(3*in  ) = FI(3*in  ) - d2(3, ld) * bc(3*k-3+ld, iej)
              end if
            end do
          else
            ! Traction prescribed on component ld of element iej
            do k = 1, nconec(iej)
              ik = kconec(k, iej)
              kode0 = 0
              do i = 1, nesim(ik)
                if (kode(3*ke(i, ik) - 3 + ld) == 0) then
                  do m = 1, nconec(ke(i, ik))
                    if (ik == kconec(m, ke(i, ik))) then
                      vp = bc(3*m-3+ld, ke(i, ik))
                      kode0 = 1
                      exit
                    end if
                  end do
                  if (kode0 /= 0) exit
                end if
              end do

              if (kode0 == 0) then
                FI(3*in-2) = FI(3*in-2) + g_elem(1, 3*k-3+ld) * bc(3*k-3+ld, iej)
                FI(3*in-1) = FI(3*in-1) + g_elem(2, 3*k-3+ld) * bc(3*k-3+ld, iej)
                FI(3*in  ) = FI(3*in  ) + g_elem(3, 3*k-3+ld) * bc(3*k-3+ld, iej)

                if (in /= ik) then
                  A(3*in-2, 3*ik-3+ld) = A(3*in-2, 3*ik-3+ld) + h_elem(1, 3*k-3+ld)
                  A(3*in-1, 3*ik-3+ld) = A(3*in-1, 3*ik-3+ld) + h_elem(2, 3*k-3+ld)
                  A(3*in  , 3*ik-3+ld) = A(3*in  , 3*ik-3+ld) + h_elem(3, 3*k-3+ld)
                else
                  A(3*in-2, 3*in-3+ld) = A(3*in-2, 3*in-3+ld) + d2(1, ld)
                  A(3*in-1, 3*in-3+ld) = A(3*in-1, 3*in-3+ld) + d2(2, ld)
                  A(3*in  , 3*in-3+ld) = A(3*in  , 3*in-3+ld) + d2(3, ld)
                end if
              else
                if (in /= ik) then
                  FI(3*in-2) = FI(3*in-2) + g_elem(1, 3*k-3+ld) * bc(3*k-3+ld, iej) - h_elem(1, 3*k-3+ld) * vp
                  FI(3*in-1) = FI(3*in-1) + g_elem(2, 3*k-3+ld) * bc(3*k-3+ld, iej) - h_elem(2, 3*k-3+ld) * vp
                  FI(3*in  ) = FI(3*in  ) + g_elem(3, 3*k-3+ld) * bc(3*k-3+ld, iej) - h_elem(3, 3*k-3+ld) * vp
                else
                  FI(3*in-2) = FI(3*in-2) + g_elem(1, 3*k-3+ld) * bc(3*k-3+ld, iej) - d2(1, ld) * vp
                  FI(3*in-1) = FI(3*in-1) + g_elem(2, 3*k-3+ld) * bc(3*k-3+ld, iej) - d2(2, ld) * vp
                  FI(3*in  ) = FI(3*in  ) + g_elem(3, 3*k-3+ld) * bc(3*k-3+ld, iej) - d2(3, ld) * vp
                end if
              end if
            end do
          end if
        end do

        if (nodo == 0) then
          do ld = 1, 3
            kode0 = 0
            do i = 1, nesim(in)
              if (kode(3*ke(i, in) - 3 + ld) == 0) then
                do m = 1, nconec(ke(i, in))
                  if (in == kconec(m, ke(i, in))) then
                    vp = bc(3*m-3+ld, ke(i, in))
                    kode0 = 1
                    exit
                  end if
                end do
                if (kode0 /= 0) exit
              end if
            end do

            if (kode0 == 0) then
              A(3*in-2, 3*in-3+ld) = A(3*in-2, 3*in-3+ld) + d1(1, ld)
              A(3*in-1, 3*in-3+ld) = A(3*in-1, 3*in-3+ld) + d1(2, ld)
              A(3*in  , 3*in-3+ld) = A(3*in  , 3*in-3+ld) + d1(3, ld)
            else
              FI(3*in-2) = FI(3*in-2) - d1(1, ld) * vp
              FI(3*in-1) = FI(3*in-1) - d1(2, ld) * vp
              FI(3*in  ) = FI(3*in  ) - d1(3, ld) * vp
            end if
          end do
        end if

      end do ! iej = 1, ne
    end do ! in = 1, n
    !$omp end parallel do
  end subroutine ghmat8

  subroutine etaxja(xj1, xj2, xj3, eta1, eta2, eta3, xja, xg1, xg2, xg3, ff, ncone)
    implicit none
    integer, intent(in) :: ncone
    real(dp), intent(in)  :: xj1(9), xj2(9), xj3(9)
    real(dp), intent(out) :: eta1(36), eta2(36), eta3(36), xja(36)
    real(dp), intent(out) :: xg1(36), xg2(36), xg3(36), ff(9, 36)

    integer :: npg, i, j, k, ik, in_idx, jn_idx, jk
    real(dp) :: g1, g2, g3
    real(dp) :: fd1(9), fd2(9), dc(9, 9), dct(9, 9)
    real(dp) :: aux(9), bux(9), cux(9), s1, s2, s3

    npg = 6
    if (ncone < 9) npg = 4

    eta1 = 0.0_dp; eta2 = 0.0_dp; eta3 = 0.0_dp; xja = 0.0_dp
    xg1 = 0.0_dp;  xg2 = 0.0_dp;  xg3 = 0.0_dp;  ff = 0.0_dp

    do i = 1, npg
      do j = 1, npg
        k = npg * (i - 1) + j

        if (ncone == 9) then
          g1 = gi_quad(i)
          g2 = gi_quad(j)
          ff(1, k) = 0.25_dp * g1 * g2 * (g1 - 1.0_dp) * (g2 - 1.0_dp)
          ff(2, k) = 0.5_dp * (1.0_dp - g1**2) * g2 * (g2 - 1.0_dp)
          ff(3, k) = 0.25_dp * g1 * g2 * (g1 + 1.0_dp) * (g2 - 1.0_dp)
          ff(4, k) = 0.5_dp * g1 * (g1 + 1.0_dp) * (1.0_dp - g2**2)
          ff(5, k) = 0.25_dp * g1 * g2 * (g1 + 1.0_dp) * (g2 + 1.0_dp)
          ff(6, k) = 0.5_dp * (1.0_dp - g1**2) * g2 * (g2 + 1.0_dp)
          ff(7, k) = 0.25_dp * g1 * g2 * (g1 - 1.0_dp) * (g2 + 1.0_dp)
          ff(8, k) = 0.5_dp * g1 * (g1 - 1.0_dp) * (1.0_dp - g2**2)
          ff(9, k) = (1.0_dp - g1**2) * (1.0_dp - g2**2)

          fd1(1) = 0.25_dp * g2 * (g2 - 1.0_dp) * (2.0_dp * g1 - 1.0_dp)
          fd1(2) = -g1 * g2 * (g2 - 1.0_dp)
          fd1(3) = 0.25_dp * g2 * (g2 - 1.0_dp) * (2.0_dp * g1 + 1.0_dp)
          fd1(4) = 0.5_dp * (1.0_dp - g2**2) * (2.0_dp * g1 + 1.0_dp)
          fd1(5) = 0.25_dp * g2 * (g2 + 1.0_dp) * (2.0_dp * g1 + 1.0_dp)
          fd1(6) = -g1 * g2 * (g2 + 1.0_dp)
          fd1(7) = 0.25_dp * g2 * (g2 + 1.0_dp) * (2.0_dp * g1 - 1.0_dp)
          fd1(8) = 0.5_dp * (1.0_dp - g2**2) * (2.0_dp * g1 - 1.0_dp)
          fd1(9) = -2.0_dp * g1 * (1.0_dp - g2**2)

          fd2(1) = 0.25_dp * g1 * (g1 - 1.0_dp) * (2.0_dp * g2 - 1.0_dp)
          fd2(2) = 0.5_dp * (1.0_dp - g1**2) * (2.0_dp * g2 - 1.0_dp)
          fd2(3) = 0.25_dp * g1 * (g1 + 1.0_dp) * (2.0_dp * g2 - 1.0_dp)
          fd2(4) = -g1 * g2 * (g1 + 1.0_dp)
          fd2(5) = 0.25_dp * g1 * (g1 + 1.0_dp) * (2.0_dp * g2 + 1.0_dp)
          fd2(6) = 0.5_dp * (1.0_dp - g1**2) * (2.0_dp * g2 + 1.0_dp)
          fd2(7) = 0.25_dp * g1 * (g1 - 1.0_dp) * (2.0_dp * g2 + 1.0_dp)
          fd2(8) = -g1 * g2 * (g1 - 1.0_dp)
          fd2(9) = -2.0_dp * g2 * (1.0_dp - g1**2)
        else
          g1 = ss_tri(j)
          g2 = rr_tri(i) * (1.0_dp - ss_tri(j))
          g3 = 1.0_dp - g1 - g2
          ff(1, k) = g1 * (2.0_dp * g1 - 1.0_dp)
          ff(2, k) = g2 * (2.0_dp * g2 - 1.0_dp)
          ff(3, k) = g3 * (2.0_dp * g3 - 1.0_dp)
          ff(4, k) = 4.0_dp * g1 * g2
          ff(5, k) = 4.0_dp * g2 * g3
          ff(6, k) = 4.0_dp * g3 * g1

          fd1(1) = 4.0_dp * g1 - 1.0_dp
          fd1(2) = 0.0_dp
          fd1(3) = 1.0_dp - 4.0_dp * g3
          fd1(4) = 4.0_dp * g2
          fd1(5) = -4.0_dp * g2
          fd1(6) = 4.0_dp * (g3 - g1)

          fd2(1) = 0.0_dp
          fd2(2) = 4.0_dp * g2 - 1.0_dp
          fd2(3) = 1.0_dp - 4.0_dp * g3
          fd2(4) = 4.0_dp * g1
          fd2(5) = 4.0_dp * (g3 - g2)
          fd2(6) = -4.0_dp * g1
        end if

        do ik = 1, ncone
          xg1(k) = xg1(k) + ff(ik, k) * xj1(ik)
          xg2(k) = xg2(k) + ff(ik, k) * xj2(ik)
          xg3(k) = xg3(k) + ff(ik, k) * xj3(ik)
        end do

        do in_idx = 1, ncone
          do jn_idx = 1, ncone
            dc(in_idx, jn_idx)  = fd1(in_idx) * fd2(jn_idx)
            dct(in_idx, jn_idx) = fd2(in_idx) * fd1(jn_idx)
          end do
        end do

        aux = 0.0_dp; bux = 0.0_dp; cux = 0.0_dp
        do jk = 1, ncone
          do ik = 1, ncone
            aux(jk) = aux(jk) + xj2(ik) * (dc(ik, jk) - dct(ik, jk))
            bux(jk) = bux(jk) + xj3(ik) * (dc(ik, jk) - dct(ik, jk))
            cux(jk) = cux(jk) + xj1(ik) * (dc(ik, jk) - dct(ik, jk))
          end do
        end do

        s1 = sum(aux(1:ncone) * xj3(1:ncone))
        s2 = sum(bux(1:ncone) * xj1(1:ncone))
        s3 = sum(cux(1:ncone) * xj2(1:ncone))

        xja(k) = sqrt(s1**2 + s2**2 + s3**2)
        eta1(k) = s1 / xja(k)
        eta2(k) = s2 / xja(k)
        eta3(k) = s3 / xja(k)
      end do
    end do
  end subroutine etaxja

  subroutine extin8(eta1, eta2, eta3, xja, xg1, xg2, xg3, ff, xi, h_elem, g_elem, ncone)
    implicit none
    integer, intent(in) :: ncone
    real(dp), intent(in) :: eta1(36), eta2(36), eta3(36), xja(36)
    real(dp), intent(in) :: xg1(36), xg2(36), xg3(36), ff(9, 36), xi(3)
    real(dp), intent(out) :: h_elem(3, 27), g_elem(3, 27)

    integer :: npg, i, j, k, ik, kk, ii, jj, ind
    real(dp) :: eta(3), p12(9)
    real(dp) :: u(3, 3), t(3, 3)

    h_elem = 0.0_dp
    g_elem = 0.0_dp

    npg = 6
    if (ncone < 9) npg = 4

    do i = 1, npg
      do j = 1, npg
        k = npg * (i - 1) + j
        eta(1) = eta1(k)
        eta(2) = eta2(k)
        eta(3) = eta3(k)

        do ik = 1, ncone
          if (ncone == 9) then
            p12(ik) = ff(ik, k) * xja(k) * ome_quad(i) * ome_quad(j)
          else
            p12(ik) = ff(ik, k) * xja(k) * omei_tri(i) * omej_tri(j)
          end if
        end do

        call funkelvin(u, t, xi(1), xi(2), xi(3), eta, xg1(k), xg2(k), xg3(k))

        do kk = 1, ncone
          do ii = 1, 3
            do jj = 1, 3
              ind = 3 * kk + ii - 3
              g_elem(jj, ind) = g_elem(jj, ind) + u(jj, ii) * p12(kk)
              h_elem(jj, ind) = h_elem(jj, ind) + t(jj, ii) * p12(kk)
            end do
          end do
        end do
      end do
    end do
  end subroutine extin8

  subroutine locin8(xj1, xj2, xj3, h_elem, g_elem, nodo, ncone)
    implicit none
    integer, intent(in) :: nodo, ncone
    real(dp), intent(in) :: xj1(9), xj2(9), xj3(9)
    real(dp), intent(out) :: h_elem(3, 27), g_elem(3, 27)

    integer :: k, m
    real(dp) :: xi(3), v1(3), v2(3), area
    real(dp) :: gi1(9), gi2(9), gii1(15), gii2(15)

    h_elem = 0.0_dp
    g_elem = 0.0_dp

    if (ncone == 9) then
      gi1(1) = -1.0_dp; gi1(2) =  0.0_dp; gi1(3) =  1.0_dp; gi1(4) =  1.0_dp
      gi1(5) =  1.0_dp; gi1(6) =  0.0_dp; gi1(7) = -1.0_dp; gi1(8) = -1.0_dp; gi1(9) = 0.0_dp

      gi2(1) = -1.0_dp; gi2(2) = -1.0_dp; gi2(3) = -1.0_dp; gi2(4) =  0.0_dp
      gi2(5) =  1.0_dp; gi2(6) =  1.0_dp; gi2(7) =  1.0_dp; gi2(8) =  0.0_dp; gi2(9) = 0.0_dp

      gii1(1:8) = gi1(1:8)
      gii2(1:8) = gi2(1:8)
      gii1(9:15) = gi1(1:7)
      gii2(9:15) = gi2(1:7)

      v1(1) = gi1(nodo); v2(1) = gi2(nodo)
      xi(1) = xj1(nodo); xi(2) = xj2(nodo); xi(3) = xj3(nodo)

      select case (nodo)
      case (1, 3, 5, 7)
        area = 2.0_dp
        do k = 1, 2
          m = 2 * k
          v1(2) = gii1(nodo + m)
          v2(2) = gii2(nodo + m)
          v1(3) = gii1(nodo + m + 2)
          v2(3) = gii2(nodo + m + 2)
          call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, nodo, 9)
        end do
      case (2, 4, 6, 8)
        do k = 1, 3
          m = 2 * k - 1
          v1(2) = gii1(nodo + m)
          v2(2) = gii2(nodo + m)
          v1(3) = gii1(nodo + m + 2)
          v2(3) = gii2(nodo + m + 2)
          area = 0.5_dp * (v1(2)*v2(3) + v1(1)*v2(2) + v1(3)*v2(1) - &
                           v1(2)*v2(1) - v1(3)*v2(2) - v1(1)*v2(3))
          call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, nodo, 9)
        end do
      case (9)
        area = 1.0_dp
        do k = 1, 4
          m = 2 * k - 1
          v1(2) = gii1(m)
          v2(2) = gii2(m)
          v1(3) = gii1(m + 2)
          v2(3) = gii2(m + 2)
          call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, nodo, 9)
        end do
      end select
    else
      gi1(1) = 1.0_dp; gi1(2) = 0.0_dp; gi1(3) = 0.0_dp; gi1(4) = 0.5_dp; gi1(5) = 0.0_dp; gi1(6) = 0.5_dp
      gi2(1) = 0.0_dp; gi2(2) = 1.0_dp; gi2(3) = 0.0_dp; gi2(4) = 0.5_dp; gi2(5) = 0.5_dp; gi2(6) = 0.0_dp

      v1(1) = gi1(nodo); v2(1) = gi2(nodo)
      xi(1) = xj1(nodo); xi(2) = xj2(nodo); xi(3) = xj3(nodo)

      select case (nodo)
      case (1)
        area = 0.5_dp
        v1(2) = gi1(2); v2(2) = gi2(2)
        v1(3) = gi1(3); v2(3) = gi2(3)
        call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, 1, 6)
      case (2)
        area = 0.5_dp
        v1(2) = gi1(3); v2(2) = gi2(3)
        v1(3) = gi1(1); v2(3) = gi2(1)
        call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, 2, 6)
      case (3)
        area = 0.5_dp
        v1(2) = gi1(1); v2(2) = gi2(1)
        v1(3) = gi1(2); v2(3) = gi2(2)
        call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, 3, 6)
      case (4)
        area = 0.25_dp
        v1(2) = gi1(2); v2(2) = gi2(2); v1(3) = gi1(3); v2(3) = gi2(3)
        call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, 4, 6)
        v1(2) = gi1(3); v2(2) = gi2(3); v1(3) = gi1(1); v2(3) = gi2(1)
        call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, 4, 6)
      case (5)
        area = 0.25_dp
        v1(2) = gi1(3); v2(2) = gi2(3); v1(3) = gi1(1); v2(3) = gi2(1)
        call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, 5, 6)
        v1(2) = gi1(1); v2(2) = gi2(1); v1(3) = gi1(2); v2(3) = gi2(2)
        call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, 5, 6)
      case (6)
        area = 0.25_dp
        v1(2) = gi1(1); v2(2) = gi2(1); v1(3) = gi1(2); v2(3) = gi2(2)
        call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, 6, 6)
        v1(2) = gi1(2); v2(2) = gi2(2); v1(3) = gi1(3); v2(3) = gi2(3)
        call triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, 6, 6)
      end select
    end if
  end subroutine locin8

  subroutine triloc8(xi, xj1, xj2, xj3, h_elem, g_elem, area, v1, v2, nodo, ncone)
    implicit none
    integer, intent(in) :: nodo, ncone
    real(dp), intent(in) :: xi(3), xj1(9), xj2(9), xj3(9), area, v1(3), v2(3)
    real(dp), intent(inout) :: h_elem(3, 27), g_elem(3, 27)

    integer :: i, j, ik, in_idx, jn_idx, jk, kk, ii, jj, ind
    real(dp) :: ss1, ss2, g1, g2, g3
    real(dp) :: f(9), fd1(9), fd2(9), dc(9, 9), dct(9, 9)
    real(dp) :: aux(9), bux(9), cux(9), s1, s2, s3, xja, eta(3)
    real(dp) :: x_pt(3), p12(9)
    real(dp) :: u(3, 3), t(3, 3)

    do i = 1, 6
      do j = 1, 6
        ss1 = 0.5_dp * (gi_quad(i) + 1.0_dp)
        ss2 = 0.5_dp * (gi_quad(j) + 1.0_dp)
        g1 = (1.0_dp - ss1)*v1(1) + ss1*(1.0_dp - ss2)*v1(2) + ss1*ss2*v1(3)
        g2 = (1.0_dp - ss1)*v2(1) + ss1*(1.0_dp - ss2)*v2(2) + ss1*ss2*v2(3)

        if (ncone == 9) then
          f(1) = 0.25_dp * g1 * g2 * (g1 - 1.0_dp) * (g2 - 1.0_dp)
          f(2) = 0.5_dp * (1.0_dp - g1**2) * g2 * (g2 - 1.0_dp)
          f(3) = 0.25_dp * g1 * g2 * (g1 + 1.0_dp) * (g2 - 1.0_dp)
          f(4) = 0.5_dp * g1 * (g1 + 1.0_dp) * (1.0_dp - g2**2)
          f(5) = 0.25_dp * g1 * g2 * (g1 + 1.0_dp) * (g2 + 1.0_dp)
          f(6) = 0.5_dp * (1.0_dp - g1**2) * g2 * (g2 + 1.0_dp)
          f(7) = 0.25_dp * g1 * g2 * (g1 - 1.0_dp) * (g2 + 1.0_dp)
          f(8) = 0.5_dp * g1 * (g1 - 1.0_dp) * (1.0_dp - g2**2)
          f(9) = (1.0_dp - g1**2) * (1.0_dp - g2**2)

          fd1(1) = 0.25_dp * g2 * (g2 - 1.0_dp) * (2.0_dp * g1 - 1.0_dp)
          fd1(2) = -g1 * g2 * (g2 - 1.0_dp)
          fd1(3) = 0.25_dp * g2 * (g2 - 1.0_dp) * (2.0_dp * g1 + 1.0_dp)
          fd1(4) = 0.5_dp * (1.0_dp - g2**2) * (2.0_dp * g1 + 1.0_dp)
          fd1(5) = 0.25_dp * g2 * (g2 + 1.0_dp) * (2.0_dp * g1 + 1.0_dp)
          fd1(6) = -g1 * g2 * (g2 + 1.0_dp)
          fd1(7) = 0.25_dp * g2 * (g2 + 1.0_dp) * (2.0_dp * g1 - 1.0_dp)
          fd1(8) = 0.5_dp * (1.0_dp - g2**2) * (2.0_dp * g1 - 1.0_dp)
          fd1(9) = -2.0_dp * g1 * (1.0_dp - g2**2)

          fd2(1) = 0.25_dp * g1 * (g1 - 1.0_dp) * (2.0_dp * g2 - 1.0_dp)
          fd2(2) = 0.5_dp * (1.0_dp - g1**2) * (2.0_dp * g2 - 1.0_dp)
          fd2(3) = 0.25_dp * g1 * (g1 + 1.0_dp) * (2.0_dp * g2 - 1.0_dp)
          fd2(4) = -g1 * g2 * (g1 + 1.0_dp)
          fd2(5) = 0.25_dp * g1 * (g1 + 1.0_dp) * (2.0_dp * g2 + 1.0_dp)
          fd2(6) = 0.5_dp * (1.0_dp - g1**2) * (2.0_dp * g2 + 1.0_dp)
          fd2(7) = 0.25_dp * g1 * (g1 - 1.0_dp) * (2.0_dp * g2 + 1.0_dp)
          fd2(8) = -g1 * g2 * (g1 - 1.0_dp)
          fd2(9) = -2.0_dp * g2 * (1.0_dp - g1**2)
        else
          g3 = 1.0_dp - g1 - g2
          f(1) = g1 * (2.0_dp * g1 - 1.0_dp)
          f(2) = g2 * (2.0_dp * g2 - 1.0_dp)
          f(3) = g3 * (2.0_dp * g3 - 1.0_dp)
          f(4) = 4.0_dp * g1 * g2
          f(5) = 4.0_dp * g2 * g3
          f(6) = 4.0_dp * g3 * g1

          fd1(1) = 4.0_dp * g1 - 1.0_dp
          fd1(2) = 0.0_dp
          fd1(3) = 1.0_dp - 4.0_dp * g3
          fd1(4) = 4.0_dp * g2
          fd1(5) = -4.0_dp * g2
          fd1(6) = 4.0_dp * (g3 - g1)

          fd2(1) = 0.0_dp
          fd2(2) = 4.0_dp * g2 - 1.0_dp
          fd2(3) = 1.0_dp - 4.0_dp * g3
          fd2(4) = 4.0_dp * g1
          fd2(5) = 4.0_dp * (g3 - g2)
          fd2(6) = -4.0_dp * g1
        end if

        do in_idx = 1, ncone
          do jn_idx = 1, ncone
            dc(in_idx, jn_idx)  = fd1(in_idx) * fd2(jn_idx)
            dct(in_idx, jn_idx) = fd2(in_idx) * fd1(jn_idx)
          end do
        end do

        aux = 0.0_dp; bux = 0.0_dp; cux = 0.0_dp
        do jk = 1, ncone
          do ik = 1, ncone
            aux(jk) = aux(jk) + xj2(ik) * (dc(ik, jk) - dct(ik, jk))
            bux(jk) = bux(jk) + xj3(ik) * (dc(ik, jk) - dct(ik, jk))
            cux(jk) = cux(jk) + xj1(ik) * (dc(ik, jk) - dct(ik, jk))
          end do
        end do

        s1 = sum(aux(1:ncone) * xj3(1:ncone))
        s2 = sum(bux(1:ncone) * xj1(1:ncone))
        s3 = sum(cux(1:ncone) * xj2(1:ncone))

        xja = sqrt(s1**2 + s2**2 + s3**2)
        eta(1) = s1 / xja
        eta(2) = s2 / xja
        eta(3) = s3 / xja

        x_pt = 0.0_dp
        do ik = 1, ncone
          x_pt(1) = x_pt(1) + f(ik) * xj1(ik)
          x_pt(2) = x_pt(2) + f(ik) * xj2(ik)
          x_pt(3) = x_pt(3) + f(ik) * xj3(ik)
        end do

        do ik = 1, ncone
          p12(ik) = f(ik) * xja * (gi_quad(i) + 1.0_dp) * area * 0.25_dp * ome_quad(i) * ome_quad(j)
        end do

        call funkelvin(u, t, xi(1), xi(2), xi(3), eta, x_pt(1), x_pt(2), x_pt(3))

        do kk = 1, ncone
          do ii = 1, 3
            do jj = 1, 3
              ind = 3 * kk + ii - 3
              g_elem(jj, ind) = g_elem(jj, ind) + u(jj, ii) * p12(kk)
              h_elem(jj, ind) = h_elem(jj, ind) + t(jj, ii) * p12(kk)
            end do
          end do
        end do

      end do
    end do
  end subroutine triloc8

  subroutine funkelvin(u, t, xp, yp, zp, eta, xg, yg, zg)
    implicit none
    real(dp), intent(in) :: xp, yp, zp, eta(3), xg, yg, zg
    real(dp), intent(out) :: u(3, 3), t(3, 3)

    real(dp) :: r1, r2, r3, r, rdn, rd(3)
    real(dp) :: c3_over_r, c5_over_r2
    integer :: ii, jj

    r1 = xg - xp
    r2 = yg - yp
    r3 = zg - zp
    r = sqrt(r1**2 + r2**2 + r3**2)

    rdn = (r1 * eta(1) + r2 * eta(2) + r3 * eta(3)) / r
    rd(1) = r1 / r
    rd(2) = r2 / r
    rd(3) = r3 / r

    c3_over_r  = c3 / r
    c5_over_r2 = c5 / (r**2)

    do ii = 1, 3
      do jj = 1, 3
        ! Kelvin displacement fundamental solution U_jj,ii
        u(jj, ii) = (c4 * delta(jj, ii) + rd(ii) * rd(jj)) * c3_over_r

        ! Kelvin traction fundamental solution T_jj,ii
        t(jj, ii) = c5_over_r2 * (rdn * (c6 * delta(jj, ii) + 3.0_dp * rd(ii) * rd(jj)) + &
                                  c6 * (eta(jj) * rd(ii) - eta(ii) * rd(jj)))
      end do
    end do
  end subroutine funkelvin

end module assembly_mod
