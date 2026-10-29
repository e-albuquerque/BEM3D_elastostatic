module io_mod
  use data_mod
  use assembly_mod, only: etaxja, extin8
  implicit none

contains

  subroutine read_bc_line(unit_num, ie, kode_val, vals)
    implicit none
    integer, intent(in) :: unit_num
    integer, intent(out) :: ie, kode_val
    real(dp), intent(out) :: vals(9)
    character(len=512) :: line
    real(dp) :: temp_pairs(2, 9)
    integer :: idx, p, n_pts

    vals = 0.0_dp
    read(unit_num, '(A)') line
    if (index(line, '(') > 0) then
      ! Handle legacy complex-formatted BC lines: (val_real, val_imag)
      do idx = 1, len_trim(line)
        if (line(idx:idx) == '(' .or. line(idx:idx) == ')') line(idx:idx) = ' '
      end do
      read(line, *) ie, kode_val
      n_pts = nconec(ie)
      read(line, *) ie, kode_val, (temp_pairs(1, p), temp_pairs(2, p), p = 1, n_pts)
      vals(1:n_pts) = temp_pairs(1, 1:n_pts)
    else
      ! Handle modern clean real BC lines: val1 val2 ...
      read(line, *) ie, kode_val
      n_pts = nconec(ie)
      read(line, *) ie, kode_val, (vals(p), p = 1, n_pts)
    end if
  end subroutine read_bc_line

  subroutine input_data()
    implicit none
    integer :: in_idx, ie, j, k, i, num_nodes, num_elems, num_pts, ios
    integer :: dummy_nfr, kode_tmp, n_pts_elem
    real(dp) :: in_x, in_y, in_z, dummy_dam, dummy_ro
    real(dp) :: bc_temp(9)
    character(len=512) :: line, dummy_freq_line

    open(unit=inp, file=trim(filein), status='old', action='read')
    open(unit=ipr, file=trim(fileout), status='replace', action='write')

    write(ipr, '(A)') ' ' // repeat('*', 79)
    read(inp, '(A)') title
    write(ipr, '(A)') trim(title)

    ! Read line 2: support either modern 5-value format (nodes, elems, pts, mu, nu)
    ! or legacy 8-value format (nodes, elems, pts, nfr, mu, dam, ro, nu)
    read(inp, '(A)') line
    read(line, *, iostat=ios) num_nodes, num_elems, num_pts, dummy_nfr, rge, dummy_dam, dummy_ro, xnu
    if (ios == 0) then
      ! Legacy format: consume frequencies line
      read(inp, '(A)') dummy_freq_line
    else
      ! Modern format
      read(line, *) num_nodes, num_elems, num_pts, rge, xnu
    end if

    call update_material_constants()
    call allocate_bem_data(num_nodes, num_elems, num_pts)

    write(ipr, 300) n, ne, l, rge, young, xnu
300 format(//2X, 'DATA'//2X, 'NUMBER OF BOUNDARY NODES    =', I5/2X, &
           'NUMBER OF BOUNDARY ELEMENTS =', I5/2X, &
           'NUMBER OF INTERNAL POINTS   =', I5/2X, &
           'SHEAR MODULUS (MU)          =', G14.7/2X, &
           'YOUNG MODULUS (E)           =', G14.7/2X, &
           'POISSON RATIO (NU)          =', G14.7)

    ! Read boundary node coordinates
    write(ipr, 500)
500 format(//2X, 'BOUNDARY NODES COORDINATES'///4X, 'NODE', 12X, 'X', 18X, 'Y', 18X, 'Z'/)

    do i = 1, n
      read(inp, *) in_idx, in_x, in_y, in_z
      x(in_idx) = in_x
      y(in_idx) = in_y
      z(in_idx) = in_z
      write(ipr, 240) in_idx, x(in_idx), y(in_idx), z(in_idx)
    end do
240 format(5X, I5, 3(5X, E14.7))

    ! Read element connectivities
    write(ipr, 600)
600 format(//2X, 'ELEMENT CONNECTIVITIES'///2X, 'ELEMENT', 1X, 'NO.CON.', 24X, 'NODES'/)
    do i = 1, ne
      read(inp, *) ie, nconec(ie), (kconec(j, ie), j = 1, nconec(ie))
      write(ipr, 340) ie, nconec(ie), (kconec(j, ie), j = 1, nconec(ie))
    end do
340 format(4X, I5, 4X, I2, 2X, 9(3X, I5))

    ! Read boundary conditions
    write(ipr, 800)
800 format(//2X, 'BOUNDARY CONDITIONS'//1X, 'ELEMEN', 1X, 'CODE', 22X, 'PRESCRIBED VALUES'/)
    do i = 1, ne
      do k = 1, 3
        call read_bc_line(inp, ie, kode_tmp, bc_temp)
        kode(3*ie - 3 + k) = kode_tmp
        do j = 1, nconec(ie)
          bc(3*j - 3 + k, ie) = bc_temp(j)
        end do
      end do
      write(ipr, 950) ie, kode(3*ie - 2), (bc(3*j - 2, ie), j = 1, nconec(ie))
      write(ipr, 960)     kode(3*ie - 1), (bc(3*j - 1, ie), j = 1, nconec(ie))
      write(ipr, 970)     kode(3*ie),     (bc(3*j,     ie), j = 1, nconec(ie))
    end do
950 format(1X, I5, 4X, I2, 1X, 'X', 9(1X, G12.5))
960 format(8X,     I2, 1X, 'Y', 9(1X, G12.5))
970 format(8X,     I2, 1X, 'Z', 9(1X, G12.5))

    ! Read internal points
    if (l > 0) then
      read(inp, *) (cx(i), cy(i), cz(i), i = 1, l)
      write(ipr, 1000)
1000  format(//2X, 'COORDINATES OF INTERNAL POINTS'//4X, 'POINT', 11X, 'X', 18X, 'Y', 17X, ' Z'/)
      do i = 1, l
        write(ipr, 240) i, cx(i), cy(i), cz(i)
      end do
    end if

    ! Determine elements sharing each node
    nesim = 0
    do ie = 1, ne
      do j = 1, nconec(ie)
        in_idx = kconec(j, ie)
        nesim(in_idx) = nesim(in_idx) + 1
        if (nesim(in_idx) <= 8) then
          ke(nesim(in_idx), in_idx) = ie
        end if
      end do
    end do

    ! Read tracked node from last line
    read(inp, *, end=190, err=190) ntrack
    goto 195
190 ntrack = 1
195 continue
    write(ipr, '(A, I6)') ' TRACKED NODE = ', ntrack
    write(*,   '(A, I6)') ' Tracking Node = ', ntrack

    close(inp)
  end subroutine input_data

  subroutine reorganize_solution()
    implicit none
    integer :: ie, ld, k, ik, in_idx, j, m

    ! Reorganize FI and DFI arrays: FI holds displacements, DFI holds tractions
    do ie = 1, ne
      do ld = 1, 3
        if (kode(3*ie - 3 + ld) == 0) then
          do k = 1, nconec(ie)
            ik = kconec(k, ie)
            DFI(3*k - 3 + ld, ie) = FI(3*ik - 3 + ld) * rge
          end do
        else
          do k = 1, nconec(ie)
            DFI(3*k - 3 + ld, ie) = bc(3*k - 3 + ld, ie)
          end do
        end if
      end do
    end do

    do in_idx = 1, n
      do ld = 1, 3
        do j = 1, nesim(in_idx)
          if (kode(3*ke(j, in_idx) - 3 + ld) == 0) then
            do m = 1, nconec(ke(j, in_idx))
              if (in_idx == kconec(m, ke(j, in_idx))) then
                FI(3*in_idx - 3 + ld) = bc(3*m - 3 + ld, ke(j, in_idx))
                exit
              end if
            end do
            exit
          end if
        end do
      end do
    end do
  end subroutine reorganize_solution

  subroutine inter8()
    implicit none
    ! Backward-compatible wrapper
    call reorganize_solution()
  end subroutine inter8

  subroutine output_data()
    implicit none
    integer :: i, ie, k, ik

    write(ipr, '(A)') ' ' // repeat('*', 79)
    write(ipr, '(A)') ' QUADESTH: 3D ELASTOSTATIC ANALYSIS RESULTS'
    write(ipr, '(A)') ' ' // repeat('*', 79)
    write(ipr, '(A)') ' BOUNDARY NODES DISPLACEMENTS'
    write(ipr, '(A)') ' NODE           DISPLACEMENT X          DISPLACEMENT Y          DISPLACEMENT Z'

    do i = 1, n
      write(ipr, 200) i, FI(3*i-2), FI(3*i-1), FI(3*i)
    end do
200 format(1X, I5, 3(5X, E18.8))

    write(ipr, 110)
110 format(//2X, 'BOUNDARY ELEMENTS TRACTIONS'//1X, 'ELEMENT', 2X, 'NODE', 8X, &
           'TRACTION X', 14X, 'TRACTION Y', 14X, 'TRACTION Z'/)

    do ie = 1, ne
      write(ipr, '(1X, I5)') ie
      do k = 1, nconec(ie)
        ik = kconec(k, ie)
        write(ipr, 130) ik, DFI(3*k-2, ie), DFI(3*k-1, ie), DFI(3*k, ie)
      end do
    end do
130 format(10X, I5, 3(5X, E18.8))

    if (l > 0) then
      write(ipr, 310)
310   format(//2X, 'INTERNAL POINTS DISPLACEMENTS'//'POINT', 7X, &
             'DISPLACEMENT X', 10X, 'DISPLACEMENT Y', 10X, 'DISPLACEMENT Z'/)
      do k = 1, l
        write(ipr, 200) k, dsol(3*k-2), dsol(3*k-1), dsol(3*k)
      end do
    end if
    write(ipr, '(A)') ' ' // repeat('*', 79)
  end subroutine output_data

  subroutine save_disp_trac()
    implicit none
    integer :: ie, kloc, k
    real(dp) :: ux, uy, uz, tx, ty, tz

    ie = ke(1, ntrack)
    kloc = 1
    do k = 1, nconec(ie)
      if (kconec(k, ie) == ntrack) then
        kloc = k
        exit
      end if
    end do

    ux = FI(3*ntrack - 2)
    uy = FI(3*ntrack - 1)
    uz = FI(3*ntrack)

    tx = DFI(3*kloc - 2, ie)
    ty = DFI(3*kloc - 1, ie)
    tz = DFI(3*kloc,     ie)

    write(iout_track, 600) ux, uy, uz, tx, ty, tz
600 format(6(2X, E18.8))
    flush(iout_track)
  end subroutine save_disp_trac

end module io_mod
