program quadesth_modern
  use data_mod
  use solver_mod
  use assembly_mod
  use stress_mod
  use io_mod
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  real(dp) :: t_start, t_end, t_assembly_start, t_assembly_end, t_solve_start, t_solve_end
  real(dp) :: t_stress_start, t_stress_end, t_int_start, t_int_end
  integer :: num_threads

  !$ integer, external :: omp_get_max_threads

  write(*, '(A)') '============================================================='
  write(*, '(A)') ' QUADESTH Modernized BEM 3D Elastostatics Solver'
  write(*, '(A)') ' Quadratic Elements - Kelvin Fundamental Solution'
  write(*, '(A)') ' OpenMP Parallelization & LAPACK Real Solver (DGESV)'
  write(*, '(A)') ' Boundary Stresses (Traction-Recovery Method)'
  write(*, '(A)') '============================================================='

  num_threads = 1
  !$ num_threads = omp_get_max_threads()
  write(*, '(A, I4)') ' OpenMP threads available: ', num_threads

  call cpu_time(t_start)

  ! 1. Read input data
  call input_data()
  write(*, '(A, I6, A, I6, A, I4, A)') ' Mesh loaded: ', n, ' nodes, ', ne, ' elements, ', l, ' internal points.'

  ! 2. Open tracking file disp_trac.dat
  open(unit=iout_track, file=trim(filetrack), status='replace', action='write')

  ! 3. Precompute element Gauss geometry once for all elements
  write(*, '(A)') ' Precomputing element Gauss point geometries...'
  call precompute_element_geometry()

  ! 4. Assemble BEM system matrix A and RHS FI in parallel
  write(*, '(A)', advance='no') ' Assembling influence matrices (G and H)... '
  call cpu_time(t_assembly_start)
  call ghmat8()
  call cpu_time(t_assembly_end)
  write(*, '(A, F6.2, A)') 'done in ', (t_assembly_end - t_assembly_start), ' s'

  ! 5. Solve linear system A * X = FI using LAPACK DGESV
  write(*, '(A)', advance='no') ' Solving linear system with LAPACK DGESV... '
  call cpu_time(t_solve_start)
  call solve_linear_system(A, FI, 3*n)
  call cpu_time(t_solve_end)
  write(*, '(A, F6.2, A)') 'done in ', (t_solve_end - t_solve_start), ' s'

  ! 6. Reorganize solution vector (displacements in FI, tractions in DFI)
  call reorganize_solution()

  ! 7. Compute boundary stresses (traction-recovery method)
  write(*, '(A)', advance='no') ' Computing boundary stresses (Traction-Recovery)... '
  call cpu_time(t_stress_start)
  call compute_boundary_stresses()
  call cpu_time(t_stress_end)
  write(*, '(A, F6.2, A)') 'done in ', (t_stress_end - t_stress_start), ' s'

  ! 8. Compute internal points displacements and stresses (Somigliana identities)
  if (l > 0) then
    write(*, '(A, I6, A)', advance='no') ' Computing internal displacements and stresses for ', l, ' points... '
    call cpu_time(t_int_start)
    call compute_internal_points()
    call cpu_time(t_int_end)
    write(*, '(A, F6.2, A)') 'done in ', (t_int_end - t_int_start), ' s'
  end if

  ! 9. Write full results to output files
  call output_data()
  call output_boundary_stresses()
  if (l > 0) call output_internal_points()

  ! 10. Write tracked node displacements and tractions
  call save_disp_trac()

  ! Close files
  close(ipr)
  close(iout_track)

  ! Clean up memory
  call deallocate_stress_data()
  call deallocate_bem_data()

  call cpu_time(t_end)
  write(*, '(A)') '============================================================='
  write(*, '(A, F8.2, A)') ' Computation completed successfully in ', (t_end - t_start), ' s'
  write(*, '(A)') ' Results exported to: ' // trim(fileout) // ', ' // trim(filestress) // ', ' // trim(fileintpts) // ' and ' // trim(filetrack)
  write(*, '(A)') '============================================================='

end program quadesth_modern
