# Modernized BEM for Three-Dimensional Elastostatics using Quadratic Elements (QUADESTH)

[![Fortran](https://img.shields.io/badge/Language-Fortran%2090-734f96?logo=fortran)](https://fortran-lang.org/)
[![Python](https://img.shields.io/badge/Language-Python%203-3776ab?logo=python)](https://python.org/)
[![OpenMP](https://img.shields.io/badge/Parallel-OpenMP-blue)](https://www.openmp.org/)
[![LAPACK](https://img.shields.io/badge/Solver-LAPACK%20DGESV-orange)](https://netlib.org/lapack/)
[![License](https://img.shields.io/badge/License-MIT%20%2F%20Academic-green)]()

This repository contains a **Fortran 90 / OpenMP / LAPACK** implementation of the **QUADESTH** program: *A Boundary Element code using (QUAD)ratic elements for the solution of (E)elasto-(S)tatic (TH)ree-dimensional problems*. 

It is an evolution, modernization, and specialization of the original **QUADTHEH** dynamic code presented by José Dominguez in his textbook:

> **Dominguez, J. (1993)**. *Boundary Elements in Dynamics*. Computational Mechanics Publications, Southampton & Boston (Section 3.5).

The formulation solves 3D linear elastostatics problems using 9-node Lagrangian quadrilateral and 6-node triangular quadratic boundary elements, based on the **3D Kelvin fundamental solutions**.

---

## 🚀 Main Features

The monolithic Fortran 77 frequency-domain program (`quadtheh.for`) was thoroughly re-engineered into modern Fortran 90 standards and tailored specifically for **pure elastostatics**, eliminating all complex arithmetic overhead:

### 1. 3D Kelvin Fundamental Solutions
Employs exact real-valued displacement ($U_{ij}^\ast$) and surface traction ($T_{ij}^\ast$) fundamental solutions for the Kelvin problem in 3D:

$$
U_{ij}^\ast(x, \xi) = \frac{1}{16\pi (1-\nu)\mu r} \left[ (3 - 4\nu)\delta_{ij} + r_{,i} r_{,j} \right]
$$

$$
T_{ij}^\ast(x, \xi) = \frac{-1}{8\pi (1-\nu) r^2} \left\{ \frac{\partial r}{\partial n}\left[ (1 - 2\nu)\delta_{ij} + 3 r_{,i} r_{,j} \right] + (1 - 2\nu)\left( n_j r_{,i} - n_i r_{,j} \right) \right\}
$$

### 2. Exact Rigid Body Translation Condition
The singular diagonal blocks of the influence matrix $H$ (containing the free-term tensor $C_{ij}(\xi)$ and the Cauchy Principal Value of $T_{ij}^\ast$) are determined directly and rigorously by applying rigid body translations ($\{u\} = \{1\}$, $\{t\} = \{0\}$):

$$
\sum_{j} H_{ij} = 0 \implies H_{ii} = - \sum_{j \neq i} H_{ij}
$$

This eliminates the need for dynamic frequency sweeps, complex exponentials, and Bessel series expansions (`tdif` / `hdif`).

### 3. Pure Real Double Precision (64-bit)
All arrays and system variables are strictly `real(dp)` (`real64`). This cuts memory consumption in half and avoids all complex number operations.

### 4. Optimized Linear System Solver (LAPACK `DGESV`)
Solves the fully populated system matrix using LAPACK's vectorized real double precision solver `DGESV` with partial pivoting.

### 5. OpenMP Parallel Collocation Assembly
The outer collocation loop in `ghmat8` is parallelized across CPU threads using OpenMP (`!$omp parallel do schedule(dynamic)`). Precomputed Gauss quadrature geometries allow assembling the entire 3D mesh in ~1 second.

### 6. Modern Modular Architecture
Decomposed into clean, maintainable Fortran modules (`data_mod.f90`, `solver_mod.f90`, `assembly_mod.f90`, `io_mod.f90`, `main.f90`) enforcing `IMPLICIT NONE` and dynamic memory allocation (`ALLOCATABLE`).

### 7. Automated Gmsh & Python Pipeline
- Seamless conversion from Gmsh quadratic quadrilateral meshes (`quad9`) via `meshio`.
- Node-ordering permutation algorithm aligning Gmsh sequencing `[0, 1, 2, 3, 4, 5, 6, 7, 8]` with Dominguez's boundary integration convention `[0, 4, 1, 5, 2, 6, 3, 7, 8]`.
- 3D physical surface visualization (`mesh_faces.png`) and automated analytical validation (`result_plot.png`).

### 8. Accurate Boundary Stress Recovery (Traction-Recovery Method)
Computes the full 3D Cauchy stress tensor ($\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{yz}, \sigma_{zx}$), principal stresses ($\sigma_1 \ge \sigma_2 \ge \sigma_3$), and Von Mises equivalent stress directly on the boundary:
- Uses a local orthonormal coordinate frame $(\vec{e}_1, \vec{e}_2, \vec{e}_3)$ on each element surface, where $\vec{e}_3 = \vec{n}$ is the unit normal.
- Out-of-plane stresses are obtained directly from surface tractions: $\hat{\sigma}_{13} = \hat{t}_1, \hat{\sigma}_{23} = \hat{t}_2, \hat{\sigma}_{33} = \hat{t}_3$.
- In-plane strains ($\hat{\varepsilon}_{11}, \hat{\varepsilon}_{22}, \hat{\varepsilon}_{12}$) are derived from surface displacement gradients using analytical derivatives of **9-node Lagrangian quadrilateral** and **6-node quadratic triangular** elements (adapted from Gao & Davies, 2001).
- In-plane stresses are recovered via isotropic Hooke's Law:

$$
\hat{\sigma}_{11} = \frac{2\mu}{1-\nu}\hat{\varepsilon}_{11} + \frac{2\mu\nu}{1-\nu}\hat{\varepsilon}_{22} + \frac{\nu}{1-\nu}\hat{\sigma}_{33}
$$

$$
\hat{\sigma}_{22} = \frac{2\mu}{1-\nu}\hat{\varepsilon}_{22} + \frac{2\mu\nu}{1-\nu}\hat{\varepsilon}_{11} + \frac{\nu}{1-\nu}\hat{\sigma}_{33}
$$

$$
\hat{\sigma}_{12} = 2\mu\hat{\varepsilon}_{12}
$$

- Stresses are transformed back to global coordinates ($\boldsymbol{\sigma} = \mathbf{L}^T \hat{\boldsymbol{\sigma}} \mathbf{L}$) and averaged across adjacent elements at boundary nodes.

---

## 📂 Repository Structure

```text
├── src / root
│   ├── main.f90              # Main driver (single-pass static solution)
│   ├── data_mod.f90          # Global data structures, dynamic arrays, and material constants
│   ├── assembly_mod.f90      # Parallel assembly (G & H), Gauss quadrature, and Kelvin solutions
│   ├── stress_mod.f90        # Boundary stress recovery (Traction-Recovery Method, 9-node Lagrange)
│   ├── solver_mod.f90        # Interface to LAPACK DGESV real linear solver
│   ├── io_mod.f90            # Flexible file readers, Somigliana internal point solver, and output
│   ├── quadtheh.for          # Original F77 reference code (frequency domain)
│   ├── Makefile              # GCC / GFortran build file
│   └── build.bat             # Windows compilation batch script
├── mesh / preproc
│   ├── placa.geo             # Gmsh script defining geometry, physical entities, and mesh
│   ├── placa.msh             # Quadratic 3D boundary element mesh (quad9)
│   └── preproc.py            # Preprocessor: Gmsh reader, real BC formatting, and 3D visualizer
├── postproc / validation
│   ├── posproc.py            # Postprocessor: 1D analytical validation and error metrics
│   ├── disp_trac.dat         # Tracked node displacement and traction output
│   ├── stresses.dat          # Complete boundary nodal stresses, Von Mises, and principal stresses
│   ├── cube8.dat             # Simulation input file generated from mesh
│   ├── cube8_out.dat         # Detailed node displacements, element tractions, and internal points
│   ├── mesh_faces.png        # 3D boundary mesh faces visualization
│   └── result_plot.png       # Numerical vs Analytical validation plot
└── README.md
```

---

## ⚙️ Compilation & Prerequisites

### Prerequisites
- **GFortran** (GCC 10 or later recommended)
- **LAPACK & BLAS** libraries
- **OpenMP** runtime library (`libgomp`)
- **Python 3.8+** with `numpy`, `matplotlib`, and `meshio`

#### Installing Prerequisites on Windows (MSYS2 UCRT64 / MinGW):
```bash
pacman -S mingw-w64-ucrt-x86_64-gcc-fortran \
          mingw-w64-ucrt-x86_64-lapack \
          mingw-w64-ucrt-x86_64-make
```

#### Installing Prerequisites on Ubuntu / Debian:
```bash
sudo apt-get update
sudo apt-get install gfortran liblapack-dev libblas-dev make python3-pip
pip install numpy matplotlib meshio
```

### Compiling the Solver

#### Using Make:
```bash
make
```
*(On Windows MSYS2, use `mingw32-make` or `make`)*

#### Using the Windows Batch Script:
```cmd
build.bat
```

This compiles all modules with `-O3 -fopenmp -fdefault-real-8 -static -llapack -lblas` and generates `quadesth_modern.exe`.

---

## 🔬 Benchmark Problem: Uniaxial Strain on a 3D Elastic Cube

The implementation is verified against the 1D uniaxial strain benchmark (oedometric condition):

* **Geometry**: $6\text{ m} \times 6\text{ m} \times 6\text{ m}$ elastic cube ($L = 6\text{ m}$).
* **Material Properties**:
  * Shear modulus: $\mu = 1.0 \times 10^6\text{ N/m}^2$
  * Poisson's ratio: $\nu = 0.25$
  * Young's modulus: $E = 2\mu(1+\nu) = 2.5 \times 10^6\text{ N/m}^2$
  * Constrained (oedometric) modulus: $M = \frac{E(1-\nu)}{(1+\nu)(1-2\nu)} = 3.0 \times 10^6\text{ N/m}^2$
* **Boundary Conditions**:
  * **Upper Face ($z = 6\text{ m}$, `sup1`)**: Uniform normal traction $p = 100\text{ N/m}^2$; zero shear tractions.
  * **Bottom Face ($z = 0\text{ m}$, `sup3`)**: Constrained normal displacement ($u_z = 0$); zero shear tractions.
  * **Lateral Faces ($x = 0, 6\text{ m}$ and $y = 0, 6\text{ m}$)**: Roller supports with constrained normal displacements ($u_x = 0$ or $u_y = 0$); zero shear tractions.
* **Discretization**:
  * 216 quadratic quadrilateral boundary elements (9-node Lagrangian quadrilaterals, 36 elements per face).
  * 866 boundary collocation nodes (2598 degrees of freedom).
  * 1 internal evaluation point at the domain center $(3, 3, 3)$.

---

## 🏃 Workflow & Execution

### 1. Preprocessing & Mesh Visualization
Generate the clean elastostatic input file and visualize the mesh:
```bash
python preproc.py
```
This reads `placa.msh`, generates `cube8.dat` (and `quad8.dat`), detects the top-center monitor node at $(3, 3, 6)$ (Node 770), and exports `mesh_faces.png`:

![3D Boundary Mesh Faces](mesh_faces.png)

### 2. Numerical Simulation
Execute the parallel Fortran 90 elastostatics solver:
```bash
./quadesth_modern.exe
```

Console output:
```text
=============================================================
 QUADESTH Modernized BEM 3D Elastostatics Solver
 Quadratic Elements - Kelvin Fundamental Solution
 OpenMP Parallelization & LAPACK Real Solver (DGESV)
=============================================================
 OpenMP threads available:    8
 Tracking Node =    770
 Mesh loaded:    866 nodes,    216 elements,    1 internal points.
 Precomputing element Gauss point geometries...
 Assembling influence matrices (G and H)... done in   1.02 s
 Solving linear system with LAPACK DGESV... done in   1.94 s
=============================================================
 Computation completed successfully in     2.98 s
 Results exported to: cube8_out.dat and disp_trac.dat
=============================================================
```

### 3. Post-Processing & Validation
Validate the results against the analytical solution ($u_z(z) = \frac{p \cdot z}{M}$):
```bash
python posproc.py
```

---

## 📈 Validation Results

Under uniaxial strain conditions, the analytical solution yields a linear vertical displacement profile:

$$
u_z(z) = \frac{p \cdot z}{M}
$$

- At the upper face ($z = 6\text{ m}$): $u_z^{\text{analyt}} = \frac{100 \times 6}{3.0 \times 10^6} = 2.00000 \times 10^{-4}\text{ m}$
- At the domain center ($z = 3\text{ m}$): $u_z^{\text{analyt}} = \frac{100 \times 3}{3.0 \times 10^6} = 1.00000 \times 10^{-4}\text{ m}$

### Comparison Summary:

#### 1. Displacements and Applied Tractions:
| Evaluation Point | Coordinate $z$ | Numerical (QUADESTH) | Analytical | Relative Error |
| :--- | :---: | :---: | :---: | :---: |
| **Top-center Node 770 ($u_z$)** | $6.0\text{ m}$ | $2.00044 \times 10^{-4}\text{ m}$ | $2.00000 \times 10^{-4}\text{ m}$ | **0.022 %** |
| **Internal Point 1 ($u_z$)** | $3.0\text{ m}$ | $1.00057 \times 10^{-4}\text{ m}$ | $1.00000 \times 10^{-4}\text{ m}$ | **0.057 %** |
| **Top Traction $t_z$** | $6.0\text{ m}$ | $100.00000\text{ N/m}^2$ | $100.00000\text{ N/m}^2$ | **0.000 %** |

#### 2. Recovered Boundary Stresses (Traction-Recovery Method):
| Stress Component | Node / Region | Numerical (QUADESTH) | Analytical | Relative Error |
| :--- | :---: | :---: | :---: | :---: |
| **Vertical Stress $\sigma_{zz}$** | Node 770 ($z=6\text{m}$) | $100.0000\text{ Pa}$ | $100.0000\text{ Pa}$ | **0.0000 %** |
| **Lateral Stress $\sigma_{xx}$** | Node 770 ($z=6\text{m}$) | $33.3490\text{ Pa}$ | $33.3333\text{ Pa}$ | **0.0469 %** |
| **Lateral Stress $\sigma_{yy}$** | Node 770 ($z=6\text{m}$) | $33.3382\text{ Pa}$ | $33.3333\text{ Pa}$ | **0.0147 %** |
| **Von Mises Stress $\sigma_{vm}$** | Node 770 ($z=6\text{m}$) | $66.6564\text{ Pa}$ | $66.6667\text{ Pa}$ | **0.0154 %** |
| **Mean $\sigma_{zz}$ (Mesh Average)** | All 866 nodes | $100.0363\text{ Pa}$ | $100.0000\text{ Pa}$ | **0.0363 %** |
| **Mean $\sigma_{xx}$ (Mesh Average)** | All 866 nodes | $33.3321\text{ Pa}$ | $33.3333\text{ Pa}$ | **0.0036 %** |
| **Mean $\sigma_{yy}$ (Mesh Average)** | All 866 nodes | $33.3369\text{ Pa}$ | $33.3333\text{ Pa}$ | **0.0108 %** |
| **Mean $\sigma_{vm}$ (Mesh Average)** | All 866 nodes | $66.7231\text{ Pa}$ | $66.6667\text{ Pa}$ | **0.0846 %** |

![Numerical vs Analytical Validation Plot](result_plot.png)

---

## 📖 References

* **Dominguez, J.** (1993). *Boundary Elements in Dynamics*. Computational Mechanics Publications, Southampton & Boston. ISBN: 1-85312-258-0.
* **Gao, X.-W., & Davies, T. G.** (2001). *Boundary Element Programming in Mechanics*. Cambridge University Press. ISBN: 978-0521772594.
* **Brebbia, C. A., & Dominguez, J.** (1992). *Boundary Elements: An Introductory Course*. Computational Mechanics Publications & McGraw-Hill.
