# nlfem — Nonlinear Finite Element Solver

**Author:** Gokul G. Anugrah  
**Contact:** gokulanugrah@gmail.com

---

## Overview

`nlfem` is a parallel MPI nonlinear structural finite element solver for large-deformation analysis of flexible structures, including parachute canopies, suspension cables, and shell structures.

### Capabilities

- **Elements:** MITC3 triangular shell, Bernoulli-Euler beam, cable elements
- **Kinematics:** Geometric and material nonlinearity (large deformation, large rotation via director vectors)
- **Linear algebra:** PETSc parallel KSP solvers (GAMG, MUMPS, etc.)
- **Domain decomposition:** Parallel partitioning via ParMETIS
- **Time integration:** Implicit Newmark-Beta (nonlinear and linear unsteady); explicit Runge-Kutta (linear)
- **FSI coupling:** Load and displacement transfer to/from an auxiliary fluid surface mesh via WLSQR interpolation
- **Restart:** Binary checkpoint/restart
- **Output:** Nodal time histories; VTK snapshots of deformed geometry and stress fields

### Validated for

- MITC3 shell elements under large deformation (parachute inflation)
- Cable element networks (suspension line dynamics)
- Fluid-Structure Interaction with external CFD solvers

---

## Dependencies

| Dependency | Notes |
|---|---|
| MPI (OpenMPI / MPICH) | `mpicxx`, `mpif90` compiler wrappers |
| PETSc | Required. Tested with PETSc 3.6.4 |
| ParMETIS | Required for parallel domain decomposition |
| LAPACK / BLAS | Required (linked through PETSc) |
| gfortran | Required for Fortran 90 MITC3 tying routines |

---

## Build

### 1. Edit `config.mk`

Set paths to match your installation:

```makefile
CXX          = mpicxx
CC           = mpicc
FC           = mpif90

PETSC_DIR    = /path/to/petsc-3.6.4
PARMETIS_DIR = /path/to/parmetis-4.0.3
LAPACK_DIR   = /path/to/lapack-3.5.0
```

On Gokul's cluster this is typically:
```makefile
PETSC_DIR    = $(HOME)/dirs/research/bitcart/dependencies/petsc-3.6.4
PARMETIS_DIR = $(HOME)/dirs/research/bitcart/dependencies/parmetis-4.0.3
LAPACK_DIR   = $(HOME)/dirs/research/bitcart/dependencies/lapack-3.5.0
```

### 2. Build

```bash
make
```

This produces the `./nlfem` executable.

For a debug build:
```bash
make DEBUG=1
```

---

## Running

```bash
mpirun -np 10 ./nlfem input.dat
```

Required input files in the working directory:
- `input.dat` — solver configuration (see `example/input.dat`)
- `fem.vtk` — structural mesh (VTK unstructured grid)
- `bc.txt` — fixed node boundary conditions

For FSI runs, also provide:
- `geo.vtk` — auxiliary aerodynamic surface mesh

---

## Input File Format

Nested bracket syntax. See the fully annotated `example/input.dat`.

```
[Structural]
  [[ElementProperties]]   # Mesh file, material, element type
  [[Solver]]              # Steady/unsteady, linear/nonlinear, FSI toggle
  [[Loading]]             # Body forces and pressure loading
  [[FSI]]                 # Auxiliary geometry and interpolation settings
  [[Output]]              # Node time-history recording
  [[Unsteady]]            # dt, Newmark-Beta constants, number of steps
  [[Restart]]             # Checkpoint interval and restart timestep
  [[Nonlinear]]           # Load steps and Newton iteration limit
[/Structural]
```

---

## Output directories (created automatically)

| Path | Contents |
|---|---|
| `post/fem/` | VTK snapshots of deformed FEM mesh |
| `post/aux/` | VTK snapshots of auxiliary (fluid) surface |
| `post/nodeData/` | Nodal time-history ASCII files |
| `post/forces/` | Transferred force distributions |
| `restart_fem/` | Binary restart checkpoint files |

---

## Pushing to GitHub

```bash
git init
git add .
git commit -m "Initial commit: nlfem nonlinear FEM solver"
git branch -M main
git remote add origin https://github.com/magadator/nlfem.git
git push -u origin main
```

---

## License

Developed as part of doctoral research in computational structural mechanics.  
All rights reserved by the author. Contact gokulanugrah@gmail.com for licensing or collaboration.
