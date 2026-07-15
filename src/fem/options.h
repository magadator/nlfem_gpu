/*
 * nlfem - Nonlinear FEM Solver
 * Author  : Gokul G. Anugrah
 * Contact : gokulanugrah@gmail.com
 *
 * options.h — Compile-time solver configuration.
 *
 * USEPETSC=1       : Use PETSc for parallel linear algebra (required for full solver)
 * DEDICATED_FEM_PROC=0 : All MPI ranks participate in FEM domain decomposition
 *
 * These are normally set via -D flags in the Makefile (from config.mk).
 * The #ifndef guards here provide safe defaults if not set on command line.
 */

#ifndef USEPETSC
#define USEPETSC 1
#endif

#ifndef DEDICATED_FEM_PROC
#define DEDICATED_FEM_PROC 0
#endif
