# ============================================================
# nlfem configuration — edit paths for your system
# ============================================================

# MPI compiler wrappers
CXX = mpicxx
CC  = mpicc
FC  = mpif90

# PETSc: set PETSC_DIR and PETSC_ARCH
# PETSC_ARCH is the build subdirectory (e.g. arch-linux2-c-opt).
# Run: ls $PETSC_DIR  to see which arch directories exist.
PETSC_DIR  = $(HOME)/dirs/research/bitcart/dependencies/petsc-3.6.4/petsc-3.6.4/install_opt
#PETSC_ARCH = arch-linux2-c-opt

# ParMETIS
PARMETIS_DIR = $(HOME)/dirs/research/bitcart/dependencies/parmetis-4.0.3/install

METIS_DIR    = $(HOME)/dirs/research/bitcart/dependencies/parmetis-4.0.3/metis/install

# LAPACK (set to empty if LAPACK is bundled inside PETSc)
LAPACK_DIR =

# Set to 1 for debug build
DEBUG = 0

# ------------------------------------------------------------
# GPU (CUDA) element assembly -- Phase 1: CABLE elements only.
# Set USEGPU=0 to build the original all-CPU codebase unchanged
# (src/gpu/gpu_assembly_stub.cc is used instead of the real kernel,
# so nothing about the CPU path or its output changes).
# ------------------------------------------------------------
USEGPU   = 1
NVCC     = nvcc
CUDA_DIR = /usr/local/cuda
# Quadro K1200 = Maxwell, compute capability 5.0 (sm_50).
# Check yours with: nvidia-smi --query-gpu=compute_cap --format=csv
CUDA_ARCH = sm_50
