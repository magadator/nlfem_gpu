# ============================================================
# nlfem — Nonlinear FEM Solver
# Author: Gokul G. Anugrah (gokulanugrah@gmail.com)
#
# Edit config.mk to set your PETSc, ParMETIS, and MPI paths.
# Then: make
# Run:  mpirun -np <N> ./nlfem example/input.dat
# ============================================================

include config.mk

BUILD_DIR = build

# -------------------------------------------------------
# Compiler flags
# -------------------------------------------------------
ifeq ($(DEBUG),1)
    OPT_FLAGS = -O0 -g -DDEBUG
else
    OPT_FLAGS = -O3 -g -DOPT
endif

CXXFLAGS = $(OPT_FLAGS) -std=c++11 -fPIC -DKARMA_MPI -DUSEPETSC=1 -DDEDICATED_FEM_PROC=0 -DUSEGPU=$(USEGPU)
FFLAGS   = $(OPT_FLAGS) -fPIC -DKARMA_MPI
NVCCFLAGS = -O3 -std=c++11 -arch=$(CUDA_ARCH) -Xcompiler -fPIC -Isrc/gpu

# -------------------------------------------------------
# PETSc paths (supports both arch-based and flat installs)
# -------------------------------------------------------
# PETSc 3.x installs headers and libs under $PETSC_DIR/$PETSC_ARCH/
# If PETSC_ARCH is empty, fall back to $PETSC_DIR directly.
ifneq ($(PETSC_ARCH),)
    PETSC_INC = $(PETSC_DIR)/include $(PETSC_DIR)/$(PETSC_ARCH)/include
    PETSC_LIB_DIR = $(PETSC_DIR)/$(PETSC_ARCH)/lib
else
    PETSC_INC = $(PETSC_DIR)/include
    PETSC_LIB_DIR = $(PETSC_DIR)/lib
endif

# -------------------------------------------------------
# Include paths
# -------------------------------------------------------
INCLUDES = \
    $(addprefix -I, $(PETSC_INC)) \
    -I$(PARMETIS_DIR)/include \
    -Isrc/fem \
    -Isrc/geometry \
    -Isrc/linear_algebra \
    -Isrc/mpi_ops \
    -Isrc/options_parser \
    -Isrc/screen \
    -Isrc/string_utilities \
    -Isrc/file_io \
    -Isrc/gpu

ifeq ($(USEGPU),1)
    INCLUDES += -I$(CUDA_DIR)/include
endif

# -------------------------------------------------------
# Link flags
# -------------------------------------------------------
# PETSc 3.6 bundles BLAS/LAPACK internally; link via petsc's own lib dir.
# Add -Wl,-rpath so the executable finds shared libs at runtime.
LDFLAGS = \
    -L$(PETSC_LIB_DIR) -Wl,-rpath,$(PETSC_LIB_DIR) -lpetsc \
    -L$(PARMETIS_DIR)/lib -L$(METIS_DIR)/lib -lparmetis -lmetis \
    -lgfortran -lm -lstdc++

# Add LAPACK only if LAPACK_DIR is set
ifneq ($(LAPACK_DIR),)
    LDFLAGS += -L$(LAPACK_DIR) -llapack -lblas
endif

ifeq ($(USEGPU),1)
    LDFLAGS += -L$(CUDA_DIR)/lib64 -Wl,-rpath,$(CUDA_DIR)/lib64 -lcudart
endif

# -------------------------------------------------------
# Source files
# -------------------------------------------------------
CXX_SRCS = \
    src/main.cc \
    src/fem/fem.cc \
    src/fem/linearTriElement.cc \
    src/geometry/surfaceTriangulation.cc \
    src/linear_algebra/linear_algebra.cc \
    src/mpi_ops/mpi_ops.cc \
    src/options_parser/options_parser.cc \
    src/options_parser/BuiltInFunctionList.cpp \
    src/options_parser/BuiltInFunctions.cpp \
    src/options_parser/BuiltInUtils.cpp \
    src/options_parser/PreprocessorContext.cpp \
    src/screen/screen.cc \
    src/string_utilities/string_utilities.cc \
    src/file_io/file_io.cc

ifeq ($(USEGPU),1)
    CU_SRCS  = src/gpu/cable_kernels.cu
else
    CXX_SRCS += src/gpu/gpu_assembly_stub.cc
endif

F90_SRCS = \
    src/fem/SM3MB.F90 \
    src/fem/SM3MH.F90 \
    src/fem/SM4M2.F90 \
    src/fem/LUFACT.F90 \
    src/fem/LUSOLV.F90

CXX_OBJS = $(patsubst src/%.cc, $(BUILD_DIR)/%.o, \
            $(patsubst src/%.cpp,$(BUILD_DIR)/%.o, $(CXX_SRCS)))
F90_OBJS = $(patsubst src/%.F90,$(BUILD_DIR)/%.o, $(F90_SRCS))
CU_OBJS  = $(patsubst src/%.cu, $(BUILD_DIR)/%.o, $(CU_SRCS))
OBJS     = $(CXX_OBJS) $(F90_OBJS) $(CU_OBJS)

# -------------------------------------------------------
# Targets
# -------------------------------------------------------
.PHONY: all clean run info

all: nlfem

nlfem: $(OBJS)
	@echo "Linking nlfem..."
	$(CXX) -o $@ $^ $(LDFLAGS)
	@echo ""
	@echo "  Build complete: ./nlfem"
	@echo "  Run: mpirun -np <N> ./nlfem example/input.dat"

# Compile .cc and .cpp
$(BUILD_DIR)/%.o: src/%.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compile Fortran 90
$(BUILD_DIR)/%.o: src/%.F90
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

# Compile CUDA .cu files (only reached when USEGPU=1, i.e. CU_SRCS is non-empty)
$(BUILD_DIR)/%.o: src/%.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# -------------------------------------------------------
# Standalone GPU validation harness (see src/gpu/validate_cable.cc).
# Deliberately does NOT depend on PETSc/MPI/ParMETIS -- just the kernel
# and a plain-CPU reference transcription of the same math, so you can
# sanity-check the kernel before trusting it inside a real nlfem run.
# -------------------------------------------------------
validate_cable: src/gpu/validate_cable.cc src/gpu/cable_kernels.cu
	@mkdir -p $(BUILD_DIR)/gpu
	$(NVCC) $(NVCCFLAGS) -c src/gpu/cable_kernels.cu -o $(BUILD_DIR)/gpu/cable_kernels.o
	$(NVCC) $(NVCCFLAGS) src/gpu/validate_cable.cc $(BUILD_DIR)/gpu/cable_kernels.o -o validate_cable
	@echo "Run ./validate_cable -- expect PASS with diff ~1e-12 or smaller."

clean:
	rm -rf $(BUILD_DIR) nlfem validate_cable

# Print resolved paths for debugging
info:
	@echo "PETSC_DIR     = $(PETSC_DIR)"
	@echo "PETSC_ARCH    = $(PETSC_ARCH)"
	@echo "PETSC_LIB_DIR = $(PETSC_LIB_DIR)"
	@echo "PARMETIS_DIR  = $(PARMETIS_DIR)"
	@echo "LAPACK_DIR    = $(LAPACK_DIR)"
	@echo ""
	@echo "Checking library files:"
	@ls $(PETSC_LIB_DIR)/libpetsc* 2>/dev/null || echo "  WARNING: no libpetsc* found in $(PETSC_LIB_DIR)"
	@ls $(PARMETIS_DIR)/lib/libparmetis* 2>/dev/null || echo "  WARNING: no libparmetis* found in $(PARMETIS_DIR)/lib"

NP ?= 4
run:
	mpirun -np $(NP) ./nlfem example/input.dat
