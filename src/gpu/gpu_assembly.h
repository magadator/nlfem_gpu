/*
 * nlfem GPU element assembly — Phase 1 (CABLE) + Phase 2 (TETRA)
 * Author: Gokul G. Anugrah (GPU port scaffolding by Claude, verified against
 *         src/fem/linearTriElement.cc comp_Ke_Cable / comp_Me_Cable and
 *         src/fem/tetra_element.cc comp_Ke_tetra / comp_Me_tetra)
 *
 * DESIGN NOTE
 * -----------
 * This module computes local element stiffness (Ke), mass (Me), and internal
 * force (Fvec) for ALL cable elements AND all tetrahedral elements in one
 * batched GPU kernel launch per type, then hands the results back to
 * fem::assembleGlobKMat() in exactly the layout it already expects.
 * PETSc assembly (MatSetValues) is UNCHANGED and still runs on the host —
 * only the per-element dense-matrix computation moves to GPU.
 *
 * CABLE (elType==5):  6 DOFs/element,  6×6  Ke/Me
 * TETRA (elType==10): 12 DOFs/element, 12×12 Ke/Me  (linear 4-node tet)
 *
 * WHY BATCH ALL ELEMENTS OF ONE TYPE AT ONCE:
 * Launch overhead on a K1200 dominates for small per-element work. Batching
 * every element into a single kernel launch amortises that overhead instead
 * of doing one kernel launch per element (which would be strictly slower than
 * the CPU loop it replaces).
 */
#ifndef NLFEM_GPU_ASSEMBLY_H
#define NLFEM_GPU_ASSEMBLY_H

#include <vector>

namespace KARMA {

/* ================================================================== */
/*  CABLE element (elType == 5)                                        */
/* ================================================================== */

// Flat, struct-of-arrays input describing every cable element in the mesh.
// Populated on the host from fem's existing member arrays (nodes_, nodes0_,
// faces_, eModulus_, density_, Area_) — no new physics, just repacking.
struct GPUCableInput {
  int nCable = 0;
  std::vector<int>    fn0, fn1;
  std::vector<double> x0_0, x1_0;   // reference coords, 3*nCable, xyz-interleaved
  std::vector<double> x0_t, x1_t;   // current  coords, 3*nCable, xyz-interleaved
  std::vector<double> eModulus;
  std::vector<double> density;
  std::vector<double> area;
};

// Flat output: one 6×6 Ke, one 6×6 Me (diagonal), one 6-vec Fvec,
// and scalar tension, per element.
struct GPUCableOutput {
  std::vector<double> Ke;      // nCable * 36
  std::vector<double> Me;      // nCable * 36
  std::vector<double> Fvec;    // nCable * 6
  std::vector<double> tension; // nCable
};

// Runs the batched CUDA kernel for all cable elements in `in` and fills `out`.
// Safe to call with nCable == 0 (no-op). Throws std::runtime_error on any
// CUDA error so the caller can fall back to CPU.
void gpuComputeCableElements(const GPUCableInput& in, GPUCableOutput& out);

/* ================================================================== */
/*  TETRA element (elType == 10) — linear 4-node tetrahedron           */
/*                                                                      */
/*  DOF ordering per element:                                           */
/*    [ux0,uy0,uz0, ux1,uy1,uz1, ux2,uy2,uz2, ux3,uy3,uz3]           */
/*    indices 0..11                                                     */
/*                                                                      */
/*  Ke: 12×12 stiffness  (V * B^T C B, B constant for linear tet)     */
/*  Me: 12×12 consistent mass  (rho*V * integral N^T N dV)            */
/*  Fe: 12    internal force   (V * B^T * C * B * u_e)                */
/*                                                                      */
/*  Voigt stress/strain order:                                          */
/*    [s_xx, s_yy, s_zz, s_yz, s_xz, s_xy]                           */
/* ================================================================== */

// Struct-of-arrays input for all tetra elements owned by this MPI rank.
//
// coords0 / coordst: reference / current node coordinates, stored as
//   coords[e*12 + n*3 + {0,1,2}] = x,y,z of node n of element e
//
// u_e: current element displacement vector (same DOF ordering as Ke),
//   u_e[e*12 + dof]
//
// nodeIDs: global node indices, nodeIDs[e*4 + n] = global ID of node n
//   of element e (used by the host to scatter Fe back into the global
//   force vector; the GPU kernel itself does not need them).
struct GPUTetraInput {
  int nTetra = 0;
  std::vector<int>    nodeIDs;   // nTetra * 4
  std::vector<double> coords0;   // nTetra * 12  (reference)
  std::vector<double> coordst;   // nTetra * 12  (current)
  std::vector<double> u_e;       // nTetra * 12  (current displacements)
  std::vector<double> eModulus;  // nTetra
  std::vector<double> nu;        // nTetra  (Poisson's ratio)
  std::vector<double> density;   // nTetra
};

// Flat output per tetra element.
struct GPUTetraOutput {
  std::vector<double> Ke;   // nTetra * 144  (row-major 12×12)
  std::vector<double> Me;   // nTetra * 144  (row-major 12×12)
  std::vector<double> Fe;   // nTetra * 12   (internal force)
};

// Runs the batched CUDA kernel for all tetra elements in `in` and fills `out`.
// Safe to call with nTetra == 0 (no-op). Throws std::runtime_error on any
// CUDA error.
void gpuComputeTetraElements(const GPUTetraInput& in, GPUTetraOutput& out);

/* ================================================================== */
/*  Utility                                                             */
/* ================================================================== */

// Returns true if a CUDA device is present and usable at runtime.
// fem.cc calls this once at startup; if false the CPU path is used for
// all element types without any further GPU calls.
bool gpuIsAvailable();

} // namespace KARMA

#endif // NLFEM_GPU_ASSEMBLY_H
