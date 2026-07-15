/*
 * nlfem GPU element assembly — Phase 1 (CABLE elements)
 * Author: Gokul G. Anugrah (GPU port scaffolding by Claude, verified against
 *         src/fem/linearTriElement.cc comp_Ke_Cable / comp_Me_Cable)
 *
 * DESIGN NOTE
 * -----------
 * This module computes local element stiffness (Ke), mass (Me), and internal
 * force (Fvec) for ALL cable elements in one batched GPU kernel launch, then
 * hands the results back to fem::assembleGlobKMat() in exactly the layout it
 * already expects. PETSc assembly (MatSetValues) is UNCHANGED and still runs
 * on the host — only the per-element dense-matrix computation moves to GPU.
 *
 * This is intentionally element-type-scoped. Only CABLE (elType==5) is
 * GPU-accelerated in this phase. BEAM and SIMPLETRI/MITC3 fall through to
 * the existing CPU path (see comp_Ke_loc dispatch in fem.cc) until they are
 * ported with the same care.
 *
 * WHY BATCH ALL CABLE ELEMENTS AT ONCE:
 * Launch overhead on a K1200 dominates for small per-element work. Batching
 * every cable element into a single kernel launch amortizes that overhead
 * instead of doing one kernel launch per element (which would be strictly
 * slower than the CPU loop it replaces).
 */
#ifndef NLFEM_GPU_ASSEMBLY_H
#define NLFEM_GPU_ASSEMBLY_H

#include <vector>

namespace KARMA {

// Flat, struct-of-arrays input describing every cable element in the mesh.
// Populated on the host from fem's existing member arrays (nodes_, nodes0_,
// faces_, eModulus_, density_, Area_) — no new physics, just repacking.
struct GPUCableInput {
  int nCable = 0;              // number of cable elements this rank owns
  std::vector<int>    fn0, fn1;     // global node indices, size nCable
  std::vector<double> x0_0, x1_0;   // reference (undeformed) coords, 3*nCable each, xyz-interleaved
  std::vector<double> x0_t, x1_t;   // current (deformed) coords,    3*nCable each, xyz-interleaved
  std::vector<double> eModulus;     // size nCable
  std::vector<double> density;      // size nCable
  std::vector<double> area;         // size nCable
};

// Flat output: one 6x6 Ke, one 6x6 (diagonal, lumped) Me, one length-6 Fvec,
// and the scalar cable tension, per element — same quantities
// fem::assembleGlobKMat() already consumes per-element from the CPU path.
struct GPUCableOutput {
  std::vector<double> Ke;      // nCable * 36, row-major 6x6 per element
  std::vector<double> Me;      // nCable * 36, row-major 6x6 per element (diagonal populated)
  std::vector<double> Fvec;    // nCable * 6
  std::vector<double> tension; // nCable
};

// Runs the batched CUDA kernel for all cable elements in `in` and fills `out`.
// Safe to call with nCable == 0 (no-op). Throws std::runtime_error on any
// CUDA error (out-of-memory, launch failure, etc.) so the caller can decide
// whether to fall back to CPU rather than silently producing wrong results.
void gpuComputeCableElements(const GPUCableInput& in, GPUCableOutput& out);

// Returns true if a CUDA device is present and usable. fem.cc should call
// this once at startup and fall back to the CPU path entirely if false,
// rather than partially running on a device that isn't there.
bool gpuIsAvailable();

} // namespace KARMA

#endif
