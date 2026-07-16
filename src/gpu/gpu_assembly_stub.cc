/*
 * Stub implementation used when nlfem is built with USEGPU=0 (no nvcc/CUDA
 * available). Lets the rest of the codebase call gpuIsAvailable(),
 * gpuComputeCableElements(), and gpuComputeTetraElements() unconditionally
 * without needing #ifdef USEGPU scattered through fem.cc.
 * gpuIsAvailable() always returns false here, so fem.cc's dispatch logic
 * falls back to the CPU path automatically.
 */
#include "gpu_assembly.h"
#include <stdexcept>

namespace KARMA {

bool gpuIsAvailable() { return false; }

void gpuComputeCableElements(const GPUCableInput&, GPUCableOutput&) {
  throw std::runtime_error(
      "gpuComputeCableElements() called but nlfem was built with USEGPU=0. "
      "This should be unreachable -- fem.cc must check gpuIsAvailable() "
      "before calling this.");
}

void gpuComputeTetraElements(const GPUTetraInput&, GPUTetraOutput&) {
  throw std::runtime_error(
      "gpuComputeTetraElements() called but nlfem was built with USEGPU=0. "
      "This should be unreachable -- fem.cc must check gpuIsAvailable() "
      "before calling this.");
}

} // namespace KARMA
