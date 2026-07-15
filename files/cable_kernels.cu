/*
 * nlfem GPU element assembly — cable element kernel
 *
 * This is a literal, one-thread-per-element port of:
 *   linearTriElement::comp_Ke_Cable()   (src/fem/linearTriElement.cc:80)
 *   linearTriElement::comp_Me_Cable()   (src/fem/linearTriElement.cc:257)
 *
 * PORTING NOTES (read before touching the math below):
 *
 * 1. The CPU version relies on execution order: comp_Ke_Cable(f) sets member
 *    scratch x0_t_/x1_t_ for element f, and the very next call in the same
 *    loop iteration, comp_Me_Cable(f), reads those same member variables
 *    without recomputing them. That's safe in the serial per-element loop
 *    but is NOT safe for independent parallel threads. This kernel
 *    recomputes x0_t_/x1_t_ (current node positions) explicitly for both
 *    the stiffness and mass calculation, per thread. This is the one
 *    semantic change from the original code, and it's required for
 *    correctness under parallel execution, not an optimization.
 *
 * 2. Everything else below — including quirks like `btb[i*6+j] =
 *    BL[j]*BL[j]*emod` (uses BL[j] twice, not BL[i]*BL[j]) and the fact that
 *    comp_Me_Cable's HrTHr_/T_ construction is dead code (only the diagonal
 *    lumped-mass line actually feeds Me_loc) — is preserved EXACTLY as
 *    written in the original, including the parts that look like they might
 *    be bugs. Changing element-formulation math is not this pass's job;
 *    only Gokul should decide whether BL[j]*BL[j] vs BL[i]*BL[j] is
 *    intentional. Flag this to him, don't silently "fix" it.
 *
 * VALIDATION: build with `make gpu-validate` (see README) to run this
 * kernel and the original CPU function side by side on the same mesh and
 * print the max abs difference in Ke/Me/Fvec. Do not trust this kernel's
 * output until that diff is at floating-point round-off level.
 */

#include <cuda_runtime.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include "gpu_assembly.h"

namespace KARMA {

#define CUDA_CHECK(call)                                                     \
  do {                                                                       \
    cudaError_t _e = (call);                                                 \
    if (_e != cudaSuccess) {                                                 \
      throw std::runtime_error(std::string("CUDA error at ") + __FILE__ +    \
                                ":" + std::to_string(__LINE__) + " - " +     \
                                cudaGetErrorString(_e));                     \
    }                                                                        \
  } while (0)

__device__ __forceinline__ double vecMag3(const double* v) {
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// One thread == one cable element. Work per element is tiny (6x6 dense),
// so thread-per-element (not block-per-element) is the right granularity
// here -- keeps occupancy high on a small part like the K1200.
__global__ void cableElementKernel(
    int nCable,
    const int*    fn0, const int* fn1,
    const double* x0_0, const double* x1_0,   // reference coords, xyz-interleaved, size 3*nCable
    const double* x0_t, const double* x1_t,   // current coords,   xyz-interleaved, size 3*nCable
    const double* eModulus, const double* density, const double* area,
    double* Ke_out, double* Me_out, double* Fvec_out, double* tension_out)
{
  int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= nCable) return;

  // ---- gather this element's data (SoA -> local registers) ----
  double x0t[3], x1t[3], x00[3], x10[3], x2t[3], x20[3];
  double xhat[6], uhat[6];
  #pragma unroll
  for (int idir = 0; idir < 3; ++idir) {
    x0t[idir] = x0_t[e*3 + idir];
    x1t[idir] = x1_t[e*3 + idir];
    x00[idir] = x0_0[e*3 + idir];
    x10[idir] = x1_0[e*3 + idir];
    x2t[idir] = x1t[idir] - x0t[idir];
    x20[idir] = x10[idir] - x00[idir];
    xhat[idir]     = x00[idir];
    xhat[idir + 3] = x10[idir];
    uhat[idir]     = x0t[idir] - x00[idir];
    uhat[idir + 3] = x1t[idir] - x10[idir];
  }

  double L0 = vecMag3(x20);
  double defaultLength = vecMag3(x2t);
  double deltaL = defaultLength - L0;
  double emod = eModulus[e];
  double erst = deltaL / L0 + 0.5 * (deltaL / L0) * (deltaL / L0);

  // HrTHr for stiffness (symmetric 6x6, row-major, M[i*6+j])
  double HrTHr[36];
  #pragma unroll
  for (int i = 0; i < 36; ++i) HrTHr[i] = 0.0;
  HrTHr[0*6+0]=1.; HrTHr[1*6+1]=1.; HrTHr[2*6+2]=1.;
  HrTHr[3*6+3]=1.; HrTHr[4*6+4]=1.; HrTHr[5*6+5]=1.;
  HrTHr[0*6+3]=-1.; HrTHr[1*6+4]=-1.; HrTHr[2*6+5]=-1.;
  HrTHr[3*6+0]=-1.; HrTHr[4*6+1]=-1.; HrTHr[5*6+2]=-1.;

  double BL[6] = {0,0,0,0,0,0};
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      BL[i] += 1.0/(L0*L0) * xhat[j] * HrTHr[j*6 + i];
    }
  }
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      BL[i] += 1.0/(L0*L0) * uhat[i] * HrTHr[j*6 + i];
    }
  }

  // direction cosines from CURRENT geometry
  double cosx = (x1t[0] - x0t[0]) / defaultLength;
  double cosy = (x1t[1] - x0t[1]) / defaultLength;
  double cosz = (x1t[2] - x0t[2]) / defaultLength;

  double T[36];
  T[0*6+0]=cosx*cosx;  T[0*6+1]=cosx*cosy;  T[0*6+2]=cosx*cosz;
  T[0*6+3]=-cosx*cosx; T[0*6+4]=-cosx*cosy; T[0*6+5]=-cosx*cosz;
  T[1*6+0]=cosx*cosy;  T[1*6+1]=cosy*cosy;  T[1*6+2]=cosy*cosz;
  T[1*6+3]=-cosx*cosy; T[1*6+4]=-cosy*cosy; T[1*6+5]=-cosy*cosz;
  T[2*6+0]=cosx*cosz;  T[2*6+1]=cosy*cosz;  T[2*6+2]=cosz*cosz;
  T[2*6+3]=-cosx*cosz; T[2*6+4]=-cosy*cosz; T[2*6+5]=-cosz*cosz;
  T[3*6+0]=-cosx*cosx; T[3*6+1]=-cosx*cosy; T[3*6+2]=-cosx*cosz;
  T[3*6+3]=cosx*cosx;  T[3*6+4]=cosx*cosy;  T[3*6+5]=cosx*cosz;
  T[4*6+0]=-cosx*cosy; T[4*6+1]=-cosy*cosy; T[4*6+2]=-cosy*cosz;
  T[4*6+3]=cosx*cosy;  T[4*6+4]=cosy*cosy;  T[4*6+5]=cosy*cosz;
  T[5*6+0]=-cosx*cosz; T[5*6+1]=-cosy*cosz; T[5*6+2]=-cosz*cosz;
  T[5*6+3]=cosx*cosz;  T[5*6+4]=cosy*cosz;  T[5*6+5]=cosz*cosz;

  // preserved as-in-original: btb[i*6+j] = BL[j]*BL[j]*emod (see file header note)
  double btb[36];
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      btb[i*6+j] = BL[j] * BL[j] * emod;

  double areaE = area[e];
  double Ke[36] = {0};
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      for (int k = 0; k < 6; ++k)
        Ke[i*6+j] += T[k*6+i] * btb[i*6+k] * T[k*6+j] * L0 * areaE;

  double Srst = emod * fmax(erst, 0.0);
  tension_out[e] = Srst * areaE;

  double ftmp[6];
  for (int i = 0; i < 6; ++i) ftmp[i] = BL[i] * Srst * L0 * areaE;

  double Fvec[6] = {0};
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      Fvec[i] += T[j*6+i] * ftmp[j];

  // ---- mass: lumped diagonal only (matches original -- see file header) ----
  double Me[36] = {0};
  double massDiag = 0.5 * density[e] * areaE * L0;
  #pragma unroll
  for (int i = 0; i < 6; ++i) Me[i*6+i] = massDiag;

  // ---- write out ----
  #pragma unroll
  for (int i = 0; i < 36; ++i) {
    Ke_out[e*36 + i] = Ke[i];
    Me_out[e*36 + i] = Me[i];
  }
  #pragma unroll
  for (int i = 0; i < 6; ++i) Fvec_out[e*6 + i] = Fvec[i];
}

bool gpuIsAvailable() {
  int count = 0;
  cudaError_t e = cudaGetDeviceCount(&count);
  return (e == cudaSuccess) && (count > 0);
}

void gpuComputeCableElements(const GPUCableInput& in, GPUCableOutput& out) {
  const int n = in.nCable;
  out.Ke.assign(n * 36, 0.0);
  out.Me.assign(n * 36, 0.0);
  out.Fvec.assign(n * 6, 0.0);
  out.tension.assign(n, 0.0);
  if (n == 0) return;

  int *d_fn0, *d_fn1;
  double *d_x0_0, *d_x1_0, *d_x0_t, *d_x1_t;
  double *d_eMod, *d_density, *d_area;
  double *d_Ke, *d_Me, *d_Fvec, *d_tension;

  CUDA_CHECK(cudaMalloc(&d_fn0, n * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&d_fn1, n * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&d_x0_0, 3 * n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_x1_0, 3 * n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_x0_t, 3 * n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_x1_t, 3 * n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_eMod, n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_density, n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_area, n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_Ke, 36 * n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_Me, 36 * n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_Fvec, 6 * n * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_tension, n * sizeof(double)));

  CUDA_CHECK(cudaMemcpy(d_fn0, in.fn0.data(), n * sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_fn1, in.fn1.data(), n * sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_x0_0, in.x0_0.data(), 3*n*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_x1_0, in.x1_0.data(), 3*n*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_x0_t, in.x0_t.data(), 3*n*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_x1_t, in.x1_t.data(), 3*n*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_eMod, in.eModulus.data(), n*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_density, in.density.data(), n*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_area, in.area.data(), n*sizeof(double), cudaMemcpyHostToDevice));

  int threads = 128;
  int blocks = (n + threads - 1) / threads;
  cableElementKernel<<<blocks, threads>>>(
      n, d_fn0, d_fn1, d_x0_0, d_x1_0, d_x0_t, d_x1_t,
      d_eMod, d_density, d_area, d_Ke, d_Me, d_Fvec, d_tension);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  CUDA_CHECK(cudaMemcpy(out.Ke.data(), d_Ke, 36*n*sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(out.Me.data(), d_Me, 36*n*sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(out.Fvec.data(), d_Fvec, 6*n*sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(out.tension.data(), d_tension, n*sizeof(double), cudaMemcpyDeviceToHost));

  cudaFree(d_fn0); cudaFree(d_fn1);
  cudaFree(d_x0_0); cudaFree(d_x1_0); cudaFree(d_x0_t); cudaFree(d_x1_t);
  cudaFree(d_eMod); cudaFree(d_density); cudaFree(d_area);
  cudaFree(d_Ke); cudaFree(d_Me); cudaFree(d_Fvec); cudaFree(d_tension);
}

} // namespace KARMA
