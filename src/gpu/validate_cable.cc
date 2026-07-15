/*
 * validate_cable.cc
 *
 * Standalone check: computes Ke/Me/Fvec for a handful of cable elements
 * using (a) the exact same formulas as linearTriElement::comp_Ke_Cable /
 * comp_Me_Cable, reproduced here as plain CPU functions with no PETSc/MPI/
 * fem-class dependency, and (b) the GPU kernel, then reports max abs diff.
 *
 * This intentionally does NOT link against fem.cc/linearTriElement.cc (that
 * would require a full PETSc+MPI+ParMETIS build just to run a unit test).
 * Instead the CPU reference below is a direct, unmodified transcription of
 * the two functions -- if you change comp_Ke_Cable/comp_Me_Cable in
 * linearTriElement.cc, update the reference here to match, and re-run this
 * before trusting the GPU path on real data.
 *
 * Build:  make -f Makefile.gpu validate
 * Run:    ./validate_cable
 * Pass criterion: max abs diff on the order of 1e-12 or smaller (double
 * precision round-off). Anything larger means the kernel and CPU code have
 * diverged and the GPU path must not be trusted until it's fixed.
 */
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "gpu_assembly.h"

using namespace std;
using namespace KARMA;

static double vecMag3(const double* v) {
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

// Direct CPU transcription of comp_Ke_Cable + comp_Me_Cable for ONE element,
// taking plain coordinates instead of fem-class member state.
static void cpuReferenceCable(
    const double x0_0[3], const double x1_0[3],
    const double x0_t[3], const double x1_t[3],
    double emod, double dens, double area,
    double Ke[36], double Me[36], double Fvec[6], double& tension)
{
  double x2t[3], x20[3], xhat[6], uhat[6];
  for (int i = 0; i < 3; ++i) {
    x2t[i] = x1_t[i] - x0_t[i];
    x20[i] = x1_0[i] - x0_0[i];
    xhat[i] = x0_0[i]; xhat[i+3] = x1_0[i];
    uhat[i] = x0_t[i]-x0_0[i]; uhat[i+3] = x1_t[i]-x1_0[i];
  }
  double L0 = vecMag3(x20);
  double defaultLength = vecMag3(x2t);
  double deltaL = defaultLength - L0;
  double erst = deltaL/L0 + 0.5*(deltaL/L0)*(deltaL/L0);

  double HrTHr[36] = {0};
  HrTHr[0*6+0]=1.; HrTHr[1*6+1]=1.; HrTHr[2*6+2]=1.;
  HrTHr[3*6+3]=1.; HrTHr[4*6+4]=1.; HrTHr[5*6+5]=1.;
  HrTHr[0*6+3]=-1.; HrTHr[1*6+4]=-1.; HrTHr[2*6+5]=-1.;
  HrTHr[3*6+0]=-1.; HrTHr[4*6+1]=-1.; HrTHr[5*6+2]=-1.;

  double BL[6] = {0};
  for (int i=0;i<6;++i) for (int j=0;j<6;++j) BL[i] += 1./(L0*L0)*xhat[j]*HrTHr[j*6+i];
  for (int i=0;i<6;++i) for (int j=0;j<6;++j) BL[i] += 1./(L0*L0)*uhat[i]*HrTHr[j*6+i];

  double cosx=(x1_t[0]-x0_t[0])/defaultLength;
  double cosy=(x1_t[1]-x0_t[1])/defaultLength;
  double cosz=(x1_t[2]-x0_t[2])/defaultLength;

  double T[36];
  T[0*6+0]=cosx*cosx; T[0*6+1]=cosx*cosy; T[0*6+2]=cosx*cosz; T[0*6+3]=-cosx*cosx; T[0*6+4]=-cosx*cosy; T[0*6+5]=-cosx*cosz;
  T[1*6+0]=cosx*cosy; T[1*6+1]=cosy*cosy; T[1*6+2]=cosy*cosz; T[1*6+3]=-cosx*cosy; T[1*6+4]=-cosy*cosy; T[1*6+5]=-cosy*cosz;
  T[2*6+0]=cosx*cosz; T[2*6+1]=cosy*cosz; T[2*6+2]=cosz*cosz; T[2*6+3]=-cosx*cosz; T[2*6+4]=-cosy*cosz; T[2*6+5]=-cosz*cosz;
  T[3*6+0]=-cosx*cosx; T[3*6+1]=-cosx*cosy; T[3*6+2]=-cosx*cosz; T[3*6+3]=cosx*cosx; T[3*6+4]=cosx*cosy; T[3*6+5]=cosx*cosz;
  T[4*6+0]=-cosx*cosy; T[4*6+1]=-cosy*cosy; T[4*6+2]=-cosy*cosz; T[4*6+3]=cosx*cosy; T[4*6+4]=cosy*cosy; T[4*6+5]=cosy*cosz;
  T[5*6+0]=-cosx*cosz; T[5*6+1]=-cosy*cosz; T[5*6+2]=-cosz*cosz; T[5*6+3]=cosx*cosz; T[5*6+4]=cosy*cosz; T[5*6+5]=cosz*cosz;

  double btb[36];
  for (int i=0;i<6;++i) for (int j=0;j<6;++j) btb[i*6+j] = BL[j]*BL[j]*emod;

  for (int i=0;i<36;++i) Ke[i]=0.;
  for (int i=0;i<6;++i) for (int j=0;j<6;++j) for (int k=0;k<6;++k)
    Ke[i*6+j] += T[k*6+i]*btb[i*6+k]*T[k*6+j]*L0*area;

  double Srst = emod*max(erst,0.);
  tension = Srst*area;
  double ftmp[6];
  for (int i=0;i<6;++i) ftmp[i]=BL[i]*Srst*L0*area;
  for (int i=0;i<6;++i) Fvec[i]=0.;
  for (int i=0;i<6;++i) for (int j=0;j<6;++j) Fvec[i]+=T[j*6+i]*ftmp[j];

  for (int i=0;i<36;++i) Me[i]=0.;
  double massDiag = 0.5*dens*area*L0;
  for (int i=0;i<6;++i) Me[i*6+i]=massDiag;
}

int main() {
  // A handful of synthetic elements: straight, angled, and pre-stretched.
  const int n = 4;
  GPUCableInput in;
  in.nCable = n;
  in.fn0 = {0,0,0,0};
  in.fn1 = {1,1,1,1};
  in.x0_0 = {0,0,0,  0,0,0,   0,0,0,   0,0,0};
  in.x1_0 = {1,0,0,  1,1,0,   0,0,2,   1,0,0};
  // element 3 (index 3) is pre-stretched: current length != reference length
  in.x0_t = {0,0,0,  0,0,0,   0,0,0,   0,0,0};
  in.x1_t = {1,0,0,  1.01,1,0.02, 0,0,2.05, 1.10,0,0};
  in.eModulus = {2.0e5, 1.5e5, 1.0e5, 3.0e5};
  in.density  = {7800, 7800, 1200, 7800};
  in.area     = {1e-4, 2e-4, 5e-5, 1e-4};

  GPUCableOutput out;
  if (!gpuIsAvailable()) {
    fprintf(stderr, "No CUDA device visible -- build/run this on a machine with the GPU enabled.\n");
    return 1;
  }
  gpuComputeCableElements(in, out);

  double maxDiff = 0.0;
  for (int e = 0; e < n; ++e) {
    double x00[3] = {in.x0_0[e*3+0], in.x0_0[e*3+1], in.x0_0[e*3+2]};
    double x10[3] = {in.x1_0[e*3+0], in.x1_0[e*3+1], in.x1_0[e*3+2]};
    double x0t[3] = {in.x0_t[e*3+0], in.x0_t[e*3+1], in.x0_t[e*3+2]};
    double x1t[3] = {in.x1_t[e*3+0], in.x1_t[e*3+1], in.x1_t[e*3+2]};
    double Ke[36], Me[36], Fvec[6], tension;
    cpuReferenceCable(x00, x10, x0t, x1t, in.eModulus[e], in.density[e], in.area[e],
                       Ke, Me, Fvec, tension);

    for (int i = 0; i < 36; ++i) maxDiff = max(maxDiff, fabs(Ke[i] - out.Ke[e*36+i]));
    for (int i = 0; i < 36; ++i) maxDiff = max(maxDiff, fabs(Me[i] - out.Me[e*36+i]));
    for (int i = 0; i < 6;  ++i) maxDiff = max(maxDiff, fabs(Fvec[i] - out.Fvec[e*6+i]));
    maxDiff = max(maxDiff, fabs(tension - out.tension[e]));

    printf("element %d: tension CPU=%.10e GPU=%.10e\n", e, tension, out.tension[e]);
  }
  printf("\nmax abs diff (Ke, Me, Fvec, tension) across %d elements: %.3e\n", n, maxDiff);
  printf(maxDiff < 1e-9 ? "PASS\n" : "FAIL -- do not trust the GPU path yet\n");
  return maxDiff < 1e-9 ? 0 : 1;
}
