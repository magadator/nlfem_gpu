/*
 * tetra_element.cc
 *
 * CPU reference implementation of the linear 4-node tetrahedral element.
 *
 * All formulas follow Zienkiewicz & Taylor Vol 1, §4.3 (constant-strain tet).
 *
 * DOF ordering per element:
 *   [ux0, uy0, uz0,  ux1, uy1, uz1,  ux2, uy2, uz2,  ux3, uy3, uz3]
 *    col:  0    1    2     3    4    5     6    7    8     9   10   11
 *
 * Voigt strain/stress order:
 *   [eps_xx, eps_yy, eps_zz, gamma_yz, gamma_xz, gamma_xy]
 *   row:  0       1       2        3        4        5
 */

#include "tetra_element.h"
#include <cstring>
#include <cmath>
#include <stdexcept>

namespace nlfem {

/* ------------------------------------------------------------------ */
/*  Constitutive matrix  C  (6×6, Voigt, isotropic linear elastic)    */
/* ------------------------------------------------------------------ */

void comp_C_iso(double E, double nu, double C[36])
{
    std::memset(C, 0, 36 * sizeof(double));

    const double lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu  = E / (2.0 * (1.0 + nu));

    /* Normal-normal block (rows/cols 0,1,2) */
    C[0*6+0] = lam + 2.0*mu;  C[0*6+1] = lam;            C[0*6+2] = lam;
    C[1*6+0] = lam;            C[1*6+1] = lam + 2.0*mu;  C[1*6+2] = lam;
    C[2*6+0] = lam;            C[2*6+1] = lam;            C[2*6+2] = lam + 2.0*mu;

    /* Shear block (rows/cols 3,4,5) — diagonal only */
    C[3*6+3] = mu;
    C[4*6+4] = mu;
    C[5*6+5] = mu;
}

/* ------------------------------------------------------------------ */
/*  Strain-displacement matrix  B  (6×12, constant over element)       */
/* ------------------------------------------------------------------ */
/*
 * For a linear tet with nodes 0,1,2