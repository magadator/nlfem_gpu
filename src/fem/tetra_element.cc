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
 * For a linear tet with nodes 0,1,2,3 at (x_i, y_i, z_i):
 *
 * The Jacobian of the mapping from barycentric to Cartesian coords is:
 *
 *   J = [x1-x0  x2-x0  x3-x0]
 *       [y1-y0  y2-y0  y3-y0]
 *       [z1-z0  z2-z0  z3-z0]
 *
 * det(J) = 6 * V  (V = element volume, positive if nodes ordered correctly)
 *
 * Shape function gradients in Cartesian coords are obtained from the
 * cofactors of J divided by det(J).  For node 0:
 *   dN0/dx = -(dN1/dx + dN2/dx + dN3/dx)   (partition of unity)
 *
 * B layout (6 rows × 12 cols), Voigt order
 * [eps_xx, eps_yy, eps_zz, gamma_yz, gamma_xz, gamma_xy]:
 *
 *   For node i (cols 3i, 3i+1, 3i+2 = ux_i, uy_i, uz_i):
 *     row 0 (eps_xx):   dNi/dx,  0,       0
 *     row 1 (eps_yy):   0,       dNi/dy,  0
 *     row 2 (eps_zz):   0,       0,       dNi/dz
 *     row 3 (gamma_yz): 0,       dNi/dz,  dNi/dy
 *     row 4 (gamma_xz): dNi/dz,  0,       dNi/dx
 *     row 5 (gamma_xy): dNi/dy,  dNi/dx,  0
 */

void comp_B_tetra(const TetraCoords& c, double B[72], double* vol)
{
    std::memset(B, 0, 72 * sizeof(double));

    /* Vectors from node 0 to nodes 1, 2, 3 — columns of J */
    const double a1 = c.x[1]-c.x[0],  b1 = c.y[1]-c.y[0],  d1 = c.z[1]-c.z[0];
    const double a2 = c.x[2]-c.x[0],  b2 = c.y[2]-c.y[0],  d2 = c.z[2]-c.z[0];
    const double a3 = c.x[3]-c.x[0],  b3 = c.y[3]-c.y[0],  d3 = c.z[3]-c.z[0];

    /* det(J) = 6V */
    const double detJ = a1*(b2*d3 - b3*d2)
                      - b1*(a2*d3 - a3*d2)
                      + d1*(a2*b3 - a3*b2);

    if (std::fabs(detJ) < 1.0e-30)
        throw std::runtime_error(
            "tetra_element::comp_B_tetra: degenerate element (|det J| < 1e-30). "
            "Check node ordering — nodes must not be coplanar.");

    *vol = detJ / 6.0;

    const double inv6V = 1.0 / detJ;   /* 1/(6V) */

    /*
     * Shape function gradients via cofactors of J.
     *
     * J = [a1 a2 a3]   (rows = x,y,z; cols = L1,L2,L3 directions)
     *     [b1 b2 b3]
     *     [d1 d2 d3]
     *
     * Cofactor expansion gives dN_i/d{x,y,z} for nodes 1,2,3.
     * Node 0 gradient follows from partition of unity: sum_i dN_i = 0.
     */

    /* Node 1 (L2 direction) */
    const double dN1dx =  (b2*d3 - b3*d2) * inv6V;
    const double dN1dy = -(a2*d3 - a3*d2) * inv6V;
    const double dN1dz =  (a2*b3 - a3*b2) * inv6V;

    /* Node 2 (L3 direction) */
    const double dN2dx = -(b1*d3 - b3*d1) * inv6V;
    const double dN2dy =  (a1*d3 - a3*d1) * inv6V;
    const double dN2dz = -(a1*b3 - a3*b1) * inv6V;

    /* Node 3 (L4 direction) */
    const double dN3dx =  (b1*d2 - b2*d1) * inv6V;
    const double dN3dy = -(a1*d2 - a2*d1) * inv6V;
    const double dN3dz =  (a1*b2 - a2*b1) * inv6V;

    /* Node 0 (L1 direction) — partition of unity */
    const double dN0dx = -(dN1dx + dN2dx + dN3dx);
    const double dN0dy = -(dN1dy + dN2dy + dN3dy);
    const double dN0dz = -(dN1dz + dN2dz + dN3dz);

    /* Fill B row by row for each node */
    auto fill_node = [&](int node, double dNdx, double dNdy, double dNdz) {
        const int c0 = 3 * node;
        /* row 0: eps_xx */   B[0*12 + c0+0] = dNdx;
        /* row 1: eps_yy */   B[1*12 + c0+1] = dNdy;
        /* row 2: eps_zz */   B[2*12 + c0+2] = dNdz;
        /* row 3: gamma_yz */ B[3*12 + c0+1] = dNdz;  B[3*12 + c0+2] = dNdy;
        /* row 4: gamma_xz */ B[4*12 + c0+0] = dNdz;  B[4*12 + c0+2] = dNdx;
        /* row 5: gamma_xy */ B[5*12 + c0+0] = dNdy;  B[5*12 + c0+1] = dNdx;
    };

    fill_node(0, dN0dx, dN0dy, dN0dz);
    fill_node(1, dN1dx, dN1dy, dN1dz);
    fill_node(2, dN2dx, dN2dy, dN2dz);
    fill_node(3, dN3dx, dN3dy, dN3dz);
}

/* ------------------------------------------------------------------ */
/*  Stiffness matrix  Ke = vol * B^T * C * B  (12×12)                 */
/* ------------------------------------------------------------------ */

void comp_Ke_tetra(const TetraCoords& coords,
                   const TetraMaterial& mat,
                   TetraKe& Ke)
{
    double B[72], C[36];
    double vol;

    comp_B_tetra(coords, B, &vol);
    comp_C_iso(mat.E, mat.nu, C);

    /* CB = C * B  (6×12) */
    double CB[72];
    std::memset(CB, 0, sizeof(CB));
    for (int i = 0; i < 6; ++i)
        for (int k = 0; k < 6; ++k)
            if (C[i*6+k] != 0.0)
                for (int j = 0; j < 12; ++j)
                    CB[i*12+j] += C[i*6+k] * B[k*12+j];

    /* Ke = vol * B^T * CB  (12×12) */
    std::memset(Ke.K, 0, sizeof(Ke.K));
    for (int i = 0; i < 12; ++i)
        for (int k = 0; k < 6; ++k)
            if (B[k*12+i] != 0.0)
                for (int j = 0; j < 12; ++j)
                    Ke.K[i*12+j] += vol * B[k*12+i] * CB[k*12+j];
}

/* ------------------------------------------------------------------ */
/*  Consistent mass matrix  Me  (12×12)                                */
/* ------------------------------------------------------------------ */
/*
 * Exact integration of rho * N^T * N over a linear tet:
 *
 *   integral_V N_i * N_j dV = V/20   if i == j   (same node)
 *                            = V/60   if i != j   (different nodes)
 *
 * The 12×12 matrix has a 4×4 block structure (one 3×3 identity-scaled
 * block per node pair):
 *
 *   Me[3*I+d1][3*J+d2] = rho * V * (1/20 if I==J, 1/60 if I!=J) * delta(d1,d2)
 *
 * where I,J in {0,1,2,3} are node indices and d1,d2 in {0,1,2} are
 * spatial DOF indices.
 */

void comp_Me_tetra(const TetraCoords& coords,
                   const TetraMaterial& mat,
                   TetraMe& Me)
{
    double B[72], vol;
    comp_B_tetra(coords, B, &vol);   /* only need vol */

    std::memset(Me.M, 0, sizeof(Me.M));

    const double diag_coef    = mat.rho * vol / 20.0;
    const double offdiag_coef = mat.rho * vol / 60.0;

    for (int I = 0; I < TETRA_NODES; ++I) {
        for (int J = 0; J < TETRA_NODES; ++J) {
            const double coef = (I == J) ? diag_coef : offdiag_coef;
            for (int d = 0; d < 3; ++d) {
                /* Only the diagonal in the spatial-DOF sub-block is non-zero
                 * (translational inertia, no coupling between x,y,z DOFs) */
                Me.M[(3*I+d)*12 + (3*J+d)] = coef;
            }
        }
    }
}

/* ------------------------------------------------------------------ */
/*  Internal force vector  Fe = vol * B^T * sigma  (12)               */
/* ------------------------------------------------------------------ */
/*
 * sigma = C * eps,  eps = B * u_e
 * Fe_i  = vol * sum_k  B[k][i] * sigma[k]
 *
 * This is the residual internal force that enters the Newton-Raphson
 * right-hand side as -Fe (work-conjugate to the virtual displacement).
 */

void comp_Fe_tetra(const TetraCoords& coords,
                   const TetraMaterial& mat,
                   const double u_e[TETRA_DOFS],
                   TetraFe& Fe)
{
    double B[72], C[36];
    double vol;

    comp_B_tetra(coords, B, &vol);
    comp_C_iso(mat.E, mat.nu, C);

    /* eps = B * u_e  (6-vector of strains, Voigt) */
    double eps[6];
    std::memset(eps, 0, sizeof(eps));
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 12; ++j)
            eps[i] += B[i*12+j] * u_e[j];

    /* sigma = C * eps  (6-vector of stresses, Voigt) */
    double sigma[6];
    std::memset(sigma, 0, sizeof(sigma));
    for (int i = 0; i < 6; ++i)
        for (int k = 0; k < 6; ++k)
            sigma[i] += C[i*6+k] * eps[k];

    /* Fe = vol * B^T * sigma  (12-vector) */
    std::memset(Fe.F, 0, sizeof(Fe.F));
    for (int i = 0; i < 12; ++i)
        for (int k = 0; k < 6; ++k)