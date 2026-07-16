#ifndef TETRA_ELEMENT_H
#define TETRA_ELEMENT_H

/*
 * tetra_element.h
 *
 * Linear 4-node tetrahedral element for 3-D solid mechanics.
 *
 * Degrees of freedom per node : 3  (ux, uy, uz)
 * Total DOFs per element      : 12
 * Element type ID (elType)    : 10  (TETRA, added alongside CABLE=5, MITC3=4, …)
 *
 * Shape functions (barycentric / volume coords L0..L3, sum = 1):
 *   N_i = L_i   (i = 0,1,2,3)
 *
 * The strain-displacement matrix B is CONSTANT over the element
 * (constant-strain tet, CST-3D), so Ke = V * B^T C B exactly with no
 * quadrature needed.
 *
 * Consistent mass matrix (exact integration over a linear tet):
 *   M_scalar_ij = rho * V / 20   (i == j, scalar node-node block)
 *               = rho * V / 60   (i != j, scalar node-node block)
 * The 12×12 matrix repeats this 4×4 scalar pattern for x, y, z DOFs.
 *
 * Material model: linear isotropic elasticity.
 * Voigt stress/strain order: [s_xx, s_yy, s_zz, s_yz, s_xz, s_xy]
 *
 * References:
 *   Zienkiewicz & Taylor, "The Finite Element Method", Vol 1, 6th ed., §4.3
 *   Hughes, "The Finite Element Method", Dover, §3.
 */

#include <cmath>
#include <stdexcept>

namespace nlfem {

/* ------------------------------------------------------------------ */
/*  Compile-time constants                                              */
/* ------------------------------------------------------------------ */
static constexpr int TETRA_NODES   = 4;
static constexpr int TETRA_DOFS    = 12;   // 4 nodes × 3 DOFs
static constexpr int TETRA_STRAINS = 6;    // Voigt notation

/* ------------------------------------------------------------------ */
/*  Plain-old-data structs (shared by CPU and GPU code)                */
/* ------------------------------------------------------------------ */

/** Physical node coordinates for one element (reference or current). */
struct TetraCoords {
    double x[4], y[4], z[4];   // node 0..3
};

/** Isotropic linear-elastic material parameters. */
struct TetraMaterial {
    double E;    ///< Young's modulus  [Pa]
    double nu;   ///< Poisson's ratio  [-]
    double rho;  ///< mass density     [kg/m^3]
};

/** 12×12 element stiffness matrix (row-major). */
struct TetraKe {
    double K[TETRA_DOFS * TETRA_DOFS];
};

/** 12×12 consistent element mass matrix (row-major). */
struct TetraMe {
    double M[TETRA_DOFS * TETRA_DOFS];
};

/** 12-component internal force vector. */
struct TetraFe {
    double F[TETRA_DOFS];
};

/* ------------------------------------------------------------------ */
/*  CPU reference functions (implemented in tetra_element.cc)          */
/* ------------------------------------------------------------------ */

/**
 * Build the 6×6 isotropic constitutive matrix C (Voigt notation).
 * Stored row-major in C[36].
 *
 * Voigt order: [eps_xx, eps_yy, eps_zz, gamma_yz, gamma_xz, gamma_xy]
 */
void comp_C_iso(double E, double nu, double C[36]);

/**
 * Build the constant 6×12 strain-displacement matrix B for a linear tet.
 * Stored row-major in B[72].
 * Also returns the signed element volume via *vol (positive when nodes are
 * ordered so that det(J) > 0).
 *
 * Throws std::runtime_error if the element is degenerate (|det J| < 1e-30).
 */
void comp_B_tetra(const TetraCoords& coords, double B[72], double* vol);

/**
 * Compute 12×12 stiffness matrix  Ke = vol * B^T * C * B.
 */
void comp_Ke_tetra(const TetraCoords& coords,
                   const TetraMaterial& mat,
                   TetraKe& Ke);

/**
 * Compute 12×12 consistent mass matrix using exact integration.
 */
void comp_Me_tetra(const TetraCoords& coords,
                   const TetraMaterial& mat,
                   TetraMe& Me);

/**
 * Compute 12-vector internal force  Fe = vol * B^T * C * B * u_e.
 * u_e[12] is the current element displacement vector in DOF order
 * [ux0,uy0,uz0, ux1,uy1,uz1, ux2,uy2,uz2, ux3,uy3,uz3].
 */
void comp_Fe_tetra(const TetraCoords& coords,
                   const TetraMaterial& mat,
                   const double u_e[TETRA_DOFS],
                   TetraFe& Fe);

} // namespace nlfem

#endif // TETRA_ELEMENT_H
