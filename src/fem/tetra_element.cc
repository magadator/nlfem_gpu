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

#include "te