/*
 * nlfem - Nonlinear FEM Solver
 * Author  : Gokul G. Anugrah
 * Contact : gokulanugrah@gmail.com
 *
 * linearTriElement.h — Concrete element class declaration.
 *   Derives from fem. Implements comp_Ke_loc() and comp_Me_loc()
 *   for MITC3 shell, beam, and cable element types.
 */
#ifndef KARMA_LINEAR_TRI_ELEMENT_H
#define KARMA_LINEAR_TRI_ELEMENT_H
//
#include "fem.h"
#include "linearTriElementFort.h"
#include "linear_algebra.h"
#include "surfaceTriangulation.h"

using namespace std;

namespace KARMA {
  
  class linearTriElement: public fem {
  public:
    void comp_Ke_loc     (vector<double>& mat,vector<double>& Fvec,const int& triID,const int& elType);
    void comp_Ke_TriShell(vector<double>& mat,vector<double>& Fvec,const int& triID);
    void comp_Ke_MITC3   (vector<double>& mat,vector<double>& Fvec,const int& triID);
    void comp_NL_Ke_MITC3(vector<double>& mat,vector<double>& Fvec,const int& triID);
    void comp_Ke_BEBeam  (vector<double>& mat,vector<double>& Fvec,const int& triID);
    void comp_Ke_Cable   (vector<double>& mat,vector<double>& Fvec,const int& triID);
    //
    void comp_Me_loc     (vector<double>& mat,const int& triID,const int& elType);
    void comp_Me_TriShell(vector<double>& mat,const int& triID);
    void comp_Me_MITC3   (vector<double>& mat,const int& triID);
    void comp_NL_Me_MITC3(vector<double>& mat,const int& triID);
    void comp_Me_BEBeam  (vector<double>& mat,const int& triID);
    void comp_Me_Cable   (vector<double>& mat,const int& triID);
    //
    void compBTCB(vector<vector<double> >& BTCBmat_,
                  const int& faceID,
                  double& rPoint,
                  double& sPoint,
                  double& tPoint,
                  const double& A,
                  vector<double> x0,
                  vector<double> x1,
                  vector<double> x2,
                  vector<double> ex,
                  vector<double> ey,   
                  vector<double> ez);   
    void comp_NL_BTCB(const int& faceID,
                      double& rPoint,
                      double& sPoint,
                      double& tPoint);
    void calculateStrainsMITC3(double& rPoint,
                               double& sPoint,
                               double& tPoint,
                               double& a1,
                               double& a2,
                               double& a3,
                               const int& faceID);
    void comp_KNL1_MITC3(const int& faceID,
                         double& rPoint,
                         double& sPoint,
                         double& tPoint);
    void comp_KNL2_MITC3(const int& faceID,
                         double& rPoint,
                         double& sPoint,
                         double& tPoint);
    void compPKST(const int& faceID);
    void calculateDeformationGradient(double& rPoint,
				      double& sPoint,
				      double& tPoint,
				      double& a0,
				      double& a1,
				      double& a2,
				      const int& faceID);
    void PKST2Cauchy(const int& faceID);
    void calculateVonMisesStress(const int& faceID);
    //create tying point vectors to fill B_rst using MITC scheme
    int faceID;
    vector<double> x0_0;
    vector<double> x1_0;
    vector<double> x2_0;
    //
    void assignExactValuesPerFace(vector<double>& u,
                                  int& fn0,
                                  int& fn1,
                                  int& fn2);
    //
    void compKloc3Dnew(vector<vector<double> >& mat,const int& triID);
    void compMloc3Dnew(vector<vector<double> >& mat,const int& triID);
    void compShapeFunction(vector<vector<double> >& Nmat,
			   double&                  area,  
			   const int&             triID,
			   const vector<double>&    LVec); 
    void compDloc(vector<vector<double> >& D,const double& A,const int& triID);
    void compddPloc(vector<vector<double> >& KKl,
		    const vector<double>&    mu,
		    const vector<double>&    L,
		    const double&            Delta,
		    const vector<double>&    a,
		    const vector<double>&    b,
		    const vector<double>&    c,
		    const vector<double>&    x0,
		    const vector<double>&    x1,
		    const vector<double>&    x2,
		    const int&               faceID);    
  };
}
#endif   
