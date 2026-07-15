#ifndef KARMA_LINEAR_BEAM_ELEMENT_H
#define KARMA_LINEAR_BEAM_ELEMENT_H
//
#include "fem.h"
#include "linear_algebra.h"
using namespace std;

namespace KARMA {
  
  class linearBeamElement: public fem {  
  public:
    void comp_Ke_loc(vector<vector<double> >& mat,const int& beamID);
    void comp_Me_loc(vector<vector<double> >& a_mat,const int& beamID);
    void setForce();
    void getElementProps();
    vector<double> thickness_;
    vector<double> poisRatio_;
    vector<double> eModulus_;
    vector<double> density_;
    vector<double> stiffness_;
    vector<double> crossSecA_;
    double defaultPoisRatio_;
    double defaultThickness_;
    double defaultEModulus_;
    double defaultDensity_;
    double defaultStiffness_;
    double defaultCrossSecA_;
    vector<int> nodeFixed_;
    vector<int> nodeFixedType_;
    int DIM = 2;
    int faceID;
    vector<vector<double> > nodeForce_;
    vector<int> forcedNodes_;
    vector<int> momentNodes_;
    int nof_nodesForces_;
    int numForcedNodes_;
    int numMomentNodes_;
    vector<vector<double> > nodeMoment_;
    
//    void compShapeFcn(vector<vector<double> >& a_Nmat,
//			 double&                  a_area,  
//			 const int&             a_beamID,
//			 const vector<double>&    a_LVec); 
//    void compDloc(vector<vector<double> >& a_D,const double& a_A,const int& a_beamID);
//    void compddPloc(vector<vector<double> >& KKl,
//		    const vector<double>&    mu,
//		    const vector<double>&    L,
//		    const double&            Delta,
//		    const vector<double>&    a,
//		    const vector<double>&    b,
//		    const vector<double>&    c,
//		    const vector<double>&    x0,
//		    const vector<double>&    x1,
//		    const vector<double>&    x2,
//		    const int&               faceID);
//    void compddPlocNew(vector<vector<double> >& a_KKl,
//		       const vector<double>&    a_mu,
//		       const vector<double>&    a_L,
//		       const double&            a_Delta,
//		       const vector<double>&    a_a,
//		       const vector<double>&    a_b,
//		       const vector<double>&    a_c,
//		       const vector<double>&          a_x0,
//		       const vector<double>&          a_x1,
//		       const vector<double>&          a_x2,
//		       const int&             a_beamID);
  };    
}

#endif   
