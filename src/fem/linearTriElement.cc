/*
 * nlfem - Nonlinear FEM Solver
 * Author  : Gokul G. Anugrah
 * Contact : gokulanugrah@gmail.com
 *
 * linearTriElement.cc — Element stiffness and mass routines.
 *   Implements local element matrices for:
 *     - MITC3 triangular shell (large deformation, Assumed Natural Strain tying)
 *     - Bernoulli-Euler beam element
 *     - Cable element
 *   Validated against benchmark MITC3 shell problems.
 */
#ifndef KARMA_LINEAR_TRI_ELEMENT
#define KARMA_LINEAR_TRI_ELEMENT
#include "linearTriElement.h"
#define SHAPEFACTOR 1.
//#define DEBUG 1
//#define DEBUG 2
//#define numericalJacobian

namespace KARMA {
  
  void linearTriElement::comp_Ke_loc(vector<double>& Ke_loc,
                                     vector<double>& Fvec,
                                     const int&      faceID,
                                     const int&      elType)
  {

    if (elType == SIMPLETRI) 
      {
        comp_Ke_TriShell(Ke_loc,Fvec,faceID);
      } 
    else if (elType == LMITC3) 
      {
        comp_Ke_MITC3(Ke_loc,Fvec,faceID);      
      } 
    else if (elType == MITC3) 
      {
        comp_NL_Ke_MITC3(Ke_loc,Fvec,faceID);      
      }
    else if (elType == BEAM) 
      {
        comp_Ke_BEBeam(Ke_loc,Fvec,faceID);
      }
    else if (elType == CABLE)
      {
        comp_Ke_Cable(Ke_loc,Fvec,faceID);
      }

  }
  
  void linearTriElement::comp_Me_loc(vector<double>& Me_loc,
                                     const int&      faceID,
                                     const int&      elType)
  {

    if (elType == SIMPLETRI)
      {
        comp_Me_TriShell(Me_loc,faceID);
      } 
    else if (elType == LMITC3)
      {
        comp_Me_MITC3(Me_loc,faceID);
      } 
    else if (elType == MITC3)
      {
        comp_Me_MITC3(Me_loc,faceID);
      }
    else if (elType == BEAM)
      {
        comp_Me_BEBeam(Me_loc,faceID);
      }
    else if (elType == CABLE)
      {
        comp_Me_Cable(Me_loc,faceID);
      }

  }
  
  void linearTriElement::comp_Ke_Cable(vector<double>& Ke_loc,
                                       vector<double>& Fvec,
                                       const int&      faceID)
  {
    //Nonlinear 3D, 2-node cable element

    int nnodes_loc   = elementData_[faceID].nnodes;
    int ndof_loc     = elementData_[faceID].ndof;
    int nnodesOffset = nnodesLocalOffset_[faceID];

    vector<int> face(nnodes_loc);
    for (int n=0;n<nnodes_loc;++n) 
      {
        int index = n + nnodesOffset;
        face[n] = faces_[index];
      }
    //face nodes
    int fn0 = face[0];
    int fn1 = face[1];

    vector<double> lastXYZ(xyz,0.);
    for (int idir=0;idir<xyz;++idir)
      {
        //current node locations
        x0_t_[idir] = nodes_ [fn0*xyz+idir];
        x1_t_[idir] = nodes_ [fn1*xyz+idir];
        //initial node locations
        x0_0_[idir] = nodes0_[fn0*xyz+idir];
        x1_0_[idir] = nodes0_[fn1*xyz+idir];
        //store element lengths in x2_t_ and x2_0_ because they are already allocated for tris
        x2_t_[idir] = x1_t_[idir] - x0_t_[idir];//current length
        x2_0_[idir] = x1_0_[idir] - x0_0_[idir];//initial length

        xhat_[idir]     = x0_0_[idir];
        xhat_[idir+xyz] = x1_0_[idir];

        uhat_[idir]     = x0_t_[idir] - x0_0_[idir];
        uhat_[idir+xyz] = x1_t_[idir] - x1_0_[idir];

        //lastXYZ[idir] = nodesM1_[fn1*xyz+idir]-nodesM1_[fn0*xyz+idir];

      }
    //compute initial element length --- Jacobian 'matrix'
    double L0 = vectorMag(x2_0_);
    //compute current element length
    defaultLength_ = vectorMag(x2_t_);

    double lastDL  = vectorMag(lastXYZ)-L0;
    //compute incremental length difference
    double deltaL  = defaultLength_ - L0;

    double emod = eModulus_[faceID];

    //calculate strains in (r,s,t)
    double erst = deltaL / L0 + 0.5 * (deltaL/L0) * (deltaL/L0);

    //displacement interpolation matrix --- store (H,r)^T * H,r
    HrTHr_[0*6+0] = 1.;
    HrTHr_[1*6+1] = 1.;
    HrTHr_[2*6+2] = 1.;
    HrTHr_[3*6+3] = 1.;
    HrTHr_[4*6+4] = 1.;
    HrTHr_[5*6+5] = 1.;

    HrTHr_[0*6+3] = -1.;
    HrTHr_[1*6+4] = -1.;
    HrTHr_[2*6+5] = -1.;
    HrTHr_[3*6+0] = -1.;
    HrTHr_[4*6+1] = -1.;
    HrTHr_[5*6+2] = -1.;

    fill(BL_.begin(),BL_.end(),0.);
    //compute BL in (r,s,t)
    //xhat contribution
    for (int i=0;i<6;++i)
      {
        for (int j=0;j<6;++j)
          {
            BL_[i] += 1./(L0*L0) * xhat_[j] * HrTHr_[j*nnodes_loc*ndof_loc + i];
          }
      }
    //uhat contribution
    for (int i=0;i<6;++i)
      {
        for (int j=0;j<6;++j)
          {
            BL_[i] += 1./(L0*L0) * uhat_[i] * HrTHr_[j*nnodes_loc*ndof_loc + i];
          }
      }
    //directional cosines to bring (r,s,t) to (x,y,z)
    double cosx = (nodes_[fn1*xyz+0] - nodes_[fn0*xyz+0])/defaultLength_;
    double cosy = (nodes_[fn1*xyz+1] - nodes_[fn0*xyz+1])/defaultLength_;
    double cosz = (nodes_[fn1*xyz+2] - nodes_[fn0*xyz+2])/defaultLength_;

    //first row
    T_[0*ndof_loc*nnodes_loc+0] =        cosx * cosx;
    T_[0*ndof_loc*nnodes_loc+1] =        cosx * cosy;
    T_[0*ndof_loc*nnodes_loc+2] =        cosx * cosz;
    T_[0*ndof_loc*nnodes_loc+3] = -1. *  cosx * cosx;
    T_[0*ndof_loc*nnodes_loc+4] = -1. *  cosx * cosy;
    T_[0*ndof_loc*nnodes_loc+5] = -1. *  cosx * cosz;

    //second row
    T_[1*ndof_loc*nnodes_loc+0] =        cosx * cosy;
    T_[1*ndof_loc*nnodes_loc+1] =        cosy * cosy;
    T_[1*ndof_loc*nnodes_loc+2] =        cosy * cosz;
    T_[1*ndof_loc*nnodes_loc+3] = -1. *  cosx * cosy;
    T_[1*ndof_loc*nnodes_loc+4] = -1. *  cosy * cosy;
    T_[1*ndof_loc*nnodes_loc+5] = -1. *  cosy * cosz;

    //third row
    T_[2*ndof_loc*nnodes_loc+0] =        cosx * cosz;
    T_[2*ndof_loc*nnodes_loc+1] =        cosy * cosz;
    T_[2*ndof_loc*nnodes_loc+2] =        cosz * cosz;
    T_[2*ndof_loc*nnodes_loc+3] = -1. *  cosx * cosz;
    T_[2*ndof_loc*nnodes_loc+4] = -1. *  cosy * cosz;
    T_[2*ndof_loc*nnodes_loc+5] = -1. *  cosz * cosz;

    //fourth row
    T_[3*ndof_loc*nnodes_loc+0] = -1. *  cosx * cosx;
    T_[3*ndof_loc*nnodes_loc+1] = -1. *  cosx * cosy;
    T_[3*ndof_loc*nnodes_loc+2] = -1. *  cosx * cosz;
    T_[3*ndof_loc*nnodes_loc+3] =        cosx * cosx;
    T_[3*ndof_loc*nnodes_loc+4] =        cosx * cosy;
    T_[3*ndof_loc*nnodes_loc+5] =        cosx * cosz;

    //fifth row
    T_[4*ndof_loc*nnodes_loc+0] = -1. *  cosx * cosy;
    T_[4*ndof_loc*nnodes_loc+1] = -1. *  cosy * cosy;
    T_[4*ndof_loc*nnodes_loc+2] = -1. *  cosy * cosz;
    T_[4*ndof_loc*nnodes_loc+3] =        cosx * cosy;
    T_[4*ndof_loc*nnodes_loc+4] =        cosy * cosy;
    T_[4*ndof_loc*nnodes_loc+5] =        cosy * cosz;

    //sixth row
    T_[5*ndof_loc*nnodes_loc+0] = -1. *  cosx * cosz;
    T_[5*ndof_loc*nnodes_loc+1] = -1. *  cosy * cosz;
    T_[5*ndof_loc*nnodes_loc+2] = -1. *  cosz * cosz;
    T_[5*ndof_loc*nnodes_loc+3] =        cosx * cosz;
    T_[5*ndof_loc*nnodes_loc+4] =        cosy * cosz;
    T_[5*ndof_loc*nnodes_loc+5] =        cosz * cosz;

    vector<double> btb(36,0.);
    for (int i=0;i<6;++i)
      {
        for (int j=0;j<6;++j)
          {
            btb[i*6+j] = BL_[j] * BL_[j] * emod;
          }
      }

    //compute Ke_loc in (x,y,z)
    for (int i=0;i<6;++i)
      {
        for (int j=0;j<6;++j)
          {
            for (int k=0;k<6;++k)
              {
                Ke_loc[i*6+j] += T_[k*6+i] *  btb[i*6+k] * T_[k*6+j] * L0 * Area_[faceID];
              }
          }
      }

    //calculate stresses in (r,s,t)
    double Srst = emod * max(erst,0.);

    myCableTensionVector_[faceID] = Srst * Area_[faceID];

    //compute work due to internal forces in (r,s,t)
    vector<double> ftmp(6,0.);
    for (int i=0;i<6;++i) ftmp[i] = BL_[i] * Srst * L0 * Area_[faceID];

    //bring internal work to (x,y,z)
    for (int i=0;i<6;++i) for (int j=0;j<6;++j) Fvec[i] += T_[j*6+i] * ftmp[j];

  }
    
  void linearTriElement::comp_Me_Cable(vector<double>& Me_loc,
                                       const int&      faceID)
  {
    //compute mass matrix for 2-node cable element

    int nnodes_loc   = elementData_[faceID].nnodes;
    int ndof_loc     = elementData_[faceID].ndof;
    int nnodesOffset = nnodesLocalOffset_[faceID];

    //lumped mass matrix
    vector<int> face(nnodes_loc);
    for (int n=0;n<nnodes_loc;++n) 
      {
        int index = n + nnodesOffset;
        face[n] = faces_[index];
      }
    //face nodes
    int fn0 = face[0];
    int fn1 = face[1];
    vector<double> x0(DIM,0.);
    vector<double> x1(DIM,0.);

    for (int idir=0;idir<xyz;++idir) 
      {
        x0_0_[idir] = nodes0_[fn0*xyz+idir];
        x1_0_[idir] = nodes0_[fn1*xyz+idir];
        x2_0_[idir] = x1_0_[idir] - x0_0_[idir];//initial length
        x2_t_[idir] = x1_t_[idir] - x0_t_[idir];//current length
      }

    double L0 = vectorMag(x2_0_);
    defaultLength_ = vectorMag(x2_t_);
    //precompute and store H^T * H evaluated at r = 1/2 in HrTHr_
    HrTHr_[0*6+0] = 0.0625;
    HrTHr_[1*6+1] = 0.0625;
    HrTHr_[2*6+2] = 0.0625;
    HrTHr_[3*6+3] = 0.5625;
    HrTHr_[4*6+4] = 0.5625;
    HrTHr_[5*6+5] = 0.5625;

    HrTHr_[0*6+3] = 0.1875;
    HrTHr_[1*6+4] = 0.1875;
    HrTHr_[2*6+5] = 0.1875;
    HrTHr_[3*6+0] = 0.1875;
    HrTHr_[4*6+1] = 0.1875;
    HrTHr_[5*6+2] = 0.1875;

    //compute transformation matrix
    //directional cosines to bring (r,s,t) to (x,y,z)
    double cosx = (nodes_[fn1*xyz+0] - nodes_[fn0*xyz+0])/defaultLength_;
    double cosy = (nodes_[fn1*xyz+1] - nodes_[fn0*xyz+1])/defaultLength_;
    double cosz = (nodes_[fn1*xyz+2] - nodes_[fn0*xyz+2])/defaultLength_;

    //first row
    T_[0*ndof_loc*nnodes_loc+0] =        cosx * cosx;
    T_[0*ndof_loc*nnodes_loc+1] =        cosx * cosy;
    T_[0*ndof_loc*nnodes_loc+2] =        cosx * cosz;
    T_[0*ndof_loc*nnodes_loc+3] = -1. *  cosx * cosx;
    T_[0*ndof_loc*nnodes_loc+4] = -1. *  cosx * cosy;
    T_[0*ndof_loc*nnodes_loc+5] = -1. *  cosx * cosz;

    //second row
    T_[1*ndof_loc*nnodes_loc+0] =        cosx * cosy;
    T_[1*ndof_loc*nnodes_loc+1] =        cosy * cosy;
    T_[1*ndof_loc*nnodes_loc+2] =        cosy * cosz;
    T_[1*ndof_loc*nnodes_loc+3] = -1. *  cosx * cosy;
    T_[1*ndof_loc*nnodes_loc+4] = -1. *  cosy * cosy;
    T_[1*ndof_loc*nnodes_loc+5] = -1. *  cosy * cosz;

    //third row
    T_[2*ndof_loc*nnodes_loc+0] =        cosx * cosz;
    T_[2*ndof_loc*nnodes_loc+1] =        cosy * cosz;
    T_[2*ndof_loc*nnodes_loc+2] =        cosz * cosz;
    T_[2*ndof_loc*nnodes_loc+3] = -1. *  cosx * cosz;
    T_[2*ndof_loc*nnodes_loc+4] = -1. *  cosy * cosz;
    T_[2*ndof_loc*nnodes_loc+5] = -1. *  cosz * cosz;

    //fourth row
    T_[3*ndof_loc*nnodes_loc+0] = -1. *  cosx * cosx;
    T_[3*ndof_loc*nnodes_loc+1] = -1. *  cosx * cosy;
    T_[3*ndof_loc*nnodes_loc+2] = -1. *  cosx * cosz;
    T_[3*ndof_loc*nnodes_loc+3] =        cosx * cosx;
    T_[3*ndof_loc*nnodes_loc+4] =        cosx * cosy;
    T_[3*ndof_loc*nnodes_loc+5] =        cosx * cosz;

    //fifth row
    T_[4*ndof_loc*nnodes_loc+0] = -1. *  cosx * cosy;
    T_[4*ndof_loc*nnodes_loc+1] = -1. *  cosy * cosy;
    T_[4*ndof_loc*nnodes_loc+2] = -1. *  cosy * cosz;
    T_[4*ndof_loc*nnodes_loc+3] =        cosx * cosy;
    T_[4*ndof_loc*nnodes_loc+4] =        cosy * cosy;
    T_[4*ndof_loc*nnodes_loc+5] =        cosy * cosz;

    //sixth row
    T_[5*ndof_loc*nnodes_loc+0] = -1. *  cosx * cosz;
    T_[5*ndof_loc*nnodes_loc+1] = -1. *  cosy * cosz;
    T_[5*ndof_loc*nnodes_loc+2] = -1. *  cosz * cosz;
    T_[5*ndof_loc*nnodes_loc+3] =        cosx * cosz;
    T_[5*ndof_loc*nnodes_loc+4] =        cosy * cosz;
    T_[5*ndof_loc*nnodes_loc+5] =        cosz * cosz;

    vector<double> me_tmp(36,0.);
    for (int i=0;i<6;++i) Me_loc[i*6+i] = 0.5 * density_[faceID] * Area_[faceID] * L0;

  }
  
  void linearTriElement::comp_Ke_BEBeam(vector<double>& Ke_loc,
                                        vector<double>& Fvec,
                                        const int&      faceID)
  {
    //2D Bernoulli-Euler beam element stiffness matrix

    int nnodes_loc = elementData_[faceID].nnodes;
    int ndof_loc   = elementData_[faceID].ndof;

    vector<int> face(nnodes_loc);
    for (int n=0;n<nnodes_loc;++n) 
      {
        int index = n + nnodes_loc*faceID;
        face[n] = faces_[index];
      }

    //face nodes
    int fn0=face[0];
    int fn1=face[1];

    //node coordinates
    vector<double> x0(DIM,0.);
    vector<double> x1(DIM,0.);
    for (int idir=0;idir<DIM;++idir) 
      {
        x0[idir] = nodes_[fn0*xyz+idir];
        x1[idir] = nodes_[fn1*xyz+idir];
      }
    //calculate length of element 
    defaultLength_ = sqrt(pow((abs(x1[0]-x0[0])),2) + pow((abs(x1[1]-x0[1])),2));
    //analytical beam stiffness matrix
    Ke_loc[0*ndof_loc*nnodes_loc+0] = eModulus_[faceID]*MoI_[faceID]/(defaultLength_*defaultLength_*defaultLength_)*12;
    Ke_loc[0*ndof_loc*nnodes_loc+1] = eModulus_[faceID]*MoI_[faceID]/(defaultLength_*defaultLength_*defaultLength_)*6*defaultLength_;
    Ke_loc[0*ndof_loc*nnodes_loc+2] = eModulus_[faceID]*MoI_[faceID]/(defaultLength_*defaultLength_*defaultLength_)*-12;
    Ke_loc[0*ndof_loc*nnodes_loc+3] = eModulus_[faceID]*MoI_[faceID]/(defaultLength_*defaultLength_*defaultLength_)*6*defaultLength_;
    Ke_loc[1*ndof_loc*nnodes_loc+0] = Ke_loc[0*ndof_loc*nnodes_loc+1];
    Ke_loc[1*ndof_loc*nnodes_loc+1] = eModulus_[faceID]*MoI_[faceID]/(defaultLength_*defaultLength_*defaultLength_)*4*defaultLength_*defaultLength_;
    Ke_loc[1*ndof_loc*nnodes_loc+2] = eModulus_[faceID]*MoI_[faceID]/(defaultLength_*defaultLength_*defaultLength_)*-6*defaultLength_;
    Ke_loc[1*ndof_loc*nnodes_loc+3] = eModulus_[faceID]*MoI_[faceID]/(defaultLength_*defaultLength_*defaultLength_)*2*defaultLength_*defaultLength_;
    Ke_loc[2*ndof_loc*nnodes_loc+0] = Ke_loc[0*ndof_loc*nnodes_loc+2];
    Ke_loc[2*ndof_loc*nnodes_loc+1] = Ke_loc[1*ndof_loc*nnodes_loc+2];
    Ke_loc[2*ndof_loc*nnodes_loc+2] = eModulus_[faceID]*MoI_[faceID]/(defaultLength_*defaultLength_*defaultLength_)*12;
    Ke_loc[2*ndof_loc*nnodes_loc+3] = eModulus_[faceID]*MoI_[faceID]/(defaultLength_*defaultLength_*defaultLength_)*-6*defaultLength_;
    Ke_loc[3*ndof_loc*nnodes_loc+0] = Ke_loc[0*ndof_loc*nnodes_loc+3];
    Ke_loc[3*ndof_loc*nnodes_loc+1] = Ke_loc[1*ndof_loc*nnodes_loc+3];
    Ke_loc[3*ndof_loc*nnodes_loc+2] = Ke_loc[2*ndof_loc*nnodes_loc+3];
    Ke_loc[3*ndof_loc*nnodes_loc+3] = eModulus_[faceID]*MoI_[faceID]/(defaultLength_*defaultLength_*defaultLength_)*4*defaultLength_*defaultLength_;

    //dumpMat(nodes_,"nodes");
    //dumpMat((*faces_),"faces");
    //dumpMat(Ke_loc,"Ke_loc");
    //wait();

  }

  void linearTriElement::comp_Me_BEBeam(vector<double>& Me_loc,
                                        const int&      faceID)
  {
    //compute local mass matrix for 2D Bernoulli-Euler beam

    int nnodes_loc = elementData_[faceID].nnodes;
    int ndof_loc   = elementData_[faceID].ndof;

    //get face ID
    vector<int> face(nnodes_loc);
    for (int n=0;n<nnodes_loc;++n) 
      {
        int index = n + nnodes_loc*faceID;
        face[n] = faces_[index];
      }
    //face nodes
    int fn0 = face[0];
    int fn1 = face[1];
    vector<double> x0(DIM,0.);
    vector<double> x1(DIM,0.);

    for (int idir=0;idir<DIM;++idir) 
      {
        x0[idir] = nodes_[fn0*xyz+idir];
        x1[idir] = nodes_[fn1*xyz+idir];      
      }
    //calculate length of element 
    defaultLength_ = sqrt(pow((abs(x1[0]-x0[0])),2) + pow((abs(x1[1]-x0[1])),2));

    //analytical form for mass matrix
    if (timeIntegrationMethod_ == "implicit")
      {
        //consistent mass matrix
        double lead = density_[faceID]*defaultLength_*beamCSA_/420.;        
        Me_loc[0*ndof_loc*nnodes_loc+0] = lead*156.;
        Me_loc[0*ndof_loc*nnodes_loc+1] = lead*22.*defaultLength_;
        Me_loc[0*ndof_loc*nnodes_loc+2] = lead*54.;
        Me_loc[0*ndof_loc*nnodes_loc+3] = lead*-13.*defaultLength_;
        Me_loc[1*ndof_loc*nnodes_loc+0] = Me_loc[0*ndof_loc*nnodes_loc+1];
        Me_loc[1*ndof_loc*nnodes_loc+1] = lead*4.*defaultLength_*defaultLength_;
        Me_loc[1*ndof_loc*nnodes_loc+2] = lead*13.*defaultLength_;
        Me_loc[1*ndof_loc*nnodes_loc+3] = lead*-3.*defaultLength_*defaultLength_;
        Me_loc[2*ndof_loc*nnodes_loc+0] = Me_loc[0*ndof_loc*nnodes_loc+2];
        Me_loc[2*ndof_loc*nnodes_loc+1] = Me_loc[1*ndof_loc*nnodes_loc+2];
        Me_loc[2*ndof_loc*nnodes_loc+2] = lead*156.;
        Me_loc[2*ndof_loc*nnodes_loc+3] = lead*-22.*defaultLength_;
        Me_loc[3*ndof_loc*nnodes_loc+0] = Me_loc[0*ndof_loc*nnodes_loc+3];
        Me_loc[3*ndof_loc*nnodes_loc+1] = Me_loc[1*ndof_loc*nnodes_loc+3];
        Me_loc[3*ndof_loc*nnodes_loc+2] = Me_loc[2*ndof_loc*nnodes_loc+3];
        Me_loc[3*ndof_loc*nnodes_loc+3] = lead*4.*defaultLength_*defaultLength_;

        //normal lumped mass matrix
        //Me_loc[0][0] = 0.5*(beamCSA_ * defaultDensity_ * defaultLength_);
        //Me_loc[1][1] = 1.E-4;
        //Me_loc[2][2] = 0.5*(beamCSA_ * defaultDensity_ * defaultLength_);
        //Me_loc[3][3] = 1.E-4;
      }
    else
      {
        //pre-inverted lumped mass matrix
        if (faceID == 0)
          { 
            Me_loc[0*ndof_loc*nnodes_loc+0] = 1.  / (beamCSA_ * density_[faceID] * defaultLength_);
            Me_loc[1*ndof_loc*nnodes_loc+1] = 1.  / 1.E-6;
            Me_loc[2*ndof_loc*nnodes_loc+2] = 0.5 / (beamCSA_ * density_[faceID] * defaultLength_);
            Me_loc[3*ndof_loc*nnodes_loc+3] = 0.5 / 1.E-6;
         }
        else if (faceID == nof_faces_glob_-1)
          {
            Me_loc[0*ndof_loc*nnodes_loc+0] = 0.5 / (beamCSA_ * density_[faceID] * defaultLength_);
            Me_loc[1*ndof_loc*nnodes_loc+1] = 0.5 / 1.E-6;
            Me_loc[2*ndof_loc*nnodes_loc+2] = 1.  / (beamCSA_ * density_[faceID] * defaultLength_);
            Me_loc[3*ndof_loc*nnodes_loc+3] = 1.  / 1.E-6;
          }
        else
          {
            Me_loc[0*ndof_loc*nnodes_loc+0] = 0.5 / (beamCSA_ * density_[faceID] * defaultLength_);
            Me_loc[1*ndof_loc*nnodes_loc+1] = 0.5 / 1.E-6;
            Me_loc[2*ndof_loc*nnodes_loc+2] = 0.5 / (beamCSA_ * density_[faceID] * defaultLength_);
            Me_loc[3*ndof_loc*nnodes_loc+3] = 0.5 / 1.E-6;
          }
        //dumpMat(Me_loc,"Me_loc");
      }

  }

  void linearTriElement::comp_Ke_TriShell(vector<double>& Ke_loc,
                                          vector<double>& Fvec,
                                          const int&      faceID)
  {

    int nnodes_loc = elementData_[faceID].nnodes;
    int ndof_loc   = elementData_[faceID].ndof;

    vector<int> face(nnodes_loc);
    for (int n=0;n<nnodes_loc;++n)
      {
        int index = n+nnodes_loc*faceID;
        face[n] = faces_[index];
      }

    // face nodes
    int fn0=face[0];
    int fn1=face[1];
    int fn2=face[2];

    // node coordinates
    vector<double> x0(DIM,0.);      
    vector<double> x1(DIM,0.);      
    vector<double> x2(DIM,0.);
    for (int idir=0;idir<DIM;++idir) {
      x0[idir] = nodes_[fn0*xyz+idir];
      x1[idir] = nodes_[fn1*xyz+idir];
      x2[idir] = nodes_[fn2*xyz+idir];
    }      
    
    //determine transformation matrix from three node locations
    vector<double> dummyVec(3,0.);    
    vector<vector<double> > transMat(3,dummyVec);    
    transMat=compTransMat(x0,x1,x2); 

    vector<vector<double> > transNodes(3,dummyVec);

    for (int d=0;d<DIM;++d) {
      int fn=face[d]; //get face node
      vector<double> xyzCoords(DIM,0.);
      for (int dir=0;dir<DIM;++dir) {
        xyzCoords[dir]=nodes_[fn*xyz+dir];
      }

      mTv_mult(transMat,xyzCoords,transNodes[d]);
    }

    vector<vector<double> > nodes(DIM,dummyVec);
    for (int d=0;d<DIM;++d) {
      nodes[1][d] = x1[d] - x0[d];
      nodes[2][d] = x2[d] - x0[d];
    }

    mmT_mult(nodes,transMat,transNodes);

    vector<double> a (DIM,0.);
    vector<double> b (DIM,0.);
    vector<double> c (DIM,0.);
    vector<double> l (DIM,0.);

    for (int d=0;d<DIM;++d) {

      int dir1=dirIJK_[d][0];
      int dir2=dirIJK_[d][1];

      vector<double> xx1 = transNodes[dir1];      
      vector<double> xx2 = transNodes[dir2];      

      a[d]=xx1[0]*xx2[1]-xx2[0]*xx1[1];
      b[d]=xx1[1]-xx2[1];
      c[d]=xx2[0]-xx1[0];
      l[d]=sqrt(c[d]*c[d]+b[d]*b[d]+(xx1[2]-xx2[2])*(xx1[2]-xx2[2]));

    } 
    vector<double> mu(DIM,0.);

    for (int d=0;d<DIM;++d) {
      int dir1=dirIJK_[d][0];
      int dir2=dirIJK_[d][1];      
      mu[d]=(l[dir2]*l[dir2]-l[dir1]*l[dir1])/(l[d]*l[d]);
    }  
    double Delta=0.5*(b[0]*c[1]-b[1]*c[0]);
    // - shape functions
    //         Nij is the shape function for node i dof j
    //         Ni is the cst shapefunction for node i
    // Shape function
    //N = [0  0  0   0   0   0  0  0   0   0   0  0  0   0   0   ;
    //     0  0  0   0   0   0  0  0   0   0   0  0  0   0   0   ;
    //     0  0  N11 N12 N13 0  0  N21 N22 N23 0  0  N31 N32 N33 ];

    double N1 = 1.;
    double N2 = 1.;
    double N3 = 1.;
    double N11 = 1.;
    double N12 = 1.;
    double N13 = 1.;
    double N21 = 1.;
    double N22 = 1.;
    double N23 = 1.;
    double N31 = 1.;
    double N32 = 1.;
    double N33 = 1.;

    vector<double> dummyVec2(15,0.);
    vector<vector<double> > N(DIM,dummyVec2);
    N[2][2 ]=N11;
    N[2][3 ]=N12;
    N[2][4 ]=N13;
    N[2][7 ]=N21;
    N[2][8 ]=N22;
    N[2][9 ]=N23;
    N[2][12]=N31;
    N[2][13]=N32;
    N[2][14]=N33;

    vector<double> L(3,1./3.);
    vector<double> dummyVec3_0pt5(3,0.5);
    vector<vector<double> > LVec(3,dummyVec3_0pt5);
    LVec[0][2]=0.;
    LVec[1][0]=0.;
    LVec[2][1]=0.;

    vector<double> weights(1,1.);
    int noWeights=weights.size();
    if (noWeights==1) {
      LVec[0][0]=1./3.;
      LVec[0][1]=1./3.;
      LVec[0][2]=1./3.;
    }

    vector<double> dummyVec3(18,0.);
    vector<vector<double> > KK(18,dummyVec3);
    for (int k=0;k<noWeights;++k)  {
      vector<vector<double> > KKloc(18,dummyVec3);
      compddPloc(KKloc,mu,LVec[k],Delta,a,b,c,transNodes[0],transNodes[1],transNodes[2],faceID);
      for (int r=0;r<18;++r) {
        for (int c=0;c<18;++c) {
          //KK[r][c]+=weights[k]*KKloc[r][c];
        }
      }
    }

    ////dumpMat(KK,"KK LTE");
    //int pop;
    //cin >> pop;

    vector<vector<double> > TT(18,dummyVec3);
    for (int n=0;n<3;++n) {
      for (int i=0;i<3;++i) {
        for (int j=0;j<3;++j) {
          int iloc=i+n*ndof_loc;
          int jloc=j+n*ndof_loc;
          TT[iloc][jloc]=transMat[i][j];
        }
      }
      for (int i=0;i<3;++i) {
        for (int j=0;j<3;++j) {
          int iloc=i+n*ndof_loc+3;
          int jloc=j+n*ndof_loc+3;
          TT[iloc][jloc]=transMat[i][j];
        }
      }
    }

    //aTma_mult(TT,KK,Ke_loc);
    ////dumpMat(transMat,"Transmat");
    ////dumpMat(Ke_loc,"Ke_loc LinearTri");
    //int stiffer;
    //cin >> stiffer;
  }

  void linearTriElement::comp_Me_TriShell(vector<double>& Me_loc,
                                          const int&      faceID)
  {
    //computing local mass matrix

    int nnodes_loc = elementData_[faceID].nnodes;
    int ndof_loc   = elementData_[faceID].ndof;

    // get face Id
    vector<int> face(nnodes_loc);
    for (int n=0;n<nnodes_loc;++n) 
      {
        int index = n + nnodes_loc*faceID;
        face[n] = faces_[index];
      }
    // face nodes
    int fn0=face[0];
    int fn1=face[1];
    int fn2=face[2];
    vector<double> x0(DIM,0.);
    vector<double> x1(DIM,0.);
    vector<double> x2(DIM,0.);

    for (int idir=0;idir<DIM;++idir) {
      x0[idir] = nodes_[fn0*xyz+idir];
      x1[idir] = nodes_[fn1*xyz+idir];      
      x2[idir] = nodes_[fn2*xyz+idir]; 
    }

    // obtain transformation matrix from local nodes
    vector<double> dummyVec3(3,0.5);
    vector<vector<double> > transMat(3,dummyVec3);
    transMat=compTransMat(x0,x1,x2);  
    vector<vector<double> > LVec(3,dummyVec3);
    LVec[0][2]=0.;
    LVec[1][0]=0.;
    LVec[2][1]=0.;
    vector<double> dummyVec2(15,0.);
    vector<double> dummyVec4(18,0.);
    vector<vector<double> > Nmatloc( 3,dummyVec2);
    vector<vector<double> > Nmat   (15,dummyVec2);
    vector<vector<double> > Mat    (18,dummyVec4);

    //    vector<double> weights(3,1./3.);
    vector<double> weights(1,1.);
    int noWeights = weights.size();
    if (noWeights == 1) { 
      LVec[0][0] = 1./3;
      LVec[0][1] = 1./3;
      LVec[0][2] = 1./3;
    }

    double area=0.;
    for (int k=0;k<noWeights;++k) {
      area=0.;
      compShapeFunction(Nmatloc,area,faceID,LVec[k]);        

      mTm_mult(Nmatloc,Nmatloc,Nmat); 

      double factor=area*thickness_[faceID]*density_[faceID];
      for (int r=0;r<3;++r) {
        for (int c=0;c<3;++c) {
          for (int rr=0;rr<5;++rr) {
            for (int cc=0;cc<5;++cc) {
              Mat[rr+r*6][cc+c*6]+=weights[k]*Nmat[rr+r*5][cc+c*5]*factor;
            }
          }
        }
      }
    }

    vector<double> dumVec(18,0.);      
    vector<vector<double> > TT(18,dumVec);
    for (int n=0;n<3;++n) {
      for (int i=0;i<3;++i) {
        for (int j=0;j<3;++j) {
          int iloc=i+n*ndof_loc;
          int jloc=j+n*ndof_loc;
          TT[iloc][jloc]=transMat[i][j];
        }
      }
      for (int i=0;i<3;++i) {
        for (int j=0;j<3;++j) {
          int iloc=i+n*ndof_loc+3;
          int jloc=j+n*ndof_loc+3;
          TT[iloc][jloc]=transMat[i][j]; 
        }
      }
    }

    //aTma_mult(TT,Mat,Me_loc);
  }
  
  void linearTriElement::compShapeFunction(vector<vector<double> >& Nmat,
                                           double&                  area,  
                                           const int&               faceID,
                                           const vector<double>&    LVec)
  {
    //comp local shape function

    int nnodes_loc = elementData_[faceID].nnodes;
    int ndof_loc   = elementData_[faceID].ndof;
    vector<int> face(nnodes_loc,0);
    for (int n=0;n<nnodes_loc;++n) 
      {
        int index = n + nnodes_loc*faceID;
        face[n] = faces_[index];
      }

    int fn0=face[0];
    int fn1=face[1];
    int fn2=face[2];

    vector<double> x0(3,0.);      
    vector<double> x1(3,0.);      
    vector<double> x2(3,0.);

    for (int idir=0;idir<DIM;++idir) 
      {
        x0[idir] = nodes0_[fn0*xyz+idir];
        x1[idir] = nodes0_[fn1*xyz+idir];
        x2[idir] = nodes0_[fn2*xyz+idir];
      }      
    vector<double> dummyVec3(3,0.);    
    vector<vector<double> > transMat(3,dummyVec3);
    transMat=compTransMat(x0,x1,x2); 

    vector<vector<double> > transNodes(3,dummyVec3);

    for (int d=0;d<DIM;++d) {
      int fn=face[d];
      vector<double> xyzvec(DIM,0.);
      for (int idir=0;idir<DIM;++idir) xyzvec[idir] = nodes0_[fn*xyz+idir];

      mTv_mult(transMat,xyzvec,transNodes[d]);  
    }

    vector<double> a(3,0.);
    vector<double> b(3,0.);
    vector<double> c(3,0.);
    vector<double> l(3,0.);

    for (int d=0; d<DIM; ++d) {

      int dir1=dirIJK_[d][0];
      int dir2=dirIJK_[d][1];
      vector<double> xx1 = transNodes[dir1];      
      vector<double> xx2 = transNodes[dir2];      

      a[d]=xx1[0]*xx2[1]-xx2[0]*xx1[1];
      b[d]=xx1[1]-xx2[1];
      c[d]=xx2[0]-xx1[0];

      l[d]=sqrt(c[d]*c[d]+b[d]*b[d]+(xx1[2]-xx2[2])*(xx1[2]-xx2[2]));    
    }

    vector<double> mu(3,0.);  
    for (int d=0;d<DIM;++d) {
      int dir1=dirIJK_[d][0];
      int dir2=dirIJK_[d][1];

      mu[d]=(l[dir2]*l[dir2]-l[dir1]*l[dir1])/(l[d]*l[d]);
    }

    double Delta=0.5*(b[0]*c[1]-b[1]*c[0]);
    area=Delta;
    vector<double> Pvec(9,0.);
    double L1=LVec[0];
    Pvec[0]=L1;
    double L2=LVec[1];
    Pvec[1]=L2;
    double L3=LVec[2];
    Pvec[2]=L3;
    double L1L2=L1*L2;
    Pvec[3]=L1L2;
    double L2L3=L2*L3;
    Pvec[4]=L2L3;
    double L3L1=L3*L1;
    Pvec[5]=L3L1;
    double L1L2L3=L1*L2*L3;
    Pvec[6]=L1*L1L2+0.5*L1L2L3*(3.*(1.-mu[2])*L1-(1.+3.*mu[2])*L2+(1.+3.*mu[2])*L3);
    Pvec[7]=L2*L2L3+0.5*L1L2L3*(3.*(1.-mu[0])*L2-(1.+3.*mu[0])*L3+(1.+3.*mu[0])*L1);
    Pvec[8]=L3*L3L1+0.5*L1L2L3*(3.*(1.-mu[1])*L3-(1.+3.*mu[1])*L1+(1.+3.*mu[1])*L2);

    // - shape functions
    //         Nij is the shape function for node i dof j
    //         Ni  is the cst shapefunction for node i

    double N1 = L1;
    double N2 = L2;
    double N3 = L3;
    double N11 = Pvec[0]-Pvec[3]+Pvec[5]+2.*(Pvec[6]-Pvec[8]);
    double N12 = -b[1]*(Pvec[8]-Pvec[5])-b[2]*Pvec[6];
    double N13 = -c[1]*(Pvec[8]-Pvec[5])-c[2]*Pvec[6];
    double N21 = Pvec[1]-Pvec[4]+Pvec[3]+2.*(Pvec[7]-Pvec[6]);
    double N22 = -b[2]*(Pvec[6]-Pvec[3])-b[0]*Pvec[7];
    double N23 = -c[2]*(Pvec[6]-Pvec[3])-c[0]*Pvec[7];
    double N31 = Pvec[2]-Pvec[5]+Pvec[4]+2.*(Pvec[8]-Pvec[7]);
    double N32 = -b[0]*(Pvec[7]-Pvec[4])-b[1]*Pvec[8];
    double N33 = -c[0]*(Pvec[7]-Pvec[4])-c[1]*Pvec[8];
    Nmat[0][ 0]=N1;
    Nmat[1][ 1]=N1;
    Nmat[0][ 5]=N2;
    Nmat[1][ 6]=N2;
    Nmat[0][10]=N3;
    Nmat[1][11]=N3;
    Nmat[2][ 2]=N11;
    Nmat[2][ 3]=N12;
    Nmat[2][ 4]=N13;
    Nmat[2][ 7]=N21;
    Nmat[2][ 8]=N22;
    Nmat[2][ 9]=N23;
    Nmat[2][12]=N31;
    Nmat[2][13]=N32;
    Nmat[2][14]=N33;        
  }
  
  void linearTriElement::compDloc(vector<vector<double> >& Dmat,
                                  const double&            A,
                                  const int&               faceID)
  {

    int nnodes_loc = elementData_[faceID].nnodes;
    int ndof_loc   = elementData_[faceID].ndof;

    Dmat.resize(6);
    for (int i=0;i<6;++i) {      
      Dmat[i].resize(6);
      for (int j=0;j<6;++j) {
	    Dmat[i][j]=0.;
	  }
    }

    double I=A*thickness_[faceID]*eModulus_[faceID]/(1.-poisRatio_[faceID]*poisRatio_[faceID]);
    //pow(thickness_[faceID],3)/12.;
    //in-plane
    Dmat[0][0]=I;
    Dmat[0][1]=I*poisRatio_[faceID];
    Dmat[1][1]=I;
    Dmat[1][0]=I*poisRatio_[faceID];
    Dmat[2][2]=I*(1.-poisRatio_[faceID])*0.5;
    //out-of-plane
    I*=pow(thickness_[faceID],2)/(12.);
    Dmat[3][3]=I;
    Dmat[3][4]=I*poisRatio_[faceID];
    Dmat[4][4]=I;
    Dmat[4][3]=I*poisRatio_[faceID];
    Dmat[5][5]=I*(1.-poisRatio_[faceID])*0.5;
  }
 
  void linearTriElement::compddPloc(vector<vector<double> >& KKl,
                                    const vector<double>&    mu,
                                    const vector<double>&    L,
                                    const double&            Delta,
                                    const vector<double>&    a,
                                    const vector<double>&    b,
                                    const vector<double>&    c,
                                    const vector<double>&    x0,
                                    const vector<double>&    x1,
                                    const vector<double>&    x2,
                                    const int&               faceID)
  {

    int nnodes_loc = elementData_[faceID].nnodes;
    int ndof_loc   = elementData_[faceID].ndof;

    vector<double> dummy6(6,0.);
    vector<vector<double> > ddP(9,dummy6);

    ddP[6][0]=2*L[1]+L[1]*L[2]*3*(1-mu[2]) ; 
    ddP[7][0]=L[1]*L[2]*(1+3*mu[0]) ; 
    ddP[8][0]=-L[1]*L[2]*(1+3*mu[1]);
    // d^2P/dL2^2
    ddP[6][1]=-L[0]*L[2]*(1+3*mu[2]) ; 
    ddP[7][1]=2*L[2]+L[0]*L[2]*3*(1-mu[0]) ; 
    ddP[8][1]=L[0]*L[2]*(1+3*mu[1]);
    // d^2P/dL3^2
    ddP[6][2]=L[0]*L[1]*(1+3*mu[2]) ; 
    ddP[7][2]=-L[0]*L[1]*(1+3*mu[0]) ; 
    ddP[8][2]=2*L[0]+L[0]*L[1]*3*(1-mu[1]);
    // d^2P/dL1dL2
    ddP[3][3]=1.;
    ddP[6][3]=2*L[0]+L[0]*L[2]*3*(1-mu[2])-L[1]*L[2]*(1+3*mu[2])+0.5*L[2]*L[2]*(1+3*mu[2]) ; 
    ddP[7][3]=L[1]*L[2]*3*(1-mu[0])-0.5*L[2]*L[2]*(1+3*mu[0])+L[0]*L[2]*(1+3*mu[0]); 
    ddP[8][3]=0.5*L[2]*L[2]*3*(1-mu[1])-L[0]*L[2]*(1+3*mu[1])+L[1]*L[2]*(1+3*mu[1]);
    // d^2P/dL1dL3
    ddP[5][4]=1.;
    ddP[6][4]=L[0]*L[1]*3*(1-mu[2])-0.5*L[1]*L[1]*(1+3*mu[2])+L[1]*L[2]*(1+3*mu[2]) ; 
    ddP[7][4]=0.5*L[1]*L[1]*3*(1-mu[0])-L[1]*L[2]*(1+3*mu[0])+L[0]*L[1]*(1+3*mu[0]); 
    ddP[8][4]=2*L[2]+L[1]*L[2]*3*(1-mu[1])-L[0]*L[1]*(1+3*mu[1])+0.5*L[1]*L[1]*(1+3*mu[1]);
    // d^2P/dL2dL3
    ddP[4][5]=1.;
    ddP[6][5]=0.5*L[0]*L[0]*3*(1-mu[2])-L[0]*L[1]*(1+3*mu[2])+L[0]*L[2]*(1+3*mu[2]) ; 
    ddP[7][5]=2*L[1]+L[0]*L[1]*3*(1-mu[0])-L[0]*L[2]*(1+3*mu[0])+0.5*L[0]*L[0]*(1+3*mu[0]); 
    ddP[8][5]=L[0]*L[2]*3*(1-mu[1])-0.5*L[0]*L[0]*(1+3*mu[1])+L[0]*L[1]*(1+3*mu[1]);

    vector<double> dummyVec3(3,0.5);
    vector<vector<double> > ddN11(3,dummyVec3);
    ddN11[0][0]=ddP[0][0]-ddP[3][0]+ddP[5][0]+2*(ddP[6][0]-ddP[8][0]); 
    ddN11[0][1]=ddP[0][3]-ddP[3][3]+ddP[5][3]+2*(ddP[6][3]-ddP[8][3]);
    ddN11[0][2]=ddP[0][4]-ddP[3][4]+ddP[5][4]+2*(ddP[6][4]-ddP[8][4]);

    ddN11[1][0]=ddP[0][3]-ddP[3][3]+ddP[5][3]+2*(ddP[6][3]-ddP[8][3]); 
    ddN11[1][1]=ddP[0][1]-ddP[3][1]+ddP[5][1]+2*(ddP[6][1]-ddP[8][1]);  
    ddN11[1][2]=ddP[0][5]-ddP[3][5]+ddP[5][5]+2*(ddP[6][5]-ddP[8][5]);

    ddN11[2][0]=ddP[0][4]-ddP[3][4]+ddP[5][4]+2*(ddP[6][4]-ddP[8][4]); 
    ddN11[2][1]=ddP[0][5]-ddP[3][5]+ddP[5][5]+2*(ddP[6][5]-ddP[8][5]); 
    ddN11[2][2]=ddP[0][2]-ddP[3][2]+ddP[5][2]+2*(ddP[6][2]-ddP[8][2]);

    vector<vector<double> > ddN12(3,dummyVec3);
    ddN12[0][0]=-b[1]*(ddP[8][0]-ddP[5][0])-b[2]*ddP[6][0];
    ddN12[0][1]=-b[1]*(ddP[8][3]-ddP[5][3])-b[2]*ddP[6][3]; 
    ddN12[0][2]=-b[1]*(ddP[8][4]-ddP[5][4])-b[2]*ddP[6][4];

    ddN12[1][0]=-b[1]*(ddP[8][3]-ddP[5][3])-b[2]*ddP[6][3]; 
    ddN12[1][1]=-b[1]*(ddP[8][1]-ddP[5][1])-b[2]*ddP[6][1];
    ddN12[1][2]=-b[1]*(ddP[8][5]-ddP[5][5])-b[2]*ddP[6][5];

    ddN12[2][0]=-b[1]*(ddP[8][4]-ddP[5][4])-b[2]*ddP[6][4]; 
    ddN12[2][1]=-b[1]*(ddP[8][5]-ddP[5][5])-b[2]*ddP[6][5];
    ddN12[2][2]=-b[1]*(ddP[8][2]-ddP[5][2])-b[2]*ddP[6][2];

    vector<vector<double> > ddN13(3,dummyVec3);
    ddN13[0][0]=-c[1]*(ddP[8][0]-ddP[5][0])-c[2]*ddP[6][0];
    ddN13[0][1]=-c[1]*(ddP[8][3]-ddP[5][3])-c[2]*ddP[6][3];
    ddN13[0][2]=-c[1]*(ddP[8][4]-ddP[5][4])-c[2]*ddP[6][4];

    ddN13[1][0]=-c[1]*(ddP[8][3]-ddP[5][3])-c[2]*ddP[6][3];
    ddN13[1][1]=-c[1]*(ddP[8][1]-ddP[5][1])-c[2]*ddP[6][1];
    ddN13[1][2]=-c[1]*(ddP[8][5]-ddP[5][5])-c[2]*ddP[6][5];

    ddN13[2][0]=-c[1]*(ddP[8][4]-ddP[5][4])-c[2]*ddP[6][4];
    ddN13[2][1]=-c[1]*(ddP[8][5]-ddP[5][5])-c[2]*ddP[6][5];
    ddN13[2][2]=-c[1]*(ddP[8][2]-ddP[5][2])-c[2]*ddP[6][2];

    vector<vector<double> > ddN21(3,dummyVec3);
    ddN21[0][0]=ddP[1][0]-ddP[4][0]+ddP[3][0]+2*(ddP[7][0]-ddP[6][0]);
    ddN21[0][1]=ddP[1][3]-ddP[4][3]+ddP[3][3]+2*(ddP[7][3]-ddP[6][3]);
    ddN21[0][2]=ddP[1][4]-ddP[4][4]+ddP[3][4]+2*(ddP[7][4]-ddP[6][4]);

    ddN21[1][0]=ddP[1][3]-ddP[4][3]+ddP[3][3]+2*(ddP[7][3]-ddP[6][3]);
    ddN21[1][1]=ddP[1][1]-ddP[4][1]+ddP[3][1]+2*(ddP[7][1]-ddP[6][1]);
    ddN21[1][2]=ddP[1][5]-ddP[4][5]+ddP[3][5]+2*(ddP[7][5]-ddP[6][5]);

    ddN21[2][0]=ddP[1][4]-ddP[4][4]+ddP[3][4]+2*(ddP[7][4]-ddP[6][4]);
    ddN21[2][1]=ddP[1][5]-ddP[4][5]+ddP[3][5]+2*(ddP[7][5]-ddP[6][5]);
    ddN21[2][2]=ddP[1][2]-ddP[4][2]+ddP[3][2]+2*(ddP[7][2]-ddP[6][2]);

    vector<vector<double> > ddN22(3,dummyVec3);
    ddN22[0][0]=-b[2]*(ddP[6][0]-ddP[3][0])-b[0]*ddP[7][0];
    ddN22[0][1]=-b[2]*(ddP[6][3]-ddP[3][3])-b[0]*ddP[7][3];
    ddN22[0][2]=-b[2]*(ddP[6][4]-ddP[3][4])-b[0]*ddP[7][4];

    ddN22[1][0]=-b[2]*(ddP[6][3]-ddP[3][3])-b[0]*ddP[7][3];    
    ddN22[1][1]=-b[2]*(ddP[6][1]-ddP[3][1])-b[0]*ddP[7][1];
    ddN22[1][2]=-b[2]*(ddP[6][5]-ddP[3][5])-b[0]*ddP[7][5];

    ddN22[2][0]=-b[2]*(ddP[6][4]-ddP[3][4])-b[0]*ddP[7][4];
    ddN22[2][1]=-b[2]*(ddP[6][5]-ddP[3][5])-b[0]*ddP[7][5];
    ddN22[2][2]=-b[2]*(ddP[6][2]-ddP[3][2])-b[0]*ddP[7][2];

    vector<vector<double> > ddN23(3,dummyVec3);
    ddN23[0][0]=-c[2]*(ddP[6][0]-ddP[3][0])-c[0]*ddP[7][0];
    ddN23[0][1]=-c[2]*(ddP[6][3]-ddP[3][3])-c[0]*ddP[7][3];
    ddN23[0][2]=-c[2]*(ddP[6][4]-ddP[3][4])-c[0]*ddP[7][4];

    ddN23[1][0]=-c[2]*(ddP[6][3]-ddP[3][3])-c[0]*ddP[7][3];
    ddN23[1][1]=-c[2]*(ddP[6][1]-ddP[3][1])-c[0]*ddP[7][1];
    ddN23[1][2]=-c[2]*(ddP[6][5]-ddP[3][5])-c[0]*ddP[7][5];

    ddN23[2][0]=-c[2]*(ddP[6][4]-ddP[3][4])-c[0]*ddP[7][4];
    ddN23[2][1]=-c[2]*(ddP[6][5]-ddP[3][5])-c[0]*ddP[7][5];
    ddN23[2][2]=-c[2]*(ddP[6][2]-ddP[3][2])-c[0]*ddP[7][2];

    vector<vector<double> > ddN31(3,dummyVec3);
    ddN31[0][0]=ddP[2][0]-ddP[5][0]+ddP[4][0]+2*(ddP[8][0]-ddP[7][0]);
    ddN31[0][1]=ddP[2][3]-ddP[5][3]+ddP[4][3]+2*(ddP[8][3]-ddP[7][3]);
    ddN31[0][2]=ddP[2][4]-ddP[5][4]+ddP[4][4]+2*(ddP[8][4]-ddP[7][4]);

    ddN31[1][0]=ddP[2][3]-ddP[5][3]+ddP[4][3]+2*(ddP[8][3]-ddP[7][3]);
    ddN31[1][1]=ddP[2][1]-ddP[5][1]+ddP[4][1]+2*(ddP[8][1]-ddP[7][1]);
    ddN31[1][2]=ddP[2][5]-ddP[5][5]+ddP[4][5]+2*(ddP[8][5]-ddP[7][5]);

    ddN31[2][0]=ddP[2][4]-ddP[5][4]+ddP[4][4]+2*(ddP[8][4]-ddP[7][4]);
    ddN31[2][1]=ddP[2][5]-ddP[5][5]+ddP[4][5]+2*(ddP[8][5]-ddP[7][5]);
    ddN31[2][2]=ddP[2][2]-ddP[5][2]+ddP[4][2]+2*(ddP[8][2]-ddP[7][2]);

    vector<vector<double> > ddN32(3,dummyVec3);
    ddN32[0][0]=-b[0]*(ddP[7][0]-ddP[4][0])-b[1]*ddP[8][0];
    ddN32[0][1]=-b[0]*(ddP[7][3]-ddP[4][3])-b[1]*ddP[8][3];
    ddN32[0][2]=-b[0]*(ddP[7][4]-ddP[4][4])-b[1]*ddP[8][4];

    ddN32[1][0]=-b[0]*(ddP[7][3]-ddP[4][3])-b[1]*ddP[8][3];
    ddN32[1][1]=-b[0]*(ddP[7][1]-ddP[4][1])-b[1]*ddP[8][1];
    ddN32[1][2]=-b[0]*(ddP[7][5]-ddP[4][5])-b[1]*ddP[8][5];

    ddN32[2][0]=-b[0]*(ddP[7][4]-ddP[4][4])-b[1]*ddP[8][4];
    ddN32[2][1]=-b[0]*(ddP[7][5]-ddP[4][5])-b[1]*ddP[8][5];
    ddN32[2][2]=-b[0]*(ddP[7][2]-ddP[4][2])-b[1]*ddP[8][2];

    vector<vector<double> > ddN33(3,dummyVec3);
    ddN33[0][0]=-c[0]*(ddP[7][0]-ddP[4][0])-c[1]*ddP[8][0];
    ddN33[0][1]=-c[0]*(ddP[7][3]-ddP[4][3])-c[1]*ddP[8][3];
    ddN33[0][2]=-c[0]*(ddP[7][4]-ddP[4][4])-c[1]*ddP[8][4];

    ddN33[1][0]=-c[0]*(ddP[7][3]-ddP[4][3])-c[1]*ddP[8][3];
    ddN33[1][1]=-c[0]*(ddP[7][1]-ddP[4][1])-c[1]*ddP[8][1];
    ddN33[1][2]=-c[0]*(ddP[7][5]-ddP[4][5])-c[1]*ddP[8][5];

    ddN33[2][0]=-c[0]*(ddP[7][4]-ddP[4][4])-c[1]*ddP[8][4];
    ddN33[2][1]=-c[0]*(ddP[7][5]-ddP[4][5])-c[1]*ddP[8][5];
    ddN33[2][2]=-c[0]*(ddP[7][2]-ddP[4][2])-c[1]*ddP[8][2];

    vector<double> dummy2(2,0.);
    vector<vector<double> > ddN11p(2,dummy2);
    vector<vector<double> > ddN12p(2,dummy2);
    vector<vector<double> > ddN13p(2,dummy2);
    vector<vector<double> > ddN21p(2,dummy2);
    vector<vector<double> > ddN22p(2,dummy2);
    vector<vector<double> > ddN23p(2,dummy2);
    vector<vector<double> > ddN31p(2,dummy2);
    vector<vector<double> > ddN32p(2,dummy2);
    vector<vector<double> > ddN33p(2,dummy2);

    //- first order derivatives

    vector<vector<double> > A(2,dummyVec3);
    double rDelta0p5=0.5/Delta;
    for (int d=0;d<DIM;++d) {
      A[0][d]=rDelta0p5*b[d];
      A[1][d]=rDelta0p5*c[d];
    }
    double dN1dx = b[0]*rDelta0p5;
    double dN1dy = c[0]*rDelta0p5;
    double dN2dx = b[1]*rDelta0p5;
    double dN2dy = c[1]*rDelta0p5;
    double dN3dx = b[2]*rDelta0p5;
    double dN3dy = c[2]*rDelta0p5;

    // - second order derivatives

    // dd = [  d^2/dL1^2 d^2/dL1dL2 d^2/dL1dL3 
    //                   d^2/dL2^2  d^2/dL2dL3 
    //         sym                  d^2/dL3^2 ]

    // ddNij is shapefunction for node i, i.e. i in eq.(4.57)
    // shapefunction row j, i.e j=1 for w, j=2 for thetx, j=3 for thety
    // d^2P/dL1^2

    // add matrix multiplications

    amaT_mult(A,ddN11,ddN11p);
    amaT_mult(A,ddN12,ddN12p);
    amaT_mult(A,ddN13,ddN13p);
    amaT_mult(A,ddN21,ddN21p);
    amaT_mult(A,ddN22,ddN22p);
    amaT_mult(A,ddN23,ddN23p);
    amaT_mult(A,ddN31,ddN31p);
    amaT_mult(A,ddN32,ddN32p);
    amaT_mult(A,ddN33,ddN33p);    

    vector<double> dummyVec2(15,0.);
    vector<vector<double> > B(6,dummyVec2);

    B[0][ 0]=dN1dx;
    B[0][ 5]=dN2dx;
    B[0][10]=dN3dx;
    B[1][1]=dN1dy;
    B[1][6]=dN2dy;
    B[1][11]=dN3dy;
    B[2][0]=dN1dy;
    B[2][1]=dN1dx;
    B[2][5]=dN2dy;
    B[2][6]=dN2dx;
    B[2][10]=dN3dy;
    B[2][11]=dN3dx;
    B[3][2 ]=ddN11p[0][0];
    B[3][3 ]=ddN12p[0][0];
    B[3][4 ]=ddN13p[0][0];
    B[3][7 ]=ddN21p[0][0];
    B[3][8 ]=ddN22p[0][0];
    B[3][9 ]=ddN23p[0][0];
    B[3][12]=ddN31p[0][0];
    B[3][13]=ddN32p[0][0];
    B[3][14]=ddN33p[0][0];
    B[4][2 ]=ddN11p[1][1];
    B[4][3 ]=ddN12p[1][1];
    B[4][4 ]=ddN13p[1][1];
    B[4][7 ]=ddN21p[1][1];
    B[4][8 ]=ddN22p[1][1];
    B[4][9 ]=ddN23p[1][1];
    B[4][12]=ddN31p[1][1];
    B[4][13]=ddN32p[1][1];
    B[4][14]=ddN33p[1][1];
    B[5][2 ]=2.*ddN11p[0][1];
    B[5][3 ]=2.*ddN12p[0][1];
    B[5][4 ]=2.*ddN13p[0][1];
    B[5][7 ]=2.*ddN21p[0][1];
    B[5][8 ]=2.*ddN22p[0][1];
    B[5][9 ]=2.*ddN23p[0][1];
    B[5][12]=2.*ddN31p[0][1];
    B[5][13]=2.*ddN32p[0][1];
    B[5][14]=2.*ddN33p[0][1];
#if (DEBUG==1)
    //dumpMat(B,"LTE Bmat");
#endif

    vector<vector<double> > D;
    compDloc(D,Delta,faceID);
    vector<vector<double> > KKlMat(15,dummyVec2);

    aTma_mult(B,D,KKlMat);

    //over-write in-plane DOFs
    vector<double> xtri(3,0.);
    xtri[0]=x0[0];
    xtri[1]=x1[0];
    xtri[2]=x2[0];
    vector<double> ytri(3,0.);
    ytri[0]=x0[1];
    ytri[1]=x1[1];
    ytri[2]=x2[1];
    vector<double> dm(9,0.);
    double nu=poisRatio_[faceID];
    double E =eModulus_[faceID];
    double d =thickness_[faceID];
    dm[0]=   d*E/(1.-nu*nu);
    dm[1]=nu*d*E/(1.-nu*nu);
    dm[3]=nu*d*E/(1.-nu*nu);
    dm[4]=   d*E/(1.-nu*nu);
    dm[8]=(1.-nu)*0.5*d*E/(1.-nu*nu);
    double  alpha=0.;
    double f=1.;
    vector<int> ls(9,0);
    ls[0]=1;
    ls[1]=2;
    ls[2]=4;
    ls[3]=5;
    ls[4]=7;
    ls[5]=8;
    ls[6]=3;
    ls[7]=6;
    ls[8]=9;
    vector<double> smOut(81,0.);
    int m=9;

    double beta=0.5;
    //Call from Fortran
    sm4m2_(&xtri[0],&ytri[0],&dm[0],&alpha,&beta,&ls[0],&smOut[0],&m);
    //    sm3mb_(&xtri[0],&ytri[0],&dm[0],&alpha,&f,&ls[0],&smOut[0],&m);
    //(x, y, dm, alpha, beta, ls, sm, m, status)
    int nby=ndof_loc;
    for (int n2=0;n2<3;++n2) {
      for (int n1=0;n1<3;++n1) {
        KKl[0+n2*nby][0+n1*nby]=smOut[ 0+3*n1+n2*27];
        KKl[0+n2*nby][1+n1*nby]=smOut[ 1+3*n1+n2*27];
        KKl[0+n2*nby][5+n1*nby]=smOut[ 2+3*n1+n2*27];
        //       2         1
        KKl[1+n2*nby][0+n1*nby]=smOut[ 9+3*n1+n2*27];
        KKl[1+n2*nby][1+n1*nby]=smOut[10+3*n1+n2*27];
        KKl[1+n2*nby][5+n1*nby]=smOut[11+3*n1+n2*27];
        //       2         1      
        KKl[5+n2*nby][0+n1*nby]=smOut[18+3*n1+n2*27];
        KKl[5+n2*nby][1+n1*nby]=smOut[19+3*n1+n2*27];
        KKl[5+n2*nby][5+n1*nby]=smOut[20+3*n1+n2*27];
      }  
    }

    for (int r=0;r<3;++r) {
      for (int c=0;c<3;++c) {
        for (int rr=0;rr<5;++rr) {
          for (int cc=0;cc<5;++cc) {
            KKl[rr+r*6][cc+c*6]=KKlMat[rr+r*5][cc+c*5];
          }
        }
      }
    }

    for (int r=0;r<3;++r) {
      for (int c=0;c<3;++c) {
        for (int rr=2;rr<5;++rr) {
          for (int cc=2;cc<5;++cc) {
            KKl[rr+r*6][cc+c*6]=KKlMat[rr+r*5][cc+c*5];
          }
        }
      }
    }
  }

  void linearTriElement::comp_Ke_MITC3(vector<double>& Ke_loc,
                                       vector<double>& Fvec,
                                       const int&      faceID)
  {

    int nnodes_loc   = elementData_[faceID].nnodes;
    int ndof_loc     = elementData_[faceID].ndof;
    int nnodesOffset = nnodesLocalOffset_[faceID];

    for (int n=0;n<nnodes_loc;++n) 
      {
        int index = n + nnodesOffset;
        face_nodes_[n] = faces_[index];
      }
    //face nodes
    int fn0 = face_nodes_[0];
    int fn1 = face_nodes_[1];
    int fn2 = face_nodes_[2];

    //global node coordinates
    vector<double> xx0(DIM,0.);      
    vector<double> xx1(DIM,0.);      
    vector<double> xx2(DIM,0.);

    //nodes_.resize(3,vector<double>(3,0.));      
    for (int idir=0;idir<DIM;++idir) {
      xx0[idir] = nodes_[fn0*xyz+idir];
      xx1[idir] = nodes_[fn1*xyz+idir];
      xx2[idir] = nodes_[fn2*xyz+idir];
    } 
    //shift global coordinates
    vector<double> xx0_tilda(DIM,0.);      
    vector<double> xx1_tilda(DIM,0.);      
    vector<double> xx2_tilda(DIM,0.);

    for (int idir=0;idir<DIM;++idir) {
      xx0_tilda[idir] = 0.;
      xx1_tilda[idir] = xx1[idir] - xx0[idir];
      xx2_tilda[idir] = xx2[idir] - xx0[idir];
    }      
    //compute ex', ey', and, ez' and form T (transformation matrix)
    vector<double> exPrime(DIM,0.);      
    vector<double> eyPrime(DIM,0.);      
    vector<double> eyTilda(DIM,0.);      
    vector<double> ezPrime(DIM,0.);

    for (int idir=0;idir<DIM;++idir) {
      exPrime[idir] = xx1_tilda[idir];
      eyTilda[idir] = xx2_tilda[idir];
    }
    normalizeVector(exPrime);
    normalizeVector(eyTilda);
    ezPrime = cross(exPrime,eyTilda);
    normalizeVector(ezPrime);
    eyPrime = cross(ezPrime,exPrime);
    normalizeVector(eyPrime);

    vector<vector<double> > dispT;
    vector<vector<double> > dispTT;
    dispTT.resize(DIM,vector<double>(DIM,0.));
    dispT.resize(DIM,vector<double>(DIM,0.));
    for (int idir=0;idir<DIM;++idir) { 
      dispTT[0][idir] = exPrime[idir];
      dispTT[1][idir] = eyPrime[idir];
      dispTT[2][idir] = ezPrime[idir];
    }    
    matrixTranspose(dispTT,dispT);

    //transform shifted global coordinates xx0_tilda,xx1_tilda,xx2_tilda to local x0,x1,x2 with T'
    vector<double> x0(DIM,0.);      
    vector<double> x1(DIM,0.);      
    vector<double> x2(DIM,0.);
    mTv_mult(dispTT,xx0_tilda,x0);
    mTv_mult(dispTT,xx1_tilda,x1);
    mTv_mult(dispTT,xx2_tilda,x2);    

    //integration points
    vector<double> rPoint(3,0.);
    vector<double> sPoint(3,0.);
    rPoint[0] = 0.1666666666667;
    rPoint[1] = 0.6666666666667;
    rPoint[2] = 0.1666666666667;
    sPoint[0] = 0.1666666666667;
    sPoint[1] = 0.1666666666667;
    sPoint[2] = 0.6666666666667;
    //vector<double> tPoint(3,0.);
    //tPoint[0] = 0.;
    //tPoint[1] = pow(3./5.,0.5);
    //tPoint[2] = -pow(3./5.,0.5);
    //vector<double> t_weights(3,0.);
    //t_weights[0] = 8./9.*0.5;
    //t_weights[1]= 5./9.*0.5;
    //t_weights[2]= 5./9.*0.5;
    vector<double> tPoint(2,0.);
    tPoint[0] =   pow(3.,-0.5);
    tPoint[1] =  -pow(3.,-0.5);
    vector<double> t_weights(2,0.5);
    vector<double> rs_weights(3,0.3333333333333);    
    int noWeights = rs_weights.size();    
    vector<vector<double> > BTCBmat;
    vector<vector<double> > KK;    
    BTCBmat.resize(18,vector<double>(18,0.));
    KK.resize(18,vector<double>(18,0.));

    vector<double> dummy1(DIM,0.);      
    vector<double> dummy2(DIM,0.);
    vector<double> dummy3(DIM,0.);
    for (int i=0;i<DIM;++i) { 
      dummy1[i] = xx1[i] - xx0[i];
      dummy2[i] = xx2[i] - xx0[i];
    }
    dummy3 = cross(dummy1,dummy2);
    double mag = vectorMag(dummy3);
    double Area = 0.5*mag;
    double a = thickness_[faceID];

    //assemble intermediate local stiffness matrix
    for (int k=0;k<noWeights;++k) {
      for (int kk=0;kk<tPoint.size();++kk) {
        compBTCB(BTCBmat,faceID,rPoint[k],sPoint[k],tPoint[kk],Area,x0,x1,x2,exPrime,eyPrime,ezPrime);
        for (int r=0;r<BTCBmat.size();++r) {
          for (int c=0;c<BTCBmat[0].size();++c) {
            KK[r][c] += t_weights[kk] * rs_weights[k] * BTCBmat[r][c]; 
          }
        }
      }
    }
    //wait();
    ////dumpMat(KK,"MITC3 KK");
    vector<vector<double> >transMat(18,vector<double>(18,0.));
    for (int n=0;n<3;++n) {
      for (int i=0;i<3;++i) { 
        for (int j=0;j<3;++j) { 
          int iloc=i+n*ndof_loc;
          int jloc=j+n*ndof_loc;
          transMat[iloc][jloc]=dispTT[i][j];
        }
      }
     for (int i=0;i<3;++i) { 
       for (int j=0;j<3;++j) { 
         int iloc=i+n*ndof_loc+3;
         int jloc=j+n*ndof_loc+3;
         transMat[iloc][jloc]=dispTT[i][j];
       }
     }
    }
    //aTma_mult(transMat,KK,Ke_loc);
    ////dumpMat(transMat,"transMat MITC3");
    //wait();
    ////dumpMat(BTCBmat,"BTCB Linear");
    //wait();
    ////dumpMat(Ke_loc,"Ke linear");
    //wait();
  }
  
  void linearTriElement::compBTCB(vector<vector<double> >& BTCBmat,
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
                                  vector<double> ez) {
    //create (B^T)CB for computation of local stiffness matrix

    int nnodes_loc = elementData_[faceID].nnodes;
    int ndof_loc   = elementData_[faceID].ndof;

    Cmat_.resize(6,vector<double>(6,0.));
    Bmat_.resize(6,vector<double>(18,0.));

    //fill Cmat_
    double k = SHAPEFACTOR;
    double I = A*thickness_[faceID]*eModulus_[faceID]/(1.-poisRatio_[faceID]*poisRatio_[faceID]);
    //in-plane
    Cmat_[0][0] = I;
    Cmat_[0][1] = I*poisRatio_[faceID];
    Cmat_[1][0] = I*poisRatio_[faceID];
    Cmat_[1][1] = I;
    //out-of-plane
    //I *=pow(thic_[faceID],2)/(12.);
    Cmat_[3][3] = I*(1-poisRatio_[faceID])*0.5;
    Cmat_[4][4] = I*k*(1-poisRatio_[faceID])*0.5; 
    Cmat_[5][5] = I*k*(1-poisRatio_[faceID])*0.5;  
    //Cmat_[0][0] = 1.;
    //Cmat_[0][1] = 0.;
    //Cmat_[1][0] = 0.;
    //Cmat_[1][1] = 1.;
    ////out-of-plane
    //Cmat_[3][3] = 1.;
    //Cmat_[4][4] = 1.;
    //Cmat_[5][5] = 1.;

    double a = thickness_[faceID];

    vector<double> g_1(DIM,0.);
    vector<double> g_2(DIM,0.);
    vector<double> g_3(DIM,0.);
    //covariant
    g_1[0] = x1[0] - x0[0];
    g_1[1] = x1[1] - x0[1];
    g_2[0] = x2[0] - x0[0];   
    g_2[1] = x2[1] - x0[1];

    g_3[2] = a/2.;

#if (DEBUG==1)
    //cout << "COVARIANT: " << endl;
    ////dumpVec(g_1,"g_1");
    ////dumpVec(g_2,"g_2");
    ////dumpVec(g_3,"g_3");
    //wait();
#endif
    //create constants for contravariant vector transforms (i.e. g^1 'dot' ex_hat)
    vector<double> g1(DIM,0.);
    vector<double> g2(DIM,0.);
    vector<double> g3(DIM,0.);
    vector<vector<double> > gVec(DIM,vector<double>(DIM,0.));
    vector<vector<vector<double> > > strain(2,vector<vector<double> >(2,vector<double>(18,0.)));
    //vector<double> ex(DIM,0.);
    //vector<double> ey(DIM,0.);
    //vector<double> ez(DIM,0.);
    //contravariant vectors
    vector<double> dumVec(3,0.);
    dumVec = cross(g_2,g_3);
    double dumDub = dprod(g_1,dumVec);
    for (int i=-0;i<3;++i) g1[i] = dumVec[i]/dumDub;
    dumVec = cross(g_3,g_1);
    dumDub = dprod(g_2,dumVec);
    for (int i=-0;i<3;++i) g2[i] = dumVec[i]/dumDub;
    dumVec = cross(g_1,g_2);
    dumDub = dprod(g_3,dumVec);
    for (int i=-0;i<3;++i) g3[i] = dumVec[i]/dumDub;

    //cout << "contravariantssss " << endl;
    ////dumpVec(g1,"g1");
    ////dumpVec(g2,"g2");
    ////dumpVec(g3,"g3");
    ////dumpVec(ex,"ex");
    ////dumpVec(ey,"ey");
    ////dumpVec(ez,"ez");

    int DDDIM=3;
    double g1_g_1=dprod(g1,g_1,DDDIM);
    //cout << "g1_g1_linear " << g1_g_1 << endl;
    g1[0] /=g1_g_1;
    g1[1] /=g1_g_1;
    g1[2] /=g1_g_1;
    double g2_g_2=dprod(g2,g_2,DDDIM);
    //cout << "g2_g2_linear " << g2_g_2 << endl;
    g2[0] /=g2_g_2;
    g2[1] /=g2_g_2;
    g2[2] /=g2_g_2;

    double g3_g_3=dprod(g3,g_3,DDDIM);
    //cout << "g3_g3_linear " << g3_g_3 << endl;
    g3[0] /= g3_g_3;
    g3[1] /= g3_g_3;
    g3[2] /= g3_g_3;

    ////dumpVec(ex,"ex");
    ////dumpVec(ey,"ey");
    ////dumpVec(ez,"ez");

    //ex[0] = 1.;
    //ey[1] = 1.;
    //ez[2] = 1.;

    double g1_ex = dprod(g1,ex,DIM);
    double g1_ey = dprod(g1,ey,DIM);
    double g2_ex = dprod(g2,ex,DIM);
    double g2_ey = dprod(g2,ey,DIM);
    double g3_ez = dprod(g3,ez,DIM);
    ////dumpVec(g1,"g1");
    ////dumpVec(g2,"g2");
    ////dumpVec(g3,"g3");
    //cout << "test covariance" << endl;
    //cout << "g1 dot g_2:" << dprod(g1,g_2) << endl;
    //cout << "g1 dot g_3:" << dprod(g1,g_3) << endl;
    //cout << "g2 dot g_1:" << dprod(g2,g_1) << endl;
    //cout << "g2 dot g_1:" << dprod(g2,g_1) << endl;
    //cout << "g3 dot g_1:" << dprod(g3,g_1) << endl;
    //cout << "g3 dot g_2:" << dprod(g3,g_2) << endl;

    //use contravariant vectors to fill Bmat_
    //exx
                           //err//                                       //esr + ers = 2ers//                             //ess//
    Bmat_[0][0]   =  g1_ex*g1_ex*(x0[0]-x1[0])               + g2_ex*g1_ex*((x0[0]-x2[0])+(x0[0]-x1[0]))        + g2_ex*g2_ex*(x0[0]-x2[0]);         //u1
    Bmat_[0][1]   =  g1_ex*g1_ex*(x0[1]-x1[1])               + g2_ex*g1_ex*((x0[1]-x1[1])+(x0[1]-x2[1]))        + g2_ex*g2_ex*(x0[1]-x2[1]);         //v1
    Bmat_[0][2]   =  0.;                                                                                                                                //w1
    Bmat_[0][3]   =  g1_ex*g1_ex*(-a/2.)*(x0[1]-x1[1])       + g2_ex*g1_ex*(a/2.)*((x1[1]-x0[1])+(x2[1]-x0[1])) + g2_ex*g2_ex*(-a/2.)*(x0[1]-x2[1]); //alpha1
    Bmat_[0][4]   =  g1_ex*g1_ex*( a/2.)*(x0[0]-x1[0])       + g2_ex*g1_ex*(a/2.)*((x0[0]-x1[0])+(x0[0]-x2[0])) + g2_ex*g2_ex*( a/2.)*(x0[0]-x2[0]);  //beta1
    Bmat_[0][5]   =  0.;                                                                                                                    //thetaZ1
    Bmat_[0][6]   = -g1_ex*g1_ex*(x0[0]-x1[0])               + g2_ex*g1_ex*(x2[0]-x0[0]);           //u2                      
    Bmat_[0][7]   = -g1_ex*g1_ex*(x0[1]-x1[1])               + g2_ex*g1_ex*(x2[1]-x0[1]);           //v2                      
    Bmat_[0][8]   =  0.;                                                                               //w2                      
    Bmat_[0][9]   =  g1_ex*g1_ex*(a/2.)*(x0[1]-x1[1])        - g2_ex*g1_ex*(a/2.)*(x2[1]-x0[1]);    //alpha2                  
    Bmat_[0][10]  = -g1_ex*g1_ex*(a/2.)*(x0[0]-x1[0])        - g2_ex*g1_ex*(a/2.)*(x0[0]-x2[0]);    //beta2                   
    Bmat_[0][11]  =  0.;                                                                 //thetaZ2                                 
    //                                                         /////helloooo//                                                                   
    Bmat_[0][12]  =                                           -g2_ex*g1_ex*(x0[0]-x1[0])                        - g2_ex*g2_ex*(x0[0]-x2[0]);         //u3
    Bmat_[0][13]  =                                           -g2_ex*g1_ex*(x0[1]-x1[1])                        - g2_ex*g2_ex*(x0[1]-x2[1]);         //v3
    Bmat_[0][14]  =  0.;                                                                                                                                    //w3
    Bmat_[0][15]  =                                           -g2_ex*g1_ex*(a/2.)*(x1[1]-x0[1])                 + g2_ex*g2_ex*(a/2.)*(x0[1]-x2[1]);  //alpha3
    Bmat_[0][16]  =                                           -g2_ex*g1_ex*(a/2.)*(x0[0]-x1[0])                 - g2_ex*g2_ex*(a/2.)*(x0[0]-x2[0]);  //beta3
    Bmat_[0][17]  =  0.;                                                                                                                              //thetaZ3
    //eyy
    //                       //err//                                       //esr + ers = 2ers//                               //ess//
    Bmat_[1][0]   =  g1_ey*g1_ey*(x0[0]-x1[0])               + g2_ey*g1_ey*((x0[0]-x2[0])+(x0[0]-x1[0]))          + g2_ey*g2_ey*(x0[0]-x2[0]);          //u1
    Bmat_[1][1]   =  g1_ey*g1_ey*(x0[1]-x1[1])               + g2_ey*g1_ey*((x0[1]-x1[1])+(x0[1]-x2[1]))          + g2_ey*g2_ey*(x0[1]-x2[1]);          //v1
    Bmat_[1][2]   =  0.;                                                                                                                            //w1
    Bmat_[1][3]   =  g1_ey*g1_ey*(-a/2.)*(x0[1]-x1[1])       + g2_ey*g1_ey*(a/2.)*((x1[1]-x0[1])+(x2[1]-x0[1]))   + g2_ey*g2_ey*(-a/2.)*(x0[1]-x2[1]);  //alpha1
    Bmat_[1][4]   =  g1_ey*g1_ey*( a/2.)*(x0[0]-x1[0])       + g2_ey*g1_ey*(a/2.)*((x0[0]-x1[0])+(x0[0]-x2[0]))   + g2_ey*g2_ey*( a/2.)*(x0[0]-x2[0]);   //beta1
    Bmat_[1][5]   =  0.;                                                                                                                            //thetaZ1

    Bmat_[1][6]   = -g1_ey*g1_ey*(x0[0]-x1[0])               - g2_ey*g1_ey*(x0[0]-x2[0]);          //u2
    Bmat_[1][7]   = -g1_ey*g1_ey*(x0[1]-x1[1])               - g2_ey*g1_ey*(x0[1]-x2[1]);          //v2
    Bmat_[1][8]   =  0.;                                                                       //w2
    Bmat_[1][9]   =  g1_ey*g1_ey*(a/2.)*(x0[1]-x1[1])        - g2_ey*g1_ey*(a/2.)*(x2[1]-x0[1]);   //alpha2
    Bmat_[1][10]  = -g1_ey*g1_ey*(a/2.)*(x0[0]-x1[0])        - g2_ey*g1_ey*(a/2.)*(x0[0]-x2[0]);   //beta2
    Bmat_[1][11]  =  0.;                                                                       //thetaZ2

    Bmat_[1][12]  =                                           -g2_ey*g1_ey*(x0[0]-x1[0])                          - g2_ey*g2_ey*(x0[0]-x2[0]);         //u3
    Bmat_[1][13]  =                                           -g2_ey*g1_ey*(x0[1]-x1[1])                          - g2_ey*g2_ey*(x0[1]-x2[1]);         //v3
    Bmat_[1][14]  =  0.;                                                                                                                           //w3
    Bmat_[1][15]  =                                           -g2_ey*g1_ey*(a/2.)*(x1[1]-x0[1])                   + g2_ey*g2_ey*(a/2.)*(x0[1]-x2[1]);  //alpha3
    Bmat_[1][16]  =                                           -g2_ey*g1_ey*(a/2.)*(x0[0]-x1[0])                   - g2_ey*g2_ey*(a/2.)*(x0[0]-x2[0]);  //beta3
    Bmat_[1][17]  =  0.;                                                                                                                               //thetaZ3
    //exy (use gamma_xy = 2*epsilon_xy)
    //                        //err//                                            //ers//                                                   //esr//                                        //ess//
    Bmat_[3][0]   =  2.*(g1_ex*g1_ey*(x0[0]-x1[0])           + 0.5*(g1_ex*g2_ey*((x0[0]-x2[0])+(x0[0]-x1[0]))         + g2_ex*g1_ey*((x0[0]-x2[0])+(x0[0]-x1[0])))        + g2_ex*g2_ey*(x0[0]-x2[0]));          //u1
    Bmat_[3][1]   =  2.*(g1_ex*g1_ey*(x0[1]-x1[1])           + 0.5*(g1_ex*g2_ey*((x0[1]-x1[1])+(x0[1]-x2[1]))         + g2_ex*g1_ey*((x0[1]-x1[1])+(x0[1]-x2[1])))        + g2_ex*g2_ey*(x0[1]-x2[1]));          //v1
    Bmat_[3][2]   =  0.;                                                                                                                                                                                         //w1
    Bmat_[3][3]   =  2.*(g1_ex*g1_ey*(-a/2.)*(x0[1]-x1[1])   + 0.5*(g1_ex*g2_ey*(a/2.)*((x1[1]-x0[1])+(x2[1]-x0[1]))  + g2_ex*g1_ey*(a/2.)*((x1[1]-x0[1])+(x2[1]-x0[1]))) + g2_ex*g2_ey*(-a/2.)*(x0[1]-x2[1]));  //alpha1
    Bmat_[3][4]   =  2.*(g1_ex*g1_ey*( a/2.)*(x0[0]-x1[0])   + 0.5*(g1_ex*g2_ey*(a/2.)*((x0[0]-x1[0])+(x0[0]-x2[0]))  + g2_ex*g1_ey*(a/2.)*((x0[0]-x1[0])+(x0[0]-x2[0]))) + g2_ex*g2_ey*( a/2.)*(x0[0]-x2[0]));  //beta1
    Bmat_[3][5]   =  0.;                                                                                                                                                                                        //thetaZ1
    //                                                                         //ers + esr//                                                                                    
    Bmat_[3][6]   =  2.*(-g1_ex*g1_ey*(x0[0]-x1[0])          + 0.5*(g2_ey*g1_ex+g1_ey*g2_ex)*(x2[0]-x0[0]));          //u2                                                     
    Bmat_[3][7]   =  2.*(-g1_ex*g1_ey*(x0[1]-x1[1])          + 0.5*(g2_ey*g1_ex+g1_ey*g2_ex)*(x2[1]-x0[1]));          //v2                                                     
    Bmat_[3][8]   =  0.;                                                                                         //w3                                                          
    Bmat_[3][9]   =  2.*( g1_ex*g1_ey*(a/2.)*(x0[1]-x1[1])   - 0.5*(g2_ey*g1_ex+g2_ex*g1_ey)*(a/2.)*(x2[1]-x0[1]));   //alpha2                                                 
    Bmat_[3][10]  =  2.*(-g1_ex*g1_ey*(a/2.)*(x0[0]-x1[0])   - 0.5*(g2_ey*g1_ex+g2_ex*g1_ey)*(a/2.)*(x0[0]-x2[0]));   //beta2                                                  
    Bmat_[3][11]  =  0.;                                                                                         //thetaZ2                                                     
    Bmat_[3][12]  =  2.*(                                     -0.5*(g2_ey*g1_ex+g2_ex*g1_ey)*(x0[0]-x1[0])                                                                - g2_ex*g2_ey*(x0[0]-x2[0]));         //u3
    Bmat_[3][13]  =  2.*(                                     -0.5*(g2_ey*g1_ex+g2_ex*g1_ey)*(x0[1]-x1[1])                                                                - g2_ex*g2_ey*(x0[1]-x2[1]));         //v3
    Bmat_[3][14]  =  0.;                                                                                                                                                                                       //w3
    Bmat_[3][15]  =  2.*(                                     -0.5*(g2_ey*g1_ex+g2_ex*g1_ey)*(a/2.)*(x1[1]-x0[1])                                                         + g2_ex*g2_ey*(a/2.)*(x0[1]-x2[1]));  //alpha3
    Bmat_[3][16]  =  2.*(                                     -0.5*(g2_ey*g1_ex+g2_ex*g1_ey)*(a/2.)*(x0[0]-x1[0])                                                         - g2_ex*g2_ey*(a/2.)*(x0[0]-x2[0]));  //beta3
    Bmat_[3][17]  =  0.;                                                                                                                                    //thetaZ3
    //eyz (use gamma_yz = 2*epsilon_yz)                                  
    //                                   //ert//                                                        //est//
    Bmat_[4][0]   =  0.;                                                                                                                               //u1
    Bmat_[4][1]   =  0.;                                                                                                                               //v1
    Bmat_[4][2]   =  a*g1_ey*g3_ez*(-0.25)                                             + a*g2_ey*g3_ez*(-0.25);                                            //w1
    Bmat_[4][3]   =  a*g1_ey*g3_ez*(sPoint*0.125*(x1[1]-x2[1])-0.125*(x1[1]-x0[1]))    - a*g2_ey*g3_ez*(rPoint*0.125*(x1[1]-x2[1])+0.125*(x2[1]-x0[1]));   //alpha1
    Bmat_[4][4]   =  a*g1_ey*g3_ez*(sPoint*0.125*(x2[0]-x1[0])+0.125*(x1[0]-x0[0]))    - a*g2_ey*g3_ez*(rPoint*0.125*(x2[0]-x1[0])-0.125*(x2[0]-x0[0]));   //beta1
    Bmat_[4][5]   =  0.;                                                                                                                              //thetaZ1
    Bmat_[4][6]   =  0.;                                                                                                                               //u2
    Bmat_[4][7]   =  0.;                                                                                                                               //v2
    Bmat_[4][8]   =  a*g1_ey*g3_ez*(0.25);                                                                                                                //w2
    Bmat_[4][9]   =  a*g1_ey*g3_ez*(sPoint*0.125*(x2[1]-x0[1])-0.125*(x1[1]-x0[1]))    - a*g2_ey*g3_ez*(rPoint*0.125*(x2[1]-x0[1]));    //alpha2
    Bmat_[4][10]  =  a*g1_ey*g3_ez*(sPoint*0.125*(x0[0]-x2[0])+0.125*(x1[0]-x0[0]))    - a*g2_ey*g3_ez*(rPoint*0.125*(x0[0]-x2[0]));    //beta2
    Bmat_[4][11]  =  0.;                                                                                                                               //thetaZ2
    Bmat_[4][12]  =  0.;                                                                                                             //u3
    Bmat_[4][13]  =  0.;                                                                                                             //v3
    Bmat_[4][14]  =                                                                      a*g2_ey*g3_ez*(0.25);                     //     new         //       //w3
    Bmat_[4][15]  =  a*g1_ey*g3_ez*(sPoint*0.125*(x0[1]-x1[1]))                        + a*g2_ey*g3_ez*(rPoint*0.125*(x1[1]-x0[1])-0.125*(x2[1]-x0[1]));      //alpha3
    Bmat_[4][16]  =  a*g1_ey*g3_ez*(sPoint*0.125*(x1[0]-x0[0]))                        - a*g2_ey*g3_ez*(rPoint*0.125*(x1[0]-x0[0])-0.125*(x2[0]-x0[0]));      //beta3  
    Bmat_[4][17]  =  0.;                                                                                                              //thetaZ3
    //exz (use gamma_xz = 2*epsilon_xz)
    //                                    //ert//                                                        //est//
    Bmat_[5][0]   =  0.;                                                                                                                                //u1 
    Bmat_[5][1]   =  0.;                                                                                                                                //v1            
    Bmat_[5][2]   =  a*g1_ex*g3_ez*(-0.25)                                             + a*g2_ex*g3_ez*(-0.25);                                             //w1
    Bmat_[5][3]   =  a*g1_ex*g3_ez*(sPoint*0.125*(x1[1]-x2[1])-0.125*(x1[1]-x0[1]))    - a*g2_ex*g3_ez*(rPoint*0.125*(x1[1]-x2[1])+0.125*(x2[1]-x0[1]));    //alpha1
    Bmat_[5][4]   =  a*g1_ex*g3_ez*(sPoint*0.125*(x2[0]-x1[0])+0.125*(x1[0]-x0[0]))    - a*g2_ex*g3_ez*(rPoint*0.125*(x2[0]-x1[0])-0.125*(x2[0]-x0[0]));    //beta1
    Bmat_[5][5]   =  0.;                                                                                                                                 //thetaZ1
    Bmat_[5][6]   =  0.;                                                                                                                                   //u2
    Bmat_[5][7]   =  0.;                                                                                                                                 //v2
    Bmat_[5][8]   =  a*g1_ex*g3_ez*(0.25);                                                                                                                //w2
    Bmat_[5][9]   =  a*g1_ex*g3_ez*(sPoint*0.125*(x2[1]-x0[1])-0.125*(x1[1]-x0[1]))    - a*g2_ex*g3_ez*(rPoint*0.125*(x2[1]-x0[1]));    //alpha2
    Bmat_[5][10]  =  a*g1_ex*g3_ez*(sPoint*0.125*(x0[0]-x2[0])+0.125*(x1[0]-x0[0]))    - a*g2_ex*g3_ez*(rPoint*0.125*(x0[0]-x2[0]));    //beta2
    Bmat_[5][11]  =  0.;                                                                                                                                //thetaZ2
    Bmat_[5][12]  =  0.;                                                                                                           //u3
    Bmat_[5][13]  =  0.;                                                                                                           //v3
    Bmat_[5][14]  =                                                                      a*g2_ex*g3_ez*(0.25);                         //w3
    Bmat_[5][15]  =  a*g1_ex*g3_ez*(sPoint*0.125*(x0[1]-x1[1]))                        + a*g2_ex*g3_ez*(rPoint*0.125*(x1[1]-x0[1])-0.125*(x2[1]-x0[1]));   //alpha3
    Bmat_[5][16]  =  a*g1_ex*g3_ez*(sPoint*0.125*(x1[0]-x0[0]))                        - a*g2_ex*g3_ez*(rPoint*0.125*(x1[0]-x0[0])-0.125*(x2[0]-x0[0]));   //beta3
    Bmat_[5][17]  =  0.;                                                                                                         //hellloooo/           //thetaZ3

    for (int n=0;n<3;++n) {
      for (int nt=3;nt<5;++nt) {
        Bmat_[0][nt+6*n] *= tPoint;
        Bmat_[1][nt+6*n] *= tPoint;
        Bmat_[3][nt+6*n] *= tPoint;
      }
    }

    for (int n=0;n<18;++n) {
      Bmat_[4][n] *= 2.;
      Bmat_[5][n] *= 2.;
    } 

    //for (int n=0;n<6;++n) {
    //  Bmat_[n][4] =-Bmat_[n][4];
    //  Bmat_[n][10]=-Bmat_[n][10];
    //  Bmat_[n][16]=-Bmat_[n][16];
    //}

    ////dumpMat(Bmat_,"Bmat linear");
    //wait();
    ////dumpMat(Bmat_,"Bmat");
    vector<vector<double> > BTmat  (18,vector<double>(6,0.));
    vector<vector<double> > BTCmat (18,vector<double>(6,0.));
    matrixTranspose(Bmat_,BTmat);
    mm_mult(BTmat,Cmat_,BTCmat);
    mm_mult(BTCmat,Bmat_,BTCBmat);
    //set drilling dof to 10E-10
    BTCBmat[5][5]   = 0.00000000001;
    BTCBmat[11][11] = 0.00000000001;
    BTCBmat[17][17] = 0.00000000001;

    //dumpMat(BTCBmat, " BTCB");
    wait();
  }  
  
  void linearTriElement::comp_NL_Ke_MITC3(vector<double>& Ke_loc,
                                          vector<double>& Fvec,
                                          const int&      faceID)
  {

    int nnodes_loc   = elementData_[faceID].nnodes;
    int ndof_loc     = elementData_[faceID].ndof;
    int nnodesOffset = nnodesLocalOffset_[faceID];

    for (int n=0;n<nnodes_loc;++n) 
      {
        int index = n + nnodesOffset;
        face_nodes_[n] = faces_[index];
      }
    //face nodes
    int fn0 = face_nodes_[0];
    int fn1 = face_nodes_[1];
    int fn2 = face_nodes_[2];

    for (int idir=0;idir<DIM;++idir) 
      {
        //current node locations
        x0_t_[idir] = nodes_[fn0*xyz+idir];
        x1_t_[idir] = nodes_[fn1*xyz+idir];
        x2_t_[idir] = nodes_[fn2*xyz+idir];
        //intial node locations
        x0_0_[idir] = nodes0_[fn0*xyz+idir];
        x1_0_[idir] = nodes0_[fn1*xyz+idir];
        x2_0_[idir] = nodes0_[fn2*xyz+idir];
      }
    //integration points
    //r-s plane
    vector<double> rPoint(3,0.);
    vector<double> sPoint(3,0.);
    rPoint[0] =  0.1666666666667;
    rPoint[1] =  0.6666666666667;
    rPoint[2] =  0.1666666666667;
    sPoint[0] =  0.1666666666667;
    sPoint[1] =  0.1666666666667;
    sPoint[2] =  0.6666666666667;
    vector<double> rs_weights(3,0.3333333333333);        
    //t-plane
    vector<double> tPoint(2,0.);
    tPoint[0] = -0.5773502691896257;
    tPoint[1] =  0.5773502691896257;
    vector<double> t_weights(2,0.5);
    //calculate area of element over reference configuration
    if (firstTimeComputingK_ == true)
      {
        for (int i=0;i<DIM;++i) 
          { 
            dummy1_[i] = x1_0_[i] - x0_0_[i];
            dummy2_[i] = x2_0_[i] - x0_0_[i];
          }
        dummy3_ = cross(dummy1_,dummy2_);
        double mag    = vectorMag(dummy3_);
        Area_[faceID] = 0.5*mag;
      }

    for (int dir=0;dir<3;++dir)
      {
        Vn0_0_[dir] = Vn0_[myLocalNodesNormalVectorStatus_[fn0]*xyz+dir];
        Vn0_1_[dir] = Vn0_[myLocalNodesNormalVectorStatus_[fn1]*xyz+dir];
        Vn0_2_[dir] = Vn0_[myLocalNodesNormalVectorStatus_[fn2]*xyz+dir];

        Vnt_0_[dir] = Vnt_[myLocalNodesNormalVectorStatus_[fn0]*xyz+dir];
        Vnt_1_[dir] = Vnt_[myLocalNodesNormalVectorStatus_[fn1]*xyz+dir];
        Vnt_2_[dir] = Vnt_[myLocalNodesNormalVectorStatus_[fn2]*xyz+dir];

        V2t_0_[dir] = V2t_[myLocalNodesNormalVectorStatus_[fn0]*xyz+dir];
        V2t_1_[dir] = V2t_[myLocalNodesNormalVectorStatus_[fn1]*xyz+dir];
        V2t_2_[dir] = V2t_[myLocalNodesNormalVectorStatus_[fn2]*xyz+dir];

        V1t_0_[dir] = V1t_[myLocalNodesNormalVectorStatus_[fn0]*xyz+dir];
        V1t_1_[dir] = V1t_[myLocalNodesNormalVectorStatus_[fn1]*xyz+dir];
        V1t_2_[dir] = V1t_[myLocalNodesNormalVectorStatus_[fn2]*xyz+dir];

        deltaVnt_0_[dir] = deltaVnt_[myLocalNodesNormalVectorStatus_[fn0]*xyz+dir];
        deltaVnt_1_[dir] = deltaVnt_[myLocalNodesNormalVectorStatus_[fn1]*xyz+dir];
        deltaVnt_2_[dir] = deltaVnt_[myLocalNodesNormalVectorStatus_[fn2]*xyz+dir];

      }

    //dumpVec(Vn0_0_,"Vn0_0");
    //dumpVec(Vn0_1_,"Vn0_1");
    //dumpVec(Vn0_2_,"Vn0_2");
    //cout << endl;
    //dumpVec(Vnt_0_,"Vnt_0");
    //dumpVec(Vnt_1_,"Vnt_1");
    //dumpVec(Vnt_2_,"Vnt_2");
    //cout << endl;
    //dumpVec(V1t_0_,"V1t_0");
    //dumpVec(V1t_1_,"V1t_1");
    //dumpVec(V1t_2_,"V1t_2");
    //cout << endl;
    //dumpVec(V2t_0_,"V2t_0");
    //dumpVec(V2t_1_,"V2t_1");
    //dumpVec(V2t_2_,"V2t_2");
    //cout << endl;
    //dumpVec(deltaVnt_0_,"deltaVnt_0");
    //dumpVec(deltaVnt_1_,"deltaVnt_1");
    //dumpVec(deltaVnt_2_,"deltaVnt_2");
    //cout << endl;
    //wait();

    //assemble intermediate local stiffness matrix
    for (int k=0;k<3;++k)//3 points in the plane
      {
	for (int kk=0;kk<2;++kk)//2 points through the thickness
	  {

	    comp_NL_BTCB(faceID,rPoint[k],sPoint[k],tPoint[kk]);

	    for (int r=0;r<nnodes_loc*ndof_loc;++r)
	      {
		for (int c=0;c<nnodes_loc*ndof_loc;++c) 
		  {
		    Ke_loc[c + r*nnodes_loc*ndof_loc] += t_weights[kk] * rs_weights[k] * BTCBmat_[r][c];
		  }
		Fvec[r] += t_weights[kk] * rs_weights[k] * Fvec_loc_[r];
	      }
	    for (int i=0;i<3;++i) for (int j=0;j<3;++j) cauchyStressTensor_[faceID][i][j] += t_weights[kk] * rs_weights[k] * integratingCauchyST_2_[i][j];
	  }
      }

    calculateVonMisesStress(faceID);

  }

  void linearTriElement::comp_NL_BTCB(const int& faceID,
                                      double& rPoint,
                                      double& sPoint,
                                      double& tPoint)
    
  {

    int nnodes_loc = elementData_[faceID].nnodes;
    int ndof_loc   = elementData_[faceID].ndof;

    //create Qmat for transforming C from r,s,t to global reference configuration
    if (firstTimeComputingK_ == true)
      {
        vector<vector<double> > Qsh(6,vector<double>(6,0.));
        vector<vector<double> > QshT(6,vector<double>(6,0.));
        vector<vector<double> > QshTC(6,vector<double>(6,0.));
        vector<double> rvec0(3,0.);
        vector<double> svec0(3,0.);
        vector<double> rBar(3,0.);
        vector<double> sBar(3,0.);
        vector<double> tBar(3,0.);

        for (int i=0;i<3;++i)
          {
            rvec0[i] = x1_0_[i] - x0_0_[i];
            svec0[i] = x2_0_[i] - x0_0_[i];
          }

        tBar = cross(rvec0,svec0);
        rBar = cross(svec0,tBar);
        sBar = cross(tBar,rBar);
        normalizeVector(tBar);
        normalizeVector(rBar);
        normalizeVector(sBar);

        double l1 = dprod(ex_,rBar);
        double l2 = dprod(ex_,sBar);
        double l3 = dprod(ex_,tBar);

        double m1 = dprod(ey_,rBar);
        double m2 = dprod(ey_,sBar);
        double m3 = dprod(ey_,tBar);

        double n1 = dprod(ez_,rBar);
        double n2 = dprod(ez_,sBar);
        double n3 = dprod(ez_,tBar);

        Qsh[0][0] = l1*l1;
        Qsh[0][1] = m1*m1;
        Qsh[0][2] = n1*n1;
        Qsh[0][3] = l1*m1;
        Qsh[0][4] = m1*n1;
        Qsh[0][5] = n1*l1;

        Qsh[1][0] = l2*l2;
        Qsh[1][1] = m2*m2;
        Qsh[1][2] = n2*n2;
        Qsh[1][3] = l2*m2;
        Qsh[1][4] = m2*n2;
        Qsh[1][5] = n2*l2;

        Qsh[2][0] = l3*l3;
        Qsh[2][1] = m3*m3;
        Qsh[2][2] = n3*n3;
        Qsh[2][3] = l3*m3;
        Qsh[2][4] = m3*n3;
        Qsh[2][5] = n3*l3;

        Qsh[3][0] = 2.*l1*l2;
        Qsh[3][1] = 2.*m1*m2;
        Qsh[3][2] = 2.*n1*n2;
        Qsh[3][3] = l1*m2 + l2*m1;
        Qsh[3][4] = m1*n2 + m2*n1;
        Qsh[3][5] = n1*l2 + n2*l1;

        Qsh[4][0] = 2.*l3*l2;
        Qsh[4][1] = 2.*m3*m2;
        Qsh[4][2] = 2.*n3*n2;
        Qsh[4][3] = l3*m2 + l2*m3;
        Qsh[4][4] = m3*n2 + m2*n3;
        Qsh[4][5] = n3*l2 + n2*l3;

        Qsh[5][0] = 2.*l3*l1;
        Qsh[5][1] = 2.*m3*m1;
        Qsh[5][2] = 2.*n3*n1;
        Qsh[5][3] = l3*m1 + l1*m3;
        Qsh[5][4] = m3*n1 + m1*n3;
        Qsh[5][5] = n3*l1 + n1*l3;

        matrixTranspose(Qsh,QshT);
        //fill Cmat_
        double k          = SHAPEFACTOR;
        double I_inplane  = Area_[faceID]*thickness_[faceID]*inPlaneEmod_[faceID] /(1.-poisRatio_[faceID]*poisRatio_[faceID]);
        double I_outplane = Area_[faceID]*thickness_[faceID]*outPlaneEmod_[faceID]/(1.-poisRatio_[faceID]*poisRatio_[faceID]);
        //in-plane
        Cmat_[0][0] = I_inplane;
        Cmat_[0][1] = I_inplane*poisRatio_[faceID];
        Cmat_[1][0] = I_inplane*poisRatio_[faceID];
        Cmat_[1][1] = I_inplane;
        Cmat_[3][3] = I_inplane*(1.-poisRatio_[faceID])*0.5;
        //out-of-plane
        Cmat_[4][4] = I_outplane*k*(1.-poisRatio_[faceID])*0.5;
        Cmat_[5][5] = I_outplane*k*(1.-poisRatio_[faceID])*0.5;

        //cout << " Area_[faceID] = " <<  Area_[faceID] << endl;
        //cout << "thickness_[faceID] = " << thickness_[faceID] << endl;
        //cout << "I_inplane =  "<<I_inplane<<" "<<inPlaneEmod_[faceID]<<endl;
        //cout << "I_outplane = "<<I_outplane<<" "<<outPlaneEmod_[faceID]<<endl;
        //transform constitutive matrix from r,s,t to global
        mm_mult(QshT,Cmat_,QshTC);
        mm_mult(QshTC,Qsh,CmatNL_[faceID]);
        //dumpMat(Qsh,"Qsh");
        //dumpMat(Cmat_,"Cmat");
        //dumpMat(CmatNL_[faceID],"cmatl");
        //wait();
      }

    int nnodesOffset = nnodesLocalOffset_[faceID];
    for (int n=0;n<nnodes_loc;++n) 
      {
        int index = n + nnodesOffset;
        face_nodes_[n] = faces_[index];
      }

    int fn0 = face_nodes_[0];
    int fn1 = face_nodes_[1];
    int fn2 = face_nodes_[2];

    //thickness 
    double a0 = thickness_[faceID];
    double a1 = thickness_[faceID];
    double a2 = thickness_[faceID];

    //2D shape functions
    double h1 = 1. - rPoint - sPoint;
    double h2 = rPoint;
    double h3 = sPoint;
    //2D shape function derivatives
    //h1 = 1-r-s
    double dh1dr = -1.;
    double dh1ds = -1.;
    double dh1dt =  0.;
    //h2 = r
    double dh2dr =  1.;
    double dh2ds =  0.;
    double dh2dt =  0.;
    //h3 = s
    double dh3dr =  0.;
    double dh3ds =  1.;
    double dh3dt =  0.;

    dhmat_[0][0] = dh1dr;
    dhmat_[0][1] = dh2dr;
    dhmat_[0][2] = dh3dr;
    dhmat_[1][0] = dh1ds;
    dhmat_[1][1] = dh2ds;
    dhmat_[1][2] = dh3ds;
    dhmat_[2][0] = dh1dt;
    dhmat_[2][1] = dh2dt;
    dhmat_[2][2] = dh3dt;
    //covariant ZERO set
    for (int i=0;i<3;++i) 
      { 
        //g0_1
        g0_1_[i] = dh1dr*x0_0_[i] + dh2dr*x1_0_[i] + dh3dr*x2_0_[i] + 0.5*tPoint * (a0*dh1dr*Vn0_0_[i] + a1*dh2dr*Vn0_1_[i] + a2*dh3dr*Vn0_2_[i]);
        //g0_2
        g0_2_[i] = dh1ds*x0_0_[i] + dh2ds*x1_0_[i] + dh3ds*x2_0_[i] + 0.5*tPoint * (a0*dh1ds*Vn0_0_[i] + a1*dh2ds*Vn0_1_[i] + a2*dh3ds*Vn0_2_[i]);
        //g0_3
        g0_3_[i] = 0.5*(a0*h1*Vn0_0_[i] + a1*h2*Vn0_1_[i] + a2*h3*Vn0_2_[i]); 
      }

    calculateStrainsMITC3(rPoint,sPoint,tPoint,a0,a1,a2,faceID);
    compPKST(faceID);
    calculateDeformationGradient(rPoint,sPoint,tPoint,a0,a1,a2,faceID);
    PKST2Cauchy(faceID);

    //evaluate strains for Bmat
    //fill in-plane strain translational DOF
    for (int i=0;i<nnodes_loc;++i)
      {
        for (int j=0;j<3;++j)//number of translational DOF per node
          { 
            int jj = j+ndof_loc*i;

            B_rst_[0][jj] = dhmat_[0][i]*gt_matT_[0][j];
            B_rst_[1][jj] = dhmat_[1][i]*gt_matT_[1][j];
            B_rst_[2][jj] = dhmat_[2][i]*gt_matT_[2][j];
            B_rst_[3][jj] = 0.5*(dhmat_[0][i]*gt_matT_[1][j] + dhmat_[1][i]*gt_matT_[0][j]);

          }
      }
    //create 'k' and 'f' vectors for alpha and beta DOF
    //let a0 be the zero-based a1, a1 be a2, a2 be a3
    for (int i=0;i<3;++i) 
      { 
        //k1 = -0.5*t*a1*h1*V2^0
        dk1dr_[i] = -0.5*tPoint*a0*dh1dr*V2t_0_[i];
        dk1ds_[i] = -0.5*tPoint*a0*dh1ds*V2t_0_[i];
        dk1dt_[i] = -0.5*a0*h1*V2t_0_[i];
        //k2 = -0.5*t*a2*h2*V2^1
        dk2dr_[i] = -0.5*tPoint*a1*dh2dr*V2t_1_[i];
        dk2ds_[i] = -0.5*tPoint*a1*dh2ds*V2t_1_[i];
        dk2dt_[i] = -0.5*a1*h2*V2t_1_[i];
        //k3 = -0.5*t*a3*h3*V2^2
        dk3dr_[i] = -0.5*tPoint*a2*dh3dr*V2t_2_[i];
        dk3ds_[i] = -0.5*tPoint*a2*dh3ds*V2t_2_[i];
        dk3dt_[i] = -0.5*a2*h3*V2t_2_[i];        
        //f1 = 0.5*t*a1*h1*V1^0
        df1dr_[i] = 0.5*tPoint*a0*dh1dr*V1t_0_[i];
        df1ds_[i] = 0.5*tPoint*a0*dh1ds*V1t_0_[i];
        df1dt_[i] = 0.5*a0*h1*V1t_0_[i];
        //f2 = 0.5*t*a2*h2*V1^1
        df2dr_[i] = 0.5*tPoint*a1*dh2dr*V1t_1_[i];
        df2ds_[i] = 0.5*tPoint*a1*dh2ds*V1t_1_[i];
        df2dt_[i] = 0.5*a1*h2*V1t_1_[i];
        //f3 = 0.5*t*a3*h3*V1^2
        df3dr_[i] = 0.5*tPoint*a2*dh3dr*V1t_2_[i];
        df3ds_[i] = 0.5*tPoint*a2*dh3ds*V1t_2_[i];
        df3dt_[i] = 0.5*a2*h3*V1t_2_[i];
      }

    //fill in rotational DOF
    //alpha 1
    B_rst_[0][3]   = dprod(dk1dr_,gt_matT_[0]);
    B_rst_[1][3]   = dprod(dk1ds_,gt_matT_[1]);
    B_rst_[2][3]   = dprod(dk1dt_,gt_matT_[2]);
    B_rst_[3][3]   = 0.5*(dprod(dk1dr_,gt_matT_[1]) + dprod(dk1ds_,gt_matT_[0]));
    //alpha2
    B_rst_[0][8]   = dprod(dk2dr_,gt_matT_[0]);
    B_rst_[1][8]   = dprod(dk2ds_,gt_matT_[1]);
    B_rst_[2][8]   = dprod(dk2dt_,gt_matT_[2]);
    B_rst_[3][8]   = 0.5*(dprod(dk2dr_,gt_matT_[1]) + dprod(dk2ds_,gt_matT_[0]));
    //aplha 3
    B_rst_[0][13]  = dprod(dk3dr_,gt_matT_[0]);
    B_rst_[1][13]  = dprod(dk3ds_,gt_matT_[1]);
    B_rst_[2][13]  = dprod(dk3dt_,gt_matT_[2]);
    B_rst_[3][13]  = 0.5*(dprod(dk3dr_,gt_matT_[1]) + dprod(dk3ds_,gt_matT_[0]));
    //beta 1
    B_rst_[0][4]   = dprod(df1dr_,gt_matT_[0]);
    B_rst_[1][4]   = dprod(df1ds_,gt_matT_[1]);
    B_rst_[2][4]   = dprod(df1dt_,gt_matT_[2]);
    B_rst_[3][4]   = 0.5*(dprod(df1dr_,gt_matT_[1]) + dprod(df1ds_,gt_matT_[0]));
    //beta 2
    B_rst_[0][9]  = dprod(df2dr_,gt_matT_[0]);
    B_rst_[1][9]  = dprod(df2ds_,gt_matT_[1]);
    B_rst_[2][9]  = dprod(df2dt_,gt_matT_[2]);
    B_rst_[3][9]  = 0.5*(dprod(df2dr_,gt_matT_[1]) + dprod(df2ds_,gt_matT_[0]));
    //beta 3
    B_rst_[0][14]  = dprod(df3dr_,gt_matT_[0]);
    B_rst_[1][14]  = dprod(df3ds_,gt_matT_[1]);
    B_rst_[2][14]  = dprod(df3dt_,gt_matT_[2]);
    B_rst_[3][14]  = 0.5*(dprod(df3dr_,gt_matT_[1]) + dprod(df3ds_,gt_matT_[0]));

    //MITC tying method
    //add transverse shear strains ert, est
    //in-plane components (u,v) are zero
    //r,s,t tying points
    double r1 = 0.5;
    double s1 = 0.;
    double r2 = 0.;
    double s2 = 0.5;
    double r3 = 0.5;
    double s3 = 0.5;
    double h1_TP1 = 1.-r1-s1;
    double h2_TP1 = r1;
    double h3_TP1 = s1;
    double h1_TP2 = 1.-r2-s2;
    double h2_TP2 = r2;
    double h3_TP2 = s2;
    double h1_TP3 = 1.-r3-s3;
    double h2_TP3 = r3;
    double h3_TP3 = s3;
    //update dk/dt and df/dt derivatives for each tying point
    //dk/dr, dkds, df/dr, df/ds are independent of r,s, and t tying points
    for (int i=0;i<3;++i) 
      { 
        dk1dt_TP1_[i] = -0.5*a0*h1_TP1*V2t_0_[i];
        dk2dt_TP1_[i] = -0.5*a1*h2_TP1*V2t_1_[i];
        dk3dt_TP1_[i] = -0.5*a2*h3_TP1*V2t_2_[i];
        df1dt_TP1_[i] =  0.5*a0*h1_TP1*V1t_0_[i];
        df2dt_TP1_[i] =  0.5*a1*h2_TP1*V1t_1_[i];
        df3dt_TP1_[i] =  0.5*a2*h3_TP1*V1t_2_[i];
        dk1dt_TP2_[i] = -0.5*a0*h1_TP2*V2t_0_[i];
        dk2dt_TP2_[i] = -0.5*a1*h2_TP2*V2t_1_[i];
        dk3dt_TP2_[i] = -0.5*a2*h3_TP2*V2t_2_[i];
        df1dt_TP2_[i] =  0.5*a0*h1_TP2*V1t_0_[i];
        df2dt_TP2_[i] =  0.5*a1*h2_TP2*V1t_1_[i];
        df3dt_TP2_[i] =  0.5*a2*h3_TP2*V1t_2_[i];
        dk1dt_TP3_[i] = -0.5*a0*h1_TP3*V2t_0_[i];
        dk2dt_TP3_[i] = -0.5*a1*h2_TP3*V2t_1_[i];
        dk3dt_TP3_[i] = -0.5*a2*h3_TP3*V2t_2_[i];
        df1dt_TP3_[i] =  0.5*a0*h1_TP3*V1t_0_[i];
        df2dt_TP3_[i] =  0.5*a1*h2_TP3*V1t_1_[i];
        df3dt_TP3_[i] =  0.5*a2*h3_TP3*V1t_2_[i];

      }

    //NOTE: ert = ert(1) + cs, est = est(2) - cr
    //this is just ert(1), ert(3), est(2), est(3)
    //translational components
    for (int i=0;i<nnodes_loc;++i)
      {
        for (int j=0;j<3;++j)//number of translation degree of freedom
          {
            int jj = j+ndof_loc*i;
            //est2
            est2_[jj] = 0.5*(dhmat_[1][i]*gt_3_TP2_[j] + dhmat_[2][i]*gt_2_[j]); 
            //ert1
            ert1_[jj] = 0.5*(dhmat_[0][i]*gt_3_TP1_[j] + dhmat_[2][i]*gt_1_[j]);         
            //est3
            est3_[jj] = 0.5*(dhmat_[1][i]*gt_3_TP3_[j] + dhmat_[2][i]*gt_2_[j]); 
            //ert3
            ert3_[jj] = 0.5*(dhmat_[0][i]*gt_3_TP3_[j] + dhmat_[2][i]*gt_1_[j]); 
          }
      }
    //rotational components
    //alpha 1
    est2_[3]  =  0.5*(dprod(dk1ds_,gt_3_TP2_) + dprod(dk1dt_TP2_,gt_2_));
    ert1_[3]  =  0.5*(dprod(dk1dr_,gt_3_TP1_) + dprod(dk1dt_TP1_,gt_1_));
    est3_[3]  =  0.5*(dprod(dk1ds_,gt_3_TP3_) + dprod(dk1dt_TP3_,gt_2_));
    ert3_[3]  =  0.5*(dprod(dk1dr_,gt_3_TP3_) + dprod(dk1dt_TP3_,gt_1_));
    //beta 1
    est2_[4]  =  0.5*(dprod(df1ds_,gt_3_TP2_) + dprod(df1dt_TP2_,gt_2_)); 
    ert1_[4]  =  0.5*(dprod(df1dr_,gt_3_TP1_) + dprod(df1dt_TP1_,gt_1_)); 
    est3_[4]  =  0.5*(dprod(df1ds_,gt_3_TP3_) + dprod(df1dt_TP3_,gt_2_)); 
    ert3_[4]  =  0.5*(dprod(df1dr_,gt_3_TP3_) + dprod(df1dt_TP3_,gt_1_)); 
    //alpha 2
    est2_[8]  =  0.5*(dprod(dk2ds_,gt_3_TP2_) + dprod(dk2dt_TP2_,gt_2_)); 
    ert1_[8]  =  0.5*(dprod(dk2dr_,gt_3_TP1_) + dprod(dk2dt_TP1_,gt_1_)); 
    est3_[8]  =  0.5*(dprod(dk2ds_,gt_3_TP3_) + dprod(dk2dt_TP3_,gt_2_)); 
    ert3_[8]  =  0.5*(dprod(dk2dr_,gt_3_TP3_) + dprod(dk2dt_TP3_,gt_1_)); 
    //beta 2
    est2_[9]  =  0.5*(dprod(df2ds_,gt_3_TP2_) + dprod(df2dt_TP2_,gt_2_)); 
    ert1_[9]  =  0.5*(dprod(df2dr_,gt_3_TP1_) + dprod(df2dt_TP1_,gt_1_)); 
    est3_[9]  =  0.5*(dprod(df2ds_,gt_3_TP3_) + dprod(df2dt_TP3_,gt_2_)); 
    ert3_[9]  =  0.5*(dprod(df2dr_,gt_3_TP3_) + dprod(df2dt_TP3_,gt_1_)); 
    //alpha 3
    est2_[13] =  0.5*(dprod(dk3ds_,gt_3_TP2_) + dprod(dk3dt_TP2_,gt_2_)); 
    ert1_[13] =  0.5*(dprod(dk3dr_,gt_3_TP1_) + dprod(dk3dt_TP1_,gt_1_)); 
    est3_[13] =  0.5*(dprod(dk3ds_,gt_3_TP3_) + dprod(dk3dt_TP3_,gt_2_)); 
    ert3_[13] =  0.5*(dprod(dk3dr_,gt_3_TP3_) + dprod(dk3dt_TP3_,gt_1_)); 
    //beta 3
    est2_[14] =  0.5*(dprod(df3ds_,gt_3_TP2_) + dprod(df3dt_TP2_,gt_2_)); 
    ert1_[14] =  0.5*(dprod(df3dr_,gt_3_TP1_) + dprod(df3dt_TP1_,gt_1_)); 
    est3_[14] =  0.5*(dprod(df3ds_,gt_3_TP3_) + dprod(df3dt_TP3_,gt_2_)); 
    ert3_[14] =  0.5*(dprod(df3dr_,gt_3_TP3_) + dprod(df3dt_TP3_,gt_1_));            

    //create c = est(2) - ert(1) - est(3) + ert(3)     
    for (int i=0;i<nnodes_loc*ndof_loc;++i) cTying_[i]  = est2_[i] - ert1_[i] - est3_[i] + ert3_[i];

    //create ert = ert(1) + cs and est = est(2) - cr
    for (int i=0;i<nnodes_loc*ndof_loc;++i)
      {

        B_rst_[4][i]  = est2_[i] - cTying_[i]*rPoint;
        B_rst_[5][i]  = ert1_[i] + cTying_[i]*sPoint;

      }

    //NO CORRECTION
    //B_rst_[4][3]   = 0.5*(dprod(dk1ds,gt_matT[2]) + dprod(dk1dt,gt_matT[1]));
    //B_rst_[5][3]   = 0.5*(dprod(dk1dr,gt_matT[2]) + dprod(dk1dt,gt_matT[0]));
    //B_rst_[4][9]   = 0.5*(dprod(dk2ds,gt_matT[2]) + dprod(dk2dt,gt_matT[1]));
    //B_rst_[5][9]   = 0.5*(dprod(dk2dr,gt_matT[2]) + dprod(dk2dt,gt_matT[0]));
    //B_rst_[4][15]  = 0.5*(dprod(dk3ds,gt_matT[2]) + dprod(dk3dt,gt_matT[1]));
    //B_rst_[5][15]  = 0.5*(dprod(dk3dr,gt_matT[2]) + dprod(dk3dt,gt_matT[0]));
    //B_rst_[4][4]   = 0.5*(dprod(df1ds,gt_matT[2]) + dprod(df1dt,gt_matT[1]));
    //B_rst_[5][4]   = 0.5*(dprod(df1dr,gt_matT[2]) + dprod(df1dt,gt_matT[0]));
    //B_rst_[4][10]  = 0.5*(dprod(df2ds,gt_matT[2]) + dprod(df2dt,gt_matT[1]));
    //B_rst_[5][10]  = 0.5*(dprod(df2dr,gt_matT[2]) + dprod(df2dt,gt_matT[0]));
    //B_rst_[4][16]  = 0.5*(dprod(df3ds,gt_matT[2]) + dprod(df3dt,gt_matT[1]));
    //B_rst_[5][16]  = 0.5*(dprod(df3dr,gt_matT[2]) + dprod(df3dt,gt_matT[0]));

    //zero out last workingvec
    for (int k=0;k<nnodes_loc*ndof_loc;++k) for (int i=0;i<3;++i) for (int j=0;j<3;++j) workingvectmp2_[k][i][j] = 0.;

    //transform B_rst_ to global system witth tensor transform
    //create emn = [err,ers,ert,esr,ess,est,etr,ets,ett]
    for (int i=0;i<nnodes_loc*ndof_loc;++i) 
      { 
        workingvectmp1_[i][0][0] = B_rst_[0][i]; //err
        workingvectmp1_[i][0][1] = B_rst_[3][i]; //ers
        workingvectmp1_[i][0][2] = B_rst_[5][i]; //ert
        workingvectmp1_[i][1][0] = B_rst_[3][i]; //esr
        workingvectmp1_[i][1][1] = B_rst_[1][i]; //ess
        workingvectmp1_[i][1][2] = B_rst_[4][i]; //est
        workingvectmp1_[i][2][0] = B_rst_[5][i]; //etr
        workingvectmp1_[i][2][1] = B_rst_[4][i]; //ets
        workingvectmp1_[i][2][2] = B_rst_[2][i]; //ett
      }

    //bring B_rst_ to global system with tensor transform
    for (int k=0;k<nnodes_loc*ndof_loc;++k) 
      {
        int count = 0;
        for (int i=0;i<3;++i)
          {
            for (int j=0;j<3;++j)
              {
                for (int m=0;m<3;++m)
                  {
                    for (int n=0;n<3;++n)
                      {
                        workingvectmp2_[k][i][j] += workingvectmp1_[k][m][n]*tensorprod_[faceID][count];
                        count += 1;
                      }
                  }
              }
          }
      }

    for (int i=0;i<nnodes_loc*ndof_loc;++i)
      {
        B_xyz_[0][i] = workingvectmp2_[i][0][0];
        B_xyz_[1][i] = workingvectmp2_[i][1][1];
        B_xyz_[2][i] = workingvectmp2_[i][2][2];
        B_xyz_[3][i] = 2.*workingvectmp2_[i][0][1];
        B_xyz_[4][i] = 2.*workingvectmp2_[i][1][2];
        B_xyz_[5][i] = 2.*workingvectmp2_[i][0][2];
      }

    //jedit. 9/5/19: i still want to unpack this multiplication and remove the transpose variable
    matrixTranspose(B_xyz_,B_xyzT_);
    mm_mult(B_xyzT_,CmatNL_[faceID],B_xyzTC_);
    mm_mult(B_xyzTC_,B_xyz_,BTCBmat_);
    //set drilling dof to 10E-10
    //BTCBmat_[5][5]   = 0.00000000001;
    //BTCBmat_[11][11] = 0.00000000001;
    //BTCBmat_[17][17] = 0.00000000001;
    //compute Fvec local for rhs Newton-Raphson
    mTv_mult(B_xyzT_,PKSTvec_,Fvec_loc_);

  }

  void linearTriElement::calculateStrainsMITC3(double& rPoint,
                                               double& sPoint,
                                               double& tPoint,
                                               double& a0,
                                               double& a1,
                                               double& a2,
                                               const int& faceID)
  {

    int nnodes_loc = elementData_[faceID].nnodes;
    int ndof_loc   = elementData_[faceID].ndof;

    //2D shape functions
    double h1 = 1. -rPoint - sPoint;
    double h2 = rPoint;
    double h3 = sPoint;
    //2D shape function derivatives
    //h1 = 1-r-s
    double dh1dr = -1.;
    double dh1ds = -1.;
    double dh1dt =  0.;
    //h2 = r
    double dh2dr =  1.;
    double dh2ds =  0.;
    double dh2dt =  0.;
    //h3 = s
    double dh3dr =  0.;
    double dh3ds =  1.;
    double dh3dt =  0.;

    //create total change in rotation vector 
    for (int idir=0;idir<3;++idir) 
      { 
        deltaVnt0_0_[idir] = Vnt_0_[idir] - Vn0_0_[idir];
        deltaVnt0_1_[idir] = Vnt_1_[idir] - Vn0_1_[idir];
        deltaVnt0_2_[idir] = Vnt_2_[idir] - Vn0_2_[idir];
      }

    //create ut,i = dut/dri for updating gti = g0_i + ut,i
    //where ut = xt - x0
    //uti = ut/dri
    for (int i=0;i<3;++i)
      {
        ut1_[i] = dh1dr*(x0_t_[i]-x0_0_[i]) + dh2dr*(x1_t_[i]-x1_0_[i]) + dh3dr*(x2_t_[i]-x2_0_[i]) + 0.5*tPoint*(a0*dh1dr*deltaVnt0_0_[i]+a1*dh2dr*deltaVnt0_1_[i]+a2*dh3dr*deltaVnt0_2_[i]);
        ut2_[i] = dh1ds*(x0_t_[i]-x0_0_[i]) + dh2ds*(x1_t_[i]-x1_0_[i]) + dh3ds*(x2_t_[i]-x2_0_[i]) + 0.5*tPoint*(a0*dh1ds*deltaVnt0_0_[i]+a1*dh2ds*deltaVnt0_1_[i]+a2*dh3ds*deltaVnt0_2_[i]);
        ut3_[i] = 0.5*(a0*h1*deltaVnt0_0_[i] + a1*h2*deltaVnt0_1_[i] + a2*h3*deltaVnt0_2_[i]);
        //covariant t set
        gt_1_[i] = g0_1_[i] + ut1_[i];
        gt_2_[i] = g0_2_[i] + ut2_[i];
        gt_3_[i] = g0_3_[i] + ut3_[i];
        //create covariant 't' matrix
        gt_matT_[0][i] = gt_1_[i];
        gt_matT_[1][i] = gt_2_[i];
        gt_matT_[2][i] = gt_3_[i];
      }

    //calculate Green-Langrange strains in (r,s,t) -> eij = gti*gtj - g0i*g0j
    //in-plane components
    GLstrains_[0] = 0.5*(dprod(gt_1_,gt_1_) - dprod(g0_1_,g0_1_)); //err
    GLstrains_[1] = 0.5*(dprod(gt_2_,gt_2_) - dprod(g0_2_,g0_2_)); //ess
    GLstrains_[2] = 0.5*(dprod(gt_3_,gt_3_) - dprod(g0_3_,g0_3_)); //ett
    GLstrains_[3] = 0.5*(dprod(gt_1_,gt_2_) - dprod(g0_1_,g0_2_)); //ers

    //no corrections
    //GLstrains[4] = 0.5*(dprod(gt_2,gt_3) - dprod(g0_2,g0_3)); //ers
    //GLstrains[5] = 0.5*(dprod(gt_1,gt_3) - dprod(g0_1,g0_3)); //ers

    //implement MITC tying scheme for est and ert transverse shear strains
    //NOTE: only ut3 and g0_3 change with r,s at the tying points
    //tying points
    double r1 = 0.5;
    double s1 = 0.;
    double r2 = 0.;
    double s2 = 0.5;
    double r3 = 0.5;
    double s3 = 0.5;
    //shape functions at tying points
    double h1_TP1 = 1.-r1-s1;
    double h2_TP1 = r1;
    double h3_TP1 = s1;
    double h1_TP2 = 1.-r2-s2;
    double h2_TP2 = r2;
    double h3_TP2 = s2;
    double h1_TP3 = 1.-r3-s3;
    double h2_TP3 = r3;
    double h3_TP3 = s3;

    for (int i=0;i<3;++i)
      {
        //ut3 at tying points
        ut3_TP1_[i] = 0.5*(a0*h1_TP1*deltaVnt0_0_[i] + a1*h2_TP1*deltaVnt0_1_[i] + a2*h3_TP1*deltaVnt0_2_[i]);
        ut3_TP2_[i] = 0.5*(a0*h1_TP2*deltaVnt0_0_[i] + a1*h2_TP2*deltaVnt0_1_[i] + a2*h3_TP2*deltaVnt0_2_[i]);
        ut3_TP3_[i] = 0.5*(a0*h1_TP3*deltaVnt0_0_[i] + a1*h2_TP3*deltaVnt0_1_[i] + a2*h3_TP3*deltaVnt0_2_[i]);
        //g0_3 at tying points
        g0_3_TP1_[i] = 0.5*(a0*h1_TP1*Vn0_0_[i] + a1*h2_TP1*Vn0_1_[i] + a2*h3_TP1*Vn0_2_[i]); 
        g0_3_TP2_[i] = 0.5*(a0*h1_TP2*Vn0_0_[i] + a1*h2_TP2*Vn0_1_[i] + a2*h3_TP2*Vn0_2_[i]); 
        g0_3_TP3_[i] = 0.5*(a0*h1_TP3*Vn0_0_[i] + a1*h2_TP3*Vn0_1_[i] + a2*h3_TP3*Vn0_2_[i]); 

      }    
    //covariant 't' set at tying points
    for (int i=0;i<3;++i) 
      {
        gt_3_TP1_[i] = g0_3_TP1_[i] + ut3_TP1_[i];
        gt_3_TP2_[i] = g0_3_TP2_[i] + ut3_TP2_[i];
        gt_3_TP3_[i] = g0_3_TP3_[i] + ut3_TP3_[i];
      }

    //transverse strains at tying points
    double est_TP1;
    double est_TP2;
    double est_TP3;
    double ert_TP1;
    double ert_TP2;
    double ert_TP3;

    est_TP1 = 0.5*(dprod(gt_2_,gt_3_TP1_) - dprod(g0_2_,g0_3_TP1_));
    est_TP2 = 0.5*(dprod(gt_2_,gt_3_TP2_) - dprod(g0_2_,g0_3_TP2_));
    est_TP3 = 0.5*(dprod(gt_2_,gt_3_TP3_) - dprod(g0_2_,g0_3_TP3_));

    ert_TP1 = 0.5*(dprod(gt_1_,gt_3_TP1_) - dprod(g0_1_,g0_3_TP1_));
    ert_TP2 = 0.5*(dprod(gt_1_,gt_3_TP2_) - dprod(g0_1_,g0_3_TP2_));
    ert_TP3 = 0.5*(dprod(gt_1_,gt_3_TP3_) - dprod(g0_1_,g0_3_TP3_));

    //fill transverse shear strains with MITC interpolation
    double cConst = est_TP2 - ert_TP1 - est_TP3 + ert_TP3;
    GLstrains_[4]  = est_TP2 - cConst*rPoint;
    GLstrains_[5]  = ert_TP1 + cConst*sPoint;

    //dumpVec(GLstrains,"GL strain vector (after correction)");

    //reorganize GLstrains_ into 3x3 matrix
    strains_rst_[0][0] = GLstrains_[0];       //err
    strains_rst_[0][1] = GLstrains_[3];       //ers
    strains_rst_[0][2] = GLstrains_[5];       //ert
    strains_rst_[1][0] = strains_rst_[0][1];  //esr
    strains_rst_[1][1] = GLstrains_[1];       //ess
    strains_rst_[1][2] = GLstrains_[4];       //est
    strains_rst_[2][0] = strains_rst_[0][2];  //etr 
    strains_rst_[2][1] = strains_rst_[1][2];  //ets
    strains_rst_[2][2] = GLstrains_[2];       //ett

    //contravariant ZERO set
    double dum;
    dummy1_ = cross(g0_2_,g0_3_);
    dum = dprod(g0_1_,dummy1_);
    for (int i=0;i<3;++i) g01_[i] = dummy1_[i]/dum;

    dummy1_ = cross(g0_3_,g0_1_);
    dum = dprod(g0_2_,dummy1_);
    for (int i=0;i<3;++i) g02_[i] = dummy1_[i]/dum;

    dummy1_ = cross(g0_1_,g0_2_);
    dum = dprod(g0_3_,dummy1_);
    for (int i=0;i<3;++i) g03_[i] = dummy1_[i]/dum;

    //create contravariant ZERO matrix
    for (int i=0;i<DIM;++i) 
      {
        g0matT_[0][i] = g01_[i];
        g0matT_[1][i] = g02_[i];
        g0matT_[2][i] = g03_[i];
      }

    int count = 0;
    //store tensor products to reference configuration once for efficiency
    if (firstTimeComputingK_ == true)
      {
        for (int i=0;i<3;++i)
          {
            for (int j=0;j<3;++j)
              {
                for (int m=0;m<3;++m)
                  {
                    for (int n=0;n<3;++n)
                      {
                        tensorprod_[faceID][count] = dprod(ematT_[i],g0matT_[m])*dprod(ematT_[j],g0matT_[n]);
                        count += 1;
                      }
                  }
              }
          }
      }
    //transform r,s,t strains into global system
    count = 0;
    for (int i=0;i<3;++i) 
      {
        for (int j=0;j<3;++j) 
          { 
            GLstrainMat_[i][j] = 0.;
            for (int m=0;m<3;++m) 
              {
                for (int n=0;n<3;++n) 
                  {
                    GLstrainMat_[i][j] += strains_rst_[m][n]*tensorprod_[faceID][count];
                    count += 1;
                  }
              }
          }
      }
    //gamma = 2*epsilon
    GLstrainMat_[0][1] = 2.*GLstrainMat_[0][1];
    GLstrainMat_[1][0] = 2.*GLstrainMat_[1][0];
    GLstrainMat_[0][2] = 2.*GLstrainMat_[0][2];
    GLstrainMat_[2][0] = 2.*GLstrainMat_[2][0];
    GLstrainMat_[1][2] = 2.*GLstrainMat_[1][2];
    GLstrainMat_[2][1] = 2.*GLstrainMat_[2][1];
    //refill GLstrains with transformed GLstrainMat_
    GLstrains_[0] = GLstrainMat_[0][0]; //err
    GLstrains_[1] = GLstrainMat_[1][1]; //ess
    GLstrains_[2] = GLstrainMat_[2][2]; //ett
    GLstrains_[3] = GLstrainMat_[0][1]; //ers
    GLstrains_[4] = GLstrainMat_[1][2]; //est
    GLstrains_[5] = GLstrainMat_[0][2]; //ert

  }
    
  void linearTriElement::compPKST(const int& faceID)
  {

    ////calculate Second Piola Kirchoff Stress Tensor from Green-Lagrange Strains Matrix
    //USING EQ 6.184 IN FEP - Jedit 11/8/17
    mTv_mult(CmatNL_[faceID],GLstrains_,PKSTvec_);

    PKSTmat_[0][0] = PKSTvec_[0];
    PKSTmat_[1][1] = PKSTvec_[1];
    PKSTmat_[2][2] = PKSTvec_[2];
    PKSTmat_[0][1] = PKSTvec_[3];
    PKSTmat_[1][2] = PKSTvec_[4];
    PKSTmat_[0][2] = PKSTvec_[5];
    PKSTmat_[2][0] = PKSTmat_[0][2];
    PKSTmat_[2][1] = PKSTmat_[1][2];
    PKSTmat_[1][0] = PKSTmat_[0][1];

  }

  void linearTriElement::calculateDeformationGradient(double& rPoint,
						      double& sPoint,
						      double& tPoint,
						      double& a0,
						      double& a1,
						      double& a2,
						      const int& faceID)
  {

    //calculate the deformation gradient for the MITC3 shell element, X = d(tx)/d(0x) = d(tu) + I

    int nnodes_loc = elementData_[faceID].nnodes;
    int ndof_loc   = elementData_[faceID].ndof;

    //2D shape functions
    double h1 = 1. -rPoint - sPoint;
    double h2 = rPoint;
    double h3 = sPoint;
    //2D shape function derivatives
    //h1 = 1-r-s
    double dh1dr = -1.;
    double dh1ds = -1.;
    double dh1dt =  0.;
    //h2 = r
    double dh2dr =  1.;
    double dh2ds =  0.;
    double dh2dt =  0.;
    //h3 = s
    double dh3dr =  0.;
    double dh3ds =  1.;
    double dh3dt =  0.;

    deformationGradient_rst_[0][0] = dh1dr*(x0_t_[0]-x0_0_[0]) + dh2dr*(x1_t_[0]-x1_0_[0]) + dh3dr*(x2_t_[0]-x2_0_[0]) + 0.5*tPoint*(a0*dh1dr*deltaVnt0_0_[0]+a1*dh2dr*deltaVnt0_1_[0]+a2*dh3dr*deltaVnt0_2_[0]) + 1.;
    deformationGradient_rst_[0][1] = dh1ds*(x0_t_[0]-x0_0_[0]) + dh2ds*(x1_t_[0]-x1_0_[0]) + dh3ds*(x2_t_[0]-x2_0_[0]) + 0.5*tPoint*(a0*dh1ds*deltaVnt0_0_[0]+a1*dh2ds*deltaVnt0_1_[0]+a2*dh3ds*deltaVnt0_2_[0]);
    deformationGradient_rst_[0][2] = 0.5*(a0*h1*deltaVnt0_0_[0] + a1*h2*deltaVnt0_1_[0] + a2*h3*deltaVnt0_2_[0]);

    deformationGradient_rst_[1][0] = dh1dr*(x0_t_[1]-x0_0_[1]) + dh2dr*(x1_t_[1]-x1_0_[1]) + dh3dr*(x2_t_[1]-x2_0_[1]) + 0.5*tPoint*(a0*dh1dr*deltaVnt0_0_[1]+a1*dh2dr*deltaVnt0_1_[1]+a2*dh3dr*deltaVnt0_2_[1]);
    deformationGradient_rst_[1][1] = dh1ds*(x0_t_[1]-x0_0_[1]) + dh2ds*(x1_t_[1]-x1_0_[1]) + dh3ds*(x2_t_[1]-x2_0_[1]) + 0.5*tPoint*(a0*dh1ds*deltaVnt0_0_[1]+a1*dh2ds*deltaVnt0_1_[1]+a2*dh3ds*deltaVnt0_2_[1]) + 1.;
    deformationGradient_rst_[1][2] = 0.5*(a0*h1*deltaVnt0_0_[1] + a1*h2*deltaVnt0_1_[1] + a2*h3*deltaVnt0_2_[1]);

    deformationGradient_rst_[2][0] = dh1dr*(x0_t_[2]-x0_0_[2]) + dh2dr*(x1_t_[2]-x1_0_[2]) + dh3dr*(x2_t_[2]-x2_0_[2]) + 0.5*tPoint*(a0*dh1dr*deltaVnt0_0_[2]+a1*dh2dr*deltaVnt0_1_[2]+a2*dh3dr*deltaVnt0_2_[2]);
    deformationGradient_rst_[2][1] = dh1ds*(x0_t_[2]-x0_0_[2]) + dh2ds*(x1_t_[2]-x1_0_[2]) + dh3ds*(x2_t_[2]-x2_0_[2]) + 0.5*tPoint*(a0*dh1ds*deltaVnt0_0_[2]+a1*dh2ds*deltaVnt0_1_[2]+a2*dh3ds*deltaVnt0_2_[2]);
    deformationGradient_rst_[2][2] = 0.5*(a0*h1*deltaVnt0_0_[2] + a1*h2*deltaVnt0_1_[2] + a2*h3*deltaVnt0_2_[2]) + 1.;

    //bring deformation gradient to global Cartesian
    int count = 0;
    for (int i=0;i<3;++i) 
      {
        for (int j=0;j<3;++j) 
          { 
            deformationGradient_xyz_[i][j] = 0.;
            for (int m=0;m<3;++m) 
              {
                for (int n=0;n<3;++n) 
                  {
                    //deformationGradient_xyz_[i][j] += deformationGradient_rst_[m][n]*tensorprod_[faceID][count];
		    deformationGradient_xyz_[i][j] = deformationGradient_rst_[i][j];
                    count += 1;
                  }
              }
          }
      }

    //jacobian is det(deformationGradient_xyz_)
    jacobian_ = deformationGradient_xyz_[0][0]*(deformationGradient_xyz_[1][1]*deformationGradient_xyz_[2][2] - deformationGradient_xyz_[2][1]*deformationGradient_xyz_[1][2]) -
                deformationGradient_xyz_[0][1]*(deformationGradient_xyz_[1][0]*deformationGradient_xyz_[2][2] - deformationGradient_xyz_[1][2]*deformationGradient_xyz_[2][0]) +
                deformationGradient_xyz_[0][2]*(deformationGradient_xyz_[1][0]*deformationGradient_xyz_[2][1] - deformationGradient_xyz_[1][1]*deformationGradient_xyz_[2][0]);

    inv_jacobian_ = 1./jacobian_;

  }
  
  void linearTriElement::PKST2Cauchy(const int& faceID)
  {

    //turn 2nd Piola-Kirchhoff stress tensor into Cauchy stress tensor via deformation gradient, Cauchy = 1/J * X * PKST * XT

    for (int i=0;i<3;++i)
      {
        for (int j=0;j<3;++j)
          {
            integratingCauchyST_2_[i][j] = 0.;
            for (int k=0;k<3;++k)
              {
                integratingCauchyST_2_[i][j] += deformationGradient_xyz_[i][k] * PKSTmat_[k][j];
              }
          }
      }

    for (int i=0;i<3;++i)
      {
        for (int j=0;j<3;++j)
          {
            integratingCauchyST_1_[i][j] = 0.;
            for (int k=0;k<3;++k)
              {
                integratingCauchyST_1_[i][j] += integratingCauchyST_2_[i][k] * deformationGradient_xyz_[j][k];
              }
          }
      }
    //bring deformation gradient to global Cartesian
    int count = 0;
    for (int i=0;i<3;++i)
      {
        for (int j=0;j<3;++j)
          {
            integratingCauchyST_2_[i][j] = 0.;
            for (int m=0;m<3;++m)
              {
                for (int n=0;n<3;++n)
                  {
                    integratingCauchyST_2_[i][j] += integratingCauchyST_1_[m][n]*tensorprod_[faceID][count];
                    count += 1;
                  }
              }
          }
      }

    for (int i=0;i<3;++i) for (int j=0;j<3;++j) integratingCauchyST_2_[i][j] *= inv_jacobian_;

  }
  
  void linearTriElement::calculateVonMisesStress(const int& faceID)
  {

    //calculate von Mises using deviatoric stress tensor VM = sqrt(3/2 devi : devi)

    double hydrostaticStress = 0.3333333333333*(cauchyStressTensor_[faceID][0][0] + cauchyStressTensor_[faceID][1][1] + cauchyStressTensor_[faceID][2][2]);

    vonMisesStress_[faceID] = 0.;
    for (int i=0;i<3;++i)
      {
	for (int j=0;j<3;++j)
	  {
	    vonMisesStress_[faceID] += 1.5 * (cauchyStressTensor_[faceID][i][j] - hydrostaticStress) * (cauchyStressTensor_[faceID][i][j] - hydrostaticStress);
	  }
      }

    vonMisesStress_[faceID] = pow(vonMisesStress_[faceID],0.5);

  }

  void linearTriElement::comp_Me_MITC3(vector<double>& Me_loc,
                                       const int&      faceID)
  {
    //computing local mass matrix for MITC3 triangular shell element

    int nnodes_loc   = elementData_[faceID].nnodes;
    int ndof_loc     = elementData_[faceID].ndof;
    int nnodesOffset = nnodesLocalOffset_[faceID];

    for (int n=0;n<nnodes_loc;++n) 
      {
        int index = n + nnodesOffset;
        face_nodes_[n] = faces_[index];
      }
    //face nodes
    int fn0 = face_nodes_[0];
    int fn1 = face_nodes_[1];
    int fn2 = face_nodes_[2];
    //nodal positions
    for (int idir=0;idir<xyz;++idir) 
      {
        x0_0_[idir] = nodes0_[fn0*xyz+idir];
        x1_0_[idir] = nodes0_[fn1*xyz+idir];
        x2_0_[idir] = nodes0_[fn2*xyz+idir];
      }
    //integration points
    //r-s plane
    vector<double> rPoint(3,0.);
    vector<double> sPoint(3,0.);
    rPoint[0] =  0.1666666666667;
    rPoint[1] =  0.6666666666667;
    rPoint[2] =  0.1666666666667;
    sPoint[0] =  0.1666666666667;
    sPoint[1] =  0.1666666666667;
    sPoint[2] =  0.6666666666667;
    vector<double> rs_weights(3,0.3333333333333);    
    //t-plane
    vector<double> tPoint(2,0.);
    tPoint[0] = -0.5773502691896257;
    tPoint[1] =  0.5773502691896257;
    vector<double> t_weights(2,0.5);
    //thickness 
    double a0 = thickness_[faceID];
    double a1 = thickness_[faceID];
    double a2 = thickness_[faceID];

    //assemble local mass matrix
    for (int k=0;k<rPoint.size();++k) 
      {
        //create 2D shape functions
        double h1 = 1. - rPoint[k] - sPoint[k];
        double h2 = rPoint[k];
        double h3 = sPoint[k];
        //fill displacement interpolation matrix for translational components
        N_[0][0] = h1;
        N_[0][3] = h2;
        N_[0][6] = h3;
        N_[1][1] = h1;
        N_[1][4] = h2;
        N_[1][7] = h3;
        N_[2][2] = h1;
        N_[2][5] = h2;
        N_[2][8] = h3;
        //transpose
        matrixTranspose(N_,NT_);
        //NT * N
        mm_mult(NT_,N_,NTN_);

        for (int kk=0;kk<tPoint.size();++kk) 
          {
            //fill displacement interpolation matrix for rotational components
            for (int i=0;i<3;++i)
              {
                N_rot_[i][0] = -0.5*tPoint[kk]*a0*h1*V20_[myLocalNodesNormalVectorStatus_[fn0]*xyz+i];
                N_rot_[i][1] =  0.5*tPoint[kk]*a0*h1*V10_[myLocalNodesNormalVectorStatus_[fn0]*xyz+i];
                N_rot_[i][2] = -0.5*tPoint[kk]*a1*h2*V20_[myLocalNodesNormalVectorStatus_[fn1]*xyz+i];
                N_rot_[i][3] =  0.5*tPoint[kk]*a1*h2*V10_[myLocalNodesNormalVectorStatus_[fn1]*xyz+i];
                N_rot_[i][4] = -0.5*tPoint[kk]*a2*h3*V20_[myLocalNodesNormalVectorStatus_[fn2]*xyz+i];
                N_rot_[i][5] =  0.5*tPoint[kk]*a2*h3*V10_[myLocalNodesNormalVectorStatus_[fn2]*xyz+i];
              }
            //transpose
            matrixTranspose(N_rot_,NT_rot_);
            //NT_rot * N_rot
            mm_mult(NT_rot_,N_rot_,NTN_rot_);

            //translational dof
            for (int r=0;r<3;++r) 
              {
                Me_loc[(r   )*ndof_loc*nnodes_loc+(r   )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_[r  ][r  ] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r   )*ndof_loc*nnodes_loc+(r+5 )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_[r  ][r+3] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r   )*ndof_loc*nnodes_loc+(r+10)] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_[r  ][r+6] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+5 )*ndof_loc*nnodes_loc+(r   )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_[r+3][r  ] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+5 )*ndof_loc*nnodes_loc+(r+5 )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_[r+3][r+3] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+5 )*ndof_loc*nnodes_loc+(r+10)] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_[r+3][r+6] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+10)*ndof_loc*nnodes_loc+(r   )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_[r+6][r  ] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+10)*ndof_loc*nnodes_loc+(r+5 )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_[r+6][r+3] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+10)*ndof_loc*nnodes_loc+(r+10)] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_[r+6][r+6] * thickness_[faceID] * Area_[faceID];
              }    
            //rotational dof
            for (int r=0;r<2;++r) 
              {
                Me_loc[(r+3 )*ndof_loc*nnodes_loc+(r+3 )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_rot_[r  ][r  ] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+3 )*ndof_loc*nnodes_loc+(r+8 )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_rot_[r  ][r+2] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+3 )*ndof_loc*nnodes_loc+(r+13)] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_rot_[r  ][r+4] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+8 )*ndof_loc*nnodes_loc+(r+3 )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_rot_[r+2][r  ] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+8 )*ndof_loc*nnodes_loc+(r+8 )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_rot_[r+2][r+2] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+8 )*ndof_loc*nnodes_loc+(r+13)] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_rot_[r+2][r+4] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+13)*ndof_loc*nnodes_loc+(r+3 )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_rot_[r+4][r  ] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+13)*ndof_loc*nnodes_loc+(r+8 )] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_rot_[r+4][r+2] * thickness_[faceID] * Area_[faceID];
                Me_loc[(r+13)*ndof_loc*nnodes_loc+(r+13)] += density_[faceID] * t_weights[kk] * rs_weights[k] * NTN_rot_[r+4][r+4] * thickness_[faceID] * Area_[faceID];
              }
          }
      }
    //fill drilling DOF
    //Me_loc[5 *ndof_loc*nnodes_loc+5 ] = 10.E-10;
    //Me_loc[11*ndof_loc*nnodes_loc+11] = 10.E-10;
    //Me_loc[17*ndof_loc*nnodes_loc+17] = 10.E-10;

  }

  void linearTriElement::comp_KNL1_MITC3(const int& faceID,
                                         double& rPoint,
                                         double& sPoint,
                                         double& tPoint)
  {
//    //KNL1 = ul/dri * ul/drj
//    //thickness 
//    double a0 = thickness_[faceID];
//    double a1 = thickness_[faceID];
//    double a2 = thickness_[faceID];
//    //create 2D shape function derivatives
//    //h1 = 1-r-s
//    double dh1dr = -1.;
//    double dh1ds = -1.;
//    double dh1dt =  0.;
//    //h2 = r
//    double dh2dr =  1.;
//    double dh2ds =  0.;
//    double dh2dt =  0.;
//    //h3 = s
//    double dh3dr =  0.;
//    double dh3ds =  1.;
//    double dh3dt =  0.;
//    //shape functions
//    double h1 = 1. - rPoint - sPoint;
//    double h2 = rPoint;
//    double h3 = sPoint;
//    vector<double> hmat(3,0.);
//    vector<double> hmat_TP1(3,0.);
//    vector<double> hmat_TP2(3,0.);
//    vector<double> hmat_TP3(3,0.);
//    //for tying points
//    double r1 = 0.5;
//    double s1 = 0.;
//    double r2 = 0.;
//    double s2 = 0.5;
//    double r3 = 0.5;
//    double s3 = 0.5;
//    //2D shape functions at tying points
//    double h1_TP1 = 1.-r1-s1;
//    double h2_TP1 = r1;
//    double h3_TP1 = s1;
//    double h1_TP2 = 1.-r2-s2;
//    double h2_TP2 = r2;
//    double h3_TP2 = s2;
//    double h1_TP3 = 1.-r3-s3;
//    double h2_TP3 = r3;
//    double h3_TP3 = s3;
//    hmat = {h1,h2,h3};
//    hmat_TP1 = {h1_TP1,h2_TP1,h3_TP1};
//    hmat_TP2 = {h1_TP2,h2_TP2,h3_TP2};
//    hmat_TP3 = {h1_TP3,h2_TP3,h3_TP3};
//    //translational DOF
//    BNL1_[0][0]  = dhmat_[0][0];
//    BNL1_[1][1]  = dhmat_[0][0];
//    BNL1_[2][2]  = dhmat_[0][0];
//    BNL1_[0][5]  = dhmat_[0][1];
//    BNL1_[1][6]  = dhmat_[0][1];
//    BNL1_[2][7]  = dhmat_[0][1];
//    BNL1_[0][10] = dhmat_[0][2];
//    BNL1_[1][11] = dhmat_[0][2];
//    BNL1_[2][12] = dhmat_[0][2];
//    BNL1_[3][0]  = dhmat_[1][0];
//    BNL1_[4][1]  = dhmat_[1][0];
//    BNL1_[5][2]  = dhmat_[1][0];
//    BNL1_[3][5]  = dhmat_[1][1];
//    BNL1_[4][6]  = dhmat_[1][1];
//    BNL1_[5][7]  = dhmat_[1][1];
//    BNL1_[3][10] = dhmat_[1][2];
//    BNL1_[4][11] = dhmat_[1][2];
//    BNL1_[5][12] = dhmat_[1][2];
//    BNL1_[6][0]  = dhmat_[2][0];
//    BNL1_[7][1]  = dhmat_[2][0];
//    BNL1_[8][2]  = dhmat_[2][0];
//    BNL1_[6][5]  = dhmat_[2][1];
//    BNL1_[7][6]  = dhmat_[2][1];
//    BNL1_[8][7]  = dhmat_[2][1];
//    BNL1_[6][10] = dhmat_[2][2];
//    BNL1_[7][11] = dhmat_[2][2];
//    BNL1_[8][12] = dhmat_[2][2];
//    vector<vector<double> > dudt_TP1(3,vector<double>(15,0.));
//    vector<vector<double> > dudt_TP2(3,vector<double>(15,0.));
//    vector<vector<double> > dudt_TP3(3,vector<double>(15,0.));
//    //rotational dof
//    for (int i=0;i<3;++i)
//      {
//        //alpha
//        BNL1_[i][3 ]     = -0.5*tPoint*a0*dhmat_[0][0]*V2t_0_[i];
//        BNL1_[i][8 ]     = -0.5*tPoint*a1*dhmat_[0][1]*V2t_1_[i];
//        BNL1_[i][13]     = -0.5*tPoint*a2*dhmat_[0][2]*V2t_2_[i];
//        BNL1_[i+3][3 ]   = -0.5*tPoint*a0*dhmat_[1][0]*V2t_0_[i];
//        BNL1_[i+3][8 ]   = -0.5*tPoint*a1*dhmat_[1][1]*V2t_1_[i];
//        BNL1_[i+3][13]   = -0.5*tPoint*a2*dhmat_[1][2]*V2t_2_[i];
//        BNL1_[i+6][3 ]   = -0.5*a0*hmat[0]*V2t_0_[i];
//        BNL1_[i+6][8 ]   = -0.5*a1*hmat[1]*V2t_1_[i];
//        BNL1_[i+6][13]   = -0.5*a2*hmat[2]*V2t_2_[i];
//        //corrections
//        dudt_TP1[i][3 ] = -0.5*a0*hmat_TP1[0]*V2t_0_[i];
//        dudt_TP1[i][8 ] = -0.5*a1*hmat_TP1[1]*V2t_1_[i];
//        dudt_TP1[i][13] = -0.5*a2*hmat_TP1[2]*V2t_2_[i];
//        dudt_TP2[i][3 ] = -0.5*a0*hmat_TP2[0]*V2t_0_[i];
//        dudt_TP2[i][8 ] = -0.5*a1*hmat_TP2[1]*V2t_1_[i];
//        dudt_TP2[i][13] = -0.5*a2*hmat_TP2[2]*V2t_2_[i];
//        dudt_TP3[i][3 ] = -0.5*a0*hmat_TP3[0]*V2t_0_[i];
//        dudt_TP3[i][8 ] = -0.5*a1*hmat_TP3[1]*V2t_1_[i];
//        dudt_TP3[i][13] = -0.5*a2*hmat_TP3[2]*V2t_2_[i];
//        //beta
//        BNL1_[i][4 ]     =  0.5*tPoint*a0*dhmat_[0][0]*V1t_0_[i];
//        BNL1_[i][9 ]     =  0.5*tPoint*a1*dhmat_[0][1]*V1t_1_[i];
//        BNL1_[i][14]     =  0.5*tPoint*a2*dhmat_[0][2]*V1t_2_[i];
//        BNL1_[i+3][4 ]   =  0.5*tPoint*a0*dhmat_[1][0]*V1t_0_[i];
//        BNL1_[i+3][9 ]   =  0.5*tPoint*a1*dhmat_[1][1]*V1t_1_[i];
//        BNL1_[i+3][14]   =  0.5*tPoint*a2*dhmat_[1][2]*V1t_2_[i];
//        BNL1_[i+6][4 ]   =  0.5*a0*hmat[0]*V1t_0_[i];
//        BNL1_[i+6][9 ]   =  0.5*a1*hmat[1]*V1t_1_[i];
//        BNL1_[i+6][14]   =  0.5*a2*hmat[2]*V1t_2_[i];
//        //corrections
//        dudt_TP1[i][4 ] =  0.5*a0*hmat_TP1[0]*V1t_0_[i];
//        dudt_TP1[i][9 ] =  0.5*a1*hmat_TP1[1]*V1t_1_[i];
//        dudt_TP1[i][14] =  0.5*a2*hmat_TP1[2]*V1t_2_[i];
//        dudt_TP2[i][4 ] =  0.5*a0*hmat_TP2[0]*V1t_0_[i];
//        dudt_TP2[i][9 ] =  0.5*a1*hmat_TP2[1]*V1t_1_[i];
//        dudt_TP2[i][14] =  0.5*a2*hmat_TP2[2]*V1t_2_[i];
//        dudt_TP3[i][4 ] =  0.5*a0*hmat_TP3[0]*V1t_0_[i];
//        dudt_TP3[i][9 ] =  0.5*a1*hmat_TP3[1]*V1t_1_[i];
//        dudt_TP3[i][14] =  0.5*a2*hmat_TP3[2]*V1t_2_[i];
//      }
//    vector<vector<double> > dudr(3,vector<double>(15,0.));
//    vector<vector<double> > duds(3,vector<double>(15,0.));
//    vector<vector<double> > dudt(3,vector<double>(15,0.));
//    vector<vector<double> > n11(15,vector<double>(15,0.));
//    vector<vector<double> > n12(15,vector<double>(15,0.));
//    vector<vector<double> > n13(15,vector<double>(15,0.));
//    vector<vector<double> > n13_TP1(15,vector<double>(15,0.));
//    vector<vector<double> > n13_TP2(15,vector<double>(15,0.));
//    vector<vector<double> > n13_TP3(15,vector<double>(15,0.));
//    vector<vector<double> > n21(15,vector<double>(15,0.));
//    vector<vector<double> > n22(15,vector<double>(15,0.));
//    vector<vector<double> > n23(15,vector<double>(15,0.));
//    vector<vector<double> > n23_TP1(15,vector<double>(15,0.));
//    vector<vector<double> > n23_TP2(15,vector<double>(15,0.));
//    vector<vector<double> > n23_TP3(15,vector<double>(15,0.));
//    vector<vector<double> > n31(15,vector<double>(15,0.));
//    vector<vector<double> > n31_TP1(15,vector<double>(15,0.));
//    vector<vector<double> > n31_TP2(15,vector<double>(15,0.));
//    vector<vector<double> > n31_TP3(15,vector<double>(15,0.));
//    vector<vector<double> > n32(15,vector<double>(15,0.));
//    vector<vector<double> > n32_TP1(15,vector<double>(15,0.));
//    vector<vector<double> > n32_TP2(15,vector<double>(15,0.));
//    vector<vector<double> > n32_TP3(15,vector<double>(15,0.));
//    vector<vector<double> > n33(15,vector<double>(15,0.));
//    vector<vector<double> > dudr_T(15,vector<double>(3,0.));
//    vector<vector<double> > duds_T(15,vector<double>(3,0.));
//    vector<vector<double> > dudt_T(15,vector<double>(3,0.));
//    vector<vector<double> > dudt_TP1_T(15,vector<double>(3,0.));
//    vector<vector<double> > dudt_TP2_T(15,vector<double>(3,0.));
//    vector<vector<double> > dudt_TP3_T(15,vector<double>(3,0.));
//    //decompose BNL1
//    for (int i=0;i<3;++i) 
//      {
//        for (int j=0;j<15;++j) 
//          {
//            dudr[i][j] = BNL1_[i  ][j];
//            duds[i][j] = BNL1_[i+3][j];
//            dudt[i][j] = BNL1_[i+6][j]; 
//          }
//      }
//    matrixTranspose(dudr,dudr_T);
//    matrixTranspose(duds,duds_T);
//    matrixTranspose(dudt,dudt_T);
//    matrixTranspose(dudt_TP1,dudt_TP1_T);
//    matrixTranspose(dudt_TP2,dudt_TP2_T);
//    matrixTranspose(dudt_TP3,dudt_TP3_T);
//    //create eta_ij
//    mm_mult(dudr_T,dudr,n11);
//    mm_mult(dudr_T,duds,n12);
//    mm_mult(duds_T,dudr,n21);
//    mm_mult(duds_T,duds,n22);
//    mm_mult(dudt_T,dudt,n33);
//    mm_mult(dudr_T,dudt_TP1,n13_TP1);
//    mm_mult(dudr_T,dudt_TP2,n13_TP2);
//    mm_mult(dudr_T,dudt_TP3,n13_TP3);
//    mm_mult(duds_T,dudt_TP1,n23_TP1);
//    mm_mult(duds_T,dudt_TP2,n23_TP2);
//    mm_mult(duds_T,dudt_TP3,n23_TP3);
//    mm_mult(dudt_TP1_T,dudr,n31_TP1);
//    mm_mult(dudt_TP2_T,dudr,n31_TP2);
//    mm_mult(dudt_TP3_T,dudr,n31_TP3);
//    mm_mult(dudt_TP1_T,duds,n32_TP1);
//    mm_mult(dudt_TP2_T,duds,n32_TP2);
//    mm_mult(dudt_TP3_T,duds,n32_TP3);
//    vector<vector<double> > cTyingF(15,vector<double>(15,0.));
//    vector<vector<double> > cTyingB(15,vector<double>(15,0.));
//    for (int i=0;i<15;++i)
//      {
//        for (int j=0;j<15;++j)
//          {
//            cTyingF[i][j] = n23_TP2[i][j] - n13_TP1[i][j] - n23_TP3[i][j] + n13_TP3[i][j];
//            cTyingB[i][j] = n32_TP2[i][j] - n31_TP1[i][j] - n32_TP3[i][j] + n31_TP3[i][j];
//            n13[i][j] = n13_TP1[i][j] + cTyingF[i][j]*sPoint;
//            n23[i][j] = n23_TP2[i][j] - cTyingF[i][j]*rPoint;
//            n31[i][j] = n31_TP1[i][j] + cTyingB[i][j]*sPoint;
//            n32[i][j] = n32_TP2[i][j] - cTyingB[i][j]*rPoint;
//          }
//      }
//    ////dumpMat(n13,"n13 with tying");
//    ////dumpMat(n23,"n23 with tying");
//    ////dumpMat(n31,"n31 with tying");
//    ////dumpMat(n32,"n32 with tying");
//    //no tying
//    //mm_mult(dudr_T,dudt,n13);
//    //mm_mult(duds_T,dudt,n23);
//    //mm_mult(dudt_T,dudr,n31);
//    //mm_mult(dudt_T,duds,n32);
//    ////dumpMat(n13,"n13");
//    ////dumpMat(n23,"n23");
//    ////dumpMat(n31,"n31");
//    ////dumpMat(n32,"n32");
//    //wait();
//    vector<vector<vector<vector<double> > > > nrs(3,vector<vector<vector<double> > >(3,vector<vector<double> >(15,vector<double>(15,0.))));
//    for (int i=0;i<15;++i) 
//      { 
//        for (int j=0;j<15;++j) 
//          { 
//            nrs[0][0][i][j] = n11[i][j]; //rr
//            nrs[0][1][i][j] = n12[i][j]; //rs
//            nrs[0][2][i][j] = n13[i][j]; //rt
//            nrs[1][0][i][j] = n21[i][j]; //sr
//            nrs[1][1][i][j] = n22[i][j]; //ss
//            nrs[1][2][i][j] = n23[i][j]; //st
//            nrs[2][0][i][j] = n31[i][j]; //tr
//            nrs[2][1][i][j] = n32[i][j]; //ts
//            nrs[2][2][i][j] = n33[i][j]; //tt
//          }
//      }
//    vector<vector<vector<vector<double> > > > nxy(3,vector<vector<vector<double> > >(3,vector<vector<double> >(15,vector<double>(15,0.))));
//    for (int i=0;i<15;++i) {
//      for (int j=0;j<15;++j) {                
//        for (int k=0;k<3;++k) {
//          for (int l=0;l<3;++l) {
//            for (int m=0;m<3;++m) {
//              for (int n=0;n<3;++n) {                
//                nxy[k][l][i][j] += nrs[m][n][i][j]*dprod(ematT_[k],g0matT_[m])*dprod(ematT_[l],g0matT_[n]);                
//              }
//            }            
//          }
//        }
//      }
//    }
//    //gamma = 2*epsilon
//    for (int i=0;i<15;++i) 
//      {
//      for (int j=0;j<15;++j) 
//        { 
//          nxy[0][1][i][j] = 1.*nxy[0][1][i][j];
//          nxy[1][0][i][j] = 1.*nxy[1][0][i][j];
//          nxy[0][2][i][j] = 1.*nxy[0][2][i][j];
//          nxy[2][0][i][j] = 1.*nxy[2][0][i][j];
//          nxy[1][2][i][j] = 1.*nxy[1][2][i][j];
//          nxy[2][1][i][j] = 1.*nxy[2][1][i][j];
//        }
//      }
//    vector<vector<double> > knl1_15(15,vector<double>(15,0.));  
//    for (int k=0;k<3;++k) {
//      for (int l=0;l<3;++l) {
//        for (int i=0;i<15;++i) {
//          for (int j=0;j<15;++j) {            
//            knl1_15[i][j] += nxy[k][l][i][j]*PKSTmat_[k][l];
//          }
//        }
//      } 
//    }
//    //fill KNL1 with zero in drilling dof
//    for (int i=0;i<5;++i) 
//      { 
//        for (int j=0;j<5;++j) 
//          { 
//            KNL1mat_[i   ][j   ] = knl1_15[i   ][j];
//            KNL1mat_[i+6 ][j   ] = knl1_15[i+5 ][j];
//            KNL1mat_[i+12][j   ] = knl1_15[i+10][j];
//            KNL1mat_[i   ][j+6 ] = knl1_15[i   ][j+5];
//            KNL1mat_[i+6 ][j+6 ] = knl1_15[i+5 ][j+5];
//            KNL1mat_[i+12][j+6 ] = knl1_15[i+10][j+5];
//            KNL1mat_[i   ][j+12] = knl1_15[i   ][j+10];
//            KNL1mat_[i+6 ][j+12] = knl1_15[i+5 ][j+10];
//            KNL1mat_[i+12][j+12] = knl1_15[i+10][j+10];
//          }
//      }
//    ////dumpMat(KNL1mat,"knl1mat me");

//    //create Jacobian
//    vector<vector<double> > Jacobian(3,vector<double>(3,0.));
//    vector<vector<double> > Jacobian_T(3,vector<double>(3,0.));
//    Jacobian_T[0][0] = x0_0[0]*dh1dr + x1_0[0]*dh2dr + x2_0[0]*dh3dr + 0.5*tPoint * (a0*dh1dr*Vn0_0[0] + a1*dh2dr*Vn0_1[0] + a2*dh3dr*Vn0_2[0]); //dx/dr
//    Jacobian_T[0][1] = x0_0[0]*dh1ds + x1_0[0]*dh2ds + x2_0[0]*dh3ds + 0.5*tPoint * (a0*dh1ds*Vn0_0[0] + a1*dh2ds*Vn0_1[0] + a2*dh3ds*Vn0_2[0]); //dx/ds
//    Jacobian_T[0][2] = 0.5 * (a0*h1*Vn0_0[0] + a1*h2*Vn0_1[0] + a2*h3*Vn0_2[0]); //dx/dt
//    Jacobian_T[1][0] = x0_0[1]*dh1dr + x1_0[1]*dh2dr + x2_0[1]*dh3dr + 0.5*tPoint * (a0*dh1dr*Vn0_0[1] + a1*dh2dr*Vn0_1[1] + a2*dh3dr*Vn0_2[1]); //dy/dr
//    Jacobian_T[1][1] = x0_0[1]*dh1ds + x1_0[1]*dh2ds + x2_0[1]*dh3ds + 0.5*tPoint * (a0*dh1ds*Vn0_0[1] + a1*dh2ds*Vn0_1[1] + a2*dh3ds*Vn0_2[1]); //dy/ds
//    Jacobian_T[1][2] = 0.5 * (a0*h1*Vn0_0[1] + a1*h2*Vn0_1[1] + a2*h3*Vn0_2[1]); //dy/dt
//    Jacobian_T[2][0] = x0_0[2]*dh1dr + x1_0[2]*dh2dr + x2_0[2]*dh3dr + 0.5*tPoint * (a0*dh1dr*Vn0_0[2] + a1*dh2dr*Vn0_1[2] + a2*dh3dr*Vn0_2[2]); //dz/dr
//    Jacobian_T[2][1] = x0_0[2]*dh1ds + x1_0[2]*dh2ds + x2_0[2]*dh3ds + 0.5*tPoint * (a0*dh1ds*Vn0_0[2] + a1*dh2ds*Vn0_1[2] + a2*dh3ds*Vn0_2[2]); //dz/ds
//    Jacobian_T[2][2] = 0.5 * (a0*h1*Vn0_0[2] + a1*h2*Vn0_1[2] + a2*h3*Vn0_2[2]); //dz/dt
//    matrixTranspose(Jacobian_T,Jacobian);
//    //invert Jacobian to bring r,s,t variables to x,y,z
//    double detJ = Jacobian[0][0]*(Jacobian[1][1]*Jacobian[2][2] - Jacobian[2][1]*Jacobian[1][2]) - Jacobian[0][1]*(Jacobian[1][0]*Jacobian[2][2] - Jacobian[1][2]*Jacobian[2][0]) + Jacobian[0][2]*(Jacobian[1][0]*Jacobian[2][1] - Jacobian[1][1]*Jacobian[2][0]);
//    double invdet = 1. / detJ;
//    vector<vector<double> > Jinv(3,vector<double>(3,0.));
//    Jinv[0][0] = (Jacobian[1][1] * Jacobian[2][2] - Jacobian[2][1] * Jacobian[1][2]) * invdet;
//    Jinv[0][1] = (Jacobian[0][2] * Jacobian[2][1] - Jacobian[0][1] * Jacobian[2][2]) * invdet;
//    Jinv[0][2] = (Jacobian[0][1] * Jacobian[1][2] - Jacobian[0][2] * Jacobian[1][1]) * invdet;
//    Jinv[1][0] = (Jacobian[1][2] * Jacobian[2][0] - Jacobian[1][0] * Jacobian[2][2]) * invdet;
//    Jinv[1][1] = (Jacobian[0][0] * Jacobian[2][2] - Jacobian[0][2] * Jacobian[2][0]) * invdet;
//    Jinv[1][2] = (Jacobian[1][0] * Jacobian[0][2] - Jacobian[0][0] * Jacobian[1][2]) * invdet;
//    Jinv[2][0] = (Jacobian[1][0] * Jacobian[2][1] - Jacobian[2][0] * Jacobian[1][1]) * invdet;
//    Jinv[2][1] = (Jacobian[2][0] * Jacobian[0][1] - Jacobian[0][0] * Jacobian[2][1]) * invdet;
//    Jinv[2][2] = (Jacobian[0][0] * Jacobian[1][1] - Jacobian[1][0] * Jacobian[0][1]) * invdet;
//    ////dumpMat(Jacobian,"Jacobain");
//    ////dumpMat(Jinv,"Jinv");
//    vector<vector<double> > dhmat0(3,vector<double>(3,0.));
//    vector<vector<double> > G0(3,vector<double>(3,0.));
//    vector<vector<double> > g1t(3,vector<double>(3,0.));
//    vector<vector<double> > g2t(3,vector<double>(3,0.));
//    hmat = {h1,h2,h3};
//    for (int i=0;i<3;++i) { 
//      for (int k=0;k<3;++k) { 
//        //0_hk,i
//        dhmat0[i][k] = Jinv[i][0]*dhmat[0][k] + Jinv[i][1]*dhmat[1][k];
//        //0_Gi^k
//        G0[i][k] = tPoint * (Jinv[i][0]*dhmat[0][k] + Jinv[i][1]*dhmat[1][k]) + Jinv[i][2]*hmat[k];
//      }
//    }
//    //t_g1i^k
//    g1t[0][0] = -0.5*a0*V2t_0[0]; 
//    g1t[0][1] = -0.5*a1*V2t_1[0]; 
//    g1t[0][2] = -0.5*a2*V2t_2[0]; 
//    g1t[1][0] = -0.5*a0*V2t_0[1]; 
//    g1t[1][1] = -0.5*a1*V2t_1[1]; 
//    g1t[1][2] = -0.5*a2*V2t_2[1]; 
//    g1t[2][0] = -0.5*a0*V2t_0[2]; 
//    g1t[2][1] = -0.5*a1*V2t_1[2]; 
//    g1t[2][2] = -0.5*a2*V2t_2[2]; 
//    //t_g2i^k   
//    g2t[0][0] =  0.5*a0*V1t_0[0]; 
//    g2t[0][1] =  0.5*a1*V1t_1[0]; 
//    g2t[0][2] =  0.5*a2*V1t_2[0]; 
//    g2t[1][0] =  0.5*a0*V1t_0[1]; 
//    g2t[1][1] =  0.5*a1*V1t_1[1]; 
//    g2t[1][2] =  0.5*a2*V1t_2[1]; 
//    g2t[2][0] =  0.5*a0*V1t_0[2]; 
//    g2t[2][1] =  0.5*a1*V1t_1[2]; 
//    g2t[2][2] =  0.5*a2*V1t_2[2]; 
//    vector<vector<double> > BNL(9,vector<double>(15,0.));
//    vector<vector<double> > BNL_T(15,vector<double>(9,0.));
//    vector<vector<double> > BNL_TS(15,vector<double>(9,0.));
//    vector<vector<double> > identPKST(9,vector<double>(9,0.));
//    for (int i=0;i<3;++i) { 
//      identPKST[i  ][i  ] = PKSTmat_[0][0];
//      identPKST[i  ][i+3] = PKSTmat_[0][1];
//      identPKST[i  ][i+6] = PKSTmat_[0][2];
//      identPKST[i+3][i  ] = PKSTmat_[1][0];
//      identPKST[i+3][i+3] = PKSTmat_[1][1];
//      identPKST[i+3][i+6] = PKSTmat_[1][2];
//      identPKST[i+6][i  ] = PKSTmat_[2][0];
//      identPKST[i+6][i+3] = PKSTmat_[2][1];
//      identPKST[i+6][i+6] = PKSTmat_[2][2];
//    }
//    //fill BNL
//    for (int i=0;i<3;++i) { 
//      for (int k=0;k<3;++k) { 
//        int kIndex = i + (ndof_loc-1)*k;
//        BNL[i  ][kIndex] = dhmat0[0][k];
//        BNL[i+3][kIndex] = dhmat0[1][k];
//        BNL[i+6][kIndex] = dhmat0[2][k];
//      }
//    }
//    for (int i=0;i<3;++i) { 
//      for (int k=0;k<3;++k) { 
//        int g1Index = 3 + (ndof_loc-1)*k;
//        int g2Index = 4 + (ndof_loc-1)*k;
//        BNL[i  ][g1Index] = g1t[i][k] * G0[0][k];
//        BNL[i  ][g2Index] = g2t[i][k] * G0[0][k];
//        BNL[i+3][g1Index] = g1t[i][k] * G0[1][k];
//        BNL[i+3][g2Index] = g2t[i][k] * G0[1][k];
//        BNL[i+6][g1Index] = g1t[i][k] * G0[2][k];
//        BNL[i+6][g2Index] = g2t[i][k] * G0[2][k];
//      }
//    }
//    matrixTranspose(BNL,BNL_T);
//    vector<vector<double> >knl1(18,vector<double>(18,0.));
//    vector<vector<double> >knl1_15B(15,vector<double>(15,0.));
//    mm_mult(BNL_T,identPKST,BNL_TS);
//    mm_mult(BNL_TS,BNL,knl1_15B);
//    //fill KNL1 with zero in drilling dof
//    for (int i=0;i<5;++i) { 
//      for (int j=0;j<5;++j) { 
//        knl1[i   ][j   ] = knl1_15B[i   ][j];
//        knl1[i+6 ][j   ] = knl1_15B[i+5 ][j];
//        knl1[i+12][j   ] = knl1_15B[i+10][j];
//        knl1[i   ][j+6 ] = knl1_15B[i   ][j+5];
//        knl1[i+6 ][j+6 ] = knl1_15B[i+5 ][j+5];
//        knl1[i+12][j+6 ] = knl1_15B[i+10][j+5];
//        knl1[i   ][j+12] = knl1_15B[i   ][j+10];
//        knl1[i+6 ][j+12] = knl1_15B[i+5 ][j+10];
//        knl1[i+12][j+12] = knl1_15B[i+10][j+10];
//      }
//    }
//    //dumpMat(knl1,"knl1 BATHE");
//    wait();

  }
  
  void linearTriElement::comp_KNL2_MITC3(const int& faceID,
                                         double& rPoint,
                                         double& sPoint,
                                         double& tPoint)
  {
//    //thickness 
//    double a0 = thickness_[faceID];
//    double a1 = thickness_[faceID];
//    double a2 = thickness_[faceID];
//    //create 2D shape function derivatives
//    //h1 = 1-r-s
//    double dh1dr = -1.;
//    double dh1ds = -1.;
//    double dh1dt =  0.;
//    //h2 = r
//    double dh2dr =  1.;
//    double dh2ds =  0.;
//    double dh2dt =  0.;
//    //h3 = s
//    double dh3dr =  0.;
//    double dh3ds =  1.;
//    double dh3dt =  0.;
//    //for tying points
//    double r1 = 0.5;
//    double s1 = 0.;
//    double r2 = 0.;
//    double s2 = 0.5;
//    double r3 = 0.5;
//    double s3 = 0.5;
//    //2D shape functions at tying points
//    double h1_TP1 = 1.-r1-s1;
//    double h2_TP1 = r1;
//    double h3_TP1 = s1;
//    double h1_TP2 = 1.-r2-s2;
//    double h2_TP2 = r2;
//    double h3_TP2 = s2;
//    double h1_TP3 = 1.-r3-s3;
//    double h2_TP3 = r3;
//    double h3_TP3 = s3;
//    //2D shape functions
//    double h1 = 1.-rPoint-sPoint;
//    double h2 = rPoint;
//    double h3 = sPoint;
//    vector<double> duq1dr(3,0.);
//    vector<double> duq2dr(3,0.);
//    vector<double> duq3dr(3,0.);
//    vector<double> duq1ds(3,0.);
//    vector<double> duq2ds(3,0.);
//    vector<double> duq3ds(3,0.);
//    vector<double> duq1dt(3,0.);
//    vector<double> duq2dt(3,0.);
//    vector<double> duq3dt(3,0.);
//    vector<double> duq1dt_TP1(3,0.);
//    vector<double> duq2dt_TP1(3,0.);
//    vector<double> duq3dt_TP1(3,0.);
//    vector<double> duq1dt_TP2(3,0.);
//    vector<double> duq2dt_TP2(3,0.);
//    vector<double> duq3dt_TP2(3,0.);
//    vector<double> duq1dt_TP3(3,0.);
//    vector<double> duq2dt_TP3(3,0.);
//    vector<double> duq3dt_TP3(3,0.);
//    for (int i=0;i<3;++i) 
//      { 
//        duq1dr[i] = -0.25*tPoint*a0*dh1dr*Vnt_0_[i];
//        duq2dr[i] = -0.25*tPoint*a1*dh2dr*Vnt_1_[i];
//        duq3dr[i] = -0.25*tPoint*a2*dh3dr*Vnt_2_[i];
//        duq1ds[i] = -0.25*tPoint*a0*dh1ds*Vnt_0_[i];
//        duq2ds[i] = -0.25*tPoint*a1*dh2ds*Vnt_1_[i];
//        duq3ds[i] = -0.25*tPoint*a2*dh3ds*Vnt_2_[i];
//        duq1dt[i] = -0.25*a0*h1*Vnt_0_[i];
//        duq2dt[i] = -0.25*a1*h2*Vnt_1_[i];
//        duq3dt[i] = -0.25*a2*h3*Vnt_2_[i];
//        duq1dt_TP1[i] = -0.25*a0*h1_TP1*Vnt_0_[i];
//        duq2dt_TP1[i] = -0.25*a1*h2_TP1*Vnt_1_[i];
//        duq3dt_TP1[i] = -0.25*a2*h3_TP1*Vnt_2_[i];
//        duq1dt_TP2[i] = -0.25*a0*h1_TP2*Vnt_0_[i];
//        duq2dt_TP2[i] = -0.25*a1*h2_TP2*Vnt_1_[i];
//        duq3dt_TP2[i] = -0.25*a2*h3_TP2*Vnt_2_[i];
//        duq1dt_TP3[i] = -0.25*a0*h1_TP3*Vnt_0_[i];
//        duq2dt_TP3[i] = -0.25*a1*h2_TP3*Vnt_1_[i];
//        duq3dt_TP3[i] = -0.25*a2*h3_TP3*Vnt_2_[i];

//      }
//    vector<vector<double> > BNL2(6,vector<double>(18,0.));
//    //alpha 1 and beta 1
//    vector<vector<double> > n11(15,vector<double>(15,0.));
//    vector<vector<double> > n12(15,vector<double>(15,0.));
//    vector<vector<double> > n13(15,vector<double>(15,0.));
//    vector<vector<double> > n13_TP1(15,vector<double>(15,0.));
//    vector<vector<double> > n13_TP2(15,vector<double>(15,0.));
//    vector<vector<double> > n13_TP3(15,vector<double>(15,0.));
//    vector<vector<double> > n21(15,vector<double>(15,0.));
//    vector<vector<double> > n22(15,vector<double>(15,0.));
//    vector<vector<double> > n23(15,vector<double>(15,0.));
//    vector<vector<double> > n23_TP1(15,vector<double>(15,0.));
//    vector<vector<double> > n23_TP2(15,vector<double>(15,0.));
//    vector<vector<double> > n23_TP3(15,vector<double>(15,0.));
//    vector<vector<double> > n31(15,vector<double>(15,0.));
//    vector<vector<double> > n31_TP1(15,vector<double>(15,0.));
//    vector<vector<double> > n31_TP2(15,vector<double>(15,0.));
//    vector<vector<double> > n31_TP3(15,vector<double>(15,0.));
//    vector<vector<double> > n32(15,vector<double>(15,0.));
//    vector<vector<double> > n32_TP1(15,vector<double>(15,0.));
//    vector<vector<double> > n32_TP2(15,vector<double>(15,0.));
//    vector<vector<double> > n32_TP3(15,vector<double>(15,0.));
//    vector<vector<double> > n33(15,vector<double>(15,0.));
//    vector<vector<double> > cTying(15,vector<double>(15,0.));
//    n11[ 3][ 3]       = dprod(duq1dr,gt_1_);
//    n11[ 4][ 4]       = dprod(duq1dr,gt_1_);
//    n11[ 8][ 8]       = dprod(duq2dr,gt_1_);
//    n11[ 9][ 9]       = dprod(duq2dr,gt_1_);
//    n11[13][13]       = dprod(duq3dr,gt_1_);
//    n11[14][14]       = dprod(duq3dr,gt_1_);
//    n22[ 3][ 3]       = dprod(duq1ds,gt_2_);
//    n22[ 4][ 4]       = dprod(duq1ds,gt_2_);
//    n22[ 8][ 8]       = dprod(duq2ds,gt_2_);
//    n22[ 9][ 9]       = dprod(duq2ds,gt_2_);
//    n22[13][13]       = dprod(duq3ds,gt_2_);
//    n22[14][14]       = dprod(duq3ds,gt_2_);
//    n33[ 3][ 3]       = dprod(duq1dt,gt_3_);
//    n33[ 4][ 4]       = dprod(duq1dt,gt_3_);
//    n33[ 8][ 8]       = dprod(duq2dt,gt_3_);
//    n33[ 9][ 9]       = dprod(duq2dt,gt_3_);
//    n33[13][13]       = dprod(duq3dt,gt_3_);
//    n33[14][14]       = dprod(duq3dt,gt_3_);
//    n12[ 3][ 3]       = dprod(duq1dr,gt_2_) + dprod(duq1ds,gt_1_);
//    n12[ 4][ 4]       = dprod(duq1dr,gt_2_) + dprod(duq1ds,gt_1_);
//    n12[ 8][ 8]       = dprod(duq2dr,gt_2_) + dprod(duq2ds,gt_1_);
//    n12[ 9][ 9]       = dprod(duq2dr,gt_2_) + dprod(duq2ds,gt_1_);
//    n12[13][13]       = dprod(duq3dr,gt_2_) + dprod(duq3ds,gt_1_);
//    n12[14][14]       = dprod(duq3dr,gt_2_) + dprod(duq3ds,gt_1_);
//    n13_TP1[ 3][ 3]   = dprod(duq1dr,gt_3_TP1_) + dprod(duq1dt_TP1,gt_1_);
//    n13_TP1[ 4][ 4]   = dprod(duq1dr,gt_3_TP1_) + dprod(duq1dt_TP1,gt_1_);
//    n13_TP1[ 8][ 8]   = dprod(duq2dr,gt_3_TP1_) + dprod(duq2dt_TP1,gt_1_);
//    n13_TP1[ 9][ 9]   = dprod(duq2dr,gt_3_TP1_) + dprod(duq2dt_TP1,gt_1_);
//    n13_TP1[13][13]   = dprod(duq3dr,gt_3_TP1_) + dprod(duq3dt_TP1,gt_1_);
//    n13_TP1[14][14]   = dprod(duq3dr,gt_3_TP1_) + dprod(duq3dt_TP1,gt_1_);
//    n13_TP2[ 3][ 3]   = dprod(duq1dr,gt_3_TP2_) + dprod(duq1dt_TP2,gt_1_);
//    n13_TP2[ 4][ 4]   = dprod(duq1dr,gt_3_TP2_) + dprod(duq1dt_TP2,gt_1_);
//    n13_TP2[ 8][ 8]   = dprod(duq2dr,gt_3_TP2_) + dprod(duq2dt_TP2,gt_1_);
//    n13_TP2[ 9][ 9]   = dprod(duq2dr,gt_3_TP2_) + dprod(duq2dt_TP2,gt_1_);
//    n13_TP2[13][13]   = dprod(duq3dr,gt_3_TP2_) + dprod(duq3dt_TP2,gt_1_);
//    n13_TP2[14][14]   = dprod(duq3dr,gt_3_TP2_) + dprod(duq3dt_TP2,gt_1_);
//    n13_TP3[ 3][ 3]   = dprod(duq1dr,gt_3_TP3_) + dprod(duq1dt_TP3,gt_1_);
//    n13_TP3[ 4][ 4]   = dprod(duq1dr,gt_3_TP3_) + dprod(duq1dt_TP3,gt_1_);
//    n13_TP3[ 8][ 8]   = dprod(duq2dr,gt_3_TP3_) + dprod(duq2dt_TP3,gt_1_);
//    n13_TP3[ 9][ 9]   = dprod(duq2dr,gt_3_TP3_) + dprod(duq2dt_TP3,gt_1_);
//    n13_TP3[13][13]   = dprod(duq3dr,gt_3_TP3_) + dprod(duq3dt_TP3,gt_1_);
//    n13_TP3[14][14]   = dprod(duq3dr,gt_3_TP3_) + dprod(duq3dt_TP3,gt_1_);
//    n23_TP1[ 3][ 3]   = dprod(duq1ds,gt_3_TP1_) + dprod(duq1dt_TP1,gt_2_);
//    n23_TP1[ 4][ 4]   = dprod(duq1ds,gt_3_TP1_) + dprod(duq1dt_TP1,gt_2_);
//    n23_TP1[ 8][ 8]   = dprod(duq2ds,gt_3_TP1_) + dprod(duq2dt_TP1,gt_2_);
//    n23_TP1[ 9][ 9]   = dprod(duq2ds,gt_3_TP1_) + dprod(duq2dt_TP1,gt_2_);
//    n23_TP1[13][13]   = dprod(duq3ds,gt_3_TP1_) + dprod(duq3dt_TP1,gt_2_);
//    n23_TP1[14][14]   = dprod(duq3ds,gt_3_TP1_) + dprod(duq3dt_TP1,gt_2_);
//    n23_TP2[ 3][ 3]   = dprod(duq1ds,gt_3_TP2_) + dprod(duq1dt_TP2,gt_2_);
//    n23_TP2[ 4][ 4]   = dprod(duq1ds,gt_3_TP2_) + dprod(duq1dt_TP2,gt_2_);
//    n23_TP2[ 8][ 8]   = dprod(duq2ds,gt_3_TP2_) + dprod(duq2dt_TP2,gt_2_);
//    n23_TP2[ 9][ 9]   = dprod(duq2ds,gt_3_TP2_) + dprod(duq2dt_TP2,gt_2_);
//    n23_TP2[13][13]   = dprod(duq3ds,gt_3_TP2_) + dprod(duq3dt_TP2,gt_2_);
//    n23_TP2[14][14]   = dprod(duq3ds,gt_3_TP2_) + dprod(duq3dt_TP2,gt_2_);
//    n23_TP3[ 3][ 3]   = dprod(duq1ds,gt_3_TP3_) + dprod(duq1dt_TP3,gt_2_);
//    n23_TP3[ 4][ 4]   = dprod(duq1ds,gt_3_TP3_) + dprod(duq1dt_TP3,gt_2_);
//    n23_TP3[ 8][ 8]   = dprod(duq2ds,gt_3_TP3_) + dprod(duq2dt_TP3,gt_2_);
//    n23_TP3[ 9][ 9]   = dprod(duq2ds,gt_3_TP3_) + dprod(duq2dt_TP3,gt_2_);
//    n23_TP3[13][13]   = dprod(duq3ds,gt_3_TP3_) + dprod(duq3dt_TP3,gt_2_);
//    n23_TP3[14][14]   = dprod(duq3ds,gt_3_TP3_) + dprod(duq3dt_TP3,gt_2_);        
//    //with corrections
//    for (int i=0;i<15;++i)
//      {
//        for (int j=0;j<15;++j)
//          {            
//            cTying[i][j] = n23_TP2[i][j] - n13_TP1[i][j] - n23_TP3[i][j] + n13_TP3[i][j];
//            n23[i][j] = n23_TP2[i][j] - cTying[i][j]*rPoint;
//            n13[i][j] = n13_TP1[i][j] + cTying[i][j]*sPoint;
//          }
//      }
//    //without corrections
//    //n13[ 3][ 3] = (dprod(duq1dr,gt_3) + dprod(duq1dt,gt_1));
//    //n13[ 4][ 4] = (dprod(duq1dr,gt_3) + dprod(duq1dt,gt_1));
//    //n13[ 8][ 8] = (dprod(duq2dr,gt_3) + dprod(duq2dt,gt_1));
//    //n13[ 9][ 9] = (dprod(duq2dr,gt_3) + dprod(duq2dt,gt_1));
//    //n13[13][13] = (dprod(duq3dr,gt_3) + dprod(duq3dt,gt_1));
//    //n13[14][14] = (dprod(duq3dr,gt_3) + dprod(duq3dt,gt_1));
//    //n23[ 3][ 3] = (dprod(duq1ds,gt_3) + dprod(duq1dt,gt_2));
//    //n23[ 4][ 4] = (dprod(duq1ds,gt_3) + dprod(duq1dt,gt_2));
//    //n23[ 8][ 8] = (dprod(duq2ds,gt_3) + dprod(duq2dt,gt_2));
//    //n23[ 9][ 9] = (dprod(duq2ds,gt_3) + dprod(duq2dt,gt_2));
//    //n23[13][13] = (dprod(duq3ds,gt_3) + dprod(duq3dt,gt_2));
//    //n23[14][14] = (dprod(duq3ds,gt_3) + dprod(duq3dt,gt_2));        
//    //    //dumpMat(n11,"n11");
//    //    //dumpMat(n12,"n12");
//    //    //dumpMat(n13,"n13");
//    //    //dumpMat(n22,"n22");
//    //    //dumpMat(n23,"n23");
//    //    //dumpMat(n33,"n33");

//    vector<vector<vector<vector<double> > > > nrs(3,vector<vector<vector<double> > >(3,vector<vector<double> >(15,vector<double>(15,0.))));
//    for (int i=0;i<15;++i) 
//      { 
//        for (int j=0;j<15;++j) 
//          { 
//            nrs[0][0][i][j] = n11[i][j]; //rr
//            nrs[0][1][i][j] = n12[i][j]; //rs
//            nrs[0][2][i][j] = n13[i][j]; //rt
//            nrs[1][0][i][j] = n12[i][j]; //sr
//            nrs[1][1][i][j] = n22[i][j]; //ss
//            nrs[1][2][i][j] = n23[i][j]; //st
//            nrs[2][0][i][j] = n13[i][j]; //tr
//            nrs[2][1][i][j] = n23[i][j]; //ts
//            nrs[2][2][i][j] = n33[i][j]; //tt
//          }
//      }
//    vector<vector<vector<vector<double> > > > nxy(3,vector<vector<vector<double> > >(3,vector<vector<double> >(15,vector<double>(15,0.))));
//    for (int i=0;i<15;++i) {
//      for (int j=0;j<15;++j) {                
//        for (int k=0;k<3;++k) {
//          for (int l=0;l<3;++l) {
//            for (int m=0;m<3;++m) {
//              for (int n=0;n<3;++n) {                
//                nxy[k][l][i][j] += nrs[m][n][i][j]*dprod(ematT_[k],g0matT_[m])*dprod(ematT_[l],g0matT_[n]);
//              }
//            }            
//          }
//        }
//      }
//    }
//    //gamma = 2*epsilon
//    for (int i=0;i<15;++i) 
//      {
//        for (int j=0;j<15;++j) 
//          {
//            nxy[0][1][i][j] = 1.*nxy[0][1][i][j];
//            nxy[1][0][i][j] = 1.*nxy[1][0][i][j];
//            nxy[0][2][i][j] = 1.*nxy[0][2][i][j];
//            nxy[2][0][i][j] = 1.*nxy[2][0][i][j];
//            nxy[1][2][i][j] = 1.*nxy[1][2][i][j];
//            nxy[2][1][i][j] = 1.*nxy[2][1][i][j];
//          }
//      }
//    vector<vector<double> > knl2_15(15,vector<double>(15,0.));  
//    for (int k=0;k<3;++k) {
//      for (int l=0;l<3;++l) {
//        for (int i=0;i<15;++i) {
//          for (int j=0;j<15;++j) {            
//            knl2_15[i][j] += nxy[k][l][i][j]*PKSTmat_[k][l];
//          }
//        }
//      } 
//    }
//    //fill KNL1 with zero in drilling dof
//    for (int i=0;i<5;++i) 
//      { 
//        for (int j=0;j<5;++j) 
//          { 
//            KNL2mat_[i   ][j   ] = knl2_15[i   ][j];
//            KNL2mat_[i+6 ][j   ] = knl2_15[i+5 ][j];
//            KNL2mat_[i+12][j   ] = knl2_15[i+10][j];
//            KNL2mat_[i   ][j+6 ] = knl2_15[i   ][j+5];
//            KNL2mat_[i+6 ][j+6 ] = knl2_15[i+5 ][j+5];
//            KNL2mat_[i+12][j+6 ] = knl2_15[i+10][j+5];
//            KNL2mat_[i   ][j+12] = knl2_15[i   ][j+10];
//            KNL2mat_[i+6 ][j+12] = knl2_15[i+5 ][j+10];
//            KNL2mat_[i+12][j+12] = knl2_15[i+10][j+10];
//          }
//      }
////    //dumpMat(KNL2mat,"knl2mat");
////    wait();
  }
  
  void linearTriElement::assignExactValuesPerFace(vector<double>& u,
                                                  int& fn0,
                                                  int& fn1,
                                                  int& fn2)
  {
//    vector<int> nodes_loc(3,0);
//    nodes_loc[0] = fn0;
//    nodes_loc[1] = fn1;
//    nodes_loc[2] = fn2;
//    //cout << endl;
//    for (int i=0;i<3;++i)
//      {
//        int n = nodes_loc[i];
//        double xcoord = nodes0_[n*xyz+0];
//        double ycoord = nodes0_[n*xyz+1];
//        double zcoord = nodes0_[n*xyz+2];
//        //translations
//        u[i*6+0] = 0.;//2. *xcoord +    ycoord;
//        u[i*6+1] = 0.;//    xcoord + 5.*ycoord;
//        u[i*6+2] = 0.;//0.5*xcoord + 3.*ycoord;
//        //rotations
//        u[i*6+3] = 2. *xcoord +    ycoord + 2.;
//        u[i*6+4] = 3. *ycoord +    zcoord + 1.;
//        u[i*6+5] = 0.;
//      }     
  }
    
  void linearTriElement::comp_NL_Me_MITC3(vector<double>& Me_loc,
                                          const int&      faceID)
  {
  }
  /*
  //computing local mass matrix
  // get face Id
  vector<int> face(nnodes_loc);
  for (int n=0;n<nnodes_loc;++n) face_nodes[n]= (*faces_)[faceID][n]; 
  // face nodes
  int fn0=face_nodes[0];
  int fn1=face_nodes[1];
  int fn2=face_nodes[2];
  vector<double> x0(DIM,0.);
  vector<double> x1(DIM,0.);
  vector<double> x2(DIM,0.);

  for (int idir=0;idir<DIM;++idir) {
  x0[idir] = nodes_[fn0][idir];
  x1[idir] = nodes_[fn1][idir];      
  x2[idir] = nodes_[fn2][idir]; 
  }

  // obtain transformation matrix from local nodes
  vector<double> myVec3(3,0.5);
  vector<vector<double> > transMat(3,dummyVec3);
  transMat=compTransMat(x0,x1,x2);  
  vector<vector<double> > LVec(3,dummyVec3);
  LVec[0][2]=0.;
  LVec[1][0]=0.;
  LVec[2][1]=0.;
  vector<double> dummyVec2(15,0.);
  vector<double> dummyVec4(18,0.);
  vector<vector<double> > Nmatloc( 3,dummyVec2);
  vector<vector<double> > Nmat   (15,dummyVec2);
  vector<vector<double> > Mat    (18,dummyVec4);

  //    vector<double> weights(3,1./3.);
  vector<double> weights(1,1.);
  int noWeights = weights.size();
  if (noWeights == 1) { 
  LVec[0][0] = 1./3;
  LVec[0][1] = 1./3;
  LVec[0][2] = 1./3;
  }

  double area=0.;
  for (int k=0;k<noWeights;++k) {
  area=0.;
  compShapeFunction(Nmatloc,area,faceID,LVec[k]);        

  mTm_mult(Nmatloc,Nmatloc,Nmat); 

  double factor=area*thickness_[faceID]*density_[faceID];
  for (int r=0;r<3;++r) {
  for (int c=0;c<3;++c) {
  for (int rr=0;rr<5;++rr) {
  for (int cc=0;cc<5;++cc) {
  Mat[rr+r*6][cc+c*6]+=weights[k]*Nmat[rr+r*5][cc+c*5]*factor;
  }
  }
  }
  }
  }

  vector<double> dumVec(18,0.);      
  vector<vector<double> > TT(18,dumVec);
  for (int n=0;n<3;++n) {
  for (int i=0;i<3;++i) {
  for (int j=0;j<3;++j) {
  int iloc=i+n*ndof_loc;
  int jloc=j+n*ndof_loc;
  TT[iloc][jloc]=transMat[i][j];
  }
  }
  for (int i=0;i<3;++i) {
  for (int j=0;j<3;++j) {
  int iloc=i+n*ndof_loc+3;
  int jloc=j+n*ndof_loc+3;
  TT[iloc][jloc]=transMat[i][j]; 
  }
  }
  }

  aTma_mult(TT,Mat,Me_loc);
  }
  */
}

//    /////////////////////Jedit 11/13/17 - adding KNL1 from Bathe paper//////////////////////////////////////////////////////
//    //create Jacobian
//    vector<vector<double> > Jacobian(3,vector<double>(3,0.));
//    vector<vector<double> > Jacobian_T(3,vector<double>(3,0.));
//    Jacobian_T[0][0] = x0_0[0]*dh1dr + x1_0[0]*dh2dr + x2_0[0]*dh3dr + 0.5*tPoint * (a0*dh1dr*Vn0_0[0] + a1*dh2dr*Vn0_1[0] + a2*dh3dr*Vn0_2[0]); //dx/dr
//    Jacobian_T[0][1] = x0_0[0]*dh1ds + x1_0[0]*dh2ds + x2_0[0]*dh3ds + 0.5*tPoint * (a0*dh1ds*Vn0_0[0] + a1*dh2ds*Vn0_1[0] + a2*dh3ds*Vn0_2[0]); //dx/ds
//    Jacobian_T[0][2] = 0.5 * (a0*h1*Vn0_0[0] + a1*h2*Vn0_1[0] + a2*h3*Vn0_2[0]); //dx/dt
//    Jacobian_T[1][0] = x0_0[1]*dh1dr + x1_0[1]*dh2dr + x2_0[1]*dh3dr + 0.5*tPoint * (a0*dh1dr*Vn0_0[1] + a1*dh2dr*Vn0_1[1] + a2*dh3dr*Vn0_2[1]); //dy/dr
//    Jacobian_T[1][1] = x0_0[1]*dh1ds + x1_0[1]*dh2ds + x2_0[1]*dh3ds + 0.5*tPoint * (a0*dh1ds*Vn0_0[1] + a1*dh2ds*Vn0_1[1] + a2*dh3ds*Vn0_2[1]); //dy/ds
//    Jacobian_T[1][2] = 0.5 * (a0*h1*Vn0_0[1] + a1*h2*Vn0_1[1] + a2*h3*Vn0_2[1]); //dy/dt
//    Jacobian_T[2][0] = x0_0[2]*dh1dr + x1_0[2]*dh2dr + x2_0[2]*dh3dr + 0.5*tPoint * (a0*dh1dr*Vn0_0[2] + a1*dh2dr*Vn0_1[2] + a2*dh3dr*Vn0_2[2]); //dz/dr
//    Jacobian_T[2][1] = x0_0[2]*dh1ds + x1_0[2]*dh2ds + x2_0[2]*dh3ds + 0.5*tPoint * (a0*dh1ds*Vn0_0[2] + a1*dh2ds*Vn0_1[2] + a2*dh3ds*Vn0_2[2]); //dz/ds
//    Jacobian_T[2][2] = 0.5 * (a0*h1*Vn0_0[2] + a1*h2*Vn0_1[2] + a2*h3*Vn0_2[2]); //dz/dt
//    matrixTranspose(Jacobian_T,Jacobian);
//    //invert Jacobian to bring r,s,t variables to x,y,z
//    double detJ = Jacobian[0][0]*(Jacobian[1][1]*Jacobian[2][2] - Jacobian[2][1]*Jacobian[1][2]) - Jacobian[0][1]*(Jacobian[1][0]*Jacobian[2][2] - Jacobian[1][2]*Jacobian[2][0]) + Jacobian[0][2]*(Jacobian[1][0]*Jacobian[2][1] - Jacobian[1][1]*Jacobian[2][0]);
//    double invdet = 1. / detJ;
//    vector<vector<double> > Jinv(3,vector<double>(3,0.));
//    Jinv[0][0] = (Jacobian[1][1] * Jacobian[2][2] - Jacobian[2][1] * Jacobian[1][2]) * invdet;
//    Jinv[0][1] = (Jacobian[0][2] * Jacobian[2][1] - Jacobian[0][1] * Jacobian[2][2]) * invdet;
//    Jinv[0][2] = (Jacobian[0][1] * Jacobian[1][2] - Jacobian[0][2] * Jacobian[1][1]) * invdet;
//    Jinv[1][0] = (Jacobian[1][2] * Jacobian[2][0] - Jacobian[1][0] * Jacobian[2][2]) * invdet;
//    Jinv[1][1] = (Jacobian[0][0] * Jacobian[2][2] - Jacobian[0][2] * Jacobian[2][0]) * invdet;
//    Jinv[1][2] = (Jacobian[1][0] * Jacobian[0][2] - Jacobian[0][0] * Jacobian[1][2]) * invdet;
//    Jinv[2][0] = (Jacobian[1][0] * Jacobian[2][1] - Jacobian[2][0] * Jacobian[1][1]) * invdet;
//    Jinv[2][1] = (Jacobian[2][0] * Jacobian[0][1] - Jacobian[0][0] * Jacobian[2][1]) * invdet;
//    Jinv[2][2] = (Jacobian[0][0] * Jacobian[1][1] - Jacobian[1][0] * Jacobian[0][1]) * invdet;
//    ////dumpMat(Jacobian,"Jacobain");
//    //dumpMat(Jinv,"Jinv");
//    vector<vector<double> > dhmat0(3,vector<double>(3,0.));
//    vector<double> hmat(3,0.);
//    vector<vector<double> > G0(3,vector<double>(3,0.));
//    vector<vector<double> > g1t(3,vector<double>(3,0.));
//    vector<vector<double> > g2t(3,vector<double>(3,0.));
//    hmat = {h1,h2,h3};
//    for (int i=0;i<3;++i) { 
//      for (int k=0;k<3;++k) { 
//        //0_hk,i
//        dhmat0[i][k] = Jinv[i][0]*dhmat[0][k] + Jinv[i][1]*dhmat[1][k];
//        //0_Gi^k
//        G0[i][k] = tPoint * (Jinv[i][0]*dhmat[0][k] + Jinv[i][1]*dhmat[1][k]) + Jinv[i][2]*hmat[k];
//      }
//    }
//    //t_g1i^k
//    g1t[0][0] = -0.5*a0*V2t_0[0]; 
//    g1t[0][1] = -0.5*a1*V2t_1[0]; 
//    g1t[0][2] = -0.5*a2*V2t_2[0]; 
//    g1t[1][0] = -0.5*a0*V2t_0[1]; 
//    g1t[1][1] = -0.5*a1*V2t_1[1]; 
//    g1t[1][2] = -0.5*a2*V2t_2[1]; 
//    g1t[2][0] = -0.5*a0*V2t_0[2]; 
//    g1t[2][1] = -0.5*a1*V2t_1[2]; 
//    g1t[2][2] = -0.5*a2*V2t_2[2]; 
//    //t_g2i^k   
//    g2t[0][0] =  0.5*a0*V1t_0[0]; 
//    g2t[0][1] =  0.5*a1*V1t_1[0]; 
//    g2t[0][2] =  0.5*a2*V1t_2[0]; 
//    g2t[1][0] =  0.5*a0*V1t_0[1]; 
//    g2t[1][1] =  0.5*a1*V1t_1[1]; 
//    g2t[1][2] =  0.5*a2*V1t_2[1]; 
//    g2t[2][0] =  0.5*a0*V1t_0[2]; 
//    g2t[2][1] =  0.5*a1*V1t_1[2]; 
//    g2t[2][2] =  0.5*a2*V1t_2[2]; 
//    vector<vector<double> > BNL(9,vector<double>(15,0.));
//    vector<vector<double> > BNL_T(15,vector<double>(9,0.));
//    vector<vector<double> > BNL_TS(15,vector<double>(9,0.));
//    vector<vector<double> > knl1_15(15,vector<double>(15,0.));
//    vector<vector<double> > identPKST(9,vector<double>(9,0.));
//    for (int i=0;i<3;++i) { 
//      identPKST[i  ][i  ] = PKSTmat_[0][0];
//      identPKST[i  ][i+3] = PKSTmat_[0][1];
//      identPKST[i  ][i+6] = PKSTmat_[0][2];
//      identPKST[i+3][i  ] = PKSTmat_[1][0];
//      identPKST[i+3][i+3] = PKSTmat_[1][1];
//      identPKST[i+3][i+6] = PKSTmat_[1][2];
//      identPKST[i+6][i  ] = PKSTmat_[2][0];
//      identPKST[i+6][i+3] = PKSTmat_[2][1];
//      identPKST[i+6][i+6] = PKSTmat_[2][2];
//    }
//    //fill BNL
//    for (int i=0;i<3;++i) { 
//      for (int k=0;k<3;++k) { 
//        int kIndex = i + (ndof_loc-1)*k;
//        BNL[i  ][kIndex] = dhmat0[0][k];
//        BNL[i+3][kIndex] = dhmat0[1][k];
//        BNL[i+6][kIndex] = dhmat0[2][k];
//      }
//    }
//    for (int i=0;i<3;++i) { 
//      for (int k=0;k<3;++k) { 
//        int g1Index = 3 + (ndof_loc-1)*k;
//        int g2Index = 4 + (ndof_loc-1)*k;
//        BNL[i  ][g1Index] = g1t[i][k] * G0[0][k];
//        BNL[i  ][g2Index] = g2t[i][k] * G0[0][k];
//        BNL[i+3][g1Index] = g1t[i][k] * G0[1][k];
//        BNL[i+3][g2Index] = g2t[i][k] * G0[1][k];
//        BNL[i+6][g1Index] = g1t[i][k] * G0[2][k];
//        BNL[i+6][g2Index] = g2t[i][k] * G0[2][k];
//      }
//    }
//    //dumpMat(BNL,"BNL");
//    vector<vector<double> > newIdent(9,vector<double>(9,0.));
//    for (int i=0;i<3;++i) 
//      { 
//        newIdent[i  ][i  ] = 1.;//1.;
//        newIdent[i  ][i+3] = 0.;//1.;
//        newIdent[i  ][i+6] = 0.;//1.;
//        newIdent[i+3][i  ] = 0.;//1.;
//        newIdent[i+3][i+3] = 1.;//1.;
//        newIdent[i+3][i+6] = 0.;//1.;
//        newIdent[i+6][i  ] = 0.;//1.;
//        newIdent[i+6][i+3] = 0.;//1.;
//        newIdent[i+6][i+6] = 1.;//1.;
//      }
//    matrixTranspose(BNL,BNL_T);
//    //dumpMat(newIdent,"newIdent");
//    vector<vector<double> >BTI(15,vector<double>(9,0.));
//    vector<vector<double> >BTB(15,vector<double>(15,0.));
//    mm_mult(BNL_T,newIdent,BTI);
//    mm_mult(BTI,BNL,BTB);
//    //dumpMat(BTB,"BTB NEW");
//    //dumpMat(BNL,"BNL");

//    mm_mult(BNL_T,identPKST,BNL_TS);
//    mm_mult(BNL_TS,BNL,knl1_15);
//    //fill KNL1 with zero in drilling dof
//    for (int i=0;i<5;++i) { 
//      for (int j=0;j<5;++j) { 
//        KNL1mat[i   ][j   ] = knl1_15[i   ][j];
//        KNL1mat[i+6 ][j   ] = knl1_15[i+5 ][j];
//        KNL1mat[i+12][j   ] = knl1_15[i+10][j];
//        KNL1mat[i   ][j+6 ] = knl1_15[i   ][j+5];
//        KNL1mat[i+6 ][j+6 ] = knl1_15[i+5 ][j+5];
//        KNL1mat[i+12][j+6 ] = knl1_15[i+10][j+5];
//        KNL1mat[i   ][j+12] = knl1_15[i   ][j+10];
//        KNL1mat[i+6 ][j+12] = knl1_15[i+5 ][j+10];
//        KNL1mat[i+12][j+12] = knl1_15[i+10][j+10];
//      }
//    }

    //dumpMat(BNL,"BNL");
    //dumpMat(knl1_15,"knl1_15");
    //dumpMat(KNL1mat,"KNL1mat new");
    //wait();

    //ADDING BL FROM BATHE - JEDIT 12/21/17

//    vector<vector<double> > BL_0(6,vector<double>(15,0.)); 
//    //create Jacobian
//    vector<vector<double> > Jacobian(3,vector<double>(3,0.));
//    Jacobian[0][0] = x0_0[0]*dh1dr + x1_0[0]*dh2dr + x2_0[0]*dh3dr + 0.5*tPoint * (a0*dh1dr*Vn0_0[0] + a1*dh2dr*Vn0_1[0] + a2*dh3dr*Vn0_2[0]); //dx/dr
//    Jacobian[1][0] = x0_0[0]*dh1ds + x1_0[0]*dh2ds + x2_0[0]*dh3ds + 0.5*tPoint * (a0*dh1ds*Vn0_0[0] + a1*dh2ds*Vn0_1[0] + a2*dh3ds*Vn0_2[0]); //dx/ds
//    Jacobian[2][0] = 0.5 * (a0*h1*Vn0_0[0] + a1*h2*Vn0_1[0] + a2*h3*Vn0_2[0]); //dx/dt
//    Jacobian[0][1] = x0_0[1]*dh1dr + x1_0[1]*dh2dr + x2_0[1]*dh3dr + 0.5*tPoint * (a0*dh1dr*Vn0_0[1] + a1*dh2dr*Vn0_1[1] + a2*dh3dr*Vn0_2[1]); //dy/dr
//    Jacobian[1][1] = x0_0[1]*dh1ds + x1_0[1]*dh2ds + x2_0[1]*dh3ds + 0.5*tPoint * (a0*dh1ds*Vn0_0[1] + a1*dh2ds*Vn0_1[1] + a2*dh3ds*Vn0_2[1]); //dy/ds
//    Jacobian[2][1] = 0.5 * (a0*h1*Vn0_0[1] + a1*h2*Vn0_1[1] + a2*h3*Vn0_2[1]); //dy/dt
//    Jacobian[0][2] = x0_0[2]*dh1dr + x1_0[2]*dh2dr + x2_0[2]*dh3dr + 0.5*tPoint * (a0*dh1dr*Vn0_0[2] + a1*dh2dr*Vn0_1[2] + a2*dh3dr*Vn0_2[2]); //dz/dr
//    Jacobian[1][2] = x0_0[2]*dh1ds + x1_0[2]*dh2ds + x2_0[2]*dh3ds + 0.5*tPoint * (a0*dh1ds*Vn0_0[2] + a1*dh2ds*Vn0_1[2] + a2*dh3ds*Vn0_2[2]); //dz/ds
//    Jacobian[2][2] = 0.5 * (a0*h1*Vn0_0[2] + a1*h2*Vn0_1[2] + a2*h3*Vn0_2[2]); //dz/dt
//    //invert Jacobian to bring r,s,t variables to x,y,z
//    double detJ = Jacobian[0][0]*(Jacobian[1][1]*Jacobian[2][2] - Jacobian[2][1]*Jacobian[1][2]) - Jacobian[0][1]*(Jacobian[1][0]*Jacobian[2][2] - Jacobian[1][2]*Jacobian[2][0]) + Jacobian[0][2]*(Jacobian[1][0]*Jacobian[2][1] - Jacobian[1][1]*Jacobian[2][0]);
//    double invdet = 1. / detJ;
//    vector<vector<double> > Jinv(3,vector<double>(3,0.));
//    Jinv[0][0] = (Jacobian[1][1] * Jacobian[2][2] - Jacobian[2][1] * Jacobian[1][2]) * invdet;
//    Jinv[0][1] = (Jacobian[0][2] * Jacobian[2][1] - Jacobian[0][1] * Jacobian[2][2]) * invdet;
//    Jinv[0][2] = (Jacobian[0][1] * Jacobian[1][2] - Jacobian[0][2] * Jacobian[1][1]) * invdet;
//    Jinv[1][0] = (Jacobian[1][2] * Jacobian[2][0] - Jacobian[1][0] * Jacobian[2][2]) * invdet;
//    Jinv[1][1] = (Jacobian[0][0] * Jacobian[2][2] - Jacobian[0][2] * Jacobian[2][0]) * invdet;
//    Jinv[1][2] = (Jacobian[1][0] * Jacobian[0][2] - Jacobian[0][0] * Jacobian[1][2]) * invdet;
//    Jinv[2][0] = (Jacobian[1][0] * Jacobian[2][1] - Jacobian[2][0] * Jacobian[1][1]) * invdet;
//    Jinv[2][1] = (Jacobian[2][0] * Jacobian[0][1] - Jacobian[0][0] * Jacobian[2][1]) * invdet;
//    Jinv[2][2] = (Jacobian[0][0] * Jacobian[1][1] - Jacobian[1][0] * Jacobian[0][1]) * invdet;
//    vector<vector<double> > dhmat0(3,vector<double>(3,0.));
//    vector<double> hmat(3,0.);
//    vector<vector<double> > G0(3,vector<double>(3,0.));
//    vector<vector<double> > g1t(3,vector<double>(3,0.));
//    vector<vector<double> > g2t(3,vector<double>(3,0.));
//    hmat = {h1,h2,h3};
//    for (int i=0;i<3;++i) { 
//      for (int k=0;k<3;++k) { 
//        //0_hk,i
//        dhmat0[i][k] = Jinv[i][0]*dhmat[0][k] + Jinv[i][1]*dhmat[1][k];
//        //0_Gi^k
//        G0[i][k] = tPoint * (Jinv[i][0]*dhmat[0][k] + Jinv[i][1]*dhmat[1][k]) + Jinv[i][2]*hmat[k];
//      }
//    }
//    //t_g1i^k
//    g1t[0][0] = -0.5*a0*V2t_0[0]; 
//    g1t[0][1] = -0.5*a1*V2t_1[0]; 
//    g1t[0][2] = -0.5*a2*V2t_2[0]; 
//    g1t[1][0] = -0.5*a0*V2t_0[1]; 
//    g1t[1][1] = -0.5*a1*V2t_1[1]; 
//    g1t[1][2] = -0.5*a2*V2t_2[1]; 
//    g1t[2][0] = -0.5*a0*V2t_0[2]; 
//    g1t[2][1] = -0.5*a1*V2t_1[2]; 
//    g1t[2][2] = -0.5*a2*V2t_2[2]; 
//    //t_g2i^k   
//    g2t[0][0] =  0.5*a0*V1t_0[0]; 
//    g2t[0][1] =  0.5*a1*V1t_1[0]; 
//    g2t[0][2] =  0.5*a2*V1t_2[0]; 
//    g2t[1][0] =  0.5*a0*V1t_0[1]; 
//    g2t[1][1] =  0.5*a1*V1t_1[1]; 
//    g2t[1][2] =  0.5*a2*V1t_2[1]; 
//    g2t[2][0] =  0.5*a0*V1t_0[2]; 
//    g2t[2][1] =  0.5*a1*V1t_1[2]; 
//    g2t[2][2] =  0.5*a2*V1t_2[2]; 
//    //translational DOF
//    for (int i=0;i<3;++i) 
//      { 
//        for (int k=0;k<3;++k) 
//          { 
//            int kIndex = i + (ndof_loc-1)*k;
//            BL_0[i][kIndex] = dhmat0[i][k];
//          }
//      }
//    for (int k=0;k<3;++k) 
//      { 
//        int kIndex0 =     (ndof_loc-1)*k;
//        int kIndex1 = 1 + (ndof_loc-1)*k;
//        int kIndex2 = 2 + (ndof_loc-1)*k;
//        BL_0[3][kIndex0] = dhmat0[1][k];
//        BL_0[3][kIndex1] = dhmat0[0][k];
//        BL_0[4][kIndex1] = dhmat0[2][k];
//        BL_0[4][kIndex2] = dhmat0[1][k];
//        BL_0[5][kIndex0] = dhmat0[2][k];
//        BL_0[5][kIndex2] = dhmat0[0][k];
//      }
//    //rotational DOF
//    for (int i=0;i<3;++i) 
//      { 
//        for (int k=0;k<3;++k) 
//          { 
//            int g1Index = 3 + (ndof_loc-1)*k;
//            int g2Index = 4 + (ndof_loc-1)*k;
//            BL_0[i][g1Index] = g1t[i][k] * G0[i][k];
//            BL_0[i][g2Index] = g2t[i][k] * G0[i][k];
//            BL_0[3][g1Index] = g1t[0][k] * G0[1][k] + g1t[1][k] * G0[0][k];
//            BL_0[3][g2Index] = g2t[0][k] * G0[1][k] + g2t[1][k] * G0[0][k];
//            BL_0[4][g1Index] = g1t[1][k] * G0[2][k] + g1t[2][k] * G0[1][k];
//            BL_0[4][g2Index] = g2t[1][k] * G0[2][k] + g2t[2][k] * G0[1][k];
//            BL_0[5][g1Index] = g1t[0][k] * G0[2][k] + g1t[2][k] * G0[0][k];
//            BL_0[5][g2Index] = g2t[0][k] * G0[2][k] + g2t[2][k] * G0[0][k];
//          }
//      }
//    //dumpMat(BL_0,"BL_0");
//    vector<vector<double> > BL_1 (6,vector<double>(15,0.));
//    vector<vector<double> > l    (3,vector<double>(3 ,0.));
//    vector<double> deltaVnt0_0(3,0.);
//    vector<double> deltaVnt0_1(3,0.);
//    vector<double> deltaVnt0_2(3,0.);
//    vector<vector<double> > deltaVnt0mat(3,vector<double>(3,0.));
//    vector<double> phi_11(3,0.);
//    vector<double> phi_12(3,0.);
//    vector<double> phi_13(3,0.);
//    vector<double> phi_21(3,0.);
//    vector<double> phi_22(3,0.);
//    vector<double> phi_23(3,0.);
//    for (int idir=0;idir<3;++idir) 
//      { 
//        deltaVnt0_0[idir] = Vnt_0[idir] - Vn0_0[idir];
//        deltaVnt0_1[idir] = Vnt_1[idir] - Vn0_1[idir];
//        deltaVnt0_2[idir] = Vnt_2[idir] - Vn0_2[idir];
//        deltaVnt0mat[idir][0] = deltaVnt0_0[idir];
//        deltaVnt0mat[idir][1] = deltaVnt0_1[idir];
//        deltaVnt0mat[idir][2] = deltaVnt0_2[idir];
//      }
//    //create ut,i = dut/dri for updating gti = g0_i + ut,i
//    //where ut = xt - x0
//    double dudx;
//    double dudy;
//    double dudz;
//    double dvdx;
//    double dvdy;
//    double dvdz;
//    double dwdx;
//    double dwdy;
//    double dwdz;
//    vector<double> u(3,0.);
//    vector<double> v(3,0.);
//    vector<double> w(3,0.);
//    u[0] = x0_t[0] - x0_0[0];
//    u[1] = x1_t[0] - x1_0[0];
//    u[2] = x2_t[0] - x2_0[0];
//    v[0] = x0_t[1] - x0_0[1];
//    v[1] = x1_t[1] - x1_0[1];
//    v[2] = x2_t[1] - x2_0[1];
//    w[0] = x0_t[2] - x0_0[2];
//    w[1] = x1_t[2] - x1_0[2];
//    w[2] = x2_t[2] - x2_0[2];
//    for (int k=0;k<3;++k)
//      {
//        dudx += dhmat0[0][k]*u[k] + 0.5*a0*(tPoint*dhmat0[0][k] + Jinv[0][2]*hmat[k])*deltaVnt0mat[0][k];
//        dudy += dhmat0[1][k]*u[k] + 0.5*a0*(tPoint*dhmat0[1][k] + Jinv[1][2]*hmat[k])*deltaVnt0mat[0][k];
//        dudz += dhmat0[2][k]*u[k] + 0.5*a0*(tPoint*dhmat0[2][k] + Jinv[2][2]*hmat[k])*deltaVnt0mat[0][k];
//        dvdx += dhmat0[0][k]*v[k] + 0.5*a0*(tPoint*dhmat0[0][k] + Jinv[0][2]*hmat[k])*deltaVnt0mat[1][k];
//        dvdy += dhmat0[1][k]*v[k] + 0.5*a0*(tPoint*dhmat0[1][k] + Jinv[1][2]*hmat[k])*deltaVnt0mat[1][k];
//        dvdz += dhmat0[2][k]*v[k] + 0.5*a0*(tPoint*dhmat0[2][k] + Jinv[2][2]*hmat[k])*deltaVnt0mat[1][k];
//        dwdx += dhmat0[0][k]*w[k] + 0.5*a0*(tPoint*dhmat0[0][k] + Jinv[0][2]*hmat[k])*deltaVnt0mat[2][k];
//        dwdy += dhmat0[1][k]*w[k] + 0.5*a0*(tPoint*dhmat0[1][k] + Jinv[1][2]*hmat[k])*deltaVnt0mat[2][k];
//        dwdz += dhmat0[2][k]*w[k] + 0.5*a0*(tPoint*dhmat0[2][k] + Jinv[2][2]*hmat[k])*deltaVnt0mat[2][k];
//      }
//    l[0][0] = dudx;
//    l[0][1] = dudy;
//    l[0][2] = dudz;
//    l[1][0] = dvdx;
//    l[1][1] = dvdy;
//    l[1][2] = dvdz;
//    l[2][0] = dwdx;
//    l[2][1] = dwdy;
//    l[2][2] = dwdz;
//    for (int k=0;k<3;++k)
//      {
//        for (int m=0;m<3;++m)
//          {
//            phi_11[k] += g1t[m][k]*l[m][0];
//            phi_12[k] += g1t[m][k]*l[m][1];
//            phi_13[k] += g1t[m][k]*l[m][2];
//            phi_21[k] += g2t[m][k]*l[m][0];
//            phi_22[k] += g2t[m][k]*l[m][1];
//            phi_23[k] += g2t[m][k]*l[m][2];
//          }
//      }
//    for (int k=0;k<3;++k)
//      {
//        int kIndex0 =     k*(ndof_loc-1);
//        int kIndex1 = 1 + k*(ndof_loc-1);
//        int kIndex2 = 2 + k*(ndof_loc-1);
//        BL_1[0][kIndex0] = l[0][0]*dhmat0[0][k];
//        BL_1[0][kIndex1] = l[1][0]*dhmat0[0][k];
//        BL_1[0][kIndex2] = l[2][0]*dhmat0[0][k];
//        BL_1[1][kIndex0] = l[0][1]*dhmat0[1][k];
//        BL_1[1][kIndex1] = l[1][1]*dhmat0[1][k];
//        BL_1[1][kIndex2] = l[2][1]*dhmat0[1][k];
//        BL_1[2][kIndex0] = l[0][2]*dhmat0[2][k];
//        BL_1[2][kIndex1] = l[1][2]*dhmat0[2][k];
//        BL_1[2][kIndex2] = l[2][2]*dhmat0[2][k];
//      }
//    for (int k=0;k<3;++k)
//      {
//        int kIndex3 = 3 + k*(ndof_loc-1);
//        int kIndex4 = 4 + k*(ndof_loc-1);
//        BL_1[0][kIndex3] = phi_11[k]*G0[0][k];
//        BL_1[0][kIndex4] = phi_21[k]*G0[0][k];
//        BL_1[1][kIndex3] = phi_12[k]*G0[1][k];
//        BL_1[1][kIndex4] = phi_22[k]*G0[1][k];
//        BL_1[2][kIndex3] = phi_13[k]*G0[2][k];
//        BL_1[2][kIndex4] = phi_23[k]*G0[2][k];
//      }            
//    for (int k=0;k<3;++k)
//      {
//        int kIndex0 =     k*(ndof_loc-1);
//        int kIndex1 = 1 + k*(ndof_loc-1);
//        int kIndex2 = 2 + k*(ndof_loc-1);
//        int kIndex3 = 3 + k*(ndof_loc-1);
//        int kIndex4 = 4 + k*(ndof_loc-1);
//        BL_1[3][kIndex0] = l[0][0]*dhmat0[1][k] + l[0][1]*dhmat0[0][k];
//        BL_1[3][kIndex1] = l[1][0]*dhmat0[1][k] + l[1][1]*dhmat0[0][k];
//        BL_1[3][kIndex2] = l[2][0]*dhmat0[1][k] + l[2][1]*dhmat0[0][k];
//        BL_1[3][kIndex3] = phi_12[k]*G0[0][k]   + phi_11[k]*G0[1][k];
//        BL_1[3][kIndex4] = phi_22[k]*G0[0][k]   + phi_21[k]*G0[1][k];
//        BL_1[4][kIndex0] = l[0][1]*dhmat0[2][k] + l[0][2]*dhmat0[1][k];
//        BL_1[4][kIndex1] = l[1][1]*dhmat0[2][k] + l[1][2]*dhmat0[1][k];
//        BL_1[4][kIndex2] = l[2][1]*dhmat0[2][k] + l[2][2]*dhmat0[1][k];
//        BL_1[4][kIndex3] = phi_13[k]*G0[1][k]   + phi_12[k]*G0[2][k];
//        BL_1[4][kIndex4] = phi_23[k]*G0[1][k]   + phi_22[k]*G0[2][k];
//        BL_1[5][kIndex0] = l[0][0]*dhmat0[2][k] + l[0][2]*dhmat0[0][k];
//        BL_1[5][kIndex1] = l[1][0]*dhmat0[2][k] + l[1][2]*dhmat0[0][k];
//        BL_1[5][kIndex2] = l[2][0]*dhmat0[2][k] + l[2][2]*dhmat0[0][k];
//        BL_1[5][kIndex3] = phi_13[k]*G0[0][k]   + phi_11[k]*G0[2][k];
//        BL_1[5][kIndex4] = phi_23[k]*G0[0][k]   + phi_21[k]*G0[2][k];
//      }            
//    //dumpMat(BL_1,"BL_1");
//    vector<vector<double> > B_xyz_tmp(6,vector<double>(15,0.));
//    vector<vector<double> > B_xyzNew (6,vector<double>(18,0.));
//    for (int i=0;i<6;++i) 
//      {
//        for (int j=0;j<15;++j)
//          {
//            B_xyz_tmp[i][j] = BL_0[i][j] + BL_1[i][j];
//          }
//      }
//    for (int i=0;i<6;++i) 
//      {
//        for (int j=0;j<5;++j)
//          {
//            B_xyzNew[i][j]    = B_xyz_tmp[i][j];
//            B_xyzNew[i][j+6]  = B_xyz_tmp[i][j+5];
//            B_xyzNew[i][j+12] = B_xyz_tmp[i][j+10];
//          }
//      }    
//    //dumpMat(B_xyz_tmp,"B_xyz_tmp");
//    //dumpMat(B_xyzNew,"B_xyzNew");
//    //wait();
//    vector<vector<double> > B_xyzNewT (18,vector<double>(6,0.));
//    vector<vector<double> > B_xyzNewTC(18,vector<double>(6,0.));
//    matrixTranspose(B_xyzNew,B_xyzNewT);
//    mm_mult(B_xyzNewT,CmatNL_,B_xyzNewTC);
//    //mm_mult(B_xyzNewTC,B_xyzNew,BTCBmat);

//    //dumpMat(dhmat0,"dhmat0");
//    //dumpMat(l,"l");
//    //dumpMat(phi_0,"phi_0");
//    //dumpMat(phi_1,"phi_1");
//    //dumpMat(phi_2,"phi_2");
//    //wait();

//    BNL2[0][3]  = dprod(duq1dr,gt_1);
//    BNL2[1][3]  = dprod(duq1ds,gt_2);
//    BNL2[2][3]  = dprod(duq1dt,gt_3);
//    BNL2[3][3]  = 0.5*(dprod(duq1dr,gt_2) + dprod(duq1ds,gt_1));
//    BNL2[4][3]  = 0.5*(dprod(duq1ds,gt_3) + dprod(duq1dt,gt_2));
//    BNL2[5][3]  = 0.5*(dprod(duq1dr,gt_3) + dprod(duq1dt,gt_1));
//    BNL2[0][4]  = dprod(duq1dr,gt_1);
//    BNL2[1][4]  = dprod(duq1ds,gt_2);
//    BNL2[2][4]  = dprod(duq1dt,gt_3);
//    BNL2[3][4]  = 0.5*(dprod(duq1dr,gt_2) + dprod(duq1ds,gt_1));
//    BNL2[4][4]  = 0.5*(dprod(duq1ds,gt_3) + dprod(duq1dt,gt_2));
//    BNL2[5][4]  = 0.5*(dprod(duq1dr,gt_3) + dprod(duq1dt,gt_1));
//    //alpha 2 and beta 2
//    BNL2[0][9 ] = dprod(duq2dr,gt_1);
//    BNL2[1][9 ] = dprod(duq2ds,gt_2);
//    BNL2[2][9 ] = dprod(duq2dt,gt_3);
//    BNL2[3][9 ] = 0.5*(dprod(duq2dr,gt_2) + dprod(duq2ds,gt_1));
//    BNL2[4][9 ] = 0.5*(dprod(duq2ds,gt_3) + dprod(duq2dt,gt_2));
//    BNL2[5][9 ] = 0.5*(dprod(duq2dr,gt_3) + dprod(duq2dt,gt_1));
//    BNL2[0][10] = dprod(duq2dr,gt_1);
//    BNL2[1][10] = dprod(duq2ds,gt_2);
//    BNL2[2][10] = dprod(duq2dt,gt_3);
//    BNL2[3][10] = 0.5*(dprod(duq2dr,gt_2) + dprod(duq2ds,gt_1));
//    BNL2[4][10] = 0.5*(dprod(duq2ds,gt_3) + dprod(duq2dt,gt_2));
//    BNL2[5][10] = 0.5*(dprod(duq2dr,gt_3) + dprod(duq2dt,gt_1));
//    //alpha 3 and beta 3
//    BNL2[0][15] = dprod(duq3dr,gt_1);
//    BNL2[1][15] = dprod(duq3ds,gt_2);
//    BNL2[2][15] = dprod(duq3dt,gt_3);
//    BNL2[3][15] = 0.5*(dprod(duq3dr,gt_2) + dprod(duq3ds,gt_1));
//    BNL2[4][15] = 0.5*(dprod(duq3ds,gt_3) + dprod(duq3dt,gt_2));
//    BNL2[5][15] = 0.5*(dprod(duq3dr,gt_3) + dprod(duq3dt,gt_1));
//    BNL2[0][16] = dprod(duq3dr,gt_1);
//    BNL2[1][16] = dprod(duq3ds,gt_2);
//    BNL2[2][16] = dprod(duq3dt,gt_3);
//    BNL2[3][16] = 0.5*(dprod(duq3dr,gt_2) + dprod(duq3ds,gt_1));
//    BNL2[4][16] = 0.5*(dprod(duq3ds,gt_3) + dprod(duq3dt,gt_2));
//    BNL2[5][16] = 0.5*(dprod(duq3dr,gt_3) + dprod(duq3dt,gt_1));
//    //CORRECTION
//    //create tying point vectors to fill BNL2 with using MITC scheme
//    vector<double> nrt1(18,0.);
//    vector<double> nrt3(18,0.);
//    vector<double> nst2(18,0.);
//    vector<double> nst3(18,0.);
//    vector<double> cTying(18,0.);
//    //alpha 1 and beta 1
//    nst2[3]   = 0.5*(dprod(duq1ds,gt_3_TP2) + dprod(duq1dt_TP2,gt_2));
//    nst2[4]   = 0.5*(dprod(duq1ds,gt_3_TP2) + dprod(duq1dt_TP2,gt_2));
//    nst3[3]   = 0.5*(dprod(duq1ds,gt_3_TP3) + dprod(duq1dt_TP3,gt_2));
//    nst3[4]   = 0.5*(dprod(duq1ds,gt_3_TP3) + dprod(duq1dt_TP3,gt_2));
//    nrt1[3]   = 0.5*(dprod(duq1dr,gt_3_TP1) + dprod(duq1dt_TP1,gt_1));
//    nrt1[4]   = 0.5*(dprod(duq1dr,gt_3_TP1) + dprod(duq1dt_TP1,gt_1));
//    nrt3[3]   = 0.5*(dprod(duq1dr,gt_3_TP3) + dprod(duq1dt_TP3,gt_1));
//    nrt3[4]   = 0.5*(dprod(duq1dr,gt_3_TP3) + dprod(duq1dt_TP3,gt_1));
//    nst2[3]   = 0.5*(dprod(duq1ds,gt_3_TP2) + dprod(duq1dt_TP2,gt_2));
//    nst2[4]   = 0.5*(dprod(duq1ds,gt_3_TP2) + dprod(duq1dt_TP2,gt_2));
//    nst3[3]   = 0.5*(dprod(duq1ds,gt_3_TP3) + dprod(duq1dt_TP3,gt_2));
//    nst3[4]   = 0.5*(dprod(duq1ds,gt_3_TP3) + dprod(duq1dt_TP3,gt_2));
//    nrt1[3]   = 0.5*(dprod(duq1dr,gt_3_TP1) + dprod(duq1dt_TP1,gt_1));
//    nrt1[4]   = 0.5*(dprod(duq1dr,gt_3_TP1) + dprod(duq1dt_TP1,gt_1));
//    nrt3[3]   = 0.5*(dprod(duq1dr,gt_3_TP3) + dprod(duq1dt_TP3,gt_1));
//    nrt3[4]   = 0.5*(dprod(duq1dr,gt_3_TP3) + dprod(duq1dt_TP3,gt_1));
//    //alpha 2 and beta 2
//    nst2[9 ]  = 0.5*(dprod(duq2ds,gt_3_TP2) + dprod(duq2dt_TP2,gt_2));
//    nst2[10]  = 0.5*(dprod(duq2ds,gt_3_TP2) + dprod(duq2dt_TP2,gt_2));
//    nst3[9 ]  = 0.5*(dprod(duq2ds,gt_3_TP3) + dprod(duq2dt_TP3,gt_2));
//    nst3[10]  = 0.5*(dprod(duq2ds,gt_3_TP3) + dprod(duq2dt_TP3,gt_2));
//    nrt1[9 ]  = 0.5*(dprod(duq2dr,gt_3_TP1) + dprod(duq2dt_TP1,gt_1));
//    nrt1[10]  = 0.5*(dprod(duq2dr,gt_3_TP1) + dprod(duq2dt_TP1,gt_1));
//    nrt3[9 ]  = 0.5*(dprod(duq2dr,gt_3_TP3) + dprod(duq2dt_TP3,gt_1));
//    nrt3[10]  = 0.5*(dprod(duq2dr,gt_3_TP3) + dprod(duq2dt_TP3,gt_1));
//    nst2[9 ]  = 0.5*(dprod(duq2ds,gt_3_TP2) + dprod(duq2dt_TP2,gt_2));
//    nst2[10]  = 0.5*(dprod(duq2ds,gt_3_TP2) + dprod(duq2dt_TP2,gt_2));
//    nst3[9 ]  = 0.5*(dprod(duq2ds,gt_3_TP3) + dprod(duq2dt_TP3,gt_2));
//    nst3[10]  = 0.5*(dprod(duq2ds,gt_3_TP3) + dprod(duq2dt_TP3,gt_2));
//    nrt1[9 ]  = 0.5*(dprod(duq2dr,gt_3_TP1) + dprod(duq2dt_TP1,gt_1));
//    nrt1[10]  = 0.5*(dprod(duq2dr,gt_3_TP1) + dprod(duq2dt_TP1,gt_1));
//    nrt3[9 ]  = 0.5*(dprod(duq2dr,gt_3_TP3) + dprod(duq2dt_TP3,gt_1));
//    nrt3[10]  = 0.5*(dprod(duq2dr,gt_3_TP3) + dprod(duq2dt_TP3,gt_1));
//    //alpha 3 and beta 3
//    nst2[15 ] = 0.5*(dprod(duq2ds,gt_3_TP2) + dprod(duq2dt_TP2,gt_2));
//    nst2[16]  = 0.5*(dprod(duq2ds,gt_3_TP2) + dprod(duq2dt_TP2,gt_2));
//    nst3[15]  = 0.5*(dprod(duq2ds,gt_3_TP3) + dprod(duq2dt_TP3,gt_2));
//    nst3[16]  = 0.5*(dprod(duq2ds,gt_3_TP3) + dprod(duq2dt_TP3,gt_2));
//    nrt1[15]  = 0.5*(dprod(duq2dr,gt_3_TP1) + dprod(duq2dt_TP1,gt_1));
//    nrt1[16]  = 0.5*(dprod(duq2dr,gt_3_TP1) + dprod(duq2dt_TP1,gt_1));
//    nrt3[15]  = 0.5*(dprod(duq2dr,gt_3_TP3) + dprod(duq2dt_TP3,gt_1));
//    nrt3[16]  = 0.5*(dprod(duq2dr,gt_3_TP3) + dprod(duq2dt_TP3,gt_1));
//    nst2[15]  = 0.5*(dprod(duq2ds,gt_3_TP2) + dprod(duq2dt_TP2,gt_2));
//    nst2[16]  = 0.5*(dprod(duq2ds,gt_3_TP2) + dprod(duq2dt_TP2,gt_2));
//    nst3[15]  = 0.5*(dprod(duq2ds,gt_3_TP3) + dprod(duq2dt_TP3,gt_2));
//    nst3[16]  = 0.5*(dprod(duq2ds,gt_3_TP3) + dprod(duq2dt_TP3,gt_2));
//    nrt1[15]  = 0.5*(dprod(duq2dr,gt_3_TP1) + dprod(duq2dt_TP1,gt_1));
//    nrt1[16]  = 0.5*(dprod(duq2dr,gt_3_TP1) + dprod(duq2dt_TP1,gt_1));
//    nrt3[15]  = 0.5*(dprod(duq2dr,gt_3_TP3) + dprod(duq2dt_TP3,gt_1));
//    nrt3[16]  = 0.5*(dprod(duq2dr,gt_3_TP3) + dprod(duq2dt_TP3,gt_1));

#endif
