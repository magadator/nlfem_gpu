/*
 * nlfem - Nonlinear FEM Solver
 * Author  : Gokul G. Anugrah
 * Contact : gokulanugrah@gmail.com
 *
 * fem.h — FEM base class declaration.
 *   Abstract base class for the parallel FEM solver. Provides all member
 *   variables, MPI decomposition interfaces, and the public API:
 *     Setup_fem()              — initialize and partition the mesh
 *     Step_structural_solver() — advance one time/load step
 *     Finalize_fem()           — clean up and finalize
 */

#ifndef KARMA_UTILS_FEM_H
#define KARMA_UTILS_FEM_H
//
#include <sstream>
#include <string>
#include <map>
#include "file_io.h"
#include "string_utilities.h"
#include <vector>
#include <cmath>
#include "options.h"
#if (USEPETSC>0)
#include "petsc.h"
#include "parmetis.h"
#endif
#include <iomanip>
#include <iostream>
#include "auxTransfer.h"

#include "mpi_ops.h"

#include <ctime>
#include "screen.h"
#include "linear_algebra.h"
#include "surfaceTriangulation.h"
#include <algorithm>
using namespace std;

//elements
#define SIMPLETRI   1
#define LMITC3      2
#define MITC3       3
#define BEAM        4
#define CABLE       5

//load/displacement transfer interpolation methods
#define PWC         1
#define WLSQR       2 
#define LSQRZ       3
#define WLSQRZ      4
#define WLSQRX      5

namespace KARMA {
  
  struct element {
    int globID;
    int ndof;
    int nnodes;
    int type;
    int compID;
  };
  
  class Mappings {
  public:
    std::map<int,int> nodeLocal2Global;
    std::map<int,int> elemLocal2Global;
    std::map<int,int> nodeLocal2Global_ext;
    std::map<int,int> elemLocal2Global_ext;
    std::map<int,int> globalElement2ProcMap;
  };
  
  class fem {
  protected:
  public:
    double R0;
    double D0;
    //workhorses
    vector<double> Ke_loc_;
    vector<double> Me_loc_;
    vector<double> Fvec_assemble_;
    vector<double> mass_vec_;
    vector<double> stiff_vec_;
    vector<double> unmat_vec_;
    vector<double> damp_vec_;
    //GPU element assembly (Phase 1: CABLE elements only, see src/gpu/)
    bool useGPUAssembly_;
    void gpuPrecomputeCableElements();
    //local-face-id -> index into GPU batch result, -1 if not a cable element
    vector<int>    gpuCableIndex_;
    vector<double> gpuCableKe_;
    vector<double> gpuCableMe_;
    vector<double> gpuCableFvec_;
    //MPI Vars
    int myElementCount_;
    int myNodeCount_;
    int myOffset_;
    vector<int> partitionOffset_;
    vector<int> partitionCellCounts_;
    int numprocs_;
    //
    double geomThickFactor_;
    double rampFactor_;
    vector<double> density0_;
    int numForceInterpPts_;
    void initializeAuxillaryBodyVelocity();
    vector<double> iniAuxNodeVel_;
    double nm_beta_;
    double nm_gamma_;
    int strucSkip_;
    //processor ID
    int mypeno_;
    bool one2one_;
    bool enforceContact_;
    //for restart
    bool do_restart_;
    bool revertNormalsByCompID_;
    vector<int> revertNormalsCompIDList_;
    string restartFilename_;
    string assistedMappingFilename_;
    int restartSkip_;
    int restartTimestep_;
    vector<vector<double> > restartVec_;
    void markNodesInBox(vector<int>& markedNodes);
    void readLinearRestart(int& restart_timestep,vector<vector<double> >& restartVec);
    void writeLinearRestart(vector<vector<double> >& restartVec, int& timeStep);
    void writeRestart(int& nt);
    void readRestart();
    //boolean
    bool read_normals_;
    bool is_defined_;
    bool is_unsteady_;
    bool is_nonlinear_;
    bool revertGeometryNormals_;
    bool MollerTrumboreIntersection(vector<double>& origin,
                                    vector<double>& normal,
                                    vector<double>& edge1,
                                    vector<double>& edge2,
                                    vector<double>& vertex);
    void dumpGeometryPartition(vector<int>& nodePart,
			       vector<int>& facePart);
    //strings
    string inputFile_;    
    string pctype_;
    string timeIntegrationMethod_;
    bool firstTimeComputingK_;
    bool firstTimeStep_;
    //functions
    fem(void);
    void assembleGlobKMat();
    virtual void comp_Ke_loc(vector<double>& mat,vector<double>& Fvec,const int& triID,const int& elType) = 0;
    virtual void comp_Me_loc(vector<double>& mat,const int& triID,const int& elType) = 0;
    void setForceSteady();
    void initializePetsc();
    void initializeDirVecs();
    void reInit();
    void solve(const string PCtype);
    void nonLinearSolve();
    //condensed main
    void Read_structural_input();
    void Step_structural_solver(int& nt);
    void Finalize_fem();
    void initializeDisplacementTransfer();
    void initializeLoadTransfer();
    void restructureFSI();
    void computeDisplacementTransferStencils();
    void computeLoadTransferStencils();
    void Setup_fem();
    void ini3Das2D();
    string surfaceFilename_;
    //parallelization
    vector<double> rootTotalNodes_;
    vector<double> rootTotalNodes0_;
    vector<int>    rootTotalFaces_;
    vector<int>    rootNdofPerNode_;
    vector<int>    rootNnodesPerEl_;
    vector<int>    rootNdofOffset_;
    bool parallelSolve_;
    Mappings elemMaps_;
    vector<element> elementData_;
    vector<element> elementData_ext_;
    void initializeDomainPartition(vector<int>& elementTypesList);
    void distributeElements();
    void partitionDomain_parmetis(vector<int>& elementTypesList_loc,
                                  vector<int>& rootComponentIDs);
    void define();
    void setupElements(int& nonUniqueNodeCount,
                       vector<int>& globalNodeList);
    void dumpPartition(vector<double>& Nodes,
                       vector<int>&    Faces);
    void finalize_fem();
    void treatSingleProcFEM(vector<int>& elementTypesList,
                            vector<int>& rootComponentIDs);
    void treatSingleProcFSI();
    void obtainUniqueAuxillaryMapping();
    void identifyAuxillaryParticipants(); 
    void initializeNormalVectors(int& nonUniqueNodeCount,
                                 vector<int>& globalNodeList,
                                 vector<int>& elementTypesList,
                                 vector<int>& temporary_nnodes_offset,
                                 vector<int>& rootComponentIDs);
    vector<int> myLocalNodesNormalVectorStatus_;//is -1 if not a normal vector node, otherwise is ID of normal vector nodes
    vector<int> myNormalVectorNodes_;//local node IDs of my nodes that need a normal vector
    int myNumberOfNormalVectorNodes_;//number of my nodes that need a normal vector
    vector<double> myDirVecs_;
    vector<double> myDirVecs0_;
    vector<double> myDirVecs_bak_;
    vector<int> myAuxParticipatingElements_;
    int myNumberOfAuxParticipatingElements_;
    vector<int> myAuxParticipatingNodes_;
    int myNumberOfAuxParticipatingNodes_;
    void identifyAuxCornerTreatmentNodes();
    bool auxCornerTreatment_;
    int myNumberOfAuxCornerTreatmentNodes_;
    vector<int> myAuxCornerTreatmentNodes_;
    void prepareVarsForParallelIO(vector<double>& localFaceForces,
				  vector<double>& localFaceDispls,
				  vector<int>&    nnodesVec);
    int ierr;
#if (USEPETSC>0)
    void identifyHaloElements(vector<int>& elementTypesList_loc,
			      vector<int>& nnodesOffset,
			      idx_t *part,
			      vector<int>& rootComponentIDs);
#endif
    vector<int> myHaloElementGlobIDs_;
    int myNumberOfHaloElements_;
    vector<int> faces_ext_;
    vector<double> nodes_ext_;
    vector<double> nodes0_ext_;
    vector<int> nnodesLocalOffset_ext_;
    int myNodeCount_ext_;
    int myElementCount_ext_;
    vector<int> ndofPerNode_ext_;
    vector<int> ndofLocalOffset_ext_;
    void updateGeometricNormals();
    //prestrained
    bool prestrained_;
    string prestrainedFilename_;
    //MITC3 tying
    vector<double>  ert1_;
    vector<double>  ert3_;
    vector<double>  est2_;
    vector<double>  est3_;
    vector<double>  cTying_;
    vector<double> deltaVnt0_0_;
    vector<double> deltaVnt0_1_;
    vector<double> deltaVnt0_2_;
    vector<double> ut1_;
    vector<double> ut2_;
    vector<double> ut3_;    
    vector<double> ut3_TP1_;
    vector<double> ut3_TP2_;
    vector<double> ut3_TP3_;
    //trapezoidal integration (unconditionally stable Newmark-Beta)
    void BetaNewmarkTrapezoidal(int& nt);
    void updateBetaNewmarkAccel();
    void compBetaNewmarkRHS();
    void updateBetaNewmarkVars();
    //updated Newmark-Beta functions have freedom to choose scehme constants
    void BetaNewmarkGeneral(int& nt);
    void updateBetaNewmarkAccel_updated();
    void compBetaNewmarkRHS_updated();
    void updateBetaNewmarkVars_updated();
    void updateBetaNewmarkVel_updated();
    //backtrack line search
    void BetaNewmarkGeneral_BacktrackLineSearch(int& nt);
    bool backtrack_;
    bool advanceSolution_;
    bool solAdvanced_;
    bool normGrew_;
    //linear unsteady functions
    void linearRungeKutta4(int& nt);
    void compRHS_RK4(int& rkss,
                     double& coef,
                     vector<vector<double> >& rhs_u,
                     vector<vector<double> >& rhs_p);
    void ini_uSolve();
    void uSolveCN2(int& nt);
    void compLoadVector(int& nt);
#if (USEPETSC>0)
    void checkNewtonConvergence(PetscReal& RHSnorm,double& Fvecnorm,ofstream& rhsInf,double& Unorm,PetscReal& ETOL);
#else
    void checkNewtonConvergence(double& RHSnorm,double& Fvecnorm,ofstream& rhsInf,double& Unorm,double& ETOL);
#endif
    void getStrains();
    void updateDirectorVecs();
    void unUpdateDirectorVecs();
    void centralDifference(int& nt);
    void updateCentralDifferenceVars();
    void compCentralDifferenceRHS();
    //loading
    void applyPressureLoading();
    void applyBodyForces();
    void getElementProps(vector<int>& markedFaces);
    void setForce(int& nonUniqueNodeCount,
                  vector<int>& globalNodeList);
    //node output
    bool recordNodeOutput_;
    vector<int> outputNodes_;
    vector<string> outputNodeDOF_;
    vector<int> outputNodeFreq_;
    //
    vector<vector<int> > elementsThatClaimMe_;
    vector<vector<int> > elementsThatClaimMe_ext_;
    vector<int> numberOfElementsThatClaimMe_;
    vector<int> numberOfElementsThatClaimMe_ext_;
    //internal force output
    vector<double> myCableTensionVector_;
    bool recordCableOutput_;
    vector<int> outputCables_;
    bool noForceFile_;
    vector<string> cablesForcesFile;
    vector<string> nodeOutFile;
    ofstream node_output[10];
    ofstream cable_output[100];
    int nIOCables_;
    void outputNodeData(int& nt);
    void outputCableData(int& nt);
    int myIOCables_[100];
    //
    void dumpFEMForces(int& nt);
    void dumpAUXForces(int& nt);
    void intermediateOutput(int& nt,vector<vector<double> >& mat1);
    void intermediateOutput(int& nt,vector<vector<double> >& mat1,vector<vector<double> >& mat2);
    //c++ objects
    vector<vector<vector<double> > > workingvectmp1_;
    vector<vector<vector<double> > > workingvectmp2_;
    vector<int> movingNodes_;
    ofstream rhsInf;
    ofstream tip;
    ofstream panel;
    int DIM;
    int xyz;
    vector<double> dirVec_;
    int dirVecIndex_;
    double newtonTolerance_;
    int it_;
    int nt_;
    double relaxationFactor_;
    vector<double> Area_;
    vector<double> nodeVel_; 
    vector<double> nodes_; 
    vector<double> nodes0_;
    //aux geometry
    vector<int> myAuxNodeNumVec_;
    vector<int> myAuxNodeNumVecXYZ_;
    vector<int> nodeIndexOffset_;
    vector<int> nodeXYZOffset_;
    int myAuxNodeCount_;
    int myFixedDOFCount_;
    vector<int> geometryMap_;
    vector<double> myAuxNodes_;
    vector<int> myAuxFaces_;
    double rayLength_;
    string auxForceInterpMethod_;
    vector<auxStencil> geom2fem;
    vector<auxStencil> fem2geom;
    string auxDispInterpMethod_;
    bool useAuxGeometry_;
    bool isFSI_;
    void threeAstwoLoadTransfer();
    void threeDLoadTransfer();
    void twoDLoadTransfer();
    void auxillaryTransfer_ini();
    void auxDispTransfer();
    void auxLoadTransfer();
    string forceInterpOrder_;
    string dispInterpOrder_;
    void getDistanceOrder(vector<int>& ordered,
                          vector<double>& dist);
    void getLSQRcoefs(vector<vector<double> >& mat,
                      int ID,
                      vector<double>& coefs);
    void getWLSQRcoefs(vector<vector<double> >& mat,
                       vector<double>& weights,
                       int ID,
                       vector<double>& coefs);
    bool shiftGeometry_;
    vector<double> shiftGeomXYZ_;
    vector<double> aux_nodes_;
    vector<double> aux_nodes0_;
    vector<double> aux_nodes_copy_;
    vector<int>    aux_faces_;    
    vector<int>    is_aux_node_;
    vector<double> auxGeomForces_;
    vector<int> myAuxLoadingNodes_;
    int myNumberOfAuxLoadingNodes_;
    string auxFilename_;
    int dimLoc_;
    int leastSQR_npnts_;
    int numGeomNodes_;
    int numGeomFaces_;
    vector<auxStencil> stencils_;
    vector<auxStencil> loadStencils_;
    //pressure loading
    bool pressureLoading_;
    int pressureDirection_;
    bool normalPressure_;
    double pressureTime_;
    double pressureMag_;
    //body forces
    bool bodyForces_;
    double bodyForceMag_;
    vector<double> gravity_;
    int bodyForceDirection_; //direction of EXTERNAL body force
    vector<double> elementBodyForce_;
    //
    int nof_effective_nodes_;
    int fsiSkip_;
    vector<vector<double> > strainGlob_; 
    vector<vector<double> > emat_rot_;
    vector<double> Fvec_;
    int ls_;
    vector<double> thickness_;
    //optimization    
    vector<vector<double> > centroidForce_;
    vector<int> faces_;
    int nof_faces_;
    vector<double> deltaUi_;
    vector<double> Uim1_;
    vector<double> Ui_;
    bool threeD_as_twoD_;
    int ntSkip_;
    //multiple elements
    int numberOfElementTypes_;      //how many types of elements are present
    vector<int> elementTypes_;      //which element types are present
    vector<int> myNumberOfEachElementType_;//how many of each type of element type do I have
    vector<int> ndofPerNode_;       //nDOF on my nodes
    vector<int> ndofGlobalOffset_;  //global node offsets on ndof
    vector<int> ndofLocalOffset_;   //local node offsets on ndof
    vector<int> nnodesGlobalOffset_;//global face node offsets --- probably never need
    vector<int> nnodesLocalOffset_; //local face node offsets
    vector<double> density_;
    int globalSystemSize_;
    int localSize_;
    int blockSize_;
    int nof_nodes_glob_;
    int nof_faces_glob_;
    vector<vector<double> > dirIJK_;
    vector<int> nodeBC_;
    vector<int> myFixedDOFIndexList_;
    vector<double> forceVec_; 
    string fem_force_file_;
    string fem_BC_file_;
    int numFixedNodes_;
    vector<int> nodeFixed_;
    vector<int> fixedDOF_;
    int numForcedNodes_;
    int numMomentNodes_;
    //regional BC
    bool loadsByCompID_;
    void applyBC(int& nonUniqueNodeCount,
                 vector<int>& globalNodeList);
    bool regionalBC_;
    string regionShape_;
    double regionDiameter_;
    string regionNormal_;
    vector<double> regionCenter_;
    int noIterations;
    //unsteady
    double dt_;
    int numTimesteps_;
    double time_;
    double time0_;
    // 
    double fillLevels; 
    int numLoadSteps_;
    int nonlinearIts_;
    //non linear director Vecs
    vector<double> Vn0_;
    vector<double> Vnt_;
    vector<double> Vnt_bak_;
    vector<double> Vnt_geom_;
    vector<double> V1t_;
    vector<double> V2t_;
    vector<double> V10_;
    vector<double> V20_;
    vector<double> deltaVnt_;
    vector<double> deltaVnt_bak_;
    //
    int itCount_;
    vector<vector<double> > GLstrainMat_;
    vector<vector<double> > PKSTmat_;
    //matrices
    vector<vector<double> > BTCBmat_;
    vector<vector<double> > KNL1mat_;
    vector<vector<double> > KNL2mat_;
    vector<vector<double> > dhmat_;
    vector<vector<double> > ematT_;
    //create displacement interpolation matrices
    vector<vector<double> > N_;
    vector<vector<double> > NT_;
    vector<vector<double> > NTN_;
    vector<vector<double> > N_rot_;
    vector<vector<double> > NT_rot_;
    vector<vector<double> > NTN_rot_;
    //covariant, contravariant vectors
    vector<double> gt_1_;
    vector<double> gt_2_;
    vector<double> gt_3_;
    vector<double> g0_1_;
    vector<double> g0_2_;
    vector<double> g0_3_;
    vector<double> g01_;
    vector<double> g02_;
    vector<double> g03_;
    vector<double> g0_3_TP1_;
    vector<double> g0_3_TP2_;
    vector<double> g0_3_TP3_;
    vector<double> gt_3_TP1_;
    vector<double> gt_3_TP2_;
    vector<double> gt_3_TP3_;
    vector<vector<double> > g0matT_;
    vector<vector<double> > gtmatT_;
    vector<vector<double> > gt_matT_;
    vector<double> GLstrains_;
    vector<vector<double> > strains_rst_;
    vector<double> PKSTvec_;
    vector<vector<double> > BNL1_;
    //internal stress vector
    vector<double> Fvec_loc_;
    //create global Catesian basis
    vector<double> ex_;
    vector<double> ey_;
    vector<double> ez_;
    //create element representations of director vectors
    vector<double> Vn0_0_;
    vector<double> Vn0_1_;
    vector<double> Vn0_2_;
    vector<double> Vnt_0_;
    vector<double> Vnt_1_;
    vector<double> Vnt_2_;
    vector<double> V1t_0_;
    vector<double> V1t_1_;
    vector<double> V1t_2_;
    vector<double> V2t_0_;
    vector<double> V2t_1_;
    vector<double> V2t_2_;
    vector<double> deltaVnt_0_;
    vector<double> deltaVnt_1_;
    vector<double> deltaVnt_2_;
    vector<int> face_nodes_;
    //current global node coordinates
    vector<double> x0_t_;
    vector<double> x1_t_;
    vector<double> x2_t_;
    //intial global node coordinates
    vector<double> x0_0_;
    vector<double> x1_0_;
    vector<double> x2_0_;
    //dummy vectors
    vector<double> dummy1_;
    vector<double> dummy2_;
    vector<double> dummy3_;
    //for rotational dof in MITC3 element
    vector<double> dk1dr_;
    vector<double> dk1ds_;
    vector<double> dk1dt_;
    vector<double> dk2dr_;
    vector<double> dk2ds_;
    vector<double> dk2dt_;
    vector<double> dk3dr_;
    vector<double> dk3ds_;
    vector<double> dk3dt_;
    vector<double> df1dr_;
    vector<double> df1ds_;
    vector<double> df1dt_;
    vector<double> df2dr_;
    vector<double> df2ds_;
    vector<double> df2dt_;
    vector<double> df3dr_;
    vector<double> df3ds_;
    vector<double> df3dt_;
    vector<double> dk1dt_TP1_;
    vector<double> dk2dt_TP1_;
    vector<double> dk3dt_TP1_;
    vector<double> df1dt_TP1_;
    vector<double> df2dt_TP1_;
    vector<double> df3dt_TP1_;
    vector<double> dk1dt_TP2_;
    vector<double> dk2dt_TP2_;
    vector<double> dk3dt_TP2_;
    vector<double> df1dt_TP2_;
    vector<double> df2dt_TP2_;
    vector<double> df3dt_TP2_;
    vector<double> dk1dt_TP3_;
    vector<double> dk2dt_TP3_;
    vector<double> dk3dt_TP3_;
    vector<double> df1dt_TP3_;
    vector<double> df2dt_TP3_;
    vector<double> df3dt_TP3_;
    vector<vector<double> > Cmat_;
    vector<vector<vector<double> > > CmatNL_;
    vector<vector<double> > deformationGradient_rst_;
    vector<vector<double> > deformationGradient_xyz_;
    vector<vector<double> > integratingCauchyST_1_;
    vector<vector<double> > integratingCauchyST_2_;
    vector<vector<vector<double> > > cauchyStressTensor_;
    vector<double> vonMisesStress_;
    double jacobian_;
    double inv_jacobian_;
    vector<vector<double> > Bmat_;
    vector<vector<double> > B_rst_;
    vector<vector<double> > B_xyz_;
    vector<vector<double> > B_xyzT_;
    vector<vector<double> > B_xyzTC_;
    vector<vector<double> > tensorprod_;
    //element properties
    double inPlaneElastic_;
    double outPlaneElastic_;
    vector<double> inPlaneEmod_;
    vector<double> outPlaneEmod_;
    vector<double> eModulus_;
    vector<double> poisRatio_;
    vector<double> MoI_;
    double defaultEModulus_;
    double defaultPoisRatio_;
    double defaultThickness_;
    vector<double> defaultDensity_;
    double defaultMoI_;
    double beamCSA_;
    double defaultWidth_;
    double defaultLength_;
    vector<double> HrTHr_;
    vector<double> xhat_;
    vector<double> uhat_;
    vector<double> BL_;
    vector<double> T_;
    double lambda;
    double lambdam1;
    //petsc objects
#if (USEPETSC>0)
    PetscScalar Tint_phi_,Tint_beta_,Tint_a1_,Tint_a2_,Tint_a3_,Tint_a4_,Tint_a5_,Tint_a6_;
    Mat Kmat_p_,Mmat_p_,mat_p_,Cmat_p_;
    Vec tempVec_p_,tempVec2_p_,tempVec3_p_;
    Vec U_p_,Up1_p_,Udot_p_,Udotp1_p_,Uddot_p_,Uddotp1_p_,Uk_p_,deltaUk_p_,deltaUkm1_p_,Udotp1_km1_p_,Uddotp1_km1_p_;
    Vec loadVector_p_,avec_p_;
    Vec rhs_p_;
    Vec UpM1_p_;
    Vec mySolutionVector_p_;
    Vec myLastSolutionVector_p_;
    Vec mySolutionVector_p_ext_;
    Vec myForceVector_p_;
    Vec forceVector_p_;
    Vec rootTotalSolution_p_;
    Vec myRestartVec_p_;
    PetscInt maxNoIts_;
    PetscScalar RelTol_,AbsTol_;
    KSP kspOp_p_;
    PC pc;
    IS myGlobalDOFIndexList_p_;
    IS myLocalDOFIndexList_p_;
    IS myGlobalDOFIndexList_p_ext_;
    IS myLocalDOFIndexList_p_ext_;
    VecScatter scatterToAll_p_; //strategic broadcast
    VecScatter scatterToAll_p_ext_; //strategic broadcast to extended map with ghosts
    VecScatter scatterToRoot_p_;//total dump to root
    PetscReal rhsNorm_,dispNorm_;
#else

    int maxNoIts_;
    double RelTol_,AbsTol_;
    double rhsNorm_,dispNorm_;

    /* double Tint_phi_,Tint_beta_,Tint_a1_,Tint_a2_,Tint_a3_,Tint_a4_,Tint_a5_,Tint_a6_; */
    /* vector<vector<double> > Kmat_p_,Mmat_p_,mat_p_,Cmat_p_; */
    /* vector<double> tempVec_p_,tempVec2_p_,tempVec3_p_; */
    /* vector<double> U_p_,Up1_p_,Udot_p_,Udotp1_p_,Uddot_p_,Uddotp1_p_,Uk_p_,deltaUk_p_,deltaUkm1_p_,Udotp1_km1_p_,Uddotp1_km1_p_; */
    /* vector<double> loadVector_p_,avec_p_; */
    /* vector<double> rhs_p_; */
    /* vector<double> UpM1_p_; */
    /* vector<double> mySolutionVector_p_; */
    /* vector<double> myLastSolutionVector_p_; */
    /* vector<double> mySolutionVector_p_ext_; */
    /* vector<double> myForceVector_p_; */
    /* vector<double> forceVector_p_; */
    /* vector<double> rootTotalSolution_p_; */
    /* vector<double> myRestartVec_p_; */

    /* KSP kspOp_p_; */
    /* PC pc; */
    /* IS myGlobalDOFIndexList_p_; */
    /* IS myLocalDOFIndexList_p_; */
    /* IS myGlobalDOFIndexList_p_ext_; */
    /* IS myLocalDOFIndexList_p_ext_; */
    /* VecScatter scatterToAll_p_; //strategic broadcast */
    /* VecScatter scatterToAll_p_ext_; //strategic broadcast to extended map with ghosts */
    /* VecScatter scatterToRoot_p_;//total dump to root */


#endif
  };
}
#endif
