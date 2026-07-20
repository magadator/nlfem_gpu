/*
 * nlfem - Nonlinear FEM Solver
 * Author  : Gokul G. Anugrah
 * Contact : gokulanugrah@gmail.com
 *
 * fem.cc — Core FEM class implementation.
 *   Implements parallel domain decomposition, global stiffness/mass matrix
 *   assembly, Newmark-Beta time integration (linear and nonlinear), PETSc
 *   linear system setup, restart I/O, FSI load/displacement transfer,
 *   and VTK output of deformed geometry and stress fields.
 */
#include <cstdlib>
#include <stdexcept>
#include "options.h"
#include "fem.h"
#include "../gpu/gpu_assembly.h"

#define OUTPUT_FORCE_POINT_CLOUD        0 //dump stencil cloud for force transfer
#define OUTPUT_DISPLACEMENT_POINT_CLOUD 0 //dump stencil cloud for displacement transfer
#define DUMP_KMAT_FROM_PETSC            0 //dump global stiffness matrix from petsc in parallel
#define TIME_ASSEMBLY                   0 //time global matrix assembly for one assembly
#define DUMP_GEOMETRY_PARTITION         1 //dump geometry pseudo-partition to vtk
#define DUMP_PRESSURE_DISTRIBUTION      0 //dump pressure distribution interpolated from geometry
#define MARK_NODES_IN_BOX               0 //1: parachute 24 gores, 2: 1-by-1 flag leading edge, 3: 80 gore parachute
#define RAY_INTERSECTION                0 //decompose geometry based on ray-triangle intersection test
#define FOUND_NODES_PENO               -1 //0 for root owning all non-identified geometry nodes --- only for non-participating geometry, i.e., static capsule etc
#define ASSISTED_MAPPING                0 //1 for using an alternate geometry to create the load and displacement mapping
#define OVERWRITE_FEM_NODES             0 //overwrite the FEM node locations after auxillary mappings are created --- to build a geometry around arbtrarily defined FEM meshes

using namespace std;
namespace KARMA {
  
  fem::fem(void) 
  {
    is_defined_ = false;
    //---- FIX: these members were declared in fem.h but never assigned.
    //     Reading them was undefined behaviour (see Setup_fem for the real values).
    DIM                   = 3;
    xyz                   = 3;
    numberOfElementTypes_ = 0;
    inPlaneElastic_       = 0.;
    outPlaneElastic_      = 0.;
    //GPU assembly on by default when built with USEGPU=1 and a device is
    //present (see gpuPrecomputeCableElements()); set NLFEM_NO_GPU=1 in the
    //environment to force the CPU path at runtime without recompiling,
    //e.g. for debugging or A/B comparison against the GPU path.
    useGPUAssembly_ = (getenv("NLFEM_NO_GPU") == nullptr);
    return;
  }
  
  void fem::define()
  {

    rhsNorm_ = 0.;
    dispNorm_ = 0.;
    //general size for myNodeCount*ndof where ndof varies with node type
    int size = 0;
    for (int n=0;n<myNodeCount_;++n) size += ndofPerNode_[n];

    myCableTensionVector_.resize(myElementCount_,0.);
    //my solution vectors
    Ui_.resize(size,0.);
    Uim1_.resize(size,0.);
    deltaUi_.resize(size,0.);
    //my nodal velocities
    nodeVel_.resize(myNodeCount_*xyz,0.);
    //area of my elements
    Area_.resize(myElementCount_,0.);
    //my interal work
    Fvec_.resize(size,0.);
    //Second Piola-Kirchoff stress and Green-Lagrange strain
    PKSTmat_    .resize(3,vector<double>(3,0.));
    GLstrainMat_.resize(3,vector<double>(3,0.));
    //optimization
    workingvectmp1_.resize(15,vector<vector<double> >(3,vector<double>(3,0.)));
    workingvectmp2_.resize(15,vector<vector<double> >(3,vector<double>(3,0.)));
    tensorprod_.resize(myElementCount_,vector<double>(81,0.));
    Cmat_.resize(6,vector<double>(6,0.));
    CmatNL_.resize(myElementCount_,vector<vector<double> >(6,vector<double>(6,0.)));
    deformationGradient_rst_.resize(3,vector<double>(3,0.));
    deformationGradient_xyz_.resize(3,vector<double>(3,0.));
    cauchyStressTensor_.resize(myElementCount_,vector<vector<double> >(3,vector<double>(3,0.)));
    integratingCauchyST_1_.resize(3,vector<double>(3,0.));
    integratingCauchyST_2_.resize(3,vector<double>(3,0.));
    vonMisesStress_.resize(myElementCount_,0.);
    //hard coded for MITC3
    Bmat_.resize   (6,vector<double>(15,0.));
    B_rst_.resize  (6,vector<double>(15,0.));
    B_xyz_.resize  (6,vector<double>(15,0.));
    B_xyzT_.resize (15,vector<double>(6,0.));
    B_xyzTC_.resize(15,vector<double>(6,0.));
    BTCBmat_.resize(15,vector<double>(15,0.));
    dhmat_.resize(3,vector<double>(3,0.));
    gt_1_.resize(3,0.);
    gt_2_.resize(3,0.);
    gt_3_.resize(3,0.);
    g0_1_.resize(3,0.);
    g0_2_.resize(3,0.);
    g0_3_.resize(3,0.);
    g01_.resize(3,0.);
    g02_.resize(3,0.);
    g03_.resize(3,0.);
    g0_3_TP1_.resize(3,0.);
    g0_3_TP2_.resize(3,0.);
    g0_3_TP3_.resize(3,0.);
    gt_3_TP1_.resize(3,0.);
    gt_3_TP2_.resize(3,0.);
    gt_3_TP3_.resize(3,0.);
    BNL1_.resize(9,vector<double>(15,0.));
    strains_rst_.resize(3,vector<double>(3,0.));
    g0matT_.resize(3,vector<double>(3,0.));
    gtmatT_.resize(3,vector<double>(3,0.));
    gt_matT_.resize(3,vector<double>(3,0.));
    //overallocated for MITC3 --- to big for cables
    //but want to avoid dynamic allocation
    Fvec_loc_.resize(15,0.);
    ex_.resize(3,0.);
    ey_.resize(3,0.);
    ez_.resize(3,0.);
    ex_[0] = 1.;
    ey_[1] = 1.;
    ez_[2] = 1.;
    ematT_.resize(3,vector<double>(3,0.));
    for (int i=0;i<3;++i) 
      { 
        ematT_[0][i] = ex_[i];
        ematT_[1][i] = ey_[i];
        ematT_[2][i] = ez_[i];
      }
    //from nonlinear LHS
    //KNL1mat_.resize(nnodes*systemDOF_,vector<double>(nnodes*systemDOF_,0.));
    //KNL2mat_.resize(nnodes*systemDOF_,vector<double>(nnodes*systemDOF_,0.));
    //for mass matrix, displacement interpolation matrices
    N_.resize(3,vector<double>(9,0.));
    NT_.resize(9,vector<double>(3,0.));
    NTN_.resize(9,vector<double>(9,0.));
    N_rot_.resize(3,vector<double>(6,0.));
    NT_rot_.resize(6,vector<double>(3,0.));
    NTN_rot_.resize(6,vector<double>(6,0.));
    Vn0_0_.resize(3,0.);
    Vn0_1_.resize(3,0.);
    Vn0_2_.resize(3,0.);
    Vnt_0_.resize(3,0.);
    Vnt_1_.resize(3,0.);
    Vnt_2_.resize(3,0.);
    V1t_0_.resize(3,0.);
    V1t_1_.resize(3,0.);
    V1t_2_.resize(3,0.);
    V2t_0_.resize(3,0.);
    V2t_1_.resize(3,0.);
    V2t_2_.resize(3,0.);
    deltaVnt_0_.resize(3,0.);
    deltaVnt_1_.resize(3,0.);
    deltaVnt_2_.resize(3,0.);    
    face_nodes_.resize(3);
    x0_t_.resize(3,0.);
    x1_t_.resize(3,0.);
    x2_t_.resize(3,0.);
    x0_0_.resize(3,0.);
    x1_0_.resize(3,0.);
    x2_0_.resize(3,0.);
    GLstrains_.resize(6,0.);
    PKSTvec_.resize(6,0.);
    dummy1_.resize(3,0.);
    dummy2_.resize(3,0.);
    dummy3_.resize(3,0.);
    dk1dr_.resize(3,0.);
    dk1ds_.resize(3,0.);
    dk1dt_.resize(3,0.);
    dk2dr_.resize(3,0.);
    dk2ds_.resize(3,0.);
    dk2dt_.resize(3,0.);
    dk3dr_.resize(3,0.);
    dk3ds_.resize(3,0.);
    dk3dt_.resize(3,0.);
    df1dr_.resize(3,0.);
    df1ds_.resize(3,0.);
    df1dt_.resize(3,0.);
    df2dr_.resize(3,0.);
    df2ds_.resize(3,0.);
    df2dt_.resize(3,0.);
    df3dr_.resize(3,0.);
    df3ds_.resize(3,0.);
    df3dt_.resize(3,0.);
    dk1dt_TP1_.resize(3,0.);
    dk2dt_TP1_.resize(3,0.);
    dk3dt_TP1_.resize(3,0.);
    df1dt_TP1_.resize(3,0.);
    df2dt_TP1_.resize(3,0.);
    df3dt_TP1_.resize(3,0.);
    dk1dt_TP2_.resize(3,0.);
    dk2dt_TP2_.resize(3,0.);
    dk3dt_TP2_.resize(3,0.);
    df1dt_TP2_.resize(3,0.);
    df2dt_TP2_.resize(3,0.);
    df3dt_TP2_.resize(3,0.);
    dk1dt_TP3_.resize(3,0.);
    dk2dt_TP3_.resize(3,0.);
    dk3dt_TP3_.resize(3,0.);
    df1dt_TP3_.resize(3,0.);
    df2dt_TP3_.resize(3,0.);
    df3dt_TP3_.resize(3,0.);
    ert1_.resize(15,0.);
    ert3_.resize(15,0.);
    est2_.resize(15,0.);
    est3_.resize(15,0.);
    cTying_.resize(15,0.);
    deltaVnt0_0_.resize(3,0.);
    deltaVnt0_1_.resize(3,0.);
    deltaVnt0_2_.resize(3,0.);
    ut1_.resize(3,0.);
    ut2_.resize(3,0.);
    ut3_.resize(3,0.);
    ut3_TP1_.resize(3,0.);
    ut3_TP2_.resize(3,0.);
    ut3_TP3_.resize(3,0.);
    elementBodyForce_.resize(size);
    //cable elements
    HrTHr_.resize(36,0.);
    xhat_.resize(6,0.);
    uhat_.resize(6,0.);
    BL_.resize(6,0.);
    T_.resize(36,0.);
    //always size local matrices for assembly to MITC3 size for now!
    Ke_loc_.resize(225,0.);
    Me_loc_.resize(225,0.);
    Fvec_assemble_.resize(15,0.);
    stiff_vec_.resize(12,0.);
    mass_vec_.resize(12,0.);
    damp_vec_.resize(12,0.);
    unmat_vec_.resize(12,0.);

    is_defined_          = true;
    firstTimeComputingK_ = true;
    firstTimeStep_       = true;

    //create a list of elements that own each of my local nodes within my partition
    //only considering tris for now! don't have a need for a more general list yet

    numberOfElementsThatClaimMe_.resize(myNodeCount_,0.);
    for (int n=0;n<myNodeCount_;++n)
      {
	elementsThatClaimMe_.push_back(vector<int>());
	for (int e=0;e<myElementCount_;++e)
	  {
	    if (elementData_[e].nnodes == 3)
	      {
		int nnodesOffset = nnodesLocalOffset_[e];

		if ( faces_[0 + nnodesOffset] == n ||
		     faces_[1 + nnodesOffset] == n || 
		     faces_[2 + nnodesOffset] == n )
		  {
		    numberOfElementsThatClaimMe_[n] += 1;
		    elementsThatClaimMe_[n].push_back(e);
		  }		  
	      }
	  }
      }

    //create a list of elements that own each of my local nodes including halo elements
    numberOfElementsThatClaimMe_ext_.resize(myNodeCount_ext_);
    for (int n=0;n<myNodeCount_ext_;++n)
      {
	elementsThatClaimMe_ext_.push_back(vector<int>());
	for (int e=0;e<myElementCount_ext_;++e)
	  {
	    if (elementData_ext_[e].nnodes == 3)
	      {
		int nnodesOffset_ext = nnodesLocalOffset_ext_[e];

		if ( faces_ext_[0 + nnodesOffset_ext] == n ||
		     faces_ext_[1 + nnodesOffset_ext] == n || 
		     faces_ext_[2 + nnodesOffset_ext] == n )
		  {
		    numberOfElementsThatClaimMe_ext_[n] += 1;
		    elementsThatClaimMe_ext_[n].push_back(e);
		  }		  
	      }
	  }
      }

  }

  void fem::initializePetsc()
  {

#if (USEPETSC>0)

    if (mypeno_ == 0) cout << "----------------------------" << endl;
    PetscInt alloc = 1000;

    if (parallelSolve_)
      {
	Screen::MasterInfo("Initializing PETSc - parallel");
#if (DEDICATED_FEM_PROC == 0)
        KSPCreate(PETSC_COMM_WORLD,&kspOp_p_);
        MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DETERMINE,PETSC_DETERMINE,globalSystemSize_,globalSystemSize_,alloc,NULL,alloc,NULL,&Kmat_p_);
#else
        KSPCreate(PETSC_COMM_SELF,&kspOp_p_);
        MatCreateAIJ(PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,globalSystemSize_,globalSystemSize_,alloc,NULL,alloc,NULL,&Kmat_p_);
#endif

        MatSetOption(Kmat_p_,MAT_SYMMETRIC,PETSC_TRUE);
        MatSetOption(Kmat_p_,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

        if (is_unsteady_)
          {
            //MatCreateBAIJ(PETSC_COMM_WORLD,blockSize_,myNodeCount_*systemDOF_,myNodeCount_*systemDOF_,nof_nodes_glob_*systemDOF_,nof_nodes_glob_*systemDOF_,alloc,NULL,alloc,NULL,& mat_p_);
            //MatCreateBAIJ(PETSC_COMM_WORLD,blockSize_,myNodeCount_*systemDOF_,myNodeCount_*systemDOF_,nof_nodes_glob_*systemDOF_,nof_nodes_glob_*systemDOF_,alloc,NULL,alloc,NULL,&Mmat_p_);
#if (DEDICATED_FEM_PROC == 0)
            MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DETERMINE,PETSC_DETERMINE,globalSystemSize_,globalSystemSize_,alloc,NULL,alloc,NULL,& mat_p_);
            MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DETERMINE,PETSC_DETERMINE,globalSystemSize_,globalSystemSize_,alloc,NULL,alloc,NULL,&Mmat_p_);
            MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DETERMINE,PETSC_DETERMINE,globalSystemSize_,globalSystemSize_,alloc,NULL,alloc,NULL,&Cmat_p_);
#else
            MatCreateAIJ(PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,globalSystemSize_,globalSystemSize_,alloc,NULL,alloc,NULL,& mat_p_);
            MatCreateAIJ(PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,globalSystemSize_,globalSystemSize_,alloc,NULL,alloc,NULL,&Mmat_p_);
            MatCreateAIJ(PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,globalSystemSize_,globalSystemSize_,alloc,NULL,alloc,NULL,&Cmat_p_);
#endif

            MatSetOption( mat_p_,MAT_SYMMETRIC,PETSC_TRUE);
            MatSetOption( mat_p_,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
            MatSetOption(Mmat_p_,MAT_SYMMETRIC,PETSC_TRUE);
            MatSetOption(Mmat_p_,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

            MatSetOption(Cmat_p_,MAT_SYMMETRIC,PETSC_TRUE);
            MatSetOption(Cmat_p_,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

          }
      }
    else 
      {

	Screen::MasterInfo("Initializing PETSc - serial");
        if (numberOfElementTypes_ > 1) Screen::MasterError("Multiple elements selected with block matrices");

        KSPCreate(PETSC_COMM_SELF,&kspOp_p_);

        MatCreateSeqBAIJ(PETSC_COMM_SELF,ndofPerNode_[0],globalSystemSize_,globalSystemSize_,alloc,NULL,&Kmat_p_);

        MatSetOption(Kmat_p_,MAT_SYMMETRIC,PETSC_TRUE);
        MatSetOption(Kmat_p_,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

        if (is_unsteady_)
          {
            MatCreateSeqBAIJ(PETSC_COMM_SELF,ndofPerNode_[0],globalSystemSize_,globalSystemSize_,alloc,NULL,&Mmat_p_);
            MatCreateSeqBAIJ(PETSC_COMM_SELF,ndofPerNode_[0],globalSystemSize_,globalSystemSize_,alloc,NULL,& mat_p_);

            MatSetOption(Mmat_p_,MAT_SYMMETRIC,PETSC_TRUE);
            MatSetOption(Mmat_p_,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

            MatSetOption( mat_p_,MAT_SYMMETRIC,PETSC_TRUE);
            MatSetOption( mat_p_,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

          }
      }

    MatCreateVecs  (Kmat_p_,&deltaUk_p_,&loadVector_p_);   //'incremental' solution vector and external-internal loading
    MatCreateVecs  (Kmat_p_,&deltaUkm1_p_,&forceVector_p_);//backup incremental solution vector, external force vector
    MatCreateVecs  (Kmat_p_,&Uk_p_     ,NULL);             //'total' solution vector

    if (is_unsteady_)
      {
        MatCreateVecs(Kmat_p_,&avec_p_        ,NULL);
        MatCreateVecs(Kmat_p_,&tempVec_p_     ,NULL);
        MatCreateVecs(Kmat_p_,&tempVec2_p_    ,NULL);
        MatCreateVecs(Kmat_p_,&tempVec3_p_    ,NULL);
        MatCreateVecs(Kmat_p_,&U_p_           ,&rhs_p_);
        MatCreateVecs(Kmat_p_,&Up1_p_         ,NULL);
        MatCreateVecs(Kmat_p_,&Udot_p_        ,NULL);
        MatCreateVecs(Kmat_p_,&Udotp1_p_      ,NULL);
        MatCreateVecs(Kmat_p_,&Udotp1_km1_p_  ,NULL);	
        MatCreateVecs(Kmat_p_,&Uddot_p_       ,NULL);
        MatCreateVecs(Kmat_p_,&Uddotp1_p_     ,NULL);
        MatCreateVecs(Kmat_p_,&Uddotp1_km1_p_ ,NULL);
      }
    //create solution vector for the root so it can hold the entire solution vector
    //aka, the root holds his local solution vector and the global solution vector
    if (mypeno_ == 0)
      {
        VecCreateSeq(MPI_COMM_SELF,globalSystemSize_,&rootTotalSolution_p_);
        VecZeroEntries(rootTotalSolution_p_);
      }

    int mysize = 0;
    for (int n=0;n<myNodeCount_;++n) mysize += ndofPerNode_[n];
    localSize_ = mysize;
    //create petsc local solution vector on each proc
    VecCreateSeq(MPI_COMM_SELF,mysize,&mySolutionVector_p_);
    VecZeroEntries(mySolutionVector_p_);
    VecCreateSeq(MPI_COMM_SELF,mysize,&myForceVector_p_);
    VecZeroEntries(myForceVector_p_);
    VecCreateSeq(MPI_COMM_SELF,mysize,&myLastSolutionVector_p_);
    VecZeroEntries(myLastSolutionVector_p_);
    //create myRestartVector for quickly dumping PETSc vars to restart
    VecCreateSeq(MPI_COMM_SELF,mysize,&myRestartVec_p_);
    VecZeroEntries(myRestartVec_p_);

#if (DEDICATED_FEM_PROC == 0)
    //fill local and global node displacement indices for petsc solution vector strategic scatter
    int myLocalDOFIndexList [mysize];
    int myGlobalDOFIndexList[mysize];
    int counter = 0;
    for (int n=0;n<myNodeCount_;++n)
      {
        int ndof_loc = ndofPerNode_[n];
        for (int dof=0;dof<ndof_loc;++dof)
          {
            myLocalDOFIndexList [counter] = counter;
            myGlobalDOFIndexList[counter] = dof + ndofGlobalOffset_[n];
            counter += 1;
          }
      }
    //petsc object index list for local solution vector
    ISCreateGeneral(MPI_COMM_SELF,mysize,myLocalDOFIndexList ,PETSC_COPY_VALUES,&myLocalDOFIndexList_p_ );
    ISCreateGeneral(MPI_COMM_SELF,mysize,myGlobalDOFIndexList,PETSC_COPY_VALUES,&myGlobalDOFIndexList_p_);

    //create scatter context for strategic broadcast of solution vector
    VecScatterCreate(deltaUk_p_,myGlobalDOFIndexList_p_,mySolutionVector_p_,myLocalDOFIndexList_p_,&scatterToAll_p_);
    //create scatter context for solution vector dumping to root
    VecScatterCreateToZero(Uk_p_,&scatterToRoot_p_,&rootTotalSolution_p_);

    //this way delayed for timing, calling this here for ghost solution vector
    VecScatterCreate(Uk_p_,myGlobalDOFIndexList_p_ext_,mySolutionVector_p_ext_,myLocalDOFIndexList_p_ext_,&scatterToAll_p_ext_);

#endif

    //VecSetBlockSize(avec_p_        ,blockSize_);
    //VecSetBlockSize(loadVector_p_  ,blockSize_);
    //VecSetBlockSize(U_p_           ,blockSize_);
    //VecSetBlockSize(Up1_p_         ,blockSize_);
    //VecSetBlockSize(Udot_p_        ,blockSize_);
    //VecSetBlockSize(Udotp1_p_      ,blockSize_);
    //VecSetBlockSize(Uddot_p_       ,blockSize_);
    //VecSetBlockSize(Uddotp1_p_     ,blockSize_);
    //VecSetBlockSize(tempVec_p_     ,blockSize_);
    //VecSetBlockSize(tempVec2_p_    ,blockSize_);
    //VecSetBlockSize(tempVec3_p_    ,blockSize_);
    //VecSetBlockSize(Uk_p_          ,blockSize_);
    //VecSetBlockSize(deltaUk_p_     ,blockSize_);
    //VecSetBlockSize(rhs_p_         ,blockSize_);    

    VecSet(loadVector_p_,0.);
    VecSet(Uk_p_        ,0.);
    VecSet(deltaUk_p_   ,0.);

    if (is_unsteady_)
      {
        VecSet(U_p_      ,0.);
        VecSet(Up1_p_    ,0.);
        VecSet(Udot_p_   ,0.);
        VecSet(Udotp1_p_ ,0.);
        VecSet(Udotp1_km1_p_ ,0.);
        VecSet(Uddot_p_  ,0.);
        VecSet(Uddotp1_p_,0.);
        VecSet(Uddotp1_km1_p_,0.);
        VecSet(avec_p_      ,0.);
        VecSet(tempVec_p_   ,0.);
        VecSet(tempVec2_p_  ,0.);
        VecSet(tempVec3_p_  ,0.);
        VecSet(rhs_p_       ,0.);
        //MatSetBlockSize(Kmat_p_,blockSize_);
        //MatSetBlockSize(Mmat_p_,blockSize_);
        //MatSetBlockSize( mat_p_,blockSize_);
        MatZeroEntries(Kmat_p_);
        MatZeroEntries(mat_p_ );
        MatZeroEntries(Mmat_p_);
        MatZeroEntries(Cmat_p_);
        KSPSetOperators(kspOp_p_,mat_p_,mat_p_);
      }
    else
      {
        //MatSetBlockSize(Kmat_p_,blockSize_);
        MatZeroEntries(Kmat_p_);
        KSPSetOperators(kspOp_p_,Kmat_p_,Kmat_p_);
      }
    //set KSP options
    if (parallelSolve_)
      {
#if (PETSC_HAVE_MUMPS)
        //PETSC with MUMPS uses parallel, direct LU to solve Ku = F --- expensive, but good for stiff FEM systems
        Screen::MasterInfo("Setting up PETSC with MUMPS for single step LU");
        PetscBool flg_mumps_lu = PETSC_TRUE; //default to lu
        PetscBool flg_mumps_ch = PETSC_FALSE;//cholesky only works with SBAIJ matrices

#if PETSC_VERSION_LT(3,7,0)
        PetscOptionsGetBool(NULL,"-use_mumps_lu",&flg_mumps_lu,NULL);
        PetscOptionsGetBool(NULL,"-use_mumps_ch",&flg_mumps_ch,NULL);
#else
        PetscOptionsGetBool(NULL,NULL,"-use_mumps_lu",&flg_mumps_lu,NULL);
        PetscOptionsGetBool(NULL,NULL,"-use_mumps_ch",&flg_mumps_ch,NULL);
#endif

	//select preconditioner and set as exact solve
	KSPSetType(kspOp_p_,KSPPREONLY);
	KSPGetPC(kspOp_p_,&pc);
	if (flg_mumps_lu)
	  {
	    PCSetType(pc,PCLU);
	  }
	else if (flg_mumps_ch)
	  {
	    MatSetOption(Kmat_p_,MAT_SPD,PETSC_TRUE);
	    PCSetType(pc,PCCHOLESKY);
	  }

	PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
        PCFactorSetUpMatSolverPackage(pc);

        PetscInt icntl, ival;
        if (is_unsteady_)
          {
	    //PCFactorGetMatrix(pc,&mat_p_);
	    //	    icntl = 7;ival = 3;
	    //            MatMumpsSetIcntl(mat_p_,icntl,ival);
	    icntl = 14;ival = 100000;
            MatMumpsSetIcntl(mat_p_,icntl,ival);
            MatMumpsSetIcntl(Kmat_p_,icntl,ival);
            MatMumpsSetIcntl(Mmat_p_,icntl,ival);
            // MUMPS memory pre-estimation: avoids INFO(1)=-9 on first factorization
            icntl = 33;ival = 1;
            MatMumpsSetIcntl(mat_p_,icntl,ival);
            MatMumpsSetIcntl(Kmat_p_,icntl,ival);
            MatMumpsSetIcntl(Mmat_p_,icntl,ival);
          }
        else
          {
	    //PCFactorGetMatrix(pc,&Kmat_p_);
	    //icntl = 7;ival = 3;
            //MatMumpsSetIcntl(Kmat_p_,icntl,ival);
	    icntl = 14;ival = 100000;
            MatMumpsSetIcntl(Kmat_p_,icntl,ival);
            // MUMPS memory pre-estimation: avoids INFO(1)=-9 on first factorization
            icntl = 33;ival = 1;
            MatMumpsSetIcntl(Kmat_p_,icntl,ival);
          }

#else
        //if PETSC with MUMPS is not installed, default to iterative schemes --- bad for stiff FEM systems!!!
        Screen::MasterWarning("MUMPS not installed, defaulting to iterative schemes");
        KSPSetTolerances  (kspOp_p_,RelTol_,AbsTol_,1.E15,maxNoIts_);
	KSPSetType        (kspOp_p_,KSPBICG);
        KSPGetPC          (kspOp_p_,&pc);
#if (DEDICATED_FEM_PROC == 1)
        PCSetType(pc,PCBJACOBI);
#else
	PCSetType(pc,PCBJACOBI);
#endif
        KSPSetInitialGuessNonzero(kspOp_p_,PETSC_TRUE);
        KSPSetFromOptions(kspOp_p_);
#endif
      }
    else
      {
        //serial runs default to PETSC LU preconditioning - default KSP combo
        KSPGetPC(kspOp_p_,&pc);
        PCSetType(pc,PCLU);
        KSPSetInitialGuessNonzero(kspOp_p_,PETSC_TRUE);
      }

#endif

  }
  
  void fem::finalize_fem()
  {

#if (USEPETSC>0)

    if (is_unsteady_)
      {
        if (mypeno_ == 0) cout << "---------------------------------" << endl;   
        Screen::MasterInfo("Structural Time Integration Complete");
        if (mypeno_ == 0) cout << "---------------------------------" << endl;   
        if (mypeno_ == 0) cout << endl;    

        if (mypeno_ == 0)
          {
            rhsInf.close();
            tip.close();
            for (int n=0;n<outputNodes_.size();++n) node_output[n].close();	    
          }
	if (nIOCables_ > 0) for (int c=0;c<nIOCables_;++c) cable_output[c].close();
      }
    else
      {
        for (int n=0;n<outputNodes_.size();++n) node_output[n].close();
	if (nIOCables_ > 0) for (int c=0;c<nIOCables_;++c) cable_output[c].close();
        if (mypeno_ == 0) cout << endl;
        if (mypeno_ == 0) cout << "---------------------------------" << endl;           
        Screen::MasterInfo("Static FEM Solution Obtained");
        if (mypeno_ == 0) cout << "---------------------------------" << endl;   
        if (mypeno_ == 0) cout << endl;
      }


#endif

  }
  
  void fem::getElementProps(vector<int>& markedFaces)
  {

#if (USEPETSC>0)

    if (threeD_as_twoD_)
      {
        //assuming z-comp is always into the plane
        defaultWidth_  = abs(nodes_[2+(nof_effective_nodes_)*xyz]-nodes_[2]);
        //for 3Das2D simulations, element length should be constant in x, y,or z (depending on orientation)
        //for 2D simuatlions, element length is automatically computing in stiffness matrix computation
        //for 3D simulations, these are not needed, the element area is used instead        
        if (dirVecIndex_ == 0)
          {
            defaultLength_ = abs(nodes_[4]-nodes_[1]);
          }
        else if (dirVecIndex_ == 1)
          {
            defaultLength_ = abs(nodes_[3]-nodes_[0]);
          }
      }

    //FIX: use assign() not resize() --- these vectors already hold the single
    //value parsed from the input deck, so resize() would leave entry 0 alone and
    //fill 1..myElementCount_-1 with the (previously uninitialized) fill value.
    eModulus_ .assign(myElementCount_,defaultEModulus_ );
    //for nonlinear MITC3 -- direction dependent elastic moduli
    inPlaneEmod_ .assign(myElementCount_, inPlaneElastic_);
    outPlaneEmod_.assign(myElementCount_,outPlaneElastic_);

    poisRatio_.assign(myElementCount_,defaultPoisRatio_);
    thickness_.assign(myElementCount_,defaultThickness_);

    vector<int> myfaces(3,0);
    vector<int> theirfaces(3,0);

    vector<int> markedNeighbors;

    for (int e=0;e<myElementCount_;++e)
      {
        if (elementData_[e].type == MITC3)
          {
            for (int i=0;i<markedFaces.size();++i)
              {
                int markedE = markedFaces[i];

                if (e == markedE)
                  {
                    thickness_[e] = 4.*defaultThickness_;
                  }
#if (MARK_NODES_IN_BOX == 1)
                else
                  {

                    int my_nnodes_offset    = nnodesLocalOffset_[e];
                    int their_nnodes_offset = nnodesLocalOffset_[markedE];
                    for (int n=0;n<3;++n)
                      {
                        myfaces[n] = faces_[my_nnodes_offset + n];
                        theirfaces[n] = faces_[their_nnodes_offset + n];
                      }

                    for (int n=0;n<3;++n)
                      {
                        int mynode = myfaces[n];
                        if ( find(theirfaces.begin(),theirfaces.end(),mynode) != theirfaces.end() )
                          {
                            if ( find(markedFaces.begin(),markedFaces.end(),e) == markedFaces.end() )
                              {
                                thickness_[e] = 1.*defaultThickness_;
                                markedNeighbors.push_back(e);
                                n = 10;
                              }
                          }
                      }
                  }
#endif
              }
          }
      }

    //Screen::MasterInfo("Outputting marked neighbors");
    //if (markedNeighbors.size() > 0)
    //  {
    //	string filename = "debugFEM/markFaces_neighbors_"+to_string(mypeno_)+".vtk";
    //	ofstream file;
    //	file.open(filename);
    //	file.precision(10);
    //	file.setf(ios::fixed);
    //	file.setf(ios::showpoint);
    //	file << "# vtk DataFile Version 3.0" << endl;
    //	file << "vtk output" << endl;;
    //	file << "ASCII" << endl;
    //	file << "DATASET POLYDATA" << endl;
    //	file << "POINTS " << markedNeighbors.size() << " double" << endl;
    //	for (int e=0;e<markedNeighbors.size();++e)
    //	  {
    //	    int element       = markedNeighbors[e];
    //	    int nnodes_loc    = elementData_[element].nnodes;
    //	    int nnodes_offset = nnodesLocalOffset_[element];
    //	    double x = 0;
    //	    double y = 0;
    //	    double z = 0;
    //	    for (int n=0;n<nnodes_loc;++n)
    //	      {
    //		int node = faces_[n + nnodes_offset];
    //		x += 0.333333333 * nodes_[node*xyz + 0];
    //		y += 0.333333333 * nodes_[node*xyz + 1];
    //		z += 0.333333333 * nodes_[node*xyz + 2];
    //	      }
    //	    file << x+1.e-10 << " " << y+1.e-10 << " " << z+1.e-10 << endl;
    //	  }
    //	file << "POINT_DATA " << markedNeighbors.size() << endl;
    //	file << "SCALARS Data double" << endl;
    //	file << "LOOKUP_TABLE default" << endl;
    //	for (int e=0;e<markedNeighbors.size();++e)
    //	  {
    //	    file << 1.000000 << endl;
    //	  }
    //	file.close();
    //  }

    density_  .resize(myElementCount_);
    density0_ .resize(myElementCount_);

    if (numberOfElementTypes_ > 1)
      {
        for (int e=0;e<myElementCount_;++e)
          {
            if (elementData_[e].type == CABLE)
              {
                density_ [e] = defaultDensity_[0];
                density0_[e] = defaultDensity_[0];
              }
            else if (elementData_[e].type == MITC3)
              {
                density_ [e] = defaultDensity_[1];
                density0_[e] = defaultDensity_[1];
              }
          }
      }
    else
      {
        for (int e=0;e<myElementCount_;++e)
          {
            density_ [e] = defaultDensity_[0];
            density0_[e] = defaultDensity_[0];
          }
      }
    MoI_      .resize(myElementCount_,defaultMoI_);
    //cables
    fill(Area_.begin(),Area_.end(),beamCSA_);

#endif

  }

  void fem::setForce(int& nonUniqueNodeCount,
                     vector<int>& globalNodeList)
  {

#if (USEPETSC>0)

    Screen::MasterInfo("Setting Forces");
    //smart allocation
    int sendcount = 0;
    for (int mn=0;mn<myNodeCount_;++mn)
      {
        sendcount += ndofPerNode_[mn];
      }
    forceVec_.resize(sendcount,0.);

#if (DEDICATED_FEM_PROC == 0)

    vector<int> nodeForceID;
    vector<int> forceDIR;
    vector<double> nodeForce;

    vector<int> nodeMomentID;
    vector<int> momentDIR;
    vector<double> nodeMoment;

    if (mypeno_ == 0)
      {
        //read loading from file
        ifstream force;
        force.open(fem_force_file_);
        noForceFile_ = force.fail();
        if (noForceFile_ && !pressureLoading_ && !bodyForces_ && !isFSI_)
          {
            if (useAuxGeometry_)
              { 
                Screen::MasterWarning("No loading specified, but useAuxGeometry_ is on");
              }
            else
              {
		
                Screen::MasterError("File not found for structural loading conditions, and no other loading options have been specified");
              }
          }
        else if (noForceFile_ && (pressureLoading_||bodyForces_||isFSI_))
          {
            Screen::MasterInfo("No structural loading from file detected, but loading is provided from");
            cout << "    Pressure Loading = " << pressureLoading_ << endl;
            cout << "    bodyForces_      = " << bodyForces_  << endl;
            cout << "    isFSI_           = " << isFSI_ << endl;
          }
        else if (!noForceFile_)
          {
            Screen::MasterInfo("Reading loading from file: " + fem_force_file_);
          }

        if (!noForceFile_)
          {            
            //read in nodeForce and nodeForceID
            force >> numForcedNodes_;
            Screen::MasterInfo("Number of Forced Nodes from File: "+to_string(numForcedNodes_));

            for (int forcedNode=0;forcedNode<numForcedNodes_;++forcedNode)
              {
                //read forced node ID
                int node;
                force >> node;
                nodeForceID.push_back(node);

                //read force DIR
                int dir;
                force >> dir;

                //allocate nodeForce
                int ndof_loc = rootNdofPerNode_[node];
                for (int d=0;d<ndof_loc;++d) nodeForce.push_back(0.);

                //read force magnitude
                double mag = 0.;
                force >> mag;            
                nodeForce[forcedNode*ndof_loc+dir] = mag;

                cout << "    Force Direction (x=0,y=1,z=2) = " << dir << endl;
                cout << "    Node Force = " << nodeForce[forcedNode*ndof_loc+dir] << endl;
                cout << "    Node Force ID = " << nodeForceID[forcedNode] << endl;

              }

            numMomentNodes_ = 0;
            force >> numMomentNodes_;
            Screen::MasterInfo("Number of Moment Nodes: "+to_string(numMomentNodes_-1));

            for (int momentNode=0;momentNode<numMomentNodes_;++momentNode) 
              {
                //read moment node ID
                int node;
                force >> node;
                nodeMomentID.push_back(node);

                //read moment DIR
                int dir;
                force >> dir;

                //allocate nodeMoment
                int ndof_loc = rootNdofPerNode_[node];
                for (int d=0;d<ndof_loc;++d) nodeMoment.push_back(0.);

                //read moment magnitude
                double mag = 0.;
                force >> mag;            
                nodeMoment[momentNode*ndof_loc+dir] = mag;

                cout << "    Moment Direction (x=0,y=1,z=2) = " << dir << endl;
                cout << "    Node Moment = " << nodeMoment[momentNode*ndof_loc+dir] << endl;
                cout << "    Node Moment ID = " << nodeMomentID[momentNode] << endl;

              }
          }
      }

    ierr = MPI_Bcast(&noForceFile_,1,MPI::BOOL,0,MPI_COMM_WORLD);

    if (!noForceFile_)
      {
        int fsize1;
        int fsize2;

        int msize1;
        int msize2;

        if (mypeno_ == 0)
          {
            fsize1 = numForcedNodes_;
            fsize2 = nodeForce.size();

            msize1 = numMomentNodes_;
            msize2 = nodeMoment.size();
          }

        ierr = MPI_Bcast(&fsize1,1,MPI::INTEGER,0,MPI_COMM_WORLD);
        ierr = MPI_Bcast(&msize1,1,MPI::INTEGER,0,MPI_COMM_WORLD);

        ierr = MPI_Bcast(&fsize2,1,MPI::INTEGER,0,MPI_COMM_WORLD);
        ierr = MPI_Bcast(&msize2,1,MPI::INTEGER,0,MPI_COMM_WORLD);

        if (mypeno_ != 0)
          {
            nodeForce   .resize(fsize2,0.);
            nodeForceID .resize(fsize1,-1);

            nodeMoment  .resize(msize2,0.);
            nodeMomentID.resize(msize1,-1);
          }      
        //broadcast the loading and forced nodes to all procs
        ierr = MPI_Bcast(&nodeForce   [0],fsize2,MPI::DOUBLE ,0,MPI_COMM_WORLD);
        ierr = MPI_Bcast(&nodeForceID [0],fsize1,MPI::INTEGER,0,MPI_COMM_WORLD);
        ierr = MPI_Bcast(&nodeMoment  [0],msize2,MPI::DOUBLE ,0,MPI_COMM_WORLD);
        ierr = MPI_Bcast(&nodeMomentID[0],msize1,MPI::INTEGER,0,MPI_COMM_WORLD);

        //to do this, temporarily broadcast rootNdofPerNode to get global context
        if (mypeno_ != 0) rootNdofPerNode_.resize(nof_nodes_glob_,-1);
        ierr = MPI_Bcast(&rootNdofPerNode_[0],nof_nodes_glob_,MPI::INT,0,MPI_COMM_WORLD);
        //apply forces loadvector
        const int& nof_nodeForces  = nodeForceID .size();
        const int& nof_nodeMoments = nodeMomentID.size();            

        int globalIndex = 0;
        for (int i=0;i<nof_nodeForces;++i)
          {
            //global forced node ID
            const int& node = nodeForceID[i];
            //ndof for global forced node ID
            int ndof_glob   = rootNdofPerNode_[node];

            int localIndex = 0;
            for (int localNode=0;localNode<myNodeCount_;++localNode)
              {
                //local forced node ID is localNode
                //ndof for local forced node ID
                int globalID = elemMaps_.nodeLocal2Global[localNode];
                int ndof_loc = rootNdofPerNode_[globalID];

                for (int idir=0;idir<DIM;++idir)
                  {
                    //check if local nodes are also forced nodes
                    if (elemMaps_.nodeLocal2Global[localNode] == node)
                      {
                        forceVec_[idir + localIndex] = nodeForce[idir + globalIndex];
                      }
                  }
                //march along index
                localIndex += ndof_loc;
              }
            //march along index
            globalIndex += ndof_glob;

          }

        globalIndex = 0;
        for (int i=0;i<nof_nodeMoments;++i)
          {
            //global moment node ID
            const int& node = nodeMomentID[i];
            //ndof for global moment node ID
            int ndof_glob   = rootNdofPerNode_[node];

            int localIndex = 0;
            for (int localNode=0;localNode<myNodeCount_;++localNode)
              {
                //local forced node ID is localNode
                //ndof for local forced node ID
                int globalID = elemMaps_.nodeLocal2Global[localNode];
                int ndof_loc = rootNdofPerNode_[globalID];

                for (int idir=0;idir<2;++idir)//only supports two moments
                  {
                    //check if local nodes are also moment nodes
                    if (elemMaps_.nodeLocal2Global[localNode] == node)
                      {
                        forceVec_[DIM + idir + localIndex] = nodeMoment[idir + globalIndex];
                      }
                  }
                //march along index
                localIndex += ndof_loc;                
              }
            //march along index
            globalIndex += ndof_glob;            
          }
      }
    //memory
    if (mypeno_ != 0) rootNdofPerNode_ = vector<int>();

#else

    //for serial FEM in parallel simualtion

    vector<int> nodeForceID;
    int forceDIR;
    vector<double> nodeForce;

    vector<int> nodeMomentID;
    int momentDIR;
    vector<double> nodeMoment;

    //read loading from file
    ifstream force;
    force.open(fem_force_file_);
    noForceFile_ = force.fail();
    if (noForceFile_ && !pressureLoading_ && !bodyForces_ && !isFSI_)
      {
        Screen::MasterError("File not found for structural loading conditions, and no other loading options have been specified");
      }
    else if (noForceFile_ && (pressureLoading_||bodyForces_||isFSI_))
      {
        Screen::MasterInfo("No structural loading from file detected, but loading is provided from another source:");
        cout << "    isFSI_           = " << isFSI_ << endl;
        cout << "    pressureLoading_ = " << pressureLoading_ << endl;
        cout << "    bodyForces_      = " << bodyForces_ << endl;
      }
    else if (!noForceFile_)
      {
        Screen::MasterInfo("Reading loading from file: " + fem_force_file_);
      }

    if (!noForceFile_)
      {
        //read in nodeForce and nodeForceID
        force >> numForcedNodes_;
        Screen::MasterInfo("Number of Forced Nodes from File: "+to_string(numForcedNodes_));

        for (int forcedNode=0;forcedNode<numForcedNodes_;++forcedNode)
          {
            //read forced node ID
            int node;
            force >> node;
            nodeForceID.push_back(node);

            //read force DIR
            force >> forceDIR;

            //allocate nodeForce
            int ndof_loc = rootNdofPerNode_[node];
            for (int d=0;d<ndof_loc;++d) nodeForce.push_back(0.);

            //read force magnitude
            double mag = 0.;
            force >> mag;            
            nodeForce[forcedNode*ndof_loc+forceDIR] = mag;

            cout << "    Force Direction (x=0,y=1,z=2) = " << forceDIR << endl;
            cout << "    Node Force = " << nodeForce[forcedNode*ndof_loc+forceDIR] << endl;
            cout << "    Node Force ID = " << nodeForceID[forcedNode] << endl;

          }

        numMomentNodes_ = 0;
        force >> numMomentNodes_;
        Screen::MasterInfo("Number of Moment Nodes: "+to_string(numMomentNodes_));

        for (int momentNode=0;momentNode<numMomentNodes_;++momentNode) 
          {
            //read moment node ID
            int node;
            force >> node;
            nodeMomentID.push_back(node);

            //read moment DIR
            force >> momentDIR;

            //allocate nodeMoment
            int ndof_loc = rootNdofPerNode_[node];
            for (int d=0;d<ndof_loc;++d) nodeMoment.push_back(0.);

            //read moment magnitude
            double mag = 0.;
            force >> mag;            
            nodeMoment[momentNode*ndof_loc+momentDIR] = mag;

            cout << "    Moment Direction (x=0,y=1,z=2) = " << momentDIR << endl;
            cout << "    Node Moment = " << nodeMoment[momentNode*ndof_loc+momentDIR] << endl;
            cout << "    Node Moment ID = " << nodeMomentID[momentNode] << endl;

          }
        //apply forces on loadvector
        const int& nof_nodeForces  = nodeForceID .size();
        const int& nof_nodeMoments = nodeMomentID.size();

        int globalIndex = 0;
        for (int i=0;i<nof_nodeForces;++i) 
          {
            //global forced node ID
            const int& node = nodeForceID[i];
            //ndof for global forced node ID
            int ndof_glob   = rootNdofPerNode_[node];

            int localIndex = 0;
            for (int localNode=0;localNode<myNodeCount_;++localNode)
              {
                //local forced node ID is localNode
                //ndof for local forced node ID
                int globalID = elemMaps_.nodeLocal2Global[localNode];
                int ndof_loc = rootNdofPerNode_[globalID];

                for (int idir=0;idir<DIM;++idir)
                  {
                    //check if local nodes are also forced nodes
                    if (elemMaps_.nodeLocal2Global[localNode] == node)
                      {
                        forceVec_[idir + localIndex] = nodeForce[idir + globalIndex];
                      }
                  }
                //march along index
                localIndex += ndof_loc;
              }
            //march along index 
            globalIndex += ndof_glob;

          }

        globalIndex = 0;
        for (int i=0;i<nof_nodeMoments;++i)
          {
            //global moment node ID
            const int& node = nodeMomentID[i];
            //ndof for global moment node ID 
            int ndof_glob   = rootNdofPerNode_[node];

            int localIndex = 0;
            for (int localNode=0;localNode<myNodeCount_;++localNode)
                  {
                    //local forced node ID is localNode
                    //ndof for local forced node ID 
                    int globalID = elemMaps_.nodeLocal2Global[localNode];
                    int ndof_loc = rootNdofPerNode_[globalID];

                    for (int idir=0;idir<2;++idir)//only supports 2 moments
                      {
                        //check if local nodes are also moment nodes
                        if (elemMaps_.nodeLocal2Global[localNode] == node)
                          {
                            forceVec_[DIM + idir + localIndex] = nodeMoment[idir + globalIndex];
                          }
                      }
                    //march along index
                    localIndex += ndof_loc;
                  }
            //march along index
            globalIndex += ndof_glob;
          }
      }
#endif

#endif
  }

  void fem::Setup_fem()
  {
#if (USEPETSC>0)
    //set num xyz dirs to 3
    xyz = 3;
    //---- FIX: DIM was declared in fem.h but never assigned anywhere ----
    //DIM is the spatial dimension of the structural problem. The MITC3 shell,
    //beam and cable elements are all embedded in 3D space, so DIM = 3.
    //(3Das2D handles the pseudo-2D case separately, via dimLoc_ in ini3Das2D.)
    DIM = 3;
    //---- FIX: these three members were never initialized anywhere ----
    //number of distinct element types present (from the ElementType option vector)
    numberOfElementTypes_ = int(elementTypes_.size());
    if (numberOfElementTypes_ < 1) Screen::MasterError("No element types specified");
    //the driver pushes the parsed in/out-of-plane moduli into element 0 of these vectors.
    //cache them as scalars so getElementProps() can broadcast them to every element.
    inPlaneElastic_  = inPlaneEmod_ .empty() ? defaultEModulus_ : inPlaneEmod_ [0];
    outPlaneElastic_ = outPlaneEmod_.empty() ? defaultEModulus_ : outPlaneEmod_[0];
    //-----------------------------------------------------------------
    noForceFile_ = true;
    surfaceTriangulation surfTri;
    //root proc creates directories
    if (mypeno_ == 0)
      {
        mkdir("post"         ,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir("post/aux"     ,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir("post/fem"     ,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir("post/timestep",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir("post/nodeData",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir("post/forces"  ,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir("post/loadstep",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir("restart_fem"  ,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir("debugFEM"     ,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      }
    //root processor reads triangulation --- keeps track of element types read in and ndof/node
    vector<int> elementTypesList_loc;
    vector<int> rootComponentIDs;
    if (mypeno_ == 0) surfTri.smartTriangulationRead(rootTotalNodes0_,rootTotalFaces_,surfaceFilename_,nof_nodes_glob_,nof_faces_glob_,elementTypes_,elementTypesList_loc,rootNdofPerNode_,rootComponentIDs,revertNormalsByCompID_);

    //initialize the current node structure on the root
    if (mypeno_ == 0 && !do_restart_)
      {
        Screen::MasterInfo("Initializing root structure via copy");
        rootTotalNodes_.resize(nof_nodes_glob_*xyz,-99.);
        rootTotalNodes_.assign(rootTotalNodes0_.begin(),rootTotalNodes0_.end());
      }
    else if (mypeno_ == 0 && do_restart_)
      {
        Screen::MasterInfo("Initializing root structure via restart");
        rootTotalNodes_.resize(nof_nodes_glob_*xyz,-99.);
        //get restart filename from the restart timestep
        string restart_filename_root;
        ostringstream OSS_root;
        OSS_root << "restart_fem/restart_root_nt_" << setw(8) << setfill('0') << restartTimestep_ << ".dat";
        restart_filename_root = OSS_root.str();

        ifstream file;
        file.open(restart_filename_root);
        file.precision(12);
        file.setf(ios::fixed);
        file.setf(ios::showpoint);

        for (int n=0;n<nof_nodes_glob_;++n)
          {
            for (int dir=0;dir<xyz;++dir)
              {
                file >> rootTotalNodes_[n*xyz + dir];
              }
          }
      }

    //                                DOMAIN PARTITION                                  //
    // This partitioning will distribute the total triangulation across all procs.      //
    // All processors, except the root, will only own a portion of the triangulation.   //
    // The root owns the entire triangulation for easy serial IO without communication. //
    // FSI is essentially N FSI problems where the geometry and FEM are partitioned by  //
    // N procs. Memory should be minimized where possible when writing code!            //

#if (DEDICATED_FEM_PROC == 0)
    PetscInitialize(NULL,NULL,NULL,NULL);
    //broadcast global node and element count
    ierr = MPI_Bcast(&nof_nodes_glob_,1,MPI::INTEGER,0,MPI_COMM_WORLD);
    ierr = MPI_Bcast(&nof_faces_glob_,1,MPI::INTEGER,0,MPI_COMM_WORLD);
    //manually initialize simple domain partition for ParMETIS
    initializeDomainPartition(elementTypesList_loc);
    //distribute elements from tmp_faces to each proc
    distributeElements();
    //use ParMETIS to re-partition the domain across procs
    //applyBC called inside
    //setForce called inside
    //setupElements called inside
    //identifyHaloElements called inside
    //initializeNormalVectors called inside
    partitionDomain_parmetis(elementTypesList_loc,rootComponentIDs);
#else
    //treat parallel simulations with serial FEM
    //setupElements called inside
    treatSingleProcFEM(elementTypesList_loc,rootComponentIDs);
#endif

    //initialize 3Das2D - introduce modified node/element counters
    ini3Das2D();

    //allocate global variables for each proc
    define();

    //initialize petsc
    initializePetsc();

    //prepare sim for restart
    if (do_restart_) readRestart();

    //set initial director vector directions
    initializeDirVecs();

    //tag nodes with special properties/BCs
    vector<int> markedFaces;
    markNodesInBox(markedFaces);

    //fill element propety vectors
    getElementProps(markedFaces);

    //compute load and displacement transfer maps
    if (useAuxGeometry_ || isFSI_) auxillaryTransfer_ini();

    Screen::MasterInfo("Prepping IO");

    //prep node IO
    if (recordNodeOutput_ && mypeno_ == 0)
      {
        nodeOutFile.resize(outputNodes_.size());
        for (int n=0;n<outputNodes_.size();++n)
          {
            nodeOutFile[n] = "post/nodeData/node_"+to_string(outputNodes_[n])+".dat";
            node_output[n].open(nodeOutFile[n],fstream::in|fstream::out|fstream::app);

	    //if file does not exist, create it
	    if (!node_output[n])
	      {
		node_output[n].open(nodeOutFile[n],fstream::in|fstream::out|fstream::trunc);
	      }
            node_output[n].precision(10);
          }
      }
    if (mypeno_ == 0) cout << "----------------------------" << endl;

    nIOCables_ = 0;
    //prep cable IO
    if (recordCableOutput_)
      {
	cablesForcesFile.resize(outputCables_.size());
	for (int f=0;f<myElementCount_;++f)
	  {
	    int globID = elemMaps_.elemLocal2Global[f];

	    if (find(outputCables_.begin(),outputCables_.end(),globID) != outputCables_.end())
	      {
		cablesForcesFile[nIOCables_] = "post/nodeData/cable_"+to_string(globID)+".dat";
		cable_output[nIOCables_].open(cablesForcesFile[nIOCables_],fstream::in|fstream::out|fstream::app);

		//if file does not exist, create it
		if (!cable_output[nIOCables_])
		  {
		    cable_output[nIOCables_].open(cablesForcesFile[nIOCables_],fstream::in|fstream::out|fstream::trunc);
		  }
		cable_output[nIOCables_].precision(10);
		myIOCables_[nIOCables_] = f;
		nIOCables_ += 1;
	      }
	  }
      }

    //initialize time integration vars
    if (!is_unsteady_)
      {
        numTimesteps_   = 0;
      }
    else if (!is_nonlinear_ and is_unsteady_)
      {
        assembleGlobKMat();
      }

    if (!do_restart_) restartTimestep_ = 0;
    itCount_ = 1;
    if (!is_unsteady_ && is_nonlinear_ && mypeno_ == 0) rhsInf.open("debugFEM/rhs_norm.dat");
    if (!is_nonlinear_ && is_unsteady_) ini_uSolve();

  }
  
  void fem::markNodesInBox(vector<int>& markedFaces)
  {

    //mark elements with a node in a box to apply special properties or BCs later

#if (MARK_NODES_IN_BOX == 1)

    Screen::MasterInfo("Marking nodes in box");    
    // --- for DGB parachute ---

    // for disk with 24 gores spaced 15degs apart starting at -7.5 CCW deg offset from y-axis in y-z plane

    vector<double> box(6,0.);
    box[0] = -1.E-6;
    box[1] =  1.E-6;
    box[2] =  0.;
    box[3] =  0.3;
    box[4] = -1.E-7;
    box[5] =  1.E-7;

    int counter = 0;
    double pi = 3.1415926535;    
    double theta = 7.5*pi/180.;

    vector<vector<double> > rotMat(3,vector<double>(3,0.));
    vector<double> xyzVec (3,0.);
    vector<double> xyz_rot(3,0.);

    //disk
    for (int gore=0;gore<24;++gore)
      {

        rotMat[0][0] =  1.;
        rotMat[1][1] =  cos(theta);
        rotMat[1][2] = -sin(theta);
        rotMat[2][1] =  sin(theta);
        rotMat[2][2] =  cos(theta);

        for (int e=0;e<myElementCount_;++e)
          {
            if (elementData_[e].type == MITC3)
              {
                int nnodes_loc    = elementData_[e].nnodes;
                int nnodes_offset = nnodesLocalOffset_[e];

                for (int n=0;n<nnodes_loc;++n)
                  {
                    bool inBox = true;
                    int node = faces_[n + nnodes_offset];

                    xyzVec[0] = nodes0_[node*xyz + 0];
                    xyzVec[1] = nodes0_[node*xyz + 1];
                    xyzVec[2] = nodes0_[node*xyz + 2];

                    mv_mult(rotMat,xyzVec,xyz_rot);

                    double x = xyz_rot[0];
                    double y = xyz_rot[1];
                    double z = xyz_rot[2];

                    if (x < box[0])
                      {
                        inBox = false;
                      }
                    else if (x > box[1])
                      {
                        inBox = false;
                      }
                    else if (y < box[2])
                      {
                        inBox = false;
                      }
                    else if (y > box[3])
                      {
                        inBox = false;
                      }
                    else if (z < box[4])
                      {
                        inBox = false;
                      }
                    else if (z > box[5])
                      {
                        inBox = false;
                      }

                    if (inBox)
                      {
                        markedFaces.push_back(e);
                        break;
                      }

                  }
              }
          }
        theta += 15.*pi/180.;
      }

    // for band with 24 gores etc...

    box[0] = -0.15;
    box[1] = -0.025;
    box[2] =  0.28;
    box[3] =  0.30;
    box[4] = -1.E-6;
    box[5] =  1.E-6;

    theta = 7.5*pi/180.;

    for (int gore=0;gore<24;++gore)
      {

        rotMat[0][0] =  1.;
        rotMat[1][1] =  cos(theta);
        rotMat[1][2] = -sin(theta);
        rotMat[2][1] =  sin(theta);
        rotMat[2][2] =  cos(theta);

        for (int e=0;e<myElementCount_;++e)
          {
            if (elementData_[e].type == MITC3)
              {
                int nnodes_loc    = elementData_[e].nnodes;
                int nnodes_offset = nnodesLocalOffset_[e];

                for (int n=0;n<nnodes_loc;++n)
                  {
                    bool inBox = true;
                    int node = faces_[n + nnodes_offset];

                    xyzVec[0] = nodes0_[node*xyz + 0];
                    xyzVec[1] = nodes0_[node*xyz + 1];
                    xyzVec[2] = nodes0_[node*xyz + 2];

                    mv_mult(rotMat,xyzVec,xyz_rot);

                    double x = xyz_rot[0];
                    double y = xyz_rot[1];
                    double z = xyz_rot[2];

                    if (x < box[0])
                      {
                        inBox = false;
                      }
                    else if (x > box[1])
                      {
                        inBox = false;
                      }
                    else if (y < box[2])
                      {
                        inBox = false;
                      }
                    else if (y > box[3])
                      {
                        inBox = false;
                      }
                    else if (z < box[4])
                      {
                        inBox = false;
                      }
                    else if (z > box[5])
                      {
                        inBox = false;
                      }

                    if (inBox)
                      {
                        markedFaces.push_back(e);
                        break;
                      }

                  }
              }
          }
        theta += 15.*pi/180.;
      }

    //thicken top and bottom of band and inner and outer disk radius

    for (int e=0;e<myElementCount_;++e)
      {
        if (elementData_[e].type == MITC3)
          {
            int nnodes_loc    = elementData_[e].nnodes;
            int nnodes_offset = nnodesLocalOffset_[e];

            for (int n=0;n<nnodes_loc;++n)
              {
                bool inBox = true;
                int node = faces_[n + nnodes_offset];

                xyzVec[0] = nodes0_[node*xyz + 0];
                xyzVec[1] = nodes0_[node*xyz + 1];
                xyzVec[2] = nodes0_[node*xyz + 2];

                inBox = false;
                if (xyzVec[0] < -0.13)
                  {
                    if ( find(markedFaces.begin(),markedFaces.end(),e) == markedFaces.end() )
                      {
                        markedFaces.push_back(e);
                        inBox = true;
                      }
                  }

                if ( (xyzVec[0] < -0.034036) && (xyzVec[0] > -0.0377236) )
                  {
                    if ( find(markedFaces.begin(),markedFaces.end(),e) == markedFaces.end() )
                      {
                        markedFaces.push_back(e);
                        inBox = true;
                      }
                  }

                if ( (pow( (xyzVec[1]*xyzVec[1]) + (xyzVec[2]*xyzVec[2]) , 0.5) > (0.29897578314 - 5.E-3)) && (xyzVec[0] > -0.0001) )
                  {
                    if ( find(markedFaces.begin(),markedFaces.end(),e) == markedFaces.end() )
                      {
                        markedFaces.push_back(e);
                        inBox = true;
                      }
                  }

                if ( (pow( (xyzVec[1]*xyzVec[1]) + (xyzVec[2]*xyzVec[2]) , 0.5) < (0.02869348316 + 5.E-3)) && (xyzVec[0] > -0.0001) )
                  {
                    if ( find(markedFaces.begin(),markedFaces.end(),e) == markedFaces.end() )
                      {
                        markedFaces.push_back(e);
                        inBox = true;
                      }
                  }

                if (inBox) break;

              }
          }
      }

    //string filename = "debugFEM/markFaces_"+to_string(mypeno_)+".vtk";
    //ofstream file;
    //file.open(filename);
    //file.precision(10);
    //file.setf(ios::fixed);
    //file.setf(ios::showpoint);
    //file << "# vtk DataFile Version 3.0" << endl;
    //file << "vtk output" << endl;;
    //file << "ASCII" << endl;
    //file << "DATASET POLYDATA" << endl;
    //file << "POINTS " << markedFaces.size() << " double" << endl;
    //for (int e=0;e<markedFaces.size();++e)
    //  {
    //    int element       = markedFaces[e];
    //    int nnodes_loc    = elementData_[element].nnodes;
    //    int nnodes_offset = nnodesLocalOffset_[element];
    //    double x = 0;
    //    double y = 0;
    //    double z = 0;
    //    for (int n=0;n<nnodes_loc;++n)
    //      {
    //        int node = faces_[n + nnodes_offset];
    //        x += 0.333333333 * nodes_[node*xyz + 0];
    //        y += 0.333333333 * nodes_[node*xyz + 1];
    //        z += 0.333333333 * nodes_[node*xyz + 2];
    //      }
    //    file << x << " " << y << " " << z << endl;
    //  }
    //file << "POINT_DATA " << markedFaces.size() << endl;
    //file << "SCALARS Data double" << endl;
    //file << "LOOKUP_TABLE default" << endl;
    //for (int e=0;e<markedFaces.size();++e)
    //  {
    //    file << 1.000000 << endl;
    //  }
    //file.close();

//#elif (MARK_NODES_IN_BOX == 2)
//    //for supersonic flag leading edge
//    Screen::MasterInfo("Marking nodes in box");
//    vector<double> box(6,0.);
//    box[0] = -1.;
//    box[1] = -0.4;
//    box[2] = -1.;
//    box[3] = 1.;
//    box[4] = -0.1;
//    box[5] = 0.1;
//    for (int e=0;e<myElementCount_;++e)
//      {
//        if (elementData_[e].type == MITC3)
//          {
//            int nnodes_loc    = elementData_[e].nnodes;
//            int nnodes_offset = nnodesLocalOffset_[e];
//            for (int n=0;n<nnodes_loc;++n)
//              {
//                bool inBox = true;
//                int node = faces_[n + nnodes_offset];
//                double x = nodes_[node*xyz+0];
//                double y = nodes_[node*xyz+1];
//                double z = nodes_[node*xyz+2];
//                if (x < box[0])
//                  {
//                    inBox = false;
//                  }
//                else if (x > box[1])
//                  {
//                    inBox = false;
//                  }
//                else if (y < box[2])
//                  {
//                    inBox = false;
//                  }
//                else if (y > box[3])
//                  {
//                    inBox = false;
//                  }
//                else if (z < box[4])
//                  {
//                    inBox = false;
//                  }
//                else if (z > box[5])
//                  {
//                    inBox = false;
//                  }
//                if (inBox)
//                  {
//                    markedFaces.push_back(e);
//                  }
//              }
//          }
//      }
#elif (MARK_NODES_IN_BOX == 3)

    Screen::MasterInfo("Marking nodes in box - 80-gore parachute");
    // --- for DGB parachute ---

    // for disk with 80 gores spaced 4.5degs apart starting at -2.25 CCW deg offset from y-axis in x-y plane

    vector<double> box(6,0.);
    box[0] = -1.E-4;
    box[1] =  1.E-4;
    box[2] =  0.;
    box[3] =  7.2;
    box[4] = -1.E-6;
    box[5] =  1.E-6;

    int counter = 0;
    double pi = 3.1415926535;    
    double theta = 2.25*pi/180.;

    vector<vector<double> > rotMat(3,vector<double>(3,0.));
    vector<double> xyzVec (3,0.);
    vector<double> xyz_rot(3,0.);

    //disk
    for (int gore=0;gore<80;++gore)
      {

        rotMat[0][0] =  cos(theta);
        rotMat[0][1] = -sin(theta); 
        rotMat[1][0] =  sin(theta);
        rotMat[1][1] =  cos(theta);
        rotMat[2][2] =  1.;

        for (int e=0;e<myElementCount_;++e)
          {
            if (elementData_[e].type == MITC3)
              {
                int nnodes_loc    = elementData_[e].nnodes;
                int nnodes_offset = nnodesLocalOffset_[e];

                for (int n=0;n<nnodes_loc;++n)
                  {
                    bool inBox = true;
                    int node = faces_[n + nnodes_offset];

                    xyzVec[0] = nodes_[node*xyz + 0];
                    xyzVec[1] = nodes_[node*xyz + 1];
                    xyzVec[2] = nodes_[node*xyz + 2];

                    mv_mult(rotMat,xyzVec,xyz_rot);

                    double x = xyz_rot[0];
                    double y = xyz_rot[1];
                    double z = xyz_rot[2];

                    if (x < box[0])
                      {
                        inBox = false;
                      }
                    else if (x > box[1])
                      {
                        inBox = false;
                      }
                    else if (y < box[2])
                      {
                        inBox = false;
                      }
                    else if (y > box[3])
                      {
                        inBox = false;
                      }
                    else if (z < box[4])
                      {
                        inBox = false;
                      }
                    else if (z > box[5])
                      {
                        inBox = false;
                      }

                    if (inBox)
                      {
                        markedFaces.push_back(e);
                        break;
                      }

                  }
              }
          }
        theta += 4.5*pi/180.;
      }

    // for band with 80 gores etc...

    box[0] = -1.E-4;
    box[1] =  1.E-4;
    box[2] =  7.18;
    box[3] =  7.2;
    box[4] = -0.904-2.580 - 1.e-4;
    box[5] = -0.904 + 1.e-4;

    theta = 2.25*pi/180.;

    for (int gore=0;gore<80;++gore)
      {

        rotMat[0][0] =  cos(theta);
        rotMat[0][1] = -sin(theta); 
        rotMat[1][0] =  sin(theta);
        rotMat[1][1] =  cos(theta);
        rotMat[2][2] =  1.;

        for (int e=0;e<myElementCount_;++e)
          {
            if (elementData_[e].type == MITC3)
              {
                int nnodes_loc    = elementData_[e].nnodes;
                int nnodes_offset = nnodesLocalOffset_[e];

                for (int n=0;n<nnodes_loc;++n)
                  {
                    bool inBox = true;
                    int node = faces_[n + nnodes_offset];

                    xyzVec[0] = nodes_[node*xyz + 0];
                    xyzVec[1] = nodes_[node*xyz + 1];
                    xyzVec[2] = nodes_[node*xyz + 2];

                    mv_mult(rotMat,xyzVec,xyz_rot);

                    double x = xyz_rot[0];
                    double y = xyz_rot[1];
                    double z = xyz_rot[2];

                    if (x < box[0])
                      {
                        inBox = false;
                      }
                    else if (x > box[1])
                      {
                        inBox = false;
                      }
                    else if (y < box[2])
                      {
                        inBox = false;
                      }
                    else if (y > box[3])
                      {
                        inBox = false;
                      }
                    else if (z < box[4])
                      {
                        inBox = false;
                      }
                    else if (z > box[5])
                      {
                        inBox = false;
                      }

                    if (inBox)
                      {
                        markedFaces.push_back(e);
                        break;
                      }

                  }
              }
          }
        theta += 4.5*pi/180.;
      }

#endif

    //if (markedFaces.size() > 0)
    //  {
    //	string filename = "debugFEM/markFaces_"+to_string(mypeno_)+".vtk";
    //	ofstream file;
    //	file.open(filename);
    //	file.precision(10);
    //	file.setf(ios::fixed);
    //	file.setf(ios::showpoint);
    //	file << "# vtk DataFile Version 3.0" << endl;
    //	file << "vtk output" << endl;;
    //	file << "ASCII" << endl;
    //	file << "DATASET POLYDATA" << endl;
    //	file << "POINTS " << markedFaces.size() << " double" << endl;
    //	for (int e=0;e<markedFaces.size();++e)
    //	  {
    //	    int element       = markedFaces[e];
    //	    int nnodes_loc    = elementData_[element].nnodes;
    //	    int nnodes_offset = nnodesLocalOffset_[element];
    //	    double x = 0;
    //	    double y = 0;
    //	    double z = 0;
    //	    for (int n=0;n<nnodes_loc;++n)
    //	      {
    //		int node = faces_[n + nnodes_offset];
    //		x += 0.333333333 * nodes_[node*xyz + 0];
    //		y += 0.333333333 * nodes_[node*xyz + 1];
    //		z += 0.333333333 * nodes_[node*xyz + 2];
    //	      }
    //	    file << x+1.e-10 << " " << y+1.e-10 << " " << z+1.e-10 << endl;
    //	  }
    //	file << "POINT_DATA " << markedFaces.size() << endl;
    //	file << "SCALARS Data double" << endl;
    //	file << "LOOKUP_TABLE default" << endl;
    //	for (int e=0;e<markedFaces.size();++e)
    //	  {
    //	    file << 1.000000 << endl;
    //	  }
    //	file.close();
    //  }

#endif
  }
  
  void fem::setupElements(int& nonUniqueNodeCount,
                          vector<int>& globalNodeList)
  {
    //distribute respective parts of ndofPerNode_loc to all procs
#if (DEDICATED_FEM_PROC == 0)    

    vector<int> nonUnique_ndofPerNodeList;
    if (mypeno_ == 0)
      {
        //perform non-unique, proc-based ordering for scattering
        nonUnique_ndofPerNodeList.resize(nonUniqueNodeCount,-1);
        for (int n=0;n<nonUniqueNodeCount;++n)
          {
            int node = globalNodeList[n];
            nonUnique_ndofPerNodeList[n] = rootNdofPerNode_[node];
          }
      }

    //recompute offsets and displs with nodes
    vector<int> partitionNodeCounts(numprocs_,-1);
    MPI_Allgather(&myNodeCount_,1,MPI_INT,&partitionNodeCounts[0],1,MPI_INT,MPI_COMM_WORLD);

    vector<int> partitionNodeOffset(numprocs_,-1);
    int myNodeOffset = myNodeCount_*mypeno_;
    partitionNodeOffset[0] = 0;
    for (int p=1;p<numprocs_;++p) partitionNodeOffset[p] = partitionNodeOffset[p-1] + partitionNodeCounts[p-1];

    //recover ndof/node from proc-ordered list
    ndofPerNode_.resize(myNodeCount_,0);
    MPI_Scatterv(&nonUnique_ndofPerNodeList[0],&partitionNodeCounts[0],&partitionNodeOffset[0],MPI::INT,&ndofPerNode_[0],partitionNodeCounts[mypeno_],MPI::INT,0,MPI_COMM_WORLD);

    //get my dof global indices to avoid recomputing
    if (mypeno_ != 0) rootNdofPerNode_.resize(nof_nodes_glob_,-1);
    ierr = MPI_Bcast(&rootNdofPerNode_[0],nof_nodes_glob_,MPI::INT,0,MPI_COMM_WORLD);

    int globalOffset = 0;
    ndofGlobalOffset_.resize(myNodeCount_,0);
    for (int nglob=0;nglob<nof_nodes_glob_;++nglob)
      {
        for (int nloc=0;nloc<myNodeCount_;++nloc)
          {
            int l2gID = elemMaps_.nodeLocal2Global[nloc];

            if (l2gID == nglob)
              {
                ndofGlobalOffset_[nloc] = globalOffset;
              }
          }
        globalOffset += rootNdofPerNode_[nglob];
      }
    //memory
    if (mypeno_ != 0) rootNdofPerNode_ = vector<int>();

    //get my dof local indices to avoid recomputing
    ndofLocalOffset_.resize(myNodeCount_,0);
    int index = 0;
    for (int nloc=0;nloc<myNodeCount_;++nloc)
      {
        ndofLocalOffset_[nloc] = index;
        index += ndofPerNode_[nloc];
      }

    //get nnodes offset in local faces_
    nnodesLocalOffset_.resize(myElementCount_,-1);
    index = 0;
    for (int e=0;e<myElementCount_;++e)
      {
        nnodesLocalOffset_[e] = index;
        index += elementData_[e].nnodes;
      }

#else

    globalSystemSize_ = 0;
    for (int n=0;n<nof_nodes_glob_;++n) globalSystemSize_ += rootNdofPerNode_[n];
    if (globalSystemSize_ > INT_MAX) Screen::MasterError("Global system size is larger than the maximum size of an integer!, System size: "+to_string(globalSystemSize_));

    ndofPerNode_.resize(nof_nodes_glob_,0);
    for (int n=0;n<nonUniqueNodeCount;++n)
      {
        int node = globalNodeList[n];
        ndofPerNode_[node] = rootNdofPerNode_[node];
      }

    int offset = 0;
    ndofLocalOffset_.resize(nof_nodes_glob_,0);
    ndofGlobalOffset_.resize(nof_nodes_glob_,0);
    for (int nglob=0;nglob<nof_nodes_glob_;++nglob)
      {
        ndofLocalOffset_[nglob] = offset;
        ndofGlobalOffset_[nglob] = offset;
        offset += rootNdofPerNode_[nglob];
      }

#endif

  }
  
  void fem::initializeDomainPartition(vector<int>& elementTypesList)
  {
    //manually initialize domain partition for ParMETIS
    Screen::MasterInfo("Initializing domain partition");

    myElementCount_ = floor(double(nof_faces_glob_)/double(numprocs_));
    myOffset_       = mypeno_*myElementCount_;
    if (mypeno_ == numprocs_-1) myElementCount_ = myElementCount_ + nof_faces_glob_ - numprocs_*myElementCount_;

    partitionCellCounts_.resize(numprocs_,0);
    ierr = MPI_Allgather(&myElementCount_,1,MPI::INT,&partitionCellCounts_[0],1,MPI::INT,MPI_COMM_WORLD);
    if (ierr != 0) Screen::MasterError("MPI_Allgather - Ierr not equal zero fem::initializeDomainPartition, ierr = "+to_string(ierr));
    partitionOffset_.resize(numprocs_);
    myOffset_ = 0;
    partitionOffset_[0] = 0;
    for (int p=1;p<numprocs_;++p) partitionOffset_[p] = partitionOffset_[p-1] + partitionCellCounts_[p-1];
    myOffset_ = partitionOffset_[mypeno_];

    //scatter element types from root
    vector<int> myElementTypes(myElementCount_,-1);
    MPI_Scatterv(&elementTypesList[0],&partitionCellCounts_[0],&partitionOffset_[0],MPI::INT,&myElementTypes[0],partitionCellCounts_[mypeno_],MPI::INT,0,MPI_COMM_WORLD);
    if (ierr != 0) Screen::MasterError("MPI_Scatterv - Ierr not equal zero fem::initializeDomainPartition, ierr = "+to_string(ierr));

    //get rootNnodesPerEl_ before clear
    if (mypeno_ == 0)
      {
        rootNnodesPerEl_.resize(nof_faces_glob_,-1);
        for (int ng=0;ng<nof_faces_glob_;++ng)
          {
            if (elementTypesList[ng] == SIMPLETRI)
              {
                rootNnodesPerEl_[ng] = 3;
              }
            if (elementTypesList[ng] == LMITC3)
              {
                rootNnodesPerEl_[ng] = 3;
              }
            if (elementTypesList[ng] == MITC3)
              {
                rootNnodesPerEl_[ng] = 3;
              }
            if (elementTypesList[ng] == BEAM)
              {
                rootNnodesPerEl_[ng] = 2;
              }
            if (elementTypesList[ng] == CABLE)
              {
                rootNnodesPerEl_[ng] = 2;
              }
          }
      }

    myNumberOfEachElementType_.resize(numberOfElementTypes_,0);
    //get element global IDs and types into element struct
    element elem;
    elementData_.resize(myElementCount_,elem);
    for (int e=0;e<myElementCount_;++e)
      {
        elementData_[e].globID = e + partitionOffset_[mypeno_];
        elementData_[e].type   = myElementTypes[e];

        //count number of each element type
        for (int ne=0;ne<numberOfElementTypes_;++ne)
          {
            if (elementData_[e].type == elementTypes_[ne])
              {
                myNumberOfEachElementType_[ne] += 1;
              }
          }
      }

  }
  
  void fem::distributeElements()
  {
    //distribute connectivity from root proc to all procs based off dumb partition
    Screen::MasterInfo("Scattering connectivity to all procs from root");

    vector<int> partitionCellCounts_with_nnodes(numprocs_,0);
    vector<int> partitionOffset_with_nnodes    (numprocs_,0);

    int sendCount = 0;
    for (int e=0;e<myElementCount_;++e)
      {
        if (elementData_[e].type == SIMPLETRI)
          {
            sendCount += 3;
          }
        if (elementData_[e].type == LMITC3)
          {
            sendCount += 3;
          }
        if (elementData_[e].type == MITC3)
          {
            sendCount += 3;
          }
        if (elementData_[e].type == BEAM)
          {
            sendCount += 2;
          }
        if (elementData_[e].type == CABLE)
          {
            sendCount += 2;
          }
      }

    MPI_Allgather(&sendCount,1,MPI::INT,&partitionCellCounts_with_nnodes[0],1,MPI::INT,MPI_COMM_WORLD);
    partitionOffset_with_nnodes[0] = 0;
    for (int np=1;np<numprocs_;++np) partitionOffset_with_nnodes[np] = partitionOffset_with_nnodes[np-1] + partitionCellCounts_with_nnodes[np-1];

    faces_.resize(partitionCellCounts_with_nnodes[mypeno_]);
    MPI_Scatterv(&rootTotalFaces_[0],&partitionCellCounts_with_nnodes[0],&partitionOffset_with_nnodes[0],MPI::INT,&faces_[0],partitionCellCounts_with_nnodes[mypeno_],MPI::INT,0,MPI_COMM_WORLD);

    //sendbuf:   rootTotalFaces_ is total connectivity from root processor
    //sendcount: partitionCellCounts_with_nnodes is myElementCount_*nnodes_ to accomodate 1D face vector
    //displs:    partitionOffset_with_nnodes provides the usual ppartitionOffset multipied by nnodes, for 'jumping' to each procs starting point in tmp_faces
    //recvbuf:   faces_ is the receive buffer, sized to &partitionCellCounts_with_nnodes[mypeno_]
    //recvcount: partitionCellCounts_with_nnodes is the size of the data being received

    //get nnodes for element struct --- for dumb partition right now
    for (int e=0;e<myElementCount_;++e)
      {
        if (elementData_[e].type == SIMPLETRI)
          {
            elementData_[e].nnodes = 3;
            elementData_[e].ndof   = 6;
          }
        if (elementData_[e].type == LMITC3)
          {
            elementData_[e].nnodes = 3;
            elementData_[e].ndof   = 6;
          }
        if (elementData_[e].type == MITC3)
          {
            elementData_[e].nnodes = 3;
            elementData_[e].ndof   = 5;
          }
        if (elementData_[e].type == BEAM)
          {
            elementData_[e].nnodes = 2;
            elementData_[e].ndof   = 2;
          }
        if (elementData_[e].type == CABLE)
          {
            elementData_[e].nnodes = 2;
            elementData_[e].ndof   = 3;
          }
      }

  }
  
  void fem::partitionDomain_parmetis(vector<int>& elementTypesList_loc,
                                     vector<int>& rootComponentIDs)
  {
#if (USEPETSC>0)

    //Partition mesh using ParMETIS
    Screen::MasterInfo("Initializing ParMETIS");

    //ParMETIS_V3_PartMeshKway(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt, int *wgtflag, int *numflag, int *ncon, int * ncommonnodes, int *nparts, float *tpwgts, float *ubvec, int *option$

    //elmdist      : partition offsets. Sized                    (np+1)
    //eptr         : Element connectivity indices                (i.e., 0, nnodes, nnodes*2, etc.) - read ParMETIS manual Sec. 4.2.3
    //eind         : Element connectivity list                   (i.e., its just 'faces_' for each proc and its local elements - the face nodes)
    //elmwgt       : Element weights                             (null)
    //wgtflag      : Can take value of 0,1,2,3                   (use 0, no weights)
    //numflag      : 0 C-style numbers, 1 Fortran-style numbers  (0, C++)
    //ncon         : Number of weights each cell has             (use 1)
    //ncommonnodes : Min number of nodes shared between elements (2 for tris,  3 for quads)
    //nparts       : Number of partitions requested              (numprocs_)
    //tpwgts       : Function of weights for each part           (set all to 1/np, Size ncon x nparts.)
    //ubvec        : Balance tolerance, 1.05 is recommended      (5% imbalance allowed)
    //options      : [0 1 15] for default
    //edgecut      : OUTPUT, number of edges
    //part         : OUTPUT, sized to myElementCount_, integer value determines which proc to send element to
    //comm         : MPI communicator                            (MPI_COMM_WORLD usually)

    //transfer partitionOffset_ to ParMETIS var elemdist
    idx_t elmdist[numprocs_+1];
    for (int p=0;p<numprocs_;++p) elmdist[p] = partitionOffset_[p];
    elmdist[numprocs_] = nof_faces_glob_;

    //fill in eptr and eind
    //eptr: indices in faces
    idx_t *eptr;
    eptr = new idx_t[myElementCount_+1];
    eptr[0] = 0;
    int index = 0;
    for (int ec=0;ec<myElementCount_;++ec)
      {
        index += elementData_[ec].nnodes;
        eptr[ec+1] = index;
      }
    //eind: node IDs
    idx_t *eind;
    eind = new idx_t[eptr[myElementCount_]];
    for (int i=0;i<faces_.size();++i) eind[i] = faces_[i];

    //rest of ParMETIS inputs
    idx_t* ncommonnodes = new idx_t; //number of shared nodes to be considered an edge in the graph
    (*ncommonnodes)     = 1;
    idx_t* part         = new idx_t[myElementCount_]; //output
    idx_t* elmwgt       = NULL;
    idx_t* wgtflag      = new idx_t; //no weights associated with elem or edges
    (*wgtflag)          = 0;
    idx_t* numflag      = new idx_t; //C-style numbering
    (*numflag)          = 0;
    idx_t* ncon         = new idx_t; //# of weights or constraints
    (*ncon)             = 1;
    real_t tpwgts[numprocs_];
    for (int p=0;p<numprocs_;++p) tpwgts[p] = 1./float(numprocs_);
    float ubvec = 1.05; //Allow up to 5% imbalance
    idx_t options[3];   //default values for timing info set 0 -> 1
    options[0] = 0;
    options[1] = 1;
    options[2] = 15;
    idx_t*  edgecut      = new idx_t; //output
    idx_t*  numprocs_tmp = new idx_t;
    real_t* ubvec_tmp    = new real_t;
    (*ubvec_tmp)         = ubvec;
    (*numprocs_tmp)      = numprocs_;
    MPI_Comm commWorld   = MPI_COMM_WORLD;

    Screen::MasterInfo("Calling ParMETIS for mesh re-partitioning");
    ierr = ParMETIS_V3_PartMeshKway(elmdist,eptr,eind, elmwgt,
                                    wgtflag, numflag, ncon, ncommonnodes,
                                    numprocs_tmp, tpwgts, ubvec_tmp, options, edgecut,
                                    part,&commWorld);
    Screen::MasterInfo("ParMETIS complete - distributing new partition");    

    vector<int> nsend         (numprocs_,0);
    vector<int> nrecv         (numprocs_,0); 
    vector<vector<int> > send (numprocs_);//send element global IDs
    vector<vector<int> > send2(numprocs_);//send element types

    //'part' details where to send the elements the procs contain from the initial partition
    for (int e=0;e<myElementCount_;++e)
      {
        nsend[part[e]]++;
        send [part[e]].push_back(elementData_[e].globID);
        send2[part[e]].push_back(elementData_[e].type);
      }
    //tell all procs how many elements everyone else is sending and receiving
    MPI_Alltoall(&nsend[0],1,MPI_INT,&nrecv[0],1,MPI_INT,MPI_COMM_WORLD);

    //prepare for 'on the fly' element list consolidation
    faces_ = vector<int>();
    elementData_.resize(0);
    vector<int> localElementList;
    vector<int> globalElementList;
    vector<int> myProcList;
    vector<int> globalProcList;
    if (mypeno_ == 0) globalElementList.resize(nof_faces_glob_,-1);
    if (mypeno_ == 0) globalProcList.resize   (nof_faces_glob_,-1);

    Screen::MasterInfo("Creating <local,global> element mappings");
    //send and receive elements, exchange elements based off part
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    for (int p=0;p<numprocs_;++p)
      {

        vector<int> recv (nrecv[p],-1);
        vector<int> recv2(nrecv[p],-1);

        MPI_Sendrecv(&send [p][0],nsend[p],MPI_INT,p,0,&recv [0],nrecv[p],MPI_INT,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Sendrecv(&send2[p][0],nsend[p],MPI_INT,p,0,&recv2[0],nrecv[p],MPI_INT,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        for (int c=0;c<recv.size();++c)
          {
            element tmp;
            tmp.globID = recv [c];//receive glob ID
            tmp.type   = recv2[c];//recieve type
            elemMaps_.elemLocal2Global.insert(pair<int,int> (elementData_.size(),recv[c]));
            elementData_.push_back(tmp);//saving globID and type in elementData_ struct
          }
        //memory
        recv.clear();
      }
    //memory
    send.clear();
    //delete[] part;
    delete[] eptr;
    delete[] eind;

    //trim padding due to push_back
    vector<element> (elementData_).swap(elementData_);

    //count element number for each proc
    myElementCount_ = elementData_.size();
    //check if any proc did not recieve any elements
    if (myElementCount_ == 0) Screen::MasterError("Reduce the number of procs (or increase the number of elements) - FEM domain is too small");
    fill(partitionCellCounts_.begin(),partitionCellCounts_.end(),0);
    MPI_Allgather(&myElementCount_,1,MPI_INT,&partitionCellCounts_[0],1,MPI_INT,MPI_COMM_WORLD);
    //recompute offsets
    partitionOffset_.resize(numprocs_);
    myOffset_ = 0;
    partitionOffset_[0] = 0;
    for (int p=1;p<numprocs_;++p) partitionOffset_[p] = partitionOffset_[p-1] + partitionCellCounts_[p-1];
    partitionOffset_.push_back(nof_faces_glob_);
    myOffset_ = partitionOffset_[mypeno_];

    //get localElementList's from the <local,global> element mapping
    localElementList.resize(myElementCount_,-1);
    for (int local=0;local<myElementCount_;++local) localElementList[local] = elemMaps_.elemLocal2Global[local];  

    //gather localElementList's into globalElementList on root
    ierr = MPI_Gatherv(&localElementList[0],partitionCellCounts_[mypeno_],MPI::INT,&globalElementList[0],&partitionCellCounts_[0],&partitionOffset_[0],MPI::INT,0,MPI_COMM_WORLD);
    //memory
    localElementList = vector<int>();
    //get nnodes for element struct
    for (int e=0;e<myElementCount_;++e)
      {
        if (elementData_[e].type == SIMPLETRI)
          {
            elementData_[e].nnodes = 3;
            elementData_[e].ndof   = 6;
          }
        if (elementData_[e].type == LMITC3)
          {
            elementData_[e].nnodes = 3;
            elementData_[e].ndof   = 6;
          }
        if (elementData_[e].type == MITC3)
          {
            elementData_[e].nnodes = 3;
            elementData_[e].ndof   = 5;
          }
        if (elementData_[e].type == BEAM)
          {
            elementData_[e].nnodes = 2;
            elementData_[e].ndof   = 2;
          }
        if (elementData_[e].type == CABLE)
          {
            elementData_[e].nnodes = 2;
            elementData_[e].ndof   = 3;
          }
      }
    //get proc list for globalElement2ProcMap
    myProcList.resize(myElementCount_,mypeno_);
    ierr = MPI_Gatherv(&myProcList[0],partitionCellCounts_[mypeno_],MPI::INT,&globalProcList[0],&partitionCellCounts_[0],&partitionOffset_[0],MPI::INT,0,MPI_COMM_WORLD);
    //memory
    myProcList = vector<int>();

    vector<int> temporary_nnodes_offset;
    //fill globalElement2ProcMap and dump CPU partition
    vector<int> allFacesListCopy;
    if (mypeno_ == 0)
      {
        allFacesListCopy.resize(rootTotalFaces_.size(),-1);
        for (int e=0;e<nof_faces_glob_;++e)
          {
            elemMaps_.globalElement2ProcMap.insert(pair<int,int> (globalElementList[e],globalProcList[e]));
          }

        dumpPartition(rootTotalNodes0_,rootTotalFaces_);

        //get global system size
        globalSystemSize_ = 0;
        for (int n=0;n<nof_nodes_glob_;++n) globalSystemSize_ += rootNdofPerNode_[n];
	if (globalSystemSize_ > INT_MAX) Screen::MasterError("Global system size is larger than the maximum size of an integer!, System size: "+to_string(globalSystemSize_));

        //create temporary nnodes offset on root
        int offset = 0;
        temporary_nnodes_offset.resize(nof_faces_glob_,-1);
        for (int e=0;e<nof_faces_glob_;++e)
          {
            temporary_nnodes_offset[e] = offset;
            offset += rootNnodesPerEl_[e];
          }
        //use proc-based element ordering in globalElementList to re-order face nodes in a proc-order inside allFacesListCopy
        // --- then scatter allFacesListCopy to each processes' faces_ --- 
        int counter = 0;
        for (int e=0;e<nof_faces_glob_;++e)
          {
            int elem         = globalElementList[e];
            int nnodes_loc   = rootNnodesPerEl_ [elem];
            int nnodesoffset = temporary_nnodes_offset[elem];
            for (int n=0;n<nnodes_loc;++n)
              {
                allFacesListCopy[counter] = rootTotalFaces_[n + nnodesoffset];
                counter += 1;
              }
          }
      }
    //broadcast global system size
    ierr = MPI_Bcast(&globalSystemSize_,1,MPI::INTEGER,0,MPI_COMM_WORLD);

    //recover local faces_ information (face nodes) from re-ordered allFacesList
    Screen::MasterInfo("Getting proc local connectivity and nodal locations");
    vector<int> partitionCellCounts_with_nnodes(numprocs_,0);
    vector<int> partitionOffset_with_nnodes    (numprocs_,0);

    int sendCount = 0;
    for (int e=0;e<myElementCount_;++e)
      {        
        if (elementData_[e].type == SIMPLETRI)
          {
            sendCount += 3;
          }
        if (elementData_[e].type == LMITC3)
          {
            sendCount += 3;
          }
        if (elementData_[e].type == MITC3)
          {
            sendCount += 3;
          }
        if (elementData_[e].type == BEAM)
          {
            sendCount += 2;
          }
        if (elementData_[e].type == CABLE)
          {
            sendCount += 2;
          }
      }

    MPI_Allgather(&sendCount,1,MPI::INT,&partitionCellCounts_with_nnodes[0],1,MPI::INT,MPI_COMM_WORLD);
    partitionOffset_with_nnodes[0] = 0;
    for (int np=1;np<numprocs_;++np) partitionOffset_with_nnodes[np] = partitionOffset_with_nnodes[np-1] + partitionCellCounts_with_nnodes[np-1];

    faces_.resize(partitionCellCounts_with_nnodes[mypeno_],-1);
    MPI_Scatterv(&allFacesListCopy[0],&partitionCellCounts_with_nnodes[0],&partitionOffset_with_nnodes[0],MPI::INT,&faces_[0],partitionCellCounts_with_nnodes[mypeno_],MPI::INT,0,MPI_COMM_WORLD);
    allFacesListCopy = vector<int>();
    //memory
    globalProcList                  = vector<int>();
    globalElementList               = vector<int>();
    partitionCellCounts_with_nnodes = vector<int>();
    partitionOffset_with_nnodes     = vector<int>();

    //get a unique <local,global> node mapping
    vector<int> uniqueFaceNodes;
    myNodeCount_ = 0;
    int offset = 0;
    for (int e=0;e<myElementCount_;++e)
      {
        int nnodes_loc = elementData_[e].nnodes;
        for (int n=0;n<nnodes_loc;++n)
          {
            int index = n + offset;
            int node  = faces_[index];
            if (find(uniqueFaceNodes.begin(),uniqueFaceNodes.end(),node) == uniqueFaceNodes.end())
              {
                uniqueFaceNodes.push_back(node);
                elemMaps_.nodeLocal2Global.insert(pair<int,int>(myNodeCount_,node));
                myNodeCount_++;
              }
          }
        offset += nnodes_loc;
      }
    //memory
    uniqueFaceNodes = vector<int>();

    //take a moment to replace global node index in faces_ with local node index
    vector<int> facesCopy;
    facesCopy = faces_;
    for (int n=0;n<myNodeCount_;++n)
      {
        int globalNode = elemMaps_.nodeLocal2Global[n];

        for (int i=0;i<faces_.size();++i)
          {
            if (globalNode == faces_[i])
              {
                facesCopy[i] = n;
              }
          }
      }
    faces_.assign(facesCopy.begin(),facesCopy.end());
    facesCopy = vector<int>();

    //proc-order allNodesList for communication so all procs can get their own nodal locs from allNodesList

    vector<int> localNodeList;
    vector<int> globalNodeList;
    localNodeList.resize(myNodeCount_,-1.);

    //compute offset and displs
    vector<int> partitionNodeCounts(numprocs_,-1);
    vector<int> partitionNodeOffset(numprocs_,-1);
    int myNodeOffset = myNodeCount_*mypeno_;

    MPI_Allgather(&myNodeCount_,1,MPI_INT,&partitionNodeCounts[0],1,MPI_INT,MPI_COMM_WORLD);

    partitionNodeOffset[0] = 0;
    for (int p=1;p<numprocs_;++p) partitionNodeOffset[p] = partitionNodeOffset[p-1] + partitionNodeCounts[p-1];

    //compute number of total non-unique nodes (i.e., shared nodes between n procs are counted n times)
    int nonUniqueNodeCount = 0;
    for (int p=0;p<numprocs_;++p)
      {
        nonUniqueNodeCount += partitionNodeCounts[p];
      }
    if (mypeno_ == 0) globalNodeList.resize(nonUniqueNodeCount,-1);    

    //fill local node lists
    for (int local=0;local<myNodeCount_;++local) localNodeList[local] = elemMaps_.nodeLocal2Global[local];
    //gatherv nodes IDs in proc order
    ierr = MPI_Gatherv(&localNodeList[0],partitionNodeCounts[mypeno_],MPI::INT,&globalNodeList[0],&partitionNodeCounts[0],&partitionNodeOffset[0],MPI::INT,0,MPI_COMM_WORLD);
    //memory
    localNodeList = vector<int>();

    //use proc-based node ordering in globalNodeList to re-order allNodesList on root
    vector<double> allNodesListNonUnique;
    vector<double> allNodesListNonUnique0;
    if (mypeno_ == 0)
      {
        allNodesListNonUnique .resize(nonUniqueNodeCount*xyz,-99.);
	allNodesListNonUnique0.resize(nonUniqueNodeCount*xyz,-99.);
        for (int n=0;n<nonUniqueNodeCount;++n)
          {
            for (int idir=0;idir<xyz;++idir)
              {
                int globalIndex                          = idir + xyz*globalNodeList[n];
                int procOrderedIndex                     = idir + xyz*n;
                allNodesListNonUnique [procOrderedIndex] = rootTotalNodes_ [globalIndex];
                allNodesListNonUnique0[procOrderedIndex] = rootTotalNodes0_[globalIndex];
              }
          }
      }
    //re-compute offset and displs to change index from node ID to (x,y,z) coords
    for (int np=0;np<numprocs_;++np)
      {
        partitionNodeCounts[np] *= xyz;
        partitionNodeOffset[np] *= xyz;
      }
    //recover local node locs from re-ordered allNodesList
    nodes_.resize  (myNodeCount_*xyz,0.);
    nodes0_.resize (myNodeCount_*xyz,0.);
    MPI_Scatterv(&allNodesListNonUnique0[0],&partitionNodeCounts[0],&partitionNodeOffset[0],MPI::DOUBLE,&nodes0_[0],partitionNodeCounts[mypeno_],MPI::DOUBLE,0,MPI_COMM_WORLD);
    if (!do_restart_)
      {
	nodes_.assign(nodes0_.begin(),nodes0_.end());
      }
    else
      {
	MPI_Scatterv(&allNodesListNonUnique[0],&partitionNodeCounts[0],&partitionNodeOffset[0],MPI::DOUBLE,&nodes_[0],partitionNodeCounts[mypeno_],MPI::DOUBLE,0,MPI_COMM_WORLD);
      }

    //memory
    partitionNodeCounts   = vector<int>();
    partitionNodeOffset   = vector<int>();
    allNodesListNonUnique = vector<double>();
    if (mypeno_ == 0) cout << "----------------------------" << endl;

    //set up element structure (nnodes, ndof, etc)
    setupElements(nonUniqueNodeCount,globalNodeList);

    //apply FEM BC while non-unique node list is available
    applyBC(nonUniqueNodeCount,globalNodeList);
    if (mypeno_ == 0) cout << "----------------------------" << endl;

    //read FEM loading from file while non-unique node list is available
    setForce(nonUniqueNodeCount,globalNodeList);
    if (mypeno_ == 0) cout << "----------------------------" << endl;

    //memory
    globalNodeList = vector<int>();

    //identify one layer of elements deep outside of my partition
    identifyHaloElements(elementTypesList_loc,temporary_nnodes_offset,part,rootComponentIDs);

    //compute normal vectors at nodes, use halo elements for parallel
    initializeNormalVectors(nonUniqueNodeCount,globalNodeList,elementTypesList_loc,temporary_nnodes_offset,rootComponentIDs);

#endif
  }
  
  void fem::treatSingleProcFEM(vector<int>& elementTypesList,
                               vector<int>& rootComponentIDs)
  {

    //true defaults to full LU with MUMPS --- good for poorly conditioned systems
    //false goes into the LU pre-conditioned GMRES
    parallelSolve_ = true;
    //when DEDICATED_FEM_PROC is set to 1, a parallel method is employed with only the root proc
    //handling the FEM solvers and load/displacement transfer. This rountine will handle the
    //initialization of the shared arrays used for parallel.

    Screen::MasterInfo("Treating dedicated FEM processor");

    myElementCount_     = nof_faces_glob_;
    myNodeCount_        = nof_nodes_glob_;
    myElementCount_ext_ = nof_faces_glob_;
    myNodeCount_ext_    = nof_nodes_glob_;

    //treat element global to local mapping
    for (int e=0;e<nof_nodes_glob_;++e)
      {
        elemMaps_.elemLocal2Global.insert(pair<int,int> (e,e));
      }
    //fill nodes_
    nodes_.resize     (rootTotalNodes_.size());
    nodes0_.resize    (rootTotalNodes_.size());
    nodes_ext_.resize (rootTotalNodes_.size());
    nodes0_ext_.resize(rootTotalNodes_.size());

    nodes_.assign     (rootTotalNodes_.begin() ,rootTotalNodes_.end());
    nodes0_.assign    (rootTotalNodes0_.begin(),rootTotalNodes0_.end());
    nodes_ext_.assign (rootTotalNodes_.begin() ,rootTotalNodes_.end());
    nodes0_ext_.assign(rootTotalNodes0_.begin(),rootTotalNodes0_.end());
    //fill faces_
    faces_    .resize(rootTotalFaces_.size());
    faces_ext_.resize(rootTotalFaces_.size());

    faces_.assign    (rootTotalFaces_.begin(),rootTotalFaces_.end());
    faces_ext_.assign(rootTotalFaces_.begin(),rootTotalFaces_.end());

    int nonUniqueNodeCount = faces_.size();    
    setupElements(nonUniqueNodeCount,faces_);

    rootNnodesPerEl_.resize(nof_faces_glob_,-1);
    element elem;
    elementData_    .resize(nof_faces_glob_,elem);
    elementData_ext_.resize(nof_faces_glob_,elem);

    int offset = 0;
    vector<int> temporary_nnodes_offset(nof_faces_glob_,-1);
    for (int ng=0;ng<nof_faces_glob_;++ng)
      {

	if (revertNormalsByCompID_)
	  {
	    elementData_    [ng].compID = rootComponentIDs[ng];
	    elementData_ext_[ng].compID = rootComponentIDs[ng];
	  }

        if (elementTypesList[ng] == SIMPLETRI)
          {
            rootNnodesPerEl_[ng]    = 3;
            elementData_[ng].nnodes = 3;
            elementData_[ng].ndof   = 6;
	    elementData_[ng].type   = SIMPLETRI;
            elementData_ext_[ng].nnodes = 3;
            elementData_ext_[ng].ndof   = 6;
          }
        if (elementTypesList[ng] == LMITC3)
          {
            rootNnodesPerEl_[ng]    = 3;
            elementData_[ng].nnodes = 3;
            elementData_[ng].ndof   = 6;
	    elementData_[ng].type   = LMITC3;
            elementData_ext_[ng].nnodes = 3;
            elementData_ext_[ng].ndof   = 6;
          }
        if (elementTypesList[ng] == MITC3)
          {
            rootNnodesPerEl_[ng]    = 3;
            elementData_[ng].nnodes = 3;
            elementData_[ng].ndof   = 5;
	    elementData_[ng].type   = MITC3;
            elementData_ext_[ng].nnodes = 3;
            elementData_ext_[ng].ndof   = 5;
          }
        if (elementTypesList[ng] == BEAM)
          {
            rootNnodesPerEl_[ng]    = 2;
            elementData_[ng].nnodes = 2;
            elementData_[ng].ndof   = 2;
	    elementData_[ng].type   = BEAM;
            elementData_ext_[ng].nnodes = 2;
            elementData_ext_[ng].ndof   = 2;
          }
        if (elementTypesList[ng] == CABLE)
          {
            rootNnodesPerEl_[ng]    = 2;
            elementData_[ng].nnodes = 2;
            elementData_[ng].ndof   = 3;
	    elementData_[ng].type   = CABLE;
            elementData_ext_[ng].nnodes = 2;
            elementData_ext_[ng].ndof   = 3;
          }

        temporary_nnodes_offset[ng] = offset;
        offset += rootNnodesPerEl_[ng];
      }

    nnodesLocalOffset_.resize(nof_faces_glob_,-1);
    nnodesLocalOffset_ext_.resize(nof_faces_glob_,-1);
    int index = 0;
    for (int e=0;e<nof_faces_glob_;++e)
      {
        nnodesLocalOffset_[e] = index;
	nnodesLocalOffset_ext_[e] = index;
        index += rootNnodesPerEl_[e];
      }

    myNumberOfEachElementType_.resize(numberOfElementTypes_,0);
    //treat node global to local mapping
    for (int e=0;e<nof_faces_glob_;++e)
      {
        int nnodes_loc = elementData_[e].nnodes;
        elementData_[e].type     = elementTypesList[e];
        elementData_ext_[e].type = elementTypesList[e];

        for (int ne=0;ne<numberOfElementTypes_;++ne)
          {
            if (elementData_[e].type == elementTypes_[ne])
              {
                myNumberOfEachElementType_[ne] += 1;
              }
          }

        for (int n=0;n<nnodes_loc;++n)
          {
            int index = n + e*nnodes_loc;
            int node  = faces_[index];
            elemMaps_.nodeLocal2Global.insert(pair<int,int>(node,node));
            elemMaps_.nodeLocal2Global_ext.insert(pair<int,int>(node,node));
          }
      }

    index = 0;
    rootNdofOffset_.resize(nof_nodes_glob_,-1);
    for (int n=0;n<nof_nodes_glob_;++n)
      {
        rootNdofOffset_[n] = index;
        index += rootNdofPerNode_[n];
      }
    //apply BC
    applyBC(nonUniqueNodeCount,faces_);
    //set forces
    setForce(nonUniqueNodeCount,faces_);
    //initialize normal vectors
    initializeNormalVectors(nonUniqueNodeCount,faces_,elementTypesList,temporary_nnodes_offset,rootComponentIDs);

    //memory
    rootComponentIDs = vector<int>();

  }

  void fem::initializeNormalVectors(int& nonUniqueNodeCount,
                                    vector<int>& globalNodeList,
                                    vector<int>& elementTypesList,
                                    vector<int>& nnodesOffset,
                                    vector<int>& rootComponentIDs)
  {
    //initialize FEM normal vectors on the root proc for general meshes
    Screen::MasterInfo("Computing initial normal vectors");

    int numSharedElements = 0;
    vector<double> vec1(xyz,0.);
    vector<double> vec2(xyz,0.);
    vector<double> vec3(xyz,0.);
    vector<double> vec4(xyz,0.);
    vector<double> nonUniqueInitNormVecList;
    vector<double> nonUniqueInitNormVecList0;
    vector<double> Vn0;
    vector<double> Vnt;
    vector<int> participantNodeList; 
    int numberOfParticipantNodes = 0;
#if (DEDICATED_FEM_PROC == 0)

    Screen::MasterInfo("Identifying normal vector participants");

    //store elements participating in normal vector computation
    // --- right now tris only ---
    vector<int> myParticipantElementList;
    int myNumberOfParticipantElements = 0;
    for (int f=0;f<myElementCount_;++f)
      {
	if (elementData_[f].type == MITC3)
	  {
	    myParticipantElementList.push_back(f);
	    myNumberOfParticipantElements += 1;
	  }
      }

    //store normal vector participating nodes
    myNumberOfNormalVectorNodes_ = 0;
    myLocalNodesNormalVectorStatus_.resize(myNodeCount_,-1);

    int offset = 0;
    for (int e=0;e<myElementCount_;++e)
      {
        int nnodes_loc = elementData_[e].nnodes;
        for (int pe=0;pe<myNumberOfParticipantElements;++pe)
          {
            int participantID = myParticipantElementList[pe];

            if (e == participantID)
              {
                for (int n=0;n<nnodes_loc;++n)
                  {
                    int index = n + offset;
                    int node  = faces_[index];
                    if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),node) == myNormalVectorNodes_.end())
                      {
                        myLocalNodesNormalVectorStatus_[node] = myNumberOfNormalVectorNodes_;
			myNormalVectorNodes_.push_back(node);
                        myNumberOfNormalVectorNodes_ += 1;                        
                      }
                  }
              }
          }
        offset += nnodes_loc;
      }

    //use ghost/halo elements to compute geometric node normals for displacement transfer
    Screen::MasterInfo("Initializing normal vectors via halo elements");
    Vn0_.resize(myNumberOfNormalVectorNodes_*xyz,0.);
    Vnt_.resize(myNumberOfNormalVectorNodes_*xyz,0.);
    Vnt_geom_.resize(myNumberOfNormalVectorNodes_*xyz,0.);

    if (revertNormalsByCompID_)
      {
	for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
	  {
	    int localNodeID = myNormalVectorNodes_[n];
	    numSharedElements = 0;

	   int globID = elemMaps_.nodeLocal2Global[localNodeID];
	    //loop over extended face list, include ghosts
	    for (int f=0;f<myElementCount_ext_;++f)
	      {
		//hard coded for triangles
		if (elementData_ext_[f].nnodes == 3)
		  {
		    int compID = elementData_ext_[f].compID;
		    int nnodesOffset_ext = nnodesLocalOffset_ext_[f];

		    //see which faces in the extended map owns this normal vector node
		    if ( faces_ext_[0 + nnodesOffset_ext] == localNodeID ||
			 faces_ext_[1 + nnodesOffset_ext] == localNodeID || 
			 faces_ext_[2 + nnodesOffset_ext] == localNodeID )
		      {
			numSharedElements += 1;

			if ( find(revertNormalsCompIDList_.begin(),revertNormalsCompIDList_.end(),abs(compID)) == revertNormalsCompIDList_.end())
			  {			   
			    for (int idir=0;idir<xyz;++idir)
			      {
				vec1[idir] = nodes0_ext_[faces_ext_[1 + nnodesOffset_ext]*xyz + idir] - nodes0_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
				vec2[idir] = nodes0_ext_[faces_ext_[2 + nnodesOffset_ext]*xyz + idir] - nodes0_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
			      }
			  }
			else
			  {
			    for (int idir=0;idir<xyz;++idir)
			      {
				vec2[idir] = nodes0_ext_[faces_ext_[1 + nnodesOffset_ext]*xyz + idir] - nodes0_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
				vec1[idir] = nodes0_ext_[faces_ext_[2 + nnodesOffset_ext]*xyz + idir] - nodes0_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
			      }
			  }

			Vnt_[n*xyz + 0] += vec1[1]*vec2[2] - vec1[2]*vec2[1];
			Vnt_[n*xyz + 1] += vec1[2]*vec2[0] - vec1[0]*vec2[2];
			Vnt_[n*xyz + 2] += vec1[0]*vec2[1] - vec1[1]*vec2[0];

			Vn0_[n*xyz + 0] += vec1[1]*vec2[2] - vec1[2]*vec2[1];
			Vn0_[n*xyz + 1] += vec1[2]*vec2[0] - vec1[0]*vec2[2];
			Vn0_[n*xyz + 2] += vec1[0]*vec2[1] - vec1[1]*vec2[0];

		      }
		  }
	      }
	    for (int idir=0;idir<xyz;++idir) Vnt_[n*xyz + idir] = Vnt_[n*xyz + idir]/float(numSharedElements);
	    normalizeVector3(Vnt_);
	    for (int idir=0;idir<xyz;++idir) Vn0_[n*xyz + idir] = Vn0_[n*xyz + idir]/float(numSharedElements);
	    normalizeVector3(Vn0_);
	  }
      }
    else
      {
	for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
	  {
	    int localNodeID = myNormalVectorNodes_[n];
	    numSharedElements = 0;

	    //loop over extended face list, include ghosts
	    for (int f=0;f<myElementCount_ext_;++f)
	      {
		//hard coded for triangles
		if (elementData_ext_[f].nnodes == 3)
		  {
		    int nnodesOffset_ext = nnodesLocalOffset_ext_[f];

		    //see which faces in the extended map owns this normal vector node
		    if ( faces_ext_[0 + nnodesOffset_ext] == localNodeID ||
			 faces_ext_[1 + nnodesOffset_ext] == localNodeID || 
			 faces_ext_[2 + nnodesOffset_ext] == localNodeID )
		      {
			numSharedElements += 1;

			for (int idir=0;idir<xyz;++idir)
			  {
			    vec1[idir] = nodes0_ext_[faces_ext_[1 + nnodesOffset_ext]*xyz + idir] - nodes0_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
			    vec2[idir] = nodes0_ext_[faces_ext_[2 + nnodesOffset_ext]*xyz + idir] - nodes0_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
			  }

			Vnt_[n*xyz + 0] += vec1[1]*vec2[2] - vec1[2]*vec2[1];
			Vnt_[n*xyz + 1] += vec1[2]*vec2[0] - vec1[0]*vec2[2];
			Vnt_[n*xyz + 2] += vec1[0]*vec2[1] - vec1[1]*vec2[0];

			Vn0_[n*xyz + 0] += vec1[1]*vec2[2] - vec1[2]*vec2[1];
			Vn0_[n*xyz + 1] += vec1[2]*vec2[0] - vec1[0]*vec2[2];
			Vn0_[n*xyz + 2] += vec1[0]*vec2[1] - vec1[1]*vec2[0];

		      }
		  }
	      }
	    for (int idir=0;idir<xyz;++idir) Vnt_[n*xyz + idir] = Vnt_[n*xyz + idir]/float(numSharedElements);
	    normalizeVector3(Vnt_);
	    for (int idir=0;idir<xyz;++idir) Vn0_[n*xyz + idir] = Vn0_[n*xyz + idir]/float(numSharedElements);
	    normalizeVector3(Vn0_);
	  }
      }
    //allocate director vectors to participating node count
    V1t_.resize          (myNumberOfNormalVectorNodes_*xyz,0.);
    V10_.resize          (myNumberOfNormalVectorNodes_*xyz,0.);
    V20_.resize          (myNumberOfNormalVectorNodes_*xyz,0.);
    V2t_.resize          (myNumberOfNormalVectorNodes_*xyz,0.);
    deltaVnt_.resize     (myNumberOfNormalVectorNodes_*xyz,0.);
    deltaVnt_bak_.resize (myNumberOfNormalVectorNodes_*xyz,0.);
    myDirVecs_.resize    (myNumberOfNormalVectorNodes_*xyz,0.);
    myDirVecs0_.resize    (myNumberOfNormalVectorNodes_*xyz,0.);
    myDirVecs_bak_.resize(myNumberOfNormalVectorNodes_*xyz,0.);
    Vnt_bak_.resize      (myNumberOfNormalVectorNodes_*xyz,0.);

#else
    //store elements participating in normal vector computation
    // --- right now tris only ---
    vector<int> participantElementList;
    int numberOfParticipantElements = 0;
    for (int f=0;f<nof_faces_glob_;++f)
      {
        if (elementTypesList[f] == MITC3)
          {
            participantElementList.push_back(f);
            numberOfParticipantElements += 1;
          }
      }
    //memory
    elementTypesList = vector<int>();

    //store participating nodes --- and on the fly, number of non-unique participating nodes
    numberOfParticipantNodes = 0;
    int nonUniquenumberOfParticipantNodes = 0;
    for (int pe=0;pe<numberOfParticipantElements;++pe)
      {
        int participantElementID = participantElementList[pe];
        int offset               = nnodesOffset[participantElementID];
        int nnodes_loc           = rootNnodesPerEl_[participantElementID];

        for (int n=0;n<nnodes_loc;++n)
          {
            int index = n + offset;
            int node  = rootTotalFaces_[index];
            if (find(participantNodeList.begin(),participantNodeList.end(),node) == participantNodeList.end())
              {
                participantNodeList.push_back(node);
                numberOfParticipantNodes += 1;
              }
            nonUniquenumberOfParticipantNodes += 1;
          }
      }

    Vn0_.resize(numberOfParticipantNodes*xyz,0.);
    Vnt_.resize(numberOfParticipantNodes*xyz,0.);
    Vnt_geom_.resize(numberOfParticipantNodes*xyz,0.);

    for (int n=0;n<numberOfParticipantNodes;++n)
      {

        int node = participantNodeList[n];
        numSharedElements = 0;

        //loop over participating elements
        for (int f=0;f<numberOfParticipantElements;++f)
          {
            int participantElementID = participantElementList[f];
            int nnodes_loc           = rootNnodesPerEl_[participantElementID];
            int nnodesOffset_loc     = nnodesOffset[participantElementID];

            //count all faces that own node n
            if ( (rootTotalFaces_[nnodesOffset_loc] == node) || (rootTotalFaces_[1 + nnodesOffset_loc] == node) || (rootTotalFaces_[2 + nnodesOffset_loc] == node) )
              {
                numSharedElements += 1;

                for (int idir=0;idir<xyz;++idir)
                  {
                    // Compute vectors along triangle face edges that own node n
                    vec1[idir] = rootTotalNodes_ [rootTotalFaces_[1 + nnodesOffset_loc]*xyz + idir] - rootTotalNodes_ [rootTotalFaces_[nnodesOffset_loc]*xyz + idir];
                    vec2[idir] = rootTotalNodes_ [rootTotalFaces_[2 + nnodesOffset_loc]*xyz + idir] - rootTotalNodes_ [rootTotalFaces_[nnodesOffset_loc]*xyz + idir];
                    vec3[idir] = rootTotalNodes0_[rootTotalFaces_[1 + nnodesOffset_loc]*xyz + idir] - rootTotalNodes0_[rootTotalFaces_[nnodesOffset_loc]*xyz + idir];
                    vec4[idir] = rootTotalNodes0_[rootTotalFaces_[2 + nnodesOffset_loc]*xyz + idir] - rootTotalNodes0_[rootTotalFaces_[nnodesOffset_loc]*xyz + idir];
                  }
                // Add vectors from all faces that share node n
                Vnt_[n*xyz + 0] += vec1[1]*vec2[2] - vec1[2]*vec2[1];
                Vnt_[n*xyz + 1] += vec1[2]*vec2[0] - vec1[0]*vec2[2];
                Vnt_[n*xyz + 2] += vec1[0]*vec2[1] - vec1[1]*vec2[0];

                Vn0_[n*xyz + 0] += vec3[1]*vec4[2] - vec3[2]*vec4[1];
                Vn0_[n*xyz + 1] += vec3[2]*vec4[0] - vec3[0]*vec4[2];
                Vn0_[n*xyz + 2] += vec3[0]*vec4[1] - vec3[1]*vec4[0];

              }

          }
        //get average normal vector from surround elements
        for (int idir=0;idir<xyz;++idir)
          {
            int index  = n*xyz + idir;
            Vn0_[index] = Vn0_[index]/float(numSharedElements);
            Vnt_[index] = Vnt_[index]/float(numSharedElements);
          }
      }

    normalizeVector3(Vn0_);
    normalizeVector3(Vnt_);

    myNumberOfNormalVectorNodes_ = 0;
    myLocalNodesNormalVectorStatus_.resize(myNodeCount_,-1);
    for (int n=0;n<myNodeCount_;++n)
      {
        int nGlob = elemMaps_.nodeLocal2Global[n];

        for (int pn=0;pn<numberOfParticipantNodes;++pn)
          {
            int participantNode = participantNodeList[pn];

            if (nGlob == participantNode)
              {
                if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),n) == myNormalVectorNodes_.end())
                  {
                    myLocalNodesNormalVectorStatus_[n] = myNumberOfNormalVectorNodes_;
                    myNormalVectorNodes_.push_back(n);
                    myNumberOfNormalVectorNodes_ += 1;
                  }
              }
          }
      }

    //allocate director vectors to participating node count
    V1t_.resize          (myNumberOfNormalVectorNodes_*xyz,0.);
    V10_.resize          (myNumberOfNormalVectorNodes_*xyz,0.);
    V20_.resize          (myNumberOfNormalVectorNodes_*xyz,0.);
    V2t_.resize          (myNumberOfNormalVectorNodes_*xyz,0.);
    deltaVnt_.resize     (myNumberOfNormalVectorNodes_*xyz,0.);
    deltaVnt_bak_.resize (myNumberOfNormalVectorNodes_*xyz,0.);
    myDirVecs_.resize    (myNumberOfNormalVectorNodes_*xyz,0.);
    myDirVecs0_.resize    (myNumberOfNormalVectorNodes_*xyz,0.);
    myDirVecs_bak_.resize(myNumberOfNormalVectorNodes_*xyz,0.);
    Vnt_bak_.resize      (myNumberOfNormalVectorNodes_*xyz,0.);

#endif
  }

#if (USEPETSC>0)
  void fem::identifyHaloElements(vector<int>& elementTypesList_loc,
				 vector<int>& nnodesOffset,
				 idx_t *part,
				 vector<int>& rootComponentIDs)
  {
    //find 1 layer deep outside my partition boundary of ghost/halo elements
    //ultimately this is great for system assembly without communication
    //it is also necessary for computing geometric node normal vectors along proc boundaries

    Screen::MasterInfo("Identifying halo elements around my partition");

    const int MAX_NEIGH = 500;

    //first step, give all procs access to relevant variables

    //element connectivity
    int facesSize;
    if (mypeno_ == 0) facesSize = rootTotalFaces_.size();
    ierr = MPI_Bcast(&facesSize,1,MPI::INTEGER,0,MPI_COMM_WORLD);
    if (mypeno_ != 0) rootTotalFaces_.resize(facesSize,-1);
    ierr = MPI_Bcast(&rootTotalFaces_[0],facesSize,MPI::INTEGER,0,MPI_COMM_WORLD);

    //element type
    if (mypeno_ != 0) elementTypesList_loc.resize(nof_faces_glob_,-1);
    ierr = MPI_Bcast(&elementTypesList_loc[0],nof_faces_glob_,MPI::INTEGER,0,MPI_COMM_WORLD);

    //nnodes per element
    if (mypeno_ != 0) rootNnodesPerEl_.resize(nof_faces_glob_,-1);
    ierr = MPI_Bcast(&rootNnodesPerEl_[0],nof_faces_glob_,MPI::INTEGER,0,MPI_COMM_WORLD);

    //nnodes offset per element
    if (mypeno_ != 0) nnodesOffset.resize(nof_faces_glob_,-1);
    ierr = MPI_Bcast(&nnodesOffset[0],nof_faces_glob_,MPI::INTEGER,0,MPI_COMM_WORLD);

    //temporary node neighbor storage
    vector<vector<int> > nodeNeigh_tmp(nof_nodes_glob_,vector<int>(MAX_NEIGH,-1));
    vector<int> nodeNeighCount_tmp(nof_nodes_glob_,0);
    vector<vector<int> > nodeNeigh(nof_faces_glob_,vector<int>(MAX_NEIGH,-1));
    vector<int> nodeNeighCount(nof_faces_glob_,0);
    vector<int> tmp(MAX_NEIGH);

    Screen::MasterInfo("Finding global node neighbors over the entire triangulation");
    //next, each process needs to determine their node neighbors
    for (int eglob=0;eglob<nof_faces_glob_;++eglob)
      {

	int eTypeGlob          = elementTypesList_loc[eglob];
	int nnodes_offset_glob = nnodesOffset[eglob];

	if (eTypeGlob == MITC3)
	  {

	    int node1 = rootTotalFaces_[0 + nnodes_offset_glob];
	    int node2 = rootTotalFaces_[1 + nnodes_offset_glob];
	    int node3 = rootTotalFaces_[2 + nnodes_offset_glob];

	    nodeNeigh_tmp[node1][nodeNeighCount_tmp[node1]] = eglob;
	    nodeNeigh_tmp[node2][nodeNeighCount_tmp[node2]] = eglob;
	    nodeNeigh_tmp[node3][nodeNeighCount_tmp[node3]] = eglob;

	    nodeNeighCount_tmp[node1] += 1;
	    nodeNeighCount_tmp[node2] += 1;
	    nodeNeighCount_tmp[node3] += 1;

	  }
	else
	  {

	    int node1 = rootTotalFaces_[0 + nnodes_offset_glob];
	    int node2 = rootTotalFaces_[1 + nnodes_offset_glob];

	    nodeNeigh_tmp[node1][nodeNeighCount_tmp[node1]] = eglob;
	    nodeNeigh_tmp[node2][nodeNeighCount_tmp[node2]] = eglob;

	    nodeNeighCount_tmp[node1] += 1;
	    nodeNeighCount_tmp[node2] += 1;

	  }
      }

    for (int eglob=0;eglob<nof_faces_glob_;++eglob)
      {

	int eTypeGlob          = elementTypesList_loc[eglob];
	int nnodesGlob         = rootNnodesPerEl_[eglob];
	int nnodes_offset_glob = nnodesOffset[eglob];
	fill(tmp.begin(),tmp.end(),-1);

	//loop over face nodes
	for (int n1=0;n1<nnodesGlob;++n1)
	  {
	    //get node global ID
	    int node = rootTotalFaces_[n1 + nnodes_offset_glob];

	    //loop over faces that own this global node ID
	    for (int n2=0;n2<nodeNeighCount_tmp[node];++n2)
	      {

		int eglob2   = nodeNeigh_tmp[node][n2];
		int shareCountLimit = 1;
		int shareCount = 0;

		for (int nn=0;nn<nodeNeighCount[eglob];++nn) tmp[nn] = nodeNeigh[eglob][nn];
		if (find(tmp.begin(),tmp.end(),eglob2) == tmp.end())
		  {
		    nodeNeigh[eglob][nodeNeighCount[eglob]] = eglob2;
		    nodeNeighCount[eglob] += 1;
		    if (nodeNeighCount[eglob] >= MAX_NEIGH)
		      {
			Screen::MasterError("node neigh count too big!");
		      }
		  }
	      }
	  }
      }
    //memory
    nodeNeigh_tmp = vector<vector<int> >();
    nodeNeighCount_tmp = vector<int>();

    vector<vector<int> > myGlobalNodeNeighbors(myElementCount_,vector<int>(MAX_NEIGH));
    vector<int> myGlobalNodeNeighborCount(myElementCount_,-1);

    Screen::MasterInfo("Finding global node neighbors on my process");
    //reduce the total map to my faces only
    for (int f=0;f<myElementCount_;++f)
      {
	int eglob = elemMaps_.elemLocal2Global[f];
	for (int count=0;count<nodeNeighCount[eglob];++count) myGlobalNodeNeighbors[f][count] = nodeNeigh[eglob][count];
	myGlobalNodeNeighborCount[f] = nodeNeighCount[eglob];
      }
    //memory
    nodeNeigh = vector<vector<int> >();
    nodeNeighCount = vector<int>();
    tmp = vector<int>();

    //store tmp global element list
    tmp.resize(myElementCount_,-1);
    for (int f=0;f<myElementCount_;++f) tmp[f] = elemMaps_.elemLocal2Global[f];

    myNumberOfHaloElements_ = 0;
    //find which of my node neighbors are on a different process
    Screen::MasterInfo("Identfying global IDs of my halo elements");
    for (int f=0;f<myElementCount_;++f)
      {
	for (int count=0;count<myGlobalNodeNeighborCount[f];++count)
	  {
	    int nodeNeighborGlobalID = myGlobalNodeNeighbors[f][count];
	    if (find(tmp.begin(),tmp.end(),nodeNeighborGlobalID) == tmp.end() &&
		find(myHaloElementGlobIDs_.begin(),myHaloElementGlobIDs_.end(),nodeNeighborGlobalID) == myHaloElementGlobIDs_.end())
	      {
		myHaloElementGlobIDs_.push_back(nodeNeighborGlobalID);
		myNumberOfHaloElements_ += 1;
	      }
	  }
      }
    //memory
    myGlobalNodeNeighbors = vector<vector<int> >();
    myGlobalNodeNeighborCount = vector<int>();

    element elem;
    elementData_ext_.resize(myElementCount_+myNumberOfHaloElements_,elem);

    //broadcast root ndof offset and ndof per node to store global dof of halo nodes
    if (mypeno_ != 0) rootNdofOffset_.resize(nof_nodes_glob_,-1);
    ierr = MPI_Bcast(&rootNdofOffset_[0],nof_nodes_glob_,MPI::INTEGER,0,MPI_COMM_WORLD);
    if (mypeno_ != 0) rootNdofPerNode_.resize(nof_nodes_glob_,-1);
    ierr = MPI_Bcast(&rootNdofPerNode_[0],nof_nodes_glob_,MPI::INTEGER,0,MPI_COMM_WORLD);

    //broadcast root comp IDs
    if (revertNormalsByCompID_)
      {
	if (mypeno_ != 0) rootComponentIDs.resize(nof_faces_glob_,-1);
	ierr = MPI_Bcast(&rootComponentIDs[0],nof_faces_glob_,MPI::INTEGER,0,MPI_COMM_WORLD);    
      }

    //keep halo elements global dof indices
    Screen::MasterInfo("Storing ghost element global connectivity");

    //building extended storage variables for ndof, nnodes, faces, nodes, etc
    //first, copy original mappings
    vector<int> uniqueFaceNodes;
    int offset = 0;
    int ndofLocalOffsetCounter = 0;
    myNodeCount_ext_ = 0;
    for (int e=0;e<myElementCount_;++e)
      {

        int nnodes_loc = elementData_[e].nnodes;
	int globE = elemMaps_.elemLocal2Global[e];

	nnodesLocalOffset_ext_.push_back(offset);
	elementData_ext_[e].nnodes = nnodes_loc;
	elementData_ext_[e].globID = globE;
	elemMaps_.elemLocal2Global_ext.insert(pair<int,int>(e,globE));
	if (revertNormalsByCompID_)
	  {
	    elementData_ext_[e].compID = rootComponentIDs[globE];
	    elementData_[e].compID = rootComponentIDs[globE];
	  }

        for (int n=0;n<nnodes_loc;++n)
          {

            int index = n + offset;
            int node  = faces_[index];
	    int nodeGlobal = elemMaps_.nodeLocal2Global[node];
	    int ndofPerNode = ndofPerNode_[node];
	    faces_ext_.push_back(nodeGlobal);

            if (find(uniqueFaceNodes.begin(),uniqueFaceNodes.end(),nodeGlobal) == uniqueFaceNodes.end())
              {

		ndofPerNode_ext_.push_back(rootNdofPerNode_[nodeGlobal]);
		ndofLocalOffset_ext_.push_back(ndofLocalOffsetCounter);
                uniqueFaceNodes.push_back(nodeGlobal);
                elemMaps_.nodeLocal2Global_ext.insert(pair<int,int>(myNodeCount_ext_,nodeGlobal));
                myNodeCount_ext_++;
		ndofLocalOffsetCounter += ndofPerNode;

              }
          }
        offset += nnodes_loc;
      }
    //getting indices for global DOF of ghost element nodes
    int counter = 0;
    vector<int> localHaloDOF;
    vector<int> globalHaloDOF;
    for (int n=0;n<myNodeCount_;++n)
      {
        int ndof_loc = ndofPerNode_[n];
        for (int dof=0;dof<ndof_loc;++dof)
          {
            localHaloDOF.push_back(counter);
            globalHaloDOF.push_back(dof + ndofGlobalOffset_[n]);
            counter += 1;
          }	
      }

    //finish local 2 global node map with ghosts and the global ghost DOF
    //next, extend mappings into the ghosts
    for (int e=0;e<myNumberOfHaloElements_;++e)
      {
	int globID_ext       = myHaloElementGlobIDs_[e];
	int nnodesOffset_ext = nnodesOffset[globID_ext];
	int nnodesPerEl_ext  = rootNnodesPerEl_[globID_ext];

	nnodesLocalOffset_ext_.push_back(offset);
	elementData_ext_[e+myElementCount_].nnodes = nnodesPerEl_ext;
	elementData_ext_[e+myElementCount_].globID = globID_ext;
	elemMaps_.elemLocal2Global_ext.insert(pair<int,int>(e+myElementCount_,globID_ext));
	if (revertNormalsByCompID_) elementData_ext_[e+myElementCount_].compID = rootComponentIDs[globID_ext];

	for (int fn=0;fn<nnodesPerEl_ext;++fn)
	  {
	    int globNodeID_ext  = rootTotalFaces_[fn + nnodesOffset_ext];
	    int ndofOffset_ext  = rootNdofOffset_[globNodeID_ext];
	    int ndofPerNode_ext = rootNdofPerNode_[globNodeID_ext];

	    faces_ext_.push_back(globNodeID_ext);
            if (find(uniqueFaceNodes.begin(),uniqueFaceNodes.end(),globNodeID_ext) == uniqueFaceNodes.end())
              {
		ndofPerNode_ext_.push_back(ndofPerNode_ext);
		ndofLocalOffset_ext_.push_back(ndofLocalOffsetCounter);
                uniqueFaceNodes.push_back(globNodeID_ext);
                elemMaps_.nodeLocal2Global_ext.insert(pair<int,int>(myNodeCount_ext_,globNodeID_ext));
		ndofLocalOffsetCounter += ndofPerNode_ext;
                myNodeCount_ext_++;

		for (int dof=0;dof<ndofPerNode_ext;++dof)
		  {
		    localHaloDOF.push_back(counter);
		    globalHaloDOF.push_back(dof + ndofOffset_ext);
		    counter += 1;
		  }

              }
	  }
	offset += nnodesPerEl_ext;
      }
    //memory
    uniqueFaceNodes = vector<int>();
    rootComponentIDs = vector<int>();

    myElementCount_ext_ = myElementCount_ + myNumberOfHaloElements_;
    //replace global node index in faces_ext_ with new local node index
    vector<int> facesCopy;
    facesCopy = faces_ext_;
    for (int n=0;n<myNodeCount_ext_;++n)
      {
        int globalNode = elemMaps_.nodeLocal2Global_ext[n];

        for (int i=0;i<faces_ext_.size();++i)
          {
            if (globalNode == faces_ext_[i])
              {
                facesCopy[i] = n;
              }
          }
      }
    faces_ext_.assign(facesCopy.begin(),facesCopy.end());
    facesCopy = vector<int>();
    //create extended scatter context
    VecCreateSeq(MPI_COMM_SELF,counter,&mySolutionVector_p_ext_);
    VecZeroEntries(mySolutionVector_p_ext_);
    ISCreateGeneral(MPI_COMM_SELF,counter,&localHaloDOF[0] ,PETSC_COPY_VALUES,&myLocalDOFIndexList_p_ext_ );
    ISCreateGeneral(MPI_COMM_SELF,counter,&globalHaloDOF[0],PETSC_COPY_VALUES,&myGlobalDOFIndexList_p_ext_);    

    //memory
    localHaloDOF  = vector<int>();
    globalHaloDOF = vector<int>();

    //memory
    if (mypeno_ != 0)
      {
	rootTotalFaces_ = vector<int>();
	elementTypesList_loc = vector<int>();
	rootNnodesPerEl_ = vector<int>();
	nnodesOffset = vector<int>();
	rootNdofOffset_ = vector<int>();
      }

    //fill nodes0_ext_ and nodes_ext_ with root information
    if (mypeno_ != 0) rootTotalNodes_.resize(nof_nodes_glob_*xyz,-1.);
    ierr = MPI_Bcast(&rootTotalNodes_[0],nof_nodes_glob_*xyz,MPI::DOUBLE,0,MPI_COMM_WORLD);
    if (mypeno_ != 0) rootTotalNodes0_.resize(nof_nodes_glob_*xyz,-1.);
    ierr = MPI_Bcast(&rootTotalNodes0_[0],nof_nodes_glob_*xyz,MPI::DOUBLE,0,MPI_COMM_WORLD);

    for (int ne=0;ne<myNodeCount_ext_;++ne)
      {
	int globalNodeID = elemMaps_.nodeLocal2Global_ext[ne];

	for (int dir=0;dir<xyz;++dir)
	  {	   
	    nodes0_ext_.push_back(rootTotalNodes0_[globalNodeID*xyz + dir]);
	    nodes_ext_ .push_back(rootTotalNodes_ [globalNodeID*xyz + dir]);
	  }
      }
    //memory
    if (mypeno_ != 0) rootTotalNodes_  = vector<double>();
    if (mypeno_ != 0) rootTotalNodes0_ = vector<double>();

    Screen::MasterInfo("Layer of halo/ghost elements identified");

    //Screen::MasterInfo("Dumping ghost partition");
    //string filename = "debugFEM/ghostPartition_mypeno_"+to_string(mypeno_)+".vtk";
    //ofstream file;
    //file.open(filename);
    //if(file.fail()) Screen::MasterError("Bad file in writeTriangulation function");
    //file.precision(10);
    //file.setf(ios::fixed);
    //file.setf(ios::showpoint);
    //file << "# ";
    //file << "vtk ";
    //file << "DataFile ";
    //file << "Version ";
    //file << "3.0";
    //file << endl;
    //file << "vtk ";
    //file << "output";
    //file << endl;
    //file << "ASCII";
    //file << endl;
    //file << "DATASET ";
    //file << "POLYDATA";
    //file << endl;
    //file << "POINTS ";
    //file << myNodeCount_ext_;
    //file << " float";
    //file << endl;
    //for (int n=0;n<myNodeCount_ext_;++n)
    //  {
    //	for (int idir=0;idir<xyz;++idir)
    //	  {
    //	    file << nodes_ext_[n*xyz+idir] << " ";
    //	  }
    //	file << endl;
    //  }
    //file << "POINT_DATA " << myNodeCount_ext_ << endl;
    //file << "SCALARS Data double" << endl;
    //file << "LOOKUP_TABLE default" << endl;
    //for (int n=0;n<myNodeCount_ext_;++n) file << 1.0 << endl;
    //file.close();

  }
#endif

void fem::applyBC(int& nonUniqueNodeCount,
                    vector<int>& globalNodeList)
  {
#if (USEPETSC>0)
    Screen::MasterInfo("Applying structural boundary conditions");
#if (DEDICATED_FEM_PROC == 0)
    vector<int> nodebc_tmp;
    vector<int> nonUniqueBoundaryConditionList;
    if (mypeno_ == 0)
      {
        int size = 0;
        for (int n=0;n<nof_nodes_glob_;++n)
          {
            size += rootNdofPerNode_[n];
          }

        nodebc_tmp.resize(size,0);

        if (!regionalBC_)
          {
            //read boundary conditions from file
            ifstream nft;
            nft.open(fem_BC_file_);
            if (nft.fail() && fem_BC_file_ != "regional") Screen::MasterError("File not found for structural boundary conditions");

            nft >> numFixedNodes_;
            Screen::MasterInfo("Number of fixed structural nodes: "+to_string(numFixedNodes_));
            Screen::MasterInfo("BC file: "+fem_BC_file_);

            vector<int> fixedNode;
            vector<int> fixedDOF;
            vector<int> ndof_fixed;

            for (int fn=0;fn<numFixedNodes_;++fn)
              {

                //read fixed node ID
                int tmp;
                nft >> tmp;
                fixedNode.push_back(tmp);

                //read number of fixed DOF
                nft >> tmp;
                ndof_fixed.push_back(tmp);

                //read fixed DOF IDs
                for (int nd=0;nd<ndof_fixed[fn];++nd)
                  {
                    nft >> tmp;
                    fixedDOF.push_back(tmp);
                  }
              }
            //run along all indices and insert fixed BC
            int offset   = 0;
            int bcoffset = 0;
            for (int n=0;n<nof_nodes_glob_;++n)
              {
                for (int fn=0;fn<numFixedNodes_;++fn)
                  {
                    if (n == fixedNode[fn])
                      {
                        int ndof_loc = rootNdofPerNode_[n];

                        for (int d=0;d<ndof_fixed[fn];++d)
                          {
                            int dof = fixedDOF[d + bcoffset];
                            nodebc_tmp[dof + offset] = 1;
                          }
                        bcoffset += ndof_fixed[fn];
                      }
                  }
                offset += rootNdofPerNode_[n];
              }
          }
        else
          {
            //apply regional BC
            if (regionShape_ == "circle")
              {
                if (regionNormal_ == "z")
                  {
                    int offset = 0;
                    for (int fn=0;fn<nof_nodes_glob_;++fn)
                      {

                        offset = fn*rootNdofPerNode_[fn];

                        int ndof_loc = rootNdofPerNode_[fn];

                        int indexX = 0 + fn*xyz;
                        int indexY = 1 + fn*xyz;
                        int indexZ = 2 + fn*xyz;

                        double xx  = nodes_[indexX]*nodes_[indexX];
                        double yy  = nodes_[indexY]*nodes_[indexY];
                        double zz  = nodes_[indexZ]*nodes_[indexZ];

                        double mag = pow(xx+yy,0.5);
                        if (mag > 0.5*regionDiameter_)
                          {                            
                            for (int d=0;d<ndof_loc;++d)
                              {
                                nodebc_tmp[d + offset] = 1;
                              }
                          }
                      }
                  }
              }
          }

        size = 0;
        for (int n=0;n<nonUniqueNodeCount;++n)
          {
            int node = globalNodeList[n];
            size    += rootNdofPerNode_[node];
          }

        int index = 0;
        rootNdofOffset_.resize(nof_nodes_glob_,-1);
        for (int n=0;n<nof_nodes_glob_;++n)
          {
            rootNdofOffset_[n] = index;
            index += rootNdofPerNode_[n];
          }
        //perform non-unique, proc-based ordering for scattering
        nonUniqueBoundaryConditionList.resize(size,-99);
        int counter = 0;
        for (int n=0;n<nonUniqueNodeCount;++n)
          {
            int node     = globalNodeList[n];
            int ndof_loc = rootNdofPerNode_[node];
            int offset   = rootNdofOffset_[node];

            for (int dof=0;dof<ndof_loc;++dof)
              {
                int globalIndex                         = dof + offset;
                nonUniqueBoundaryConditionList[counter] = nodebc_tmp[globalIndex];
                counter += 1;
              }
          }
      }
    //memory
    nodebc_tmp = vector<int>();

    //compute offset and displs
    vector<int> partitionNodeCounts_with_ndof(numprocs_,-1);
    vector<int> partitionNodeOffset_with_ndof(numprocs_,-1);
    int myNodeOffset = myNodeCount_*mypeno_;

    MPI_Allgather(&myNodeCount_,1,MPI_INT,&partitionNodeCounts_with_ndof[0],1,MPI_INT,MPI_COMM_WORLD);

    partitionNodeOffset_with_ndof[0] = 0;
    for (int p=1;p<numprocs_;++p) partitionNodeOffset_with_ndof[p] = partitionNodeOffset_with_ndof[p-1] + partitionNodeCounts_with_ndof[p-1];
    //re-compute offset and displs to change index from node ID to number of DOF
    int sendcount = 0;
    for (int mn=0;mn<myNodeCount_;++mn)
      {
        sendcount += ndofPerNode_[mn];
      }

    MPI_Allgather(&sendcount,1,MPI_INT,&partitionNodeCounts_with_ndof[0],1,MPI::INT,MPI_COMM_WORLD);
    partitionNodeOffset_with_ndof[0] = 0;
    for (int np=1;np<numprocs_;++np) partitionNodeOffset_with_ndof[np] = partitionNodeOffset_with_ndof[np-1] + partitionNodeCounts_with_ndof[np-1];    
    //recover nodal BC from proc-ordered BC list
    nodeBC_.resize(sendcount,0);
    MPI_Scatterv(&nonUniqueBoundaryConditionList[0],&partitionNodeCounts_with_ndof[0],&partitionNodeOffset_with_ndof[0],MPI::INT,&nodeBC_[0],partitionNodeCounts_with_ndof[mypeno_],MPI::INT,0,MPI_COMM_WORLD);

    //get global DOF indices from nodeBC
    int local_offset  = 0;
    for (int n=0;n<myNodeCount_;++n)
      {
        int ndof_loc      = ndofPerNode_[n];
        int global_offset = ndofGlobalOffset_[n];
        int local_offset  = ndofLocalOffset_[n];

        for (int d=0;d<ndof_loc;++d)
          {
            if (nodeBC_[d + local_offset] == 1)
              {
                myFixedDOFIndexList_.push_back(d + global_offset);
              }
          }
      }

    myFixedDOFCount_ = myFixedDOFIndexList_.size();

#else
    //for serial FEM in parallel simulation
    int size = 0;
    for (int n=0;n<nof_nodes_glob_;++n)
      {
        size += rootNdofPerNode_[n];
      }

    nodeBC_.resize(size,0);

    if (!regionalBC_)
      {
        //read boundary conditions from file
        ifstream nft;
        nft.open(fem_BC_file_);
        if (nft.fail() && fem_BC_file_ != "regional") Screen::MasterError("File not found for structural boundary conditions");

        nft >> numFixedNodes_;
        Screen::MasterInfo("Number of fixed structural nodes: "+to_string(numFixedNodes_));
        Screen::MasterInfo("BC file: "+fem_BC_file_);

        vector<int> fixedNode;
        vector<int> fixedDOF;
        vector<int> ndof_fixed;

        for (int fn=0;fn<numFixedNodes_;++fn) 
          {

            //read fixed node ID
            int tmp;
            nft >> tmp;
            fixedNode.push_back(tmp);

            //read number of fixed DOF
            nft >> tmp;
            ndof_fixed.push_back(tmp);

            //read fixed DOF IDs
            for (int nd=0;nd<ndof_fixed[fn];++nd)
              {
                nft >> tmp;
                fixedDOF.push_back(tmp);
              }
          }
        //run along all indices and insert fixed BC
        int offset = 0;
        int bcoffset = 0;
        for (int n=0;n<nof_nodes_glob_;++n)
          {
            for (int fn=0;fn<numFixedNodes_;++fn)
              {
                if (n == fixedNode[fn])
                  {
                    int ndof_loc = rootNdofPerNode_[n];

                    for (int d=0;d<ndof_fixed[fn];++d)
                      {
                        int dof = fixedDOF[d + bcoffset];
                        nodeBC_[dof + offset] = 1;                        
                      }
                    bcoffset += ndof_fixed[fn];
                  }
              }
            offset += rootNdofPerNode_[n];
          }
      }
    else
      {
        //apply regional BC
        int count = 0;
        if (regionShape_ == "circle")
          {
            if (regionNormal_ == "z")
              {
                int offset = 0;
                for (int fn=0;fn<nof_nodes_glob_;++fn)
                  {

                    offset = fn*rootNdofPerNode_[fn];

                    int ndof_loc = rootNdofPerNode_[fn];

                    int indexX = 0 + fn*xyz;
                    int indexY = 1 + fn*xyz;
                    int indexZ = 2 + fn*xyz;

                    double xx  = nodes_[indexX]*nodes_[indexX];
                    double yy  = nodes_[indexY]*nodes_[indexY];
                    double zz  = nodes_[indexZ]*nodes_[indexZ];

                    double mag = pow(xx+yy,0.5);
                    if (mag > 0.5*regionDiameter_)
                      {
                        for (int d=0;d<ndof_loc;++d)
                          {
                            nodeBC_[d + offset] = 1;
                          }
                      }
                  }
              }
          }
        //dump regional BC to file for testing
        ofstream regional;
        regional.open("post/debugFEM/regionalFEMBC.vtk");
        if(regional.fail()) Screen::MasterError("Bad file in regional writer");
        regional.precision(10);
        regional.setf(ios::fixed);
        regional.setf(ios::showpoint);

        regional << "# vtk DataFile Version 3.0" << endl;
        regional << "vtk output" << endl;
        regional << "ASCII" << endl;
        regional << "DATASET POLYDATA" << endl;
        regional << "POINTS ";
        regional << count;
        regional << " float" << endl;

        for (int fn=0;fn<nof_nodes_glob_;++fn)
          {
            int ndof_loc = rootNdofPerNode_[fn];
            if (nodeBC_[fn*ndof_loc] == 1)
              {
                //x,y,z
                regional << nodes_[fn*xyz+0] << " " << nodes_[fn*xyz+1] << " " << nodes_[fn*xyz+2] << endl;
              }
          }

        regional << "POINT_DATA ";
        regional << count << endl;
        regional << "SCALARS RegionalNodes float" << endl;
        regional << "LOOKUP_TABLE default" << endl;
        for (int fn=0;fn<nof_nodes_glob_;++fn)
          {
            int ndof_loc = rootNdofPerNode_[fn];
            if (nodeBC_[fn*ndof_loc] == 1)
              {
                regional << 1. << endl;
              }
          }

      }    
    //get global DOF indices from nodeBC    
    for (int n=0;n<myNodeCount_;++n)
      {
        int ndof_loc = rootNdofPerNode_[n];
        int nGlob = elemMaps_.nodeLocal2Global[n];

        for (int d=0;d<ndof_loc;++d)
          {
            if (nodeBC_[d+ndof_loc*n] == 1)
              {
                myFixedDOFIndexList_.push_back(nGlob*ndof_loc + d);
              }
          }
      }

    myFixedDOFCount_ = myFixedDOFIndexList_.size();

#endif
  }
  
  void fem::ini3Das2D()
  {
    //set defaults
    dimLoc_ = DIM;
    nof_effective_nodes_ = nof_nodes_glob_;

    if (threeD_as_twoD_)
      {
        //introduce effective number of nodes per face
        dimLoc_ = 2;

        if (nof_nodes_glob_%2 != 0)
          {
            Screen::MasterError("Number of nodes in FEM mesh is odd, cannot use 3Das2D");
          }
        else
          {
            //introduce effective number of nodes in mesh
            nof_effective_nodes_ = nof_nodes_glob_/2;
          }
      }

#endif
  }
    
  void fem::initializeDirVecs()
  {
    dirVec_.resize(xyz,0);    
    if (mypeno_ == 0) cout << "----------------------------" << endl;
    Screen::MasterInfo("Initializing director vectors");

    int ddir = dirVecIndex_;
    for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
      {
        //find suitable dirVec
        if (dirVecIndex_ < 0)
          {
            int dd = 0;
            for (dd=0;dd<xyz;++dd)
              {
                dirVec_[dd] = 1.;
                double product = 0.;
                for (int idir=0;idir<xyz;++idir) product += Vn0_[n*xyz + idir]*dirVec_[idir];

                //if vectors are not aligned, good to go
                if ( abs(1. - abs(product)) > 1.E-4 )
                  {
                    ddir = dd;
                    break;
                  }
                else
                  {
                    dirVec_[dd] = 0.;
                  }
              }
            //check if no direction was ever found
            if ( (dd == 2) && (dirVec_[dd] < 1.E-8) ) Screen::MasterError("Could not find available dir vec index");

          }

	if (!do_restart_) myDirVecs_[n*xyz + ddir] = 1.;
	myDirVecs0_[n*xyz + ddir] = 1.;

        if (!do_restart_)
          {

            V1t_[n*xyz    ] = myDirVecs_[n*xyz+1]*Vnt_[n*xyz+2] - myDirVecs_[n*xyz+2]*Vnt_[n*xyz+1];
            V1t_[n*xyz + 1] = myDirVecs_[n*xyz+2]*Vnt_[n*xyz+0] - myDirVecs_[n*xyz+0]*Vnt_[n*xyz+2];
            V1t_[n*xyz + 2] = myDirVecs_[n*xyz+0]*Vnt_[n*xyz+1] - myDirVecs_[n*xyz+1]*Vnt_[n*xyz+0];

            V2t_[n*xyz    ] = Vnt_[n*xyz + 1]*V1t_[n*xyz + 2] - Vnt_[n*xyz + 2]*V1t_[n*xyz + 1];
            V2t_[n*xyz + 1] = Vnt_[n*xyz + 2]*V1t_[n*xyz + 0] - Vnt_[n*xyz + 0]*V1t_[n*xyz + 2];
            V2t_[n*xyz + 2] = Vnt_[n*xyz + 0]*V1t_[n*xyz + 1] - Vnt_[n*xyz + 1]*V1t_[n*xyz + 0];

          }

        V10_[n*xyz    ] = myDirVecs0_[n*xyz+1]*Vn0_[n*xyz+2] - myDirVecs0_[n*xyz+2]*Vn0_[n*xyz+1];
        V10_[n*xyz + 1] = myDirVecs0_[n*xyz+2]*Vn0_[n*xyz+0] - myDirVecs0_[n*xyz+0]*Vn0_[n*xyz+2];
        V10_[n*xyz + 2] = myDirVecs0_[n*xyz+0]*Vn0_[n*xyz+1] - myDirVecs0_[n*xyz+1]*Vn0_[n*xyz+0];

        V20_[n*xyz    ] = Vn0_[n*xyz + 1]*V10_[n*xyz + 2] - Vn0_[n*xyz + 2]*V10_[n*xyz + 1];
        V20_[n*xyz + 1] = Vn0_[n*xyz + 2]*V10_[n*xyz + 0] - Vn0_[n*xyz + 0]*V10_[n*xyz + 2];
        V20_[n*xyz + 2] = Vn0_[n*xyz + 0]*V10_[n*xyz + 1] - Vn0_[n*xyz + 1]*V10_[n*xyz + 0];

      }

    normalizeVector3(V1t_);
    normalizeVector3(V2t_);
    normalizeVector3(V10_);
    normalizeVector3(V20_);

    //removing normal vector calculation for beams for now
//    //beams - used for displacement transfer only
//    else if (elementType_ == BEAM)
//      {
//        vector<double> tangent(xyz,0.);
//        vector<double> intoDomain(xyz,0.);
//        intoDomain[2] = 1.;
//        for (int n=0;n<myNodeCount_-1;++n)
//          {
//            for (int idir=0;idir<xyz;++idir) tangent[idir] = nodes_[(n+1)*xyz+idir] - nodes_[n*xyz+idir];
//            Vnt_[n] = cross(intoDomain,tangent);
//          }
//        //treat end node
//        for (int idir=0;idir<xyz;++idir) tangent[idir] = nodes_[(myNodeCount_-1)*xyz+idir] - nodes_[(myNodeCount_-2)*xyz+idir];
//        Vnt_[myNodeCount_-1] = cross(intoDomain,tangent);
//        //normalize
//        for (int n=0;n<myNodeCount_;++n) normalizeVector(Vnt_[n]);
//      }
  }
  
  void fem::reInit()
  {
#if (USEPETSC>0)

    VecAssemblyBegin(loadVector_p_);
    VecAssemblyEnd  (loadVector_p_);
    VecAssemblyBegin(forceVector_p_);
    VecAssemblyEnd  (forceVector_p_);
    VecZeroEntries  (loadVector_p_);
    VecZeroEntries  (forceVector_p_);
    MatZeroEntries  (Kmat_p_);    
    if (is_unsteady_)
      {
        MatZeroEntries(mat_p_ );
        MatZeroEntries(Mmat_p_);
        MatZeroEntries(Cmat_p_);
      }

    fill(Fvec_.begin(),Fvec_.end(),0.);   
    for (int i=0;i<3;++i) 
      {
        fill(PKSTmat_    [i].begin(),PKSTmat_    [i].end(),0.);
        fill(GLstrainMat_[i].begin(),GLstrainMat_[i].end(),0.);
      }

#endif
  }

  void fem::auxillaryTransfer_ini()
  {    

    if (useAuxGeometry_)
      {
        if (do_restart_)
          {
            surfaceTriangulation surfTri;
            surfTri.readTriangulation_AUX_FSI(aux_nodes0_,auxFilename_);            
          }
        else
          {
            //read auxillary geometry
            //is_aux_node_ set inside
            surfaceTriangulation surfTri;
            int node = dimLoc_;
            surfTri.readTriangulation_AUX(aux_nodes_,aux_faces_,auxFilename_,numGeomNodes_,numGeomFaces_,node,is_aux_node_);

            aux_nodes0_.resize(numGeomNodes_*xyz,0.);
            aux_nodes0_.assign(aux_nodes_.begin(),aux_nodes_.end());

          }
      }
    else if (isFSI_)
      {
        assistedMappingFilename_ = "ultra2_remesh.vtk";
        if (do_restart_)
          {
#if (ASSISTED_MAPPING == 1)
            //set initial geometry nodal positions for mapping
            Screen::MasterInfo("Assisted mapping is being generated from "+assistedMappingFilename_);
            surfaceTriangulation surfTri;
            surfTri.readTriangulation_AUX_FSI(aux_nodes0_,assistedMappingFilename_);
#else
            surfaceTriangulation surfTri;
            surfTri.readTriangulation_AUX_FSI(aux_nodes0_,auxFilename_);
#endif

            if (shiftGeometry_)
              {
                for (int n=0;n<numGeomNodes_;++n)
                  {
                    for (int dir=0;dir<xyz;++dir)
                      {
                        aux_nodes0_[n*xyz+dir] += shiftGeomXYZ_[dir];
                      }
                  }
              }
          }
        else
          {
#if (ASSISTED_MAPPING == 1)
            Screen::MasterInfo("Assisted mapping is being generated from "+assistedMappingFilename_);
            vector<int> dummy;
            vector<int> dummy2;
            int dummyint1;
            int dummyint2;
            surfaceTriangulation surfTri;
            int node = dimLoc_;
            aux_nodes_.resize(0);
            surfTri.readTriangulation_AUX(aux_nodes_,dummy,assistedMappingFilename_,dummyint1,dummyint2,node,dummy2);
            dummy  = vector<int>();
            dummy2 = vector<int>();
#endif
            aux_nodes0_.resize(numGeomNodes_*xyz,0.);
            aux_nodes0_.assign(aux_nodes_.begin(),aux_nodes_.end());
          }
      }
    //handle FORTRAN indexing
    if (isFSI_) for (int i=0;i<aux_faces_.size();++i) aux_faces_[i] -= 1;
    if (mypeno_ == 0) cout << "----------------------------" << endl;

    //identify FEM elements participating in the auxillary transfer
    identifyAuxillaryParticipants();

#if (DEDICATED_FEM_PROC == 0)
    //obtain unqiue geometry-FEM element mapping for each processor
    obtainUniqueAuxillaryMapping();
#else
    treatSingleProcFSI();
#endif

    if (mypeno_ == 0) cout << "----------------------------" << endl;
    //create local-global mappings between geometry variables to save memory and facilitate transfer to fluid
    restructureFSI();
    ierr=MPI_Barrier(MPI_COMM_WORLD);

    //identify and store mappings from geometry nodes to FEM nodes
    initializeDisplacementTransfer();
    ierr=MPI_Barrier(MPI_COMM_WORLD);

    //identify and store mappings from FEM to geometry nodes
    initializeLoadTransfer();
    ierr=MPI_Barrier(MPI_COMM_WORLD);

    //compute stencils for interpolated displacement transfer
    computeDisplacementTransferStencils();
    ierr=MPI_Barrier(MPI_COMM_WORLD);

    //compute stencils for interpolated load transfer
    computeLoadTransferStencils();
    ierr=MPI_Barrier(MPI_COMM_WORLD);

    //memory
    aux_nodes0_ = vector<double>();

    //identify nodes needing extra treatment
    if (auxCornerTreatment_) identifyAuxCornerTreatmentNodes();

#if (ASSISTED_MAPPING == 1)
    Screen::MasterInfo("Overwriting aux_nodes0_ from assisted mapping");
    vector<int> dummy;
    vector<int> dummy2;
    surfaceTriangulation surfTri;
    int node = dimLoc_; 
    int dummyint1;
    int dummyint2;
    aux_nodes_.resize(0);
    surfTri.readTriangulation_AUX(aux_nodes0_,dummy,auxFilename_,dummyint1,dummyint2,node,dummy2);
    dummy  = vector<int>();
    dummy2 = vector<int>();
    if (!do_restart_) aux_nodes_.assign(aux_nodes0_.begin(),aux_nodes0_.end());
#endif

    //delay allocation of the copy until memory intensive routines have done their work
    aux_nodes_copy_.resize(numGeomNodes_*xyz,0.);

    //initialize nodal velocities on the geometry
    if (do_restart_ && isFSI_) initializeAuxillaryBodyVelocity();

  }
  
  void fem::identifyAuxillaryParticipants()
  {
    Screen::MasterInfo("Identifying FEM elements on each process that are participating in auxillary transfers");
    //identify FEM elements that will be used in auxillary routines
    // --- triangles only right now --- 

    //store participating elements
    myNumberOfAuxParticipatingElements_ = 0;
    for (int e=0;e<myElementCount_;++e)
      {
        if (elementData_[e].type == MITC3)
          {
            myAuxParticipatingElements_.push_back(e);
            myNumberOfAuxParticipatingElements_ += 1;
          }
      }

    //store participating nodes
    myNumberOfAuxParticipatingNodes_ = 0;
    int offset = 0;
    for (int e=0;e<myElementCount_;++e)
      {
        int nnodes_loc = elementData_[e].nnodes;
        for (int pe=0;pe<myNumberOfAuxParticipatingElements_;++pe)
          {
            int participantID = myAuxParticipatingElements_[pe];

            if (e == participantID)
              {
                for (int n=0;n<nnodes_loc;++n)
                  {
                    int index = n + offset;
                    int node  = faces_[index];
#if (DEDICATED_FEM_PROC == 1)
                    if ((find(myAuxParticipatingNodes_.begin(),myAuxParticipatingNodes_.end(),node) == myAuxParticipatingNodes_.end()) && (node < nof_effective_nodes_))
#else
                    if (find(myAuxParticipatingNodes_.begin(),myAuxParticipatingNodes_.end(),node) == myAuxParticipatingNodes_.end())
#endif
                      {
                        myAuxParticipatingNodes_.push_back(node);
                        myNumberOfAuxParticipatingNodes_ += 1;                        
                      }
                  }
              }
          }
        offset += nnodes_loc;
      }

    centroidForce_.resize(myNumberOfAuxParticipatingElements_,vector<double>(xyz,0.));

  }
  
  void fem::obtainUniqueAuxillaryMapping()
  {

    //push back on moving nodes that are also intersecting FEM elements in my process

#if (ASSISTED_MAPPING == 1) 
    double distance = 0.47249;
#else
    double distance = geomThickFactor_;
#endif

    vector<int> foundNodes(numGeomNodes_,FOUND_NODES_PENO);
#if (RAY_INTERSECTION == 0)

    Screen::MasterInfo("Obtaining unique geometry node - FEM element mapping for each process via distance");

    vector<int> ordered         (numGeomNodes_,0);
    vector<vector<double> > diff(numGeomNodes_,vector<double>(xyz,99.));
    vector<double> dist         (numGeomNodes_,99.);

    for (int fb=-1;fb<=1;fb+=2)
      {
        for (int fn=0;fn<myNumberOfAuxParticipatingNodes_;++fn)
          {
            int participantNodeID = myAuxParticipatingNodes_[fn];//local FEM node ID
            int vectorNode = myLocalNodesNormalVectorStatus_[participantNodeID];

            for (int mn=0;mn<numGeomNodes_;++mn)
              {
                if ( abs(is_aux_node_[mn]) == 99 )
                  {
                    dist[mn] = 0.;
                    for (int idir=0;idir<xyz;++idir)
                      {
                        int vecIndex   = mn*xyz + idir;
                        diff[mn][idir] = aux_nodes0_[vecIndex] - (nodes0_[participantNodeID*xyz+idir] + float(fb)*0.5*distance*defaultThickness_*Vn0_[vectorNode*xyz+idir]);
                        //diff[mn][idir] = aux_nodes0_[vecIndex] - nodes0_[participantNodeID*xyz+idir];
                        dist[mn]      += diff[mn][idir]*diff[mn][idir];
                      }
                  }
              }
            //order geometry nodes closest to farthest from FEM node fn
            getDistanceOrder(ordered,dist);

            for (int i=0;i<1;++i)
              {
                int node = ordered[i];
                foundNodes[node] = mypeno_;
              }

          }
      }

#else

    Screen::MasterInfo("Obtaining unique geometry node - FEM element mapping for each process using ray-triangle intersection");

    vector<double> normal(xyz,0.);   //normal vector at geometry element centroid
    vector<double> gedge1(xyz,0.);   //geometry element edge 1
    vector<double> gedge2(xyz,0.);   //geometry element edge 2
    vector<double> fedge1(xyz,0.);   //FEM element edge 1
    vector<double> fedge2(xyz,0.);   //FEM element edge 2
    vector<double> origin(xyz,0.);   //origin of ray
    vector<double> vertex(xyz,0.);   //reference node on FEM element
    vector<double> centroid(xyz,0.); //geometry element centroid
    bool intersectedElement;

    for (int af=0;af<numGeomFaces_;++af)
      {
        //get geometry element face nodes
        int gnode1 = aux_faces_[af*xyz+0];
        int gnode2 = aux_faces_[af*xyz+1];
        int gnode3 = aux_faces_[af*xyz+2];
        //set director --- does not matter if CULL is on
        double towardsFEM = 1.;
        if (revertGeometryNormals_) towardsFEM = -1.;
        //find geometry element edges
        for (int idir=0;idir<xyz;++idir)
          {
            gedge1  [idir] = towardsFEM  * (aux_nodes0_[gnode2*xyz+idir] - aux_nodes0_[gnode1*xyz+idir]);
            gedge2  [idir] =               (aux_nodes0_[gnode3*xyz+idir] - aux_nodes0_[gnode1*xyz+idir]);
            centroid[idir] = 0.333333333 * (aux_nodes0_[gnode1*xyz+idir] + aux_nodes0_[gnode2*xyz+idir] + aux_nodes0_[gnode3*xyz+idir]);
          }
        //compute normal vector pointing towards FEM mesh
        normal = cross(gedge2,gedge1);
        normalizeVector(normal);
        //loop over geometry element face nodes
        for (int fn=0;fn<3;++fn)
          {
            int node = aux_faces_[af*xyz+fn];
            if ( (abs(is_aux_node_[node]) == 99) && (foundNodes[node] == FOUND_NODES_PENO) )
              {

                double sign = -1.;
                intersectedElement = false;
                //node perturbation loop
                for (int perturb=0;perturb<3;++perturb)
                  {
                    //fill centroid-shifted ray origin
                    if (perturb == 2)
                      {
                        for (int idir=0;idir<xyz;++idir) origin[idir] = centroid[idir];
                      }
                    else
                      {
                        for (int idir=0;idir<xyz;++idir) origin[idir] = aux_nodes0_[node*xyz+idir] + sign*(1.E-3)*(aux_nodes0_[node*xyz+idir]-centroid[idir]);
                      }
                    //test intersection of moving node and an FEM element in my process
                    for (int el=0;el<myNumberOfAuxParticipatingElements_;++el)
                      {
                        int participantID    = myAuxParticipatingElements_[el];
                        int nnodes_loc       = elementData_[participantID].nnodes;
                        int nnodesOffset_loc = nnodesLocalOffset_[participantID];
                        //get FEM element face nodes
                        int fnode1 = faces_[nnodesOffset_loc + 0];
                        int fnode2 = faces_[nnodesOffset_loc + 1];
                        int fnode3 = faces_[nnodesOffset_loc + 2];

                        for (int idir=0;idir<xyz;++idir)
                          {
                            //find FEM element edges
                            fedge1[idir] = nodes0_[fnode2*xyz+idir] - nodes0_[fnode1*xyz+idir];
                            fedge2[idir] = nodes0_[fnode3*xyz+idir] - nodes0_[fnode1*xyz+idir];
                            //fill vertex
                            vertex[idir] = nodes0_[fnode1*xyz+idir];

                          }
                        //test intersection
                        intersectedElement = MollerTrumboreIntersection(origin,normal,fedge1,fedge2,vertex);
                        //store node if it intersects the FEM element
                        if (intersectedElement)
                          {
                            foundNodes[node] = mypeno_;
                            perturb = 3;//so we exit the perturbation loop
                            break;
                          }
                      }
                    sign += 2.;
                  }
              }
          }
      }
#endif

    Screen::MasterInfo("Checking completeness of mapping in parallel");
    ierr = MPI_Allreduce(MPI_IN_PLACE,&foundNodes[0],numGeomNodes_,MPI::INT,MPI_MAX,MPI_COMM_WORLD);

#if (FOUND_NODES_PENO == -1)
    for (int n=0;n<numGeomNodes_;++n)
      {
        if ( (foundNodes[n] == FOUND_NODES_PENO) && (abs(is_aux_node_[n]) == 99) )
          {
            cout << "mypeno  " << mypeno_ << " | geometry node  " << n << " with (x,y,z) = (" << aux_nodes0_[n*xyz+0] << "," << aux_nodes0_[n*xyz+1] << "," << aux_nodes0_[n*xyz+2] << ")" << endl;
            Screen::MasterError("Mapping is not complete");
          }
      }
#endif

    Screen::MasterInfo("Distrubuting nodes uniquely to each process");
    for (int n=0;n<numGeomNodes_;++n)
      {
        if (foundNodes[n] == mypeno_)
          {
            movingNodes_.push_back(n);
          }
      }

    myAuxNodeCount_ = movingNodes_.size();
    Screen::MasterInfo("Gathering geometry partition");

    //gather number of aux nodes on each process
    int myAuxNodeNum = movingNodes_.size();    
    myAuxNodeNumVec_   .resize(numprocs_,0);
    myAuxNodeNumVecXYZ_.resize(numprocs_,0);
    MPI_Allgather(&myAuxNodeNum,1,MPI_INT,&myAuxNodeNumVec_[0],1,MPI_INT,MPI_COMM_WORLD);
    for (int p=0;p<numprocs_;++p) myAuxNodeNumVecXYZ_[p] = myAuxNodeNumVec_[p]*xyz;

    //compute offsets
    nodeIndexOffset_.resize(numprocs_);
    nodeXYZOffset_  .resize(numprocs_);
    nodeIndexOffset_[0] = 0;
    nodeXYZOffset_  [0] = 0;

    for (int p=1;p<numprocs_;++p)
      {
        nodeIndexOffset_[p] = nodeIndexOffset_[p-1] + myAuxNodeNumVec_[p-1];
        nodeXYZOffset_  [p] = nodeIndexOffset_[p]*xyz;
      }

    Screen::MasterInfo("Partitioning geometry faces");
    vector<int> foundFaces(numGeomFaces_, -1);
    for (int n = 0; n < myAuxNodeCount_; ++n)
      {
	int node = movingNodes_[n];
	for (int f = 0; f < numGeomFaces_; ++f)
	  {
	    int nnodes_offset = xyz * f;
	    if (aux_faces_[0 + nnodes_offset] == node || aux_faces_[1 + nnodes_offset] == node ||
		aux_faces_[2 + nnodes_offset] == node)
	      {
		foundFaces[f] = mypeno_;
	      }
	  }
      }

    //tie-breaker goes to higher rank process
    ierr = MPI_Allreduce(MPI_IN_PLACE, &foundFaces[0], numGeomFaces_, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    for (int f = 0; f < numGeomFaces_; ++f)
      {
	//check for unidentified 'moving' faces
	bool isMovingFace = false;
	int nnodes_offset = xyz * f;
	int nMovingNodes = 0;
	for (int n=0;n<xyz;++n)
	  {
	    if (is_aux_node_[aux_faces_[n+nnodes_offset]] == 99) nMovingNodes++;
	  }
	if (nMovingNodes == xyz) isMovingFace = true;
	if (foundFaces[f] == -1 && isMovingFace)
	  {
	    cout << f << " | " << 0+nnodes_offset << " , " << 1+nnodes_offset << " , " << 2+nnodes_offset
		 << " | " << is_aux_node_[0+nnodes_offset] << " , " << is_aux_node_[1+nnodes_offset]
		 << " , " << is_aux_node_[2+nnodes_offset] << endl;
	    Screen::MasterError("Geometry face mapping is not complete");
	  }
      }

    Screen::MasterInfo("Distrubuting faces uniquely to each process");
    for (int f = 0; f < numGeomFaces_; ++f)
      {
	if (foundFaces[f] == mypeno_)
	  {
	    myAuxFaces_.push_back(f);
	  }
      }

#if (DEDICATED_FEM_PROC==0)
#if (DUMP_GEOMETRY_PARTITION == 1)
    Screen::MasterInfo("Dumping geometry partition");
    //dump geometry partition
    dumpGeometryPartition(foundNodes,foundFaces);
#endif
#endif

  }
  
  bool fem::MollerTrumboreIntersection(vector<double>& origin,
                                       vector<double>& normal,
                                       vector<double>& edge1,
                                       vector<double>& edge2,
                                       vector<double>& vertex)
  {

    //test intersection of a ray and a triangle    
#define CULL 1

#if (CULL == 0)
    //only consider 'front facing' triangles

    //initialization
    vector<double> h(xyz,0.);
    vector<double> s(xyz,0.);
    vector<double> q(xyz,0.);
    double a,f,u,v,t;
    double epsilon = 1.E-10;

    h = cross(normal,edge2);
    a = dprod(edge1,h);

    if (a < epsilon) return false;//parallel to triangle
    for (int idir=0;idir<xyz;++idir) s[idir] = origin[idir] - vertex[idir];
    u = dprod(s,h);

    if (u < -epsilon || u > a+epsilon) return false;

    q = cross(s,edge1);
    v = dprod(normal,q);

    if (v < -epsilon || u + v > a+epsilon) return false;

    t = dprod(edge2,q);

    if (t > epsilon)
      {
        if (t > defaultThickness_*geomThickFactor_)
          {
            return false;//intersected an element but it is too far away
          }
        else
          {
            return true;//realistic intersection
          }
      }
    else
      {
        return false;
      }
#else
    //consider two-sided triangles

    //initialization
    vector<double> h(xyz,0.);
    vector<double> s(xyz,0.);
    vector<double> q(xyz,0.);
    double a,f,u,v,t;
    double epsilon = 1.E-8;

    h = cross(normal,edge2);
    a = dprod(edge1,h);

    if (a > -epsilon && a < epsilon) return false;//parallel to triangle

    f = 1./a;
    for (int idir=0;idir<xyz;++idir) s[idir] = origin[idir] - vertex[idir];
    u = f * dprod(s,h);

    if (u < -epsilon || u > 1.+epsilon) return false;

    q = cross(s,edge1);
    v = f * dprod(normal,q);

    if (v < -epsilon || u + v > 1.+epsilon) return false;

    t = f * dprod(edge2,q);

    if (t > epsilon)
      {
        if (t > 0.1)//for wiggle room on parachute
          {
            return false;//intersected an element but it is too far away
          }
        else
          {
            return true;//realistic intersection
          }
      }
    else
      {
        return false;
      }
#endif

  }
  
  void fem::treatSingleProcFSI()
  {

    Screen::MasterInfo("Treating FSI with a dedicated FEM proc");
    //trick parallel layout for a single proc FSI simulation

    for (int f=0;f<numGeomFaces_;++f)
      {
        for (int n=0;n<dimLoc_;++n)
          {
            int index = f*dimLoc_ + n;
            int node  = aux_faces_[index];

            if (abs(is_aux_node_[node]) == 99)
              {
                if ( find(movingNodes_.begin(),movingNodes_.end(),node) == movingNodes_.end() )
                  {
                    movingNodes_.push_back(node);
                  }
              }
          }
	myAuxFaces_.push_back(f);
      }

    myAuxNodeCount_     = movingNodes_.size();

    myAuxNodeNumVec_   .resize(numprocs_,-1);
    myAuxNodeNumVecXYZ_.resize(numprocs_,-1);
    nodeIndexOffset_   .resize(numprocs_,-1);
    nodeXYZOffset_     .resize(numprocs_,-1);

    myAuxNodeNumVec_   [0] = myAuxNodeCount_;
    myAuxNodeNumVecXYZ_[0] = myAuxNodeCount_*xyz;
    nodeIndexOffset_   [0] = 0;
    nodeXYZOffset_     [0] = 0;

  }
  
  void fem::restructureFSI()
  {

    Screen::MasterInfo("Restructuring geometry variables");
    //order geometry nodes contiguously across processes

    myAuxNodes_.resize(myAuxNodeCount_*xyz,0.);
    for (int gn=0;gn<myAuxNodeCount_;++gn)
      { 
        //save my portion of the total geometry
        for (int idir=0;idir<xyz;++idir)
          {
            myAuxNodes_[gn*xyz+idir]  = aux_nodes0_[movingNodes_[gn]*xyz+idir];
          }
      }
    //create mapping to unscramble aux nodes after solve
#if (DEDICATED_FEM_PROC == 0)

    geometryMap_.resize(numGeomNodes_,-1);
    ierr = MPI_Allgatherv(&movingNodes_[0],myAuxNodeNumVec_[mypeno_],MPI::INT,&geometryMap_[0],&myAuxNodeNumVec_[0],&nodeIndexOffset_[0],MPI::INT,MPI_COMM_WORLD);
#else
    geometryMap_.assign(movingNodes_.begin(),movingNodes_.end());
#endif

  }
  
  void fem::initializeDisplacementTransfer()
  {

    Screen::MasterInfo("Identifying stencils for displacement transfer");
    //storage for displacement transfer

    //initialize storage
    fem2geom.resize(myAuxNodeCount_);

    //store k closest fem nodes to geometry node 'gn'

    //calculate distance between geometry nodes in my process and all fem nodes
    for (int mn=0;mn<myAuxNodeCount_;++mn)
      {

#if (DEDICATED_FEM_PROC == 1)

        int gn = mn;//movingNodes_[mn];

        vector<int> ordered         (nof_effective_nodes_,0);
        vector<double> dist         (nof_effective_nodes_,0.);
        vector<vector<double> > diff(nof_effective_nodes_,vector<double>(xyz,0.));

        for (int fn=0;fn<nof_effective_nodes_;++fn)
          {
            int participantNodeID = myAuxParticipatingNodes_[fn];
            int vectorNode = myLocalNodesNormalVectorStatus_[participantNodeID];
#else

        int gn = mn;

        vector<int> ordered         (myNumberOfAuxParticipatingNodes_,0);
        vector<double> dist         (myNumberOfAuxParticipatingNodes_,99.);
        vector<vector<double> > diff(myNumberOfAuxParticipatingNodes_,vector<double>(xyz,99.));

        for (int fn=0;fn<myNumberOfAuxParticipatingNodes_;++fn)
          {
            dist[fn] = 0.;
            int participantNodeID = myAuxParticipatingNodes_[fn];
            int vectorNode = myLocalNodesNormalVectorStatus_[participantNodeID];
#endif
            for (int idir=0;idir<xyz;++idir)
              {
                int vecIndex   = gn*xyz + idir;
                diff[fn][idir] = myAuxNodes_[vecIndex] - nodes0_[participantNodeID*xyz+idir];
                dist[fn]      += diff[fn][idir]*diff[fn][idir];		
              }
          }

        if (myNumberOfAuxParticipatingNodes_ > 0)
          {
            //order fem nodes closest to farthest from geom node 'gn' and store k closest	    
            getDistanceOrder(ordered,dist);
            for (int fn=0;fn<leastSQR_npnts_;++fn)
              {
                int closestParticipatingFEMNodeLocalID = ordered[fn];
                int FEMNodeLocalID = myAuxParticipatingNodes_[closestParticipatingFEMNodeLocalID];
                fem2geom[gn].closestFEMnode.push_back(FEMNodeLocalID);
              }
          }
        else
          {
            fem2geom[gn].closestFEMnode.push_back(0);
          }

#if (OUTPUT_DISPLACEMENT_POINT_CLOUD == 1)
        for (int p=0;p<numprocs_;++p)
          {
            if (mypeno_ == p)
              {
                cout << "mypeno = " <<  mypeno_ << endl;
                cout << "geometry node " << gn << "(x,y,z) = (" << myAuxNodes_[gn*xyz]<<","<<myAuxNodes_[gn*xyz+1]<<","<<myAuxNodes_[gn*xyz+2]<<")" << endl;
                for (int fn=0;fn<leastSQR_npnts_;++fn)
                  {
                    cout << " identified FEM node " <<  fem2geom[gn].closestFEMnode[fn] << " (x,y,z) = (" << nodes0_[fem2geom[gn].closestFEMnode[fn]*xyz+0]<<","<<nodes0_[fem2geom[gn].closestFEMnode[fn]*xyz+1]<<","<<nodes0_[fem2geom[gn].closestFEMnode[fn]*xyz+2]<<")"<<endl;
                  }
              }
            //ierr = MPI_Barrier(MPI_COMM_WORLD);
            //wait();
          }
#endif
      }

   }   
    
  void fem::initializeLoadTransfer()
   {

     Screen::MasterInfo("Identifying stencils for load transfer");
     //storage for load transfer method with single point, one-to-one mesh with a 'single slice' of geometry laying on top on FEM
     //store fem node closest to geom node 'gn'

     //initialize storage    
     geom2fem.resize(myNumberOfAuxParticipatingNodes_);

     if (auxForceInterpMethod_ == "single")
       {
         for (int fn=0;fn<myNumberOfAuxParticipatingNodes_;++fn)
           {
             geom2fem[fn].closestGeomNode.push_back(vector<int>());
             int participantNodeID = myAuxParticipatingNodes_[fn];//local FEM node ID
             int vectorNode = myLocalNodesNormalVectorStatus_[participantNodeID];

             vector<int> ordered(movingNodes_.size(),0);
             vector<vector<double> > diff(movingNodes_.size(),vector<double>(xyz,99.));
             vector<double> dist(movingNodes_.size(),99.);

             for (int mn=0;mn<movingNodes_.size();++mn)
               {
#if (DEDICATED_FEM_PROC == 0)
                 int gn = mn;
#else
                 int gn = movingNodes_[mn];
#endif
                 for (int idir=0;idir<xyz;++idir)
                   {
                     int vecIndex   = gn*xyz + idir;
                     diff[mn][idir] = myAuxNodes_[vecIndex] - (nodes0_[participantNodeID*xyz+idir] + 0.5*geomThickFactor_*defaultThickness_*Vn0_[vectorNode*xyz+idir]);
                     dist[mn]      += diff[mn][idir]*diff[mn][idir];
                   }
               }
             //dummy, do nothing
             if (movingNodes_.size() == 0)
               {
                 for (int k=0;k<numForceInterpPts_;++k) geom2fem[fn].closestGeomNode[0].push_back(0);
               }
             else
               {
                 //order fem nodes closest to farthest from geom node 'gn'
                 getDistanceOrder(ordered,dist);
                 for (int k=0;k<numForceInterpPts_;++k) geom2fem[fn].closestGeomNode[0].push_back(movingNodes_[ordered[k]]);
               }

#if (OUTPUT_FORCE_POINT_CLOUD == 1)
             cout << "fem node " << fn << " (x,y,z) = (" << nodes0_[fn*xyz+0] << "," << nodes0_[fn*xyz+1] << "," << nodes0_[fn*xyz+2] << ")" << endl;
             cout << "is taking loading from geometry nodes: " << endl;
             cout << geom2fem[fn].closestGeomNode[0][0] << " (x,y,z) = (" << aux_nodes0_[geom2fem[fn].closestGeomNode[0][0]*xyz]<<","<<aux_nodes0_[geom2fem[fn].closestGeomNode[0][0]*xyz+1]<<","<<aux_nodes0_[geom2fem[fn].closestGeomNode[0][0]*xyz+2]<<")"<<endl;
             wait();
#endif
           }
       }
     else if (auxForceInterpMethod_ == "double")
       {
         //initialize storage

#if (ASSISTED_MAPPING == 1) 
         double distance = 0.47249;
#else
         double distance = geomThickFactor_;
#endif

         int faceID = 0;
         for (int fb=-1;fb<=1;fb+=2)
           {
             for (int fn=0;fn<myNumberOfAuxParticipatingNodes_;++fn)
               {
                 geom2fem[fn].closestGeomNode.push_back(vector<int>());
                 int participantNodeID = myAuxParticipatingNodes_[fn];//local FEM node ID
                 int vectorNode = myLocalNodesNormalVectorStatus_[participantNodeID];

                 vector<int> ordered         (numGeomNodes_,0);
                 vector<vector<double> > diff(numGeomNodes_,vector<double>(xyz,99.));
                 vector<double> dist         (numGeomNodes_,99.);

                 for (int mn=0;mn<numGeomNodes_;++mn)
                   {
                     if (is_aux_node_[mn] == 99)
                       {
                         dist[mn] = 0.;
                         for (int idir=0;idir<xyz;++idir)
                           {
                             int vecIndex   = mn*xyz + idir;
                             diff[mn][idir] = aux_nodes0_[vecIndex] - (nodes0_[participantNodeID*xyz+idir] + float(fb)*0.5*distance*defaultThickness_*Vn0_[vectorNode*xyz+idir]);
                             dist[mn]      += diff[mn][idir]*diff[mn][idir];
                           }
                       }
                   }
                 //order geometry nodes closest to farthest from FEM node fn
                 getDistanceOrder(ordered,dist);

                 int globalGeomNodeID = ordered[0];
                 geom2fem[fn].closestGeomNode[faceID].push_back(globalGeomNodeID);

#if (OUTPUT_FORCE_POINT_CLOUD == 1)
                 cout << "mypeno " << mypeno_ << " | fem node " << fn << " (x,y,z) = (" << nodes0_[fn*xyz+0] << "," << nodes0_[fn*xyz+1] << "," << nodes0_[fn*xyz+2] << ") found the following cloud on the " << faceID << " face: " << endl;
                 for (int k=0;k<numForceInterpPts_;++k)
                   {
                     cout << geom2fem[fn].closestGeomNode[faceID][k] << " (x,y,z) = (" << aux_nodes0_[geom2fem[fn].closestGeomNode[faceID][k]*xyz]<<","<<aux_nodes0_[geom2fem[fn].closestGeomNode[faceID][k]*xyz+1]<<","<<aux_nodes0_[geom2fem[fn].closestGeomNode[faceID][k]*xyz+2]<<")"<<endl;
                   }
                 wait();
#endif
               }
             faceID += 1;
           }
       }
     else
       {
         //storage for load transfer method with point clouds on 'plus-minus' geometry faces, via compID or normal vector
         int faceID = 0;
         int count  = 0;
         int ID     = -99;

         for (int fb=-1;fb<=1;fb+=2)
           {
             for (int fn=0;fn<myNumberOfAuxParticipatingNodes_;++fn)
               {
                 int participantNodeID = myAuxParticipatingNodes_[fn];//local FEM node ID
		 int vectorNode = myLocalNodesNormalVectorStatus_[participantNodeID];
                 geom2fem[fn].closestGeomNode.push_back(vector<int>());
                 vector<double> xyzvec(xyz,0.);
                 vector<double> dist         (numGeomNodes_,99.);
                 vector<int> ordered         (numGeomNodes_,0);
                 vector<vector<double> > diff(numGeomNodes_,vector<double>(xyz,99.));

                 if (loadsByCompID_ == true)
                   {
                     for (int mn=0;mn<myAuxNodeCount_;++mn)
                       {
#if (DEDICATED_FEM_PROC == 1)
                         int gn = movingNodes_[mn];
#else
                         int gn = mn;
#endif

                         dist[mn] = 0.;
                         for (int idir=0;idir<xyz;++idir)
                           {
                             int vecIndex   = gn*xyz + idir;
                             diff[mn][idir] = aux_nodes0_[vecIndex] - (nodes0_[participantNodeID*xyz+idir] + float(fb)*0.5*geomThickFactor_*defaultThickness_*Vn0_[vectorNode*xyz+idir]);
                             dist[mn]      += diff[mn][idir]*diff[mn][idir];
                           }                        
                       }
                   }
                 else
                   {
                     for (int mn=0;mn<numGeomNodes_;++mn)
                       {
			 if (is_aux_node_[mn] == 99)
			   {
			     dist[mn] = 0.;
#if (DEDICATED_FEM_PROC == 1)
			     int gn = movingNodes_[mn];
#else
			     int gn = mn;
#endif
			     for (int idir=0;idir<xyz;++idir)
			       {
				 int vecIndex   = gn*xyz + idir;
				 diff[mn][idir] = aux_nodes0_[vecIndex] - (nodes0_[participantNodeID*xyz+idir] + float(fb)*0.5*rayLength_*defaultThickness_*Vn0_[vectorNode*xyz+idir]);
				 dist[mn]      += diff[mn][idir]*diff[mn][idir];
			       }
			   }
		       }
		     //order fem nodes closest to farthest from geom node 'gn'
		     getDistanceOrder(ordered,dist);

#if (DEDICATED_FEM_PROC == 1)
                     for (int k=0;k<numForceInterpPts_;++k) geom2fem[fn].closestGeomNode[faceID].push_back(movingNodes_[ordered[k]]);
#else
                     for (int k=0;k<numForceInterpPts_;++k)
                       {
                         int node = ordered[k];
                         geom2fem[fn].closestGeomNode[faceID].push_back(node);

                       }
#endif

                   }

#if (OUTPUT_FORCE_POINT_CLOUD == 1)
                 for (int p=0;p<numprocs_;++p)
                   {
                     if (mypeno_ == p)
                       {
                         cout << "mypeno = " <<  mypeno_ << endl;
                         cout << "fem node local ID " << participantNodeID << " (x,y,z) = (" << nodes0_[participantNodeID*xyz+0] << "," << nodes0_[participantNodeID*xyz+1] << "," << nodes0_[participantNodeID*xyz+2] << ") found the following cloud on the " << faceID << " face: " << endl;
                         for (int k=0;k<numForceInterpPts_;++k)
                           {
                             int gn = geom2fem[fn].closestGeomNode[faceID][k];
                             cout << "geom node: " << gn << " (x,y,z) = (" << aux_nodes0_[gn*xyz+0] << "," << aux_nodes0_[gn*xyz+1] << "," << aux_nodes0_[gn*xyz+2] << ")" << endl;
                           }                        
                       }
                     //ierr = MPI_Barrier(MPI_COMM_WORLD);
                   }
                 //wait();
#endif
               }
             ID = 99;
             faceID += 1;
           }
       }

     auxGeomForces_.resize(5*numGeomNodes_,0.);

   }
  
  void fem::computeDisplacementTransferStencils()
  {
    Screen::MasterInfo("Computing stencil coefficients for auxilliary displacement transfer");
    auxStencil sten;
    sten.nodes.resize(leastSQR_npnts_,-1);
    sten.coefs.resize(leastSQR_npnts_,1.0E20);
#if (DEDICATED_FEM_PROC == 1)
    stencils_.resize(numGeomNodes_,sten);
#else
    stencils_.resize(myAuxNodeCount_,sten);
#endif

    int numEqs = 0;
    //for general dx,dy,dz
    if (auxDispInterpMethod_ == "LSQR") 
      {
        if (dispInterpOrder_ == "first")
          {
            numEqs = 3;
          }
        else if (dispInterpOrder_ == "second")
          {
            numEqs = 6;
          }
        else
          {
            Screen::MasterError("dispInterpOrder_ must be 'first' or 'second' :: fem.cc");
          }
      }
    //for directional dependent interpolation, i.e., leaving out x,z
    else
      {
        if (dispInterpOrder_ == "first")
          {
            numEqs = 3;
          }
        else if (dispInterpOrder_ == "second")
          {
            numEqs = 6;
          }
        else
          {
            Screen::MasterError("dispInterpOrder_ must be 'first' or 'second' :: fem.cc");
          }
      }

    vector<vector<double> > tempm(leastSQR_npnts_,vector<double>(numEqs,0.));
    vector<vector<double> > tempmT(numEqs,vector<double>(leastSQR_npnts_,0.));

    for (int mn=0;mn<movingNodes_.size();++mn)
      {
#if (DEDICATED_FEM_PROC == 1)
        int gn = mn;//movingNodes_[mn]; cylfil case for JCP
#else
        int gn = mn;
#endif
        if (auxDispInterpMethod_ == "single")
          {
            //store FEM node directly opposite geometry
            int index = fem2geom[gn].closestFEMnode[0];
            stencils_[gn].nodes[0] = index;
            stencils_[gn].coefs[0] = 1.0;

            if (myNumberOfAuxParticipatingNodes_ > 0)
              {		   
                //compute direction geometry node traveled to FEM nodes
                vector<double> xyzvec(xyz,0.);
                stencils_[gn].sign = 1.;
                double product = 0.;               
                int vectorNode = myLocalNodesNormalVectorStatus_[index];
                for (int idir=0;idir<xyz;++idir)
                  {
                    xyzvec[idir] = myAuxNodes_[gn*xyz+idir] - nodes0_[index*xyz+idir];
                    product     += Vn0_[vectorNode*xyz+idir]*xyzvec[idir];
                  }
                stencils_[gn].sign = checksign(product);
              }
            else
              {
                stencils_[gn].sign = 1.;
              }
          }

        if (auxDispInterpMethod_ == "LSQR")
          {
            //default sign to positive
            stencils_[gn].sign = 1.;
            //normal vector shoots from closest FEM node
            int vectorNode = myLocalNodesNormalVectorStatus_[fem2geom[gn].closestFEMnode[0]];
            vector<double> normal(3,0.);
            normal[0] = Vn0_[vectorNode*xyz + 0];
            normal[1] = Vn0_[vectorNode*xyz + 1];
            normal[2] = Vn0_[vectorNode*xyz + 2];
            normalizeVector(normal);

            //basis vector allocation
            vector<double> svec(3,0.);
            vector<double> tvec(3,0.);

            //find first suitable basis vector
            for (int k=0;k<leastSQR_npnts_-1;++k)
              {
                tvec[0] = nodes0_[fem2geom[gn].closestFEMnode[k+1]*xyz + 0] - nodes0_[fem2geom[gn].closestFEMnode[k]*xyz + 0];
                tvec[1] = nodes0_[fem2geom[gn].closestFEMnode[k+1]*xyz + 1] - nodes0_[fem2geom[gn].closestFEMnode[k]*xyz + 1];
                tvec[2] = nodes0_[fem2geom[gn].closestFEMnode[k+1]*xyz + 2] - nodes0_[fem2geom[gn].closestFEMnode[k]*xyz + 2];

                //check alignment with normal vector

                normalizeVector(tvec);
                double product = dprod(normal,tvec);
                if ( abs(1. - product) > 1.E-4) break;//if not aligned, good to go

              }

            //find second suitable basis vector

            svec = cross(normal,tvec);

            //re-compute orthogonal first basis vector

            //tvec = cross(svec,normal);

            //orthonormal basis vectors
            normalizeVector(tvec);
            normalizeVector(svec);

            //save geometry node position in global Cartesian basis
            vector<double> geomtryNodeXYZ(3,0.);
            geomtryNodeXYZ[0] = myAuxNodes_[gn*xyz + 0];
            geomtryNodeXYZ[1] = myAuxNodes_[gn*xyz + 1];
            geomtryNodeXYZ[2] = myAuxNodes_[gn*xyz + 2];

            //transform geometry node into new basis as basis origin
            vector<double> origin(2,0.);
            origin[0] = dprod(geomtryNodeXYZ,svec);
            origin[1] = dprod(geomtryNodeXYZ,tvec);

            //global Cartesian coordinates of the FEM nodes
            vector<double> xyzVec(3,0.);            

            //transformed coordinates of the FEM nodes
            vector<double> xyzTilda(2,0.);

            //pseudo-2D stencil generation in the computed basis
            for (int k=0;k<leastSQR_npnts_;++k)
              {
                //get closest FEM node local ID
                int index = fem2geom[gn].closestFEMnode[k];
                stencils_[gn].nodes[k] = index;

                xyzVec[0] = nodes0_[index*xyz + 0];
                xyzVec[1] = nodes0_[index*xyz + 1];
                xyzVec[2] = nodes0_[index*xyz + 2];

                xyzTilda[0] = dprod(xyzVec,svec) - origin[0];
                xyzTilda[1] = dprod(xyzVec,tvec) - origin[1];

                //assemble coefficient matrix
                tempm[k][0] = 1.;
                tempm[k][1] = xyzTilda[0];
                tempm[k][2] = xyzTilda[1];
                //second order interp.
                if (dispInterpOrder_ == "second")
                  {
                    tempm[k][3] = xyzTilda[0]*xyzTilda[0];
                    tempm[k][4] = xyzTilda[1]*xyzTilda[1];
                    tempm[k][5] = xyzTilda[0]*xyzTilda[1];
                  }
              }
            //take Penrose inverse of coefficient matrix
            getLSQRcoefs(tempm,0,stencils_[gn].coefs);

            //determine which geometry face we are shooting to (plus or minus)
            vector<double> nvec(xyz,0.);
            //compute if xyz lies with or against Vn0
            for (int p=0;p<leastSQR_npnts_;++p)
              {
                int vectorNode = myLocalNodesNormalVectorStatus_[stencils_[gn].nodes[p]];
                for (int idir=0;idir<xyz;++idir)
                  {
                    nvec[idir] += stencils_[gn].coefs[p]*Vn0_[vectorNode*xyz + idir];
                  }
              }
            //overwrite normal with difference in geometry node and closest FEM node
            normal[0] = myAuxNodes_[gn*xyz + 0] - nodes0_[fem2geom[gn].closestFEMnode[0]*xyz + 0];
            normal[1] = myAuxNodes_[gn*xyz + 1] - nodes0_[fem2geom[gn].closestFEMnode[0]*xyz + 1];
            normal[2] = myAuxNodes_[gn*xyz + 2] - nodes0_[fem2geom[gn].closestFEMnode[0]*xyz + 2];

            double dummy = dprod(normal,nvec);
            stencils_[gn].sign = checksign(dummy);

          }
        else if (auxDispInterpMethod_ == "LSQRZ")
          {
            for (int n=0;n<leastSQR_npnts_;++n)
              {
                int index = fem2geom[gn].closestFEMnode[n];
                stencils_[gn].nodes[n] = index;

                vector<double> xyzvec(xyz,0.);
                for (int idir=0;idir<xyz;++idir)
                  {
                    int vecIndex = gn*xyz + idir;
                    xyzvec[idir] = myAuxNodes_[vecIndex] - nodes0_[index*xyz+idir];
                  }
                //first-order interp
                tempm[n][0] = 1.;
                tempm[n][1] = xyzvec[0];
                tempm[n][2] = xyzvec[1];
                //second order interp
                if (dispInterpOrder_ == "second")
                  {
                    tempm[n][3] = xyzvec[0]*xyzvec[0];
                    tempm[n][4] = xyzvec[0]*xyzvec[1];
                    tempm[n][5] = xyzvec[0]*xyzvec[1];
                  }
              }

            getLSQRcoefs(tempm,0,stencils_[gn].coefs);

          }
        else if (auxDispInterpMethod_ == "WLSQRZ")
          {
            stencils_[gn].sign = 1.;
            vector<vector<double> > xyzmat(leastSQR_npnts_,vector<double>(xyz,0.));
            vector<double> weights(leastSQR_npnts_,1.);
            for (int n=0;n<leastSQR_npnts_;++n)
              {
                int index = fem2geom[gn].closestFEMnode[n];
                stencils_[gn].nodes[n] = index;

                vector<double> xyzvec(xyz,0.);
                for (int idir=0;idir<xyz;++idir)
                  {
                    int vecIndex    = gn*xyz + idir;
                    xyzvec[idir]    = myAuxNodes_[vecIndex] - nodes0_[index*xyz+idir];
                    xyzmat[n][idir] = xyzvec[idir];
                  }

                double mag = vectorMag(xyzvec);
                weights[n] = 1.;

                //first-order interp
                tempm[n][0] = 1.;
                tempm[n][1] = xyzvec[0];
                tempm[n][2] = xyzvec[1];
                //second-order interp
                if (dispInterpOrder_ == "second")
                  {
                    tempm[n][3] = xyzvec[0]*xyzvec[0];
                    tempm[n][4] = xyzvec[1]*xyzvec[1];
                    tempm[n][5] = xyzvec[0]*xyzvec[1];
                  }                
              }

            getWLSQRcoefs(tempm,weights,0,stencils_[gn].coefs);

            vector<double> normal(xyz,0.);
            vector<double> xyzsum(xyz,0.);
            //compute if xyz lies with or against Vnt
            for (int p=0;p<leastSQR_npnts_;++p)
              {
                int vectorNode = myLocalNodesNormalVectorStatus_[stencils_[gn].nodes[p]];
                for (int idir=0;idir<xyz;++idir)
                  {
                    normal[idir] += stencils_[gn].coefs[p]*Vn0_  [vectorNode*xyz + idir];
                    xyzsum[idir] += stencils_[gn].coefs[p]*xyzmat[p][idir];
                  }
              }
            double dummy = dprod(normal,xyzsum);
            stencils_[gn].sign = checksign(dummy);
          }          
        else if (auxDispInterpMethod_ == "WLSQRX")
          {
            stencils_[gn].sign = 1.;
            vector<vector<double> > xyzmat(leastSQR_npnts_,vector<double>(xyz,0.));
            vector<double> weights(leastSQR_npnts_,1.);
            for (int n=0;n<leastSQR_npnts_;++n)
              {
                int index = fem2geom[gn].closestFEMnode[n];
                stencils_[gn].nodes[n] = index;

                vector<double> xyzvec(xyz,0.);
                for (int idir=0;idir<xyz;++idir)
                  {
                    int vecIndex    = gn*xyz + idir;
                    xyzvec[idir]    = myAuxNodes_[vecIndex] - nodes0_[index*xyz+idir];
                    xyzmat[n][idir] = xyzvec[idir];
                  }

                double mag = vectorMag(xyzvec);
                weights[n] = 1.;
                //first-order interp
                tempm[n][0] = 1.;
                tempm[n][1] = xyzvec[1];
                tempm[n][2] = xyzvec[2];
                //second-order interp
                if (dispInterpOrder_ == "second")
                  {
                    tempm[n][3] = xyzvec[1]*xyzvec[1];
                    tempm[n][4] = xyzvec[2]*xyzvec[2];
                    tempm[n][5] = xyzvec[1]*xyzvec[2];
                  }
              }

            getWLSQRcoefs(tempm,weights,0,stencils_[gn].coefs);

            vector<double> normal(xyz,0.);
            vector<double> xyzsum(xyz,0.);
            //compute if xyz lies with or against Vnt
            for (int p=0;p<leastSQR_npnts_;++p)
              {
                int vectorNode = myLocalNodesNormalVectorStatus_[stencils_[gn].nodes[p]];
                for (int idir=0;idir<xyz;++idir)
                  {
                    normal[idir] += stencils_[gn].coefs[p]*Vn0_[vectorNode*xyz + idir];
                    xyzsum[idir] += stencils_[gn].coefs[p]*xyzmat[p][idir];
                  }
              }
            double dummy = dprod(normal,xyzsum);
            stencils_[gn].sign = checksign(dummy);
          }
        //test for correct coefficients
        double sum = -1.;
        for (int i=0;i<leastSQR_npnts_;++i)
          {
            sum += stencils_[gn].coefs[i]; 
            if (stencils_[gn].coefs[i] != stencils_[gn].coefs[i])               
              {
                for (int n=0;n<leastSQR_npnts_;++n)
                  {
                    cout << tempm[n][0] << " "<< tempm[n][1] << " "<< tempm[n][2] << " "<< tempm[n][3] << " "<< tempm[n][4] << " "<< tempm[n][5]<<"  coef = "<<stencils_[gn].coefs[i]<<endl;
                  }
                Screen::MasterError("displacement interpolation stencils did not pass sum test at geom node "+to_string(gn));
              }
          }

        if (abs(sum) > 1.E-3)
          {
            cout << "[I] mypeno =  " << mypeno_ << ":" << endl;
            cout << "[I] at geometry node " << gn << ":" << endl;
            cout << "[I] with position (x,y,z) = (" << aux_nodes0_[gn*xyz] << ", " << aux_nodes0_[gn*xyz + 1] << ", " << aux_nodes0_[gn*xyz + 2] << ")" << endl;
            cout << "[I] Stencil sum was " << sum << endl;
            cout << "[I] Stencil FEM nodes are " << endl;
            for (int n=0;n<leastSQR_npnts_;++n)
              {
                int index = fem2geom[gn].closestFEMnode[n];
                cout << "node " << index << ": (x,y,z) = ";
                for (int idir=0;idir<xyz;++idir)
                  {
                    cout << nodes0_[index*xyz+idir] << "  ";
                  }
                cout << endl;
              }
            Screen::MasterError("Stencil coefs are not correct in auxStencilGen. Try another interpolation method");
          }
      }

  }
 
  void fem::computeLoadTransferStencils()
  {
    Screen::MasterInfo("Computing stencil coefficients for auxilliary load transfer");
    auxStencil sten;
    sten.loadCoefs.resize(2,vector<double>(numForceInterpPts_,0.));
    loadStencils_.resize(myNumberOfAuxParticipatingNodes_,sten);
    //perform interpolation to move pressure at the 'n' closest geometry nodes to FEM node fn
    int faceID = 0;
    int numEqs;

    if (auxForceInterpMethod_ == "LSQR") 
      {
        if (forceInterpOrder_ == "first")
          {
            numEqs = 3;
          }
        else if (forceInterpOrder_ == "second")
          {
            numEqs = 6;
          }
        else
          {
            Screen::MasterError("forceInterpOrder_ must be 'first' or 'second' :: fem.cc");
          }
      }
    else
      {      
        if (forceInterpOrder_ == "first")
          {
            numEqs = 3;
          }
        else if (forceInterpOrder_ == "second")
          {
            numEqs = 6;
          }
        else
          {
            Screen::MasterError("forceInterpOrder_ must be 'first' or 'second' :: fem.cc");
          }
      }

    if (auxForceInterpMethod_ == "single")
      {
        if (myAuxNodeCount_ > 0)
          {
            for (int fn=0;fn<myNumberOfAuxParticipatingNodes_;++fn)
              {
                int gn = geom2fem[fn].closestGeomNode[0][0];
                loadStencils_[fn].loadCoefs[0][0]  = 1.;
              }
          }
      }
    else if (auxForceInterpMethod_ == "double")
      {
        int faceID = 0;
        for (int fb=-1;fb<=1;fb+=2)
          {
            for (int fn=0;fn<myNumberOfAuxParticipatingNodes_;++fn)
              {
                loadStencils_[fn].loadCoefs[faceID][0] = 1.;
              }            
            faceID += 1;
          }
      }
    else if (auxForceInterpMethod_ == "WLSQRZ" || auxForceInterpMethod_ == "WLSQRX" || auxForceInterpMethod_ == "LSQR")
      {
        vector<double> weights       (numForceInterpPts_,1.);
        vector<vector<double> > tempm(numForceInterpPts_,vector<double>(numEqs,0.));

        for (int fb=-1;fb<=1;fb+=2)
          {
            for (int fn=0;fn<myNumberOfAuxParticipatingNodes_;++fn)
              {
                int participantNodeID = myAuxParticipatingNodes_[fn];//local FEM node ID

                if (auxForceInterpMethod_ == "WLSQRZ")
                  {                    
                    for (int n=0;n<numForceInterpPts_;++n)
                      {
                        int gn = geom2fem[fn].closestGeomNode[faceID][n];

                        vector<double> xyzvec(xyz,0.);
                        int vectorNode = myLocalNodesNormalVectorStatus_[participantNodeID];
                        for (int idir=0;idir<xyz;++idir)
                          {
                            int vecIndex = gn*xyz + idir;
                            xyzvec[idir] = (nodes0_[participantNodeID*xyz+idir] + 0.5*float(faceID)*Vn0_[vectorNode*xyz+idir]*defaultThickness_) - aux_nodes0_[vecIndex];
                          }
                        loadStencils_[fn].loadCoefs[faceID][n] = 0.;
                        weights[n]       = 1.;
                        for (int k=0;k<numEqs;++k) tempm[n][k] = 0.;

                        //first-order interp
                        tempm[n][0] = 1.;
                        tempm[n][1] = xyzvec[0];
                        tempm[n][2] = xyzvec[1];
                        //second-order interp
                        if (forceInterpOrder_ == "second")
                          {
                            tempm[n][3] = xyzvec[0]*xyzvec[0];
                            tempm[n][4] = xyzvec[1]*xyzvec[1];
                            tempm[n][5] = xyzvec[0]*xyzvec[1];
                          }
                      }
                  }
                else if (auxForceInterpMethod_ == "WLSQRX")
                  {
                    for (int n=0;n<numForceInterpPts_;++n)
                      {
                        int gn = geom2fem[fn].closestGeomNode[faceID][n];

                        vector<double> xyzvec(xyz,0.);
                        int vectorNode = myLocalNodesNormalVectorStatus_[participantNodeID];
                        for (int idir=0;idir<xyz;++idir)
                          {
                            int vecIndex = gn*xyz + idir;
                            xyzvec[idir] = (nodes0_[participantNodeID*xyz+idir] + 0.5*float(faceID)*Vn0_[vectorNode*xyz+idir]*defaultThickness_) - aux_nodes0_[vecIndex];
                          }
                        loadStencils_[fn].loadCoefs[faceID][n] = 0.;
                        weights[n]       = 1.;
                        for (int k=0;k<numEqs;++k) tempm[n][k] = 0.;

                        //first-order interp
                        tempm[n][0] = 1.;
                        tempm[n][1] = xyzvec[1];
                        tempm[n][2] = xyzvec[2];
                        //second-order interp
                        if (forceInterpOrder_ == "second")
                          {
                            tempm[n][3] = xyzvec[1]*xyzvec[1];
                            tempm[n][4] = xyzvec[2]*xyzvec[2];
                            tempm[n][5] = xyzvec[1]*xyzvec[2];
                          }
                      }
                  }
                else if (auxForceInterpMethod_ == "LSQR")
                  {
                    //normal vector from FEM node
                    int vectorNode = myLocalNodesNormalVectorStatus_[participantNodeID];
                    vector<double> normal(3,0.);
                    normal[0] = Vn0_[vectorNode*xyz + 0];
                    normal[1] = Vn0_[vectorNode*xyz + 1];
                    normal[2] = Vn0_[vectorNode*xyz + 2];

                    //basis vector allocation
                    vector<double> svec(3,0.);
                    vector<double> tvec(3,0.);

                    //find first suitable basis vector
                    for (int k=0;k<numForceInterpPts_;++k)
                      {

                        tvec[0] = aux_nodes0_[geom2fem[fn].closestGeomNode[faceID][k]*xyz + 0];
                        tvec[1] = aux_nodes0_[geom2fem[fn].closestGeomNode[faceID][k]*xyz + 1];
                        tvec[2] = aux_nodes0_[geom2fem[fn].closestGeomNode[faceID][k]*xyz + 2];

                        //check alignment with normal vector

                        normalizeVector(tvec);
                        double product = dprod(normal,tvec);
                        if ( abs(1. - product) > 1.E-4) break;//if not aligned, good to go

                      }

                    //find second suitable basis vector

                    svec = cross(normal,tvec);

                    //re-compute orthogonal first basis vector

                    tvec = cross(svec,normal);

                    //orthonormal basis vectors
                    normalizeVector(tvec);
                    normalizeVector(svec);                    

                    //save FEM node XYZ in global Cartesian basis
                    vector<double> femNodeXYZ(3,0.);
                    femNodeXYZ[0] = nodes0_[participantNodeID*xyz + 0];
                    femNodeXYZ[1] = nodes0_[participantNodeID*xyz + 1];
                    femNodeXYZ[2] = nodes0_[participantNodeID*xyz + 2];

                    //transform FEM node into new basis as basis origin

                    vector<double> origin(2,0.);
                    origin[0] = dprod(femNodeXYZ,svec);
                    origin[1] = dprod(femNodeXYZ,tvec);

                    //global Cartesian coordinates of geometry nodes

                    vector<double> xyzVec(3,0.);

                    //transformed coordinates of geometry nodes

                    vector<double> xyzTilda(2,0.);

                    //pseudo-2D stencil generation in the computed basis 
                    for (int n=0;n<numForceInterpPts_;++n)
                      {
                        //get closest geometry node ID
                        int gn = geom2fem[fn].closestGeomNode[faceID][n];

                        xyzVec[0] = aux_nodes0_[gn*xyz + 0];
                        xyzVec[1] = aux_nodes0_[gn*xyz + 1];
                        xyzVec[2] = aux_nodes0_[gn*xyz + 2];

                        xyzTilda[0] = dprod(xyzVec,svec) - origin[0];
                        xyzTilda[1] = dprod(xyzVec,tvec) - origin[1];

                        loadStencils_[fn].loadCoefs[faceID][n] = 0.;
                        weights[n]       = 1.;
                        for (int k=0;k<numEqs;++k) tempm[n][k] = 0.;

                        tempm[n][0] = 1.;
                        tempm[n][1] = xyzTilda[0];
                        tempm[n][2] = xyzTilda[1];

                        if (forceInterpOrder_ == "second")
                          {
                            tempm[n][3] = xyzTilda[0]*xyzTilda[0];
                            tempm[n][4] = xyzTilda[1]*xyzTilda[1];
                            tempm[n][5] = xyzTilda[0]*xyzTilda[1];
                          }                        

                      }
                  }
                //generate coefficients
                if (auxForceInterpMethod_ == "LSQR") getLSQRcoefs(tempm,0,loadStencils_[fn].loadCoefs[faceID]);
                if (auxForceInterpMethod_ == "WLSQRX" || auxForceInterpMethod_ == "WLSQRZ") getWLSQRcoefs(tempm,weights,0,loadStencils_[fn].loadCoefs[faceID]);
                //sum test
                double sum = -1.;
                for (int n=0;n<numForceInterpPts_;++n)
                  {
                    sum += loadStencils_[fn].loadCoefs[faceID][n];
                    //nan check
                    if (loadStencils_[fn].loadCoefs[faceID][n] != loadStencils_[fn].loadCoefs[faceID][n]) Screen::MasterError("NaN in force transfer coefs!");
                  }
                if (sum > 1.E-3)
                  {
                    cout << "Computed interpolation coefficients do not pass sum test for FEM node " << participantNodeID << ", (x,y,z) = (" << nodes0_[participantNodeID*xyz+0] << "," << nodes0_[participantNodeID*xyz+1] << "," << nodes0_[participantNodeID*xyz+2] << ")" <<endl;
                    cout << "On geometry face " << faceID << endl;
                    cout << "Sum = " << sum << endl;
                    for (int n=0;n<numForceInterpPts_;++n)
                      {
                        int gn = geom2fem[fn].closestGeomNode[faceID][n];
                        cout << "geometry node " << gn << setw(15) << ", coef = " << loadStencils_[fn].loadCoefs[faceID][n] << setw(15) << ", (x,y,z) = " << aux_nodes0_[gn*xyz+0] << "," << aux_nodes0_[gn*xyz+1] << "," <<aux_nodes0_[gn*xyz+2] << ")" << endl;
                      }
                    Screen::MasterError("stopping the code");
                  }
              }
            faceID += 1;            
          }
      }
    else if (auxForceInterpMethod_ == "WLSQR")
      {
        int nEqs = 10;
        int order = 2;
        double weight = 1.;

        vector<vector<double> > tempm (numForceInterpPts_,vector<double>(nEqs,0.));
        vector<vector<double> > tempmT(nEqs,vector<double>(numForceInterpPts_,0.));

        for (int fb=-1;fb<=1;fb+=2)
          {
            for (int fn=0;fn<myNumberOfAuxParticipatingNodes_;++fn)
              {
                int participantNodeID = myAuxParticipatingNodes_[fn];
                for (int n=0;n<numForceInterpPts_;++n)
                  {
                    int gn = geom2fem[fn].closestGeomNode[faceID][n];

                    vector<double> xyzvec(xyz,0.);
                    for (int idir=0;idir<xyz;++idir)
                      {
                        int vecIndex = gn*xyz + idir;
                        xyzvec[idir]    = nodes0_[participantNodeID*xyz+idir] - aux_nodes0_[vecIndex];
                      }
                    //2nd-order
                    int p = 0;
                    for (int p0=0;p0<=order;++p0)
                      {
                        for (int p1=0;p1<=order;++p1)
                          {
                            for (int p2=0;p2<=order;++p2)
                              {
                                if (p0+p1+p2 <= order)
                                  {
                                    const double coef = pow(xyzvec[0],p0)*pow(xyzvec[1],p1)*pow(xyzvec[2],p2);
                                    tempm [n][p] = coef;
                                    tempmT[p][n] = weight*coef;
                                    p++;
                                  }
                              }
                          }
                      }
                  }
                vector<vector<double> > matTmat(nEqs,vector<double>(nEqs,0.));
                for (int k1=0;k1<nEqs;++k1)
                  {
                    for (int k2=0;k2<nEqs;++k2)
                      { 
                        for (int k3=0;k3<numForceInterpPts_;++k3)
                          { 
                            matTmat[k1][k2] += tempmT[k1][k3]*tempm[k3][k2];                        
                          }
                      }
                  }
                mat_inverse(matTmat);
                vector<vector<double> > mTmmT(nEqs,vector<double>(numForceInterpPts_,0.));
                for (int k1=0;k1<nEqs;++k1)
                  {
                    for (int k2=0;k2<numForceInterpPts_;++k2)
                      { 
                        for (int k3=0;k3<nEqs;++k3)
                          { 
                            mTmmT[k1][k2] += matTmat[k1][k3]*tempmT[k3][k2];          
                          }
                      }
                  } 

                for (int n=0;n<numForceInterpPts_;++n)
                  {
                    loadStencils_[fn].loadCoefs[faceID][n] = mTmmT[0][n];
                  }
                //sum test
                double sum = -1.;
                for (int n=0;n<numForceInterpPts_;++n)
                  {
                    sum += loadStencils_[fn].loadCoefs[faceID][n];
                    if (loadStencils_[fn].loadCoefs[faceID][n] != loadStencils_[fn].loadCoefs[faceID][n]) Screen::MasterError("NaN in force transfer coefs!");
                  }

                if (sum > 1.E-3)
                  {
                    cout << "Computed interpolation coefficients do not pass sum test for FEM node " << participantNodeID << ", (x,y,z) = " << nodes_[participantNodeID*xyz+0] << "," << nodes_[participantNodeID*xyz+1] << "," << nodes_[participantNodeID*xyz+2] << ")" <<endl;
                    cout << "On geometry face " << faceID << endl;
                    cout << "Sum = " << sum << endl;
                    for (int n=0;n<numForceInterpPts_;++n)
                      {
                        int gn = geom2fem[fn].closestGeomNode[faceID][n];
                        cout << "geometry node " << gn << setw(15) << ", coef = " << loadStencils_[fn].loadCoefs[faceID][n] << setw(15) << ", (x,y,z) = " << aux_nodes0_[gn*xyz+0] << "," << aux_nodes0_[gn*xyz+1] << "," <<aux_nodes0_[gn*xyz+2] << ")" << endl;
                      }
                    Screen::MasterError("stopping the code");                                      
                  }
              }
          }
      }

  }
 
  void fem::getDistanceOrder(vector<int>& ordered,
                             vector<double>& dist)
  {
    const int num = dist.size();
    vector<pair<double,int> > dist2_index(num);

    for (int i=0;i<ordered.size();++i)
      {
        dist2_index[i].first  = dist[i];
        dist2_index[i].second = i; 
      }

    stable_sort(dist2_index.begin(),dist2_index.end());

    for (int i=0;i<num;++i)
      {
        ordered[i] = dist2_index[i].second;
      }

  }
 
  void fem::auxLoadTransfer()
  {

    Screen::MasterInfo("Performing auxilliary load transfer");
    //nonlinear 2D with tri elements
    if (threeD_as_twoD_)
      {
        threeAstwoLoadTransfer();
      }
    //3D linear/nonlinear with tri elements
    else if (!threeD_as_twoD_ && DIM == 3)
      {
        threeDLoadTransfer();
      }
    //2D linear beam elements
    else if (!threeD_as_twoD_ && DIM == 2)
      {
        twoDLoadTransfer();
      }

  }
    
  void fem::twoDLoadTransfer()
  {
    for (int fn=0;fn<nof_nodes_glob_;++fn)
      {
        int ndof_loc = ndofPerNode_[fn];
        for (int n=0;n<numForceInterpPts_;++n)
          {
            int gn            = geom2fem[fn].closestGeomNode[0][n];
            int vecIndex      = gn*5 + dirVecIndex_;
            forceVec_[fn*ndof_loc+0] += auxGeomForces_[vecIndex]*defaultWidth_;
          }
      }

  }
    
  void fem::threeAstwoLoadTransfer()
  {

    if (auxForceInterpMethod_ == "double")
      {
	for (int fn=0;fn<nof_effective_nodes_;++fn)
	  {
	    int ndofOffset_loc = ndofLocalOffset_[fn];
	    //plus-minus geometry node
	    for (int faceID=0;faceID<2;++faceID)
	      {
		//in xy-plane only
		for (int dir=0;dir<2;++dir)
		  {
		    //closest geometry node on plus-minus face
		    int gn = geom2fem[fn].closestGeomNode[faceID][0];
		    int vecIndex = gn*5 + dir;
		    forceVec_[ndofOffset_loc + dir] += 0.5*auxGeomForces_[vecIndex]*defaultWidth_;
		  }
	      }
	  }
      }
    else if (auxForceInterpMethod_ == "single")
      {
	for (int fn=0;fn<nof_effective_nodes_;++fn)
	  {
	    int ndofOffset_loc = ndofLocalOffset_[fn];
	    //plus-minus geometry node
	    //in xy-plane only
	    for (int dir=0;dir<2;++dir)
	      {
		//closest geometry node on plus-minus face
		int gn = geom2fem[fn].closestGeomNode[0][0];
		int vecIndex = gn*5 + dir;
		forceVec_[ndofOffset_loc + dir] += 0.5*auxGeomForces_[vecIndex]*defaultWidth_;
	      }
	  }
      }
    else
      {
	Screen::MasterError("Cannot use anything except single or double for 3das2d");
      }
    //fill into-domain side
    for (int fn=0;fn<nof_effective_nodes_;++fn)
      {
        int node = fn + nof_effective_nodes_;
	int ndofOffset_loc  = ndofLocalOffset_[fn];
        int ndofOffset_loc2 = ndofLocalOffset_[node];
        for (int dir=0;dir<2;++dir)
          {
            forceVec_[ndofOffset_loc2 + dir] = forceVec_[ndofOffset_loc + dir];
          }
      }

  }
  
  void fem::threeDLoadTransfer()
  {
    //apply interpolation
    int faceID = 0;
    if (myAuxNodeCount_ > 0)
      {
	if (auxForceInterpMethod_ == "single")
	  {
	    for (int fn=0;fn<myNumberOfAuxParticipatingNodes_;++fn)
	      {
		int participantNodeID = myAuxParticipatingNodes_[fn];
		int ndof_loc          = ndofPerNode_[participantNodeID];
		int ndofOffset_loc    = ndofLocalOffset_[participantNodeID];

		for (int n=0;n<numForceInterpPts_;++n)
		  {                
		    int gn = geom2fem[fn].closestGeomNode[faceID][n];
		    if ( abs(is_aux_node_[gn]) == 99)
		      {
			for (int d=0;d<xyz;++d) 
			  {
			    int vecIndex = gn*5 + d;
			    forceVec_[d + ndofOffset_loc] += loadStencils_[fn].loadCoefs[faceID][n]*auxGeomForces_[vecIndex];
			  }
		      }
		  }
	      }
	  }
	else
	  {
	    for (int fb=-1;fb<=1;fb+=2)
	      {        
		for (int fn=0;fn<myNumberOfAuxParticipatingNodes_;++fn)
		  {
		    int participantNodeID = myAuxParticipatingNodes_[fn];
		    int ndof_loc          = ndofPerNode_[participantNodeID];
		    int ndofOffset_loc    = ndofLocalOffset_[participantNodeID];

		    for (int n=0;n<numForceInterpPts_;++n)
		      {                
			int gn = geom2fem[fn].closestGeomNode[faceID][n];
			if ( abs(is_aux_node_[gn]) == 99)
			  {
			    for (int d=0;d<xyz;++d) 
			      {
				int vecIndex = gn*5 + d;
				forceVec_[d + ndofOffset_loc] += loadStencils_[fn].loadCoefs[faceID][n]*auxGeomForces_[vecIndex];
			      }
			  }
		      }
		  }
		faceID += 1;
	      }
	  }
      }

#if (DUMP_PRESSURE_DISTRIBUTION == 1)

    vector<double> rootForces;
    vector<int> orderedList;

    if (it_ == 1)
      {

        if (nt_%ntSkip_ == 0)
          {

            Screen::MasterInfo("Dumping pressure distribution that was interpolated from geometry");
            vector<double> faceForces(3*myElementCount_,0.);
            for (int e=0;e<myElementCount_;++e)
              {
                int nnodes_loc    = elementData_[e].nnodes;
                int nnodes_offset = nnodesLocalOffset_[e];
                for (int n=0;n<nnodes_loc;++n)
                  {
                    int node           = faces_[n + nnodes_offset];
                    int ndofOffset_loc = ndofLocalOffset_[node];
                    faceForces[e*xyz + 0] += 0.33333333333333*forceVec_[ndofOffset_loc + 0];
                    faceForces[e*xyz + 1] += 0.33333333333333*forceVec_[ndofOffset_loc + 1];
                    faceForces[e*xyz + 2] += 0.33333333333333*forceVec_[ndofOffset_loc + 2];
                  }
              }

            int myNumber = 3*myElementCount_;
            int myOffset = mypeno_*myNumber;
            vector<int> partitionCellCounts(numprocs_,0);
            vector<int> partitionOffset    (numprocs_,0);

            MPI_Allgather(&myNumber,1,MPI::INT,&partitionCellCounts[0],1,MPI::INT,MPI_COMM_WORLD);
            partitionOffset[0] = 0;
            for (int p=1;p<numprocs_;++p) partitionOffset[p] = partitionOffset[p-1] + partitionCellCounts[p-1];
            //gather all forces on root for output
            if (mypeno_ == 0) rootForces.resize(3*nof_faces_glob_,0.);
            ierr = MPI_Gatherv(&faceForces[0],partitionCellCounts[mypeno_],MPI::DOUBLE,&rootForces[0],&partitionCellCounts[0],&partitionOffset[0],MPI::DOUBLE,0,MPI_COMM_WORLD);

            vector<int> myGlobalIDs(myElementCount_,-1);               
            //re-order proc-ordered gather to the global ID order
            for (int e=0;e<myElementCount_;++e)
              {
                myGlobalIDs[e] = elemMaps_.elemLocal2Global[e];
              }

            //gather global ID ordering on root
            if (mypeno_ == 0) orderedList.resize(nof_faces_glob_,-1);
            ierr = MPI_Gatherv(&myGlobalIDs[0],partitionCellCounts_[mypeno_],MPI::INT,&orderedList[0],&partitionCellCounts_[0],&partitionOffset_[0],MPI::INT,0,MPI_COMM_WORLD);

          }

        if (mypeno_ == 0 && nt_%ntSkip_ == 0)
          {

            Screen::MasterInfo("Writing interpolated pressures from root");
            ofstream surf_out;
            string surface_filename;
            ostringstream OSS;
            surfaceTriangulation surfTri;
            OSS << "debugFEM/interpolatedPressureDistribution_" << setw(8) << setfill('0') << nt_ << ".vtk";
            surface_filename = OSS.str();
            surf_out.open(surface_filename);
            vector<double> rootDisps;
            rootDisps.resize(3*nof_faces_glob_,0.);
            surfTri.smartWriteTriangulation_withForces(surface_filename,rootTotalNodes_,rootTotalFaces_,xyz,nof_nodes_glob_,nof_faces_glob_,numberOfElementTypes_,rootNnodesPerEl_,rootForces,orderedList,rootDisps);
            surf_out.close();

          }
      }
#endif

    if (myAuxNodeCount_ > 0)
      {
        //express interpolated nodal pressures as nodal force
        //1) move nodal pressure to centroid
        if (DIM == 3 && !threeD_as_twoD_)
          {
            for (int f=0;f<myNumberOfAuxParticipatingElements_;++f)
              {
                int participantElementID = myAuxParticipatingElements_[f];
                int nnodes_loc           = elementData_[participantElementID].nnodes;
                int nnodesOffset_loc     = nnodesLocalOffset_[participantElementID];
                vector<double> centroidPressure(xyz,0.);

                for (int idir=0;idir<xyz;++idir)
                  {
                    for (int n=0;n<nnodes_loc;++n)
                      {
                        int index               = n + nnodesOffset_loc;
                        int node                = faces_[index];
                        int ndofOffset_loc      = ndofLocalOffset_[node];
                        centroidPressure[idir] += forceVec_[idir + ndofOffset_loc]/float(nnodes_loc);
                      }
                    //2) calculate force at element centroid
                    centroidForce_[f][idir] = Area_[participantElementID] * centroidPressure[idir];
                  }
              }
            //3) prepare forceVec for nodal forces
            fill(forceVec_.begin(),forceVec_.end(),0.);
            //4) redistribute face force to nodes
            for (int f=0;f<myNumberOfAuxParticipatingElements_;++f)
              {
                int participantElementID = myAuxParticipatingElements_[f];
                int nnodes_loc           = elementData_[participantElementID].nnodes;
                int nnodesOffset_loc     = nnodesLocalOffset_[f];

                for (int n=0;n<nnodes_loc;++n)
                  {
                    int index          = n + nnodesOffset_loc;
                    int node           = faces_[index];
                    int ndofOffset_loc = ndofLocalOffset_[node];

                    for (int idir=0;idir<xyz;++idir)
                      {                    
                        forceVec_[idir + ndofOffset_loc] += centroidForce_[f][idir]/float(nnodes_loc);
                      }
                  }
              }
          }
      }

  }
  
  void fem::initializeAuxillaryBodyVelocity()
  {

    //intialize geometry triangulation for velocity computation at restart
    Screen::MasterInfo("Initializing geometry body velocity");

    vector<double> myFEMNodeVel(myNodeCount_*xyz,0.);
    vector<double> myAuxNodeVel(myAuxNodeCount_*xyz,0.);
    vector<double> iniAuxNodeVel_loc(numGeomNodes_*xyz,0.);
    iniAuxNodeVel_.resize           (numGeomNodes_*xyz,0.);

    //compute FEM nodal velocities
    double invdt = 1./dt_;
    for (int n=0;n<myNodeCount_;++n)
      {
        int ndof_offset = ndofLocalOffset_[n];
        for (int idir=0;idir<xyz;++idir)
          {
            myFEMNodeVel[n*xyz + idir] = invdt * deltaUi_[ndof_offset + idir];
          }
      }

    if (myNumberOfAuxParticipatingNodes_ > 0)
      {	
	//interpolate to geometry nodes
	for (int mn=0;mn<movingNodes_.size();++mn)
	  {
#if (DEDICATED_FEM_PROC == 1)
	    int gn = mn;
#else
	    int gn = mn;
#endif

	    if ( abs(is_aux_node_[movingNodes_[gn]]) == 99 )
	      {
		for (int fn=0;fn<leastSQR_npnts_;++fn)
		  {
		    for (int idir=0;idir<xyz;++idir)
		      {
			myAuxNodeVel[gn*xyz+idir] += stencils_[gn].coefs[fn] * myFEMNodeVel[stencils_[gn].nodes[fn]];
		      }
		  }
	      }
	  }
      }
    //consolidate all procs aux node vel
#if (DEDICATED_FEM_PROC == 0)
    int send = myAuxNodeCount_*xyz;
    ierr = MPI_Allgatherv(&myAuxNodeVel[0],send,MPI::DOUBLE,&iniAuxNodeVel_loc[0],&myAuxNodeNumVecXYZ_[0],&nodeXYZOffset_[0],MPI::DOUBLE,MPI_COMM_WORLD);
#else
    iniAuxNodeVel_loc.assign(myAuxNodeVel.begin(),myAuxNodeVel.end());
#endif

    //memory
    myAuxNodeVel = vector<double>();

    ierr = MPI_Barrier(MPI_COMM_WORLD);
    Screen::MasterInfo("before");

    //use stored mapping to unscramble geometry node velocity
    for (int gn=0;gn<numGeomNodes_;++gn)
      {
#if (DEDICATED_FEM_PROC == 0)
        int node  = geometryMap_[gn];
#else
        int node = gn;
#endif

	if ( abs(is_aux_node_[node]) == 99 )
	  {

	    for (int idir=0;idir<xyz;++idir)
	      {
		iniAuxNodeVel_[node*xyz+idir] = iniAuxNodeVel_loc[gn*xyz+idir];
	      }
	  }
      }
    
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    Screen::MasterInfo("final");

  }
  
  void fem::auxDispTransfer()
  {

#if (OVERWRITE_FEM_NODES==1)

    string filename = "320k_aspireSR01_msl_80gore_21m_dgb_parachute_FEM_Umbrella_IC_0p25_shrink_ASPIRE_SR01_Riser.vtk";
    Screen::MasterWarning("Reading file: "+filename+" to overwrite FEM node locations");
    ifstream ovwrite;
    ovwrite.open(filename);
    if (ovwrite.fail()) Screen::MasterError("FEM node overwrite file: "+filename+" not found!");

    vector<double> globalNodes;
    surfaceTriangulation surfTri;
    surfTri.readTriangulation_AUX_FSI(globalNodes,filename);

    //overwrite nodal locations in normal map
    for (int n=0;n<myNodeCount_;++n)
      {
	int globID = elemMaps_.nodeLocal2Global[n];
	for (int idir=0;idir<xyz;++idir)
	  {
	    nodes_[n*xyz+idir] = globalNodes[globID*xyz+idir];
	  }
      }

    //overwrite nodal locations in extended map
    for (int n=0;n<myNodeCount_ext_;++n)
      {
	int globID = elemMaps_.nodeLocal2Global_ext[n];
	for (int idir=0;idir<xyz;++idir)
	  {
	    nodes_ext_[n*xyz+idir] = globalNodes[globID*xyz+idir];
	  }
      }

    globalNodes = vector<double>();

    //update geometric normal vectors for displacement transfer
    updateGeometricNormals();

#endif

    Screen::MasterInfo("Transferring displacements from FEM to geometry");
    //replace movingNode positions with FEM nodes pushed out by plus-minusnormal vector
    for (int mn=0;mn<movingNodes_.size();++mn)
      {
#if (DEDICATED_FEM_PROC == 1)
	int gn = mn;
#else
	int gn = mn;
#endif
	for (int idir=0;idir<xyz;++idir)
	  {                
	    myAuxNodes_[gn*xyz+idir] = 0.;
	  }
      }

    if (threeD_as_twoD_)
      {
	for (int mn=0;mn<movingNodes_.size();++mn)
	  {
	    int gn = mn;
	    for (int fn=0;fn<leastSQR_npnts_;++fn)
	      {
		int vectorNode = myLocalNodesNormalVectorStatus_[stencils_[gn].nodes[fn]];
		for (int idir=0;idir<dimLoc_;++idir)
		  {
		    int vecIndex           = gn*xyz + idir;
		    myAuxNodes_[vecIndex] += stencils_[gn].coefs[fn] * (nodes_[stencils_[gn].nodes[fn]*xyz+idir] + 0.5*geomThickFactor_*defaultThickness_*stencils_[gn].sign*Vnt_[vectorNode*xyz+idir]);
		  }
	      }
	  }
      }
    else
      {
	//update all three spatial coordinates
	//(one spatial coordinate is zero for 2D sims anyway)
	if (myNumberOfAuxParticipatingNodes_ > 0)
	  {
	    for (int mn=0;mn<movingNodes_.size();++mn)
	      {
#if (DEDICATED_FEM_PROC == 1)
		int gn = mn;
#else
		int gn = mn;
#endif

		if ( abs(is_aux_node_[movingNodes_[gn]]) == 99 )
		  {
		    for (int fn=0;fn<leastSQR_npnts_;++fn)
		      {
			int vectorNode = myLocalNodesNormalVectorStatus_[stencils_[gn].nodes[fn]];
			for (int idir=0;idir<xyz;++idir)
			  {
			    int vecIndex           = gn*xyz + idir;
			    myAuxNodes_[vecIndex] += stencils_[gn].coefs[fn] * (nodes_[stencils_[gn].nodes[fn]*xyz+idir] + 0.5*geomThickFactor_*defaultThickness_*stencils_[gn].sign*Vnt_geom_[vectorNode*xyz+idir]);
			  }
		      }
		  }	       
	      }
	  }
      }
    //consolidate all procs aux nodes into aux_nodes_copy_
#if (DEDICATED_FEM_PROC == 0)
    int send = myAuxNodeCount_*xyz;
    ierr = MPI_Allgatherv(&myAuxNodes_[0],send,MPI::DOUBLE,&aux_nodes_copy_[0],&myAuxNodeNumVecXYZ_[0],&nodeXYZOffset_[0],MPI::DOUBLE,MPI_COMM_WORLD);
    //use stored mapping to unscramble aux_nodes_copy_ into aux_nodes_
    for (int gn=0;gn<numGeomNodes_;++gn)
      {
	int node  = geometryMap_[gn];
	//only update movingNodes that I am responsible for
	if ( abs(is_aux_node_[node]) == 99 )
	  {
	    for (int idir=0;idir<xyz;++idir)
	      {
		aux_nodes_[node*xyz+idir] = aux_nodes_copy_[gn*xyz+idir];
	      }
	  }
      }     
#else
    for (int gn=0;gn<movingNodes_.size();++gn)
      {
	for (int idir=0;idir<xyz;++idir)
	  {
	    int node = movingNodes_[gn];
	    aux_nodes_[node*xyz+idir] = myAuxNodes_[gn*xyz+idir];
	  }
      }
#endif
    //dump geometry
    if (useAuxGeometry_ && nt_%ntSkip_ == 0 && mypeno_ == 0)
      //if (mypeno_ == 0)
      {
	string name;
	ostringstream OSS;
	OSS << "post/aux/displacedAux_" << setw(8) << setfill('0') << nt_ << ".vtk";
	name = OSS.str();

	surfaceTriangulation surfTri;
	surfTri.writeTriangulation_aux(name,aux_nodes_,aux_faces_,dimLoc_);
      }

  }
  
 void fem::getLSQRcoefs(vector<vector<double> >& mat,
                        int ID,
                        vector<double>& coefs)
 {
   const int rows = mat.size();
   const int cols = mat[0].size();
   vector<double> dummyCol(cols,0.);
   vector<vector<double> > a(cols,dummyCol);
   vector<double> dummyRow(rows,0.);
   vector<vector<double> > b(cols,dummyRow);

   mTm_mult(mat,mat,a);
   mat_inverse(a);
   mmT_mult(a,mat,b);

   for (int n=0;n<rows;++n)
     {
       coefs[n] = b[ID][n];
     }

 }

 void fem::getWLSQRcoefs(vector<vector<double> >& mat,
                         vector<double>& weights,
                         int ID,
                         vector<double>& coefs)
 {
   const int rows = mat.size();
   const int cols = mat[0].size();
   vector<double> dummyCol(cols,0.);
   vector<vector<double> > a(cols,dummyCol);
   vector<double> dummyRow(rows,0.);
   vector<vector<double> > b(cols,dummyRow);
   vector<vector<double> > WmatT(cols,dummyRow);

   for (int m=0;m<cols;++m)
     {
       for (int n=0;n<rows;++n)
         {
           WmatT[m][n] = weights[n] * mat[n][m];
         }
     }

   mm_mult(WmatT,mat,a);
   mat_inverse(a);
   mm_mult(a,WmatT,b);

   for (int n=0;n<rows;++n)
     {
       coefs[n] = b[ID][n];
     }

 }

 void fem::updateDirectorVecs()
 {
   //update Vnt,V1t,V2t after each step
   //nonlinear - second order effects
   //grab alpha and beta displacements from deltaUi_
   int alphaIndex = 3;
   int betaIndex  = 4;

   if (dirVecIndex_ < 0)
     {
       for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
         {

           int node       = myNormalVectorNodes_[n];
           int ndofOffset = ndofLocalOffset_[node];
           double alpha   = deltaUi_[alphaIndex + ndofOffset];
           double beta    = deltaUi_[betaIndex  + ndofOffset];

           for (int idir=0;idir<xyz;++idir)
             { 

               deltaVnt_[n*xyz+idir]     = -alpha*V2t_[n*xyz+idir] + beta*V1t_[n*xyz+idir] - 0.5*(alpha*alpha+beta*beta)*Vnt_[n*xyz+idir];

               Vnt_[n*xyz+idir]         += deltaVnt_[n*xyz+idir];

               myDirVecs_[n*xyz + idir] += deltaVnt_[n*xyz+idir];

             }

           //update V1t and V2t using eq (11) from MITC3+ geometric nonlinear
           V1t_[n*xyz    ] = myDirVecs_[n*xyz+1]*Vnt_[n*xyz+2] - myDirVecs_[n*xyz+2]*Vnt_[n*xyz+1];
           V1t_[n*xyz + 1] = myDirVecs_[n*xyz+2]*Vnt_[n*xyz+0] - myDirVecs_[n*xyz+0]*Vnt_[n*xyz+2];
           V1t_[n*xyz + 2] = myDirVecs_[n*xyz+0]*Vnt_[n*xyz+1] - myDirVecs_[n*xyz+1]*Vnt_[n*xyz+0];

           V2t_[n*xyz    ] = Vnt_[n*xyz + 1]*V1t_[n*xyz + 2] - Vnt_[n*xyz + 2]*V1t_[n*xyz+1];
           V2t_[n*xyz + 1] = Vnt_[n*xyz + 2]*V1t_[n*xyz + 0] - Vnt_[n*xyz + 0]*V1t_[n*xyz+2];
           V2t_[n*xyz + 2] = Vnt_[n*xyz    ]*V1t_[n*xyz + 1] - Vnt_[n*xyz + 1]*V1t_[n*xyz  ];

         }

       normalizeVector3(Vnt_);
       normalizeVector3(myDirVecs_);
       normalizeVector3(V1t_);
       normalizeVector3(V2t_);

     }
   else
     {
       for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
         {
           int node       = myNormalVectorNodes_[n];
           int ndofOffset = ndofLocalOffset_[node];
           double alpha   = deltaUi_[alphaIndex + ndofOffset];
           double beta    = deltaUi_[betaIndex  + ndofOffset];

           for (int idir=0;idir<xyz;++idir) 
             { 

               deltaVnt_[n*xyz+idir] = -alpha*V2t_[n*xyz+idir] + beta*V1t_[n*xyz+idir] - 0.5*(alpha*alpha+beta*beta)*Vnt_[n*xyz+idir];

               Vnt_[n*xyz+idir]  += deltaVnt_[n*xyz+idir];

             }
           //update V1t and V2t using eq (11) from MITC3+ geometric nonlinear
           V1t_[n*xyz    ] = myDirVecs_[n*xyz+1]*Vnt_[n*xyz+2] - myDirVecs_[n*xyz+2]*Vnt_[n*xyz+1];
           V1t_[n*xyz + 1] = myDirVecs_[n*xyz+2]*Vnt_[n*xyz+0] - myDirVecs_[n*xyz+0]*Vnt_[n*xyz+2];
           V1t_[n*xyz + 2] = myDirVecs_[n*xyz+0]*Vnt_[n*xyz+1] - myDirVecs_[n*xyz+1]*Vnt_[n*xyz+0];

           V2t_[n*xyz    ] = Vnt_[n*xyz + 1]*V1t_[n*xyz + 2] - Vnt_[n*xyz + 2]*V1t_[n*xyz+1];
           V2t_[n*xyz + 1] = Vnt_[n*xyz + 2]*V1t_[n*xyz + 0] - Vnt_[n*xyz + 0]*V1t_[n*xyz+2];
           V2t_[n*xyz + 2] = Vnt_[n*xyz    ]*V1t_[n*xyz + 1] - Vnt_[n*xyz + 1]*V1t_[n*xyz  ]; 

         }

       normalizeVector3(Vnt_);
       normalizeVector3(V1t_);
       normalizeVector3(V2t_);
     }

   //removing normal vector calculation for beams for now
   //    //beams   
   //    else
   //      {
   //        vector<double> tangent(xyz,0.);
   //        vector<double> intoDomain(xyz,0.);
   //        intoDomain[2] = 1.;
   //        for (int n=0;n<myNodeCount_-1;++n)
   //          {
   //            for (int idir=0;idir<xyz;++idir) tangent[idir] = nodes_[(n+1)*xyz+idir] - nodes_[n*xyz+idir];
   //            Vnt_[n] = cross(intoDomain,tangent);
   //          }
   //        //treat end node
   //        for (int idir=0;idir<xyz;++idir) tangent[idir] = nodes_[(myNodeCount_-1)*xyz+idir] - nodes_[(myNodeCount_-2)*xyz+idir];
   //        Vnt_[nof_nodes_glob_-1] = cross(intoDomain,tangent);
   //        //normalize
   //        for (int n=0;n<myNodeCount_;++n) normalizeVector(Vnt_[n]);
   //      }
   //  }

 }
 
 void fem::unUpdateDirectorVecs()
 {
   //update Vnt,V1t,V2t after each step
   //nonlinear - second order effects
   //grab alpha and beta displacements from deltaUi_
   int alphaIndex = 3;
   int betaIndex  = 4;

   if (dirVecIndex_ < 0)
     {
       for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
         {

           int node       = myNormalVectorNodes_[n];
           int ndofOffset = ndofLocalOffset_[node];
           double alpha   = deltaUi_[alphaIndex + ndofOffset];
           double beta    = deltaUi_[betaIndex  + ndofOffset];

           for (int idir=0;idir<xyz;++idir)
             {

               Vnt_[n*xyz+idir] = Vnt_bak_[n*xyz+idir];

               myDirVecs_[n*xyz + idir] = myDirVecs_bak_[n*xyz + idir];
             }

           //update V1t and V2t using eq (11) from MITC3+ geometric nonlinear
           V1t_[n*xyz    ] = myDirVecs_[n*xyz+1]*Vnt_[n*xyz+2] - myDirVecs_[n*xyz+2]*Vnt_[n*xyz+1];
           V1t_[n*xyz + 1] = myDirVecs_[n*xyz+2]*Vnt_[n*xyz+0] - myDirVecs_[n*xyz+0]*Vnt_[n*xyz+2];
           V1t_[n*xyz + 2] = myDirVecs_[n*xyz+0]*Vnt_[n*xyz+1] - myDirVecs_[n*xyz+1]*Vnt_[n*xyz+0];

           V2t_[n*xyz    ] = Vnt_[n*xyz + 1]*V1t_[n*xyz + 2] - Vnt_[n*xyz + 2]*V1t_[n*xyz+1];
           V2t_[n*xyz + 1] = Vnt_[n*xyz + 2]*V1t_[n*xyz + 0] - Vnt_[n*xyz + 0]*V1t_[n*xyz+2];
           V2t_[n*xyz + 2] = Vnt_[n*xyz    ]*V1t_[n*xyz + 1] - Vnt_[n*xyz + 1]*V1t_[n*xyz  ];

         }

       normalizeVector3(Vnt_);
       normalizeVector3(myDirVecs_);
       normalizeVector3(V1t_);
       normalizeVector3(V2t_);

     }
   else
     {
       for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
         {
           for (int idir=0;idir<xyz;++idir) Vnt_[n*xyz+idir] = Vnt_bak_[n*xyz+idir];

           //update V1t and V2t using eq (11) from MITC3+ geometric nonlinear
           V1t_[n*xyz    ] = myDirVecs_[n*xyz+1]*Vnt_[n*xyz+2] - myDirVecs_[n*xyz+2]*Vnt_[n*xyz+1];
           V1t_[n*xyz + 1] = myDirVecs_[n*xyz+2]*Vnt_[n*xyz+0] - myDirVecs_[n*xyz+0]*Vnt_[n*xyz+2];
           V1t_[n*xyz + 2] = myDirVecs_[n*xyz+0]*Vnt_[n*xyz+1] - myDirVecs_[n*xyz+1]*Vnt_[n*xyz+0];

           V2t_[n*xyz    ] = Vnt_[n*xyz + 1]*V1t_[n*xyz + 2] - Vnt_[n*xyz + 2]*V1t_[n*xyz+1];
           V2t_[n*xyz + 1] = Vnt_[n*xyz + 2]*V1t_[n*xyz + 0] - Vnt_[n*xyz + 0]*V1t_[n*xyz+2];
           V2t_[n*xyz + 2] = Vnt_[n*xyz    ]*V1t_[n*xyz + 1] - Vnt_[n*xyz + 1]*V1t_[n*xyz  ]; 

         }

       normalizeVector3(Vnt_);
       normalizeVector3(V1t_);
       normalizeVector3(V2t_);
     }

   //removing normal vector calculation for beams for now
   //    //beams   
   //    else
   //      {
   //        vector<double> tangent(xyz,0.);
   //        vector<double> intoDomain(xyz,0.);
   //        intoDomain[2] = 1.;
   //        for (int n=0;n<myNodeCount_-1;++n)
   //          {
   //            for (int idir=0;idir<xyz;++idir) tangent[idir] = nodes_[(n+1)*xyz+idir] - nodes_[n*xyz+idir];
   //            Vnt_[n] = cross(intoDomain,tangent);
   //          }
   //        //treat end node
   //        for (int idir=0;idir<xyz;++idir) tangent[idir] = nodes_[(myNodeCount_-1)*xyz+idir] - nodes_[(myNodeCount_-2)*xyz+idir];
   //        Vnt_[nof_nodes_glob_-1] = cross(intoDomain,tangent);
   //        //normalize
   //        for (int n=0;n<myNodeCount_;++n) normalizeVector(Vnt_[n]);
   //      }
   //  }

 }

 void fem::updateGeometricNormals()
 {   
   //use ghost/halo elements to compute geometric node normals for displacement transfer
   Screen::MasterInfo("Updating geometric nodal normal vectors for displacement transfer");
   vector<double> vec1(xyz,0.);
   vector<double> vec2(xyz,0.);
   fill(Vnt_geom_.begin(),Vnt_geom_.end(),0.);

   if (revertNormalsByCompID_)
     {
       for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
	 {
	   int localNodeID = myNormalVectorNodes_[n];
	   int numSharedElements = 0;

	   //loop over extended face list, include ghosts
	   for (int f=0;f<myElementCount_ext_;++f)
	     {
	       //hard coded for triangles
	       if (elementData_ext_[f].nnodes == 3)
		 {
		   int compID = elementData_ext_[f].compID;
		   int nnodesOffset_ext = nnodesLocalOffset_ext_[f];

		   //see which faces in the extended map owns this normal vector node
		   if ( faces_ext_[0 + nnodesOffset_ext] == localNodeID ||
			faces_ext_[1 + nnodesOffset_ext] == localNodeID || 
			faces_ext_[2 + nnodesOffset_ext] == localNodeID )
		     {
		       numSharedElements += 1;

		       if ( find(revertNormalsCompIDList_.begin(),revertNormalsCompIDList_.end(),abs(compID)) == revertNormalsCompIDList_.end())
			 {
			   for (int idir=0;idir<xyz;++idir)
			     {
			       vec1[idir] = nodes_ext_[faces_ext_[1 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
			       vec2[idir] = nodes_ext_[faces_ext_[2 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
			     }
			 }
		       else
			 {
			   for (int idir=0;idir<xyz;++idir)
			     {
			       vec2[idir] = nodes_ext_[faces_ext_[1 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
			       vec1[idir] = nodes_ext_[faces_ext_[2 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
			     }
			 }

		       Vnt_geom_[n*xyz + 0] += vec1[1]*vec2[2] - vec1[2]*vec2[1];
		       Vnt_geom_[n*xyz + 1] += vec1[2]*vec2[0] - vec1[0]*vec2[2];
		       Vnt_geom_[n*xyz + 2] += vec1[0]*vec2[1] - vec1[1]*vec2[0];

		     }
		 }
	     }
	   for (int idir=0;idir<xyz;++idir) Vnt_geom_[n*xyz + idir] = Vnt_geom_[n*xyz + idir]/float(numSharedElements);
	   normalizeVector3(Vnt_geom_);
	 }
     }
   else
     {
       for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
	 {
	   int localNodeID = myNormalVectorNodes_[n];
	   int numSharedElements = 0;

	   //loop over extended face list, include ghosts
	   for (int f=0;f<myElementCount_ext_;++f)
	     {
	       //hard coded for triangles
	       if (elementData_ext_[f].nnodes == 3)
		 {
		   int nnodesOffset_ext = nnodesLocalOffset_ext_[f];

		   //see which faces in the extended map owns this normal vector node
		   if ( faces_ext_[0 + nnodesOffset_ext] == localNodeID ||
			faces_ext_[1 + nnodesOffset_ext] == localNodeID || 
			faces_ext_[2 + nnodesOffset_ext] == localNodeID )
		     {
		       numSharedElements += 1;

		       for (int idir=0;idir<xyz;++idir)
			 {
			   vec1[idir] = nodes_ext_[faces_ext_[1 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
			   vec2[idir] = nodes_ext_[faces_ext_[2 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
			 }

		       Vnt_geom_[n*xyz + 0] += vec1[1]*vec2[2] - vec1[2]*vec2[1];
		       Vnt_geom_[n*xyz + 1] += vec1[2]*vec2[0] - vec1[0]*vec2[2];
		       Vnt_geom_[n*xyz + 2] += vec1[0]*vec2[1] - vec1[1]*vec2[0];

		     }
		 }
	     }
	   for (int idir=0;idir<xyz;++idir) Vnt_geom_[n*xyz + idir] = Vnt_geom_[n*xyz + idir]/float(numSharedElements);
	   normalizeVector3(Vnt_geom_);
	 }
     }

   //treat special corner nodes
   //normal vector is the average of the 'plus' and 'minus' averages
   //this assumes aux participating nodes are also normal vector nodes

   if (auxCornerTreatment_)
     {

       Screen::MasterInfo("Treating special corner nodes in the geometric normal update");
       vector<int> numSharedElementsVec(2,0);
       vector<vector<double> > nVec(2,vector<double>(xyz,0.));

       if (revertNormalsByCompID_)
	 {	   
	   for (int n=0;n<myNumberOfAuxCornerTreatmentNodes_;++n)
	     {
	       int localNodeID  = myAuxCornerTreatmentNodes_[n];
	       int globalNodeID = elemMaps_.nodeLocal2Global[localNodeID];
	       int vectorNode   = myLocalNodesNormalVectorStatus_[localNodeID];
	       int sign_i = -1;

	       for (int pmFace=0;pmFace<2;++pmFace)
		 {

		   nVec[pmFace][0] = 0.;
		   nVec[pmFace][1] = 0.;
		   nVec[pmFace][2] = 0.;
		   numSharedElementsVec[pmFace] = 0;

		   //loop over extended face list, include ghosts
		   for (int f=0;f<myElementCount_ext_;++f)
		     {
		       //hard coded for triangles, and test plus-minus compID
		       if ( (elementData_ext_[f].nnodes == 3) && ((elementData_ext_[f].compID*sign_i)>1) )
			 {
			   int compID = elementData_ext_[f].compID;
			   int nnodesOffset_ext = nnodesLocalOffset_ext_[f];

			   //see which faces in the extended map own this normal vector/auxillary participating node
			   if ( faces_ext_[0 + nnodesOffset_ext] == localNodeID ||
				faces_ext_[1 + nnodesOffset_ext] == localNodeID || 
				faces_ext_[2 + nnodesOffset_ext] == localNodeID )
			     {
			       numSharedElementsVec[pmFace] += 1;

			       if ( find(revertNormalsCompIDList_.begin(),revertNormalsCompIDList_.end(),abs(compID)) == revertNormalsCompIDList_.end())
				 {
				   for (int idir=0;idir<xyz;++idir)
				     {
				       vec1[idir] = nodes_ext_[faces_ext_[1 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
				       vec2[idir] = nodes_ext_[faces_ext_[2 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
				     }
				 }
			       else
				 {
				   for (int idir=0;idir<xyz;++idir)
				     {
				       vec2[idir] = nodes_ext_[faces_ext_[1 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
				       vec1[idir] = nodes_ext_[faces_ext_[2 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
				     }
				 }

			       nVec[pmFace][0] += vec1[1]*vec2[2] - vec1[2]*vec2[1];
			       nVec[pmFace][1] += vec1[2]*vec2[0] - vec1[0]*vec2[2];
			       nVec[pmFace][2] += vec1[0]*vec2[1] - vec1[1]*vec2[0];

			     }
			 }
		     }

		   sign_i += 2;
		   for (int idir=0;idir<xyz;++idir) nVec[pmFace][idir] = nVec[pmFace][idir]/float(numSharedElementsVec[pmFace]);
		   normalizeVector(nVec[pmFace]);

		 }

	       for (int idir=0;idir<xyz;++idir) Vnt_geom_[vectorNode*xyz + idir] = 0.5 * (nVec[0][idir] + nVec[1][idir]);
	       normalizeVector3(Vnt_geom_);

	     }
	 }
       else
	 {
	   for (int n=0;n<myNumberOfAuxCornerTreatmentNodes_;++n)
	     {
	       int localNodeID = myAuxCornerTreatmentNodes_[n];
	       int sign_i = -1;

	       for (int pmFace=0;pmFace<2;++pmFace)
		 {

		   nVec[pmFace][0] = 0.;
		   nVec[pmFace][1] = 0.;
		   nVec[pmFace][2] = 0.;
		   numSharedElementsVec[pmFace] = 0;

		   //loop over extended face list, include ghosts
		   for (int f=0;f<myElementCount_ext_;++f)
		     {
		       //hard coded for triangles, and test plus-minus compID
		       if ( (elementData_ext_[f].nnodes == 3) && ((elementData_ext_[f].compID*sign_i)>1) )
			 {
			   int nnodesOffset_ext = nnodesLocalOffset_ext_[f];

			   //see which faces in the extended map own this normal vector/auxillary participating node
			   if ( faces_ext_[0 + nnodesOffset_ext] == localNodeID ||
				faces_ext_[1 + nnodesOffset_ext] == localNodeID || 
				faces_ext_[2 + nnodesOffset_ext] == localNodeID )
			     {
			       numSharedElementsVec[pmFace] += 1;

			       for (int idir=0;idir<xyz;++idir)
				 {
				   vec1[idir] = nodes_ext_[faces_ext_[1 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
				   vec2[idir] = nodes_ext_[faces_ext_[2 + nnodesOffset_ext]*xyz + idir] - nodes_ext_[faces_ext_[0 + nnodesOffset_ext]*xyz + idir];
				 }				 

			       nVec[pmFace][0] += vec1[1]*vec2[2] - vec1[2]*vec2[1];
			       nVec[pmFace][1] += vec1[2]*vec2[0] - vec1[0]*vec2[2];
			       nVec[pmFace][2] += vec1[0]*vec2[1] - vec1[1]*vec2[0];

			     }
			 }
		     }

		   sign_i += 2;
		   for (int idir=0;idir<xyz;++idir) nVec[pmFace][idir] = nVec[pmFace][idir]/float(numSharedElementsVec[pmFace]);
		   normalizeVector(nVec[pmFace]);		   

		 }

	       for (int idir=0;idir<xyz;++idir) Vnt_geom_[localNodeID*xyz + idir] = 0.5 * (nVec[0][idir] + nVec[1][idir]);

	     }
	   normalizeVector3(Vnt_geom_);
	 }
     }

 }

 void fem::identifyAuxCornerTreatmentNodes()
 {
   //seek for nodes on the intersection of 'plus-minus' compID lines
   //this routine was written assuming that 'aux participating' nodes and elements are shells
   Screen::MasterInfo("Identifying nodes at specially marked corners");

   int compMax;
   int compMin;
   myNumberOfAuxCornerTreatmentNodes_ = 0;

   //only loop over aux participating nodes
   for (int n=0;n<myNumberOfAuxParticipatingNodes_;++n)
     {
       int localNodeID  = myAuxParticipatingNodes_[n];
       int globalNodeID = elemMaps_.nodeLocal2Global[localNodeID];
       compMax = -10000;
       compMin =  10000;

       //loop over extended face list so we can idenpendently determine these corner nodes
       for (int e=0;e<myElementCount_ext_;++e)
	 {
	   //hard coded for triangles
	   if (elementData_ext_[e].nnodes == 3)
	     {	       
	       int nnodesOffset_ext = nnodesLocalOffset_ext_[e];

	       //see which faces in the extended map owns this normal aux node
	       if ( faces_ext_[0 + nnodesOffset_ext] == localNodeID ||
		    faces_ext_[1 + nnodesOffset_ext] == localNodeID || 
		    faces_ext_[2 + nnodesOffset_ext] == localNodeID )
		 {

		   //store max and min of compIDs of all elements in the extended map that touch my aux participating nodes
		   int compID = elementData_ext_[e].compID;
		   if (compID > compMax) compMax = compID;
		   if (compID < compMin) compMin = compID;

		 }
	     }
	 }

       //test if the product of the max and min the element compIDs have the same sign
       int product = compMin*compMax;

       //if the product is negative, the node is on a special boundary
       if (product < 0)
	 {
	   //store the problem node
	   myAuxCornerTreatmentNodes_.push_back(localNodeID);
	   myNumberOfAuxCornerTreatmentNodes_ += 1;
	 }
     }

 }
 
 void fem::gpuPrecomputeCableElements()
 {
   // Batch-compute Ke/Me/Fvec for every CABLE element this rank owns in one
   // GPU launch, cached in gpuCable{Ke,Me,Fvec}_. Called once per assembly
   // pass, BEFORE the per-element loop in assembleGlobKMat(). Everything
   // else (BEAM, SIMPLETRI, MITC3/LMITC3) is untouched and still computed
   // on CPU via comp_Ke_loc/comp_Me_loc inside the loop -- see the dispatch
   // in assembleGlobKMat() below.
   using namespace KARMA;

   gpuCableIndex_.assign(myElementCount_, -1);
   gpuCableKe_.clear();
   gpuCableMe_.clear();
   gpuCableFvec_.clear();

   if (!useGPUAssembly_ || !gpuIsAvailable()) return;

   GPUCableInput in;
   vector<int> localFaceOfCable; // batch index -> local face id f
   for (int f = 0; f < myElementCount_; ++f)
     {
       if (elementData_[f].type != CABLE) continue;

       int nnodesOffset = nnodesLocalOffset_[f];
       int fn0 = faces_[nnodesOffset + 0];
       int fn1 = faces_[nnodesOffset + 1];

       in.fn0.push_back(fn0);
       in.fn1.push_back(fn1);
       for (int idir = 0; idir < xyz; ++idir)
         {
           in.x0_0.push_back(nodes0_[fn0*xyz+idir]);
           in.x1_0.push_back(nodes0_[fn1*xyz+idir]);
           in.x0_t.push_back(nodes_ [fn0*xyz+idir]);
           in.x1_t.push_back(nodes_ [fn1*xyz+idir]);
         }
       in.eModulus.push_back(eModulus_[f]);
       in.density.push_back(density_[f]);
       in.area.push_back(Area_[f]);

       gpuCableIndex_[f] = (int)localFaceOfCable.size();
       localFaceOfCable.push_back(f);
     }
   in.nCable = (int)localFaceOfCable.size();
   if (in.nCable == 0) return;

   GPUCableOutput out;
   try
     {
       gpuComputeCableElements(in, out);
     }
   catch (const std::exception& ex)
     {
       // Fail loud, fail to CPU. Never silently produce a partially-GPU,
       // partially-garbage assembly.
       Screen::MasterInfo("GPU cable assembly failed (" + string(ex.what()) +
                           "); falling back to CPU for this assembly pass.");
       gpuCableIndex_.assign(myElementCount_, -1);
       return;
     }

   gpuCableKe_   = out.Ke;
   gpuCableMe_   = out.Me;
   gpuCableFvec_ = out.Fvec;
   // out.tension currently unused by the caller; wire into
   // myCableTensionVector_[f] here if/when downstream code needs it read
   // back per-element -- same indexing as gpuCableIndex_.
 }

 void fem::assembleGlobKMat()
 {
#if (USEPETSC>0)
   //assemble global stiffness, mass, and LHS matrices
   Screen::MasterInfo("Assembling global FEM matrices");

   gpuPrecomputeCableElements();

#if (TIME_ASSEMBLY == 1)
   ierr = MPI_Barrier(MPI_COMM_WORLD);
   clock_t start;
   double duration;
   if (mypeno_ == 0) start = clock();
#endif

   double invdt2      = 1./(dt_*dt_);
   double epsilon     = 1.E-10;
   double nmb_factor  = 1./nm_beta_*invdt2;
   double nmb_factor2 = nm_gamma_/(nm_beta_*dt_);
   fill(vonMisesStress_.begin(),vonMisesStress_.end(),0.);

   for (int f=0;f<myElementCount_;++f)
     {
       int fGlob        = elemMaps_.elemLocal2Global[f];
       int nnodes_loc   = elementData_[f].nnodes;
       int elType       = elementData_[f].type;
       int ndof_elem    = elementData_[f].ndof;
       int nnodesOffset = nnodesLocalOffset_[f];
       //zero out last data
       fill(Ke_loc_.begin(),Ke_loc_.end(),0.);
       fill(Me_loc_.begin(),Me_loc_.end(),0.);
       fill(Fvec_assemble_.begin(),Fvec_assemble_.end(),0.);
       for (int i=0;i<3;++i) for (int j=0;j<3;++j) cauchyStressTensor_[f][i][j] = 0.;
       //compute local stiffness and mass matrices --
       //GPU path for CABLE elements when available (see gpuPrecomputeCableElements()),
       //CPU path (original, unchanged) for everything else.
       if (elType == CABLE && gpuCableIndex_[f] >= 0)
         {
           int gi = gpuCableIndex_[f];
           for (int i=0;i<36;++i) Ke_loc_[i] = gpuCableKe_[gi*36+i];
           for (int i=0;i<36;++i) Me_loc_[i] = gpuCableMe_[gi*36+i];
           for (int i=0;i<6; ++i) Fvec_assemble_[i] = gpuCableFvec_[gi*6+i];
         }
       else
         {
           comp_Ke_loc(Ke_loc_,Fvec_assemble_,f,elType);
           comp_Me_loc(Me_loc_,f,elType);
         }

       //damping
       double alphadamp = 0.0;
       double betadamp  = 0.0;
       if (elType == 5)
         {
           for (int i=0;i<6;++i) Ke_loc_[i+i*6] += 1.E-10;
           alphadamp = 0.0;
           betadamp  = 0.0;
         }

       vector<int> nodes(nnodes_loc,0);
       for (int n=0;n<nnodes_loc;++n) 
         {
           int index = n + nnodesOffset;
           nodes[n]  = faces_[index];
         }
       //assemble total Fvec_ for each proc
       for (int nn=0;nn<nnodes_loc;++nn)
         {
           int nodeID = nodes[nn];               //local node ID
           int index  = ndofLocalOffset_[nodeID];//local ndof offset of local node
           for (int d=0;d<ndof_elem;++d)
             {
               Fvec_[d + index] += Fvec_assemble_[d + nn*ndof_elem];
             }
         }
       //linear/nonlinear Newmark-Beta time integration
       if (is_unsteady_ && timeIntegrationMethod_ == "implicit")
         {
           for (int n1=0;n1<nnodes_loc;++n1)
             {
               for (int n2=0;n2<nnodes_loc;++n2)
                 {

                   PetscInt rowJumper = 0;
                   int indx2[ndof_elem];
                   int loc1 = elemMaps_.nodeLocal2Global[nodes[n1]];
                   int loc2 = elemMaps_.nodeLocal2Global[nodes[n2]];

                   //fill vector representations
                   for (int d1=0;d1<ndof_elem;++d1)
                     {
                       for (int d2=0;d2<ndof_elem;++d2)
                         {

                           int index = d2 + n2*ndof_elem + (nnodes_loc*ndof_elem)*d1 + n1*nnodes_loc*ndof_elem*ndof_elem;
                           indx2     [d2] = d2 + ndofGlobalOffset_[nodes[n2]];
                           stiff_vec_[d2] = Ke_loc_[index];
                           mass_vec_ [d2] = Me_loc_[index];
                           damp_vec_ [d2] = alphadamp*Me_loc_[index] + betadamp*Ke_loc_[index];
                           unmat_vec_[d2] = stiff_vec_[d2] + nmb_factor*mass_vec_[d2] + nmb_factor2*damp_vec_ [d2];

                         }
                       rowJumper = d1 + ndofGlobalOffset_[nodes[n1]];
                       MatSetValues(Kmat_p_,1,&rowJumper,ndof_elem,indx2,&stiff_vec_[0],ADD_VALUES);
                       MatSetValues(Mmat_p_,1,&rowJumper,ndof_elem,indx2,& mass_vec_[0],ADD_VALUES);
                       MatSetValues(Cmat_p_,1,&rowJumper,ndof_elem,indx2,& damp_vec_[0],ADD_VALUES);
                       MatSetValues( mat_p_,1,&rowJumper,ndof_elem,indx2,&unmat_vec_[0],ADD_VALUES);
                     }
                 }
             }
         }
       //linear/nonlinear steady
       else
         {
           for (int n1=0;n1<nnodes_loc;++n1)
             {
               for (int n2=0;n2<nnodes_loc;++n2)
                 {

                   PetscInt rowJumper = 0;
                   int indx2[ndof_elem];
                   int loc1 = elemMaps_.nodeLocal2Global[nodes[n1]];
                   int loc2 = elemMaps_.nodeLocal2Global[nodes[n2]];

                   //fill vector representations
                   for (int d1=0;d1<ndof_elem;++d1)
                     {
                       for (int d2=0;d2<ndof_elem;++d2)
                         {

                           int index = d2 + n2*ndof_elem + (nnodes_loc*ndof_elem)*d1 + n1*nnodes_loc*ndof_elem*ndof_elem;
                           indx2     [d2] = d2 + ndofGlobalOffset_[nodes[n2]];
                           stiff_vec_[d2] = Ke_loc_[index];

                         }
                       rowJumper = d1 + ndofGlobalOffset_[nodes[n1]];
                       MatSetValues(Kmat_p_,1,&rowJumper,ndof_elem,indx2,&stiff_vec_[0],ADD_VALUES);
                     }
                 }
             }
         }
     }
   //apply BC to global system
   Screen::MasterInfo("Applying BC to global system");
   MatAssemblyBegin(Kmat_p_,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd  (Kmat_p_,MAT_FINAL_ASSEMBLY);

   MatZeroRows(Kmat_p_,myFixedDOFCount_,&myFixedDOFIndexList_[0],1.0,NULL,NULL);
   if (is_unsteady_)
     {

       MatAssemblyBegin(Mmat_p_,MAT_FINAL_ASSEMBLY);
       MatAssemblyEnd  (Mmat_p_,MAT_FINAL_ASSEMBLY);
       MatAssemblyBegin(Cmat_p_,MAT_FINAL_ASSEMBLY);
       MatAssemblyEnd  (Cmat_p_,MAT_FINAL_ASSEMBLY);
       MatAssemblyBegin( mat_p_,MAT_FINAL_ASSEMBLY);
       MatAssemblyEnd  ( mat_p_,MAT_FINAL_ASSEMBLY);

       MatZeroRows(Mmat_p_,myFixedDOFCount_,&myFixedDOFIndexList_[0],1.0,NULL,NULL);
       MatZeroRows(Cmat_p_,myFixedDOFCount_,&myFixedDOFIndexList_[0],1.0,NULL,NULL);
       MatZeroRows(mat_p_ ,myFixedDOFCount_,&myFixedDOFIndexList_[0],1.0,NULL,NULL);

     }

#if (DUMP_KMAT_FROM_PETSC == 1)
#if (DEDICATED_FEM_PROC == 0)
   if (is_unsteady_)
     {
       if (do_restart_)
	 {
	   ierr = MPI_Barrier(MPI_COMM_WORLD);
	   Screen::MasterInfo("Dumping global LHS matrix --- restart");
	   PetscViewer viewer;
	   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"kmat_restart.dat",&viewer);
	   MatView(Kmat_p_,viewer);
	   wait("dumped");
	 }
       else
	 {
	   ierr = MPI_Barrier(MPI_COMM_WORLD);
	   Screen::MasterInfo("Dumping global mats");
	   PetscViewer viewer;
	   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"kmat_orig.dat",&viewer);
	   MatView(Kmat_p_,viewer);
	   PetscViewer viewer2;
	   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat_orig.dat",&viewer2);
	   MatView(mat_p_,viewer2);
	   wait("dumped");
	 }
     }
   else
     {
       ierr = MPI_Barrier(MPI_COMM_WORLD);
       Screen::MasterInfo("Dumping global stiffness matrix");
       PetscViewer viewer;
       PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.m",&viewer);
       PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
       MatView(Kmat_p_,viewer);

       wait();
     }
#else
   if (is_unsteady_)
     {
       Screen::MasterInfo("Dumping global LHS matrix");
       PetscViewer viewer;
       PetscViewerASCIIOpen(PETSC_COMM_SELF,"kmat.dat",&viewer);
       MatView(Kmat_p_,viewer);
       wait();     
     }
   else
     {
       Screen::MasterInfo("Dumping global stiffness matrix");
       PetscViewer viewer;
       PetscViewerASCIIOpen(PETSC_COMM_SELF,"kmat.dat",&viewer);
       MatView(Kmat_p_,viewer);
       wait();
     }
#endif
#endif

#if (TIME_ASSEMBLY == 1)
   ierr = MPI_Barrier(MPI_COMM_WORLD);
   if (mypeno_ == 0)
     {
       duration = (clock() - start)/(double) CLOCKS_PER_SEC;
       Screen::MasterInfo("Global matrix assembly took "+to_string(duration)+" seconds.");
     }
#endif

   firstTimeComputingK_ = false;

#endif
 }
  
 void fem::setForceSteady()
 {
#if (USEPETSC>0)
   PetscScalar rhs = 0.;
   for (int nn=0;nn<myNodeCount_;++nn)
     {
       int globalNode = elemMaps_.nodeLocal2Global[nn];
       int ndof_loc = ndofPerNode_[nn];

       for (int dof=0;dof<ndof_loc;++dof)
         {
           int localIndex  = ndof_loc*nn + dof;
           int globalIndex = dof + ndof_loc*globalNode;
           rhs = forceVec_[localIndex];
           VecSetValues(loadVector_p_,1,&globalIndex,&rhs,ADD_VALUES);                
         }
     }

#endif
 }
  
 void fem::compLoadVector(int& nt)
 {
#if (USEPETSC>0)
   //if there is no forcing from file, we do not need to keep a history of forceVec --- clear it
   if (noForceFile_) fill(forceVec_.begin(),forceVec_.end(),0.);
   //apply surface forces from geometry
   if (isFSI_ || useAuxGeometry_) auxLoadTransfer();
   //apply pressure and body forces --- these are strictly += routines on forceVec
   if (pressureLoading_) applyPressureLoading();
   if (bodyForces_) applyBodyForces();

   //apply BCs -- no BCs are applied before this
   for (int n=0;n<myNodeCount_;++n)
     {
       int ndofOffset = ndofLocalOffset_[n];
       int ndof_loc   = ndofPerNode_    [n];

       for (int dof=0;dof<ndof_loc;++dof)
         {
           if (nodeBC_[dof + ndofOffset] == 1)
             {
               int index = dof + ndofOffset;
               forceVec_[index] = 0.;
               Fvec_    [index] = 0.;
             }
         }
     }
   //treat nonlinear rhs = R-F out of balance load vector
   if (is_nonlinear_)
     {
       PetscScalar rhs = 0.;
       PetscScalar f   = 0.;
       //load stepping
       if (!is_unsteady_)
         {
           int counter = 0;
           for (int nn=0;nn<myNodeCount_;++nn)
             {
               int ndofOffset = ndofGlobalOffset_[nn];
               int ndof_loc   = ndofPerNode_    [nn];

               for (int dof=0;dof<ndof_loc;++dof)
                 {
                   int globalIndex = dof + ndofOffset;
                   rhs = (forceVec_[counter] * (double(ls_)) / double(numLoadSteps_)) - Fvec_[counter];
		   f = (forceVec_[counter] * (double(ls_)) / double(numLoadSteps_));
                   VecSetValues(loadVector_p_,1,&globalIndex,&rhs,ADD_VALUES);
		   VecSetValues(forceVector_p_,1,&globalIndex,&f,ADD_VALUES);
                   counter += 1;
                 }
             }
         }
       else
         {
           int counter = 0;
           //unsteady out of balance load vector
           for (int nn=0;nn<myNodeCount_;++nn)
             {
               int ndofOffset = ndofGlobalOffset_[nn];
               int ndof_loc   = ndofPerNode_    [nn];                

               for (int dof=0;dof<ndof_loc;++dof)
                 {
                   int globalIndex = dof + ndofOffset;
                   rhs = forceVec_[counter] - Fvec_[counter];
		   f = forceVec_[counter];
                   VecSetValues(loadVector_p_,1,&globalIndex,&rhs,ADD_VALUES);
		   VecSetValues(forceVector_p_,1,&globalIndex,&f,ADD_VALUES);
                   counter += 1;
		   //cout<<rhs<<" ";
                 }
	       //cout<<endl;
             }
         }
     }
   else 
     {
       setForceSteady();
     }

#endif
 }

 void fem::solve(const string PCtype) 
 {
#if (USEPETSC>0)

   surfaceTriangulation surfTri;

   MatAssemblyBegin(mat_p_,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd  (mat_p_,MAT_FINAL_ASSEMBLY);

   MatAssemblyBegin(Kmat_p_,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd  (Kmat_p_,MAT_FINAL_ASSEMBLY);

   MatAssemblyBegin(Mmat_p_,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd  (Mmat_p_,MAT_FINAL_ASSEMBLY);

   int staticdummy = 99;
   compLoadVector(staticdummy);

   VecAssemblyBegin(loadVector_p_);
   VecAssemblyEnd  (loadVector_p_);

#if PETSC_VERSION_LT(3,5,0)
   KSPSetOperators(kspOp_p_,Kmat_p_,Kmat_p_,SAME_NONZERO_PATTERN);
#else
   KSPSetOperators(kspOp_p_,Kmat_p_,Kmat_p_);
#endif
   Screen::MasterInfo("Solving linear system");
   KSPSolve(kspOp_p_,loadVector_p_,deltaUk_p_);

   Screen::MasterInfo("Updating node locations");
   //update node locations
   for (int nn=0;nn<nof_nodes_glob_;++nn) 
     {
       int ndof_loc = ndofPerNode_[nn];
       if (DIM == 3)
         {
           for (int ndof=0;ndof<ndof_loc;++ndof) 
             {                
               PetscScalar disp;
               int dispIndex = ndof_loc*nn + ndof;
               VecGetValues(deltaUk_p_,1,&dispIndex,&disp);
               deltaUi_[nn*ndof_loc+ndof] = disp;
             }
           for (int dir=0;dir<xyz;++dir)
             {
               nodes_[nn*xyz+dir] = deltaUi_[nn*ndof_loc+dir] + nodes0_[nn*xyz+dir];
             }
         }
       else if (DIM == 2)
         {
           for (int dof=0;dof<ndof_loc;++dof) 
             {                 
               PetscScalar disp;
               int dispIndex = ndof_loc*nn + dof;
               VecGetValues(deltaUk_p_,1,&dispIndex,&disp);
               deltaUi_[nn*ndof_loc+dof] = disp;
             }
           nodes_[nn*xyz+1]  = deltaUi_[nn*ndof_loc+0] + nodes0_[nn*xyz+1];
         }        
     }
   //update director vectors for outputting to vtk
   updateDirectorVecs();

   vector<string> varnames(6);
   varnames[0]="u";
   varnames[1]="v";
   varnames[2]="w";
   varnames[3]="Rx";
   varnames[4]="Ry";
   varnames[5]="Rz";

   //surfTri.writeTriangulation("Sol.vtk",nodes_,faces_,DIM);
   PetscViewer viewer;
   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"sol.output",&viewer);
   VecView(deltaUk_p_,viewer);
   MatView(Kmat_p_,viewer);

   if (recordNodeOutput_ && mypeno_ == 0) outputNodeData(staticdummy);
   dumpFEMForces(staticdummy);

#endif
 }
  
 void fem::nonLinearSolve()
 {
#if (USEPETSC>0)

   surfaceTriangulation surfTri;
   itCount_      = 1;
   time_         = 1000.;
   pressureTime_ = 1000.;

   if (mypeno_ == 0) cout << endl;
   Screen::MasterInfo("Beginning Load stepping loop");
   //loop over load steps
   for (ls_=1;ls_<=numLoadSteps_;++ls_) 
     { 

       if (mypeno_ == 0) cout << endl;
       Screen::MasterInfo("Load Step: "+to_string(ls_));
       if (mypeno_ == 0) cout << "==========================================" << endl;
       //loop over nonlinear iterations per load step
       for (it_=0;it_<nonlinearIts_;++it_) 
         {

           if (mypeno_ == 0) cout << endl;
           if (mypeno_ == 0) cout << "Sub Iteration: " << it_+1 << endl;
           //calculate global stiffness matrix K and rhs
           reInit();
           assembleGlobKMat();
           compLoadVector(ls_);

           //test RHS convergence
           VecAssemblyBegin(loadVector_p_);
           VecAssemblyEnd  (loadVector_p_);
           PetscReal rhsNorm = 0.;
           VecNorm(loadVector_p_,NORM_INFINITY,&rhsNorm);
           if (mypeno_ == 0) cout  << "\t" << "RHS inf norm    = " << setprecision(15) << rhsNorm  << endl;
           if (mypeno_ == 0) rhsInf << itCount_ << "\t" << rhsNorm << endl;

           //test displacement convergence
#if (DEDICATED_FEM_PROC == 0)
           VecAssemblyBegin(deltaUk_p_);
           VecAssemblyEnd  (deltaUk_p_);
           PetscReal dispNorm = 0.;
           VecNorm(deltaUk_p_,NORM_INFINITY,&dispNorm);
#else
           VecAssemblyBegin(mySolutionVector_p_);
           VecAssemblyEnd  (mySolutionVector_p_);
           PetscReal dispNorm = 0.;
           VecNorm(mySolutionVector_p_,NORM_INFINITY,&dispNorm);
#endif
           if (mypeno_ == 0) cout << "\t" << "Disp inf norm   = " << setprecision(15) << dispNorm << endl;
           if (mypeno_ == 0) cout << "------------------------------" << endl;
           if (mypeno_ == 0) cout << endl;

           if (rhsNorm > 1.E20 || dispNorm > 1.E20) Screen::MasterError("Norms exceeded limits");
           if (rhsNorm < newtonTolerance_ && dispNorm < newtonTolerance_ && it_ > 0)
             //if (rhsNorm < newtonTolerance_ && dispNorm < newtonTolerance_ && it_ > 0)
             {
               MatAssemblyBegin(Kmat_p_,MAT_FINAL_ASSEMBLY);
               MatAssemblyEnd  (Kmat_p_,MAT_FINAL_ASSEMBLY);
               VecAssemblyBegin(loadVector_p_);
               VecAssemblyEnd  (loadVector_p_);
               //converged
               cout << endl;
               Screen::MasterInfo("Newton-Raphson Iterations Converged After " + to_string(it_+1) + " Iteration(s)");

               if (ls_%ntSkip_ == 0 || recordNodeOutput_)
                 {

                   VecAssemblyBegin(Uk_p_);
                   VecAssemblyEnd  (Uk_p_);
#if (DEDICATED_FEM_PROC == 0)
                   VecScatterBegin(scatterToRoot_p_,Uk_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);
                   VecScatterEnd  (scatterToRoot_p_,Uk_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);
#else
                   VecCopy(Uk_p_,rootTotalSolution_p_);
#endif
                   if (mypeno_ == 0)
                     {
                       PetscScalar disp;
                       double maxDisp = 0.;
                       for (int nn=0;nn<nof_nodes_glob_;++nn)
                         {
                           int ndof_loc = rootNdofPerNode_[nn];
                           for (int idir=0;idir<xyz;++idir)
                             {
                               int dispIndex = idir + rootNdofOffset_[nn];
                               int xyzIndex  = nn*xyz      + idir;
                               VecGetValues(rootTotalSolution_p_,1,&dispIndex,&disp);
                               rootTotalNodes_[xyzIndex] = disp + rootTotalNodes0_[xyzIndex];
                               if (fabs(disp) > maxDisp) maxDisp = fabs(disp);
                             }
                         }
                       cout << "  [nlfem] Max nodal displacement at load step " << ls_ << " = " << maxDisp << endl;
                     }
                 }

               if (mypeno_ == 0 && ls_%ntSkip_ == 0)
                 {
                   Screen::MasterInfo("Writing solution to file from root");

                   ofstream surf_out;
                   string surface_filename;
                   ostringstream OSS;
                   OSS << "post/loadstep/loadStep_" << setw(8) << setfill('0') << ls_ << ".vtk";
                   surface_filename = OSS.str();
                   surf_out.open(surface_filename);
                   surfTri.smartWriteTriangulation(surface_filename,rootTotalNodes_,rootTotalFaces_,xyz,nof_nodes_glob_,nof_faces_glob_,numberOfElementTypes_,rootNnodesPerEl_);
                   surf_out.close();

                 }
               break;
             }
           else
             {
               //not converged
#if (DEDICATED_FEM_PROC == 0)
               //solve for iteration k
               KSPSolve(kspOp_p_,loadVector_p_,deltaUk_p_);
               VecScale(deltaUk_p_,relaxationFactor_);
               VecAXPY(Uk_p_,1.,deltaUk_p_);
               //strategic broadcast --- broadcast selected components of solution vector to all procs so they can grab updates from shared nodes
               VecScatterBegin(scatterToAll_p_,deltaUk_p_,mySolutionVector_p_,INSERT_VALUES,SCATTER_FORWARD);
               VecScatterEnd  (scatterToAll_p_,deltaUk_p_,mySolutionVector_p_,INSERT_VALUES,SCATTER_FORWARD);
#else
               //solve for iteration k
               KSPSolve(kspOp_p_,loadVector_p_,mySolutionVector_p_);                
               VecScale(mySolutionVector_p_,relaxationFactor_);
               //update total displacement vector
               VecAXPY(Uk_p_,1.,mySolutionVector_p_);
#endif
               //update Ui with deltaUi
               for (int nn=0;nn<myNodeCount_;++nn)
                 {
                   int ndofOffset = ndofLocalOffset_[nn];
                   int ndof_loc   = ndofPerNode_    [nn];

                   for (int ndof=0;ndof<ndof_loc;++ndof) 
                     {
                       PetscScalar disp;
                       int localIndex  = ndof + ndofOffset;
                       VecGetValues(mySolutionVector_p_,1,&localIndex,&disp);

                       deltaUi_[localIndex]  = disp;
                       Ui_     [localIndex] += deltaUi_[localIndex];
                       Uim1_   [localIndex]  = Ui_[localIndex] - deltaUi_[localIndex];
                     }
                   //update node locations with each sub-it solution
                   for (int idir=0;idir<xyz;++idir) 
                     {
                       nodes_[nn*xyz+idir] = Ui_[idir + ndofOffset] + nodes0_[nn*xyz+idir];
                     }
                 }
             }
           itCount_ += 1;
           updateDirectorVecs();
           if (it_ == nonlinearIts_-1) Screen::MasterError("Newton method did not converge in "+to_string(nonlinearIts_)+" iterations");
         }//end nonlinear iterations
       if (recordNodeOutput_ && mypeno_ == 0) outputNodeData(ls_);
       it_ = 0;
       if (mypeno_ == 0) cout << "==========================================" << endl;
     }//end loadsteps

#endif
 }

 void fem::compBetaNewmarkRHS()
 {
#if (USEPETSC>0)
   PetscReal factor_p  = 4./(dt_*dt_);
   PetscReal factor2_p = 4./(dt_);
   //create rhs for Beta-Newmark time integration in "parts"
   Screen::MasterInfo("Creating Newmark-Beta RHS Vector");

   //create rhs for Beta-Newmark time integration
   //(R-F) - M *Uddot
   MatMult(Mmat_p_,Uddotp1_p_,tempVec2_p_);
   VecWAXPY(rhs_p_,-1.,tempVec2_p_,loadVector_p_);

   //finish assembling Petsc vars before solve
   VecAssemblyBegin(rhs_p_);
   VecAssemblyEnd  (rhs_p_);
   VecAssemblyBegin(loadVector_p_);
   VecAssemblyEnd  (loadVector_p_); 

   MatAssemblyBegin(Kmat_p_,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd  (Kmat_p_,MAT_FINAL_ASSEMBLY);
   MatAssemblyBegin(Mmat_p_,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd  (Mmat_p_,MAT_FINAL_ASSEMBLY);
   MatAssemblyBegin(mat_p_ ,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd  (mat_p_ ,MAT_FINAL_ASSEMBLY);

#endif
 }

 void fem::updateBetaNewmarkAccel()
 {
#if (USEPETSC>0)
   PetscReal factor_p  = 4./(dt_*dt_);
   PetscReal factor2_p = 4./(dt_);
   //create rhs for Beta-Newmark time integration in "parts"

   //4/dt^2*(Uk - tU + dektaUk) -- (hold tempVec)
   VecWAXPY(tempVec_p_,-1.,U_p_,Uk_p_);
   VecSet(avec_p_,factor_p);
   VecPointwiseMult(tempVec2_p_,avec_p_,tempVec_p_);

   //tempVec - 4/dt * tUdot -- (hold tempVec3)
   VecSet(avec_p_,factor2_p);
   VecPointwiseMult(tempVec_p_,avec_p_,Udot_p_);
   VecWAXPY(tempVec3_p_,-1.,tempVec_p_,tempVec2_p_);

   //tempVec3 - tUddot -- (hold Uddot)
   VecWAXPY(Uddotp1_p_,-1.,Uddot_p_,tempVec3_p_);

#endif
 }
  
 void fem::compBetaNewmarkRHS_updated()
 {
#if (USEPETSC>0)
   //(R-F) - M*Uddot - C*udot
   VecAssemblyBegin(loadVector_p_);
   VecAssemblyEnd  (loadVector_p_);
   MatMult(Mmat_p_,Uddotp1_p_,tempVec2_p_);
   MatMult(Cmat_p_,Udotp1_p_  ,tempVec_p_ );
   VecWAXPY(tempVec3_p_,-1.,tempVec2_p_,loadVector_p_);
   VecWAXPY(rhs_p_,-1.,tempVec_p_,tempVec3_p_);
   VecAssemblyBegin(rhs_p_);
   VecAssemblyEnd  (rhs_p_);

//   for (int nn=0;nn<myNodeCount_;++nn)
//     {
//       int ndofOffset = ndofLocalOffset_[nn];
//       int ndof_loc   = ndofPerNode_    [nn];
//       for (int ndof=0;ndof<ndof_loc;++ndof) 
//       	 {
//       	   PetscScalar udot;
//       	   int dispIndex = ndof + ndofOffset;
//       	   VecGetValues(Uk_p_,1,&dispIndex,&udot);
//       	   cout << udot << " ";
//       	 }
//       cout << " ||| ";
//       for (int ndof=0;ndof<ndof_loc;++ndof) 
//	 {
//	   PetscScalar udot;
//	   int dispIndex = ndof + ndofOffset;
//	   VecGetValues(Udotp1_p_,1,&dispIndex,&udot);
//	   cout << udot << " ";
//	 }
//       cout << endl;
//     }

#endif
 }
  
 void fem::updateBetaNewmarkAccel_updated()
 {
#if (USEPETSC>0)
   //compute acceleration at iteration k
   PetscReal factor1_p = 1./(dt_*dt_*nm_beta_);
   PetscReal factor2_p = 1./(dt_*nm_beta_);
   PetscReal factor3_p = 1./(2.*nm_beta_)-1.;

   //factor1_p*(Uk - tU) -- (hold tempVec2)
   VecWAXPY(tempVec_p_,-1.,U_p_,Uk_p_);
   VecSet(avec_p_,factor1_p);
   VecPointwiseMult(tempVec2_p_,avec_p_,tempVec_p_);

   //tempVec - factor2_p * tUdo -- (hold tempVec3)
   VecSet(avec_p_,factor2_p);
   VecPointwiseMult(tempVec_p_,avec_p_,Udot_p_);
   VecWAXPY(tempVec3_p_,-1.,tempVec_p_,tempVec2_p_);

   //tempVec3 - factor3_p * tUddot -- (hold Uddotp1)
   VecSet(avec_p_,factor3_p);
   VecPointwiseMult(tempVec2_p_,avec_p_,Uddot_p_);
   VecWAXPY(Uddotp1_p_,-1.,tempVec2_p_,tempVec3_p_);

#endif
 }

 void fem::updateBetaNewmarkVel_updated()
 {
#if (USEPETSC>0)
   //compute velocity at iteration k
   PetscReal factor1_p = 1. - nm_gamma_/nm_beta_;
   PetscReal factor2_p = nm_gamma_/(nm_beta_*dt_);

   //factor2_p*(Uk - tU) -- (hold tempVec2)
   VecWAXPY(tempVec_p_,-1.,U_p_,Uk_p_);
   VecSet(avec_p_,factor2_p);
   VecPointwiseMult(tempVec2_p_,avec_p_,tempVec_p_);

   //(1-factor1_p)* tUdot -- (hold tempVec)
   VecSet(avec_p_,factor1_p);
   VecPointwiseMult(tempVec_p_,avec_p_,Udot_p_);

   //add them together -- (hold Udotp1_p_)
   VecWAXPY(Udotp1_p_,1.,tempVec_p_,tempVec2_p_);

#endif
 }

 void fem::BetaNewmarkTrapezoidal(int& nt)
 {
#if (USEPETSC>0)
   surfaceTriangulation surfTri;

   nt_ = nt;
   cout << endl;
   time_ = nt_*dt_;
   //Screen::MasterInfo("SIMULATION TIME, t = "+to_string(time_));
   Screen::MasterInfo("Timestep:   "+to_string(nt+1)+"/"+to_string(numTimesteps_));
   cout << "==========================================" << endl;

   //loop over nonlinear iterations per timestep
   for (it_=0;it_<nonlinearIts_;++it_)         
     { 
       cout << endl;
       cout << "Sub Iteration: " << it_ << endl;
       //test Newton-Raphson convergence
       PetscReal RHSnorm    = 10.;
       PetscReal ETOL       = 10.;
       double Fvecnorm      = 10.;
       double Unorm         = 10.;
       //reset varibles to zero
       reInit();
       //create global stiffness and mass matrices
       assembleGlobKMat(); 
       //compute out-of-balance load vector R-F
       compLoadVector(nt);
       //get Newton-Raphson rhs infinity norms 
       checkNewtonConvergence(RHSnorm,Fvecnorm,rhsInf,Unorm,ETOL);

       if (RHSnorm < newtonTolerance_ && Unorm < newtonTolerance_ && it_ > 0)
	 {
	   //converged
	   MatAssemblyBegin(Kmat_p_,MAT_FINAL_ASSEMBLY);
	   MatAssemblyEnd  (Kmat_p_,MAT_FINAL_ASSEMBLY);
	   MatAssemblyBegin(Mmat_p_,MAT_FINAL_ASSEMBLY);
	   MatAssemblyEnd  (Mmat_p_,MAT_FINAL_ASSEMBLY);
	   MatAssemblyBegin(mat_p_ ,MAT_FINAL_ASSEMBLY);
	   MatAssemblyEnd  (mat_p_ ,MAT_FINAL_ASSEMBLY);
	   VecAssemblyBegin(loadVector_p_);
	   VecAssemblyEnd  (loadVector_p_);

	   cout << endl;
	   Screen::MasterInfo("Newton-Raphson Iterations Converged After " + to_string(it_+1) + " Iteration(s)");

	   if (nt % ntSkip_ == 0)
	     {
	       //output load vector to file
	       dumpFEMForces(nt);
	       //write updated nodal positons to VTK every timestep at equilibrium
	       ofstream surf_out;      
	       string surface_filename;
	       ostringstream OSS;
	       OSS << "post/fem/surf_" << setw(8) << setfill('0') << nt << ".vtk";
	       surface_filename = OSS.str();
	       surf_out.open(surface_filename);
	       //surfTri.writeTriangulation_data(surface_filename,nodes_,faces_,DIM,forceVec_,Ui_,nodeVel_,Area_,Vnt_,V1t_,V2t_);
	       surf_out.close();       
	     }
	   //update acceleration to compute RHS
	   updateBetaNewmarkAccel();
	   //update up1 and udotp1
	   updateBetaNewmarkVars();
	   //finish assembling u vectors, nullvecs, and avec
	   VecAssemblyBegin(U_p_);
	   VecAssemblyEnd(U_p_);
	   VecAssemblyBegin(Udot_p_);
	   VecAssemblyEnd(Udot_p_);
	   VecAssemblyBegin(Uddot_p_);
	   VecAssemblyEnd(Uddot_p_); 
	   VecAssemblyBegin(tempVec_p_);
	   VecAssemblyEnd(tempVec_p_);
	   VecAssemblyBegin(tempVec2_p_);
	   VecAssemblyEnd(tempVec2_p_);
	   VecAssemblyBegin(tempVec3_p_);
	   VecAssemblyEnd(tempVec3_p_);
	   VecAssemblyBegin(avec_p_);
	   VecAssemblyEnd(avec_p_);
	   break;
	 } 
       else 
	 {
	   //update acceleration to compute RHS
	   updateBetaNewmarkAccel();
	   //comp RHS for beta newmark
	   compBetaNewmarkRHS();
	   //not converged
	   Screen::MasterInfo("Newton-Raphson Not Converged");
	   //solve for iteration k (or i)
	   KSPSolve(kspOp_p_,rhs_p_,deltaUk_p_);
	   //update Petsc Uk_p_
	   VecScale(deltaUk_p_,relaxationFactor_);
	   VecAXPY(Uk_p_,1.,deltaUk_p_);

	   //update C++ Ui and deltaUi for nodal update
	   for (int nn=0;nn<nof_nodes_glob_;++nn) 
	     {
	       int ndof_loc = ndofPerNode_[nn];
	       for (int ndof=0;ndof<ndof_loc;++ndof) 
		 {
		   PetscScalar disp;
		   int dispIndex = ndof_loc*nn + ndof;
		   VecGetValues(deltaUk_p_,1,&dispIndex,&disp);

		   deltaUi_[nn*ndof_loc+ndof]  = disp;
		   Ui_     [nn*ndof_loc+ndof] += deltaUi_[nn*xyz+ndof];
		   Uim1_   [nn*ndof_loc+ndof]  = Ui_[nn*ndof_loc+ndof] - deltaUi_[nn*ndof_loc+ndof];
		 }
	       //update node locations with each sub-it solution
	       for (int idir=0;idir<xyz;++idir) 
		 { 
		   nodes_[nn*xyz+idir] = Ui_[nn*ndof_loc+idir] + nodes0_[nn*xyz+idir];
		 }
	     }
	   //dumpMat(nodes_,"nodes");
	   //wait();
	   //dumpMat(Ui_,"Ui_");
	   //wait();
	 }
       //update director vectors
       updateDirectorVecs();
       itCount_ += 1;
       if (it_ == nonlinearIts_-1) Screen::MasterError("Newton method did not converge in "+to_string(nonlinearIts_)+" iterations");
     }//end nonlinear iterations

   //output node data
   if (recordNodeOutput_ && mypeno_ == 0) outputNodeData(nt);
   //transfer displacements to geometry
   if (isFSI_ || useAuxGeometry_) auxDispTransfer();
   //dumpMat(Ui_,"Ui_");
   //wait();    
   //reset counters
   it_ = 0;
   cout << "==========================================" << endl;

#endif
 }

#if (USEPETSC>0)
 
 void fem::checkNewtonConvergence(PetscReal& RHSnorm,
                                  double& Fvecnorm,
                                  ofstream& rhsInf,
                                  double& Unorm,
                                  PetscReal& ETOL)
 {
   //compute norm of RHS and Fvec
   VecNorm(rhs_p_,NORM_INFINITY,&RHSnorm);
   Fvecnorm = abs(Fvec_[0]);
   for (int i=1;i<Fvec_.size();++i) if (abs(Fvec_[i]) > Fvecnorm) Fvecnorm = abs(Fvec_[i]);

   //test displacement convergence
   Unorm = abs(Ui_[0] - Uim1_[0]);
   for (int i=0;i<myNodeCount_;++i)
     {
       for (int j=0;j<ndofPerNode_[i];++j)
         {                    
           if (abs(Ui_[i*ndofPerNode_[i]+j]-Uim1_[i*ndofPerNode_[i]+j]) > Unorm) 
             {
               Unorm = abs(Ui_[i*ndofPerNode_[i]+j]-Uim1_[i*ndofPerNode_[i]+j]);
             }
         }
     }

   cout << "\t" << "RHS inf norm    = " << setprecision(15) << RHSnorm  << endl;
   cout << "\t" << "Fvec            = " << setprecision(15) << Fvecnorm << endl;
   cout << "\t" << "Disp inf norm   = " << setprecision(15) << Unorm << endl;

 }

#else

 void fem::checkNewtonConvergence(double& RHSnorm,
                                  double& Fvecnorm,
                                  ofstream& rhsInf,
                                  double& Unorm,
                                  double& ETOL)
 {
   // DUMMY ROUTINE
 }

#endif
 
 void fem::BetaNewmarkGeneral(int& nt)
 {
#if (USEPETSC>0)

   double invdt = 1./dt_;
   surfaceTriangulation surfTri;
   //double relaxationFactor_local = relaxationFactor_;

   if (firstTimeStep_)
     {
       for (int f=0;f<myElementCount_;++f) if (elementData_[f].type == MITC3) density_[f] *= rampFactor_;
     }
   else
     {
       if (rampFactor_ > 1.)
         {
           for (int f=0;f<myElementCount_;++f)
             {
               if (elementData_[f].type == MITC3)
                 {
                   density_[f] *= 0.95;
                   density_[f]  = max(density0_[f],density_[f]);
                 }
             }
         }
     }

   if (mypeno_ == 0) cout << endl;
   Screen::MasterInfo("Timestep:   "+to_string(nt)+"/"+to_string(numTimesteps_+restartTimestep_)+", dt = "+to_string(dt_)+", time = "+to_string(time_));
   if (mypeno_ == 0) cout << "==========================================" << endl;

   //loop over nonlinear iterations per timestep
   for (it_=1;it_<=nonlinearIts_+1;++it_)         
     { 

       if (mypeno_ == 0) cout << "Sub Iteration: " << it_ << endl;
       //reset varibles to zero
       reInit();
       //create global stiffness and mass matrices
       assembleGlobKMat();
       //compute RHS of system
       compLoadVector(nt);
       compBetaNewmarkRHS_updated();

       //initialize terms for relative error norms
       if (it_ == 1)
         {
           R0 = 1.;
           D0 = 1.;
         }

       //back up norms for relaxation
       double rhsNorm_m1 = rhsNorm_;
       double dispNorm_m1 = dispNorm_;

       //test RHS convergence
       VecAssemblyBegin(rhs_p_);
       VecAssemblyEnd  (rhs_p_);
       VecNorm(rhs_p_,NORM_INFINITY,&rhsNorm_);
       //test displacement convergence
#if (DEDICATED_FEM_PROC == 0)
       VecAssemblyBegin(deltaUk_p_);
       VecAssemblyEnd  (deltaUk_p_);
       VecNorm(deltaUk_p_,NORM_INFINITY,&dispNorm_);
#else
       VecAssemblyBegin(mySolutionVector_p_);
       VecAssemblyEnd  (mySolutionVector_p_);
       VecNorm(mySolutionVector_p_,NORM_INFINITY,&dispNorm_);
#endif

       if (it_ == nonlinearIts_ + 1)
         {
           rhsNorm_ = 0.;
           dispNorm_ = 0.;
           R0 =  1.;
           D0 =  1.;
         }

       if (it_ == nonlinearIts_)
	 {
	   rhsNorm_ = 0.;
	   dispNorm_ = 0.;
	   Screen::MasterWarning("Newton method did not converge in "+to_string(nonlinearIts_)+" iterations --- Pushing forward");
	 }

       if (it_ == 2)
         {
           R0 = rhsNorm_;
           D0 = dispNorm_;
         }

       if (it_ > 2)
         {
           if ( (rhsNorm_/R0 > rhsNorm_m1/R0) || (dispNorm_ > dispNorm_m1) )
             {
               relaxationFactor_ = max(0.1,relaxationFactor_*0.85);
             }
           else
             {
               relaxationFactor_ = min(1.0,relaxationFactor_*1.05);
             }
         }

       if (mypeno_ == 0)
         {
           cout << "\t" << "Abs. RHS inf norm        = " << setprecision(15) << rhsNorm_     << endl;
           cout << "\t" << "Abs. RHS inf norm (k-1)  = " << setprecision(15) << rhsNorm_m1   << endl;
           cout << "\t" << "Rel. RHS inf norm        = " << setprecision(15) << rhsNorm_/R0  << endl;
           cout << "\t" << "Abs. Disp inf norm       = " << setprecision(15) << dispNorm_    << endl;
           cout << "\t" << "Abs. Disp inf norm (k-1) = " << setprecision(15) << dispNorm_m1  << endl;
           cout << "\t" << "Rel. Disp inf norm       = " << setprecision(15) << dispNorm_/D0 << endl;
           cout << "\t" << "Relaxation Factor        = " << setprecision(15) << relaxationFactor_ << endl;
           cout << "------------------------------" << endl;
           cout << endl;
         }

       if (rhsNorm_ > 1.E20 || dispNorm_ > 1.E20 || dispNorm_ != dispNorm_ || rhsNorm_!= rhsNorm_ || D0 != D0 || R0 != R0) Screen::MasterError("Norms exceeded limits or NaN");

       if (rhsNorm_/R0 < newtonTolerance_ && dispNorm_ < 5.E-4 && it_ > 1)
         {
           //converged
           if (mypeno_ == 0) cout << endl;
           Screen::MasterInfo("Newton-Raphson Iterations Converged After " + to_string(it_) + " Iteration(s)");

           if (nt % ntSkip_ == 0 || recordNodeOutput_)
             {
               VecAssemblyBegin(Uk_p_);
               VecAssemblyEnd  (Uk_p_);
#if (DEDICATED_FEM_PROC == 1)
               VecCopy(Uk_p_,rootTotalSolution_p_);
#else
               VecScatterBegin(scatterToRoot_p_,Uk_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);
               VecScatterEnd  (scatterToRoot_p_,Uk_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);
#endif

               if (mypeno_ == 0)
                 {
                   PetscScalar disp;
                   for (int nn=0;nn<nof_nodes_glob_;++nn)
                     {
                       int ndof_loc_root = rootNdofPerNode_[nn];
                       for (int idir=0;idir<xyz;++idir)
                         {
                           int dispIndex = idir + rootNdofOffset_[nn];
                           int xyzIndex  = nn*xyz + idir;
                           VecGetValues(rootTotalSolution_p_,1,&dispIndex,&disp);
                           rootTotalNodes_[xyzIndex] = disp + rootTotalNodes0_[xyzIndex];
                         }
                     }
                 }
             }

#if (DEDICATED_FEM_PROC == 0)

           vector<double> rootVonMises;
           vector<double> rootForces;
           vector<double> rootDisps;
           vector<int> orderedList;
	   vector<double> faceDisps(xyz*nof_faces_glob_,0.);
	   vector<double> faceVonMises(nof_faces_glob_,0.);
	   vector<double> faceDisps2(xyz*nof_faces_glob_,0.);
	   vector<double> faceVonMises2(nof_faces_glob_,0.);

           if (nt % ntSkip_ == 0)
             {

	       Screen::MasterInfo("Gathering on root for IO");

               //move displacements to element centers
               for (int e=0;e<myElementCount_;++e)
                 {
                   int nnodes_loc    = elementData_[e].nnodes;
                   int nnodes_offset = nnodesLocalOffset_[e];
		   int eGlob         = elemMaps_.elemLocal2Global[e];
                   for (int n=0;n<nnodes_loc;++n)
                     {
                       int node           = faces_[n + nnodes_offset];
                       int ndofOffset_loc = ndofLocalOffset_[node];
                       faceDisps[eGlob*xyz + 0] += 0.33333333333333*Ui_[ndofOffset_loc + 0];
                       faceDisps[eGlob*xyz + 1] += 0.33333333333333*Ui_[ndofOffset_loc + 1];
                       faceDisps[eGlob*xyz + 2] += 0.33333333333333*Ui_[ndofOffset_loc + 2];
                     }
                 }

               //gather all disps on root for output
	       ierr = MPI_Reduce(&faceDisps[0],&faceDisps2[0],xyz*nof_faces_glob_,MPI::DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
               if (mypeno_ != 0) faceDisps = vector<double>();

	       //gather von Mises on root
	       for (int e=0;e<myElementCount_;++e)
                 {
                   int nnodes_loc    = elementData_[e].nnodes;
                   int nnodes_offset = nnodesLocalOffset_[e];
		   int eGlob         = elemMaps_.elemLocal2Global[e];
		   faceVonMises[eGlob] = vonMisesStress_[e];
		 }
	       ierr = MPI_Reduce(&faceVonMises[0],&faceVonMises2[0],nof_faces_glob_,MPI::DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	       if (mypeno_ != 0) faceVonMises = vector<double>();

             }
#endif

           if (mypeno_ == 0 && nt%ntSkip_ == 0)
             {

	       Screen::MasterInfo("Writing solution to file from root");

               ofstream surf_out;
               string surface_filename;
               ostringstream OSS;
               OSS << "post/fem/surf_" << setw(8) << setfill('0') << nt << ".vtk";
               surface_filename = OSS.str();
               surf_out.open(surface_filename);
#if (DEDICATED_FEM_PROC == 1)
               surfTri.smartWriteTriangulation(surface_filename,rootTotalNodes_,rootTotalFaces_,xyz,nof_nodes_glob_,nof_faces_glob_,numberOfElementTypes_,rootNnodesPerEl_);
#else
               surfTri.smartWriteTriangulation_withForces(surface_filename,rootTotalNodes_,rootTotalFaces_,xyz,nof_nodes_glob_,nof_faces_glob_,numberOfElementTypes_,rootNnodesPerEl_,rootForces,orderedList,faceDisps2,faceVonMises2);
#endif
               surf_out.close();

             }

           //update acceleration to compute RHS
           updateBetaNewmarkAccel_updated();
           //update up1 and udotp1
           updateBetaNewmarkVars_updated();
           //write restart of converged solution on the restartSkip
           if (nt % restartSkip_ == 0) writeRestart(nt);
           break;

         }
       else
         {
           //update acceleration to compute RHS
           updateBetaNewmarkAccel_updated();
           //update velocity to compute RHS
           updateBetaNewmarkVel_updated();
           //comp RHS for beta newmark
           compBetaNewmarkRHS_updated();
           //back up incremental solution vector
#if (DEDICATED_FEM_PROC == 0)
           //solve for iteration k
           KSPSolve(kspOp_p_,rhs_p_,deltaUk_p_);
           //update petsc Uk_p_
           VecScale(deltaUk_p_,relaxationFactor_);
           VecAXPY(Uk_p_,1.,deltaUk_p_);
           //strategic broadcast --- broadcast selected components of solution vector to all procs so they can grab updates from shared nodes
           VecScatterBegin(scatterToAll_p_,deltaUk_p_,mySolutionVector_p_,INSERT_VALUES,SCATTER_FORWARD);
           VecScatterEnd  (scatterToAll_p_,deltaUk_p_,mySolutionVector_p_,INSERT_VALUES,SCATTER_FORWARD);
#else
           //solve for iteration k
           KSPSolve(kspOp_p_,rhs_p_,mySolutionVector_p_);
           //update total displacement vector
           VecScale(mySolutionVector_p_,relaxationFactor_);
           VecAXPY(Uk_p_,1.,mySolutionVector_p_);
#endif
           //fill C++ vars for updates
           for (int nn=0;nn<myNodeCount_;++nn)
             {
	       int ndofOffset = ndofLocalOffset_[nn];
	       int ndof_loc   = ndofPerNode_    [nn];

	       for (int ndof=0;ndof<ndof_loc;++ndof) 
		 {
		   PetscScalar disp;
		   int dispIndex = ndof + ndofOffset;
		   VecGetValues(mySolutionVector_p_,1,&dispIndex,&disp);
		   deltaUi_[dispIndex]  = disp;
		   Ui_     [dispIndex] += deltaUi_[dispIndex];
		   Uim1_   [dispIndex]  = Ui_[dispIndex] - deltaUi_[dispIndex];

		 }
	       //update node locations with each sub-it solution
	       for (int idir=0;idir<xyz;++idir)
		 {
		   nodes_[nn*xyz+idir] = Ui_[idir + ndofOffset] + nodes0_[nn*xyz+idir];
		 }
	     }
         }
       itCount_ +=1;
       //update director vectors
       updateDirectorVecs();

     }//end nonlinear iterations

#if (DEDICATED_FEM_PROC == 0)
   //scatter the converged, extended solution vector
   VecScatterBegin(scatterToAll_p_ext_,Uk_p_,mySolutionVector_p_ext_,INSERT_VALUES,SCATTER_FORWARD);
   VecScatterEnd  (scatterToAll_p_ext_,Uk_p_,mySolutionVector_p_ext_,INSERT_VALUES,SCATTER_FORWARD);
   //update extended node vector
   for (int nn=0;nn<myNodeCount_ext_;++nn)
     {
       int ndofOffset_ext = ndofLocalOffset_ext_[nn];

       for (int idir=0;idir<xyz;++idir) 
	 {
	   PetscScalar disp_ext;
	   int dispIndex_ext = idir + ndofOffset_ext;
	   VecGetValues(mySolutionVector_p_ext_,1,&dispIndex_ext,&disp_ext);
	   nodes_ext_[nn*xyz + idir] = disp_ext + nodes0_ext_[nn*xyz + idir];
	 }
     }
#else
   nodes_ext_.assign(nodes_.begin(),nodes_.end());
#endif
   //update the node normal vectors based on the geometry of the FEM mesh
   if ( (isFSI_ || useAuxGeometry_) && !threeD_as_twoD_) updateGeometricNormals();
   //output node data
   if (recordNodeOutput_ && mypeno_ == 0) outputNodeData(nt);
   //output cable data
   if (nIOCables_ > 0) outputCableData(nt);
   //transfer displacements to geometry
   if (isFSI_ || useAuxGeometry_) auxDispTransfer();
   //reset counters
   it_ = 0;
   if (mypeno_ == 0) cout << "==========================================" << endl;

#endif
 }
 
 void fem::BetaNewmarkGeneral_BacktrackLineSearch(int& nt)
 {
#if (USEPETSC>0)

   //implement a backtrack line search. do not advance the solution if the 
   //relative RHS residual or incremental update increase in INF norm.
   //if an increase is detected, cut the step size and solve again.
   //repeat until the residual drops or a iteration limit is reached.

   double invdt = 1./dt_;

   if (rampFactor_ > 1.)
     {
       if (firstTimeStep_)
	 {
	   for (int f=0;f<myElementCount_;++f) if (elementData_[f].type == MITC3) density_[f] *= rampFactor_;
	 }
       else
	 {
           for (int f=0;f<myElementCount_;++f)
             {
               if (elementData_[f].type == MITC3)
                 {
                   density_[f] *= 0.95;
                   density_[f]  = max(density0_[f],density_[f]);
                 }
             }
         }
     }

   if (mypeno_ == 0) cout << endl;
   Screen::MasterInfo("Timestep:   "+to_string(nt)+"/"+to_string(numTimesteps_+restartTimestep_)+", dt = "+to_string(dt_)+", time = "+to_string(time_));
   if (mypeno_ == 0) cout << "==========================================" << endl;

   relaxationFactor_ = 1.;
   backtrack_ = false;
   advanceSolution_ = true;
   normGrew_ = false;
   int numTries = 1;
   double rhsNorm_m1  = 1.e10;
   double dispNorm_m1 = 1.e10;
   bool hardAdvance = false;

   //loop over nonlinear iterations per timestep
   for (it_=1;it_<=nonlinearIts_+1;++it_)
     { 

       if (!backtrack_) solAdvanced_ = false;

       //reset varibles to zero
       reInit();
       //create global stiffness and mass matrices
       assembleGlobKMat();
       //compute RHS of system
       compLoadVector(nt);
       compBetaNewmarkRHS_updated();

       //initialize terms for relative error norms
       if (it_ == 1)
         {
           R0 = 1.;
           D0 = 1.;
         }

       //test RHS convergence
       VecAssemblyBegin(rhs_p_);
       VecAssemblyEnd  (rhs_p_);
       VecNorm(rhs_p_,NORM_INFINITY,&rhsNorm_);
       //test displacement convergence
#if (DEDICATED_FEM_PROC == 0)
       VecAssemblyBegin(deltaUk_p_);
       VecAssemblyEnd  (deltaUk_p_);
       VecNorm(deltaUk_p_,NORM_INFINITY,&dispNorm_);
#else
       VecAssemblyBegin(mySolutionVector_p_);
       VecAssemblyEnd  (mySolutionVector_p_);
       VecNorm(mySolutionVector_p_,NORM_INFINITY,&dispNorm_);
#endif

       if (it_ == nonlinearIts_ + 1)
	 {
	   rhsNorm_ = 0.;
	   dispNorm_ = 0.;
	   Screen::MasterWarning("Newton method did not converge in "+to_string(nonlinearIts_)+" iterations --- Pushing forward");
	 }

       if (it_ == 2)
         {
           R0 = rhsNorm_;
           D0 = dispNorm_;
         }

       normGrew_ = false;
       if (it_ > 2)
         {
           if ( (rhsNorm_/R0 > rhsNorm_m1/R0) || (dispNorm_ > dispNorm_m1) )
             {
	       normGrew_ = true;
	       if (!backtrack_) it_ -= 1;
             }
           else
             {
	       normGrew_ = false;
             }
	 }
       if (hardAdvance) normGrew_ = false;

       bool backupNorms = false;
       if (backtrack_)
	 {
	   Screen::MasterInfo("Backtrack status - on");
	   //already cut relxation, check norms
	   if (solAdvanced_)
	     {
	       Screen::MasterInfo("Backtrack candidate solution advanced");
	       //if norm grew again, cut relaxation
	       if (normGrew_)
		 {
		   Screen::MasterInfo("Backtrack search grew residual, cutting step size again");
		   relaxationFactor_ = relaxationFactor_*0.5;
		   advanceSolution_ = false;
		   numTries += 1;
		   it_ -= 1;
		 }
	       //cutting relaxation worked, turn off backtrack
	       else
		 {
		   numTries = 1;
		   Screen::MasterInfo("Residual dropped, turning off backtrack, advancing solution");
		   backtrack_ = false;
		   relaxationFactor_ = 1.;
		   advanceSolution_ = true;
		   backupNorms = true;
		 }
	     }
	   else
	     {
	       it_ -= 1;
	       Screen::MasterInfo("Preparing for backtrack candidate solution");
	       advanceSolution_ = true;
	     }
	 }
       else
	 {
	   Screen::MasterInfo("Backtrack status - off");
	   //norm grew, turn on backtrack
	   if (normGrew_)
	     {
	       Screen::MasterInfo("Divergence detected, backtrack search is on");
	       relaxationFactor_ = max(0.01,relaxationFactor_*0.5);
	       backtrack_ = true;
	       advanceSolution_ = false;
	     }
	   else
	     {
	       Screen::MasterInfo("Solution is converging");
	       backupNorms = true;
	     }
	 }

       if (mypeno_ == 0)
         {
	   cout << "Sub Iteration: " << it_ << endl;
	   if (normGrew_) cout << "\t" << "[W] Norm grew!" << endl;
	   if (backtrack_) cout << "\t" << "[W] Number of attemps at iteration " << it_ << ": " << numTries << endl;
           cout << "\t" << "Abs. RHS inf norm        = " << setprecision(15) << rhsNorm_     << endl;
           cout << "\t" << "Abs. RHS inf norm (k-1)  = " << setprecision(15) << rhsNorm_m1   << endl;
           cout << "\t" << "Rel. RHS inf norm        = " << setprecision(15) << rhsNorm_/R0  << endl;
           cout << "\t" << "Abs. Disp inf norm       = " << setprecision(15) << dispNorm_    << endl;
           cout << "\t" << "Abs. Disp inf norm (k-1) = " << setprecision(15) << dispNorm_m1  << endl;
           cout << "\t" << "Rel. Disp inf norm       = " << setprecision(15) << dispNorm_/D0 << endl;
           cout << "\t" << "Relaxation Factor        = " << setprecision(15) << relaxationFactor_ << endl;
	 }

       if (backupNorms)
	 {
	   deltaVnt_bak_.assign(deltaVnt_.begin(),deltaVnt_.end());
	   myDirVecs_bak_.assign(myDirVecs_.begin(),myDirVecs_.end());
	   Vnt_bak_.assign(Vnt_.begin(),Vnt_.end());

	   rhsNorm_m1 = rhsNorm_;
	   dispNorm_m1 = dispNorm_;
	 }

       if (numTries >= 11)
	 {
	   backtrack_ = false;
	   normGrew_ = false;
	   advanceSolution_ = true;
	   Screen::MasterWarning("Could not drop residual after "+to_string(numTries-1)+" attempts, advancing anyway");
	   numTries = 0;
	   hardAdvance = true;
	   relaxationFactor_ = 1.;
	 }
       else
	 {
	   hardAdvance = false;
	 }

       if (rhsNorm_ > 1.E20 || dispNorm_ > 1.E20 || dispNorm_ != dispNorm_ || rhsNorm_!= rhsNorm_ || D0 != D0 || R0 != R0) Screen::MasterError("Norms exceeded limits or NaN");

       if (rhsNorm_/R0 < newtonTolerance_ && dispNorm_/D0 < 100*newtonTolerance_ && it_ > 1)
         {
           //converged
           if (mypeno_ == 0) cout << endl;
           Screen::MasterInfo("Newton-Raphson Iterations Converged After " + to_string(it_) + " Iteration(s)");

	   //gather solution vector on root
           if (nt % ntSkip_ == 0 || recordNodeOutput_)
             {

	       VecAssemblyBegin(forceVector_p_);
	       VecAssemblyEnd  (forceVector_p_);	       
               VecAssemblyBegin(Uk_p_);
               VecAssemblyEnd  (Uk_p_);
#if (DEDICATED_FEM_PROC == 1)
               VecCopy(Uk_p_,rootTotalSolution_p_);
#else
               VecScatterBegin(scatterToRoot_p_,Uk_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);
               VecScatterEnd  (scatterToRoot_p_,Uk_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);
#endif

               if (mypeno_ == 0)
                 {
                   PetscScalar disp;
                   for (int nn=0;nn<nof_nodes_glob_;++nn)
                     {
                       int ndof_loc_root = rootNdofPerNode_[nn];
                       for (int idir=0;idir<xyz;++idir)
                         {
                           int dispIndex = idir + rootNdofOffset_[nn];
                           int xyzIndex  = nn*xyz + idir;
                           VecGetValues(rootTotalSolution_p_,1,&dispIndex,&disp);
                           rootTotalNodes_[xyzIndex] = disp + rootTotalNodes0_[xyzIndex];
                         }
                     }
                 }
             }

	   //dump to vtk or unformatted file
	   if (nt % ntSkip_ == 0)
	     {

	       surfaceTriangulation surfTri;
               ofstream surf_out;
               string surface_filename;
               ostringstream OSS;
#if (DEDICATED_FEM_PROC == 1)
	       //normal .vtk gather and write from a single processor
               OSS << "post/fem/surf_" << setw(8) << setfill('0') << nt << ".vtk";
               surface_filename = OSS.str();	       
               Screen::MasterInfo("Writing solution to file from root");	       
               surf_out.open(surface_filename);
               surfTri.smartWriteTriangulation(surface_filename,rootTotalNodes_,rootTotalFaces_,xyz,nof_nodes_glob_,nof_faces_glob_,numberOfElementTypes_,rootNnodesPerEl_);
	       surf_out.close();
#else
	       //parallel, unformatted IO, hopefully fast
               OSS << "post/fem/surf_unformatted_" << setw(8) << setfill('0') << nt << ".vtk";
               surface_filename = OSS.str();

	       //prepare output vars for parallel IO
	       vector<double> localFaceForces(myElementCount_*xyz,0.);
	       vector<double> localFaceDispls(myElementCount_*xyz,0.);
	       vector<int>    nnodesVec      (myElementCount_*xyz,0);
	       prepareVarsForParallelIO(localFaceForces,localFaceDispls,nnodesVec);
               surfTri.writeTriangulation_parallel_unformatted(surface_filename,nof_nodes_glob_,nof_faces_glob_,localFaceForces,localFaceDispls,elemMaps_.nodeLocal2Global,elemMaps_.elemLocal2Global,myNodeCount_,myElementCount_,nodes_,faces_,nnodesVec,vonMisesStress_);

#endif
	     }

           //update acceleration to compute RHS
           updateBetaNewmarkAccel_updated();
           //update up1 and udotp1
           updateBetaNewmarkVars_updated();
           //write restart of converged solution on the restartSkip
           if (nt % restartSkip_ == 0) writeRestart(nt);
           break;

         }
       else
         {

	   if (advanceSolution_ == false)
	     {

	       Screen::MasterInfo("Stripping update");

	       //residual grew, we need to strip the update and solve again
	       solAdvanced_ = false;

	       //strip incremental update from total displacement vector and normal vector update
	       VecAXPY(Uk_p_,-1.,deltaUk_p_);
	       VecCopy(deltaUkm1_p_,deltaUk_p_);
	       deltaVnt_.assign(deltaVnt_bak_.begin(),deltaVnt_bak_.end());

	       //strip the update from the director vectors
	       unUpdateDirectorVecs();	       

	       //strip C++ updates
	       for (int nn=0;nn<myNodeCount_;++nn)
		 {
		   int ndofOffset = ndofLocalOffset_[nn];
		   int ndof_loc   = ndofPerNode_    [nn];

		   for (int ndof=0;ndof<ndof_loc;++ndof) 
		     {
		       PetscScalar disp;
		       PetscScalar disp_last;
		       int dispIndex = ndof + ndofOffset;
		       VecGetValues(mySolutionVector_p_,1,&dispIndex,&disp);
		       VecGetValues(myLastSolutionVector_p_,1,&dispIndex,&disp_last);
		       deltaUi_[dispIndex]  = disp_last;//revert to last 'good' incremental solution
		       Ui_     [dispIndex] -= disp;//strip off most recent update
		       Uim1_   [dispIndex]  = Ui_[dispIndex] - deltaUi_[dispIndex];//recompute Uim1

		     }
		   //update node locations with each sub-it solution
		   for (int idir=0;idir<xyz;++idir)
		     {
		       nodes_[nn*xyz+idir] = Ui_[idir + ndofOffset] + nodes0_[nn*xyz+idir];
		     }
		 }
	       itCount_ -=1;

	       VecCopy(Udotp1_km1_p_ ,Udotp1_p_);
	       VecCopy(Uddotp1_km1_p_,Uddotp1_p_);

	     }
	   else
	     {

	       //residual is converging, advance solution
	       Screen::MasterInfo("Solving system");

	       //back up last computed velocity and acceleration that was dropping the residual
	       if (!backtrack_)
		 {
		   VecCopy(Udotp1_p_ ,Udotp1_km1_p_);
		   VecCopy(Uddotp1_p_,Uddotp1_km1_p_);
		 }
	       //update acceleration to compute RHS
	       updateBetaNewmarkAccel_updated();
	       //update velocity to compute RHS
	       updateBetaNewmarkVel_updated();
	       //comp RHS for beta newmark
	       compBetaNewmarkRHS_updated();

#if (DEDICATED_FEM_PROC == 0)
	       //back up solution vectors if we are not in a backtrack loop and incremental normal update
	       if (!backtrack_)
		 {
		   VecCopy(deltaUk_p_,deltaUkm1_p_);
		   VecCopy(mySolutionVector_p_,myLastSolutionVector_p_);
		 }
	       //solve for iteration k
	       KSPSolve(kspOp_p_,rhs_p_,deltaUk_p_);
	       //update petsc Uk_p_
	       VecScale(deltaUk_p_,relaxationFactor_);
	       VecAXPY(Uk_p_,1.,deltaUk_p_);
	       //strategic broadcast --- broadcast selected components of solution vector to all procs so they can grab updates from shared nodes
	       VecScatterBegin(scatterToAll_p_,deltaUk_p_,mySolutionVector_p_,INSERT_VALUES,SCATTER_FORWARD);
	       VecScatterEnd  (scatterToAll_p_,deltaUk_p_,mySolutionVector_p_,INSERT_VALUES,SCATTER_FORWARD);
#else
	       //solve for iteration k
	       KSPSolve(kspOp_p_,rhs_p_,mySolutionVector_p_);
	       //update total displacement vector
	       VecScale(mySolutionVector_p_,relaxationFactor_);
	       VecAXPY(Uk_p_,1.,mySolutionVector_p_);
#endif
	       //fill C++ vars for updates
	       for (int nn=0;nn<myNodeCount_;++nn)
		 {
		   int ndofOffset = ndofLocalOffset_[nn];
		   int ndof_loc   = ndofPerNode_    [nn];

		   for (int ndof=0;ndof<ndof_loc;++ndof) 
		     {
		       PetscScalar disp;
		       int dispIndex = ndof + ndofOffset;
		       VecGetValues(mySolutionVector_p_,1,&dispIndex,&disp);
		       deltaUi_[dispIndex]  = disp;
		       Ui_     [dispIndex] += deltaUi_[dispIndex];
		       Uim1_   [dispIndex]  = Ui_[dispIndex] - deltaUi_[dispIndex];

		     }
		   //update node locations with each sub-it solution
		   for (int idir=0;idir<xyz;++idir)
		     {
		       nodes_[nn*xyz+idir] = Ui_[idir + ndofOffset] + nodes0_[nn*xyz+idir];
		     }
		 }
	       itCount_ +=1;

	       //update director vectors
	       updateDirectorVecs();
	       solAdvanced_ = true;
	     }
	 }
       if (mypeno_ == 0)
	 {
           cout << "------------------------------" << endl;
           cout << endl;
	 }

     }//end nonlinear iterations

#if (DEDICATED_FEM_PROC == 0)
   //scatter the converged, extended solution vector
   VecScatterBegin(scatterToAll_p_ext_,Uk_p_,mySolutionVector_p_ext_,INSERT_VALUES,SCATTER_FORWARD);
   VecScatterEnd  (scatterToAll_p_ext_,Uk_p_,mySolutionVector_p_ext_,INSERT_VALUES,SCATTER_FORWARD);
   //update extended node vector
   for (int nn=0;nn<myNodeCount_ext_;++nn)
     {
       int ndofOffset_ext = ndofLocalOffset_ext_[nn];

       for (int idir=0;idir<xyz;++idir) 
	 {
	   PetscScalar disp_ext;
	   int dispIndex_ext = idir + ndofOffset_ext;
	   VecGetValues(mySolutionVector_p_ext_,1,&dispIndex_ext,&disp_ext);
	   nodes_ext_[nn*xyz + idir] = disp_ext + nodes0_ext_[nn*xyz + idir];
	 }
     }
#else
   nodes_ext_.assign(nodes_.begin(),nodes_.end());
#endif
   //update the node normal vectors based on the geometry of the FEM mesh
   if ( (isFSI_ || useAuxGeometry_) && !threeD_as_twoD_) updateGeometricNormals();
   //output node data
   if (recordNodeOutput_ && mypeno_ == 0) outputNodeData(nt);
   //output cable data
   if (nIOCables_ > 0) outputCableData(nt);
   //transfer displacements to geometry
   if (isFSI_ || useAuxGeometry_) auxDispTransfer();
   //reset counters
   it_ = 0;
   if (mypeno_ == 0) cout << "==========================================" << endl;

#endif
 }
 
 void fem::updateBetaNewmarkVars()
 {
#if (USEPETSC>0)
   //calculate udot+1 and u+1
   PetscReal factor  = dt_/2.;

   VecSet(avec_p_,factor);
   //for Up1_p_
   //Udot + 4/dt^2 * Uddot + Uddotp1 -- (hold Udotp1)
   VecWAXPY(tempVec_p_,1.,Uddot_p_,Uddotp1_p_);
   VecPointwiseMult(tempVec2_p_,avec_p_,tempVec_p_);
   VecWAXPY(Udotp1_p_,1.,Udot_p_,tempVec2_p_);

   //for Up1_p_
   VecCopy(Uk_p_,Up1_p_);
   VecCopy(Up1_p_,Uk_p_);

   //finish assembling n+1 vectors before writing to file
   VecAssemblyBegin(Udotp1_p_);
   VecAssemblyEnd  (Udotp1_p_);
   VecAssemblyBegin(Uddotp1_p_);
   VecAssemblyEnd  (Uddotp1_p_);

   //write to data files
   ofstream sol_out;
   if (itCount_% ntSkip_ == 0) 
     {
       string filename_output;
       ostringstream sol;
       sol << "post/unsteady/sol_" << setw(8) << setfill('0') << itCount_ << ".dat";
       filename_output = sol.str();
       sol_out.open(filename_output);
     }
   //set n+1 values to n
   PetscScalar u1=0.,u2=0.,u3=0.,u1p1=0.,u2p1=0.,u3p1=0.;
   //this isnt sized right
   vector<vector<double> > restartVec(nof_nodes_glob_*6,vector<double>(6,0.));

   for (int nn=0;nn<nof_nodes_glob_;++nn) 
     {
       int ndof_loc = ndofPerNode_[nn];
       if (itCount_%ntSkip_ == 0) sol_out << setw(10) << nodes_[nn*xyz+0];

       for (int nd=0;nd<ndof_loc;++nd) 
	 {
	   int index=ndof_loc*nn+nd;
	   //get u+1 values
	   VecGetValues(Up1_p_    ,1,&index,&u1p1);
	   VecGetValues(Udotp1_p_ ,1,&index,&u2p1);
	   VecGetValues(Uddotp1_p_,1,&index,&u3p1);
	   //get u values
	   VecGetValues(U_p_    ,1,&index,&u1);
	   VecGetValues(Udot_p_ ,1,&index,&u2);
	   VecGetValues(Uddot_p_,1,&index,&u3);
	   //set n+1 values to n
	   VecSetValues(U_p_    ,1,&index,&u1p1,INSERT_VALUES);
	   VecSetValues(Udot_p_ ,1,&index,&u2p1,INSERT_VALUES);
	   VecSetValues(Uddot_p_,1,&index,&u3p1,INSERT_VALUES);                

	   if (itCount_%restartSkip_ == 0)
	     {
	       //fill restart vector with n and n+1 values
	       restartVec[index][0] = u1;
	       restartVec[index][1] = u2;
	       restartVec[index][2] = u3;
	       restartVec[index][3] = u1p1;
	       restartVec[index][4] = u2p1;
	       restartVec[index][5] = u3p1;              
	     } 
	   //write displacements to solution file
	   if (itCount_%ntSkip_ == 0) sol_out <<  setw(10) << u1p1 << setw(10) << u2p1 << setw(10) << u3p1 << setw(10);
	 }
       if (itCount_%ntSkip_ == 0) sol_out << endl;     
     }
   //write to restart file
   if (itCount_%restartSkip_ == 0) 
     {
       Screen::MasterInfo("Writing to restart file");
       writeLinearRestart(restartVec,nt_);
     }
   //write updated nodal positons
   ofstream surf_out;      
   if (itCount_%ntSkip_ == 0) 
     {
       sol_out.close();              
     }
#endif
 }

 void fem::updateBetaNewmarkVars_updated()
 {
#if (USEPETSC>0)
   //calculate udot+1 and u+1
   //udot+1 factors
   PetscReal factor4_p = nm_gamma_/(dt_*nm_beta_);
   PetscReal factor5_p = nm_gamma_/nm_beta_;
   PetscReal factor6_p = dt_*((1.-nm_gamma_)-(nm_gamma_*(0.5-nm_beta_)/nm_beta_));

   //for (t+1)U, Uk_p_ is continuously updated with deltaUk_

   VecCopy(Uk_p_,Up1_p_);

   //for Udotp1 = tUdot + factor4 * (t+1U -tU) - factor5 * tUdot + factor6 * tUddot

   //dt*gamma *(t+1)Uddot -- hold tempVec
   VecSet(avec_p_,nm_gamma_*dt_);
   VecPointwiseMult(tempVec_p_,avec_p_,Uddotp1_p_);

   //dt*(1-gamma)*tUddot -- hold tempVec2
   VecSet(avec_p_,dt_*(1.-nm_gamma_));
   VecPointwiseMult(tempVec2_p_,avec_p_,Uddot_p_);

   //add them -- hold tempVec3
   VecWAXPY(tempVec3_p_,1.,tempVec_p_,tempVec2_p_);
   //VecWAXPY(Udotp1_p_,1.,Udot_p_,tempVec3_p_); 

   //assemble p+1 level vectors
   VecAssemblyBegin(tempVec_p_);
   VecAssemblyEnd(tempVec_p_);
   VecAssemblyBegin(tempVec2_p_);
   VecAssemblyEnd(tempVec2_p_);
   VecAssemblyBegin(tempVec3_p_);
   VecAssemblyEnd(tempVec3_p_);
   VecAssemblyBegin(avec_p_);
   VecAssemblyEnd(avec_p_);
   VecAssemblyBegin(Up1_p_);
   VecAssemblyEnd  (Up1_p_);
   VecAssemblyBegin(Udotp1_p_);
   VecAssemblyEnd  (Udotp1_p_);
   VecAssemblyBegin(Uddotp1_p_);
   VecAssemblyEnd  (Uddotp1_p_);

   //set n+1 values to n
   PetscScalar u1=0.,u2=0.,u3=0.,u1p1=0.,u2p1=0.,u3p1=0.;
   VecCopy(Up1_p_,U_p_);
   VecCopy(Udotp1_p_,Udot_p_);
   VecCopy(Uddotp1_p_,Uddot_p_);

   //assemble p level vectors
   VecAssemblyBegin(U_p_);
   VecAssemblyEnd(U_p_);
   VecAssemblyBegin(Udot_p_);
   VecAssemblyEnd(Udot_p_);
   VecAssemblyBegin(Uddot_p_);
   VecAssemblyEnd(Uddot_p_); 

#endif
 }
 
 void fem::prepareVarsForParallelIO(vector<double>& localFaceForces,
				    vector<double>& localFaceDispls,
				    vector<int>&    nnodesVec)
 {
#if (USEPETSC>0)

   //gather node (XYZ comp only) forces and displacements across partition boundaries and prepare for fast parallel IO dump
   Screen::MasterInfo("Preparing for parallel IO");

   //scatter the external force vector
   VecScatterBegin(scatterToAll_p_,forceVector_p_,myForceVector_p_,INSERT_VALUES,SCATTER_FORWARD);
   VecScatterEnd  (scatterToAll_p_,forceVector_p_,myForceVector_p_,INSERT_VALUES,SCATTER_FORWARD);
   //scatter the converged  solution vector
   VecScatterBegin(scatterToAll_p_,Uk_p_,mySolutionVector_p_,INSERT_VALUES,SCATTER_FORWARD);
   VecScatterEnd  (scatterToAll_p_,Uk_p_,mySolutionVector_p_,INSERT_VALUES,SCATTER_FORWARD);

   for (int e=0;e<myElementCount_;++e)
     {

       int nnodes_loc   = elementData_[e].nnodes;
       int nnodesOffset = nnodesLocalOffset_[e];
       double invnn     = 1./float(nnodes_loc);
       nnodesVec[e]     = nnodes_loc;

       for (int n=0;n<nnodes_loc;++n)
	 {

	   PetscReal fx;PetscReal fy;PetscReal fz;
	   int node       = faces_[n + nnodesOffset];
	   int ndofOffset = ndofLocalOffset_[node];
	   int indexX     = ndofOffset + 0;
	   int indexY     = ndofOffset + 1;
	   int indexZ     = ndofOffset + 2;

	   VecGetValues(myForceVector_p_,1,&indexX,&fx);
	   VecGetValues(myForceVector_p_,1,&indexY,&fy);
	   VecGetValues(myForceVector_p_,1,&indexZ,&fz);

	   localFaceForces[e*xyz+0] += invnn*fx;
	   localFaceForces[e*xyz+1] += invnn*fy;
	   localFaceForces[e*xyz+2] += invnn*fz;

	 }
     }

   for (int e=0;e<myElementCount_;++e)
     {

       int nnodes_loc   = elementData_[e].nnodes;
       int nnodesOffset = nnodesLocalOffset_[e];
       double invnn     = 1./float(nnodes_loc);

       for (int n=0;n<nnodes_loc;++n)
	 {

	   PetscReal dx;PetscReal dy;PetscReal dz;
	   int node       = faces_[n + nnodesOffset];
	   int ndofOffset = ndofLocalOffset_[node];
	   int indexX     = ndofOffset + 0;
	   int indexY     = ndofOffset + 1;
	   int indexZ     = ndofOffset + 2;

	   VecGetValues(mySolutionVector_p_,1,&indexX,&dx);
	   VecGetValues(mySolutionVector_p_,1,&indexY,&dy);
	   VecGetValues(mySolutionVector_p_,1,&indexZ,&dz);

	   localFaceDispls[e*xyz+0] += invnn*dx;
	   localFaceDispls[e*xyz+1] += invnn*dy;
	   localFaceDispls[e*xyz+2] += invnn*dz;

	 }
     }

  }
 
  void fem::ini_uSolve() 
 {

   Screen::MasterInfo("Creating Newmark-Beta Coefficients");
   Tint_phi_  = nm_gamma_;
   Tint_beta_ = nm_beta_;//0.25*pow((Tint_phi_+0.5),2.);
   Tint_a1_   = 1./(Tint_beta_*dt_*dt_);
   Tint_a2_   = Tint_phi_/(Tint_beta_*dt_);
   Tint_a3_   = 1./(Tint_beta_*dt_);
   Tint_a4_   = 1.-Tint_phi_/Tint_beta_;
   Tint_a5_   = 1./(2.*Tint_beta_)-1.;
   Tint_a6_   = (1.-Tint_phi_/(2*Tint_beta_))*dt_;

   Screen::MasterInfo("Phi:  "   + to_string(Tint_phi_ ));
   Screen::MasterInfo("Beta: "   + to_string(Tint_beta_));
   Screen::MasterInfo("a1:   "   + to_string(Tint_a1_  ));
   Screen::MasterInfo("a2:   "   + to_string(Tint_a2_  ));
   Screen::MasterInfo("a3:   "   + to_string(Tint_a3_  ));
   Screen::MasterInfo("a4:   "   + to_string(Tint_a4_  ));
   Screen::MasterInfo("a5:   "   + to_string(Tint_a5_  ));
   Screen::MasterInfo("a6:   "   + to_string(Tint_a6_  ));

#endif
 }
  
 void fem::uSolveCN2(int& nt)
 {
#if (USEPETSC>0)

   double invdt = 1./dt_;
   surfaceTriangulation surfTri;

   Screen::MasterInfo("Timestep:   "+to_string(nt)+"/"+to_string(numTimesteps_+restartTimestep_));
   cout << "==================================================" << endl;
   //applying FSI/special loading
   compLoadVector(nt);
   //output load vector to file
   if (nt % ntSkip_ == 0 && nt != restartTimestep_) dumpFEMForces(nt);

   //create rhs for u_p1_ in "parts"

   Screen::MasterInfo("Creating RHS");
   //"part 1" (hold tempVec)
   //massmat_ * un * a1 + fn+1
   MatMult(Mmat_p_,U_p_, tempVec_p_);
   VecSet(avec_p_,Tint_a1_);
   VecPointwiseMult(tempVec2_p_,avec_p_,tempVec_p_);
   VecWAXPY(tempVec_p_,1.,loadVector_p_,tempVec2_p_);

   //"part 2" (hold tempVec2)
   //Mmat_p_ * udot * a3 + part 1
   MatMult(Mmat_p_,Udot_p_,tempVec2_p_);
   VecSet(avec_p_,Tint_a3_);
   VecPointwiseMult(tempVec3_p_,avec_p_,tempVec2_p_);

   //"part 2 + part 1"
   VecWAXPY(tempVec2_p_,1.,tempVec_p_,tempVec3_p_);

   //"part 3" (rhsU_ output)
   //Mmat_p_ * uddot * a5
   MatMult(Mmat_p_,Uddot_p_,tempVec_p_);
   VecSet(avec_p_,Tint_a5_);
   VecPointwiseMult(tempVec3_p_,avec_p_,tempVec_p_);

   //"part 3" + "part 2" + "part 1"
   VecWAXPY(rhs_p_,1.,tempVec2_p_,tempVec3_p_);

   //finish assembling rhs and sol vectors for timestep n
   VecAssemblyBegin(rhs_p_);
   VecAssemblyEnd  (rhs_p_);
   VecAssemblyBegin(loadVector_p_);
   VecAssemblyEnd(loadVector_p_); 
   VecAssemblyBegin(Up1_p_);
   VecAssemblyEnd(Up1_p_);

   //solve
   Screen::MasterInfo("Solving Linear System");
   KSPSetInitialGuessNonzero (kspOp_p_,PETSC_TRUE);      
   KSPSolve(kspOp_p_,rhs_p_,Up1_p_);

#ifdef fileOutput
   PetscViewer solViewer;
   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"sol.output",&solViewer);
   PetscViewer rhsViewer;
   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"rhs.output",&rhsViewer);
   VecView(Up1_p_,solViewer);
   VecView(rhs_p_,rhsViewer);
   MatMult(mat_p_,Up1_p_,rhs_p_);
   VecView(rhs_p_,solViewer);
   MatView(mat_p_,solViewer);
#endif
   //calculate udot+1 and uddot+1 from u_p1_ value in "parts"
   //for rhsUddot_
   //"part 1" = a1* [un+1 - un]   (hold tempVec2)
   VecWAXPY(tempVec_p_,-1,U_p_,Up1_p_);
   VecSet(avec_p_,Tint_a1_);
   VecPointwiseMult(tempVec2_p_,avec_p_,tempVec_p_);

   //"part 2" = a3 * udot (hold tempVec)
   VecSet(avec_p_,Tint_a3_);
   VecPointwiseMult(tempVec_p_,avec_p_,Udot_p_);

   //"part 1" - "part 2" (hold tempVec3)
   VecWAXPY(tempVec3_p_,-1.,tempVec_p_,tempVec2_p_);

   //"part 3" = a5 * uddot (hold tempVec)
   VecSet(avec_p_,Tint_a5_);
   VecPointwiseMult(tempVec_p_,avec_p_,Uddot_p_);

   //"part 1&2" - "part 3" 
   VecWAXPY(Uddotp1_p_,-1.,tempVec_p_,tempVec3_p_);

   //for rhsUdot_
   //"part 1" = a2*[un+1 - un] (hold tempVec2)
   VecWAXPY(tempVec_p_,-1.,U_p_,Up1_p_);
   VecSet(avec_p_,Tint_a2_);
   VecPointwiseMult(tempVec2_p_,avec_p_,tempVec_p_);

   //"part 2" = a4 * udot (hold tempVec)
   VecSet(avec_p_,Tint_a4_);
   VecPointwiseMult(tempVec_p_,avec_p_,Udot_p_);

   //"part 1" + "part 2" (hold tempVec3)
   VecWAXPY(tempVec3_p_,1.,tempVec_p_,tempVec2_p_);

   //"part 3" = a6 * uddot (hold tempVec)
   VecSet(avec_p_,Tint_a6_);
   VecPointwiseMult(tempVec_p_,avec_p_,Uddot_p_);

   //"part 1&2" + "part 3" 
   VecWAXPY(Udotp1_p_,1.,tempVec_p_,tempVec3_p_);

   //finsih assembling n+1 vectors before writing to file
   VecAssemblyBegin(Udotp1_p_);
   VecAssemblyEnd  (Udotp1_p_);
   VecAssemblyBegin(Uddotp1_p_);
   VecAssemblyEnd  (Uddotp1_p_);

   //set n+1 values to n
   PetscScalar u1,u2,u3,u1p1,u2p1,u3p1;
   //this isnt sized right
   vector<vector<double> > restartVec(nof_nodes_glob_,vector<double>(6+6+6*6+6+9,0.));

   for (int nn=0;nn<nof_nodes_glob_;++nn) 
     {
       int ndof_loc = ndofPerNode_[nn];
       //update node locations
       if (DIM == 3)
         {
           //fill incremental solution vector
           for (int nd=0;nd<ndof_loc;++nd)
             {
               PetscScalar disp;
               int dispIndex = ndof_loc*nn + nd;
               VecGetValues(Up1_p_,1,&dispIndex,&disp);
               deltaUi_[nn*xyz+nd] = disp;
             }
           //update translational DOF
           for (int idir=0;idir<DIM;++idir) 
             {                     
               //store last nodal pos in Uim1_
               Uim1_   [nn*xyz+idir] = nodes_  [nn*xyz+idir];
               nodes_  [nn*xyz+idir] = deltaUi_[nn*xyz+idir] + nodes0_[nn*xyz+idir];
               nodeVel_[nn*xyz+idir] = (nodes_ [nn*xyz+idir] - Uim1_[nn*xyz+idir])*invdt; 
             }
         }
       else
         {
           for (int idir=0;idir<DIM;++idir) 
             {                 
               PetscScalar disp;
               int dispIndex = ndof_loc*nn + idir;
               VecGetValues(Up1_p_,1,&dispIndex,&disp);
               deltaUi_[nn*xyz+idir] = disp;
             }
           //use Uim1_ to store node pos from last timestep
           Uim1_   [nn*xyz+dirVecIndex_] = nodes_  [nn*xyz+dirVecIndex_];
           nodes_  [nn*xyz+dirVecIndex_] = deltaUi_[nn*xyz+0] + nodes0_[nn*xyz+dirVecIndex_];
           nodeVel_[nn*xyz+dirVecIndex_] = (nodes_ [nn*xyz+dirVecIndex_] - Uim1_[nn*xyz+dirVecIndex_])*invdt;
         }

       int index;
       for (int nd=0;nd<ndof_loc;++nd)
         {
           int petscindex = nn*ndof_loc + nd;
           //get u+1 values
           VecGetValues(Up1_p_    ,1,&petscindex,&u1p1);
           VecGetValues(Udotp1_p_ ,1,&petscindex,&u2p1);
           VecGetValues(Uddotp1_p_,1,&petscindex,&u3p1);
           //get u values
           VecGetValues(U_p_    ,1,&petscindex,&u1);
           VecGetValues(Udot_p_ ,1,&petscindex,&u2);
           VecGetValues(Uddot_p_,1,&petscindex,&u3);
           //set n+1 values to n
           VecSetValues(U_p_    ,1,&petscindex,&u1p1,INSERT_VALUES);
           VecSetValues(Udot_p_ ,1,&petscindex,&u2p1,INSERT_VALUES);
           VecSetValues(Uddot_p_,1,&petscindex,&u3p1,INSERT_VALUES);                

           if (nt % restartSkip_ == 0)
             {
               index = 6+nd+ndof_loc;
               restartVec[nn][6+nd         ] = forceVec_[nn*ndof_loc+nd];
               restartVec[nn][index        ] = deltaUi_ [nn*xyz+nd];
               restartVec[nn][index+ndof_loc  ] = u1;
               restartVec[nn][index+ndof_loc*2] = u2;
               restartVec[nn][index+ndof_loc*3] = u3;
               restartVec[nn][index+ndof_loc*4] = u1p1;
               restartVec[nn][index+ndof_loc*5] = u2p1;
               restartVec[nn][index+ndof_loc*6] = u3p1;

             }
         }

       if (nt % restartSkip_ == 0)
         {
           //current and last nodal pos
           restartVec[nn][0] = nodes_  [nn*xyz+0];
           restartVec[nn][1] = nodes_  [nn*xyz+1];
           restartVec[nn][2] = nodes_  [nn*xyz+2];
           restartVec[nn][3] = Uim1_   [nn*xyz+0];
           restartVec[nn][4] = Uim1_   [nn*xyz+1];
           restartVec[nn][5] = Uim1_   [nn*xyz+2];
           //director vectors
           for (int i=0;i<3;++i)
             {
               //restartVec[nn][index+ndof_loc*6+(i+1)  ] = Vnt_[nn][i];
               //restartVec[nn][index+ndof_loc*6+(i+1)+3] = V1t_[nn][i];
               //restartVec[nn][index+ndof_loc*6+(i+1)+6] = V2t_[nn][i];
             }
         }
     }
   //write to restart file
   if (nt % restartSkip_ == 0) writeLinearRestart(restartVec,nt);

   if (nt % ntSkip_ == 0 && nt != restartTimestep_)
     {
       //output updated triangulation
       ofstream surf_out;
       string surface_filename;
       ostringstream OSS;
       OSS << "post/fem/surf_" << setw(8) << setfill('0') << nt << ".vtk";
       surface_filename = OSS.str();
       surf_out.open(surface_filename);
       //surfTri.writeTriangulation_data(surface_filename,nodes_,faces_,DIM,forceVec_,deltaUi_,nodeVel_,Area_,Vnt_,V1t_,V2t_);
       surf_out.close();
     }
   //dumpMat(deltaUi_,"deltaUi_");
   //wait();
   //output node data
   if (recordNodeOutput_ && mypeno_ == 0) outputNodeData(nt);
   updateDirectorVecs();
   //finish assembling u vectors, nullvecs, and avec
   VecAssemblyBegin(U_p_);
   VecAssemblyEnd(U_p_);
   VecAssemblyBegin(Udot_p_);
   VecAssemblyEnd(Udot_p_);
   VecAssemblyBegin(Uddot_p_);
   VecAssemblyEnd(Uddot_p_); 
   VecAssemblyBegin(tempVec_p_);
   VecAssemblyEnd(tempVec_p_);
   VecAssemblyBegin(tempVec2_p_);
   VecAssemblyEnd(tempVec2_p_);
   VecAssemblyBegin(tempVec3_p_);
   VecAssemblyEnd(tempVec3_p_);
   VecAssemblyBegin(avec_p_);
   VecAssemblyEnd(avec_p_);

   if (useAuxGeometry_ || isFSI_) auxDispTransfer();

   cout << "==================================================" << endl;
   cout << endl;

#endif
 }
  
 void fem::linearRungeKutta4(int& nt)
 {
   //4th-order RK scheme for geometrically linear FEM
   double invdt = 1./dt_;
   surfaceTriangulation surfTri;

   compLoadVector(nt);

   Screen::MasterInfo("Timestep:   "+to_string(nt)+"/"+to_string(numTimesteps_));
   cout << "==================================================" << endl;

   Screen::MasterInfo("Creating RK4 RHS");
   vector<double> coefs(4,0.);
   coefs = {0.,0.5,0.5,1.};
   //this isnt sized right
   vector<vector<double> > rhs_u(4,vector<double>(nof_nodes_glob_*6,0.));
   vector<vector<double> > rhs_p(4,vector<double>(nof_nodes_glob_*6,0.));
   //RK4 substeps
   for (int rkss=0;rkss<4;++rkss)
     {
       compRHS_RK4(rkss,coefs[rkss],rhs_u,rhs_p);
     }
   //update u and p = udot
   for (int nn=0;nn<nof_nodes_glob_;++nn) 
     {
       int ndof_loc = ndofPerNode_[nn];
       for (int dof=0;dof<ndof_loc;++dof)
         {
           int index = nn*ndof_loc + dof;
           //store u update in deltaUi_
           deltaUi_[nn*xyz+dof] += 0.16666666666*dt_*( rhs_u[0][index] + 2.*rhs_u[1][index] + 2.*rhs_u[2][index] + rhs_u[3][index] );
           //store p = udot in Uim1_
           Uim1_   [nn*xyz+dof] += 0.16666666666*dt_*( rhs_p[0][index] + 2.*rhs_p[1][index] + 2.*rhs_p[2][index] + rhs_p[3][index] );
         }
     }
   //dumpMat(deltaUi_,"u");
   //dumpMat(Uim1_,"udot");
   //wait();
   //update node locations
   if (DIM == 3)
     {        
       for (int nn=0;nn<nof_nodes_glob_;++nn) 
         {
           for (int idir=0;idir<DIM;++idir) 
             { 
               nodes_[nn*xyz+idir] = deltaUi_[nn*xyz+idir] + nodes0_[nn*xyz+idir];
             }
         }
     }
   else if (DIM == 2)
     {
       for (int nn=0;nn<nof_nodes_glob_;++nn) 
         {            
           nodes_[nn*xyz+1] = deltaUi_[nn*xyz+0] + nodes0_[nn*xyz+1];
         }
     }

   //dumpMat(forceVec_,"forcevec");
   //dumpMat(deltaUi_,"deltaUi");
   //wait("delta");
   //write to restart file

   if (nt % ntSkip_ == 0)
     {
       //output updated triangulation
       surfaceTriangulation surfTri;
       ofstream surf_out;
       string surface_filename;
       ostringstream OSS;
       OSS << "post/fem/surf_" << setw(8) << setfill('0') << nt << ".vtk";
       surface_filename = OSS.str();
       surf_out.open(surface_filename);
       surfTri.writeTriangulation(surface_filename,rootTotalNodes_,rootTotalFaces_,xyz,nof_nodes_glob_,nof_faces_glob_);
       surf_out.close();
     }
   //output node data
   if (recordNodeOutput_ && mypeno_ == 0) outputNodeData(nt);

   cout << "==================================================" << endl;
   cout << endl;

   if (useAuxGeometry_ || isFSI_) auxDispTransfer();

 }
    
 void fem::compRHS_RK4(int& rkss,
                       double& coef,
                       vector<vector<double> >& rhs_u,
                       vector<vector<double> >& rhs_p)
 {
#if (USEPETSC>0)
   int RK0 = 0;
   if (rkss > 0) RK0 = 1;

   //1st rhs: udot = p
   for (int nn=0;nn<nof_nodes_glob_;++nn)
     {
       int ndof_loc = ndofPerNode_[nn];
       for (int dof=0;dof<ndof_loc;++dof)
         {
           int index = nn*ndof_loc + dof;
           rhs_u[rkss][index] = Uim1_   [nn*xyz+dof] + dt_*coef*rhs_p[rkss-1*RK0][index];

           PetscScalar value  = deltaUi_[nn*xyz+dof] + dt_*coef*rhs_u[rkss-1*RK0][index];
           VecSetValues(U_p_,1,&index,&value,INSERT_VALUES);
         }
     }
   //2nd rhs: pdot = M*(f - Ku)
   MatMult (Kmat_p_,U_p_,tempVec_p_);
   VecWAXPY(tempVec2_p_,-1.,tempVec_p_,loadVector_p_);
   MatMult (Mmat_p_,tempVec2_p_,tempVec_p_);

   for (int nn=0;nn<nof_nodes_glob_;++nn)
     {
       int ndof_loc = ndofPerNode_[nn];
       for (int dof=0;dof<ndof_loc;++dof)
         {
           int index = nn*ndof_loc + dof;
           PetscScalar value;
           VecGetValues(tempVec_p_,1,&index,&value);

           rhs_p[rkss][index] = value;
         }
     }

#endif
 }
  
 void fem::outputNodeData(int& nt)
 {
   double timer = time_;
   if (!is_unsteady_) timer = float(nt)/float(numLoadSteps_);
   //track FEM node total displacement
   for (int n=0;n<outputNodes_.size();++n)
     {
       int node = outputNodes_[n];
       int ndof_loc_root = rootNdofPerNode_[node];
       if (nt%outputNodeFreq_[n] == 0)
         {
           if (outputNodeDOF_[n] == "x")
             {
               node_output[n] << nt << "\t" << timer << "\t" << rootTotalNodes_[node*xyz+0]-rootTotalNodes0_[node*xyz+0] << endl;
             }
           else if (outputNodeDOF_[n] == "y")
             {
               node_output[n] << nt << "\t" << timer << "\t" << rootTotalNodes_[node*xyz+1]-rootTotalNodes0_[node*xyz+1] << endl;
             }
           else if (outputNodeDOF_[n] == "z")
             {
               node_output[n] << nt << "\t" << timer << "\t" << rootTotalNodes_[node*xyz+2]-rootTotalNodes0_[node*xyz+2] << endl;
             }
           else if (outputNodeDOF_[n] == "alpha")
             {
             }
           else if (outputNodeDOF_[n] == "beta")
             {
             }
           else if (outputNodeDOF_[n] == "all")
             {
               node_output[n] << nt << "\t" << timer << "\t";
               for (int dir=0;dir<3;++dir)
                 {
                   node_output[n] << rootTotalNodes_[node*xyz+dir]-rootTotalNodes0_[node*xyz+dir] << "\t";
                 }
               node_output[n] << endl;
             }
         }        
     }

 }
 
 void fem::outputCableData(int& nt)
 {
   double timer = time_;
   if (!is_unsteady_) timer = float(nt)/float(numLoadSteps_);

   //track FEM cable internal forces
   for (int c=0;c<nIOCables_;++c)
     {

       int cableLocID   = myIOCables_[c];
       int cableGlobID  = elemMaps_.elemLocal2Global[cableLocID];
       int nnodesOffset = nnodesLocalOffset_[cableLocID];
       int face_node1   = faces_[0 + nnodesOffset];
       int face_node2   = faces_[1 + nnodesOffset];
       int ndofOffset1  = ndofLocalOffset_[face_node1];
       int ndofOffset2  = ndofLocalOffset_[face_node2];

       cable_output[c] << nt << "\t" << time_ << "\t" << myCableTensionVector_[cableLocID] << endl;

     }

 }
 
 void fem::applyPressureLoading()
 {

   double pressureLoading_local = pressureMag_;
   if (time_ < pressureTime_ || !is_unsteady_)
     {
       int elemCounter = myElementCount_;
       if (normalPressure_)
         {
           for(int fg=0;fg<myElementCount_;++fg)
             {
               //hard coded for MITC3 only right now
               if (elementData_[fg].type == MITC3)
                 {
                   int nnodes_loc   = elementData_[fg].nnodes;
                   int nnodesOffset = nnodesLocalOffset_[fg];
                   for (int n=0;n<nnodes_loc;++n)
                     {
                       int node         = faces_            [n + nnodesOffset];
                       int ndofOffset   = ndofLocalOffset_  [node];
                       int vectorNode   = myLocalNodesNormalVectorStatus_[node];
                       for (int idir=0;idir<DIM;++idir)
                         {
                           forceVec_[ndofOffset + idir] += Vnt_[vectorNode*xyz+idir]*pressureLoading_local*Area_[fg]/float(numberOfElementsThatClaimMe_ext_[node]);
                         }
                     }
                 }
             }
         }
       else
         {
           for(int fg=0;fg<myElementCount_;++fg)
             {
               if (elementData_[fg].type == MITC3)
                 {
                   int nnodes_loc = elementData_[fg].nnodes;
                   for (int n=0;n<nnodes_loc;++n)
                     {
                       int nnodesOffset = nnodesLocalOffset_[fg];
                       int node         = faces_            [n + nnodesOffset];
                       int ndofOffset   = ndofLocalOffset_  [node];
                       forceVec_[pressureDirection_ + ndofOffset] += pressureLoading_local*Area_[fg]/float(numberOfElementsThatClaimMe_ext_[node]);
                     }
                 }
             }
         }
     }

 }
  
 void fem::applyBodyForces()
 {
   for(int fg=0;fg<myElementCount_;++fg)
     {
       int nnodes_loc   = elementData_[fg].nnodes;
       int nnodesOffset = nnodesLocalOffset_[fg];
       for (int n=0;n<nnodes_loc;++n)
	 {
	   int node         = faces_            [n + nnodesOffset];
	   int ndofOffset   = ndofLocalOffset_  [node];
	   for (int idir=0;idir<DIM;++idir)
	     {/*                                        rho      *                V                * g    */
	       forceVec_[ndofOffset + idir] += defaultDensity_[fg]*defaultThickness_*Area_[fg]*gravity_[idir]/float(numberOfElementsThatClaimMe_ext_[node]);
	     }
	 }
     }

 }
  
 void fem::readLinearRestart(int& restart_timestep,vector<vector<double> >& restartVecIn)
 {
   //get restart filename from the restart timestep
   string restart_filename;
   ostringstream OSS;
   OSS << "restart_nt_" << setw(8) << setfill('0') << restart_timestep << ".dat";
   restart_filename = OSS.str();

   ifstream file;
   file.open("restart_fem/"+restart_filename);
   if (file.fail()) Screen::MasterError("Bad file in [fem.cc] readLinearRestart");

   for (int nn=0;nn<nof_nodes_glob_;++nn)
     {
       for (int i=0;i<restartVecIn[0].size();++i)
         {
           file >> restartVecIn[nn][i];
         }
     }

   file.close();
 }
  
 void fem::writeLinearRestart(vector<vector<double> >& restartVec, int& timeStep)
 {

   ierr = system("rm restart_fem/*");
   Screen::MasterInfo("Writing linear restart file");
   string restart_filename;
   ostringstream OSS;
   OSS << "restart_nt_" << setw(8) << setfill('0') << timeStep << ".dat";
   restart_filename = OSS.str();
   //write restart timestep to file
   ofstream restart_timestep;
   restart_timestep.open("restart_fem/restart_data.dat");
   if (restart_timestep.fail()) Screen::MasterError("Bad file in [fem.cc] writeLinearRestart | restart_timestep");
   restart_timestep.precision(12);
   restart_timestep.setf(ios::fixed);
   restart_timestep.setf(ios::showpoint);    
   restart_timestep << timeStep;
   //write restart vector to file
   ofstream restart;
   restart.open("restart_fem/"+restart_filename);
   if (restart.fail()) Screen::MasterError("Bad file in [fem.cc] writeLinearRestart | restart");
   restart.precision(12);
   restart.setf(ios::fixed);
   restart.setf(ios::showpoint);    

   restart << setw(20);
   for (int nn=0;nn<nof_nodes_glob_;++nn)
     {
       for (int i=0;i<restartVec[0].size();++i)
         {
           restart << restartVec[nn][i] << setw(20);
         }
       restart << endl;
     }
   restart.close();

 }
  
 void fem::dumpFEMForces(int& nt)
 {
   ofstream load_out;
   string filename_output;
   ostringstream sol;
   sol << "post/forces/force_" << setw(8) << setfill('0') << nt << ".dat";
   filename_output = sol.str();
   load_out.open(filename_output);
   if(load_out.fail()) Screen::MasterError("Could not open dumpFEMForces file");
   load_out.precision(10);

   for (int i=0;i<nof_nodes_glob_;++i)
     {
       int ndof_loc = ndofPerNode_[i];

       for (int idir=0;idir<xyz;++idir)
         {
           load_out << setw(20) << nodes_[i*xyz+idir];
         }

       for (int dof=0;dof<ndof_loc;++dof)
         {            
           load_out << setw(20) << forceVec_[i*ndof_loc+dof];
         }
       load_out << endl;
     }
 }
  
 void fem::intermediateOutput(int& nt,
                              vector<vector<double> >& mat)
 {
   ofstream output;
   string filename_output;
   ostringstream sol;
   sol << "intermediate_output_FEM_" << setw(8) << setfill('0') << nt << ".dat";
   filename_output = sol.str();
   output.open(filename_output);
   if(output.fail()) Screen::MasterError("Could not open dumpFEMForces file");
   output.precision(10);

   for (int r=0;r<mat.size();++r)
     {
       for (int c=0;c<mat[0].size();++c)
         {
           output << setw(20) << mat[r][c];
         }
       output << endl;
     }
 }

 void fem::intermediateOutput(int& nt,
                              vector<vector<double> >& mat1,
                              vector<vector<double> >& mat2)
 {
   ofstream output;
   string filename_output;
   ostringstream sol;
   sol << "intermediate_output_FEM_" << setw(8) << setfill('0') << nt << ".dat";
   filename_output = sol.str();
   output.open(filename_output);
   if(output.fail()) Screen::MasterError("Could not open dumpFEMForces file");
   output.precision(10);

   for (int r=0;r<mat1.size();++r)
     {
       for (int c=0;c<mat1[0].size();++c)
         {
           output << setw(20) << mat1[r][c];
         }

       for (int c=0;c<mat2[0].size();++c)
         {
           output << setw(20) << mat2[r][c];
         }
       output << endl;
     }
 }

 void fem::dumpAUXForces(int& nt)
 {
   if (mypeno_ == 0)
     {
       ofstream load_out;
       string filename_output;
       ostringstream sol;
       sol << "post/forces/aux_force_" << setw(8) << setfill('0') << nt << ".dat";
       filename_output = sol.str();
       load_out.open(filename_output);
       if(load_out.fail()) Screen::MasterError("Could not open dumpAUXForces file");
       load_out.precision(10);

       for (int mn=0;mn<numGeomNodes_;++mn)
         {

           int gn = movingNodes_[mn];

           for (int idir=0;idir<xyz;++idir)
             {
               int vecIndex = gn*xyz + idir;
               load_out << setw(20) << aux_nodes_[vecIndex];
             }

           for (int dof=0;dof<5;++dof)
             {
               int vecIndex = gn*5 + dof;
               load_out << setw(20) << auxGeomForces_[vecIndex];
             }

           load_out << endl;
         }
     }

 }
  
 void fem::dumpPartition(vector<double>& Nodes,
                         vector<int>&    Faces)
 {
   //dump FEM mesh parallel partition
   Screen::MasterInfo("Dumping domain partition");
   ofstream file;
   string filename = "debugFEM/domainPartition.vtk";
   file.open(filename);
   if(file.fail()) Screen::MasterError("Bad file in writeTriangulation function");
   file.precision(10);
   file.setf(ios::fixed);
   file.setf(ios::showpoint);

   int nodeNum = nof_nodes_glob_;
   int faceNum = nof_faces_glob_;
   vector<int> newElementOrder;

   file << "# ";
   file << "vtk ";
   file << "DataFile ";
   file << "Version ";
   file << "3.0";
   file << endl;
   file << "vtk ";
   file << "output";
   file << endl;
   file << "ASCII";
   file << endl;
   file << "DATASET ";
   file << "POLYDATA";
   file << endl;
   file << "POINTS ";
   file << nodeNum;
   file << " float";
   file << endl;

   int DIM = 3;
   //writing nodal locations
   for (int n=0;n<nodeNum;++n)
     {
       for (int idir=0;idir<xyz;++idir)
         {
           int index = idir + xyz*n;
           file << Nodes[index];
           file << " ";
         }
       file << endl;
     }

   if (numberOfElementTypes_ == 1)
     {
       if (rootNnodesPerEl_[0] == 2)
         {
           file << "LINES ";
           file << faceNum;
           file << " ";
           int idata;
           file << faceNum*3;
           file << endl;
           idata = 2;
         }
       else
         {
           file << "POLYGONS ";
           file << faceNum;
           file << " ";
           int idata;
           file << faceNum*4;
           file << endl;
           idata = 3;
         }

       for (int f=0;f<faceNum;++f)
         {
           newElementOrder.push_back(f);
           file << rootNnodesPerEl_[0] << " ";
           for (int idir=0;idir<rootNnodesPerEl_[0];++idir)
             {            
               int index = idir + rootNnodesPerEl_[0]*f;
               file << Faces[index];
               file << " ";
             }
           file << endl;
         }
     }
   //write connectivity one element type at a time
   else
     {

       int numfaces1 = 0;
       int numfaces2 = 0;
       for (int n=0;n<nof_faces_glob_;++n)
         {
           if (rootNnodesPerEl_[n] == 2)
             {
               numfaces1 +=1;
             }
           if (rootNnodesPerEl_[n] == 3)
             {
               numfaces2 +=1;
             }
         }

       file << "LINES ";
       file << numfaces1;
       file << " ";  
       int idata;
       file << numfaces1*3;
       file << endl;
       idata = 2;
       int offset = 0;
       for (int f=0;f<faceNum;++f)
         {
           if (rootNnodesPerEl_[f] == 2)
             {
               newElementOrder.push_back(f);
               int nodesPerElement = 2;
               file << nodesPerElement << " ";
               for (int idir=0;idir<nodesPerElement;++idir)
                 {
                   int index = idir + offset;
                   file << Faces[index];
                   file << " ";
                 }
               file << endl;
             }
           offset += rootNnodesPerEl_[f];
         }

       file << "POLYGONS ";
       file << numfaces2;
       file << " ";  
       idata;
       file << numfaces2*4;
       file << endl;
       idata = 3;
       offset = 0;

       for (int f=0;f<faceNum;++f)
         {
           if (rootNnodesPerEl_[f] == 3)
             {
               newElementOrder.push_back(f);
               int nodesPerElement = 3;
               file << nodesPerElement << " ";
               for (int idir=0;idir<nodesPerElement;++idir)
                 {
                   int index = idir + offset;
                   file << Faces[index];
                   file << " ";
                 }
               file << endl;
             }
           offset += rootNnodesPerEl_[f];
         }
     }

   file << "CELL_DATA ";
   file << faceNum;
   file << endl;
   file << "SCALARS ";
   file << "CPU ";
   file << "int";
   file << endl;
   file << "LOOKUP_TABLE ";
   file << "default";
   file << endl;

   for (int e=0;e<nof_faces_glob_;++e)
     {
       int map = newElementOrder[e];
       file << elemMaps_.globalElement2ProcMap[map];
       file << endl;
     }

   file.close();    

 }
  
 void fem::dumpGeometryPartition(vector<int>& nodePart, vector<int>& facePart)
 {

   if (mypeno_ == 0) cout << "-------------------------" << endl;
   if (mypeno_ == 0)
     {
       string filename = "debugFEM/geometryFacePartition.vtk";
       Screen::MasterInfo("Dumping geometry face partition to VTK file: "+filename);
       ofstream file;
       file.open(filename);
       if(file.fail()) Screen::MasterError("Bad file in writeTriangulation function");
       file.precision(10);
       file.setf(ios::fixed);
       file.setf(ios::showpoint);

       file << "# ";
       file << "vtk ";
       file << "DataFile ";
       file << "Version ";
       file << "3.0";
       file << endl;
       file << "vtk ";
       file << "output";
       file << endl;
       file << "ASCII";
       file << endl;
       file << "DATASET ";
       file << "POLYDATA";
       file << endl;
       file << "POINTS ";
       file << numGeomNodes_;
       file << " float";
       file << endl;
       //writing nodal locations
       for (int n=0;n<numGeomNodes_;++n) 
         {
           for (int idir=0;idir<xyz;++idir) 
             {
               file << aux_nodes0_[n*xyz + idir];
               file << " ";
             }
           file << endl;
         }

      //write connectivity
       file << "POLYGONS " << numGeomFaces_ << " " << numGeomFaces_ * 4 << endl;
      for (int f = 0; f < numGeomFaces_; ++f)
        {
          file << dimLoc_ << " ";
          for (int n = 0; n < dimLoc_; ++n)
            {
              file << aux_faces_[f * dimLoc_ + n] << " ";
            }
          file << endl;
        }

       //write face cpu affiliation
       file << "CELL_DATA ";
       file << numGeomFaces_;
       file << endl;
       file << "SCALARS ";
       file << "CPU_FACES ";
       file << "float";
       file << endl;
       file << "LOOKUP_TABLE ";
       file << "default";
       file << endl;

       for (int n=0;n<numGeomFaces_;++n)
         {
           file << facePart[n] << endl;;
         }

       file.close();

       //dump geomtry node partitioning
       filename = "";
       filename = "debugFEM/geometryNodePartition.vtk";
       Screen::MasterInfo("Dumping geometry node partition to vtk: " + filename);
       file.open(filename);
       if (file.fail()) Screen::MasterError("Bad node file in fem::dumpGeometryPartition");
       file.precision(12);
       file.setf(ios::fixed);
       file.setf(ios::showpoint);

       file << "# ";
       file << "vtk ";
       file << "DataFile ";
       file << "Version ";
       file << "3.0";
       file << endl;
       file << "vtk ";
       file << "output";
       file << endl;
       file << "ASCII";
       file << endl;
       file << "DATASET ";
       file << "POLYDATA";
       file << endl;
       file << "POINTS ";
       file << numGeomNodes_;
       file << " float";
       file << endl;

       //writing nodal locations
       for (int n = 0; n < numGeomNodes_; ++n)
	 {
	   for (int idir = 0; idir < xyz; ++idir)
	     {
	       file << aux_nodes0_[n * xyz + idir] << " ";
	     }
	   file << endl;
	 }

       //write node cpu affiliation
       file << "POINT_DATA ";
       file << numGeomNodes_;
       file << endl;
       file << "SCALARS ";
       file << "CPU_NODES ";
       file << "float";
       file << endl;
       file << "LOOKUP_TABLE ";
       file << "default";
       file << endl;

       for (int n = 0; n < numGeomNodes_; ++n)
	 {
	   file << nodePart[n];
	   file << endl;
	 }

       file.close();      

     }

 }
       
 void fem::writeRestart(int& nt)
 {
#if (USEPETSC>0)

   //write general restart for linear/nonlinear steady/unsteady parallel, multiple elements

#if (DEDICATED_FEM_PROC == 0)
   Screen::MasterInfo("Writing FEM restart -- MPI");

   // --- write nodal values, velocities, and accelerations from all partitions ---

   //get restart filename from the restart timestep
   ofstream restartstream;
   if (mypeno_==0)
     {
       string restart_filename;
       ostringstream OSS;
       OSS << "restart_fem/restart_nt_" << setw(8) << setfill('0') << nt << ".dat";
       restart_filename = OSS.str();
       restartstream.open(restart_filename);
       if (restartstream.fail()) Screen::MasterError("Not able to open restartstream on file "+restart_filename);
     }

   //write total displacement
   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       int dispIndex = idir + rootNdofOffset_[nn];
	       VecGetValues(rootTotalSolution_p_,1,&dispIndex,&disp);
	       restartstream << disp << " ";
	     }
	   restartstream << endl;
	 }
     }
   //write incremental displacement

   //gather incremental displacements to root, overwrite rootTotalSolution_p_
   VecScatter scatterToRoot_loc;
   VecScatterCreateToZero(deltaUk_p_,&scatterToRoot_loc,&rootTotalSolution_p_);
   VecScatterBegin(scatterToRoot_loc,deltaUk_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);
   VecScatterEnd  (scatterToRoot_loc,deltaUk_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);

   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       int dispIndex = idir + rootNdofOffset_[nn];
	       VecGetValues(rootTotalSolution_p_,1,&dispIndex,&disp);
	       restartstream << disp << " ";
	     }
	   restartstream << endl;
	 }
     }

   //write velocity at n

   //gather velocity at n to root, overwrite rootTotalSolution_p_
   VecScatterBegin(scatterToRoot_loc,Udot_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);
   VecScatterEnd  (scatterToRoot_loc,Udot_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);

   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       int dispIndex = idir + rootNdofOffset_[nn];
	       VecGetValues(rootTotalSolution_p_,1,&dispIndex,&disp);
	       restartstream << disp << " ";
	     }
	   restartstream << endl;
	 }
     }

   //write acceleration at n

   //gather acceleration at n to root, overwrite rootTotalSolution_p_
   VecScatterBegin(scatterToRoot_loc,Uddot_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);
   VecScatterEnd  (scatterToRoot_loc,Uddot_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);

   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       int dispIndex = idir + rootNdofOffset_[nn];
	       VecGetValues(rootTotalSolution_p_,1,&dispIndex,&disp);
	       restartstream << disp << " ";
	     }
	   restartstream << endl;
	 }
     }

   //write velocity at n+1

   //gather velocity at n to root, overwrite rootTotalSolution_p_
   VecScatterBegin(scatterToRoot_loc,Udotp1_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);
   VecScatterEnd  (scatterToRoot_loc,Udotp1_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);

   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       int dispIndex = idir + rootNdofOffset_[nn];
	       VecGetValues(rootTotalSolution_p_,1,&dispIndex,&disp);
	       restartstream << disp << " ";
	     }
	   restartstream << endl;
	 }
     }

   //write acceleration at n+1

   //gather acceleration at n to root, overwrite rootTotalSolution_p_
   VecScatterBegin(scatterToRoot_loc,Uddotp1_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);
   VecScatterEnd  (scatterToRoot_loc,Uddotp1_p_,rootTotalSolution_p_,INSERT_VALUES,SCATTER_FORWARD);

   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       int dispIndex = idir + rootNdofOffset_[nn];
	       VecGetValues(rootTotalSolution_p_,1,&dispIndex,&disp);
	       restartstream << disp << " ";
	     }
	   restartstream << endl;
	 }
     }

   //write Vnt, V1, V2

   //create a var and gather all Vnt's on root
   vector<double> rootVnt;
   if (mypeno_==0) rootVnt.resize(nof_nodes_glob_*xyz,-99.);

   //create a local Vnt for every node, not just localNormalVectorNodes
   vector<double> Vnt_all(nof_nodes_glob_*xyz,-999.);
   for (int n=0;n<myNodeCount_;++n)
     {
       int nglob = elemMaps_.nodeLocal2Global[n];
       int normalVectorNode = myLocalNodesNormalVectorStatus_[n];
       if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),n) != myNormalVectorNodes_.end())
	 {
	   Vnt_all[nglob*xyz + 0] = Vnt_[normalVectorNode*xyz + 0];
	   Vnt_all[nglob*xyz + 1] = Vnt_[normalVectorNode*xyz + 1];
	   Vnt_all[nglob*xyz + 2] = Vnt_[normalVectorNode*xyz + 2];
	 }
     }

   //reduce normals onto root
   ierr = MPI_Reduce(&Vnt_all[0],&rootVnt[0],nof_nodes_glob_*xyz,MPI::DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   if (ierr != 0) cout <<  ierr << " " << " help! " << endl;

   if (mypeno_ == 0)
     {
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   for (int idir=0;idir<xyz;++idir)
	     {
	       restartstream << rootVnt[nn*xyz+idir] << " ";
	     }
	   restartstream << endl;
	 }
     }

   fill(Vnt_all.begin(),Vnt_all.end() ,-999.);
   //re-do for Vnt1
   for (int n=0;n<myNodeCount_;++n)
     {
       int nglob = elemMaps_.nodeLocal2Global[n];
       int normalVectorNode = myLocalNodesNormalVectorStatus_[n];
       if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),n) != myNormalVectorNodes_.end())
	 {
	   Vnt_all[nglob*xyz + 0] = V1t_[normalVectorNode*xyz + 0];
	   Vnt_all[nglob*xyz + 1] = V1t_[normalVectorNode*xyz + 1];
	   Vnt_all[nglob*xyz + 2] = V1t_[normalVectorNode*xyz + 2];
	 }
     }

   //reduce normals onto root
   ierr = MPI_Reduce(&Vnt_all[0],&rootVnt[0],nof_nodes_glob_*xyz,MPI::DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);   

   if (mypeno_ == 0)
     {
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   for (int idir=0;idir<xyz;++idir)
	     {
	       restartstream << rootVnt[nn*xyz+idir] << " ";
	     }
	   restartstream << endl;
	 }
     }

   fill(Vnt_all.begin(),Vnt_all.end() ,-999.);
   //re-do for Vnt2
   for (int n=0;n<myNodeCount_;++n)
     {
       int nglob = elemMaps_.nodeLocal2Global[n];
       int normalVectorNode = myLocalNodesNormalVectorStatus_[n];
       if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),n) != myNormalVectorNodes_.end())
	 {
	   Vnt_all[nglob*xyz + 0] = V2t_[normalVectorNode*xyz + 0];
	   Vnt_all[nglob*xyz + 1] = V2t_[normalVectorNode*xyz + 1];
	   Vnt_all[nglob*xyz + 2] = V2t_[normalVectorNode*xyz + 2];
	 }
     }

   //reduce normals onto root
   ierr = MPI_Reduce(&Vnt_all[0],&rootVnt[0],nof_nodes_glob_*xyz,MPI::DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);   
   if (mypeno_ == 0)
     {
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   for (int idir=0;idir<xyz;++idir)
	     {
	       restartstream << rootVnt[nn*xyz+idir] << " ";
	     }
	   restartstream << endl;
	 }
     }

   fill(Vnt_all.begin(),Vnt_all.end() ,-999.);
   //re-do for myDirVecs
   for (int n=0;n<myNodeCount_;++n)
     {
       int nglob = elemMaps_.nodeLocal2Global[n];
       int normalVectorNode = myLocalNodesNormalVectorStatus_[n];
       if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),n) != myNormalVectorNodes_.end())
	 {
	   Vnt_all[nglob*xyz + 0] = myDirVecs_[normalVectorNode*xyz + 0];
	   Vnt_all[nglob*xyz + 1] = myDirVecs_[normalVectorNode*xyz + 1];
	   Vnt_all[nglob*xyz + 2] = myDirVecs_[normalVectorNode*xyz + 2];
	 }
     }

   //reduce normals onto root
   ierr = MPI_Reduce(&Vnt_all[0],&rootVnt[0],nof_nodes_glob_*xyz,MPI::DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);   
   if (mypeno_ == 0)
     {
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   for (int idir=0;idir<xyz;++idir)
	     {
	       restartstream << rootVnt[nn*xyz+idir] << " ";
	     }
	   restartstream << endl;
	 }
     }

   //write time and close
   if (mypeno_==0)
     {
       restartstream << time_ << endl;
       restartstream.close();   
     }

   // --- write root nodal locations ---

   if (mypeno_ == 0)
     {
       //get restart filename from the restart timestep
       string restart_filename_root;
       ostringstream OSS_root;
       OSS_root << "restart_fem/restart_root_nt_" << setw(8) << setfill('0') << nt << ".dat";
       restart_filename_root = OSS_root.str();

       ofstream file;
       file.open(restart_filename_root);
       file.precision(12);
       file.setf(ios::fixed);
       file.setf(ios::showpoint);

       for (int n=0;n<nof_nodes_glob_;++n)
         {
           for (int dir=0;dir<xyz;++dir)
             {
               file << rootTotalNodes_[n*xyz + dir] << " ";
             }
           file << endl;
         }
       file.close();
     }

   ierr = MPI_Barrier(MPI_COMM_WORLD);
   Screen::MasterInfo("Restart Done");

#else
   //write general restart for linear/nonlinear steady/unsteady parallel, multiple elements
   Screen::MasterInfo("Writing FEM restart -- serial");

   // --- write nodal values, velocities, and accelerations ---

   //get restart filename from the restart timestep
   string restart_filename;
   ostringstream OSS;
   OSS << "restart_fem/restart_nt_" << setw(8) << setfill('0') << nt << ".dat";
   restart_filename = OSS.str();
   ofstream restart;
   restart.open(restart_filename);
   restart.precision(10);
   restart.setf(ios::fixed);
   restart.setf(ios::showpoint);

   //loop over nodes --- write nodal positions
   vector<double> pos(3,0.);
   for (int e=0;e<myElementCount_;++e)
     {

       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {

           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];

           //get nodal positions in smaller form
           pos[0] = nodes_[locN*xyz + 0];
           pos[1] = nodes_[locN*xyz + 1];
           pos[2] = nodes_[locN*xyz + 2];
           for (int idir=0;idir<xyz;++idir)
             {
               restart << pos[idir] << " ";
             }
           restart << endl;
         }
     }
   //loop over nodes --- write nodal velocity and acceleration from PETSc
   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           vector<double> vel_n(ndof_loc,0.);

           //get nodal velocity and acceleration in smaller form
           for (int ndof=0;ndof<ndof_loc;++ndof)
             {

               PetscScalar vel_n_petsc;

               int dispIndex = offset + ndof;
               VecGetValues(Udot_p_,1,&dispIndex,&vel_n_petsc);

               restart << vel_n_petsc << " ";

             }
           restart << endl;
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           vector<double> vel_n(ndof_loc,0.);

           //get nodal velocity and acceleration in smaller form
           for (int ndof=0;ndof<ndof_loc;++ndof)
             {

               PetscScalar acc_n_petsc;

               int dispIndex = offset + ndof;
               VecGetValues(Uddot_p_,1,&dispIndex,&acc_n_petsc);

               restart << acc_n_petsc << " ";

             }
           restart << endl;
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //get nodal velocity and acceleration in smaller form
           for (int ndof=0;ndof<ndof_loc;++ndof)
             {

               PetscScalar vel_np1_petsc;

               int dispIndex = offset + ndof;
               VecGetValues(Udotp1_p_,1,&dispIndex,&vel_np1_petsc);

               restart << vel_np1_petsc << " ";

             }
           restart << endl;
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           vector<double> vel_n(ndof_loc,0.);

           //get nodal velocity and acceleration in smaller form
           for (int ndof=0;ndof<ndof_loc;++ndof)
             {

               PetscScalar acc_np1_petsc;

               int dispIndex = offset + ndof;
               VecGetValues(Uddotp1_p_,1,&dispIndex,&acc_np1_petsc);

               restart << acc_np1_petsc << " ";

             }
           restart << endl;
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //get nodal velocity and acceleration in smaller form
           for (int ndof=0;ndof<ndof_loc;++ndof)
             {
               int dispIndex = offset + ndof;
               deltaUi_[dispIndex];
               restart << deltaUi_[dispIndex] << " ";
             }

           restart << endl;

         }
     }

   //write total nodal displacements
   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           for (int ndof=0;ndof<ndof_loc;++ndof)
             {
               int dispIndex = offset + ndof;
               Ui_[dispIndex];
               restart << Ui_[dispIndex] << " ";
             }
           restart << endl;

         }
     }

   vector<double> vec(xyz,0.);

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //if this node is a normal vector node, write it, otherwise, write -99
           if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),locN) != myNormalVectorNodes_.end())
             {
               int vectorNode = myLocalNodesNormalVectorStatus_[locN];
               //get normal in smaller form
               for (int idir=0;idir<xyz;++idir)
                 {
                   vec[idir] = Vnt_[vectorNode*xyz+idir];
                 }
             }
           else
             {
               //get normal in smaller form
               for (int idir=0;idir<xyz;++idir)
                 {                
                   vec[idir] = -99.;
                 }
             }
           for (int idir=0;idir<xyz;++idir) 
             {
               restart << vec[idir] << " ";
             }
           restart << endl;
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //if this node is a normal vector node, write it, otherwise, write -99
           if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),locN) != myNormalVectorNodes_.end())
             {
               int vectorNode = myLocalNodesNormalVectorStatus_[locN];
               //get normal in smaller form
               for (int idir=0;idir<xyz;++idir)
                 {
                   vec[idir] = V1t_[vectorNode*xyz+idir];
                 }
             }
           else
             {
               //get normal in smaller form
               for (int idir=0;idir<xyz;++idir)
                 {                
                   vec[idir] = -99.;
                 }
             }
           for (int idir=0;idir<xyz;++idir) 
             {
               restart << vec[idir] << " ";
             }
           restart << endl;
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //if this node is a normal vector node, write it, otherwise, write -99
           if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),locN) != myNormalVectorNodes_.end())
             {
               int vectorNode = myLocalNodesNormalVectorStatus_[locN];
               //get normal in smaller form
               for (int idir=0;idir<xyz;++idir)
                 {
                   vec[idir] = V2t_[vectorNode*xyz+idir];
                 }
             }
           else
             {
               //get normal in smaller form
               for (int idir=0;idir<xyz;++idir)
                 {                
                   vec[idir] = -99.;
                 }
             }
           for (int idir=0;idir<xyz;++idir) 
             {
               restart << vec[idir] << " ";
             }
           restart << endl;
         }
     }

   restart << time_ << endl;
   restart.close();

   // --- write root nodal locations ---

   //get restart filename from the restart timestep
   string restart_filename_root;
   ostringstream OSS_root;
   OSS_root << "restart_fem/restart_root_nt_" << setw(8) << setfill('0') << nt << ".dat";
   restart_filename_root = OSS_root.str();

   ofstream file;
   file.open(restart_filename_root);
   file.precision(12);
   file.setf(ios::fixed);
   file.setf(ios::showpoint);

   for (int n=0;n<nof_nodes_glob_;++n)
	 {
	   for (int dir=0;dir<xyz;++dir)
	     {
	       file << rootTotalNodes_[n*xyz + dir] << " ";
	     }
	   file << endl;
	 }
   file.close();

#endif

#endif
 }
  
 void fem::readRestart()
 {
#if (USEPETSC>0)

#define WRITE_RESTART 1

   //read general restart for linear/nonlinear steady/unsteady parallel, multiple elements

#if (WRITE_RESTART == 1)

   ofstream restart_check;   
   if (mypeno_==0)
     {
       //get restart filename from the restart timestep
       string restart_filename2;
       ostringstream OSS2;
       OSS2 << "restart_fem/restart_nt_" << setw(8) << setfill('0') << restartTimestep_ << "_check.dat";
       restart_filename2 = OSS2.str();
       restart_check.open(restart_filename2);
     }

#endif

#if (DEDICATED_FEM_PROC == 0)
   Screen::MasterInfo("Reading FEM restart --- MPI");

   // --- read nodal values, velocities, and accelerations for all partitions ---

   //get restart filename from the restart timestep
   ifstream restartstream;
   if (mypeno_==0)
     {
       string restart_filename;
       ostringstream OSS;
       OSS << "restart_fem/restart_nt_" << setw(8) << setfill('0') << restartTimestep_ << ".dat";
       restart_filename = OSS.str();
       restartstream.open(restart_filename);
       if (restartstream.fail()) Screen::MasterError("Not able to open restartstream on file "+restart_filename);
     }

   vector<double> bigData(nof_nodes_glob_*5,0.);

   //read total displacement
   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       restartstream >> disp;	       
	       int dispIndex = idir + rootNdofOffset_[nn];
	       VecSetValues(rootTotalSolution_p_,1,&dispIndex,&disp,INSERT_VALUES);
	       bigData[dispIndex] = disp;
#if (WRITE_RESTART == 1)
	       if (mypeno_==0) restart_check << disp << " ";
#endif
	     }
#if (WRITE_RESTART == 1)
	   if (mypeno_==0) restart_check << endl;
#endif
	 }
     }

   //broadcast total displacements 
   ierr = MPI_Bcast(&bigData[0],nof_nodes_glob_*5,MPI::DOUBLE,0,MPI_COMM_WORLD);

   //fill U_p_ glob
   PetscScalar disp;
   for (int nn=0;nn<myNodeCount_;++nn)
     {

       int ndof_loc = ndofPerNode_[nn];
       int nGlob = elemMaps_.nodeLocal2Global[nn];
       int local_offset  = ndofLocalOffset_[nn];
       int global_offset = ndofGlobalOffset_[nn];

       for (int idir=0;idir<ndof_loc;++idir)
	 {
	   int localIndex  = idir + local_offset;
	   int globalIndex = idir + global_offset;
	   disp = bigData[globalIndex];
	   VecSetValues(U_p_,1,&globalIndex,&disp,INSERT_VALUES);
	   VecSetValues(Uk_p_,1,&globalIndex,&disp,INSERT_VALUES);
	   Ui_[localIndex] = disp;
	 }
     }

   fill(bigData.begin(),bigData.end(),0.);
   //read incremental displacement

   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       restartstream >> disp;	       
	       int dispIndex = idir + rootNdofOffset_[nn];
	       bigData[dispIndex] = disp;
#if (WRITE_RESTART == 1)
	       if (mypeno_==0) restart_check << disp << " ";
#endif
	     }
#if (WRITE_RESTART == 1)
           if (mypeno_==0) restart_check << endl;
#endif
	 }
     }

   //broadcast incremental displacements 
   ierr = MPI_Bcast(&bigData[0],nof_nodes_glob_*5,MPI::DOUBLE,0,MPI_COMM_WORLD);

   for (int nn=0;nn<myNodeCount_;++nn)
     {

       int ndof_loc = ndofPerNode_[nn];
       int nGlob = elemMaps_.nodeLocal2Global[nn];
       int local_offset  = ndofLocalOffset_[nn];
       int global_offset = ndofGlobalOffset_[nn];

       for (int idir=0;idir<ndof_loc;++idir)
	 {
	   int localIndex  = idir + local_offset;
	   int globalIndex = idir + global_offset;

	   disp = bigData[globalIndex];
	   VecSetValues(deltaUk_p_,1,&globalIndex,&disp,INSERT_VALUES);
	   deltaUi_[localIndex] = disp;
	   Uim1_[localIndex] = Ui_[localIndex] - deltaUi_[localIndex];
	 }
     }

   fill(bigData.begin(),bigData.end(),0.);
   //read velocity at n

   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       restartstream >> disp;	       
	       int dispIndex = idir + rootNdofOffset_[nn];
	       bigData[dispIndex] = disp;
#if (WRITE_RESTART == 1)
	       if (mypeno_==0) restart_check << disp << " ";
#endif
	     }
#if (WRITE_RESTART == 1)
	     if (mypeno_==0) restart_check << endl;
#endif
	 }
     }

   //broadcast velocity
   ierr = MPI_Bcast(&bigData[0],nof_nodes_glob_*5,MPI::DOUBLE,0,MPI_COMM_WORLD);

   for (int nn=0;nn<myNodeCount_;++nn)
     {

       int ndof_loc = ndofPerNode_[nn];
       int nGlob = elemMaps_.nodeLocal2Global[nn];
       int local_offset  = ndofLocalOffset_[nn];
       int global_offset = ndofGlobalOffset_[nn];

       for (int idir=0;idir<ndof_loc;++idir)
	 {
	   int localIndex  = idir + local_offset;
	   int globalIndex = idir + global_offset;

	   disp = bigData[globalIndex];
	   VecSetValues(Udot_p_,1,&globalIndex,&disp,INSERT_VALUES);	  
	 }
     }

   fill(bigData.begin(),bigData.end(),0.);
   //read acceleration at n

   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       restartstream >> disp;	       
	       int dispIndex = idir + rootNdofOffset_[nn];
	       bigData[dispIndex] = disp;
#if (WRITE_RESTART == 1)
	       if (mypeno_==0) restart_check << disp << " ";
#endif
	     }
#if (WRITE_RESTART == 1)
	     if (mypeno_==0) restart_check << endl;
#endif
	 }
     }

   //broadcast acceleration
   ierr = MPI_Bcast(&bigData[0],nof_nodes_glob_*5,MPI::DOUBLE,0,MPI_COMM_WORLD);

   for (int nn=0;nn<myNodeCount_;++nn)
     {

       int ndof_loc = ndofPerNode_[nn];
       int nGlob = elemMaps_.nodeLocal2Global[nn];
       int local_offset  = ndofLocalOffset_[nn];
       int global_offset = ndofGlobalOffset_[nn];

       for (int idir=0;idir<ndof_loc;++idir)
	 {
	   int localIndex  = idir + local_offset;
	   int globalIndex = idir + global_offset;

	   disp = bigData[globalIndex];	   
	   VecSetValues(Uddot_p_,1,&globalIndex,&disp,INSERT_VALUES);
	 }
     }

   fill(bigData.begin(),bigData.end(),0.);
   //read velocity at n+1

   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       restartstream >> disp;	       
	       int dispIndex = idir + rootNdofOffset_[nn];
	       bigData[dispIndex] = disp;
#if (WRITE_RESTART == 1)
	       if (mypeno_==0) restart_check << disp << " ";
#endif
	     }
#if (WRITE_RESTART == 1)
	   if (mypeno_==0) restart_check << endl;
#endif
	 }
     }

   //broadcast velocity
   ierr = MPI_Bcast(&bigData[0],nof_nodes_glob_*5,MPI::DOUBLE,0,MPI_COMM_WORLD);

   for (int nn=0;nn<myNodeCount_;++nn)
     {

       int ndof_loc = ndofPerNode_[nn];
       int nGlob = elemMaps_.nodeLocal2Global[nn];
       int local_offset  = ndofLocalOffset_[nn];
       int global_offset = ndofGlobalOffset_[nn];

       for (int idir=0;idir<ndof_loc;++idir)
	 {
	   int localIndex  = idir + local_offset;
	   int globalIndex = idir + global_offset;

	   disp = bigData[globalIndex];
	   VecSetValues(Udotp1_p_,1,&globalIndex,&disp,INSERT_VALUES);	  
	 }
     }

   fill(bigData.begin(),bigData.end(),0.);
   //read acceleration at n+1

   if (mypeno_ == 0)
     {
       PetscScalar disp;
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   int ndof_loc_root = rootNdofPerNode_[nn];
	   for (int idir=0;idir<ndof_loc_root;++idir)
	     {
	       restartstream >> disp;	       
	       int dispIndex = idir + rootNdofOffset_[nn];
	       bigData[dispIndex] = disp;
#if (WRITE_RESTART == 1)
	       if (mypeno_==0) restart_check << disp << " ";
#endif
	     }
#if (WRITE_RESTART == 1)
	     if (mypeno_==0) restart_check << endl;
#endif
	 }
     }

   //broadcast acceleration
   ierr = MPI_Bcast(&bigData[0],nof_nodes_glob_*5,MPI::DOUBLE,0,MPI_COMM_WORLD);

   for (int nn=0;nn<myNodeCount_;++nn)
     {

       int ndof_loc = ndofPerNode_[nn];
       int nGlob = elemMaps_.nodeLocal2Global[nn];
       int local_offset  = ndofLocalOffset_[nn];
       int global_offset = ndofGlobalOffset_[nn];

       for (int idir=0;idir<ndof_loc;++idir)
	 {
	   int localIndex  = idir + local_offset;
	   int globalIndex = idir + global_offset;

	   disp = bigData[globalIndex];	   
	   VecSetValues(Uddotp1_p_,1,&globalIndex,&disp,INSERT_VALUES);
	 }
     }

   VecAssemblyBegin(U_p_);
   VecAssemblyEnd(U_p_);
   VecAssemblyBegin(Uk_p_);
   VecAssemblyEnd(Uk_p_);
   VecAssemblyBegin(Udot_p_);
   VecAssemblyEnd(Udot_p_);
   VecAssemblyBegin(Uddot_p_);
   VecAssemblyEnd(Uddot_p_);
   VecAssemblyBegin(Udotp1_p_);
   VecAssemblyEnd(Udotp1_p_);
   VecAssemblyBegin(Uddotp1_p_);
   VecAssemblyEnd(Uddotp1_p_);

   //read Vnt

   vector<double> VecStore(nof_nodes_glob_*xyz,-99.);
   if (mypeno_==0)
     {
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   for (int idir=0;idir<xyz;++idir)
	     {	       
	       restartstream >> VecStore[nn*xyz+idir];
#if (WRITE_RESTART == 1)
	       if (mypeno_==0) restart_check << VecStore[nn*xyz+idir] << " ";
#endif
	     }
#if (WRITE_RESTART == 1)
	     if (mypeno_==0) restart_check << endl;
#endif
	 }
     }

   ierr = MPI_Bcast(&VecStore[0],nof_nodes_glob_*xyz,MPI::DOUBLE,0,MPI_COMM_WORLD);

   for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
     {
       int localNodeID = myNormalVectorNodes_[n];
       int nglob = elemMaps_.nodeLocal2Global[localNodeID];

       Vnt_[n*xyz + 0] = VecStore[nglob*xyz + 0];
       Vnt_[n*xyz + 1] = VecStore[nglob*xyz + 1];
       Vnt_[n*xyz + 2] = VecStore[nglob*xyz + 2];
     } 

   //read V1t

   if (mypeno_==0)
     {
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   for (int idir=0;idir<xyz;++idir)
	     {	       
	       restartstream >> VecStore[nn*xyz+idir];
#if (WRITE_RESTART == 1)
               if (mypeno_==0) restart_check << VecStore[nn*xyz+idir] << " ";
#endif
	     }
#if (WRITE_RESTART == 1)
	   if (mypeno_==0) restart_check << endl;
#endif
	 }
     }

   ierr = MPI_Bcast(&VecStore[0],nof_nodes_glob_*xyz,MPI::DOUBLE,0,MPI_COMM_WORLD);

   for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
     {
       int localNodeID = myNormalVectorNodes_[n];
       int nglob = elemMaps_.nodeLocal2Global[localNodeID];
       V1t_[n*xyz + 0] = VecStore[nglob*xyz + 0];
       V1t_[n*xyz + 1] = VecStore[nglob*xyz + 1];
       V1t_[n*xyz + 2] = VecStore[nglob*xyz + 2];     
     }

   //read V2t

   if (mypeno_==0)
     {
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   for (int idir=0;idir<xyz;++idir)
	     {	       
	       restartstream >> VecStore[nn*xyz+idir];
#if (WRITE_RESTART == 1)
               if (mypeno_==0) restart_check << VecStore[nn*xyz+idir] << " ";
#endif
	     }
#if (WRITE_RESTART == 1)
	   if (mypeno_==0) restart_check << endl;
#endif
	 }
     }

   ierr = MPI_Bcast(&VecStore[0],nof_nodes_glob_*xyz,MPI::DOUBLE,0,MPI_COMM_WORLD);

   for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
     {
       int localNodeID = myNormalVectorNodes_[n];
       int nglob = elemMaps_.nodeLocal2Global[localNodeID];
       V2t_[n*xyz + 0] = VecStore[nglob*xyz + 0];
       V2t_[n*xyz + 1] = VecStore[nglob*xyz + 1];
       V2t_[n*xyz + 2] = VecStore[nglob*xyz + 2];     
     }

   //read myDirVecs

   if (mypeno_==0)
     {
       for (int nn=0;nn<nof_nodes_glob_;++nn)
	 {
	   for (int idir=0;idir<xyz;++idir)
	     {	       
	       restartstream >> VecStore[nn*xyz+idir];
#if (WRITE_RESTART == 1)
               if (mypeno_==0) restart_check << VecStore[nn*xyz+idir] << " ";
#endif
	     }
#if (WRITE_RESTART == 1)
	   if (mypeno_==0) restart_check << endl;
#endif
	 }
     }

   ierr = MPI_Bcast(&VecStore[0],nof_nodes_glob_*xyz,MPI::DOUBLE,0,MPI_COMM_WORLD);

   for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
     {
       int localNodeID = myNormalVectorNodes_[n];
       int nglob = elemMaps_.nodeLocal2Global[localNodeID];
       myDirVecs_[n*xyz + 0] = VecStore[nglob*xyz + 0];
       myDirVecs_[n*xyz + 1] = VecStore[nglob*xyz + 1];
       myDirVecs_[n*xyz + 2] = VecStore[nglob*xyz + 2];     
     }

   int alphaIndex = 3;
   int betaIndex  = 4;
   //initialize delta director vector
   for (int n=0;n<myNumberOfNormalVectorNodes_;++n)
     {
       int node       = myNormalVectorNodes_[n];
       int ndofOffset = ndofLocalOffset_[node];
       double alpha   = deltaUi_[alphaIndex + ndofOffset];
       double beta    = deltaUi_[betaIndex  + ndofOffset];

       for (int idir=0;idir<xyz;++idir) 
         { 
           deltaVnt_[n*xyz+idir] = -alpha*V2t_[n*xyz+idir] + beta*V1t_[n*xyz+idir] - 0.5*(alpha*alpha+beta*beta)*Vnt_[n*xyz+idir];
         }
     }

   //read time

   if (mypeno_==0) restartstream >> time_;   

   ierr = MPI_Bcast(&time_,1,MPI::DOUBLE,0,MPI_COMM_WORLD);

#if (WRITE_RESTART == 1)
   if (mypeno_==0) restart_check << time_ << endl;
   if (mypeno_==0) restart_check.close();
#endif

#else

#if (WRITE_RESTART == 1)

   //get restart filename from the restart timestep
   string restart_filename2;
   ostringstream OSS2;
   OSS2 << "restart_fem/restart_nt_" << setw(8) << setfill('0') << restartTimestep_ << "_check.dat";
   restart_filename2 = OSS2.str();
   //ofstream restart_check;
   restart_check.open(restart_filename2);
   restart_check.precision(10);
   restart_check.setf(ios::fixed);
   restart_check.setf(ios::showpoint);

#endif

   // --- read nodal values, velocities, and accelerations ---

   //get restart filename from the restart timestep
   string restart_filename;
   ostringstream OSS;
   OSS << "restart_fem/restart_nt_" << setw(8) << setfill('0') << restartTimestep_ << ".dat";
   restart_filename = OSS.str();
   ifstream restart;
   restart.open(restart_filename);
   restart.precision(10);
   restart.setf(ios::fixed);
   restart.setf(ios::showpoint);

   //loop over nodes --- write nodal positions
   vector<double> pos(3,0.);
   for (int e=0;e<myElementCount_;++e)
     {

       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {

           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];

           for (int idir=0;idir<xyz;++idir) restart >> pos[idir];
           //get nodal positions in smaller form
           nodes_[locN*xyz + 0] = pos[0];
           nodes_[locN*xyz + 1] = pos[1];
           nodes_[locN*xyz + 2] = pos[2];

#if (WRITE_RESTART == 1)
           for (int idir=0;idir<xyz;++idir)
             {
               restart_check << pos[idir] << " ";
             }
           restart_check << endl;
#endif
         }
     }

   //loop over nodes --- write nodal velocity and acceleration from PETSc
   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           vector<double> vel_n(ndof_loc,0.);

           //get nodal velocity and acceleration in smaller form
           for (int ndof=0;ndof<ndof_loc;++ndof)
             {

               PetscScalar vel_n_petsc;
               restart >> vel_n_petsc;

               int dispIndex = offset + ndof;
               VecSetValues(Udot_p_,1,&dispIndex,&vel_n_petsc,INSERT_VALUES);

#if (WRITE_RESTART == 1)
               VecGetValues(Udot_p_,1,&dispIndex,&vel_n_petsc);
               restart_check << vel_n_petsc << " ";
#endif
             }
#if (WRITE_RESTART == 1)
           restart_check << endl;
#endif
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //get nodal velocity and acceleration in smaller form
           for (int ndof=0;ndof<ndof_loc;++ndof)
             {

               PetscScalar acc_n_petsc;
               restart >> acc_n_petsc;

               int dispIndex = offset + ndof;
               VecSetValues(Uddot_p_,1,&dispIndex,&acc_n_petsc,INSERT_VALUES);

#if (WRITE_RESTART == 1)
               VecGetValues(Uddot_p_,1,&dispIndex,&acc_n_petsc);
               restart_check << acc_n_petsc << " ";
#endif
             }
#if (WRITE_RESTART == 1)
           restart_check << endl;
#endif
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //get nodal velocity and acceleration in smaller form
           for (int ndof=0;ndof<ndof_loc;++ndof)
             {

               PetscScalar vel_np1_petsc;
               restart >> vel_np1_petsc;

               int dispIndex = offset + ndof;
               VecSetValues(Udotp1_p_,1,&dispIndex,&vel_np1_petsc,INSERT_VALUES);

#if (WRITE_RESTART == 1)
               VecGetValues(Udotp1_p_,1,&dispIndex,&vel_np1_petsc);
               restart_check << vel_np1_petsc << " ";
#endif
             }
#if (WRITE_RESTART == 1)
           restart_check << endl;
#endif
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //get nodal velocity and acceleration in smaller form
           for (int ndof=0;ndof<ndof_loc;++ndof)
             {

               PetscScalar acc_np1_petsc;
               restart >> acc_np1_petsc;

               int dispIndex = offset + ndof;
               VecSetValues(Uddotp1_p_,1,&dispIndex,&acc_np1_petsc,INSERT_VALUES);

#if (WRITE_RESTART == 1)
               VecGetValues(Uddotp1_p_,1,&dispIndex,&acc_np1_petsc);
               restart_check << acc_np1_petsc << " ";
#endif               
             }
#if (WRITE_RESTART == 1)
           restart_check << endl;
#endif
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //get nodal velocity and acceleration in smaller form
           for (int ndof=0;ndof<ndof_loc;++ndof)
             {
               int dispIndex = offset + ndof;
               restart >> deltaUi_[dispIndex];
#if (WRITE_RESTART == 1)
               restart_check << deltaUi_[dispIndex] << " ";               
#endif
             }
#if (WRITE_RESTART == 1)
           restart_check << endl;
#endif
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //get nodal velocity and acceleration in smaller form
           for (int ndof=0;ndof<ndof_loc;++ndof)
             {
               int dispIndex = offset + ndof;
               restart >> Ui_[dispIndex];
               PetscScalar ui = Ui_[dispIndex];
               VecSetValues(Uk_p_,1,&dispIndex,&ui,INSERT_VALUES);

               VecSetValues(U_p_ ,1,&dispIndex,&ui,INSERT_VALUES);
               VecSetValues(Uk_p_,1,&dispIndex,&ui,INSERT_VALUES);

#if (WRITE_RESTART == 1)
               restart_check << Ui_[dispIndex] << " ";
#endif
             }
#if (WRITE_RESTART == 1)
           restart_check << endl;
#endif
         }
     }

   vector<double> vec(xyz,0.);
   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //if this node is a normal vector node, write it, otherwise, write -99
           if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),locN) != myNormalVectorNodes_.end())
             {
               int vectorNode = myLocalNodesNormalVectorStatus_[locN];
               for (int idir=0;idir<xyz;++idir)
                 {
                   restart >> Vnt_[vectorNode*xyz+idir];
#if (WRITE_RESTART == 1)
                   vec[idir] = Vnt_[vectorNode*xyz+idir];
#endif                   
                 }
             }
#if (WRITE_RESTART == 1)
           else
             {
               for (int idir=0;idir<xyz;++idir)
                 {                
                   vec[idir] = -99.;
                 }
             }

           for (int idir=0;idir<xyz;++idir) 
             {
               restart_check << vec[idir] << " ";
             }
           restart_check << endl;
#endif
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //if this node is a normal vector node, write it, otherwise, write -99
           if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),locN) != myNormalVectorNodes_.end())
             {
               int vectorNode = myLocalNodesNormalVectorStatus_[locN];
               for (int idir=0;idir<xyz;++idir)
                 {
                   restart >> V1t_[vectorNode*xyz+idir];
#if (WRITE_RESTART == 1)
                   vec[idir] = V1t_[vectorNode*xyz+idir];
#endif                   
                 }
             }
#if (WRITE_RESTART == 1)
           else
             {
               for (int idir=0;idir<xyz;++idir)
                 {                
                   vec[idir] = -99.;
                 }
             }

           for (int idir=0;idir<xyz;++idir) 
             {
               restart_check << vec[idir] << " ";
             }
           restart_check << endl;
#endif
         }
     }

   for (int e=0;e<myElementCount_;++e)
     {
       int globE         = elemMaps_.elemLocal2Global[e];
       int nnodes_loc    = elementData_[e].nnodes;
       int nnodes_offset = nnodesLocalOffset_[e];

       for (int n=0;n<nnodes_loc;++n)
         {
           int locN     = faces_[n + nnodes_offset];
           int globN    = elemMaps_.nodeLocal2Global[locN];
           int ndof_loc = ndofPerNode_[locN];
           int offset   = ndofLocalOffset_[locN];

           //if this node is a normal vector node, write it, otherwise, write -99
           if (find(myNormalVectorNodes_.begin(),myNormalVectorNodes_.end(),locN) != myNormalVectorNodes_.end())
             {
               int vectorNode = myLocalNodesNormalVectorStatus_[locN];
               for (int idir=0;idir<xyz;++idir)
                 {
                   restart >> V2t_[vectorNode*xyz+idir];
#if (WRITE_RESTART == 1)
                   vec[idir] = V2t_[vectorNode*xyz+idir];
#endif                   
                 }
             }
#if (WRITE_RESTART == 1)
           else
             {
               for (int idir=0;idir<xyz;++idir)
                 {                
                   vec[idir] = -99.;
                 }
             }

           for (int idir=0;idir<xyz;++idir) 
             {
               restart_check << vec[idir] << " ";
             }
           restart_check << endl;
#endif
         }
     }

   restart >> time_;
#if (WRITE_RESTART == 1)
   restart_check << time_;
#endif
   restart.close();

#endif

#endif
 }
  
 void fem::Step_structural_solver(int& nt)
 {
   //timer
   nt_ = nt;
   if (!isFSI_) time_ = nt_*dt_;

   //unsteady sims
   if (is_unsteady_)
     {
       //nonlinear, unsteady 
       if (is_nonlinear_) 
         {
           if (timeIntegrationMethod_ == "implicit")
             {
	       BetaNewmarkGeneral(nt);
             }
         }
       //linear, unsteady
       else 
         {
           if (timeIntegrationMethod_ == "implicit")
             {
               uSolveCN2(nt);
             }
           else if (timeIntegrationMethod_ == "explicit")
             {
               linearRungeKutta4(nt);
             }
         }
       firstTimeStep_ = false;
     }
   //steady sims
   else 
     {
       //nonlinear, steady
       if (is_nonlinear_) 
         {
           nonLinearSolve();
         } 
       //linear, steady
       else 
         {
           assembleGlobKMat();
           setForceSteady();
           solve(pctype_);
         }
     }

 }
  void fem::Finalize_fem()
  {
#if (USEPETSC>0)
#if (DEDICATED_FEM_PROC == 0)
    VecScatterDestroy(&scatterToAll_p_);
    VecScatterDestroy(&scatterToRoot_p_);
#endif
    PetscFinalize();
#endif
  }

}//end namespace KARMA
