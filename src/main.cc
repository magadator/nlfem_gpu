/*
 * nlfem - Nonlinear FEM Solver
 * Author  : Gokul G. Anugrah
 * Contact : gokulanugrah@gmail.com
 *
 * main.cc — Main driver.
 *   Parses input.dat, populates the linearTriElement solver object,
 *   calls Setup_fem(), runs the time/load-step loop, and finalizes.
 *
 *   Requires PETSc (set PETSC_DIR in config.mk) and ParMETIS.
 *
 * Usage:
 *   mpirun -np <N> ./nlfem input.dat
 *   ./nlfem input.dat                  (single process)
 */

#include <mpi.h>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "options_parser.h"
#include "PreprocessorContext.h"
#include "linearTriElement.h"
#include "screen.h"

using namespace std;
using namespace KARMA;

// Explicitly invoke conversion operators to avoid ambiguous C-style casts
static vector<int>    toIntVec(OptionEntry e) { return e.operator vector<int>();    }
static vector<double> toDblVec(OptionEntry e) { return e.operator vector<double>(); }
static string         toStr   (OptionEntry e) { return e.operator string();         }
static int            toInt   (OptionEntry e) { return e.operator int();            }
static double         toDbl   (OptionEntry e) { return e.operator double();         }
static bool           toBool  (OptionEntry e) { return e.operator bool();           }

// Safe getters: return a default if the key is absent from the input file
static int    getInt (OptionsSection& s, const string& k, int    def=0)    { auto e=s.Get(k); return e.is_found_ ? toInt(e)    : def; }
static double getDbl (OptionsSection& s, const string& k, double def=0.)   { auto e=s.Get(k); return e.is_found_ ? toDbl(e)    : def; }
static bool   getBool(OptionsSection& s, const string& k, bool   def=false){ auto e=s.Get(k); return e.is_found_ ? toBool(e)   : def; }
static string getStr (OptionsSection& s, const string& k, const string& def=""){ auto e=s.Get(k); return e.is_found_ ? toStr(e) : def; }

static vector<double> getDblVec(OptionsSection& s, const string& k, vector<double> def={}){
    auto e = s.Get(k); return e.is_found_ ? toDblVec(e) : def;
}
static vector<int> getIntVec(OptionsSection& s, const string& k, vector<int> def={}){
    auto e = s.Get(k); return e.is_found_ ? toIntVec(e) : def;
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int mypeno, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &mypeno);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    if (argc < 2) {
        if (mypeno == 0)
            cerr << "Usage: mpirun -np <N> ./nlfem input.dat" << endl;
        MPI_Finalize();
        return 1;
    }

    string inputFile = argv[1];

    if (mypeno == 0) {
        cout << "========================================" << endl;
        cout << "  nlfem - Nonlinear FEM Solver"          << endl;
        cout << "  Author : Gokul G. Anugrah"             << endl;
        cout << "  Contact: gokulanugrah@gmail.com"        << endl;
        cout << "========================================" << endl;
        cout << "  MPI processes : " << numprocs           << endl;
        cout << "  Input file    : " << inputFile           << endl;
        cout << "========================================" << endl;
    }

    PreprocessorContext pContext(&argc, &argv, mypeno);
    OptionsParser parser(&pContext);
    parser.SetFile(inputFile);

    OptionsSection& strucSec    = parser.Section("Structural");
    OptionsSection& elemProps   = strucSec.Section("ElementProperties");
    OptionsSection& solverSec   = strucSec.Section("Solver");
    OptionsSection& loadSec     = strucSec.Section("Loading");
    OptionsSection& fsiSec      = strucSec.Section("FSI");
    OptionsSection& outSec      = strucSec.Section("Output");
    OptionsSection& unsteadySec = strucSec.Section("Unsteady");
    OptionsSection& restartSec  = strucSec.Section("Restart");
    OptionsSection& nlSec       = strucSec.Section("Nonlinear");

    linearTriElement solver;

    // MPI
    solver.mypeno_   = mypeno;
    solver.numprocs_ = numprocs;

    // ---- Element properties ----
    // ElementType uses vector syntax e.g. "(3)"
    {
        auto e = elemProps.Get("ElementType");
        if (e.is_found_) solver.elementTypes_ = toIntVec(e);
        else             solver.elementTypes_ = {3};  // default: MITC3
    }
    solver.defaultEModulus_  = getDbl(elemProps, "defaultEmod",      1.0);
    solver.defaultPoisRatio_ = getDbl(elemProps, "defaultPoisRatio", 0.3);
    // defaultDensity uses vector syntax e.g. "(0.001)"
    {
        auto e = elemProps.Get("defaultDensity");
        if (e.is_found_) solver.defaultDensity_ = toDblVec(e);
        else             solver.defaultDensity_ = {1.0};
    }
    solver.defaultThickness_ = getDbl(elemProps, "defaultThick",  0.001);
    solver.beamCSA_          = getDbl(elemProps, "defaultCSA",    0.0);   // optional for shells
    solver.surfaceFilename_  = getStr(elemProps, "surfaceFilename");
    solver.fem_BC_file_      = getStr(elemProps, "nodeFixedFile");
    solver.dirVecIndex_      = getInt(elemProps, "dirVecIndex",  -1);

    // Optional anisotropic moduli
    solver.inPlaneEmod_.push_back (getDbl(elemProps, "inPlaneEmod",  solver.defaultEModulus_));
    solver.outPlaneEmod_.push_back(getDbl(elemProps, "outPlaneEmod", solver.defaultEModulus_));

    // Optional forces file
    solver.noForceFile_ = true;
    {
        string ff = getStr(elemProps, "forcesFile");
        if (!ff.empty()) {
            solver.fem_force_file_ = ff;
            solver.noForceFile_ = false;
        }
    }

    // ---- Solver settings ----
    solver.isFSI_                 = getBool(solverSec, "isFSI",              false);
    solver.useAuxGeometry_        = getBool(solverSec, "useAuxGeometry",     false);
    solver.threeD_as_twoD_        = getBool(solverSec, "3Das2D",             false);
    solver.recordNodeOutput_      = getBool(solverSec, "recordNodeOutput",   false);
    solver.is_unsteady_           = getBool(solverSec, "isUnsteady",         false);
    solver.do_restart_            = getBool(solverSec, "doRestart",          false);
    solver.is_nonlinear_          = getBool(solverSec, "isNonLinear",        false);
    solver.read_normals_          = getBool(solverSec, "readNormals",        false);
    solver.prestrained_           = getBool(solverSec, "prestrained",        false);
    solver.relaxationFactor_      = getDbl (solverSec, "relaxFactor",        1.0);
    solver.newtonTolerance_       = getDbl (solverSec, "newTol",             1e-6);
    solver.timeIntegrationMethod_ = getStr (solverSec, "timeIntegrationMethod", "implicit");

    // ---- Loading (all optional) ----
    solver.bodyForces_      = getBool(loadSec, "BodyForces",      false);
    solver.pressureLoading_ = getBool(loadSec, "PressureLoading", false);
    solver.bodyForceMag_    = getDbl (loadSec, "ExternalBodyForceMag", 0.0);
    solver.gravity_         = getDblVec(loadSec, "Gravity", {0., 0., 0.});
    solver.pressureMag_     = getDbl (loadSec, "PressureMagnitude",    0.0);
    solver.pressureTime_    = getDbl (loadSec, "PressureLoadingTime",  0.0);

    // ---- FSI / auxiliary geometry (all optional) ----
    if (solver.isFSI_ || solver.useAuxGeometry_) {
        OptionsSection& auxGeo    = fsiSec.Section("AuxGeometry");
        OptionsSection& fsiSolver = fsiSec.Section("Solver");
        solver.auxFilename_           = getStr (auxGeo, "auxFilename");
        solver.rayLength_             = getDbl (auxGeo, "rayLength",          1.0);
        solver.geomThickFactor_       = getDbl (auxGeo, "geomThickFactor",    1.0);
        solver.revertGeometryNormals_ = getBool(auxGeo, "revertGeometryNormals", false);
        solver.enforceContact_        = getInt (auxGeo, "enforceContact",     0);
        solver.numForceInterpPts_     = getInt (auxGeo, "numForceInterpPts",  1);
        solver.auxForceInterpMethod_  = getStr (auxGeo, "auxForceInterpMethod",  "cloud");
        solver.auxDispInterpMethod_   = getStr (auxGeo, "geom2femInterpMethod",  "single");
        solver.forceInterpOrder_      = getStr (auxGeo, "auxForceInterpOrder",   "first");
        solver.dispInterpOrder_       = getStr (auxGeo, "auxDispInterpOrder",    "first");
        solver.fsiSkip_               = getInt (fsiSolver, "fsiSkip", 1);
    }

    // ---- Node output (optional) ----
    if (solver.recordNodeOutput_) {
        int nNodes = getInt(outSec, "numberOfNodes", 0);
        for (int i = 1; i <= nNodes; ++i) {
            ostringstream nodeKey, dofKey, freqKey;
            nodeKey << "node" << i;
            dofKey  << "DOF"  << i;
            freqKey << "freq" << i;
            solver.outputNodes_.push_back   (getInt(outSec, nodeKey.str(), 0));
            solver.outputNodeDOF_.push_back (getStr(outSec, dofKey.str(),  "all"));
            solver.outputNodeFreq_.push_back(getInt(outSec, freqKey.str(), 1));
        }
    }

    // ---- Unsteady (optional) ----
    solver.nm_gamma_     = getDbl(unsteadySec, "gamma",     0.5);
    solver.nm_beta_      = getDbl(unsteadySec, "beta",      0.25);
    solver.rampFactor_   = getDbl(unsteadySec, "rampFactor",1.0);
    solver.dt_           = getDbl(unsteadySec, "dt",        1e-3);
    solver.numTimesteps_ = getInt(unsteadySec, "timeSteps", 0);
    solver.ntSkip_       = getInt(unsteadySec, "ntSkip",    1);

    // ---- Restart (optional) ----
    solver.restartSkip_     = getInt(restartSec, "restartSkip", 0);
    solver.restartTimestep_ = getInt(restartSec, "restartTS",   0);

    // ---- Nonlinear (optional) ----
    solver.numLoadSteps_ = getInt(nlSec, "numLoadSteps", 1);
    solver.nonlinearIts_ = getInt(nlSec, "nonlinearIts", 100);

    // ---- Internal defaults ----
    solver.parallelSolve_         = true;
    solver.pctype_                = "gamg";
    solver.shiftGeometry_         = false;
    solver.loadsByCompID_         = false;
    solver.regionalBC_            = false;
    solver.auxCornerTreatment_    = false;
    solver.recordCableOutput_     = false;
    solver.one2one_               = false;
    solver.revertNormalsByCompID_ = false;
    solver.backtrack_             = false;
    solver.firstTimeStep_         = true;
    solver.firstTimeComputingK_   = true;

    // ---- Run ----
    solver.Setup_fem();

    if (solver.is_unsteady_) {
        for (int nt = 1; nt <= solver.numTimesteps_; ++nt) {
            solver.Step_structural_solver(nt);
            if (mypeno == 0 && nt % solver.ntSkip_ == 0)
                cout << "  Step " << nt << "/" << solver.numTimesteps_
                     << "  t = " << solver.dt_ * nt << endl;
        }
    } else {
        int nt = 0;
        solver.Step_structural_solver(nt);
    }

    solver.Finalize_fem();

    if (mypeno == 0) {
        cout << "========================================" << endl;
        cout << "  nlfem: Simulation complete."           << endl;
        cout << "========================================" << endl;
    }

    MPI_Finalize();
    return 0;
}
