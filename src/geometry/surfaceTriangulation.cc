/*
 * nlfem - Nonlinear FEM Solver
 * Author  : Gokul G. Anugrah
 * Contact : gokulanugrah@gmail.com
 *
 * surfaceTriangulation.cc — VTK mesh reader and surface geometry utilities.
 *   Reads unstructured VTK meshes, handles FSI auxiliary geometry,
 *   and provides ray-triangle intersection tests for geometry mapping.
 */
#include "surfaceTriangulation.h"
#include "screen.h"
#include "string_utilities.h"
#include "linear_algebra.h"
using namespace std;
namespace KARMA { 
  
  void surfaceTriangulation::smartTriangulationRead(vector<double>& nodes,               //xyz
                                                    vector<int>&    faces,               //connectivity
                                                    string&         filename,            //VTK filename
                                                    int&            globalNodeNum,       //global number of unqiue nodes
                                                    int&            globalFaceNum,       //global number of elements
                                                    vector<int>&    elementTypes,        //which element types were specified in input file
                                                    vector<int>&    elementTypeVector,   //tag each global element ID with its element type
                                                    vector<int>&    ndofPerNode,         //unique list of ndof per node
                                                    vector<int>&    componentIDs,        //element component IDs
                                                    bool&           revertNormalsByComp) //should we revert normals on the listed comp IDs
  {

    //this reader will identify the types of elements being read and sort accordingly.
    //i.e., a vector will be returned that has the ndof/node for a node on each elemenet type.
    //if a node is shared between 2 types of elements, the one with higher dof will win.

    cout << "-----------------------------------" << endl;
    Screen::MasterInfo("Reading FEM mesh from VTK file --- Smart Reader");
    Screen::MasterInfo("VTK filename is " + filename);
    ifstream file;
    file.open(filename);
    if(file.fail()) Screen::MasterError("File not found in surfaceTriangulation::readTriangulation");

    string cdummy;
    file >> cdummy;//#
    file >> cdummy;//vtk
    file >> cdummy;//DataFile
    file >> cdummy;//Version
    file >> cdummy;//3.0
    file >> cdummy;//vtk
    file >> cdummy;//output
    file >> cdummy;//ASCII
    file >> cdummy;//DATASET
    file >> cdummy;//POLYDATA
    file >> cdummy;//POINTS
    file >> globalNodeNum;
    Screen::MasterInfo("FEM mesh has " + to_string(globalNodeNum) + " nodes");
    file >> cdummy;//float

    //allocate nodes based off read value of globalNodeNum
    nodes.resize(globalNodeNum*3,0.);

    //read nodal locations
    for (int n=0;n<globalNodeNum;++n) 
      {
        for (int idir=0;idir<3;++idir) 
          {
            int index = idir + 3*n;
            file >> nodes[index];
          }
      }

    int numberOfElements_type1 = 0;
    int numberOfElements_type2 = 0;
    file >> cdummy;//LINES or POLYGONS
    file >> numberOfElements_type1;//number of elements    

    Screen::MasterInfo("FEM mesh has " + to_string(numberOfElements_type1) + " " + cdummy);

    int nodesPerElement1;
    int ndofPerElement1;
    if (cdummy=="LINES")
      {
        //determine if beams or cables
        bool beamElements = false;
        for (int i=0;i<elementTypes.size();++i)
          {
            if (elementTypes[i] == 4)
              {
                beamElements = true;
              }
          }

        if (beamElements)
          {
            nodesPerElement1 = 2;
            ndofPerElement1  = 2;
          }
        else//cable elements
          {
            nodesPerElement1 = 2;
            ndofPerElement1  = 3;
          }
      }
    else if (cdummy=="POLYGONS")
      {
        //treat all triangular elements same
        nodesPerElement1 = 3;
        ndofPerElement1  = 5;
      }
    else
      {
        Screen::MasterError("LINES or POLYGONS not specified in VTK");
      }
    file >> cdummy;//numberOfElements_type1*3 for LINES and numberOfElements_type1*4 for POLYGONS
    ndofPerNode.resize(globalNodeNum,-1);
    //allocation for the connectivity corresponding to element type 1
    vector<int> faces_element1(nodesPerElement1*numberOfElements_type1,-1);
    int idata;
    globalFaceNum = 0;
    //reading face connectivity
    for (int f=0;f<numberOfElements_type1;++f) 
      {
        globalFaceNum += 1;
        file >> idata;//dummy 2 for LINES or a 3 for POLYGONS
        for (int n=0;n<nodesPerElement1;++n)
          {
            int index = n + nodesPerElement1*f;
            file >> faces_element1[index];
            faces.push_back(faces_element1[index]);
          }
        //sort out elementTypeVector and ndofPerNode vector
        //set element types in vector
        elementTypeVector.push_back(elementTypes[0]);

        for (int n1=0;n1<nodesPerElement1;++n1)
          {
            int node1 = faces_element1[n1 + nodesPerElement1*f];
            ndofPerNode[node1] = ndofPerElement1;
          }
      }
    //test if there is anything left to read --- i.e., any other elements in the simulation
    if (elementTypes.size() > 1)
      {
        Screen::MasterInfo("Multiple elements specified in input file, continuing to read");

        file >> cdummy;//read LINES or POLYGONS      
        int nodesPerElement2;
        int ndofPerElement2;
        //final check
        if (cdummy != "LINES" && cdummy != "POLYGONS")
          {
            Screen::MasterWarning("Something is left in the FEM mesh VTK, but its not LINES or POLYGONS --- leaving routine");
            return;
          }
        if (cdummy == "POLYGONS")
          {
            nodesPerElement2 = 3;
            ndofPerElement2 = 5;
          }
        else if (cdummy == "LINES")
          {
            //assuming cables and beams are never run togther
            Screen::MasterError("LINES has been specified twice in the VTK");
          }

        file >> numberOfElements_type2;//number of elements with type 2
        Screen::MasterInfo("FEM mesh has " + to_string(numberOfElements_type2) + " " + cdummy);

        file >> cdummy;//dummy number

        //reading face connectivity
        vector<int> faces_element2(nodesPerElement2*numberOfElements_type2,-1);
        for (int f=0;f<numberOfElements_type2;++f)
          {
            globalFaceNum += 1;
            file >> idata;//dummy number
            for (int n=0;n<nodesPerElement2;++n)
              {
                int index = n + nodesPerElement2*f;
                file >> faces_element2[index];
                faces.push_back(faces_element2[index]);
              }
          }
        //the assumption here is that the element that has nodes with more DOF comes second
        //this allows the second loop to overwrite the ndofPerNode for shared nodes
        for (int e2=0;e2<numberOfElements_type2;++e2)
          {
            //set element types in vector with offset
            elementTypeVector.push_back(elementTypes[1]);

            for (int n2=0;n2<nodesPerElement2;++n2)
              {
                int node2 = faces_element2[n2 + nodesPerElement2*e2];
                ndofPerNode[node2] = ndofPerElement2;
              }
          }
      }

    if (revertNormalsByComp)
      {

        Screen::MasterInfo("Reading in component IDs for Reverting Normals on Marked Faces");
        componentIDs.resize(numberOfElements_type1+numberOfElements_type2,-200);

        file >> cdummy;//CELL_DATA
        file >> cdummy;//total number of elements
        file >> cdummy;//SCALARS
        file >> cdummy;//Components
        file >> cdummy;//int
        file >> cdummy;//LOOKUP_TABLE
        file >> cdummy;//default

        for (int f=0;f<(numberOfElements_type1+numberOfElements_type2);++f)
          {
            file >> componentIDs[f];
          }

      }

    Screen::MasterInfo("FEM mesh reading completed");
    cout << "-----------------------------------" << endl;
    file.close();

  } 
  
  void surfaceTriangulation::readTriangulation(vector<double>& nodes,
                                               vector<int>&    faces,
                                               string&         filename,
                                               int&            globalNodeNum,
                                               int&            globalFaceNum,
                                               int&            nodesPerElement)
  {
    cout << "-----------------------------------" << endl;
    Screen::MasterInfo("Reading FEM mesh");
    Screen::MasterInfo("Surface filename is " + filename);
    ifstream file;
    file.open(filename);
    if(file.fail()) Screen::MasterError("File not found in surfaceTriangulation::readTriangulation");

    char cdummy[100];
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> globalNodeNum;
    Screen::MasterInfo("FEM mesh has " + to_string(globalNodeNum) + " nodes");
    file >> cdummy;

    nodes.resize(globalNodeNum*3,0.);
    //reading nodal locations
    for (int n=0;n<globalNodeNum;++n) 
      {
        for (int idir=0;idir<3;++idir) 
          {
            int index = idir + 3*n;
            file >> nodes[index];
          }
      }

    file >> cdummy;
    file >> globalFaceNum;
    Screen::MasterInfo("FEM mesh has " + to_string(globalFaceNum) + " elements");
    file >> cdummy;    

    faces.resize(globalFaceNum*nodesPerElement,0);
    int idata;

    //reading face connectivity
    for (int f=0;f<globalFaceNum;++f) 
      {
        file >> idata;     
        for (int idir=0;idir<nodesPerElement;++idir)
          {
            int index = idir + nodesPerElement*f;
            file >> faces[index];
          }
      }

    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;

    components_.resize(globalFaceNum,0);

    //read components
    for (int f=0;f<globalFaceNum;++f) 
      {
        file >> components_[f];
      }

    Screen::MasterInfo("FEM mesh reading completed");
    cout << "-----------------------------------" << endl;
    file.close();

  } 

  void surfaceTriangulation::read_preStrained_nodeLocs(vector<vector<double> >& nodes,
                                                       string filename)
  {
    cout << "-----------------------------------" << endl;
    Screen::MasterInfo("Reading pre-strained FEM node locations");
    Screen::MasterInfo("Pre-strained node locations are in " + filename);
    ifstream file;
    file.open(filename);
    if(file.fail()) Screen::MasterError("Bad file for pre-strained FEM");

    char cdummy[100];
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> globalNodeNum_;
    Screen::MasterInfo("VTK file contains " + to_string(globalNodeNum_) + " nodes");
    file >> cdummy;

    vector<double> threeVec(3,0.);
    nodes.resize(globalNodeNum_,threeVec);
    //reading nodal locations
    for (int n=0;n<globalNodeNum_;++n) 
      {
        for (int idir=0;idir<3;++idir) 
          {
            file >> nodes[n][idir];
          }
      }

    Screen::MasterInfo("Pre-strained nodal locations read");
    cout << "-----------------------------------" << endl;
    file.close();

  } 

  void surfaceTriangulation::readTriangulation_AUX(vector<double>& auxnodes,
                                                   vector<int>&    auxfaces,
                                                   string&         filename,
                                                   int&            globalNodeNum,
                                                   int&            globalFaceNum,
                                                   int&            nodesPerElement,
                                                   vector<int>&    isauxnode)
  
  {
    Screen::MasterInfo("Reading Geometry File: "+filename);
    ifstream file;
    file.open(filename);
    if(file.fail()) Screen::MasterError("Bad file in geometry reader :: readTriangulation_AUX");

    char cdummy[100];
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> globalNodeNum;
    Screen::MasterInfo("AUX VTK file contains " + to_string(globalNodeNum) + " nodes");
    file >> cdummy;

    auxnodes.resize(globalNodeNum*3,0.);
    //reading nodal locations
    for (int n=0;n<globalNodeNum;++n)
      {
        for (int idir=0;idir<3;++idir)
          {
            int index = idir + 3*n;
            file >> auxnodes[index];
          }
      }

    file >> cdummy;
    file >> globalFaceNum; 
    Screen::MasterInfo("AUX VTK file contains " + to_string(globalFaceNum) + " faces");  
    file >> cdummy;  

    //reading face connectivity
    auxfaces.resize(globalFaceNum*nodesPerElement,0);
    int idata;
    for (int f=0;f<globalFaceNum;++f) 
      {
        file >> idata;     
        for (int idir=0;idir<nodesPerElement;++idir)
          {
            int index = idir + nodesPerElement*f;
            file >> auxfaces[index];
          }
      }

    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;

    //read components
    vector<int> auxcomponents(globalFaceNum,0);
    for (int f=0;f<globalFaceNum;++f) 
      {
        file >> auxcomponents[f];
      }
    //identify moving nodes and fill aux nodes
    Screen::MasterInfo("Identifying auxillary nodes");
    isauxnode.resize(globalNodeNum,0);
    for (int f=0;f<globalFaceNum;++f)
      {
        int store99 = auxcomponents[f];
        if (abs(store99) == 99)
          {
            for (int n=0;n<nodesPerElement;++n)
              {
                int node = auxfaces[n+nodesPerElement*f];
                isauxnode[node] = store99;
              }
          }
      }

    file.close();

  } 
  
  void surfaceTriangulation::readTriangulation_AUX_FSI(vector<double>& auxnodes,
                                                       string&         filename)
  {

    Screen::MasterInfo("Reading duplicate node locations from vtk: "+filename);
    ifstream file;
    file.open(filename);
    if(file.fail()) Screen::MasterError("Bad file in geometry reader :: readTriangulation_AUX_FSI");

    int globalNodeNum_loc;
    int globalFaceNum_loc;

    char cdummy[100];
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> cdummy;
    file >> globalNodeNum_loc;
    Screen::MasterInfo("VTK file contains " + to_string(globalNodeNum_loc) + " nodes");
    file >> cdummy;

    auxnodes.resize(globalNodeNum_loc*3,0.);
    //reading nodal locations
    for (int n=0;n<globalNodeNum_loc;++n)
      {
        for (int idir=0;idir<3;++idir)
          {
            int index = idir + 3*n;
            file >> auxnodes[index];
          }
      }

    file.close();

  } 

  void surfaceTriangulation::writeTriangulation(string filename,
                                                vector<double>& Nodes,
                                                vector<int>&    Faces,
                                                int& dimLoc,
                                                int& nodeNum,
                                                int& faceNum)
  {
    Screen::MasterInfo("Writing nodal positions to VTK");
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
    file << nodeNum;
    file << " float";
    file << endl;
    int DIM = 3;
    //writing nodal locations
    for (int n=0;n<nodeNum;++n) 
      {
        for (int idir=0;idir<DIM;++idir) 
          {
            file << Nodes[n*DIM+idir];
            file << " ";
          }
        file << endl;
      }

    if (dimLoc == 2)
      {
        file << "LINES ";
      }
    else
      {
        file << "POLYGONS "; 
      }
    file << faceNum;
    file << " ";  
    int idata;
    if (dimLoc == 2)
      {
        file << faceNum*3;
        idata = 2;
      }
    else
      {
        file << faceNum*4;
        idata = 3;
      }
    file << endl;    

    //writing face connectivity
    for (int f=0;f<faceNum;++f) 
      {
        file << idata;
        file << " ";    
        for (int idir=0;idir<3;++idir)
          {
            int index = idir + 3*f;
            file << Faces[index];
            file << " ";
          }
        file << endl;
      }

    file.close();

  }
  
  void surfaceTriangulation::smartWriteTriangulation(string filename,
                                                     vector<double>& Nodes,
                                                     vector<int>&    Faces,
                                                     int& dimLoc,
                                                     int& nodeNum,
                                                     int& faceNum,
                                                     int& numberOfElementTypes,
                                                     vector<int>& nnodesPerElement)
  {
    Screen::MasterInfo("Writing nodal positions to VTK --- Smart Writer");
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
    file << nodeNum;
    file << " float";
    file << endl;
    int DIM = 3;
    //writing nodal locations
    for (int n=0;n<nodeNum;++n) 
      {
        for (int idir=0;idir<DIM;++idir) 
          {
            file << Nodes[n*DIM+idir];
            file << " ";
          }
        file << endl;
      }
    //write connectivity
    int idata;
    if (numberOfElementTypes == 1)
      {
        if (nnodesPerElement[0] == 2)
          {
            file << "LINES ";
            file << faceNum;
            file << " ";
            file << faceNum*3;
            file << endl;
            idata = 2;
          }
        else
          {
            file << "POLYGONS ";
            file << faceNum;
            file << " ";
            file << faceNum*4;
            file << endl;
            idata = 3;
          }

        for (int f=0;f<faceNum;++f)
          {
            file << nnodesPerElement[f] << " ";                           
            for (int idir=0;idir<nnodesPerElement[f];++idir)
              {            
                int index = idir + nnodesPerElement[f]*f;
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
        for (int n=0;n<faceNum;++n)
          {
            if (nnodesPerElement[n] == 2)
              {
                numfaces1 +=1;
              }
            if (nnodesPerElement[n] == 3)
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
            if (nnodesPerElement[f] == 2)
              {
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
            offset += nnodesPerElement[f];
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
            if (nnodesPerElement[f] == 3)
              {
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
            offset += nnodesPerElement[f];
          }
      }

    file.close();

  }
  
  void surfaceTriangulation::smartWriteTriangulation_withForces(string filename,
                                                                vector<double>& Nodes,
                                                                vector<int>&    Faces,
                                                                int& dimLoc,
                                                                int& nodeNum,
                                                                int& faceNum,
                                                                int& numberOfElementTypes,
                                                                vector<int>& nnodesPerElement,
                                                                vector<double>& forces,
                                                                vector<int>& orderedList,
                                                                vector<double>& disps,
								vector<double>& vonMises)
  {
    Screen::MasterInfo("Writing nodal positions to VTK --- Smart Writer");
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
    file << nodeNum;
    file << " float";
    file << endl;
    int DIM = 3;
    //writing nodal locations
    for (int n=0;n<nodeNum;++n) 
      {
        for (int idir=0;idir<DIM;++idir) 
          {
            file << Nodes[n*DIM+idir];
            file << " ";
          }
        file << endl;
      }
    //write connectivity
    int idata;
    if (numberOfElementTypes == 1)
      {
        if (nnodesPerElement[0] == 2)
          {
            file << "LINES ";
            file << faceNum;
            file << " ";
            file << faceNum*3;
            file << endl;
            idata = 2;
          }
        else
          {
            file << "POLYGONS ";
            file << faceNum;
            file << " ";
            file << faceNum*4;
            file << endl;
            idata = 3;
          }

        for (int f=0;f<faceNum;++f)
          {
            file << nnodesPerElement[f] << " ";                           
            for (int idir=0;idir<nnodesPerElement[f];++idir)
              {            
                int index = idir + nnodesPerElement[f]*f;
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
        for (int n=0;n<faceNum;++n)
          {
            if (nnodesPerElement[n] == 2)
              {
                numfaces1 +=1;
              }
            if (nnodesPerElement[n] == 3)
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
            if (nnodesPerElement[f] == 2)
              {
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
            offset += nnodesPerElement[f];
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
            if (nnodesPerElement[f] == 3)
              {
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
            offset += nnodesPerElement[f];
          }
      }
    //write forces at the centroid for each element
    file << "CELL_DATA ";
    file << faceNum;
    file << endl;
//    file << "SCALARS ";
//    file << "xForces ";
//    file << "float";
//    file << endl;
//    file << "LOOKUP_TABLE ";
//    file << "default";
//    file << endl;
//    for (int e=0;e<faceNum;++e)
//      {
//        for (int oe=0;oe<faceNum;++oe)
//          {
//            int elementID = orderedList[oe];
//            if (e == elementID)
//              {
//                file << forces[oe*3 + 0];
//                file << endl;
//              }
//          }
//      }
//    file << "SCALARS ";
//    file << "yForces ";
//    file << "float";
//    file << endl;
//    file << "LOOKUP_TABLE ";
//    file << "default";
//    file << endl;
//    for (int e=0;e<faceNum;++e)
//      {
//        for (int oe=0;oe<faceNum;++oe)
//          {
//            int elementID = orderedList[oe];
//            if (e == elementID)
//              {
//                file << forces[oe*3 + 1];
//                file << endl;
//              }
//          }
//      }
//    file << "SCALARS ";
//    file << "zForces ";
//    file << "float";
//    file << endl;
//    file << "LOOKUP_TABLE ";
//    file << "default";
//    file << endl;
//    for (int e=0;e<faceNum;++e)
//      {
//        for (int oe=0;oe<faceNum;++oe)
//          {
//            int elementID = orderedList[oe];
//            if (e == elementID)
//              {
//                file << forces[oe*3 + 2];
//                file << endl;
//              }
//          }
//      }
    file << "SCALARS ";
    file << "xTotDisp ";
    file << "float";
    file << endl;
    file << "LOOKUP_TABLE ";
    file << "default";
    file << endl;

    for (int e=0;e<faceNum;++e) file << disps[e*3 + 0] << endl;

    file << "SCALARS ";
    file << "yTotDisp ";
    file << "float";
    file << endl;
    file << "LOOKUP_TABLE ";
    file << "default";
    file << endl;

    for (int e=0;e<faceNum;++e) file << disps[e*3 + 1] << endl;

    file << "SCALARS ";
    file << "zTotDisp ";
    file << "float";
    file << endl;
    file << "LOOKUP_TABLE ";
    file << "default";
    file << endl;

    for (int e=0;e<faceNum;++e) file << disps[e*3 + 2] << endl;

    file << "SCALARS ";
    file << "vonMises ";
    file << "float";
    file << endl;
    file << "LOOKUP_TABLE ";
    file << "default";
    file << endl;

    for (int e=0;e<faceNum;++e) file << vonMises[e] << endl;

    file.close();

  }
  
  void surfaceTriangulation::writeTriangulation_parallel_unformatted(string& filename,
								     int& globalNodeCount,
								     int& globalFaceCount,
								     vector<double>& localFaceForces,
								     vector<double>& localFaceDispls,
								     map<int,int>& nodeLocal2Global,
								     map<int,int>& elemLocal2Global,
								     int& localNodeCount,
								     int& localFaceCount,
								     vector<double>& nodes,
								     vector<int>& faces,
								     vector<int>& nnodesPerElement,
								     vector<double>& vonMises)
  {

    //this will write an unformatted file that needs a converter to turn into a vtk
    //only writing xyz forces and displacements
    Screen::MasterInfo("Writing unformatted file: "+filename);

    int ierr;
    //all processes open file
    MPI_File fh;
    ierr = MPI_File_open(MPI_COMM_WORLD,&filename[0],MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
    if (ierr != 0) Screen::MasterError("Error in MPI_File_open, writeTriangulation_parallel_unformatted::surfaceTriangulation.cc, ierr = "+to_string(ierr));

    //set view
    MPI_Offset disp = 0;
    ierr = MPI_File_set_view(fh,disp,MPI::DOUBLE,MPI::DOUBLE,"native",MPI_INFO_NULL);
    if (ierr != 0) Screen::MasterError("Error in MPI_File_set_view, writeTriangulation_parallel_unformatted::surfaceTriangulation.cc, ierr = "+to_string(ierr));
    MPI_Status status;

    int line_offset = 1;
    MPI_Offset line = line_offset;

    //write global node and face count
    double num1 = float(globalNodeCount);
    ierr = MPI_File_write_at(fh,line,&num1,1,MPI::DOUBLE,&status);
    if (ierr != 0) Screen::MasterError("Error in MPI_File_write_at, writeTriangulation_parallel_unformatted::surfaceTriangulation.cc::write 1, ierr = "+to_string(ierr));//write 1
    line_offset += 1;

    line = line_offset;
    double num2 = float(globalFaceCount);
    ierr = MPI_File_write_at(fh,line,&num2,1,MPI::DOUBLE,&status);
    if (ierr != 0) Screen::MasterError("Error in MPI_File_write_at, writeTriangulation_parallel_unformatted::surfaceTriangulation.cc::write 2, ierr = "+to_string(ierr));//write 2
    line_offset += 1;

    vector<double> pos(3,0.);
    //write node xyz locations
    for (int n=0;n<localNodeCount;++n)
      {	

	pos[0] = nodes[n*3+0];
	pos[1] = nodes[n*3+1];
	pos[2] = nodes[n*3+2];

	int globID = nodeLocal2Global[n];
	line = globID*3 + line_offset;
	ierr = MPI_File_write_at(fh,line,&pos[0],3,MPI::DOUBLE,&status);
	if (ierr != 0) Screen::MasterError("Error in MPI_File_write_at, writeTriangulation_parallel_unformatted::surfaceTriangulation.cc::write 3, ierr = "+to_string(ierr));//write 3

      }
    line_offset += globalNodeCount*3;

    vector<double> nodeVec(4,0.);
    int offset = 0;
    //write nnodesPerElement,face nodes
    for (int e=0;e<localFaceCount;++e)
      {

	int nnodes_loc = nnodesPerElement[e];
	int globID     = elemLocal2Global[e];
	nodeVec[0] = float(nnodes_loc);
	for (int n=0;n<nnodes_loc;++n) nodeVec[n+1] = float(nodeLocal2Global[faces[n + offset]]);
	offset += nnodes_loc;

	line = globID*4 + line_offset;
	ierr = MPI_File_write_at(fh,line,&nodeVec[0],(nnodes_loc+1),MPI::DOUBLE,&status);
	if (ierr != 0) Screen::MasterError("Error in MPI_File_write_at, writeTriangulation_parallel_unformatted::surfaceTriangulation.cc::write 4, ierr = "+to_string(ierr));//write 4

      }
    line_offset += globalFaceCount*4;

    //write face forces
    for (int e=0;e<localFaceCount;++e)
      {
	int globID = elemLocal2Global[e];
	pos[0] = localFaceForces[e*3 + 0];
	pos[1] = localFaceForces[e*3 + 1];
	pos[2] = localFaceForces[e*3 + 2];

	line = globID*3 + line_offset;
	ierr = MPI_File_write_at(fh,line,&pos[0],3,MPI::DOUBLE,&status);
	if (ierr != 0) Screen::MasterError("Error in MPI_File_write_at, writeTriangulation_parallel_unformatted::surfaceTriangulation.cc::write 5, ierr = "+to_string(ierr));//write 5

      }
    line_offset += globalFaceCount*3;

    //write face disps
    for (int e=0;e<localFaceCount;++e)
      {
	int globID = elemLocal2Global[e];
	pos[0] = localFaceDispls[e*3 + 0];
	pos[1] = localFaceDispls[e*3 + 1];
	pos[2] = localFaceDispls[e*3 + 2];

	line = globID*3 + line_offset;
	ierr = MPI_File_write_at(fh,line,&pos[0],3,MPI::DOUBLE,&status);
	if (ierr != 0) Screen::MasterError("Error in MPI_File_write_at, writeTriangulation_parallel_unformatted::surfaceTriangulation.cc::write 6, ierr = "+to_string(ierr));//write 6

      }

    line_offset += globalFaceCount*3;

    //write von Mises stresses
    for (int e=0;e<localFaceCount;++e)
      {
	int globID = elemLocal2Global[e];
	line = globID + line_offset;
	double vonmises_loc = vonMises[e];
	ierr = MPI_File_write_at(fh,line,&vonmises_loc,1,MPI::DOUBLE,&status);
	if (ierr != 0) Screen::MasterError("Error in MPI_File_write_at, writeTriangulation_parallel_unformatted::surfaceTriangulation.cc::write 8, ierr = "+to_string(ierr));//write 8

      }

    //close file
    ierr = MPI_File_close(&fh);
    if (ierr != 0) Screen::MasterError("Error in MPI_File_close, writeRestart::fem.cc, ierr = "+to_string(ierr));

  }
  
  void surfaceTriangulation::writeTriangulation_data(string filename,
                                                     vector<vector<double> >& Nodes,
                                                     vector<int>& Faces,
                                                     int& dimLoc,
                                                     vector<vector<double> >& forces,
                                                     vector<vector<double> >& disp,
                                                     vector<vector<double> >& vel,
                                                     vector<double> & area,
                                                     vector<vector<double> >& Vnt,
                                                     vector<vector<double> >& V1t,
                                                     vector<vector<double> >& V2t)
  {
    //dumpMat(forces,"forces");
    //wait();
    Screen::MasterInfo("Writing FEM data to vtk");
    ofstream file;
    file.open(filename);
    if(file.fail()) Screen::MasterError("Bad file in writeTriangulation function");
    file.precision(10);
    file.setf(ios::fixed);
    file.setf(ios::showpoint);

    int nodeNum = Nodes.size();
    int faceNum = Faces.size();
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
    file << "POLYDATA ";
    file << endl;
    file << "POINTS ";
    file << nodeNum;
    file << " float";
    file << endl;
    int DIM = 3;
    //writing nodal locations
    for (int n=0;n<nodeNum;++n) 
      {
        for (int idir=0;idir<DIM;++idir) 
          {
            file << Nodes[n][idir];
            file << " ";
          }
        file << endl;
      }

    if (dimLoc == 2)
      {
        file << "LINES ";
      }
    else
      {
        file << "POLYGONS "; 
      }
    file << faceNum;
    file << " ";  
    int idata;
    if (dimLoc == 2)
      {
        file << faceNum*3;
        idata = 2;
      }
    else
      {
        file << faceNum*4;
        idata = 3;
      }
    file << endl;    

    //writing face connectivity
    for (int f=0;f<faceNum;++f) 
      {
        file << idata;     
        file << " ";    
        for (int idir=0;idir<dimLoc;++idir)
          {
            int index = idir + dimLoc*f;
            file << Faces[index];
            file << " ";
          }
        file << endl;
      }

    if (dimLoc == 2)
      {
        //write force at face
        file << "CELL_DATA ";
        file << faceNum;
        file << endl;
        file << "SCALARS ";
        file << "yForces ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f) 
          {     
            double force = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //force += forces[Faces[f][idir]][0]/dimLoc;
              }
            file << force;
            file << endl;
          }

        file << "SCALARS ";
        file << "Moments ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f) 
          {     
            double force = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //force += forces[Faces[f][idir]][1]/dimLoc;
              }
            file << force;
            file << endl;
          }
      }
    else if (dimLoc == 3)
      {
        file << "CELL_DATA ";
        file << faceNum;
        file << endl;
        file << "SCALARS ";
        file << "xForces ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f) 
          {     
            double force = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //force += forces[Faces[f][idir]][0]/dimLoc;
              }
            file << force;
            file << endl;
          }      

        file << "SCALARS ";
        file << "yForces ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f) 
          {     
            double force = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //force += forces[Faces[f][idir]][1]/dimLoc;
              }
            file << force;
            file << endl;
          }

        file << "SCALARS ";
        file << "zForces ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f) 
          {     
            double force = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //force += forces[Faces[f][idir]][2]/dimLoc;
              }
            file << force;
            file << endl;
          }

        file << "SCALARS ";
        file << "alphaMoment ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f) 
          {     
            double force = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //force += forces[Faces[f][idir]][3]/dimLoc;
              }
            file << force;
            file << endl;
          }

        file << "SCALARS ";
        file << "betaMoment ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f) 
          {     
            double force = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //force += forces[Faces[f][idir]][4]/dimLoc;
              }
            file << force;
            file << endl;
          }

        //file << "SCALARS ";
        //file << "Pressure ";
        //file << "float";
        //file << endl;
        //file << "LOOKUP_TABLE ";
        //file << "default";
        //file << endl;
        ////faces
        //for (int f=0;f<faceNum;++f) 
        //  {     
        //    double pres = 0.;
        //    //face nodes
        //    for (int fn=0;fn<dimLoc;++fn)
        //      {
        //        //dir
        //        for (int idir=0;idir<dimLoc;++idir)
        //          {
        //            pres += forces[ Faces[f][fn] ][idir] * Vnt[ Faces[f][fn] ][idir];
        //          }
        //        pres *= 0.33333333333333;
        //      }
        //    pres *= area[f];
        //    file << pres;
        //    file << endl;
        //  }
        file << "SCALARS ";
        file << "fx ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        //faces
        for (int f=0;f<faceNum;++f) 
          {     
            double pres = 0.;
            //face nodes
            for (int fn=0;fn<dimLoc;++fn)
              {
                //dir
                //pres += forces[ Faces[f][fn] ][0];
              }
            pres *= 0.3333333333333333/area[f];

            file << pres;
            file << endl;
          }
        file << "SCALARS ";
        file << "fy ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        //faces
        for (int f=0;f<faceNum;++f) 
          {     
            double pres = 0.;
            //face nodes
            for (int fn=0;fn<dimLoc;++fn)
              {
                //dir
                //pres += forces[ Faces[f][fn] ][1];
              }
            pres *= 0.3333333333333333/area[f];

            file << pres;
            file << endl;
          }
        file << "SCALARS ";
        file << "fz ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        //faces
        for (int f=0;f<faceNum;++f) 
          {     
            double pres = 0.;
            //face nodes
            for (int fn=0;fn<dimLoc;++fn)
              {
                //dir
                //pres += forces[ Faces[f][fn] ][2];
              }
            pres *= 0.3333333333333333/area[f];

            file << pres;
            file << endl;
          }

      }
    //write total displacements
    if (dimLoc == 3)
      {
        file << "SCALARS ";
        file << "udisp ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f)
          {
            double Displ = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //Displ += disp[Faces[f][idir]][0]/3.;
              }
            file << Displ;
            file << endl;
          }

        file << "SCALARS ";
        file << "vdisp ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f)
          {
            double Displ = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //Displ += disp[Faces[f][idir]][1]/3.;
              }
            file << Displ;
            file << endl;
          }

        file << "SCALARS ";
        file << "wdisp ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f)
          {
            double Displ = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //Displ += disp[Faces[f][idir]][2]/3.;
              }
            file << Displ;
            file << endl;
          }

        file << "SCALARS ";
        file << "alpha ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f)
          {
            double Displ = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //Displ += disp[Faces[f][idir]][3]/3.;
              }
            file << Displ;
            file << endl;
          }

        file << "SCALARS ";
        file << "beta ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f)
          {
            double Displ = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //Displ += disp[Faces[f][idir]][4]/3.;
              }
            file << Displ;
            file << endl;
          }
      }
    else if (dimLoc == 2)
      {

        file << "SCALARS ";
        file << "vdisp ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f)
          {
            double Displ = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //Displ += disp[Faces[f][idir]][0]/2.;
              }
            file << Displ;
            file << endl;
          }

        file << "SCALARS ";
        file << "phi ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f)
          {
            double Displ = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //Displ += disp[Faces[f][idir]][1]/2.;
              }
            file << Displ;
            file << endl;
          }
      }

    //face velocities
    if (dimLoc == 3)
      {
        file << "SCALARS ";
        file << "uVel ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f)
          {
            double Velo  = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //Velo += vel[Faces[f][idir]][0]/2.;
              }
            file << Velo;
            file << endl;
          }

        file << "SCALARS ";
        file << "vVel ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f)
          {
            double Velo  = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //Velo += vel[Faces[f][idir]][1]/2.;
              }
            file << Velo;
            file << endl;
          }

        file << "SCALARS ";
        file << "wVel ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f)
          {
            double Velo  = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //Velo += vel[Faces[f][idir]][2]/2.;
              }
            file << Velo;
            file << endl;
          }
      }
    else if (dimLoc == 2)
      {
        file << "SCALARS ";
        file << "vVel ";
        file << "float";
        file << endl;
        file << "LOOKUP_TABLE ";
        file << "default";
        file << endl;

        for (int f=0;f<faceNum;++f)
          {
            double Velo  = 0.;
            for (int idir=0;idir<dimLoc;++idir)
              {
                //Velo += vel[Faces[f][idir]][1]/2.;
              }
            file << Velo;
            file << endl;
          }
      }

    //write normal vectors
    file << "POINT_DATA ";
    file << nodeNum;
    file << endl;
    file << "VECTORS ";
    file << "normalVecs ";
    file << "double";
    file << endl;

    for (int f=0;f<nodeNum;++f) 
      {
            file << Vnt[f][0] << "\t" << Vnt[f][1] <<  "\t" << Vnt[f][2];
            file << endl;
      }
    if (dimLoc == 3)
      {        

        //write V1 vectors
        file << "VECTORS ";
        file << "v1Vecs ";
        file << "double";
        file << endl;

        for (int f=0;f<nodeNum;++f) 
          {
            file << V1t[f][0] << "\t" << V1t[f][1] << "\t" << V1t[f][2];
            file << endl;
          }

        //write V2 vectors
        file << "VECTORS ";
        file << "v2Vecs ";
        file << "double";
        file << endl;

        for (int f=0;f<nodeNum;++f) 
          {
            file << V2t[f][0] << "\t" << V2t[f][1] << "\t" << V2t[f][2];
            file << endl;
          }
      }

    file.close();    
  }
  
  void surfaceTriangulation::writeTriangulation_aux(string filename,
                                                    vector<double>& Nodes,
                                                    vector<int>& Faces,
                                                    int& dimLoc)
  {
    Screen::MasterInfo("Writing geometry to VTK");
    ofstream file;
    file.open(filename);
    if(file.fail()) Screen::MasterError("Bad file in writeTriangulation function");
    file.precision(20);
    file.setf(ios::fixed);
    file.setf(ios::showpoint);

    int nodeNum = Nodes.size()/3;
    int faceNum = Faces.size()/dimLoc;

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
    file << "POLYDATA ";
    file << endl;
    file << "POINTS ";
    file << nodeNum;
    file << " float";
    file << endl;
    int DIM = 3;
    //writing nodal locations
    for (int n=0;n<nodeNum;++n)
      {
        for (int idir=0;idir<DIM;++idir)
          {
            file << Nodes[n*DIM+idir];
            file << " ";
          }
        file << endl;
      }

    if (dimLoc == 2)
      {
        file << "LINES ";
      }
    else
      {
        file << "POLYGONS "; 
      }
    file << faceNum;
    file << " ";  
    int idata;
    if (dimLoc == 2)
      {
        file << faceNum*3;
        idata = 2;
      }
    else
      {
        file << faceNum*4;
        idata = 3;
      }
    file << endl;    

    //writing face connectivity
    for (int f=0;f<faceNum;++f) 
      {
        file << idata;     
        file << " ";    
        for (int idir=0;idir<dimLoc;++idir)
          {
            file << Faces[f*dimLoc+idir];
            file << " ";
          }
        file << endl;
      }    

    file.close();

  }
}
